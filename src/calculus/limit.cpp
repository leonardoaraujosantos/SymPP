#include <sympp/calculus/limit.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/calculus/series.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/infinity.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

namespace {

// Split a Mul expression of the form num * Pow(den, -1) * ... into
// (numerator, denominator) where denominator is the product of all
// negative-exponent-1 factors. Used after together() so the input is
// already a single fraction.
struct NumDen { Expr num; Expr den; };

[[nodiscard]] NumDen split_after_together(const Expr& e) {
    Expr t = together(e);
    if (t->type_id() != TypeId::Mul) return {t, S::One()};
    std::vector<Expr> num_factors;
    std::vector<Expr> den_factors;
    for (const auto& f : t->args()) {
        // Match any factor with a negative-integer exponent; collect its
        // base raised to the absolute exponent into the denominator. This
        // covers both Pow(d, -1) and Pow(d, -k) cases (e.g. Pow(x, -2)
        // after the canonical (a^m)^n → a^(m·n) fold).
        if (f->type_id() == TypeId::Pow
            && f->args()[1]->type_id() == TypeId::Integer) {
            const auto& z = static_cast<const Integer&>(*f->args()[1]);
            if (z.is_negative()) {
                Expr abs_exp = mul(S::NegativeOne(), f->args()[1]);
                den_factors.push_back(pow(f->args()[0], abs_exp));
                continue;
            }
        }
        num_factors.push_back(f);
    }
    if (den_factors.empty()) return {t, S::One()};
    return {mul(num_factors), mul(den_factors)};
}

[[nodiscard]] bool is_nan(const Expr& e) noexcept {
    return e->type_id() == TypeId::NaN;
}

// True if `e` contains sin/cos/tan of an argument that diverges to infinity —
// the residue of substituting x = ∞ into an unresolved oscillation, e.g. sin(∞),
// cos(∞) + 1, x·sin(∞). Such an expression has no determinate limit (it
// oscillates over a bounded or unbounded range), so the engine reports nan
// rather than leaking the meaningless f(∞). (Genuinely convergent cases —
// sin(x)/x → 0, (x+sin x)/(x+cos x) → 1 — are resolved by earlier rules and
// never reach this guard.)
[[nodiscard]] bool has_oscillating_infinity(const Expr& e) {
    if (!e) return false;
    if (e->type_id() == TypeId::Function) {
        const auto id = static_cast<const Function&>(*e).function_id();
        if ((id == FunctionId::Sin || id == FunctionId::Cos
             || id == FunctionId::Tan)
            && e->args().size() == 1) {
            bool inf = false;
            auto scan = [&](auto&& self, const Expr& x) -> void {
                if (inf) return;
                if (is_infinity(x)) { inf = true; return; }
                for (const auto& a : x->args()) self(self, a);
            };
            scan(scan, e->args()[0]);
            if (inf) return true;
        }
    }
    for (const auto& a : e->args()) {
        if (has_oscillating_infinity(a)) return true;
    }
    return false;
}

// An ∞ appearing as a *proper subexpression* of an arithmetic node — e.g. the
// ∞·(∞+∞·log 4)⁻¹ that L'Hôpital leaves when it differentiates loggamma into
// digamma on a divergent gamma ratio. A genuine limit value is finite, or a bare
// ±∞ / zoo singleton; an ∞ buried inside a Mul/Add/Pow means the form is still
// indeterminate (∞/∞, ∞−∞) and must not be returned as the answer.
[[nodiscard]] bool has_buried_infinity(const Expr& e) {
    if (!e) return false;
    const auto t = e->type_id();
    if (t == TypeId::Infinity || t == TypeId::NegativeInfinity
        || t == TypeId::ComplexInfinity) {
        return false;  // a bare infinity is a valid limit value
    }
    auto contains = [](auto&& self, const Expr& x) -> bool {
        const auto xt = x->type_id();
        if (xt == TypeId::Infinity || xt == TypeId::NegativeInfinity
            || xt == TypeId::ComplexInfinity) {
            return true;
        }
        for (const auto& a : x->args()) {
            if (self(self, a)) return true;
        }
        return false;
    };
    return contains(contains, e);
}

// Recursively combine a rational expression into a single numerator/denominator,
// descending into the bases of integer powers so nested reciprocals flatten —
// (p/q)^(−k) = q^k/p^k. together() stops at a Pow base, so a compound fraction
// like (−x⁻²)⁻¹ stays un-flattened; this clears it (used to rationalise the
// L'Hôpital ratio when a derivative is itself a fraction).
[[nodiscard]] NumDen flatten_fraction(const Expr& e) {
    switch (e->type_id()) {
        case TypeId::Mul: {
            std::vector<Expr> nums, dens;
            for (const auto& a : e->args()) {
                auto r = flatten_fraction(a);
                nums.push_back(r.num);
                dens.push_back(r.den);
            }
            return {mul(std::move(nums)), mul(std::move(dens))};
        }
        case TypeId::Add: {
            std::vector<NumDen> parts;
            parts.reserve(e->args().size());
            for (const auto& a : e->args()) parts.push_back(flatten_fraction(a));
            Expr den = S::One();
            for (const auto& p : parts) den = mul(den, p.den);
            std::vector<Expr> terms;
            terms.reserve(parts.size());
            for (std::size_t i = 0; i < parts.size(); ++i) {
                Expr cof = S::One();
                for (std::size_t j = 0; j < parts.size(); ++j) {
                    if (j != i) cof = mul(cof, parts[j].den);
                }
                terms.push_back(mul(parts[i].num, cof));
            }
            return {add(std::move(terms)), den};
        }
        case TypeId::Pow: {
            const Expr& base = e->args()[0];
            const Expr& ex = e->args()[1];
            if (ex->type_id() == TypeId::Integer) {
                const auto& z = static_cast<const Integer&>(*ex);
                if (z.fits_long()) {
                    auto br = flatten_fraction(base);
                    const long k = z.to_long();
                    if (k < 0) {
                        return {pow(br.den, integer(-k)), pow(br.num, integer(-k))};
                    }
                    if (k > 0) return {pow(br.num, ex), pow(br.den, ex)};
                }
            }
            return {e, S::One()};
        }
        default:
            return {e, S::One()};
    }
}

// L'Hôpital on the fraction num/den. Resolves the 0/0 and ∞/∞ indeterminate
// forms by differentiating top and bottom, re-rationalising the ratio each
// step (so algebraic cancellation, e.g. x²/(x²+x) → x/(x+1), keeps the limit
// tractable). Returns the value, or std::nullopt when it cannot progress.
[[nodiscard]] std::optional<Expr> lhopital_nd(Expr num, Expr den,
                                              const Expr& var,
                                              const Expr& target) {
    // "Infinite-like": a genuine infinity, or nan arising from oo±oo when a
    // polynomial is substituted at the target (its true magnitude is ∞).
    auto infinite_like = [](const Expr& v) noexcept {
        return is_infinity(v) || v->type_id() == TypeId::NaN;
    };
    // Size budget: each L'Hôpital step differentiates num/den, and for a radical
    // integrand the nested radicals grow without bound (the ratio never
    // stabilises). Bail when the expression balloons so limit() returns the
    // unevaluated nan instead of hanging. (sqrt(x²+x) − x needs asymptotic-series
    // / Gruntz machinery, which is deferred; the guard just keeps it terminating.)
    auto node_count = [](const Expr& e) {
        std::size_t n = 0;
        auto rec = [&](auto&& self, const Expr& x) -> void {
            ++n;
            for (const auto& a : x->args()) self(self, a);
        };
        rec(rec, e);
        return n;
    };
    for (int iter = 0; iter < 16; ++iter) {
        if (node_count(num) + node_count(den) > 400) return std::nullopt;
        Expr num_at = simplify(subs(num, var, target));
        Expr den_at = simplify(subs(den, var, target));
        const bool zero_zero = (num_at == S::Zero() && den_at == S::Zero());
        const bool inf_inf = infinite_like(num_at) && infinite_like(den_at);
        if (zero_zero || inf_inf) {
            Expr num_d = diff(num, var);
            Expr den_d = diff(den, var);
            if (den_d == S::Zero()) return std::nullopt;
            // Re-rationalise the new ratio before the next step. together() does
            // not flatten a nested reciprocal like (−x⁻²)⁻¹, so when den_d is itself
            // a fraction (e.g. d/dx(1/x) = −x⁻²) the naive num_d/den_d leaves a
            // negative power in the denominator and the next substitution is nan.
            // flatten_fraction recursively combines into a single numer/denom first;
            // split_after_together then consolidates (cheaper than full simplify).
            // together() cancels common factors (e.g. x²·cos(1/x)/x² → cos(1/x))
            // but can leave a residual nested power when num and den have unrelated
            // denominators (e.g. the atan case: a stray x⁻² stuck in the
            // denominator); flatten_fraction then clears that nesting. Doing only
            // one of the two breaks the other family, so apply both.
            auto ff = flatten_fraction(together(num_d / den_d));
            auto nd = split_after_together(ff.num / ff.den);
            num = std::move(nd.num);
            den = std::move(nd.den);
            continue;
        }
        // Determinate denominator (finite non-zero or an infinity): the value
        // is num/den, which the infinity arithmetic folds (finite/∞ = 0, …).
        if (!is_nan(den_at) && !(den_at == S::Zero())) {
            Expr val = simplify(num_at / den_at);
            if (!is_nan(val)) return val;
        }
        return std::nullopt;
    }
    return std::nullopt;
}

// L'Hôpital on `expr` first split into a single fraction.
[[nodiscard]] std::optional<Expr> lhopital(const Expr& expr, const Expr& var,
                                           const Expr& target) {
    auto nd = split_after_together(expr);
    return lhopital_nd(nd.num, nd.den, var, target);
}

[[nodiscard]] Expr limit_impl(const Expr& expr, const Expr& var,
                              const Expr& target, int depth);

// exp() of a (possibly infinite) limit value: exp(+oo)=oo, exp(-oo)=0.
[[nodiscard]] Expr exp_of_limit(const Expr& v) {
    if (v->type_id() == TypeId::Infinity) return S::Infinity();
    if (v->type_id() == TypeId::NegativeInfinity) return S::Zero();
    return simplify(exp(v));
}

// Power indeterminate forms 1^∞, ∞^0, 0^0 via a^b = exp(b·log a).
[[nodiscard]] std::optional<Expr> try_power_form(const Expr& expr,
                                                 const Expr& var,
                                                 const Expr& target,
                                                 int depth) {
    if (expr->type_id() != TypeId::Pow) return std::nullopt;
    const Expr& base = expr->args()[0];
    const Expr& ex = expr->args()[1];
    Expr bl = limit_impl(base, var, target, depth + 1);
    Expr el = limit_impl(ex, var, target, depth + 1);
    const bool one_inf = (bl == S::One() && is_infinity(el));
    const bool inf_zero = (is_infinity(bl) && el == S::Zero());
    const bool zero_zero = (bl == S::Zero() && el == S::Zero());
    if (!(one_inf || inf_zero || zero_zero)) return std::nullopt;
    // For the 1^∞ form use log(base) ~ (base − 1) as base → 1, so
    // lim ex·log(base) = lim ex·(base − 1). This avoids needing the Taylor series
    // of log of a composite base — the series engine leaves log(sin x / x)
    // unexpanded — letting e.g. (sin x / x)^(1/x²) resolve to e^(−1/6). The other
    // forms (∞^0, 0^0) genuinely need log(base) and keep it.
    Expr rate = one_inf ? mul(ex, add(base, S::NegativeOne()))
                        : mul(ex, log(base));
    Expr inner = limit_impl(rate, var, target, depth + 1);
    if (is_nan(inner)) return std::nullopt;
    return exp_of_limit(inner);
}

// Product indeterminate form 0·∞: rewrite the factors that tend to ∞ into the
// denominator (as reciprocals) and resolve the resulting quotient.
[[nodiscard]] std::optional<Expr> try_product_form(const Expr& expr,
                                                   const Expr& var,
                                                   const Expr& target,
                                                   int depth) {
    if (expr->type_id() != TypeId::Mul) return std::nullopt;
    // A vanishing factor log(g) with g → 1 is the 0-side of a 0·∞ product, but
    // its naive endpoint subs gives nan when g is a ratio like (x+1)/x (∞/∞) —
    // which derails L'Hôpital's 0/0 test. Replace it by its leading asymptotic
    // log(g) ~ (g − 1): log(g)/(g − 1) → 1 as g → 1, so the product's limit is
    // preserved. Mirrors the 1^∞ power-form trick (log(base) ~ base − 1).
    // Combine a sum Σ cᵢ·log(gᵢ) into the single argument ∏ gᵢ^cᵢ (so that the
    // whole sum equals log(∏ gᵢ^cᵢ)). Returns nullopt unless every term is a
    // numeric-scaled log — simplify() does not fold a bare log difference.
    auto combine_log_args = [](const Expr& sum) -> std::optional<Expr> {
        std::vector<Expr> factors;
        for (const auto& term : sum->args()) {
            Expr coeff = S::One();
            Expr core = term;
            if (term->type_id() == TypeId::Mul) {
                std::vector<Expr> rest;
                for (const auto& g : term->args()) {
                    if (is_number(g)) coeff = mul(coeff, g);
                    else rest.push_back(g);
                }
                if (rest.size() != 1) return std::nullopt;
                core = rest[0];
            }
            if (core->type_id() != TypeId::Function) return std::nullopt;
            const auto& fn = static_cast<const Function&>(*core);
            if (fn.function_id() != FunctionId::Log || fn.args().size() != 1)
                return std::nullopt;
            factors.push_back(pow(fn.args()[0], coeff));
        }
        if (factors.empty()) return std::nullopt;
        return mul(std::move(factors));
    };
    // A vanishing log factor with argument g → 1, possibly written as a sum of
    // logs (log(x+1) − log(x) = log((x+1)/x)). Replace it by log(g) ~ (g − 1):
    // the leading asymptotic preserves the 0·∞ product's limit while avoiding
    // the nan that g's ∞/∞ endpoint subs would feed L'Hôpital.
    auto soften_log_to_one = [&](const Expr& f) -> Expr {
        std::optional<Expr> arg;
        if (f->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f);
            if (fn.function_id() == FunctionId::Log && fn.args().size() == 1)
                arg = fn.args()[0];
        } else if (f->type_id() == TypeId::Add) {
            arg = combine_log_args(f);
        }
        if (arg && limit_impl(*arg, var, target, depth + 1) == S::One()) {
            // Reduce (g − 1) to lowest terms (e.g. (x+1)/x − 1 → 1/x) so its
            // endpoint subs is clean rather than another ∞/∞ nan.
            return simplify(add(*arg, S::NegativeOne()));
        }
        return f;
    };
    std::vector<Expr> zero_f, inf_f, finite_f;
    for (const auto& f : expr->args()) {
        Expr lf = limit_impl(f, var, target, depth + 1);
        if (is_nan(lf)) return std::nullopt;
        if (is_infinity(lf)) inf_f.push_back(f);
        else if (lf == S::Zero()) zero_f.push_back(soften_log_to_one(f));
        else finite_f.push_back(f);
    }
    if (zero_f.empty() || inf_f.empty()) return std::nullopt;

    auto prod = [](std::vector<Expr> v) -> Expr {
        return v.empty() ? Expr{S::One()} : mul(std::move(v));
    };
    // Reciprocal of a factor, but exp-aware: 1/exp(g) = exp(−g). Keeping the
    // exponential as a single Function (rather than exp(g)^(−1)) means L'Hôpital
    // can hold it in the denominator across iterations instead of flipping it
    // back into the numerator each step (which stalls x^n·e^(−x)).
    auto recip_one = [](const Expr& f) -> Expr {
        if (f->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f);
            if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1) {
                return exp(mul(S::NegativeOne(), fn.args()[0]));
            }
        }
        return pow(f, S::NegativeOne());
    };
    auto recips = [&](const std::vector<Expr>& v) {
        std::vector<Expr> r;
        r.reserve(v.size());
        for (const auto& f : v) r.push_back(recip_one(f));
        return r;
    };
    Expr finite = prod(finite_f);

    // An f·g with f→0, g→∞ is an indeterminate 0·∞. Recast it as a quotient and
    // apply L'Hôpital. Two arrangements; try the 0/0 one first (it preserves the
    // existing behaviour, e.g. x·sin(1/x)), then the ∞/∞ one — needed for
    // x·e^(-x), where differentiating the 0/0 form only raises the polynomial
    // degree, whereas ∞/∞ (x / e^x) collapses in one step.
    {
        // Arrangement A — 0/0: (Πzeros · finite) / (Π 1/∞-factors).
        std::vector<Expr> nf = zero_f;
        if (!(finite == S::One())) nf.push_back(finite);
        if (auto r = lhopital_nd(prod(std::move(nf)), prod(recips(inf_f)), var,
                                 target);
            r && !is_nan(*r)) {
            return r;
        }
    }
    // Arrangement B — ∞/∞: (Π∞-factors · finite) / (Π 1/zeros).
    std::vector<Expr> nf = inf_f;
    if (!(finite == S::One())) nf.push_back(finite);
    return lhopital_nd(prod(std::move(nf)), prod(recips(zero_f)), var, target);
}

// A product/ratio of constant-base exponentials sharing one var-monomial in the
// exponent — e.g. 2^x/3^x = 2^x·3^(−x), or exp(x)/exp(2x) = exp(x)·exp(−2x).
// Evaluated factor by factor these are an ∞·0 (or 0·∞) indeterminate that
// L'Hôpital cannot crack (differentiating reproduces the same form), so the
// product path returns nan. But the combined rate is a single exponential:
// ∏ bᵢ^(cᵢ·m) · ∏ exp(dⱼ·m) = exp(m·log B), with B = ∏ bᵢ^cᵢ · e^(Σdⱼ) a concrete
// positive constant. The single-exponential path then resolves it (B<1 → 0,
// B>1 → ∞, B=1 → 1 as m → +∞). Fixes lim 2^x/3^x = 0 (was nan).
[[nodiscard]] std::optional<Expr> try_exponential_product(const Expr& expr,
                                                          const Expr& var,
                                                          const Expr& target,
                                                          int depth) {
    if (!is_infinity(target)) return std::nullopt;
    if (expr->type_id() != TypeId::Mul) return std::nullopt;

    // Split an exponent e into (coeff, monomial) with e = coeff·monomial, coeff a
    // number and monomial var-dependent; nullopt if e has no such factorization.
    auto split = [&](const Expr& e) -> std::optional<std::pair<Expr, Expr>> {
        if (!has(e, var)) return std::nullopt;
        if (e->type_id() == TypeId::Mul) {
            std::vector<Expr> nums, rest_f;
            for (const auto& f : e->args()) {
                if (f->type_id() == TypeId::Integer
                    || f->type_id() == TypeId::Rational) {
                    nums.push_back(f);
                } else {
                    rest_f.push_back(f);
                }
            }
            Expr m = rest_f.empty() ? Expr{S::One()} : mul(rest_f);
            if (!has(m, var)) return std::nullopt;
            Expr c = nums.empty() ? Expr{S::One()} : mul(nums);
            return std::make_pair(c, m);
        }
        return std::make_pair(Expr{S::One()}, e);
    };

    auto is_exp = [](const Expr& f) -> const Function* {
        if (f->type_id() != TypeId::Function) return nullptr;
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1) {
            return &fn;
        }
        return nullptr;
    };

    Expr shared_mono;             // the common var-monomial in every exponent
    Expr base_prod = S::One();    // ∏ bᵢ^cᵢ for the constant-base power factors
    Expr exp_csum = S::Zero();    // Σ dⱼ for the exp(dⱼ·m) factors
    std::vector<Expr> rest;       // remaining (must be polynomial in var) factors
    int combined = 0;
    auto take_exp_arg = [&](const Expr& arg) -> bool {
        auto sp = split(arg);
        if (!sp || (shared_mono && !(shared_mono == sp->second))) return false;
        shared_mono = sp->second;
        exp_csum = add(exp_csum, sp->first);
        ++combined;
        return true;
    };
    for (const auto& f : expr->args()) {
        if (f->type_id() == TypeId::Pow) {
            const Expr& b = f->args()[0];
            const Expr& e = f->args()[1];
            // Constant positive base ^ var-dependent exponent: bᵢ^(cᵢ·m).
            if (!has(b, var) && is_positive(b) == true && has(e, var)) {
                if (auto sp = split(e);
                    sp && (!shared_mono || shared_mono == sp->second)) {
                    shared_mono = sp->second;
                    base_prod = mul(base_prod, pow(b, sp->first));
                    ++combined;
                    continue;
                }
            }
            // exp(g)^k with numeric k = exp(k·g) — this is how a/exp(g) and
            // exp(g)^2 surface after canonicalization (1/exp(2x) is Pow(exp,−1)).
            if (const Function* fn = is_exp(b);
                fn && (e->type_id() == TypeId::Integer
                       || e->type_id() == TypeId::Rational)) {
                if (take_exp_arg(mul(e, fn->args()[0]))) continue;
            }
        }
        if (const Function* fn = is_exp(f); fn && has(fn->args()[0], var)) {
            if (take_exp_arg(fn->args()[0])) continue;
        }
        // Anything else is kept as a residual factor; it must be polynomial in var
        // (checked below) so growth dominance against the exponential is decidable.
        rest.push_back(f);
    }
    if (combined < 2 && !(combined >= 1 && !rest.empty())) return std::nullopt;
    if (!shared_mono) return std::nullopt;

    // The residual factors must form a polynomial in var: the exponential's growth
    // class strictly dominates any polynomial, so a decaying exponential drives
    // the whole product to 0 and a growing one to ±∞ regardless of polynomial
    // degree. A non-polynomial residual (another transcendental) would break that.
    Expr rest_prod = rest.empty() ? Expr{S::One()} : mul(rest);
    if (!rest.empty()) {
        try {
            Poly rp(expand(rest_prod), var);
            for (const auto& cc : rp.coeffs()) {
                if (has(cc, var)) return std::nullopt;
            }
        } catch (const std::exception&) {
            return std::nullopt;  // residual not polynomial in var
        }
    }

    // The exponential part is B^m with B = ∏ bᵢ^cᵢ · e^(Σdⱼ) a concrete positive
    // constant and m = shared_mono. As m → ±∞ its behaviour is set by B vs 1 —
    // decide it from the sign of B−1 directly. (Going through exp(m·log B) instead
    // lets simplify distribute log(2·e⁻³) → log 2 − 3, which the engine then can't
    // sign — so we avoid symbolic logs here.)
    Expr B = simplify(mul(base_prod, exp(exp_csum)));
    // Sign of B−1, falling back to a numeric evaluation: the assumption system
    // can sign a rational like 2/3−1 but not exp(−1)−1, which surfaces for the
    // exp(·)/exp(·) and mixed base·exp cases.
    auto sign_of = [&](const Expr& e) -> int {
        if (is_positive(e) == true) return 1;
        if (is_negative(e) == true) return -1;
        Expr ef = evalf(e, 30);
        if (is_positive(ef) == true) return 1;
        if (is_negative(ef) == true) return -1;
        return 0;
    };
    Expr Lm = limit_impl(shared_mono, var, target, depth + 1);
    int dir;
    if (Lm->type_id() == TypeId::Infinity) {
        dir = 1;
    } else if (Lm->type_id() == TypeId::NegativeInfinity) {
        dir = -1;
    } else {
        return std::nullopt;  // monomial does not diverge — leave to other paths
    }

    if (B == S::One()) {
        // B^m ≡ 1: the limit is just the residual's. (Reachable only with a
        // residual, since a bare B=1 product collapses to 1 upstream.)
        return rest.empty() ? Expr{S::One()}
                            : limit_impl(rest_prod, var, target, depth + 1);
    }
    const int s = sign_of(simplify(add(B, S::NegativeOne())));
    if (s == 0) return std::nullopt;  // B vs 1 undecidable

    const bool exp_decays = (dir * s < 0);  // B^m → 0
    if (exp_decays) {
        // Geometric decay dominates any polynomial residual → 0.
        return S::Zero();
    }
    // B^m → +∞ (B > 0 so the exponential itself is positive); the sign of the
    // whole product is the sign of the polynomial residual's divergence.
    Expr Lrest = limit_impl(rest_prod, var, target, depth + 1);
    if (Lrest == S::Zero()) return S::Zero();  // residual ≡ 0
    if (Lrest->type_id() == TypeId::Infinity || is_positive(Lrest) == true) {
        return S::Infinity();
    }
    if (Lrest->type_id() == TypeId::NegativeInfinity
        || is_negative(Lrest) == true) {
        return S::NegativeInfinity();
    }
    return std::nullopt;  // residual sign undecidable
}

// Merge every exponential factor of a product into a single exp:
// rest · ∏ exp(gᵢ)^{kᵢ} = rest · exp(Σ kᵢ gᵢ). A product/ratio of exponentials
// with *different* exponent monomials — exp(x²)·exp(−2x), exp(x²)/exp(x)²,
// x²·exp(x)/exp(x²) — substitutes factor-by-factor to ∞·0 = nan and is not
// handled by try_exponential_product, which only combines exponentials sharing
// one monomial. Merging the exponents leaves a single exp(g) whose limit the
// continuity rule (and the poly×exp machinery) resolves on the re-take. Requires
// ≥2 exponential factors so the merge strictly reduces their count and the
// recursion terminates. The identity exp(a)·exp(b) = exp(a+b) holds everywhere,
// so this is valid at any target.
[[nodiscard]] std::optional<Expr> try_exp_combine(const Expr& expr,
                                                  const Expr& var,
                                                  const Expr& target,
                                                  int depth) {
    if (expr->type_id() != TypeId::Mul) return std::nullopt;
    // (g, k) such that the factor equals exp(g)^k, else nullopt.
    auto as_exp_pow =
        [](const Expr& f) -> std::optional<std::pair<Expr, Expr>> {
        if (f->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f);
            if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1) {
                return std::make_pair(fn.args()[0], Expr{S::One()});
            }
        }
        if (f->type_id() == TypeId::Pow && is_number(f->args()[1])) {
            const Expr& b = f->args()[0];
            if (b->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*b);
                if (fn.function_id() == FunctionId::Exp
                    && fn.args().size() == 1) {
                    return std::make_pair(fn.args()[0], f->args()[1]);
                }
            }
        }
        return std::nullopt;
    };
    Expr exponent = S::Zero();
    std::vector<Expr> rest;
    int exp_count = 0;
    for (const auto& f : expr->args()) {
        if (auto ge = as_exp_pow(f)) {
            exponent = add(exponent, mul(ge->second, ge->first));
            ++exp_count;
        } else {
            rest.push_back(f);
        }
    }
    if (exp_count < 2) return std::nullopt;
    rest.push_back(exp(exponent));
    return limit_impl(mul(std::move(rest)), var, target, depth + 1);
}

// Sign of `expr` evaluated at `point`: +1 / −1, or 0 when indeterminate (a
// non-real value, an unsigned/indeterminate infinity, or no definite sign).
[[nodiscard]] int side_sign_at(const Expr& expr, const Expr& var,
                               const Expr& point) {
    Expr v = simplify(subs(expr, var, point));
    if (v->type_id() == TypeId::Infinity) return 1;
    if (v->type_id() == TypeId::NegativeInfinity) return -1;
    if (is_infinity(v) || v->type_id() == TypeId::NaN) return 0;
    Expr f = evalf(v, 30);
    if (is_positive(f) == std::optional<bool>{true}) return 1;
    if (is_negative(f) == std::optional<bool>{true}) return -1;
    return 0;
}

// Direct substitution at a finite pole yields zoo (the unsigned 1/0). Recover
// the signed infinity by sampling the sign of the expression just either side of
// the target: equal signs ⇒ ±oo (1/x² → +oo, −1/x² → −oo), opposite signs ⇒ the
// limit is genuinely two-sided and stays zoo (1/x). Inconclusive samples (a
// non-real value, or no definite sign) also keep zoo.
[[nodiscard]] std::optional<Expr> signed_pole(const Expr& expr, const Expr& var,
                                              const Expr& target) {
    if (!is_number(target)) return std::nullopt;
    Expr h = pow(integer(10), integer(-6));  // 1e-6
    int s_right = side_sign_at(expr, var, add(target, h));
    int s_left = side_sign_at(expr, var, add(target, mul(S::NegativeOne(), h)));
    if (s_right == 0 || s_left == 0 || s_right != s_left) return std::nullopt;
    return s_right > 0 ? S::Infinity() : S::NegativeInfinity();
}

// An expression containing sign(g) or abs(g) with g → 0 at the target is
// discontinuous there: direct substitution uses sign(0)=0 (and the abs is hidden
// inside a quotient like |x|/x), which is wrong. Resolve it from the actual
// one-sided behavior: in the right- and left-neighborhoods replace each sign(g)
// by its sampled sign s∈{±1} and each abs(g) by g·s (since |g| = g·sign g), then
// limit each now-continuous version. Agreeing sides give the limit; disagreeing
// sides mean the two-sided limit does not exist (nan). Returns nullopt when there
// is no such vanishing node or a sample is inconclusive, leaving the normal path
// untouched. This fixes sign(x) → nan, sign(x²) → 1, and |x|/x → nan.
[[nodiscard]] std::optional<Expr> resolve_sign_limit(const Expr& expr,
                                                     const Expr& var,
                                                     const Expr& target,
                                                     int depth) {
    if (!is_number(target)) return std::nullopt;
    std::vector<Expr> nodes;  // sign(g) / abs(g) with g → 0
    auto collect = [&](auto&& self, const Expr& e) -> void {
        if (e->type_id() == TypeId::Function) {
            const FunctionId id = static_cast<const Function&>(*e).function_id();
            if ((id == FunctionId::Sign || id == FunctionId::Abs)
                && limit_impl(e->args()[0], var, target, depth + 1)
                       == S::Zero()) {
                if (std::none_of(nodes.begin(), nodes.end(),
                                 [&](const Expr& s) { return s == e; })) {
                    nodes.push_back(e);
                }
            }
        }
        for (const auto& a : e->args()) self(self, a);
    };
    collect(collect, expr);
    if (nodes.empty()) return std::nullopt;
    const Expr h = pow(integer(10), integer(-6));
    auto neighborhood_limit = [&](const Expr& point) -> std::optional<Expr> {
        ExprMap<Expr> m;
        for (const auto& node : nodes) {
            const Expr& g = node->args()[0];
            int sg = side_sign_at(g, var, point);
            if (sg == 0) return std::nullopt;  // inconclusive sample
            const FunctionId id =
                static_cast<const Function&>(*node).function_id();
            // sign(g) → s;  abs(g) = g·sign(g) → g·s.
            m.emplace(node, id == FunctionId::Sign
                                ? Expr{integer(sg)}
                                : mul(g, integer(sg)));
        }
        return limit_impl(xreplace(expr, m), var, target, depth + 1);
    };
    auto right = neighborhood_limit(add(target, h));
    auto left = neighborhood_limit(add(target, mul(S::NegativeOne(), h)));
    if (!right || !left) return std::nullopt;
    if (*right == *left) return *right;
    return S::NaN();
}

// Limit of a rational function N(x)/D(x) at ±∞ by degree comparison and the
// leading-coefficient ratio. Direct substitution is unreliable at ∞ (∞·0, ∞/∞),
// and L'Hôpital mishandles symbolic (var-free) coefficients, so resolve these
// before either: (x+a)/x → 1, a·x/(x+1) → a, (x²+a)/(x+1) → ∞. Returns nullopt
// unless `expr` is N/D with a var-dependent denominator and polynomial, var-free
// coefficients on both.
[[nodiscard]] std::optional<Expr> rational_limit_at_infinity(
    const Expr& expr, const Expr& var, const Expr& target) {
    Expr tg = together(expr);
    Expr num = S::One();
    Expr den;
    auto den_of = [](const Expr& f) -> std::optional<Expr> {
        if (f->type_id() == TypeId::Pow
            && f->args()[1]->type_id() == TypeId::Integer) {
            const auto& z = static_cast<const Integer&>(*f->args()[1]);
            if (z.is_negative() && z.fits_long()) {
                return pow(f->args()[0], integer(-z.to_long()));
            }
        }
        return std::nullopt;
    };
    if (auto d0 = den_of(tg)) {
        den = *d0;
    } else if (tg->type_id() == TypeId::Mul) {
        std::vector<Expr> nf, df;
        for (const auto& f : tg->args()) {
            if (auto d = den_of(f)) df.push_back(*d);
            else nf.push_back(f);
        }
        if (!df.empty()) {
            den = mul(std::move(df));
            num = nf.empty() ? Expr{S::One()} : mul(std::move(nf));
        }
    }
    if (!den || !has(den, var)) return std::nullopt;
    try {
        Poly np(expand(num), var);
        Poly dp(expand(den), var);
        for (const auto& cc : np.coeffs()) {
            if (has(cc, var)) return std::nullopt;
        }
        for (const auto& cc : dp.coeffs()) {
            if (has(cc, var)) return std::nullopt;
        }
        const long dn = static_cast<long>(np.degree());
        const long dd = static_cast<long>(dp.degree());
        if (dn < dd) return S::Zero();
        Expr ratio = simplify(np.leading_coeff() / dp.leading_coeff());
        if (dn == dd) return ratio;
        Expr lim =
            simplify(subs(ratio * pow(var, integer(dn - dd)), var, target));
        if (!is_nan(lim)) return lim;
    } catch (const std::exception&) {
    }
    return std::nullopt;
}

// Rewrite factorial(u) → gamma(u+1) throughout, so the gamma-ratio
// simplification can collapse mixed forms like gamma(x)/factorial(x) = 1/x.
[[nodiscard]] Expr normalize_factorial(const Expr& e) {
    ExprMap<Expr> m;
    auto scan = [&](auto&& self, const Expr& x) -> void {
        if (x->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*x);
            if (fn.function_id() == FunctionId::Factorial
                && fn.args().size() == 1) {
                m.emplace(x, gamma(add(fn.args()[0], S::One())));
            }
        }
        for (const auto& a : x->args()) self(self, a);
    };
    scan(scan, e);
    if (m.empty()) return e;
    return xreplace(e, m);
}

// Rewrite reciprocal trig/hyperbolic functions as sin/cos ratios:
//   cot→cos/sin, csc→1/sin, sec→1/cos, coth→cosh/sinh, csch→1/sinh,
//   sech→1/cosh. The limit machinery (direct substitution, L'Hôpital) handles
//   sin/cos but treats cot/csc/… as opaque, so e.g. x·cot(x) → nan instead of 1.
[[nodiscard]] Expr rewrite_reciprocal_trig(const Expr& e) {
    ExprMap<Expr> m;
    auto scan = [&](auto&& self, const Expr& x) -> void {
        if (x->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*x);
            if (fn.args().size() == 1) {
                const Expr& u = fn.args()[0];
                const Expr inv_sin = pow(sin(u), S::NegativeOne());
                const Expr inv_cos = pow(cos(u), S::NegativeOne());
                const Expr inv_sinh = pow(sinh(u), S::NegativeOne());
                const Expr inv_cosh = pow(cosh(u), S::NegativeOne());
                switch (fn.function_id()) {
                    case FunctionId::Cot:  m.emplace(x, mul(cos(u), inv_sin)); break;
                    case FunctionId::Csc:  m.emplace(x, inv_sin); break;
                    case FunctionId::Sec:  m.emplace(x, inv_cos); break;
                    case FunctionId::Coth: m.emplace(x, mul(cosh(u), inv_sinh)); break;
                    case FunctionId::Csch: m.emplace(x, inv_sinh); break;
                    case FunctionId::Sech: m.emplace(x, inv_cosh); break;
                    default: break;
                }
            }
        }
        for (const auto& a : x->args()) self(self, a);
    };
    scan(scan, e);
    if (m.empty()) return e;
    return xreplace(e, m);
}

// Rewrite hyperbolic functions of a var-dependent argument to their exponential
// definitions, so combinations the closed forms hide — sinh u + cosh u = eᵘ,
// cosh u − sinh u = e⁻ᵘ, sinh u/cosh u = tanh u → 1 — collapse to elementary
// exponentials the asymptotic machinery already resolves. Used only at ±∞, where
// this is the right canonicalization (eᵘ vs e⁻ᵘ separate cleanly by growth).
[[nodiscard]] Expr rewrite_hyperbolic_exp(const Expr& e, const Expr& var,
                                          const Expr& target, int depth) {
    ExprMap<Expr> m;
    auto scan = [&](auto&& self, const Expr& x) -> void {
        if (x->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*x);
            // Only rewrite a hyperbolic whose argument itself diverges: there the
            // eᵘ / e⁻ᵘ split is the right asymptotics. For a vanishing argument
            // (e.g. tanh(e⁻ˣ), arg → 0) the leading-term/small-angle rule applies,
            // and the exponential form would defeat it.
            if (fn.args().size() == 1 && has(fn.args()[0], var)
                && is_infinity(
                    limit_impl(fn.args()[0], var, target, depth + 1))) {
                const Expr& u = fn.args()[0];
                const Expr eu = exp(u);
                const Expr enu = exp(mul(S::NegativeOne(), u));
                const Expr sum = add(eu, enu);                       // eᵘ + e⁻ᵘ
                const Expr diff = add(eu, mul(S::NegativeOne(), enu));  // eᵘ − e⁻ᵘ
                switch (fn.function_id()) {
                    case FunctionId::Sinh:
                        m.emplace(x, mul(rational(1, 2), diff)); break;
                    case FunctionId::Cosh:
                        m.emplace(x, mul(rational(1, 2), sum)); break;
                    case FunctionId::Tanh:
                        m.emplace(x, mul(diff, pow(sum, S::NegativeOne())));
                        break;
                    case FunctionId::Coth:
                        m.emplace(x, mul(sum, pow(diff, S::NegativeOne())));
                        break;
                    case FunctionId::Sech:
                        m.emplace(x, mul(integer(2),
                                         pow(sum, S::NegativeOne())));
                        break;
                    case FunctionId::Csch:
                        m.emplace(x, mul(integer(2),
                                         pow(diff, S::NegativeOne())));
                        break;
                    default: break;
                }
            }
        }
        for (const auto& a : x->args()) self(self, a);
    };
    scan(scan, e);
    if (m.empty()) return e;
    return xreplace(e, m);
}

// Flatten (eᵃ)ⁿ → eⁿᵃ everywhere. After a hyperbolic ratio is rewritten and the
// denominator expanded to a single e⁻ᵘ, its reciprocal is (e⁻ᵘ)⁻¹, which the limit
// engine otherwise reads as 0⁻¹ = zoo rather than eᵘ.
[[nodiscard]] Expr flatten_exp_powers(const Expr& e) {
    ExprMap<Expr> m;
    auto scan = [&](auto&& self, const Expr& x) -> void {
        if (x->type_id() == TypeId::Pow && !m.count(x)
            && x->args()[0]->type_id() == TypeId::Function
            && static_cast<const Function&>(*x->args()[0]).function_id()
                   == FunctionId::Exp) {
            m.emplace(x, exp(mul(x->args()[1], x->args()[0]->args()[0])));
        }
        for (const auto& a : x->args()) self(self, a);
    };
    scan(scan, e);
    if (m.empty()) return e;
    return xreplace(e, m);
}

// Number of gamma/factorial applications anywhere in e.
[[nodiscard]] int count_gamma_factorial(const Expr& e) {
    int n = 0;
    if (e->type_id() == TypeId::Function) {
        const auto id = static_cast<const Function&>(*e).function_id();
        if (id == FunctionId::Gamma || id == FunctionId::Factorial) ++n;
    }
    for (const auto& a : e->args()) n += count_gamma_factorial(a);
    return n;
}

// Asymptotic growth at +∞ obeys the strict hierarchy
//   factorial/gamma  ≫  exponential  ≫  polynomial  ≫  logarithm.
// classify_growth maps one multiplicative factor f^s (s = ±1 the side it sits
// on: +1 numerator, −1 denominator) to {rank, direction}, where direction is the
// sign of f's contribution to log|expr| (a factor → 0 contributes −). Returns
// nullopt for a factor that does not fit the hierarchy, so the caller bails.
struct Growth {
    int rank;  // 4 gamma/factorial, 3 exp, 2 polynomial, 1 log, 0 bounded
    int dir;   // +1 pushes magnitude up, −1 down (toward 0), 0 bounded
};
[[nodiscard]] std::optional<Growth> classify_growth(const Expr& factor,
                                                    const Expr& var,
                                                    int depth) {
    // Peel a numeric exponent: base^e contributes sign(e) to the side.
    Expr base = factor;
    int side = 1;
    if (factor->type_id() == TypeId::Pow) {
        const Expr& b = factor->args()[0];
        const Expr& e = factor->args()[1];
        if (e->type_id() == TypeId::Integer) {
            const auto& z = static_cast<const Integer&>(*e);
            if (!z.fits_long()) return std::nullopt;
            const long ev = z.to_long();
            if (ev == 0) return Growth{0, 0};
            side = (ev > 0) ? 1 : -1;
            base = b;
        } else if (!has(b, var) && is_positive(b) == true && has(e, var)) {
            // B^E with positive constant base B and a var-dependent exponent: an
            // exponential a^x. Direction = sign(E)·sign(log B).
            Expr le = limit_impl(e, var, S::Infinity(), depth + 1);
            const int sE = (le == S::Infinity()) ? 1
                           : (le == S::NegativeInfinity()) ? -1
                                                           : 0;
            if (sE == 0) return Growth{0, 0};  // E bounded ⇒ B^E bounded
            int sB;
            if (is_positive(simplify(add(b, S::NegativeOne()))) == true) {
                sB = 1;  // B > 1
            } else if (is_positive(simplify(add(S::One(),
                                                mul(S::NegativeOne(), b))))
                       == true) {
                sB = -1;  // 0 < B < 1
            } else {
                return std::nullopt;  // B vs 1 undecidable
            }
            return Growth{3, sB * sE};
        } else {
            return std::nullopt;  // x^x and friends — not in the hierarchy
        }
    }
    if (!has(base, var)) return Growth{0, 0};  // bounded constant factor
    auto inner_to_pos_inf = [&](const Expr& u) {
        return is_infinity(limit_impl(u, var, S::Infinity(), depth + 1))
               && limit_impl(u, var, S::Infinity(), depth + 1) == S::Infinity();
    };
    if (base->type_id() == TypeId::Function) {
        const auto id = static_cast<const Function&>(*base).function_id();
        const Expr& u = base->args()[0];
        if (id == FunctionId::Gamma || id == FunctionId::Factorial) {
            if (!inner_to_pos_inf(u)) return std::nullopt;  // gamma(−∞): wild
            return Growth{4, side};
        }
        if (id == FunctionId::Exp) {
            Expr lu = limit_impl(u, var, S::Infinity(), depth + 1);
            if (lu == S::Infinity()) return Growth{3, side};
            if (lu == S::NegativeInfinity()) return Growth{3, -side};
            return std::nullopt;
        }
        if (id == FunctionId::Log) {
            if (!inner_to_pos_inf(u)) return std::nullopt;
            return Growth{1, side};
        }
        return std::nullopt;
    }
    // var itself, or any polynomial in var that diverges: rank 2.
    if (base == var) return Growth{2, side};
    try {
        Poly p(expand(base), var);
        if (p.degree() >= 1) {
            for (const auto& cc : p.coeffs())
                if (has(cc, var)) return std::nullopt;
            return Growth{2, side};
        }
    } catch (const std::exception&) {
    }
    return std::nullopt;
}

// Limit at +∞ decided purely by the dominant growth class: if a single factor
// has a strictly higher rank than every other (no same-rank competitor with the
// opposite direction), it sets the limit to +∞ or 0. Mixed top ranks (e.g.
// gamma(2x)/gamma(x)) are left to other methods. Only attempted when a
// gamma/factorial is present, so ordinary limits are unaffected.
[[nodiscard]] std::optional<Expr> gamma_growth_limit(const Expr& expr,
                                                     const Expr& var,
                                                     int depth) {
    Expr t = together(expr);
    std::vector<Expr> factors;
    if (t->type_id() == TypeId::Mul) {
        for (const auto& f : t->args()) factors.push_back(f);
    } else {
        factors.push_back(t);
    }
    int top_rank = 0;
    int top_dir_sum = 0;
    bool top_conflict = false;
    for (const auto& f : factors) {
        auto g = classify_growth(f, var, depth);
        if (!g) return std::nullopt;  // unrecognized factor — cannot conclude
        if (g->rank == 0) continue;
        if (g->rank > top_rank) {
            top_rank = g->rank;
            top_dir_sum = g->dir;
            top_conflict = false;
        } else if (g->rank == top_rank) {
            if (g->dir != 0 && top_dir_sum != 0
                && ((g->dir > 0) != (top_dir_sum > 0))) {
                top_conflict = true;
            }
            top_dir_sum += g->dir;
        }
    }
    if (top_rank == 0 || top_conflict || top_dir_sum == 0) return std::nullopt;
    return top_dir_sum > 0 ? S::Infinity() : S::Zero();
}

// Balanced gamma-ratio asymptotic at +∞. By Stirling, log Γ(x+a) =
// (x+a−½)·log x − x + ½·log 2π + o(1), so for a product C · ∏ᵢ Γ(x+aᵢ)^{pᵢ} · x^q
// with Σ pᵢ = 0 (a ratio of gammas) the x·log x and x terms cancel and
//   ∏ᵢ Γ(x+aᵢ)^{pᵢ} ~ x^{Σ pᵢ aᵢ},
// hence the whole expression ~ C · x^{q + Σ pᵢ aᵢ}. The limit is C (exponent 0),
// sign(C)·∞ (positive), or 0 (negative). This resolves Γ(x+½)/Γ(x)/√x → 1, on
// which the growth-rank test must abstain (same rank, opposite direction) and a
// later path returns a wrong 0. Each gamma argument must be a pure shift x + aᵢ
// (so Γ(2x)/Γ(x), where the argument is not a shift, is left to other methods)
// and every pᵢ, aᵢ, q numeric; otherwise nullopt.
[[nodiscard]] std::optional<Expr> gamma_ratio_asymptotic(const Expr& expr,
                                                         const Expr& var) {
    // Flatten into multiplicative factors, distributing a power of a product
    // (∏ gᵢ)^p → ∏ gᵢ^p with p numeric. `together` leaves a ratio's denominator
    // as one combined Pow — Γ(n+½)/(√n·Γ(n)) becomes (√n·Γ(n))⁻¹·Γ(n+½) — which
    // would otherwise read as a single unrecognized factor and abort the rule.
    std::vector<Expr> factors;
    auto add_factor = [&](auto&& self, const Expr& f) -> void {
        if (f->type_id() == TypeId::Mul) {
            for (const auto& g : f->args()) self(self, g);
            return;
        }
        if (f->type_id() == TypeId::Pow && is_number(f->args()[1])
            && f->args()[0]->type_id() == TypeId::Mul) {
            for (const auto& g : f->args()[0]->args()) {
                self(self, pow(g, f->args()[1]));
            }
            return;
        }
        factors.push_back(f);
    };
    add_factor(add_factor, expr);
    Expr C = S::One(), q_total = S::Zero();
    Expr sum_p = S::Zero(), sum_pa = S::Zero();
    Expr log_rate = S::Zero();  // ρ = Σ kᵢ·log cᵢ from constant-base c^(kᵢ·x)
    int gamma_count = 0;
    for (const auto& f : factors) {
        if (!has(f, var)) {  // a var-free constant factor
            C = mul(C, f);
            continue;
        }
        Expr base = f, ex = S::One();
        if (f->type_id() == TypeId::Pow) {
            base = f->args()[0];
            ex = f->args()[1];
        }
        if (base == var && is_number(ex)) {  // x^q
            q_total = add(q_total, ex);
            continue;
        }
        if (base->type_id() == TypeId::Function && is_number(ex)) {
            const auto& fn = static_cast<const Function&>(*base);
            if (fn.function_id() == FunctionId::Gamma && fn.args().size() == 1) {
                Expr a = simplify(add(fn.args()[0], mul(S::NegativeOne(), var)));
                if (!has(a, var) && is_number(a)) {  // Γ(x + a)^p, a numeric
                    sum_p = add(sum_p, ex);
                    sum_pa = add(sum_pa, mul(ex, a));
                    ++gamma_count;
                    continue;
                }
            }
        }
        // Constant-base exponential c^(k·x + d) with c > 0: contributes k·log c to
        // the growth rate ρ, and the var-free part c^d folds into C. An exponential
        // (ρ ≠ 0) dominates any power of x, so it decides the limit below.
        if (is_number(base) && is_positive(base) == std::optional<bool>{true}) {
            Expr k = simplify(diff(ex, var));
            if (is_number(k)) {
                Expr d = simplify(add(ex, mul(S::NegativeOne(), mul(k, var))));
                if (!has(d, var)) {
                    log_rate = add(log_rate, mul(k, log(base)));
                    C = mul(C, pow(base, d));
                    continue;
                }
            }
        }
        return std::nullopt;  // unrecognized var-dependent factor
    }
    if (gamma_count == 0) return std::nullopt;
    Expr Cs = simplify(C);
    // Sign of a constant expression: assumption query first, then a numeric probe.
    auto sign_of = [](const Expr& e) -> int {
        if (is_positive(e) == std::optional<bool>{true}) return 1;
        if (is_negative(e) == std::optional<bool>{true}) return -1;
        Expr fv = evalf(e, 30);
        if (fv->type_id() == TypeId::Float) {
            try {
                const double d = std::stod(fv->str());
                if (d > 1e-12) return 1;
                if (d < -1e-12) return -1;
            } catch (...) {
            }
        }
        return 0;
    };
    // Unbalanced gamma power: with every gamma a same-rate shift Γ(x+aᵢ), a net
    // exponent Σpᵢ ≠ 0 makes the gamma part ~ (x^x)^{Σpᵢ} by Stirling, which grows
    // (Σpᵢ > 0) or decays (Σpᵢ < 0) faster than any c^{ρx}·x^E. So Σpᵢ < 0 → 0 and
    // Σpᵢ > 0 → sign(C)·∞. Resolves Γ(2x)/Γ(x) → ∞ (which duplication turns into
    // 4ˣ·½·π^(−½)·Γ(x+½)) and avoids the slow gamma-growth fallback.
    Expr sp = simplify(sum_p);
    if (!(sp == S::Zero())) {
        const int sps = sign_of(sp);
        if (sps < 0) return S::Zero();
        if (sps > 0) {
            const int cs = sign_of(Cs);
            if (cs > 0) return S::Infinity();
            if (cs < 0) return S::NegativeInfinity();
        }
        return std::nullopt;
    }
    Expr rho = simplify(log_rate);
    if (!(rho == S::Zero())) {
        // e^{ρx}·x^E: the exponential dominates. ρ < 0 → 0; ρ > 0 → sign(C)·∞.
        const int rs = sign_of(rho);
        if (rs < 0) return S::Zero();
        if (rs > 0) {
            const int cs = sign_of(Cs);
            if (cs > 0) return S::Infinity();
            if (cs < 0) return S::NegativeInfinity();
        }
        return std::nullopt;
    }
    Expr E = simplify(add(q_total, sum_pa));
    if (E == S::Zero()) return Cs;  // limit is the leading constant
    if (is_negative(E) == std::optional<bool>{true}) return S::Zero();
    if (is_positive(E) == std::optional<bool>{true}) {
        if (is_positive(Cs) == std::optional<bool>{true}) return S::Infinity();
        if (is_negative(Cs) == std::optional<bool>{true})
            return S::NegativeInfinity();
    }
    return std::nullopt;
}

// Gauss multiplication for a Γ whose argument grows at an integer multiple of the
// rate, Γ(kx+b) with k ≥ 2, which gamma_ratio_asymptotic (slope-1 only) cannot
// rank — e.g. the central binomial Γ(2x+1)/Γ(x+1)²/4ˣ → 0, Γ(2x)/Γ(x)² → ∞,
// Γ(3x)/Γ(x)³ → ∞. The identity (Legendre duplication generalized)
//   Γ(kz) = (2π)^((1−k)/2)·k^(kz−1/2)·∏_{j=0}^{k−1} Γ(z + j/k),   z = x + b/k
// rewrites the slope-k gamma into k slope-1 gammas plus a constant-base
// exponential k^(kx+…), whose rate the exponential-rate branch of
// gamma_ratio_asymptotic absorbs (cancelling an explicit kᵏ^(−x) via log sums).
// Resolved through that asymptotic directly — never a recursive limit — so it
// succeeds exactly when the rewrite is gamma-ratio-shaped and abstains otherwise.
[[nodiscard]] std::optional<Expr> try_gamma_multiplication(const Expr& expr,
                                                           const Expr& var,
                                                           const Expr& target,
                                                           int depth) {
    if (depth >= 10 || target->type_id() != TypeId::Infinity) {
        return std::nullopt;
    }
    if (count_gamma_factorial(expr) == 0) return std::nullopt;
    Expr work = normalize_factorial(expr);  // factorial → gamma
    ExprMap<Expr> m;
    auto scan = [&](auto&& self, const Expr& e) -> void {
        if (e->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*e);
            if (fn.function_id() == FunctionId::Gamma && fn.args().size() == 1
                && !m.count(e)) {
                const Expr& g = fn.args()[0];
                if (has(g, var)) {
                    Expr slope = simplify(diff(g, var));
                    // Integer rate k in [2, 6] (the product has k gammas).
                    if (slope->type_id() == TypeId::Integer
                        && static_cast<const Integer&>(*slope).fits_long()) {
                        const long k = static_cast<const Integer&>(*slope).to_long();
                        Expr b = simplify(add(
                            g, mul(integer(-k), var)));  // intercept g − k·x
                        if (k >= 2 && k <= 6 && !has(b, var)) {
                            // ∏_{j=0}^{k−1} Γ(x + (b+j)/k)
                            std::vector<Expr> factors;
                            factors.reserve(static_cast<std::size_t>(k) + 2);
                            factors.push_back(
                                pow(mul(integer(2), S::Pi()),
                                    rational(1 - k, 2)));  // (2π)^((1−k)/2)
                            factors.push_back(
                                pow(integer(k),
                                    add(mul(integer(k), var),
                                        add(b, rational(-1, 2)))));  // k^(kx+b−½)
                            for (long j = 0; j < k; ++j) {
                                factors.push_back(gamma(add(
                                    var, mul(add(b, integer(j)),
                                             pow(integer(k), S::NegativeOne())))));
                            }
                            m.emplace(e, mul(std::move(factors)));
                        }
                    }
                }
            }
        }
        for (const auto& a : e->args()) self(self, a);
    };
    scan(scan, work);
    if (m.empty()) return std::nullopt;
    // Resolve the rewritten form *only* through the gamma-ratio asymptotic, not a
    // full recursive limit. That keeps the rule sound and bounded: it succeeds
    // exactly when the rewrite produced a gamma-ratio-shaped form (gammas, x^q,
    // and constant-base exponentials whose rates it now absorbs) and abstains
    // otherwise — e.g. Γ(2n)·n^(−n), whose surviving n^(−n) is not gamma-ratio
    // shaped, is left to the super-power / power-as-exp paths rather than looping.
    return gamma_ratio_asymptotic(xreplace(work, m), var);
}

// Bounded × vanishing at ±∞: a product carrying a bounded oscillating factor
// (sin/cos of a real argument, value in [−1, 1]) times a remaining product that
// → 0 has limit 0 by the squeeze theorem (|sin(g)·rest| ≤ |rest| → 0). Resolves
// x·cos(x)·e^(−x) → 0, x²·sin(x)·e^(−x) → 0 and the like, where the oscillating
// factor otherwise leaves the growth/substitution machinery at nan. Returns
// nullopt unless an oscillating factor is present and the rest provably vanishes,
// so non-vanishing or non-oscillating products are untouched.
[[nodiscard]] std::optional<Expr> bounded_times_vanishing(const Expr& expr,
                                                          const Expr& var,
                                                          const Expr& target,
                                                          int depth) {
    if (!is_infinity(target) || expr->type_id() != TypeId::Mul) {
        return std::nullopt;
    }
    std::vector<Expr> rest;
    bool has_bounded = false;
    for (const auto& f : expr->args()) {
        bool bounded = false;
        if (f->type_id() == TypeId::Function) {
            const auto id = static_cast<const Function&>(*f).function_id();
            if ((id == FunctionId::Sin || id == FunctionId::Cos)
                && is_real(f->args()[0]) != std::optional<bool>{false}) {
                bounded = true;  // sin/cos of a real argument ∈ [−1, 1]
            }
        }
        if (bounded) {
            has_bounded = true;
        } else {
            rest.push_back(f);
        }
    }
    if (!has_bounded) return std::nullopt;
    Expr rest_prod = rest.empty() ? Expr{S::One()} : mul(std::move(rest));
    if (limit_impl(rest_prod, var, target, depth + 1) == S::Zero()) {
        return S::Zero();  // bounded · 0 = 0
    }
    return std::nullopt;
}

// A ratio of sums mixing polynomial and bounded-oscillating terms at ±∞ —
// (x + sin x)/(x + cos x) → 1. Direct substitution leaks sin(∞)/cos(∞) garbage,
// and L'Hôpital oscillates (the derivative of the bounded part has no limit).
// Splitting numerator and denominator into a polynomial skeleton plus a bounded
// O(1) remainder, the bounded remainder is negligible whenever the polynomial
// skeleton → ±∞, so the limit is the rational ratio of those skeletons. Fires
// only when an oscillating sin/cos of the variable is actually present.
[[nodiscard]] std::optional<Expr> try_oscillating_rational(const Expr& expr,
                                                           const Expr& var,
                                                           const Expr& target,
                                                           int depth) {
    if (!is_infinity(target) || expr->type_id() != TypeId::Mul) {
        return std::nullopt;
    }
    auto is_bounded_fn = [&](const Expr& f) {
        if (f->type_id() != TypeId::Function) return false;
        const auto id = static_cast<const Function&>(*f).function_id();
        return (id == FunctionId::Sin || id == FunctionId::Cos)
               && f->args().size() == 1
               && is_real(f->args()[0]) != std::optional<bool>{false};
    };
    // A term that stays bounded as var → ∞: every var-dependent factor is a
    // bounded sin/cos, with at least one such factor.
    auto is_bounded_term = [&](const Expr& t) {
        std::vector<Expr> fs;
        if (t->type_id() == TypeId::Mul) {
            for (const auto& f : t->args()) fs.push_back(f);
        } else {
            fs.push_back(t);
        }
        bool any = false;
        for (const auto& f : fs) {
            if (!has(f, var)) continue;
            if (is_bounded_fn(f)) { any = true; continue; }
            return false;  // a growing var-dependent factor (e.g. x·sin x)
        }
        return any;
    };
    // Polynomial part of a sum, dropping bounded terms; nullopt if a term is
    // neither a clean polynomial in var nor bounded.
    auto poly_part = [&](const Expr& e) -> std::optional<Expr> {
        std::vector<Expr> terms;
        if (e->type_id() == TypeId::Add) {
            for (const auto& t : e->args()) terms.push_back(t);
        } else {
            terms.push_back(e);
        }
        std::vector<Expr> poly;
        for (const auto& t : terms) {
            bool is_poly = false;
            try {
                Poly probe(t, var);
                is_poly = true;
                for (const auto& cc : probe.coeffs()) {
                    if (has(cc, var)) { is_poly = false; break; }
                }
            } catch (const std::exception&) {
                is_poly = false;
            }
            if (is_poly) {
                poly.push_back(t);
            } else if (!is_bounded_term(t)) {
                return std::nullopt;
            }
        }
        return poly.empty() ? Expr{S::Zero()} : add(std::move(poly));
    };
    // Split expr into numerator / denominator (negative integer powers → denom).
    Expr num = S::One(), den = S::One();
    for (const auto& f : expr->args()) {
        if (f->type_id() == TypeId::Pow
            && f->args()[1]->type_id() == TypeId::Integer
            && static_cast<const Integer&>(*f->args()[1]).is_negative()
            && static_cast<const Integer&>(*f->args()[1]).fits_long()) {
            den = mul(den, pow(f->args()[0],
                               integer(-static_cast<const Integer&>(*f->args()[1])
                                            .to_long())));
        } else {
            num = mul(num, f);
        }
    }
    // Only engage when a genuine oscillation of the variable is present.
    bool has_osc = false;
    auto scan = [&](auto&& self, const Expr& e) -> void {
        if (has_osc) return;
        if (is_bounded_fn(e) && has(e, var)) { has_osc = true; return; }
        for (const auto& a : e->args()) self(self, a);
    };
    scan(scan, expr);
    if (!has_osc) return std::nullopt;

    auto np = poly_part(num);
    auto dp = poly_part(den);
    if (!np || !dp) return std::nullopt;
    Expr lnp = limit_impl(*np, var, target, depth + 1);
    Expr ldp = limit_impl(*dp, var, target, depth + 1);
    if (!is_infinity(lnp) || !is_infinity(ldp)) return std::nullopt;
    return limit_impl(mul(*np, pow(*dp, S::NegativeOne())), var, target,
                      depth + 1);
}

// A super-power var^(c·var) (c a nonzero rational) dominates the factorial/gamma
// class: n^n ≫ n!. When the expression is exactly such a super-power times a
// single factorial(var) / gamma(var+1) raised to ±1 (plus constants), the
// super-power's sign decides the limit at +∞ — ∞ when it sits in the numerator,
// 0 in the denominator. Restricted to a lone factorial of the matching variable,
// so gamma(2n) (which outgrows n^n) is left to other methods rather than given a
// wrong answer. Closes n!/n^n → 0, n^n/n! → ∞.
[[nodiscard]] std::optional<Expr> superpow_vs_factorial(const Expr& expr,
                                                        const Expr& var) {
    Expr t = together(expr);
    std::vector<Expr> factors;
    if (t->type_id() == TypeId::Mul) {
        for (const auto& f : t->args()) factors.push_back(f);
    } else {
        factors.push_back(t);
    }
    int s_pow = 0;  // sign of the super-power's exponent coefficient
    bool has_factorial = false;
    for (const auto& f : factors) {
        if (!has(f, var)) continue;  // constant factor — dominated
        // Super-power var^(c·var).
        if (f->type_id() == TypeId::Pow && f->args()[0] == var
            && has(f->args()[1], var)) {
            Expr coef = simplify(f->args()[1] / var);  // c when exp = c·var
            if (has(coef, var) || !is_number(coef)) return std::nullopt;
            const int sc = is_positive(coef) == std::optional<bool>{true}  ? 1
                           : is_negative(coef) == std::optional<bool>{true} ? -1
                                                                            : 0;
            if (sc == 0 || s_pow != 0) return std::nullopt;
            s_pow = sc;
            continue;
        }
        // factorial(var) or gamma(var+1), optionally inverted.
        Expr base = f;
        if (f->type_id() == TypeId::Pow && f->args()[1] == S::NegativeOne()) {
            base = f->args()[0];
        }
        if (base->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*base);
            const auto id = fn.function_id();
            if (id == FunctionId::Factorial && fn.args().size() == 1
                && fn.args()[0] == var) {
                has_factorial = true;
                continue;
            }
            if (id == FunctionId::Gamma && fn.args().size() == 1
                && simplify(fn.args()[0] - integer(1)) == var) {
                has_factorial = true;
                continue;
            }
        }
        return std::nullopt;  // some other var-dependent factor — bail
    }
    if (s_pow == 0 || !has_factorial) return std::nullopt;
    return s_pow > 0 ? Expr{S::Infinity()} : Expr{S::Zero()};
}

// Stirling-asymptotic limits of n-th roots of factorial/gamma: (n!)^(1/n)/n →
// 1/e, (n!)^(1/n) → ∞, n/(n!)^(1/n) → e, (n!/nⁿ)^(1/n) → 1/e. These need the
// growth rate of n! that the elementary growth hierarchy lacks. Substitute the
// variable with a *positive* symbol (valid at +∞, and it lets the powdenest
// rules collapse ((m/e)^m)^(1/m) → m/e) and replace each factorial/gamma by its
// leading Stirling form g! ~ (g/e)^g. The candidate is accepted only after a
// numeric check against the original at large n, so the dropped subleading
// Stirling terms cannot produce a wrong answer.
[[nodiscard]] std::optional<Expr> try_stirling_limit(const Expr& expr,
                                                     const Expr& var,
                                                     const Expr& target,
                                                     int depth) {
    if (target != S::Infinity() || count_gamma_factorial(expr) == 0) {
        return std::nullopt;
    }
    AssumptionMask pm;
    pm.set_positive(true);
    Expr m = symbol("__stirling_m", pm);
    Expr em = subs(expr, var, m);  // recast over a positive variable
    // Build the factorial/gamma → Stirling rewrite map. `full` selects the form:
    //   leading  g! ~ (g/e)^g — its subleading prefactor cancels under an n-th root
    //            (the canonical (n!)^(1/n)/n → 1/e case);
    //   full     g! = Γ(g+1) ~ √(2π)·g^(g+1/2)·e⁻ᵍ, Γ(g) ~ √(2π)·g^(g−1/2)·e⁻ᵍ —
    //            the prefactor is the *entire* limit once (g/e)^g cancels, e.g.
    //            n!·n⁻ⁿ·eⁿ ~ √(2πn) → ∞.
    auto build_map = [&](bool full) -> std::optional<ExprMap<Expr>> {
        ExprMap<Expr> map;
        bool any = false;
        auto collect = [&](auto&& self, const Expr& e) -> void {
            if (e->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*e);
                const auto id = fn.function_id();
                if ((id == FunctionId::Factorial || id == FunctionId::Gamma)
                    && fn.args().size() == 1 && has(fn.args()[0], m)
                    && map.find(e) == map.end()) {
                    const Expr& g = fn.args()[0];
                    if (limit_impl(g, m, S::Infinity(), depth + 1)
                        == S::Infinity()) {
                        Expr rep = pow(mul(g, exp(S::NegativeOne())), g);
                        if (full) {
                            const Expr half = id == FunctionId::Factorial
                                                  ? Expr{rational(1, 2)}
                                                  : Expr{rational(-1, 2)};
                            rep = mul(mul(pow(mul(integer(2), S::Pi()),
                                              rational(1, 2)),
                                          pow(g, half)),
                                      rep);
                        }
                        map.emplace(e, std::move(rep));
                        any = true;
                    }
                }
            }
            for (const auto& a : e->args()) self(self, a);
        };
        collect(collect, em);
        if (!any) return std::nullopt;
        return map;
    };
    // Evaluate a rewrite map to a candidate and validate it numerically: sample the
    // ORIGINAL expression at increasing n and require the values to track the
    // candidate (converge to a finite value, or grow/decay consistently with ±∞).
    // A wrong candidate keeps an O(1) offset.
    auto evaluate = [&](const ExprMap<Expr>& map) -> std::optional<Expr> {
        Expr cand = limit_impl(simplify(xreplace(em, map)), m, S::Infinity(),
                               depth + 1);
        if (is_nan(cand) || has(cand, m)
            || cand->type_id() == TypeId::ComplexInfinity) {
            return std::nullopt;
        }
        const bool cand_pinf = cand->type_id() == TypeId::Infinity;
        const bool cand_ninf = cand->type_id() == TypeId::NegativeInfinity;
        double cv = 0.0;
        if (!cand_pinf && !cand_ninf) {
            Expr cvf = evalf(cand, 30);
            if (cvf->type_id() != TypeId::Float) return std::nullopt;
            try {
                cv = std::stod(cvf->str());
            } catch (...) {
                return std::nullopt;
            }
        }
        double prev_diff = 1e300, prev_val = 0.0;
        int checks = 0;
        for (long nv : {300L, 1000L, 3000L}) {
            Expr at = evalf(simplify(subs(expr, var, integer(nv))), 40);
            if (at->type_id() != TypeId::Float) return std::nullopt;
            double v = 0.0;
            try {
                v = std::stod(at->str());
            } catch (...) {
                return std::nullopt;
            }
            if (cand_pinf) {
                // The Stirling rewrite is asymptotically exact, so a divergent
                // candidate only needs the samples to climb (a slow sqrt(n)/e or
                // log(n) limit never reaches a fixed threshold).
                if (checks > 0 && !(v > prev_val)) return std::nullopt;
                if (nv == 3000L && !(v > 1.0)) return std::nullopt;
            } else if (cand_ninf) {
                if (checks > 0 && !(v < prev_val)) return std::nullopt;
                if (nv == 3000L && !(v < -1.0)) return std::nullopt;
            } else {
                const double d = std::fabs(v - cv);
                if (!(d <= prev_diff + 1e-12)) return std::nullopt;  // not diverge
                prev_diff = d;
                if (nv == 3000L && !(d < 1e-2)) return std::nullopt;  // slow ~log n/n
            }
            prev_val = v;
            ++checks;
        }
        return checks == 3 ? std::optional<Expr>{cand} : std::nullopt;
    };
    // Leading form first (keeps the n-th-root path exact); the full form rescues
    // products where (g/e)^g cancels and the prefactor carries the limit.
    for (bool full : {false, true}) {
        auto map = build_map(full);
        if (!map) return std::nullopt;  // no divergent factorial/gamma present
        if (auto r = evaluate(*map)) return *r;
    }
    return std::nullopt;
}

// True if e contains a non-integer power of a var-dependent base (a radical such
// as √(x²+1) = (x²+1)^(1/2)), which the conjugate rationalization targets.
[[nodiscard]] bool has_var_radical(const Expr& e, const Expr& var) {
    if (e->type_id() == TypeId::Pow
        && e->args()[1]->type_id() != TypeId::Integer
        && has(e->args()[0], var)) {
        return true;
    }
    for (const auto& a : e->args()) {
        if (has_var_radical(a, var)) return true;
    }
    return false;
}

// LCM of the denominators of the fractional exponents over var-dependent bases —
// the root order to raise a difference to so every radical clears (√ → 2, ∛ → 3,
// a mix of √ and ∛ → 6). Returns 1 when there is no var-radical.
[[nodiscard]] long radical_order(const Expr& e, const Expr& var) {
    long n = 1;
    auto rec = [&](auto&& self, const Expr& x) -> void {
        if (x->type_id() == TypeId::Pow
            && x->args()[1]->type_id() == TypeId::Rational
            && has(x->args()[0], var)) {
            const mpz_class den =
                static_cast<const Rational&>(*x->args()[1]).value().get_den();
            if (den.fits_slong_p()) n = std::lcm(n, den.get_si());
        }
        for (const auto& a : x->args()) self(self, a);
    };
    rec(rec, e);
    return n;
}

// ∞ − ∞ involving a radical, resolved by the conjugate. For a two-term sum
// t₁ + t₂ (one term a square root) whose direct limit is the indeterminate
// ∞ − ∞: t₁ + t₂ = (t₁² − t₂²)/(t₁ − t₂). Squaring clears the radical and the
// ratio resolves — e.g. x − √(x²+1) → 0, √(x²+x) − x → 1/2, √(x+1) − √x → 0.
[[nodiscard]] std::optional<Expr> try_conjugate_difference(
    const Expr& expr, const Expr& var, const Expr& target, int depth) {
    if (expr->type_id() != TypeId::Add || expr->args().size() != 2) {
        return std::nullopt;
    }
    if (!has_var_radical(expr, var)) return std::nullopt;
    const Expr& t1 = expr->args()[0];
    const Expr& t2 = expr->args()[1];
    Expr num = simplify(mul(t1, t1) + mul(S::NegativeOne(), mul(t2, t2)));
    Expr den = simplify(t1 + mul(S::NegativeOne(), t2));  // t₁ − t₂
    if (den == S::Zero()) return std::nullopt;
    // The conjugate t₁²−t₂² normally clears the radical; when it does not (a
    // nested radical like √(x+√x)−√x leaves num = √x), the rationalised ratio
    // num/den is still a determinate lower-order form (√x/(√(x+√x)+√x) → 1/2), so
    // attempt it rather than abstaining — but only when the residual num is a
    // *simple* var radical c·x^q. That is exactly the nested-difference shape and
    // keeps the engine from recursing on every other radical ∞−∞. Failure falls
    // through on the nan check below, so it cannot loop.
    if (has_var_radical(num, var)) {
        Expr core = num;
        if (num->type_id() == TypeId::Mul) {
            int var_factors = 0;
            for (const auto& f : num->args()) {
                if (has(f, var)) {
                    ++var_factors;
                    core = f;
                }
            }
            if (var_factors != 1) core = S::NaN();  // not a single var factor
        }
        const bool simple_radical =
            core->type_id() == TypeId::Pow && core->args()[0] == var
            && core->args()[1]->type_id() == TypeId::Rational;
        if (!simple_radical) return std::nullopt;
    }
    // Pass the ratio UNSIMPLIFIED: simplify() would rationalize the denominator
    // straight back to the original ∞ − ∞ form and loop. limit_impl substitutes
    // before simplifying, so the genuine pole collapses first.
    Expr val = limit_impl(mul(num, pow(den, S::NegativeOne())), var, target,
                          depth + 1);
    if (is_nan(val)) return std::nullopt;
    return val;
}

// Leading asymptotic term of an algebraic expression as var → +∞: returns
// (coeff, degree) with e ~ coeff · var^degree (coeff ≠ 0), or nullopt when the
// behavior is not a single algebraic power (a non-algebraic function, a symbolic
// exponent, or a leading cancellation the conjugate can't clear). Degrees may be
// rational (√ halves them). This is the leading-order slice of the Gruntz/MRV
// algorithm restricted to polynomials and their roots — enough to resolve the
// √-difference limits (x − √(x²−x), √(x²+x) − √(x²−x), x·(√(x²+1) − x)) that
// L'Hôpital can't (the radical never stabilises under repeated differentiation).
[[nodiscard]] std::optional<std::pair<Expr, Expr>> leading_pos_inf(
    const Expr& e, const Expr& var, int depth) {
    // Bound the recursion: the conjugate-on-cancellation branch can re-enter a few
    // levels deep on a single expression, but never unboundedly (each conjugate
    // clears a leading cancellation). 24 covers realistic nesting and still
    // terminates on a pathological input.
    if (depth > 24) return std::nullopt;
    if (!has(e, var)) {
        if (e == S::Zero()) return std::nullopt;  // zero has no leading term
        return std::make_pair(e, Expr{S::Zero()});
    }
    if (e == var) return std::make_pair(Expr{S::One()}, Expr{S::One()});
    switch (e->type_id()) {
        case TypeId::Pow: {
            const Expr& b = e->args()[0];
            const Expr& p = e->args()[1];
            if (has(p, var)) return std::nullopt;  // var in exponent: not algebraic
            auto lb = leading_pos_inf(b, var, depth + 1);
            if (!lb) return std::nullopt;
            // A fractional power needs a positive leading coefficient to stay real.
            if (p->type_id() != TypeId::Integer
                && is_positive(lb->first) != std::optional<bool>{true}) {
                return std::nullopt;
            }
            return std::make_pair(simplify(pow(lb->first, p)),
                                  simplify(mul(lb->second, p)));
        }
        case TypeId::Mul: {
            Expr c = S::One();
            Expr d = S::Zero();
            for (const auto& f : e->args()) {
                auto lf = leading_pos_inf(f, var, depth + 1);
                if (!lf) return std::nullopt;
                c = mul(c, lf->first);
                d = add(d, lf->second);
            }
            return std::make_pair(simplify(c), simplify(d));
        }
        case TypeId::Add: {
            std::optional<Expr> maxdeg;
            Expr csum = S::Zero();
            for (const auto& t : e->args()) {
                auto lt = leading_pos_inf(t, var, depth + 1);
                if (!lt) return std::nullopt;
                if (!maxdeg) {
                    maxdeg = lt->second;
                    csum = lt->first;
                    continue;
                }
                Expr dd = simplify(lt->second - *maxdeg);
                if (is_positive(dd) == std::optional<bool>{true}) {
                    maxdeg = lt->second;
                    csum = lt->first;
                } else if (dd == S::Zero()) {
                    csum = simplify(csum + lt->first);
                }
            }
            if (!maxdeg) return std::nullopt;
            if (simplify(csum) != S::Zero()) {
                return std::make_pair(simplify(csum), *maxdeg);
            }
            // Leading terms cancel. For a two-term sum u − v with a radical of order
            // n, the n-th-root conjugate clears it:
            //   u − v = (uⁿ − vⁿ) / Σ_{i=0}^{n-1} u^(n−1−i)·vⁱ.
            // uⁿ, vⁿ raise the radicals to integer powers, so the numerator becomes
            // radical-free; the denominator has no leading cancellation. Recurse on
            // the quotient (no simplify — that would rationalize straight back).
            // n = 2 is the classic √ conjugate; n = 3 handles (x³+x²)^(1/3) − x → 1/3.
            if (e->args().size() == 2 && has_var_radical(e, var)) {
                const Expr u = e->args()[0];
                const Expr v = mul(S::NegativeOne(), e->args()[1]);  // e = u − v
                const long n = radical_order(e, var);
                if (n >= 2 && n <= 6) {
                    Expr num = simplify(pow(u, integer(n))
                                        + mul(S::NegativeOne(), pow(v, integer(n))));
                    if (!has_var_radical(num, var)) {
                        std::vector<Expr> dterms;
                        dterms.reserve(static_cast<std::size_t>(n));
                        for (long i = 0; i < n; ++i) {
                            dterms.push_back(
                                mul(pow(u, integer(n - 1 - i)), pow(v, integer(i))));
                        }
                        Expr den = add(std::move(dterms));
                        if (!(simplify(den) == S::Zero())) {
                            return leading_pos_inf(
                                mul(num, pow(den, S::NegativeOne())), var,
                                depth + 1);
                        }
                    }
                }
            }
            return std::nullopt;
        }
        default:
            return std::nullopt;  // trig/exp/log/… are not handled here
    }
}

// Limit at +∞ of an algebraic (radical) expression via its leading term. Resolves
// the ∞ − ∞ / 0·∞ forms that L'Hôpital abandons on radicals: e ~ c·x^d gives a
// finite limit c when d = 0, ±∞ when d > 0, and 0 when d < 0.
[[nodiscard]] std::optional<Expr> try_algebraic_inf(const Expr& expr,
                                                    const Expr& var,
                                                    const Expr& target) {
    if (target->type_id() != TypeId::Infinity) return std::nullopt;  // +∞ only
    if (!has_var_radical(expr, var)) return std::nullopt;
    auto lead = leading_pos_inf(expr, var, 0);
    if (!lead) return std::nullopt;
    const Expr& c = lead->first;
    const Expr& d = lead->second;
    if (d == S::Zero()) return simplify(c);
    if (is_positive(d) == std::optional<bool>{true}) {
        if (is_positive(c) == std::optional<bool>{true}) return S::Infinity();
        if (is_negative(c) == std::optional<bool>{true}) return S::NegativeInfinity();
        return std::nullopt;
    }
    if (is_negative(d) == std::optional<bool>{true}) return S::Zero();
    return std::nullopt;
}

// Asymptotic of log(g) at +∞ when g is dominated by an exponential — the
// "log-sum-exp → max" identity. For g = Σ cᵢ·exp(eᵢ) (after rewriting cosh/sinh
// and a^x into exp), factoring out the fastest-growing exp(e_dom) gives
//   log(g) = e_dom + log(g·exp(−e_dom)),   g·exp(−e_dom) → (finite) > 0,
// so the residual log has a finite limit. Rewrites the first such log(g) inside
// `expr` and re-takes the limit. Resolves x − log(cosh x) → log 2,
// log(2^x+3^x)/x → log 3, (2^x+3^x)^(1/x) → 3.
[[nodiscard]] std::optional<Expr> try_log_exp_asymptotic(const Expr& expr,
                                                         const Expr& var,
                                                         const Expr& target,
                                                         int depth) {
    if (target->type_id() != TypeId::Infinity) return std::nullopt;  // +∞ only

    auto is_log = [](const Expr& f) {
        return f->type_id() == TypeId::Function
               && static_cast<const Function&>(*f).function_id()
                      == FunctionId::Log
               && f->args().size() == 1;
    };
    // Find the first log(g) with g depending on var.
    Expr the_log;
    auto find = [&](auto&& self, const Expr& e) -> void {
        if (the_log) return;
        if (is_log(e) && has(e->args()[0], var)) { the_log = e; return; }
        for (const auto& a : e->args()) {
            if (the_log) return;
            self(self, a);
        }
    };
    find(find, expr);
    if (!the_log) return std::nullopt;

    // Rewrite cosh(u)/sinh(u) → (e^u ± e^(−u))/2 inside g, then expand to a sum.
    ExprMap<Expr> hyp;
    auto scan = [&](auto&& self, const Expr& e) -> void {
        if (e->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*e);
            if (fn.args().size() == 1) {
                const Expr& u = fn.args()[0];
                if (fn.function_id() == FunctionId::Cosh) {
                    hyp.emplace(e, (exp(u) + exp(mul(S::NegativeOne(), u)))
                                       / integer(2));
                } else if (fn.function_id() == FunctionId::Sinh) {
                    hyp.emplace(e, (exp(u) - exp(mul(S::NegativeOne(), u)))
                                       / integer(2));
                }
            }
        }
        for (const auto& a : e->args()) self(self, a);
    };
    const Expr& g = the_log->args()[0];
    scan(scan, g);
    Expr g_exp = expand(hyp.empty() ? g : xreplace(g, hyp));

    // The exponent e such that `term = const·exp(e)` (a^h ↦ e = h·log a; a bare
    // constant ↦ e = 0). Bails on a non-exponential var factor (e.g. x·eˣ).
    auto exp_exponent = [&](const Expr& term) -> std::optional<Expr> {
        if (!has(term, var)) return Expr{S::Zero()};
        std::vector<Expr> facs;
        if (term->type_id() == TypeId::Mul) {
            for (const auto& f : term->args()) facs.push_back(f);
        } else {
            facs.push_back(term);
        }
        Expr exponent;
        bool found = false;
        for (const auto& f : facs) {
            if (!has(f, var)) continue;
            if (f->type_id() == TypeId::Function
                && static_cast<const Function&>(*f).function_id()
                       == FunctionId::Exp) {
                if (found) return std::nullopt;
                exponent = f->args()[0];
                found = true;
            } else if (f->type_id() == TypeId::Pow && !has(f->args()[0], var)
                       && is_positive(f->args()[0]) == std::optional<bool>{true}
                       && has(f->args()[1], var)) {
                if (found) return std::nullopt;
                exponent = mul(f->args()[1], log(f->args()[0]));
                found = true;
            } else {
                return std::nullopt;
            }
        }
        return found ? exponent : Expr{S::Zero()};
    };

    // Pick the term with the fastest-growing exponent (max coefficient of x).
    std::vector<Expr> terms;
    if (g_exp->type_id() == TypeId::Add) {
        for (const auto& t : g_exp->args()) terms.push_back(t);
    } else {
        terms.push_back(g_exp);
    }
    Expr e_dom;
    double best_rate = 0.0;
    for (const auto& t : terms) {
        auto e = exp_exponent(t);
        if (!e) return std::nullopt;
        // Rate = coefficient of x in the exponent (must be affine in var).
        Expr rate_coeff;
        try {
            Poly pe(expand(*e), var);
            if (pe.degree() < 1) {
                rate_coeff = S::Zero();  // constant exponent
            } else if (pe.degree() == 1) {
                rate_coeff = pe.coeffs()[1];
            } else {
                return std::nullopt;  // super-exponential — out of scope
            }
        } catch (const std::exception&) {
            return std::nullopt;
        }
        Expr rate_e = evalf(rate_coeff, 30);
        double rate = 0.0;
        try {
            rate = std::stod(rate_e->str());
        } catch (...) {
            return std::nullopt;  // non-numeric rate (e.g. a symbolic coefficient)
        }
        if (!e_dom || rate > best_rate) {
            best_rate = rate;
            e_dom = *e;
        }
    }
    if (!e_dom || best_rate <= 0.0) return std::nullopt;  // not exp-growing

    // log(g) → e_dom + log(g·exp(−e_dom)); replace and re-limit. Use the
    // exp-rewritten g_exp (cosh/sinh already in exponential form) so the residual
    // collapses — e.g. cosh(x)·e^(−x) = (1 + e^(−2x))/2 → 1/2 rather than the
    // unevaluated cosh(x)·e^(−x) the limit engine can't fold directly.
    Expr residual = simplify(expand(mul(g_exp, exp(mul(S::NegativeOne(), e_dom)))));
    Expr new_log = add(e_dom, log(residual));
    ExprMap<Expr> repl;
    repl.emplace(the_log, new_log);
    // expand so a wrapping product distributes — log(g)/x must become
    // e_dom/x + log(residual)/x for the per-term limits to resolve, rather than
    // staying an (∞)·(0) Mul.
    Expr expr2 = expand(xreplace(expr, repl));
    Expr r = limit_impl(expr2, var, target, depth + 1);
    if (is_nan(r)) return std::nullopt;
    return r;
}

// Resolve logarithms at the target.
//   (a) limit(log(g)) = log(limit g)   for a positive-finite or +∞ inner limit;
//   (b) Σ cᵢ·log(gᵢ) + rest: factor a common κ so every cᵢ/κ is an integer,
//       giving κ·log(∏ gᵢ^(cᵢ/κ)) with a single (rational) argument — this
//       resolves the ∞ − ∞ between logs, e.g. log(x+1) − log(x) → 0 and the
//       log terms of the ∫1/(1+x⁴) antiderivative at ∞.
[[nodiscard]] std::optional<Expr> try_log_limit(const Expr& expr,
                                                const Expr& var,
                                                const Expr& target, int depth) {
    auto is_log = [](const Expr& f) {
        return f->type_id() == TypeId::Function
               && static_cast<const Function&>(*f).function_id()
                      == FunctionId::Log
               && f->args().size() == 1;
    };
    if (is_log(expr)) {
        Expr gl = limit_impl(expr->args()[0], var, target, depth + 1);
        if (gl == S::Infinity()) return S::Infinity();
        if (!is_nan(gl) && !is_infinity(gl)
            && is_positive(gl) == std::optional<bool>{true}) {
            return simplify(log(gl));
        }
        return std::nullopt;
    }
    // atan is continuous on the extended reals (atan(±∞) = ±π/2), so
    // limit(atan(g)) = atan(limit g) whenever the inner limit is determinate.
    if (expr->type_id() == TypeId::Function
        && static_cast<const Function&>(*expr).function_id() == FunctionId::Atan
        && expr->args().size() == 1) {
        Expr gl = limit_impl(expr->args()[0], var, target, depth + 1);
        if (!is_nan(gl)) return simplify(atan(gl));
        return std::nullopt;
    }
    // exp is continuous with exp(+∞) = ∞, exp(−∞) = 0, so limit(exp(g)) =
    // exp(limit g) whenever the inner limit is determinate. Resolves cases where
    // the exponent only substitutes to an indeterminate ∞ − ∞ but has a definite
    // limit — e.g. exp(x²−2x) → ∞, exp(2x−x²) → 0, exp(x−√x) → ∞.
    if (expr->type_id() == TypeId::Function
        && static_cast<const Function&>(*expr).function_id() == FunctionId::Exp
        && expr->args().size() == 1) {
        Expr gl = limit_impl(expr->args()[0], var, target, depth + 1);
        if (gl == S::Infinity()) return S::Infinity();
        if (gl == S::NegativeInfinity()) return S::Zero();
        if (!is_nan(gl) && !is_infinity(gl)) return simplify(exp(gl));
        return std::nullopt;
    }
    if (expr->type_id() != TypeId::Add) return std::nullopt;
    auto as_clog =
        [&](const Expr& t) -> std::optional<std::pair<Expr, Expr>> {  // (g, c)
        if (is_log(t)) return std::make_pair(t->args()[0], Expr{S::One()});
        if (t->type_id() == TypeId::Mul) {
            Expr inner;
            std::vector<Expr> coeff;
            for (const auto& f : t->args()) {
                if (!inner && is_log(f)) {
                    inner = f->args()[0];
                    continue;
                }
                if (has(f, var)) return std::nullopt;
                coeff.push_back(f);
            }
            if (inner) return std::make_pair(inner, mul(std::move(coeff)));
        }
        return std::nullopt;
    };
    std::vector<std::pair<Expr, Expr>> logs;
    std::vector<Expr> rest;
    for (const auto& t : expr->args()) {
        if (auto lp = as_clog(t)) logs.push_back(std::move(*lp));
        else rest.push_back(t);
    }
    if (logs.size() < 2) return std::nullopt;
    const Expr& kappa = logs[0].second;
    if (kappa == S::Zero()) return std::nullopt;
    std::vector<Expr> factors;
    for (const auto& [g, c] : logs) {
        Expr ratio = simplify(c * pow(kappa, integer(-1)));
        if (ratio->type_id() != TypeId::Integer
            || !static_cast<const Integer&>(*ratio).fits_long()) {
            return std::nullopt;  // coefficients not commensurate
        }
        factors.push_back(pow(g, ratio));
    }
    Expr combined = mul(kappa, log(mul(std::move(factors))));
    rest.push_back(combined);
    Expr val = limit_impl(add(std::move(rest)), var, target, depth + 1);
    if (is_nan(val)) return std::nullopt;
    return val;
}

// Gruntz preprocessing: rewrite every general power f(x)^g(x) — where both the
// base and the exponent depend on var — as exp(g·log f). This is the canonical
// first step of the most-rapidly-varying algorithm: with all such powers on a
// common exp footing, the exp/gamma growth machinery can compare and resolve
// shapes the bare-power paths abandon, e.g. n^n against Γ(2n). Constant-base
// exponentials (2^x, exp(x)) are left untouched — they have dedicated handlers
// and a var-free base is not a "general" power.
[[nodiscard]] std::optional<Expr> try_power_as_exp(const Expr& expr,
                                                   const Expr& var,
                                                   const Expr& target,
                                                   int depth) {
    if (depth >= 10) return std::nullopt;
    ExprMap<Expr> m;
    auto collect = [&](auto&& self, const Expr& e) -> void {
        if (e->type_id() == TypeId::Pow) {
            const Expr& base = e->args()[0];
            const Expr& ex = e->args()[1];
            if (has(base, var) && has(ex, var) && !m.count(e)) {
                m.emplace(e, exp(mul(ex, log(base))));
            }
        }
        for (const auto& a : e->args()) self(self, a);
    };
    collect(collect, expr);
    if (m.empty()) return std::nullopt;
    Expr r = limit_impl(xreplace(expr, m), var, target, depth + 1);
    if (is_nan(r)) return std::nullopt;
    return r;
}

// A power f(x)^{g(x)} that tends to 1 — its exponent t = g·log f → 0 — expanded by
// exp(t) − 1 = t + t²/2 + t³/6 + …. The dropped tail is o(t), so for the unit-power
// difference forms (f^g − 1)/h → 1 the leading t carries the limit exactly; the few
// extra terms cover an h ∼ t² or t³. Resolves (x^x − 1)/(x·log x) → 1 and
// (x^(1/x) − 1)/(log x/x) → 1, where the base → 0 or ∞ defeats the series stage
// (its exponent has a log singularity). Numerically verified, so a too-short
// truncation cannot pass.
[[nodiscard]] std::optional<Expr> try_unit_power_expansion(const Expr& expr,
                                                           const Expr& var,
                                                           const Expr& target,
                                                           int depth) {
    if (depth >= 8) return std::nullopt;
    ExprMap<Expr> m;
    // Match a unit-tending power under a sum (the f^g − 1 shape), as a raw
    // f(x)^g(x) or — once the power-as-exp rewrite has fired — as exp(t).
    auto scan = [&](auto&& self, const Expr& e, bool under_add) -> void {
        Expr t;
        bool match = false;
        if (under_add && !m.count(e)) {
            if (e->type_id() == TypeId::Pow && has(e->args()[0], var)
                && has(e->args()[1], var)) {
                t = mul(e->args()[1], log(e->args()[0]));  // g·log f
                match = true;
            } else if (e->type_id() == TypeId::Function
                       && static_cast<const Function&>(*e).function_id()
                              == FunctionId::Exp
                       && has(e->args()[0], var)) {
                t = e->args()[0];
                match = true;
            }
        }
        if (match
            && limit_impl(t, var, target, depth + 1) == S::Zero()) {
            Expr poly = S::One(), tk = S::One();
            long fk = 1;
            for (int kk = 1; kk <= 3; ++kk) {
                tk = mul(tk, t);
                fk *= kk;
                poly = add(poly, mul(tk, pow(integer(fk), integer(-1))));
            }
            m.emplace(e, poly);
            return;  // expanded — do not descend further into this node
        }
        const bool c = under_add || e->type_id() == TypeId::Add;
        for (const auto& a : e->args()) self(self, a, c);
    };
    scan(scan, expr, false);
    if (m.empty()) return std::nullopt;
    // Skip when a divergent exponential exp(h → +∞) multiplies the unit difference
    // (e.g. eˣ·(e^{1/x} − 1)): substituting the expansion leaves an undistributed
    // eˣ·(1 + …) − eˣ — an ∞ − ∞ the reciprocal substitution spins on. Those 0·∞
    // exp products are the MRV-rewrite cases, left to other paths (honest nan), not
    // this elementary expansion. The matched unit exp(t → 0) never trips this.
    bool divergent_exp = false;
    auto scan2 = [&](auto&& self, const Expr& e) -> void {
        if (divergent_exp) return;
        if (e->type_id() == TypeId::Function
            && static_cast<const Function&>(*e).function_id() == FunctionId::Exp
            && has(e->args()[0], var)
            && is_infinity(limit_impl(e->args()[0], var, target, depth + 1))) {
            divergent_exp = true;
            return;
        }
        for (const auto& a : e->args()) self(self, a);
    };
    scan2(scan2, expr);
    if (divergent_exp) return std::nullopt;
    Expr cand = limit_impl(xreplace(expr, m), var, target, depth + 1);
    if (is_nan(cand) || has(cand, var)
        || cand->type_id() == TypeId::ComplexInfinity) {
        return std::nullopt;
    }
    // Numeric guard against a too-short truncation. Sample the ORIGINAL approaching
    // the target: large |x| at ∞, x = a + 1/N near a finite a.
    const bool at_inf = is_infinity(target);
    const Expr sgn = target->type_id() == TypeId::NegativeInfinity
                         ? S::NegativeOne()
                         : S::One();
    const bool c_pinf = cand->type_id() == TypeId::Infinity;
    const bool c_ninf = cand->type_id() == TypeId::NegativeInfinity;
    double cvd = 0.0;
    if (!c_pinf && !c_ninf) {
        Expr cvf = evalf(cand, 30);
        if (cvf->type_id() != TypeId::Float) return std::nullopt;
        try {
            cvd = std::stod(cvf->str());
        } catch (...) {
            return std::nullopt;
        }
    }
    double prev_diff = 1e300, prev_val = 0.0;
    int checks = 0;
    for (long gp : {300L, 3000L, 30000L}) {
        const Expr pt = at_inf
                            ? Expr{mul(sgn, integer(gp))}
                            : Expr{add(target, pow(integer(gp), integer(-1)))};
        Expr at = evalf(subs(expr, var, pt), 60);
        if (at->type_id() != TypeId::Float) return std::nullopt;
        double v = 0.0;
        try {
            v = std::stod(at->str());
        } catch (...) {
            return std::nullopt;
        }
        if (c_pinf) {
            if (checks > 0 && !(v > prev_val)) return std::nullopt;
            if (gp == 30000L && !(v > 1.0)) return std::nullopt;
        } else if (c_ninf) {
            if (checks > 0 && !(v < prev_val)) return std::nullopt;
            if (gp == 30000L && !(v < -1.0)) return std::nullopt;
        } else {
            const double d = std::fabs(v - cvd);
            if (!(d <= prev_diff + 1e-9)) return std::nullopt;
            prev_diff = d;
            if (gp == 30000L && !(d < 1e-2)) return std::nullopt;
        }
        prev_val = v;
        ++checks;
    }
    return cand;
}

// Gruntz MRV rewrite for a difference of exponentials. A sum Σ cᵢ·exp(aᵢ) of
// asymptotically-equal exponentials (every aᵢ − aⱼ → finite) is the canonical
// most-rapidly-varying case: e^{x+e⁻ˣ} − e^x → 1 (Gruntz's flagship example).
// Factor out the common exp(b) for a reference exponent b — exact, since
// exp(b)·Σ cᵢ·exp(aᵢ − b) = Σ cᵢ·exp(aᵢ) — so the difference collapses to a unit
// term the existing machinery resolves: e^x·(e^{e⁻ˣ} − 1) → 1. Each candidate b is
// tried; the rewrite is an identity, so any determinate re-take is exact.
[[nodiscard]] std::optional<Expr> try_exp_difference(const Expr& expr,
                                                     const Expr& var,
                                                     const Expr& target,
                                                     int depth) {
    if (depth >= 8 || !is_infinity(target)
        || expr->type_id() != TypeId::Add) {
        return std::nullopt;
    }
    struct Term {
        Expr coeff;
        Expr expo;
    };
    std::vector<Term> terms;
    for (const auto& t : expr->args()) {
        Expr coeff = S::One(), expo;
        bool found = false;
        std::vector<Expr> facs;
        if (t->type_id() == TypeId::Mul) {
            for (const auto& f : t->args()) facs.push_back(f);
        } else {
            facs.push_back(t);
        }
        for (const auto& f : facs) {
            if (!found && f->type_id() == TypeId::Function
                && static_cast<const Function&>(*f).function_id()
                       == FunctionId::Exp
                && has(f->args()[0], var)) {
                expo = f->args()[0];
                found = true;
            } else {
                coeff = mul(coeff, f);
            }
        }
        if (!found) return std::nullopt;  // a non-exponential term — not this shape
        terms.push_back({coeff, expo});
    }
    if (terms.size() < 2) return std::nullopt;
    for (const auto& ref : terms) {
        bool ok = true;
        Expr inner = S::Zero();
        for (const auto& tm : terms) {
            Expr d = simplify(add(tm.expo, mul(S::NegativeOne(), ref.expo)));
            Expr ld = limit_impl(d, var, target, depth + 1);
            if (is_infinity(ld) || is_nan(ld)
                || ld->type_id() == TypeId::ComplexInfinity) {
                ok = false;
                break;
            }
            inner = add(inner, mul(tm.coeff, exp(d)));  // cᵢ·exp(aᵢ − b)
        }
        if (!ok) continue;
        Expr r = limit_impl(mul(exp(ref.expo), inner), var, target, depth + 1);
        if (!is_nan(r) && r->type_id() != TypeId::ComplexInfinity) return r;
    }
    return std::nullopt;
}

// Gruntz leading-term via series. For an indeterminate form built from a variable
// power inside a sum — g(x)·(f(x)^{h(x)} − c), the 1^∞ minus its value — substitute
// x = ±1/u (so x → ±∞ becomes u → 0⁺), rewrite each u-dependent power to exp(·log·)
// so series() can expand it, factor out the explicit u-power pole, expand the
// regular part, recombine, and read the limit off the leading term. Resolves
// x·((1+1/x)^x − e) → −e/2, which the power-as-exp + reciprocal paths spin on.
[[nodiscard]] std::optional<Expr> try_series_limit(const Expr& expr,
                                                   const Expr& var,
                                                   const Expr& target,
                                                   int depth) {
    const bool at_inf = is_infinity(target);
    // Handles both x → ±∞ and a finite real point a (the same 1^∞ correction
    // appears at x → 0, e.g. ((1+x)^(1/x) − e)/x → −e/2).
    if (depth >= 8 || (!at_inf && !is_number(target))) return std::nullopt;
    // Gate to a variable power f(x)^{g(x)} sitting inside a sum (the f^g − c shape
    // the other asymptotic paths leave hanging), AND whose base tends to a finite
    // nonzero limit. The base condition is essential: x^(1/x) (base → ∞) becomes
    // log(1/u) = −log u after x = 1/u, a singular log series cannot expand — that is
    // the LIMIT-NOHANG-1 form, which must stay short-circuited. (1+1/x)^x has base
    // → 1, giving the analytic log(1+u). Other forms are left untouched.
    bool gated = false;
    auto g = [&](auto&& self, const Expr& e, bool under_add) -> void {
        if (gated) return;
        if (under_add && e->type_id() == TypeId::Pow
            && has(e->args()[0], var) && has(e->args()[1], var)) {
            Expr lb = limit_impl(e->args()[0], var, target, depth + 1);
            if (!is_nan(lb) && !is_infinity(lb) && !(lb == S::Zero())
                && lb->type_id() != TypeId::ComplexInfinity) {
                gated = true;
            }
            return;
        }
        const bool c = under_add || e->type_id() == TypeId::Add;
        for (const auto& a : e->args()) self(self, a, c);
    };
    g(g, expr, false);
    if (!gated) return std::nullopt;

    const Expr sgn = target->type_id() == TypeId::NegativeInfinity
                         ? S::NegativeOne()
                         : S::One();
    // A plain symbol (no positivity): with u marked positive the exp(g·log f)
    // rewrite below canonicalizes straight back to f^g, defeating the series step.
    // Local coordinate u → 0: x = ±1/u at ∞, x = a + u at a finite point a.
    Expr u = symbol("__series_u");
    Expr eu = at_inf ? subs(expr, var, mul(sgn, pow(u, integer(-1))))
                     : subs(expr, var, add(target, u));
    // Factor the explicit u-power pole eu = u^p · rest_raw. Pulling the pole out
    // keeps series(rest) regular — it stalls on a raw (1/u)·(nested-exp) form.
    Expr p = S::Zero();
    std::vector<Expr> rest_factors;
    if (eu->type_id() == TypeId::Mul) {
        for (const auto& f : eu->args()) {
            if (f == u) {
                p = add(p, S::One());
            } else if (f->type_id() == TypeId::Pow && f->args()[0] == u
                       && is_number(f->args()[1])) {
                p = add(p, f->args()[1]);
            } else {
                rest_factors.push_back(f);
            }
        }
    } else {
        rest_factors.push_back(eu);
    }
    Expr rest_raw = rest_factors.empty() ? Expr{S::One()} : mul(rest_factors);
    // Expansion order: enough to recover the u^{-p} coefficient (the term the pole
    // shifts to u⁰), plus a couple of spares. Capped — high orders both cost and
    // make the nested-exp expansions give up.
    std::size_t order = 4;
    {
        Expr pf = evalf(p, 20);
        if (pf->type_id() == TypeId::Float) {
            try {
                const long o = static_cast<long>(std::ceil(-std::stod(pf->str())))
                               + 3;
                order = static_cast<std::size_t>(std::max(3L, std::min(7L, o)));
            } catch (...) {
            }
        }
    }
    // Rewrite each u-dependent power f^g to exp(series(g·log f)). Pre-expanding the
    // exponent is essential: series can expand exp(polynomial) but stalls on
    // exp(log(cos u)/u) and similar nested-transcendental exponents.
    ExprMap<Expr> rwmap;
    bool rw_fail = false;
    auto rw = [&](auto&& self, const Expr& e) -> void {
        if (rw_fail) return;
        if (e->type_id() == TypeId::Pow && !rwmap.count(e)
            && has(e->args()[1], u)) {
            Expr hs;
            try {
                hs = series(mul(e->args()[1], log(e->args()[0])), u, S::Zero(),
                            order);
            } catch (...) {
                rw_fail = true;
                return;
            }
            bool bad = false;
            auto ck = [&](auto&& s, const Expr& x) -> void {
                if (bad) return;
                if (x->type_id() == TypeId::Function && has(x, u)) {
                    bad = true;
                    return;
                }
                for (const auto& a : x->args()) s(s, a);
            };
            ck(ck, hs);
            if (bad) {  // exponent did not reduce to a polynomial — abandon
                rw_fail = true;
                return;
            }
            rwmap.emplace(e, exp(hs));
        }
        for (const auto& a : e->args()) self(self, a);
    };
    rw(rw, rest_raw);
    if (rw_fail) return std::nullopt;
    Expr rest = rwmap.empty() ? rest_raw : xreplace(rest_raw, rwmap);
    Expr s;
    try {
        s = series(rest, u, S::Zero(), order);
    } catch (...) {
        return std::nullopt;
    }
    // Unify the constant's representations: series emits the lifted base as exp(k)
    // while the original carries e^k = E^k, and SymPP does not equate exp(2) with
    // E². A stray exp(k) − E^k would otherwise survive as a spurious leading term.
    {
        ExprMap<Expr> em;
        auto norm = [&](auto&& self, const Expr& e) -> void {
            if (e->type_id() == TypeId::Function && !em.count(e)
                && static_cast<const Function&>(*e).function_id()
                       == FunctionId::Exp
                && !has(e->args()[0], u)) {
                em.emplace(e, pow(S::E(), e->args()[0]));
            }
            for (const auto& a : e->args()) self(self, a);
        };
        norm(norm, s);
        if (!em.empty()) s = xreplace(s, em);
    }
    // The limit is the leading term of u^p · s. Read it directly off the truncated
    // series rather than calling limit() on a Laurent form (which a stray exp(k) −
    // E^k pole would hang, and whose one-sided u → 0⁺ sign limit() would lose).
    auto term_degree = [&](const Expr& t, Expr& coeff_out)
        -> std::optional<Expr> {
        Expr deg = S::Zero();
        std::vector<Expr> cfac;
        std::vector<Expr> facs;
        if (t->type_id() == TypeId::Mul) {
            for (const auto& f : t->args()) facs.push_back(f);
        } else {
            facs.push_back(t);
        }
        for (const auto& f : facs) {
            if (f == u) {
                deg = add(deg, S::One());
            } else if (f->type_id() == TypeId::Pow && f->args()[0] == u
                       && is_number(f->args()[1])) {
                deg = add(deg, f->args()[1]);
            } else if (has(f, u)) {
                return std::nullopt;  // non-monomial u dependence (e.g. log u)
            } else {
                cfac.push_back(f);
            }
        }
        coeff_out = cfac.empty() ? Expr{S::One()} : mul(cfac);
        return simplify(deg);
    };
    std::vector<Expr> terms;
    if (s->type_id() == TypeId::Add) {
        for (const auto& t : s->args()) terms.push_back(t);
    } else {
        terms.push_back(s);
    }
    Expr min_deg, min_coeff;
    double min_dv = 1e300;
    bool found = false;
    for (const auto& t : terms) {
        Expr c;
        auto d = term_degree(t, c);
        if (!d) return std::nullopt;
        Expr df = evalf(*d, 20);
        if (df->type_id() != TypeId::Float) return std::nullopt;
        double dv = 0.0;
        try {
            dv = std::stod(df->str());
        } catch (...) {
            return std::nullopt;
        }
        if (!found || dv < min_dv) {
            min_dv = dv;
            min_deg = *d;
            min_coeff = c;
            found = true;
        }
    }
    if (!found) return std::nullopt;
    auto sign_of = [](const Expr& e) -> int {
        if (is_positive(e) == std::optional<bool>{true}) return 1;
        if (is_negative(e) == std::optional<bool>{true}) return -1;
        Expr f = evalf(e, 30);
        if (f->type_id() == TypeId::Float) {
            try {
                const double d = std::stod(f->str());
                if (d > 1e-12) return 1;
                if (d < -1e-12) return -1;
            } catch (...) {
            }
        }
        return 0;
    };
    Expr order_total = simplify(add(p, min_deg));  // leading power of u^p · s
    Expr of = evalf(order_total, 20);
    if (of->type_id() != TypeId::Float) return std::nullopt;
    double ordv = 0.0;
    try {
        ordv = std::stod(of->str());
    } catch (...) {
        return std::nullopt;
    }
    Expr cand;
    if (ordv > 1e-9) {
        cand = S::Zero();  // vanishes
    } else if (std::fabs(ordv) < 1e-9) {
        cand = min_coeff;  // finite leading constant
    } else if (at_inf) {
        const int sc = sign_of(min_coeff);  // u^{negative} → +∞ as u → 0⁺
        if (sc > 0) cand = S::Infinity();
        else if (sc < 0) cand = S::NegativeInfinity();
        else return std::nullopt;
    } else {
        // A pole at a finite point is two-sided (the sign depends on the approach
        // direction and the power's parity); leave it to the other machinery.
        return std::nullopt;
    }
    if (has(cand, u) || is_nan(cand)) return std::nullopt;
    // Numeric guard: verify against the original at large |x|. High precision
    // absorbs the catastrophic cancellation in f(x) − c; the corrections converge
    // as ∼1/x, so the far point must be large.
    const bool cand_pinf = cand->type_id() == TypeId::Infinity;
    const bool cand_ninf = cand->type_id() == TypeId::NegativeInfinity;
    double cvd = 0.0;
    if (!cand_pinf && !cand_ninf) {
        Expr cvf = evalf(cand, 30);
        if (cvf->type_id() != TypeId::Float) return std::nullopt;
        try {
            cvd = std::stod(cvf->str());
        } catch (...) {
            return std::nullopt;
        }
    }
    // Sample points approaching the target: large |x| at ∞ (the corrections decay
    // as ∼1/x, so the far point must be large), or x = a + 1/N near a finite a.
    const std::array<long, 3> grid{300L, 3000L, 30000L};
    double prev_diff = 1e300, prev_val = 0.0;
    int checks = 0;
    for (long gp : grid) {
        const Expr pt = at_inf ? Expr{mul(sgn, integer(gp))}
                               : Expr{add(target, pow(integer(gp), integer(-1)))};
        Expr at = evalf(subs(expr, var, pt), 60);
        if (at->type_id() != TypeId::Float) return std::nullopt;
        double v = 0.0;
        try {
            v = std::stod(at->str());
        } catch (...) {
            return std::nullopt;
        }
        if (cand_pinf) {
            if (checks > 0 && !(v > prev_val)) return std::nullopt;
            if (gp == 30000L && !(v > 1.0)) return std::nullopt;
        } else if (cand_ninf) {
            if (checks > 0 && !(v < prev_val)) return std::nullopt;
            if (gp == 30000L && !(v < -1.0)) return std::nullopt;
        } else {
            const double d = std::fabs(v - cvd);
            if (!(d <= prev_diff + 1e-9)) return std::nullopt;
            prev_diff = d;
            if (gp == 30000L && !(d < 1e-2)) return std::nullopt;
        }
        prev_val = v;
        ++checks;
    }
    return cand;
}

// Continuity of a power with a constant exponent: lim base^r = (lim base)^r when
// the inner limit is determinate. Direct substitution computes the base
// pointwise and a non-rational base (e.g. log(log x)/log x) substitutes to an
// indeterminate ∞/∞ → nan, so √(…) of it surfaces as nan even though the base
// has a clean limit. Taking the base's limit first resolves it: √(log log x/log x)
// → √0 = 0. Guards skip a pole (0^negative) and a complex root of a negative base.
[[nodiscard]] std::optional<Expr> try_power_continuity(const Expr& expr,
                                                       const Expr& var,
                                                       const Expr& target,
                                                       int depth) {
    if (depth >= 12 || expr->type_id() != TypeId::Pow) return std::nullopt;
    const Expr& base = expr->args()[0];
    const Expr& ex = expr->args()[1];
    if (has(ex, var)) return std::nullopt;  // var exponent → not this rule
    Expr Lb = limit_impl(base, var, target, depth + 1);
    if (is_nan(Lb)) return std::nullopt;
    // A negative base with a non-integer exponent is complex near the limit; do
    // not fabricate a complex value for a real-limit query.
    if (ex->type_id() != TypeId::Integer
        && is_negative(Lb) == std::optional<bool>{true}) {
        return std::nullopt;
    }
    Expr r = simplify(pow(Lb, ex));
    if (is_nan(r) || r->type_id() == TypeId::ComplexInfinity) {
        return std::nullopt;
    }
    return r;
}

// True when f evaluates to a positive real at a large sample point — used to
// certify that a log-factor split log(∏fᵢ)=Σlog(fᵢ) is valid (every factor
// positive) for the asymptotic x→+∞ regime.
[[nodiscard]] bool sample_positive_inf(const Expr& f, const Expr& var) {
    Expr at = evalf(simplify(subs(f, var, integer(1000000))), 30);
    return at->type_id() == TypeId::Float
           && is_positive(at) == std::optional<bool>{true};
}

[[nodiscard]] std::optional<Expr> expand_log_positive(const Expr& arg,
                                                      const Expr& var,
                                                      int depth);

// Value of log(f) with f a single factor: log(exp g)=g, log(bᵉ)=e·log(b)
// (b expanded recursively), else log(f) when f is verified positive.
[[nodiscard]] std::optional<Expr> log_factor_positive(const Expr& f,
                                                      const Expr& var,
                                                      int depth) {
    if (depth > 6) return std::nullopt;
    if (f->type_id() == TypeId::Function
        && static_cast<const Function&>(*f).function_id() == FunctionId::Exp
        && f->args().size() == 1) {
        return f->args()[0];
    }
    if (f->type_id() == TypeId::Pow) {
        const Expr& base = f->args()[0];
        const Expr& e = f->args()[1];
        std::optional<Expr> lb = expand_log_positive(base, var, depth + 1);
        if (!lb) {
            if (!sample_positive_inf(base, var)) return std::nullopt;
            lb = log(base);
        }
        return mul(e, *lb);
    }
    // A sum factor (e.g. eᵍ − 1, which vanishes) is not a clean log-factor: its
    // log neither expands nor stays bounded, and taking log of a → 0 factor turns
    // the reduction into an expensive ∞ − ∞. Defer such products to the
    // small-angle / product paths rather than expanding here.
    if (f->type_id() == TypeId::Add) return std::nullopt;
    if (!sample_positive_inf(f, var)) return std::nullopt;
    return log(f);
}

// Expand log(arg) into a sum of logs, splitting only factors verified positive at
// +∞ (so the identity log(∏)=Σlog holds). Returns nullopt for an atomic arg with
// no expansion leverage.
[[nodiscard]] std::optional<Expr> expand_log_positive(const Expr& arg,
                                                      const Expr& var,
                                                      int depth) {
    if (depth > 6) return std::nullopt;
    if (arg->type_id() == TypeId::Function
        && static_cast<const Function&>(*arg).function_id() == FunctionId::Exp
        && arg->args().size() == 1) {
        return arg->args()[0];
    }
    if (arg->type_id() == TypeId::Mul) {
        std::vector<Expr> terms;
        terms.reserve(arg->args().size());
        for (const auto& f : arg->args()) {
            auto t = log_factor_positive(f, var, depth + 1);
            if (!t) return std::nullopt;
            terms.push_back(std::move(*t));
        }
        return add(std::move(terms));
    }
    if (arg->type_id() == TypeId::Pow) {
        return log_factor_positive(arg, var, depth);
    }
    return std::nullopt;
}

// Gruntz log-exp reduction: for a positive product/power/quotient at +∞ whose
// pointwise value is indeterminate, lim e = exp(lim log e). log e is expanded into
// a sum of logs (each split factor certified positive at large x, so the identity
// is exact in the asymptotic regime), whose limit the dominant-term / continuity
// machinery resolves; exp of that is the answer. Closes nested-transcendental
// ratios such as log x / exp(√(log x·log log x)) → 0 that the growth ranking
// misses. The positivity certificate keeps it sound, so it is strictly additive.
[[nodiscard]] std::optional<Expr> try_log_exp_reduction(const Expr& expr,
                                                        const Expr& var,
                                                        const Expr& target,
                                                        int depth) {
    if (depth >= 8 || target->type_id() != TypeId::Infinity) return std::nullopt;
    if (expr->type_id() != TypeId::Mul && expr->type_id() != TypeId::Pow) {
        return std::nullopt;
    }
    if (!has(expr, var) || !sample_positive_inf(expr, var)) return std::nullopt;
    auto le = expand_log_positive(expr, var, 0);
    if (!le) return std::nullopt;
    Expr m = limit_impl(*le, var, target, depth + 1);
    if (is_nan(m)) return std::nullopt;
    return exp_of_limit(m);
}

// Expand log(c·x), log(xᵏ), … subterms (argument a positive product/power) into a
// sum of logs and re-take, so an ∞−∞ of same-rank but differently-scaled log terms
// combines: 2x·log(2x) − x·log x = 2x·log 2 + x·log x → ∞ once log(2x) splits.
// (try_common_log_combine handles only matching coefficients — 2x vs x do not.)
[[nodiscard]] std::optional<Expr> rewrite_expand_logs(const Expr& expr,
                                                      const Expr& var,
                                                      const Expr& target,
                                                      int depth) {
    if (depth >= 10 || target->type_id() != TypeId::Infinity) {
        return std::nullopt;
    }
    ExprMap<Expr> m;
    auto scan = [&](auto&& self, const Expr& e) -> void {
        if (e->type_id() == TypeId::Function && !m.count(e)
            && static_cast<const Function&>(*e).function_id() == FunctionId::Log
            && e->args().size() == 1) {
            const Expr& arg = e->args()[0];
            if ((arg->type_id() == TypeId::Mul || arg->type_id() == TypeId::Pow)
                && has(arg, var)) {
                if (auto le = expand_log_positive(arg, var, 0)) {
                    if (!(*le == e)) m.emplace(e, *le);
                }
            }
        }
        for (const auto& a : e->args()) self(self, a);
    };
    scan(scan, expr);
    if (m.empty()) return std::nullopt;
    Expr r = limit_impl(expand(xreplace(expr, m)), var, target, depth + 1);
    if (is_nan(r)) return std::nullopt;
    return r;
}

// Gruntz dominant-term rule for an indeterminate ∞−∞ sum at ±∞: if one term
// strictly outgrows every other (each other term divided by it tends to 0), the
// sum is asymptotic to that term, so the limit equals the dominant term's limit.
// Resolves x − x·log x → −∞ (the −x·log x term dominates x), which closes the
// exp(x)/x^x chain once the powers are on an exp footing. Only fires for a unique
// strict dominator that itself diverges to a signed ∞, so genuinely cancelling
// differences (x − x, √(x²+x) − x) are left to the conjugate/algebraic machinery.
// Combine an ∞−∞ sum whose terms share a common (var-dependent) coefficient on
// a log: c·log(p) − c·log(q) = c·log(p/q). The standard log-combine only handles
// var-free coefficients, so a distributed exponent like x·log(x+1) − x·log(x)
// (which the power-as-exp/exp-combine path produces for x^x/(x+1)^x) was left as
// an unresolved ∞−∞ and sent the engine into the reciprocal-substitution loop.
// Combined, it is x·log((x+1)/x) → 1, resolved by the existing 0·∞ machinery.
[[nodiscard]] std::optional<Expr> try_common_log_combine(const Expr& expr,
                                                         const Expr& var,
                                                         const Expr& target,
                                                         int depth) {
    if (depth >= 10 || expr->type_id() != TypeId::Add
        || expr->args().size() < 2) {
        return std::nullopt;
    }
    // term = coeff · log(arg); returns nullopt if no single log factor.
    auto split = [](const Expr& term) -> std::optional<std::pair<Expr, Expr>> {
        std::vector<Expr> factors;
        if (term->type_id() == TypeId::Mul) {
            for (const auto& f : term->args()) factors.push_back(f);
        } else {
            factors.push_back(term);
        }
        Expr log_arg;
        std::vector<Expr> rest;
        for (const auto& f : factors) {
            if (!log_arg && f->type_id() == TypeId::Function
                && static_cast<const Function&>(*f).function_id()
                       == FunctionId::Log
                && f->args().size() == 1) {
                log_arg = f->args()[0];
            } else {
                rest.push_back(f);
            }
        }
        if (!log_arg) return std::nullopt;
        return std::make_pair(rest.empty() ? Expr{S::One()} : mul(rest), log_arg);
    };
    Expr ref_coeff;
    std::vector<std::pair<Expr, int>> logs;  // (arg, ±1 relative to ref_coeff)
    for (const auto& term : expr->args()) {
        auto sp = split(term);
        if (!sp) return std::nullopt;
        Expr coeff = simplify(sp->first);
        if (!has(coeff, var)) return std::nullopt;  // var-free → ordinary combine
        if (!ref_coeff) {
            ref_coeff = coeff;
            logs.emplace_back(sp->second, 1);
        } else if (coeff == ref_coeff) {
            logs.emplace_back(sp->second, 1);
        } else if (simplify(add(coeff, ref_coeff)) == S::Zero()) {
            logs.emplace_back(sp->second, -1);
        } else {
            return std::nullopt;  // coefficients not equal up to sign
        }
    }
    Expr prod = S::One();
    for (const auto& [arg, sign] : logs) prod = mul(prod, pow(arg, integer(sign)));
    Expr r = limit_impl(mul(ref_coeff, log(prod)), var, target, depth + 1);
    if (is_nan(r)) return std::nullopt;
    return r;
}

// Rational form of divergent terms at ±∞: P/Q where the denominator's dominant
// term D diverges. Divide numerator and denominator through by D, so each term
// tends to a finite value (0 for the sub-dominant ones, a constant for D's), and
// the quotient becomes finite/finite. This resolves ratios of exponentials —
// (eᵏˣ+1)/(eᵏˣ−1) → 1, (2ˣ+3ˣ)/3ˣ → 1 (after cˣ→exp rewrite) — that L'Hôpital
// loops on (the differentiated ratio never collapses) and that per-factor folding
// mis-resolves to 0. Dominance is decided by the limit of each pairwise ratio.
// Polynomial in var (var-free subterms are constant coefficients): the admissible
// exponent of a constant-base exponential c^g — excludes log/general transcendentals.
[[nodiscard]] bool is_poly_in_var(const Expr& e, const Expr& var) {
    if (!has(e, var)) return true;
    switch (e->type_id()) {
        case TypeId::Symbol:
            return true;
        case TypeId::Add:
        case TypeId::Mul:
            for (const auto& a : e->args()) {
                if (!is_poly_in_var(a, var)) return false;
            }
            return true;
        case TypeId::Pow:
            return is_poly_in_var(e->args()[0], var)
                   && !has(e->args()[1], var)
                   && e->args()[1]->type_id() == TypeId::Integer;
        default:
            return false;
    }
}

// A rational function of constant-base exponentials cˣ (c a positive number) and
// polynomials in var — the only shape try_dominant_ratio is sound and terminating
// on. Anything with a logarithm, a gamma, a general power f(x)^g(x), or an exp of a
// transcendental (e.g. exp(log x / x), into which x^(1/x) rewrites) is excluded, so
// the divide/distribute step never reopens a hard Gruntz form.
[[nodiscard]] bool is_exp_rational(const Expr& e, const Expr& var) {
    if (!has(e, var)) return true;
    switch (e->type_id()) {
        case TypeId::Symbol:
            return true;
        case TypeId::Add:
        case TypeId::Mul:
            for (const auto& a : e->args()) {
                if (!is_exp_rational(a, var)) return false;
            }
            return true;
        case TypeId::Pow: {
            const Expr& b = e->args()[0];
            const Expr& p = e->args()[1];
            if (!has(p, var) && p->type_id() == TypeId::Integer) {
                return is_exp_rational(b, var);  // f(x)ⁿ
            }
            const auto bt = b->type_id();
            const bool numeric_base = bt == TypeId::Integer
                                      || bt == TypeId::Rational
                                      || bt == TypeId::Float;
            return numeric_base
                   && is_positive(b) == std::optional<bool>{true}
                   && is_poly_in_var(p, var);  // c^poly, c > 0
        }
        case TypeId::Function:
            return static_cast<const Function&>(e.operator*()).function_id()
                       == FunctionId::Exp
                   && is_poly_in_var(e->args()[0], var);  // exp(poly)
        default:
            return false;
    }
}

[[nodiscard]] std::optional<Expr> try_dominant_ratio(const Expr& expr,
                                                     const Expr& var,
                                                     const Expr& target,
                                                     int depth) {
    if (depth >= 8 || !is_infinity(target)) return std::nullopt;
    if (!is_exp_rational(expr, var)) return std::nullopt;
    NumDen nd = split_after_together(expr);
    if (has(nd.den, var) && has(nd.num, var)) {
        std::vector<Expr> dterms;
        if (nd.den->type_id() == TypeId::Add) {
            for (const auto& t : nd.den->args()) dterms.push_back(t);
        } else {
            dterms.push_back(nd.den);
        }
        // Dominant term D of the denominator: all other terms / D → 0.
        Expr dominant;
        for (const auto& ti : dterms) {
            if (!has(ti, var)) continue;
            bool dom = true;
            for (const auto& tj : dterms) {
                if (tj.get() == ti.get()) continue;
                Expr ratio = limit_impl(simplify(mul(tj, pow(ti, integer(-1)))),
                                        var, target, depth + 1);
                if (!(ratio == S::Zero())) { dom = false; break; }
            }
            if (dom) { dominant = ti; break; }
        }
        if (dominant
            && is_infinity(limit_impl(dominant, var, target, depth + 1))) {
            Expr inv = pow(dominant, integer(-1));
            Expr lnum = limit_impl(simplify(mul(nd.num, inv)), var, target,
                                   depth + 1);
            Expr lden = limit_impl(simplify(mul(nd.den, inv)), var, target,
                                   depth + 1);
            if (!is_nan(lnum) && !is_nan(lden) && !is_infinity(lnum)
                && !is_infinity(lden) && !(lden == S::Zero())) {
                return simplify(mul(lnum, pow(lden, integer(-1))));
            }
        }
    }
    // A sum multiplied by a vanishing constant-base factor that `together` will
    // not pull into a denominator — (2ˣ+3ˣ)·3⁻ˣ — folds to 0·∞ as a product but
    // distributes to 2ˣ3⁻ˣ + 1 → 1, which per-term linearity resolves.
    if (expr->type_id() == TypeId::Mul) {
        bool has_add = false;
        for (const auto& f : expr->args()) {
            if (f->type_id() == TypeId::Add) { has_add = true; break; }
        }
        if (has_add) {
            Expr ex = expand(expr);
            if (!(ex == expr)) {
                Expr r = limit_impl(ex, var, target, depth + 1);
                if (!is_nan(r)) return r;
            }
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> try_dominant_term_sum(const Expr& expr,
                                                        const Expr& var,
                                                        const Expr& target,
                                                        int depth) {
    if (depth >= 10) return std::nullopt;
    if (expr->type_id() != TypeId::Add || !is_infinity(target)) {
        return std::nullopt;
    }
    const auto& terms = expr->args();
    if (terms.size() < 2 || terms.size() > 6) return std::nullopt;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        if (!has(terms[i], var)) continue;
        bool dominates = true;
        for (std::size_t j = 0; j < terms.size() && dominates; ++j) {
            if (j == i) continue;
            Expr ratio = simplify(mul(terms[j], pow(terms[i], integer(-1))));
            if (!(limit_impl(ratio, var, target, depth + 1) == S::Zero())) {
                dominates = false;
            }
        }
        if (!dominates) continue;
        Expr cl = limit_impl(terms[i], var, target, depth + 1);
        if (cl->type_id() == TypeId::Infinity
            || cl->type_id() == TypeId::NegativeInfinity) {
            return cl;
        }
    }
    return std::nullopt;
}

// Replace every harmonic number H(g) whose argument g → +∞ by its asymptotic
// expansion H(g) ~ log g + γ + 1/(2g) − 1/(12g²). The harmonic number is opaque
// to the limit machinery (H(∞) does not fold), which left H(n)/log n → 0 (wrong)
// and H(n) → nan. Four terms suffice for the usual limits — H(n) → ∞,
// H(n)/log n → 1, H(n) − log n → γ, n·(H(n) − log n − γ) → 1/2 — and the
// substitution is gated on g provably diverging, so a finite-argument harmonic is
// untouched.
[[nodiscard]] std::optional<Expr> rewrite_harmonic_asymptotic(const Expr& expr,
                                                              const Expr& var,
                                                              const Expr& target,
                                                              int depth) {
    if (depth >= 10 || target->type_id() != TypeId::Infinity) {
        return std::nullopt;
    }
    ExprMap<Expr> m;
    auto scan = [&](auto&& self, const Expr& e) -> void {
        if (e->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*e);
            if (fn.function_id() == FunctionId::Harmonic && fn.args().size() == 1
                && !m.count(e)) {
                const Expr& g = fn.args()[0];
                if (has(g, var)
                    && limit_impl(g, var, target, depth + 1) == S::Infinity()) {
                    m.emplace(
                        e, add({log(g), S::EulerGamma(),
                                mul(rational(1, 2), pow(g, integer(-1))),
                                mul(rational(-1, 12), pow(g, integer(-2)))}));
                }
            }
        }
        for (const auto& a : e->args()) self(self, a);
    };
    scan(scan, expr);
    if (m.empty()) return std::nullopt;
    // Expand the substituted form so a difference of harmonics flattens — without
    // it, H(2n) − H(n) leaves an undistributed −1·(log n + γ + …), where the γ does
    // not cancel and the log(2n) − log(n) sits behind a product, sending the limit
    // engine down a non-terminating path. Expanded, it is a flat sum the
    // log-combine / vanishing-tail rules resolve (H(2n) − H(n) → log 2).
    return limit_impl(expand(xreplace(expr, m)), var, target, depth + 1);
}

// Log-Stirling asymptotic at +∞. log Γ(z) = (z−½)·log z − z + ½·log 2π + o(1), so
//   log(g!)   = log Γ(g+1) ~ (g+½)·log g − g + ½·log 2π   (the log(g+1) folds away),
//   log Γ(g)               ~ (g−½)·log g − g + ½·log 2π.
// Replacing each log of a divergent factorial/gamma by its expansion turns the
// expression elementary (the dropped o(1) vanishes, so the limit is exact). This
// resolves log(n!)/(n log n) → 1, log(n!)/n → ∞ (was a wrong 0), log(n!) − n log n
// → −∞, and log(n!)/log(nⁿ) → 1 (which previously hung).
[[nodiscard]] std::optional<Expr> rewrite_loggamma_asymptotic(const Expr& expr,
                                                              const Expr& var,
                                                              const Expr& target,
                                                              int depth) {
    if (depth >= 10 || target->type_id() != TypeId::Infinity) {
        return std::nullopt;
    }
    const Expr half_log2pi =
        mul(rational(1, 2), log(mul(integer(2), S::Pi())));
    ExprMap<Expr> m;
    auto scan = [&](auto&& self, const Expr& e) -> void {
        if (e->type_id() == TypeId::Function && !m.count(e)) {
            const auto& fn = static_cast<const Function&>(*e);
            if (fn.function_id() == FunctionId::Log && fn.args().size() == 1
                && fn.args()[0]->type_id() == TypeId::Function) {
                const auto& inner =
                    static_cast<const Function&>(*fn.args()[0]);
                const auto id = inner.function_id();
                if ((id == FunctionId::Factorial || id == FunctionId::Gamma)
                    && inner.args().size() == 1) {
                    const Expr& g = inner.args()[0];
                    if (has(g, var)
                        && limit_impl(g, var, target, depth + 1)
                               == S::Infinity()) {
                        // z·log z − z form. factorial(g) = Γ(g+1) carries the +½
                        // coefficient on log g; Γ(g) carries −½. A gamma with a
                        // positive-integer shift, Γ(var+k) = (var+k−1)!, is recast as
                        // the factorial of var+k−1 so its log argument stays var-clean
                        // (Γ(n+1) ⇒ log n, matching n!) rather than log(n+1) — on which
                        // the engine's polynomial·log-ratio handling stalls.
                        Expr z = g, coeff = rational(1, 2);
                        if (id == FunctionId::Gamma) {
                            Expr s =
                                simplify(add(g, mul(S::NegativeOne(), var)));
                            if (s->type_id() == TypeId::Integer
                                && is_positive(s) == std::optional<bool>{true}) {
                                z = simplify(add(g, S::NegativeOne()));
                            } else {
                                coeff = rational(-1, 2);
                            }
                        }
                        m.emplace(e, add({mul(add(z, coeff), log(z)),
                                          mul(S::NegativeOne(), z),
                                          half_log2pi}));
                    }
                }
            }
        }
        for (const auto& a : e->args()) self(self, a);
    };
    scan(scan, expr);
    if (m.empty()) return std::nullopt;
    // Expand so a difference such as log(n!) − n·log n flattens to a sum the
    // log-combine / dominant-term rules resolve, rather than hiding behind a
    // product as it would after a bare substitution.
    return limit_impl(expand(xreplace(expr, m)), var, target, depth + 1);
}

// Verify a candidate limit numerically: sample the ORIGINAL at large +x and require
// convergence to a finite candidate / consistent divergence for ±∞.
[[nodiscard]] bool numeric_verify_at_inf(const Expr& expr, const Expr& var,
                                         const Expr& cand) {
    const bool pinf = cand->type_id() == TypeId::Infinity;
    const bool ninf = cand->type_id() == TypeId::NegativeInfinity;
    double cvd = 0.0;
    if (!pinf && !ninf) {
        Expr cvf = evalf(cand, 30);
        if (cvf->type_id() != TypeId::Float) return false;
        try {
            cvd = std::stod(cvf->str());
        } catch (...) {
            return false;
        }
    }
    // Wide sample span: asinh/log-type corrections converge only as ∼1/log x, far
    // too slow for an absolute tolerance, so a finite candidate is accepted on a
    // strictly *shrinking* gap that decays substantially across the span — a trend a
    // wrong (e.g. precision-truncated) candidate, whose gap stays O(1), fails.
    double first_diff = -1.0, prev_diff = 1e300, prev_val = 0.0;
    int checks = 0;
    for (long g : {10L, 1000L, 100000L}) {
        Expr at = evalf(subs(expr, var, integer(g)), 50);
        if (at->type_id() != TypeId::Float) return false;
        double v = 0.0;
        try {
            v = std::stod(at->str());
        } catch (...) {
            return false;
        }
        if (pinf) {
            if (checks > 0 && !(v > prev_val)) return false;
            if (g == 100000L && !(v > 1.0)) return false;
        } else if (ninf) {
            if (checks > 0 && !(v < prev_val)) return false;
            if (g == 100000L && !(v < -1.0)) return false;
        } else {
            const double d = std::fabs(v - cvd);
            if (checks == 0) first_diff = d;
            if (checks > 0 && !(d < prev_diff)) return false;  // strictly shrinking
            prev_diff = d;
            if (g == 100000L && !(d < 0.15 && d < 0.5 * first_diff)) return false;
        }
        prev_val = v;
        ++checks;
    }
    return checks == 3;
}

// Inverse-hyperbolic asymptotics at +∞: asinh(g), acosh(g) ~ log(2g) as g → +∞
// (both equal log(g + √(g²±1)), and √(g²±1) ~ g). The engine hangs on the exact
// log-of-a-radical form, so the leading log(2g) is substituted and the re-taken
// limit numerically verified — resolving asinh(x)/log x → 1, acosh(x)/log x → 1.
[[nodiscard]] std::optional<Expr> try_inverse_hyperbolic_asymptotic(
    const Expr& expr, const Expr& var, const Expr& target, int depth) {
    if (depth >= 10 || target->type_id() != TypeId::Infinity) {
        return std::nullopt;
    }
    ExprMap<Expr> m;
    auto scan = [&](auto&& self, const Expr& e) -> void {
        if (e->type_id() == TypeId::Function && !m.count(e)) {
            const auto id = static_cast<const Function&>(*e).function_id();
            if ((id == FunctionId::Asinh || id == FunctionId::Acosh)
                && e->args().size() == 1 && has(e->args()[0], var)
                && limit_impl(e->args()[0], var, target, depth + 1)
                       == S::Infinity()) {
                // Two-term asymptotic, so a difference that cancels the leading
                // log(2g) still resolves: asinh(g) = log(2g) + 1/(4g²) + O(g⁻⁴),
                // acosh(g) = log(2g) − 1/(4g²) + O(g⁻⁴), giving
                // (acosh − asinh)·g² → −1/2. The numeric guard rejects a case that
                // needs the dropped O(g⁻⁴) term.
                const Expr& g = e->args()[0];
                const Expr corr = id == FunctionId::Asinh ? Expr{rational(1, 4)}
                                                          : Expr{rational(-1, 4)};
                m.emplace(e, add(log(mul(integer(2), g)),
                                 mul(corr, pow(g, integer(-2)))));
            }
        }
        for (const auto& a : e->args()) self(self, a);
    };
    scan(scan, expr);
    if (m.empty()) return std::nullopt;
    Expr cand = limit_impl(xreplace(expr, m), var, target, depth + 1);
    if (is_nan(cand) || has(cand, var)
        || cand->type_id() == TypeId::ComplexInfinity) {
        return std::nullopt;
    }
    return numeric_verify_at_inf(expr, var, cand)
               ? std::optional<Expr>{cand}
               : std::nullopt;
}

// Leading-term (small-angle) substitution: every f(g) with f(t) = t + O(t³) at 0
// — sin, tan, sinh, tanh, asin, atan, asinh, atanh — and g → 0 is replaced by its
// argument g, after which the limit is re-taken. This resolves the 0·∞ forms the
// heuristic engine abandons, e.g. eˣ·sin(e⁻ˣ) → 1 and the canonical Gruntz
// oscillation eˣ·(sin(1/x + e⁻ˣ) − sin(1/x)) → 1. The substitution drops the
// cubic tail, so the candidate is accepted only after a numeric check against the
// original at large |x| (a wrong value keeps an O(1) offset and is rejected; the
// non-increasing ladder tolerates the ∼1/x² approach of these limits).
[[nodiscard]] std::optional<Expr> try_small_angle(const Expr& expr,
                                                  const Expr& var,
                                                  const Expr& target,
                                                  int depth) {
    // Only a product can give the 0·∞ these leading terms resolve (a diverging
    // factor times a vanishing f(g)). Restricting to Mul keeps every target while
    // skipping the many other nan-branch expressions that merely *contain* an
    // f(g → 0) — important now that eᵍ and cos g (common subterms) are expanded.
    if (depth >= 8 || !is_infinity(target) || expr->type_id() != TypeId::Mul
        || !has(expr, var)) {
        return std::nullopt;
    }
    ExprMap<Expr> m;
    auto scan = [&](auto&& self, const Expr& e) -> void {
        if (e->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*e);
            const FunctionId id = fn.function_id();
            const bool linear_at_zero =
                id == FunctionId::Sin || id == FunctionId::Tan
                || id == FunctionId::Sinh || id == FunctionId::Tanh
                || id == FunctionId::Asin || id == FunctionId::Atan
                || id == FunctionId::Asinh || id == FunctionId::Atanh;
            // Functions that are 1 + O(g) at 0 are replaced by their Maclaurin
            // head (through g⁴), so eᵍ − 1 ~ g and cos g − 1 ~ −g²/2 resolve too.
            const bool unit_at_zero =
                id == FunctionId::Exp || id == FunctionId::Cos
                || id == FunctionId::Cosh;
            if ((linear_at_zero || unit_at_zero) && fn.args().size() == 1
                && !m.count(e)) {
                const Expr& g = fn.args()[0];
                if (has(g, var)
                    && limit_impl(g, var, target, depth + 1) == S::Zero()) {
                    if (linear_at_zero) {
                        m.emplace(e, g);
                    } else if (id == FunctionId::Exp) {
                        // 1 + g + g²/2 (eᵍ − 1 ~ g); the numeric check rejects any
                        // case that would need a higher term.
                        m.emplace(e, add({S::One(), g,
                                          mul(rational(1, 2),
                                              pow(g, integer(2)))}));
                    } else {  // Cos / Cosh: 1 ∓ g²/2 (cos g − 1 ~ −g²/2). Kept to
                        // the quadratic head so a sum-valued g does not blow up
                        // expand() with a g⁴ term.
                        const Expr s2 = id == FunctionId::Cos ? rational(-1, 2)
                                                             : rational(1, 2);
                        m.emplace(e, add({S::One(),
                                          mul(s2, pow(g, integer(2)))}));
                    }
                }
            }
        }
        for (const auto& a : e->args()) self(self, a);
    };
    scan(scan, expr);
    if (m.empty()) return std::nullopt;
    Expr cand = limit_impl(expand(xreplace(expr, m)), var, target, depth + 1);
    if (is_nan(cand) || is_infinity(cand) || has(cand, var)
        || cand->type_id() == TypeId::ComplexInfinity) {
        return std::nullopt;
    }
    Expr cv = evalf(cand, 30);
    if (!is_number(cv)) return std::nullopt;  // a real numeric candidate
    double cvd = 0.0;
    try {
        cvd = std::stod(cv->str());  // "1" (int) or "0.5" (float) both parse
    } catch (...) {
        return std::nullopt;
    }
    const Expr sgn = target->type_id() == TypeId::Infinity ? S::One()
                                                           : S::NegativeOne();
    double prev = 1e300;
    int checks = 0;
    // Sample points stay moderate so a big exponential (e^x) does not overflow and
    // the working precision (which scales with the point to resolve a difference
    // sin(a+h)−sin(a) with h ∼ e^{−x} that would otherwise be lost to catastrophic
    // cancellation) stays affordable, yet are large enough that an ∼1/x² approach
    // is already within tolerance (at x = 300, 1/(2x²) ≈ 6e−6 < 1e−4).
    for (long xv : {100L, 200L, 300L}) {
        const int prec = static_cast<int>(static_cast<double>(xv) * 0.45) + 40;
        Expr at = evalf(simplify(subs(expr, var, mul(sgn, integer(xv)))), prec);
        if (!is_number(at)) return std::nullopt;
        try {
            const double d = std::fabs(std::stod(at->str()) - cvd);
            if (!(d <= prev + 1e-12)) return std::nullopt;  // must not diverge
            prev = d;
            ++checks;
        } catch (...) {
            return std::nullopt;
        }
    }
    if (checks == 3 && prev < 1e-4) return cand;
    return std::nullopt;
}

Expr limit_impl(const Expr& expr, const Expr& var, const Expr& target,
                int depth) {
    // Canonicalize (eᵃ)ⁿ → eⁿᵃ up front: otherwise a reciprocal of a vanishing
    // exponential, (e⁻ˣ)⁻¹, substitutes to 0⁻¹ = zoo instead of resolving to eˣ.
    if (depth < 14) {
        Expr fe = flatten_exp_powers(expr);
        if (!(fe == expr)) return limit_impl(fe, var, target, depth + 1);
    }
    // Reciprocal trig/hyperbolic functions are opaque to the limit machinery;
    // rewrite them as sin/cos ratios (cot x → cos x/sin x, …) and retry, so
    // forms like x·cot(x) resolve via L'Hôpital instead of returning nan.
    if (depth < 12) {
        Expr rw = rewrite_reciprocal_trig(expr);
        if (!(rw == expr)) return limit_impl(rw, var, target, depth + 1);
    }
    // Hyperbolic functions at ±∞ → exponential form, so sinh/cosh combinations
    // (sinh u + cosh u → eᵘ, sinh u/cosh u → 1) resolve via the exp machinery.
    // Expand so the eᵘ terms cancel where they should — cosh u − sinh u collapses
    // to e⁻ᵘ. The (eᵃ)ⁿ flattening at the top of limit_impl then turns the
    // reciprocal of a single-exponential denominator into a clean eⁿᵃ.
    if (depth < 12 && is_infinity(target)) {
        Expr rw = rewrite_hyperbolic_exp(expr, var, target, depth);
        if (!(rw == expr)) {
            return limit_impl(expand(rw), var, target, depth + 1);
        }
    }
    // Harmonic numbers H(g), g → +∞, expand to log g + γ + 1/(2g) − …; opaque
    // otherwise (H(n)/log n → 0 wrong, H(n) → nan).
    if (auto v = rewrite_harmonic_asymptotic(expr, var, target, depth)) {
        return *v;
    }
    // log of a divergent factorial/gamma → its log-Stirling expansion; opaque
    // otherwise (log(n!)/n → 0 wrong, log(n!)/log(nⁿ) hung).
    if (auto v = rewrite_loggamma_asymptotic(expr, var, target, depth)) {
        return *v;
    }
    // asinh(g), acosh(g) ~ log(2g) at +∞ — the exact log-of-radical form hangs.
    if (auto v = try_inverse_hyperbolic_asymptotic(expr, var, target, depth)) {
        return *v;
    }
    // An expression with sign(g), g → 0, is discontinuous at the target; resolve
    // it from the one-sided signs rather than the point value sign(0)=0.
    if (depth < 12) {
        if (auto s = resolve_sign_limit(expr, var, target, depth)) return *s;
    }
    // A rational function at ±∞: degree/leading-coefficient comparison, done
    // before the unreliable ∞ substitution and the symbolic-coefficient-blind
    // L'Hôpital path.
    if (is_infinity(target)) {
        if (auto r = rational_limit_at_infinity(expr, var, target)) return *r;
        if (depth < 12) {
            if (auto r = bounded_times_vanishing(expr, var, target, depth)) {
                return *r;
            }
            if (auto r = try_oscillating_rational(expr, var, target, depth)) {
                return *r;
            }
        }
    }

    // A balanced gamma ratio (Σ exponents = 0) resolves by the Stirling
    // asymptotic Γ(x+a)/Γ(x) ~ x^a. Tried first: it is exact and cheap, and for
    // half-integer shifts it also pre-empts the Stirling-root numeric guard
    // below, which does not terminate on some such shapes (e.g.
    // Γ(x+3/2)/Γ(x)/x^(3/2)). The growth-rank test must abstain here (same rank,
    // opposite direction), and a later path returns a wrong 0.
    if (target == S::Infinity() && count_gamma_factorial(expr) > 0) {
        if (auto r = gamma_ratio_asymptotic(normalize_factorial(expr), var))
            return *r;
        // A multiple-rate Γ(kx+b) is opaque to the slope-1 asymptotic above;
        // Gauss multiplication splits it into k slope-1 gammas plus a
        // constant-base exponential (e.g. Γ(2x+1)/Γ(x+1)²/4ˣ → 0, Γ(3x)/Γ(x)³ → ∞).
        if (auto r = try_gamma_multiplication(expr, var, target, depth)) return *r;
    }

    // A super-power n^(c·n) against a single factorial: n!/n^n → 0, n^n/n! → ∞.
    if (target == S::Infinity() && count_gamma_factorial(expr) > 0) {
        if (auto r = superpow_vs_factorial(expr, var)) return *r;
    }

    // n-th roots of factorial/gamma via Stirling (numerically guarded):
    // (n!)^(1/n)/n → 1/e, (n!)^(1/n) → ∞, n/(n!)^(1/n) → e. Skipped on a top-level
    // sum: the Stirling rewrite of a difference of comparable divergent terms
    // (n! − nⁿ) just yields another ∞−∞ form the substitution spins on — that is a
    // dominant-term case, handled below.
    if (depth < 10 && target == S::Infinity()
        && count_gamma_factorial(expr) > 0 && has_var_radical(expr, var)
        && expr->type_id() != TypeId::Add) {
        if (auto r = try_stirling_limit(expr, var, target, depth)) return *r;
    }

    // Gamma/factorial at +∞. Direct substitution gives gamma(∞)/gamma(∞), which
    // simplify wrongly cancels to 1; handle these before the substitution below.
    if (depth < 12 && target == S::Infinity()
        && count_gamma_factorial(expr) > 0) {
        // Normalize factorial → gamma so gamma(x)/factorial(x) can collapse.
        Expr g = normalize_factorial(expr);
        // (B) Collapse gamma ratios first — simplify reduces gamma(x+1)/gamma(x)
        // to x, after which the rational-at-∞ machinery applies. Only recurse
        // when simplification actually removed gamma content.
        Expr simp = simplify(g);
        if (count_gamma_factorial(simp) < count_gamma_factorial(g)) {
            return limit_impl(simp, var, target, depth + 1);
        }
        // (C) Otherwise decide by the dominant growth class (gamma ≫ exp ≫
        // poly ≫ log): e.g. exp(x)/gamma(x) → 0, gamma(x)/exp(x) → ∞.
        if (auto r = gamma_growth_limit(g, var, depth)) return *r;
        // Fall back to the gamma form (if normalization changed anything) so the
        // remaining methods never differentiate factorial into a Derivative node.
        if (!(g == expr)) return limit_impl(g, var, target, depth + 1);
    }

    // Indeterminate power forms (1^∞, 0^0, ∞^0) must be resolved before direct
    // substitution, which would collapse 1^∞ to 1 — e.g. (1+x)^(1/x) → E (not 1)
    // as x → 0. try_power_form returns nullopt for any determinate power.
    if (depth < 12 && expr->type_id() == TypeId::Pow) {
        if (auto v = try_power_form(expr, var, target, depth)) return *v;
    }

    // Continuity of atan/log: limit(f(g)) = f(limit g). Done before direct
    // substitution, which can leave atan(<unevaluated ∞ expression>) — e.g. the
    // atan terms of the ∫1/(1+x⁴) antiderivative at ∞.
    if (depth < 12 && expr->type_id() == TypeId::Function) {
        if (auto v = try_log_limit(expr, var, target, depth)) return *v;
    }

    // log of an exponential-dominated sum (log(2^x+3^x)/x → log 3): resolve before
    // direct substitution, which folds the inner log(∞) into an unevaluated
    // ∞-arithmetic mess rather than nan, so the post-substitution path can't catch
    // it. The handler bails fast unless such a log is present.
    if (depth < 12) {
        if (auto v = try_log_exp_asymptotic(expr, var, target, depth)) return *v;
    }

    // Linearity over a sum, attempted before direct substitution: when every
    // term has a determinate finite limit, the limit is their sum. Hoisted ahead
    // of substitution and L'Hôpital so a convergent special function plus a
    // vanishing term — e.g. the Si(2x) + cos(2x)/(2x) − 1/(2x) antiderivative of
    // sin²x/x² — isn't mangled: substitution would fold Si(∞) to a wrong value
    // and L'Hôpital would differentiate Si into sin(x)/x and substitute ∞.
    if (depth < 12 && expr->type_id() == TypeId::Add) {
        std::vector<Expr> finite_limits;
        std::vector<Expr> divergent_terms;
        bool any_nan = false;
        for (const auto& t : expr->args()) {
            Expr tl = limit_impl(t, var, target, depth + 1);
            if (is_nan(tl)) { any_nan = true; break; }
            if (is_infinity(tl)) {
                divergent_terms.push_back(t);
            } else {
                finite_limits.push_back(std::move(tl));
            }
        }
        if (!any_nan) {
            if (divergent_terms.empty()) {
                return simplify(add(std::move(finite_limits)));
            }
            // Mixed: finite terms alongside divergent ones whose ∞ − ∞ may
            // cancel. Peel the finite terms off and resolve the divergent
            // remainder on its own; if it has a determinate limit, add them.
            // (A lone divergent term, or an all-divergent sum, is left to the
            // ∞ − ∞ machinery below.)
            if (!finite_limits.empty() && divergent_terms.size() >= 2) {
                Expr rem = limit_impl(add(std::move(divergent_terms)), var,
                                      target, depth + 1);
                if (!is_nan(rem)) {
                    finite_limits.push_back(std::move(rem));
                    return simplify(add(std::move(finite_limits)));
                }
            }
        }
    }

    // Pull var-free constant factors out of a product before substitution:
    // lim c·g(x) = c·lim g(x). Done early so c·(convergent special-function sum)
    // isn't collapsed by the fragile substitution / L'Hôpital paths — e.g.
    // −1·(−Si(x) − cos(x)/x) → π/2 (the antiderivative of (1−cos x)/x²). Only
    // taken when the inner limit is determinate (not nan).
    if (depth < 12 && expr->type_id() == TypeId::Mul) {
        std::vector<Expr> const_factors;
        std::vector<Expr> var_factors;
        for (const auto& f : expr->args()) {
            (has(f, var) ? var_factors : const_factors).push_back(f);
        }
        if (!const_factors.empty() && !var_factors.empty()) {
            Expr inner = limit_impl(mul(std::move(var_factors)), var, target,
                                    depth + 1);
            if (!is_nan(inner)) {
                return simplify(mul(mul(std::move(const_factors)), inner));
            }
        }
    }

    Expr direct = simplify(subs(expr, var, target));

    // A finite-target pole surfaces as zoo; resolve its sign when both sides
    // agree, otherwise leave the unsigned zoo.
    if (direct->type_id() == TypeId::ComplexInfinity) {
        if (auto s = signed_pole(expr, var, target)) return *s;
    }

    if (depth < 12) {
        // Indeterminate forms surface as nan after substitution.
        if (is_nan(direct)) {
            // Polynomial at ±∞: the highest-degree term dominates, resolving the
            // ∞−∞ that direct substitution leaves as nan (x²−x → +∞, x−x² → −∞).
            if (is_infinity(target)) {
                try {
                    Poly p(expand(expr), var);
                    bool poly_ok = p.degree() >= 1;
                    for (const auto& cc : p.coeffs()) {
                        if (has(cc, var)) { poly_ok = false; break; }
                    }
                    if (poly_ok) {
                        Expr lead = mul(
                            p.leading_coeff(),
                            pow(var, integer(static_cast<long>(p.degree()))));
                        Expr lim = simplify(subs(lead, var, target));
                        if (!is_nan(lim)) return lim;
                    }
                } catch (const std::exception&) {
                    // not a polynomial in var — fall through to other methods
                }
                // ∞ − ∞ with a radical: rationalize via the conjugate.
                if (auto v = try_conjugate_difference(expr, var, target, depth)) {
                    return *v;
                }
                // Algebraic ∞−∞ / 0·∞ via the leading asymptotic term — resolves the
                // radical ratios the conjugate produces and L'Hôpital abandons.
                if (auto v = try_algebraic_inf(expr, var, target)) {
                    return *v;
                }
                // Gruntz dominant-term rule: a ∞−∞ sum with a unique strictly
                // dominant divergent term (e.g. x − x·log x → −∞).
                if (auto v = try_dominant_term_sum(expr, var, target, depth)) {
                    return *v;
                }
            }
            // Logarithms: log(g) → log(lim g), and combine a ∞ − ∞ between logs.
            if (auto v = try_log_limit(expr, var, target, depth)) return *v;
            // c·log(p) − c·log(q) with a var coefficient c → c·log(p/q); resolves
            // the distributed exponent x·log(x+1) − x·log(x) (from x^x/(x+1)^x).
            if (auto v = try_common_log_combine(expr, var, target, depth)) {
                return *v;
            }
            // Expand log(c·x), log(xᵏ) subterms so same-rank, differently-scaled log
            // terms combine: 2x·log(2x) − x·log x → ∞.
            if (auto v = rewrite_expand_logs(expr, var, target, depth)) return *v;
            // Rational form of divergent terms P/Q: divide through by Q's dominant
            // term so the quotient is finite/finite — (2ˣ+1)/(2ˣ−1) → 1,
            // (2ˣ+3ˣ)/3ˣ → 1 — which L'Hôpital loops on and folding mis-reads as 0.
            if (auto v = try_dominant_ratio(expr, var, target, depth)) {
                return *v;
            }
            // (Sum linearity for the all-finite case is handled before direct
            // substitution above; a genuine ∞ − ∞ falls through to L'Hôpital on
            // the combined fraction.)
            // Determinate product: fold the per-factor limits when there is no
            // genuine 0·∞ conflict (e.g. 2·(…→0) → 0). Direct substitution gives
            // nan when a nested factor's substitution does, even though the
            // product is determinate. The 0·∞ case is left to try_product_form.
            if (expr->type_id() == TypeId::Mul) {
                std::vector<Expr> factor_limits;
                bool ok = true;
                int zeros = 0, infs = 0;
                for (const auto& f : expr->args()) {
                    Expr fl = limit_impl(f, var, target, depth + 1);
                    if (is_nan(fl)) { ok = false; break; }
                    if (fl == S::Zero()) ++zeros;
                    else if (is_infinity(fl)) ++infs;
                    factor_limits.push_back(std::move(fl));
                }
                if (ok && !(zeros > 0 && infs > 0)) {
                    Expr p = simplify(mul(std::move(factor_limits)));
                    if (!is_nan(p)) return p;
                }
            }
            if (auto v = try_power_form(expr, var, target, depth)) return *v;
            // Merge exponentials with *different* exponent monomials
            // (exp(x²)·exp(−2x), exp(x²)/exp(x)²) into one exp(Σ) and re-take —
            // try_exponential_product below only combines a shared monomial.
            if (auto v = try_exp_combine(expr, var, target, depth)) return *v;
            // Merge a product of constant-base exponentials (2^x/3^x,
            // exp(x)/exp(2x)) into one exp(rate) before the generic product path,
            // which would otherwise see ∞·0 and stall in L'Hôpital → nan.
            if (auto v = try_exponential_product(expr, var, target, depth)) {
                return *v;
            }
            if (auto v = try_product_form(expr, var, target, depth)) return *v;
            // Gruntz MRV rewrite: a difference of asymptotically-equal exponentials,
            // e^{x+e⁻ˣ} − e^x → 1, factored to e^x·(e^{e⁻ˣ} − 1).
            if (auto v = try_exp_difference(expr, var, target, depth)) return *v;
            // A unit-tending power f(x)^g(x) → 1 in a difference: expand exp(g·log f)
            // − 1 = g·log f + …, resolving (x^x − 1)/(x·log x) → 1 where the log
            // singularity in the exponent defeats the polynomial series stage.
            if (auto v = try_unit_power_expansion(expr, var, target, depth)) {
                return *v;
            }
            // Gruntz leading-term via series for an f(x)^g(x) − c difference, before
            // the power-as-exp rewrite (which would otherwise spin on the resulting
            // ∞·0): x·((1+1/x)^x − e) → −e/2.
            if (auto v = try_series_limit(expr, var, target, depth)) return *v;
            // Gruntz: rewrite general powers f(x)^g(x) as exp(g·log f) and
            // re-take, so the exp/gamma growth machinery can resolve forms the
            // bare-power paths leave as nan — e.g. Γ(2n)/n^n → ∞.
            if (auto v = try_power_as_exp(expr, var, target, depth)) return *v;
            // Power with a constant exponent: lim base^r = (lim base)^r — resolves
            // √(non-rational base) whose pointwise substitution is nan.
            if (auto v = try_power_continuity(expr, var, target, depth)) {
                return *v;
            }
            // Gruntz log-exp reduction: lim e = exp(lim log e) for a positive
            // product/power at +∞ — closes nested-transcendental ratios like
            // log x / exp(√(log x·log log x)) → 0.
            if (auto v = try_log_exp_reduction(expr, var, target, depth)) {
                return *v;
            }
            // Small-angle: replace sin/tan/…(g) by g when g → 0 (numerically
            // verified) — closes eˣ·sin(e⁻ˣ) → 1 and the Gruntz oscillation.
            if (auto v = try_small_angle(expr, var, target, depth)) return *v;
        }
        // 0/0 and ∞/∞ quotients (also recovers finite 0/0 where direct
        // substitution collapses to 0 or nan). A nan result is not an answer —
        // fall through to the reciprocal substitution rather than returning it.
        // An oscillating f(∞) (lhopital's determinate-denominator branch divides
        // sin(x)/1 → sin(∞)) is likewise not an answer — let the end guard nan it.
        if (auto v = lhopital(expr, var, target);
            v && !is_nan(*v) && !has_oscillating_infinity(*v)
            && !has_buried_infinity(*v)) {
            return *v;
        }

        // Reciprocal substitution for a limit at ±∞: x = ±1/t maps it to t → 0,
        // where the series / L'Hôpital machinery resolves the asymptotic-expansion
        // forms the direct ∞ methods abandon — cube-root (and general n-th-root)
        // differences (x³+x²)^(1/3)−x → 1/3, and 0·∞ / ∞−∞ with a transcendental
        // subleading term: x²(1−cos 1/x) → 1/2, x − x²·log(1+1/x) → 1/2. The
        // candidate is accepted only after a numeric check against the original at
        // large x, so a one-sided/two-sided mismatch (or a wrong t-limit) cannot
        // slip through.
        // Skip the substitution when a general power f(x)^g(x) (var in both base
        // and exponent, e.g. x^(1/x)) sits inside a *sum* — the f(x)^g(x) − c
        // difference shape. x = 1/t only turns it into an equally hard t→0 form
        // (x^(1/x) ↦ exp(−t·log t)) on which the t→0 machinery can spin. A
        // standalone such power (e.g. the Stirling root (n!)^(1/n)) is *not*
        // skipped — the substitution resolves it cleanly.
        bool var_power_in_sum = false;
        {
            auto vp = [&](auto&& self, const Expr& e, bool under_add) -> void {
                if (var_power_in_sum) return;
                if (under_add && e->type_id() == TypeId::Pow
                    && has(e->args()[0], var) && has(e->args()[1], var)) {
                    var_power_in_sum = true;
                    return;
                }
                const bool child_under_add =
                    under_add || e->type_id() == TypeId::Add;
                for (const auto& a : e->args()) self(self, a, child_under_add);
            };
            vp(vp, expr, false);
        }
        if (is_infinity(target) && depth < 8 && has(expr, var)
            && !var_power_in_sum) {
            const Expr sgn = target->type_id() == TypeId::Infinity
                                 ? S::One()
                                 : S::NegativeOne();
            Expr tt = symbol("__lim_recip_t");
            Expr sub = subs(expr, var, mul(sgn, pow(tt, integer(-1))));
            Expr cand = limit_impl(sub, tt, S::Zero(), depth + 1);
            if (!is_nan(cand) && !is_infinity(cand) && !has(cand, tt)
                && cand->type_id() != TypeId::ComplexInfinity) {
                Expr cv = evalf(cand, 30);
                if (cv->type_id() == TypeId::Float) {
                    double cvd = 0.0;
                    bool parsed = false;
                    try {
                        cvd = std::stod(cv->str());
                        parsed = true;
                    } catch (...) {
                    }
                    // Sample at increasing magnitudes and require convergence to
                    // cand: the diff must shrink and the largest sample land very
                    // close. A wrong candidate keeps an O(1) offset and is rejected;
                    // the loose-to-tight ladder tolerates slow (∼1/x) convergence.
                    bool ok = parsed;
                    double prev = 1e300;
                    int checks = 0;
                    for (long xv : {1000L, 1000000L, 1000000000L}) {
                        if (!ok) break;
                        Expr at = evalf(
                            simplify(subs(expr, var, mul(sgn, integer(xv)))), 40);
                        if (at->type_id() != TypeId::Float) { ok = false; break; }
                        try {
                            const double d = std::fabs(std::stod(at->str()) - cvd);
                            if (!(d <= prev + 1e-12)) ok = false;  // must not diverge
                            prev = d;
                            ++checks;
                        } catch (...) {
                            ok = false;
                        }
                    }
                    if (ok && checks == 3 && prev < 1e-4) return cand;
                }
            }
        }
    }
    // Essential singularity at a finite point (exp(1/x), 1/x² at 0): substitute
    // u = 1/(x − a) and take u → ±∞; the two one-sided limits agree iff the
    // two-sided limit exists. Only at a finite target with a non-finite direct
    // value and a reciprocal singularity (a negative power of a factor that
    // vanishes at the target), so ordinary limits are untouched. Resolves
    // exp(−1/x²) → 0 and x/(exp(1/x)−1) → 0.
    if (depth < 12 && is_number(target)
        && (is_nan(direct) || direct->type_id() == TypeId::ComplexInfinity)) {
        bool has_sing = false;
        auto scan = [&](auto&& self, const Expr& e) -> void {
            if (has_sing) return;
            if (e->type_id() == TypeId::Pow
                && e->args()[1]->type_id() == TypeId::Integer
                && static_cast<const Integer&>(*e->args()[1]).is_negative()
                && has(e->args()[0], var)
                && limit_impl(e->args()[0], var, target, depth + 1)
                       == S::Zero()) {
                has_sing = true;
                return;
            }
            for (const auto& a : e->args()) self(self, a);
        };
        scan(scan, expr);
        if (has_sing) {
            Expr u = symbol("__recip_pt_u");
            Expr sub = subs(expr, var, add(target, pow(u, S::NegativeOne())));
            Expr lr = limit_impl(sub, u, S::Infinity(), depth + 1);
            Expr ll = limit_impl(sub, u, S::NegativeInfinity(), depth + 1);
            if (!is_nan(lr) && !is_nan(ll) && lr == ll) return lr;
        }
    }
    // An unresolved oscillation — direct substitution left sin/cos/tan of ∞
    // (e.g. 1 + sin(∞)). The limit does not exist as a single value; report nan
    // instead of the meaningless f(∞). Convergent oscillation cases are handled
    // by earlier rules and never reach here.
    if (has_oscillating_infinity(direct)) return S::NaN();
    return direct;
}

}  // namespace

Expr limit(const Expr& expr, const Expr& var, const Expr& target) {
    Expr r = limit_impl(expr, var, target, 0);
    // Final sanity: an ∞ buried inside arithmetic (e.g. ∞·(∞+∞·log 4)⁻¹ from a
    // divergent gamma ratio) is an unresolved indeterminate, not an answer.
    return has_buried_infinity(r) ? Expr{S::NaN()} : r;
}

Expr limit(const Expr& expr, const Expr& var, const Expr& target, int dir) {
    // dir 0 (two-sided) and limits at ±∞ (already one-sided by nature) use the
    // standard engine.
    if (dir == 0 || !is_number(target)) return limit_impl(expr, var, target, 0);
    // Reduce a one-sided finite limit to a limit at infinity via the
    // substitution x = target + 1/u, taking u → +∞ for the right limit and
    // u → −∞ for the left (1/u → 0 from the matching side). Using +1/u keeps
    // reciprocals un-nested (1/x ↦ u, not (±1/u)⁻¹), which the engine resolves
    // cleanly. u carries the sign of its target (positive for +∞, negative for
    // −∞) so simplify resolves the sign-dependent nodes: Abs(1/u) → ±1/u,
    // log(1/u) → −log(u). This makes the approach direction explicit and reuses
    // the well-tested ±∞ machinery, handling poles (1/x → ±∞), sign/abs
    // discontinuities (|x|/x → ±1) and essential singularities (exp(1/x) → ∞
    // from the right, 0 from the left).
    AssumptionMask um;
    um.set_positive(dir > 0).set_negative(dir < 0);
    Expr u = symbol("__onesided_u", um);
    Expr shifted = add(target, pow(u, S::NegativeOne()));
    Expr utarget = dir > 0 ? Expr{S::Infinity()} : Expr{S::NegativeInfinity()};
    Expr r = limit_impl(simplify(subs(expr, var, shifted)), u, utarget, 0);
    // If the transformed problem is intractable but the ordinary two-sided limit
    // is a determinate finite value (a continuous point), prefer that.
    if (is_nan(r) || r->type_id() == TypeId::ComplexInfinity) {
        Expr two = limit_impl(expr, var, target, 0);
        if (!is_nan(two) && two->type_id() != TypeId::ComplexInfinity) return two;
    }
    return has_buried_infinity(r) ? Expr{S::NaN()} : r;
}

}  // namespace sympp
