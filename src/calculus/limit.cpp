#include <sympp/calculus/limit.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/infinity.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
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
    std::vector<Expr> zero_f, inf_f, finite_f;
    for (const auto& f : expr->args()) {
        Expr lf = limit_impl(f, var, target, depth + 1);
        if (is_nan(lf)) return std::nullopt;
        if (is_infinity(lf)) inf_f.push_back(f);
        else if (lf == S::Zero()) zero_f.push_back(f);
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
    if (den == S::Zero() || has_var_radical(num, var)) return std::nullopt;
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
            // Leading terms cancel. For a two-term sum with a radical, the conjugate
            // t₁ + t₂ = (t₁² − t₂²)/(t₁ − t₂) clears it; recurse on that form (no
            // simplify on the quotient, which would rationalize straight back).
            if (e->args().size() == 2 && has_var_radical(e, var)) {
                const Expr& t1 = e->args()[0];
                const Expr& t2 = e->args()[1];
                Expr num = simplify(mul(t1, t1) + mul(S::NegativeOne(), mul(t2, t2)));
                Expr den = simplify(t1 + mul(S::NegativeOne(), t2));
                if (!(den == S::Zero()) && !has_var_radical(num, var)) {
                    return leading_pos_inf(
                        mul(num, pow(den, S::NegativeOne())), var, depth + 1);
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
    std::fprintf(stderr, "DBG alginf expr=%s lead=%s\n", expr->str().c_str(),
                 lead ? (lead->first->str() + " | " + lead->second->str()).c_str()
                      : "NULLOPT");
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

Expr limit_impl(const Expr& expr, const Expr& var, const Expr& target,
                int depth) {
    // Reciprocal trig/hyperbolic functions are opaque to the limit machinery;
    // rewrite them as sin/cos ratios (cot x → cos x/sin x, …) and retry, so
    // forms like x·cot(x) resolve via L'Hôpital instead of returning nan.
    if (depth < 12) {
        Expr rw = rewrite_reciprocal_trig(expr);
        if (!(rw == expr)) return limit_impl(rw, var, target, depth + 1);
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
            }
            // Logarithms: log(g) → log(lim g), and combine a ∞ − ∞ between logs.
            if (auto v = try_log_limit(expr, var, target, depth)) return *v;
            // Linearity over a sum: when every term has a determinate finite
            // limit, the limit is their sum. Direct substitution gives nan when
            // a single term is an ∞·0 product (e.g. the antiderivative
            // −x·e^(−x) − e^(−x) at +∞, whose terms each → 0). Bail if any term
            // diverges, so an ∞ − ∞ cancellation still falls through to L'Hôpital
            // on the combined fraction.
            if (expr->type_id() == TypeId::Add) {
                std::vector<Expr> term_limits;
                bool all_finite = true;
                for (const auto& t : expr->args()) {
                    Expr tl = limit_impl(t, var, target, depth + 1);
                    if (is_nan(tl) || is_infinity(tl)) { all_finite = false; break; }
                    term_limits.push_back(std::move(tl));
                }
                if (all_finite) return add(std::move(term_limits));
            }
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
            // Merge a product of constant-base exponentials (2^x/3^x,
            // exp(x)/exp(2x)) into one exp(rate) before the generic product path,
            // which would otherwise see ∞·0 and stall in L'Hôpital → nan.
            if (auto v = try_exponential_product(expr, var, target, depth)) {
                return *v;
            }
            if (auto v = try_product_form(expr, var, target, depth)) return *v;
        }
        // 0/0 and ∞/∞ quotients (also recovers finite 0/0 where direct
        // substitution collapses to 0 or nan).
        if (auto v = lhopital(expr, var, target)) return *v;
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
    return direct;
}

}  // namespace

Expr limit(const Expr& expr, const Expr& var, const Expr& target) {
    return limit_impl(expr, var, target, 0);
}

}  // namespace sympp
