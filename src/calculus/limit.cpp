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
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/exponential.hpp>
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
            // Re-rationalise the new ratio with together() (cheap) before the
            // next step; full simplify() here is too slow over many iterations.
            auto nd = split_after_together(num_d / den_d);
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
    Expr inner = limit_impl(mul(ex, log(base)), var, target, depth + 1);
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

Expr limit_impl(const Expr& expr, const Expr& var, const Expr& target,
                int depth) {
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
            }
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
            if (auto v = try_product_form(expr, var, target, depth)) return *v;
        }
        // 0/0 and ∞/∞ quotients (also recovers finite 0/0 where direct
        // substitution collapses to 0 or nan).
        if (auto v = lhopital(expr, var, target)) return *v;
    }
    return direct;
}

}  // namespace

Expr limit(const Expr& expr, const Expr& var, const Expr& target) {
    return limit_impl(expr, var, target, 0);
}

}  // namespace sympp
