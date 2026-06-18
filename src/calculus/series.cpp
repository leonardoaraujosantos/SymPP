#include <sympp/calculus/series.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/calculus/limit.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/expr_collections.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

namespace {
// A non-finite Taylor coefficient — the function is singular at the expansion
// point (e.g. log(x), 1/x at 0), or direct substitution hit an indeterminate
// 0/0 form (a removable singularity such as sin(x)/x).
[[nodiscard]] bool is_non_finite(const Expr& e) {
    switch (e->type_id()) {
        case TypeId::ComplexInfinity:
        case TypeId::Infinity:
        case TypeId::NegativeInfinity:
        case TypeId::NaN:
            return true;
        default:
            return false;
    }
}

// Ordinary Taylor expansion of an analytic function: returns the degree-(n−1)
// polynomial Σ f⁽ᵏ⁾(x0)/k!·(x−x0)ᵏ, or nullopt if a coefficient is non-finite
// (a pole or an unresolved indeterminate), in which case the caller may try a
// Laurent expansion instead.
[[nodiscard]] std::optional<Expr> taylor_series(const Expr& expr,
                                                const Expr& var, const Expr& x0,
                                                std::size_t n) {
    std::vector<Expr> terms;
    terms.reserve(n);
    Expr current_deriv = expr;
    Expr dx = var - x0;
    for (std::size_t k = 0; k < n; ++k) {
        Expr value_at_x0 = subs(current_deriv, var, x0);
        if (is_non_finite(value_at_x0)) {
            // A removable singularity (0/0 like sin(x)/x) has a finite limit;
            // a genuine singularity (log, 1/x) does not.
            value_at_x0 = limit(current_deriv, var, x0);
            if (is_non_finite(value_at_x0)) return std::nullopt;
        }
        Expr fact = factorial(integer(static_cast<long>(k)));
        terms.push_back(mul({pow(fact, S::NegativeOne()), value_at_x0,
                             pow(dx, integer(static_cast<long>(k)))}));
        if (k + 1 < n) current_deriv = diff(current_deriv, var);
    }
    return add(std::move(terms));
}

// Rewrite reciprocal trig/hyperbolic functions as sin/cos ratios, so a pole
// surfaces as an explicit vanishing denominator the Laurent path can divide
// through: cot→cos/sin, csc→1/sin, sec→1/cos, coth→cosh/sinh, csch→1/sinh,
// sech→1/cosh.
[[nodiscard]] Expr rewrite_reciprocal_trig(const Expr& e) {
    ExprMap<Expr> m;
    auto scan = [&](auto&& self, const Expr& x) -> void {
        if (x->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*x);
            if (fn.args().size() == 1) {
                const Expr& u = fn.args()[0];
                switch (fn.function_id()) {
                    case FunctionId::Cot:
                        m.emplace(x, mul(cos(u), pow(sin(u), S::NegativeOne())));
                        break;
                    case FunctionId::Csc:
                        m.emplace(x, pow(sin(u), S::NegativeOne()));
                        break;
                    case FunctionId::Sec:
                        m.emplace(x, pow(cos(u), S::NegativeOne()));
                        break;
                    case FunctionId::Coth:
                        m.emplace(x,
                                  mul(cosh(u), pow(sinh(u), S::NegativeOne())));
                        break;
                    case FunctionId::Csch:
                        m.emplace(x, pow(sinh(u), S::NegativeOne()));
                        break;
                    case FunctionId::Sech:
                        m.emplace(x, pow(cosh(u), S::NegativeOne()));
                        break;
                    default:
                        break;
                }
            }
        }
        for (const auto& a : x->args()) self(self, a);
    };
    scan(scan, e);
    if (m.empty()) return e;
    return xreplace(e, m);
}

[[nodiscard]] long valuation(const std::vector<Expr>& v) {
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (!(simplify(v[i]) == S::Zero())) return static_cast<long>(i);
    }
    return -1;  // identically zero
}

// Laurent series at 0 of a ratio N(x)/D(x) with a denominator that vanishes at
// 0. Expands N and D as Taylor polynomials, then divides the power series:
// f = x^(v_N − v_D) · (Ñ/D̃) with Ñ(0), D̃(0) ≠ 0. Returns nullopt when the
// expression is not such a ratio.
[[nodiscard]] std::optional<Expr> try_laurent_series(const Expr& expr,
                                                     const Expr& var,
                                                     const Expr& x0,
                                                     std::size_t n,
                                                     bool only_if_pole = false) {
    if (!(x0 == S::Zero())) return std::nullopt;  // implemented at 0 only
    Expr e = rewrite_reciprocal_trig(expr);

    // Split into numerator / denominator (factors with a negative integer power
    // contribute to the denominator).
    Expr num = S::One();
    Expr den = S::One();
    auto add_factor = [&](const Expr& f) {
        if (f->type_id() == TypeId::Pow
            && f->args()[1]->type_id() == TypeId::Integer) {
            const auto& z = static_cast<const Integer&>(*f->args()[1]);
            if (z.is_negative() && z.fits_long()) {
                den = mul(den, pow(f->args()[0], integer(-z.to_long())));
                return;
            }
        }
        num = mul(num, f);
    };
    if (e->type_id() == TypeId::Mul) {
        for (const auto& f : e->args()) add_factor(f);
    } else {
        add_factor(e);
    }
    if (!has(den, var)) return std::nullopt;  // no var-dependent pole

    const long W = static_cast<long>(n) + 8;  // working order with headroom
    auto ns = taylor_series(num, var, x0, static_cast<std::size_t>(W));
    auto ds = taylor_series(den, var, x0, static_cast<std::size_t>(W));
    if (!ns || !ds) return std::nullopt;

    std::vector<Expr> a;
    std::vector<Expr> b;
    try {
        a = Poly(expand(*ns), var).coeffs();
        b = Poly(expand(*ds), var).coeffs();
    } catch (const std::exception&) {
        return std::nullopt;  // numerator/denominator not polynomial in var
    }
    auto pad = [&](std::vector<Expr>& v) {
        while (static_cast<long>(v.size()) < W) v.push_back(S::Zero());
    };
    pad(a);
    pad(b);

    const long vN = valuation(a);
    const long vD = valuation(b);
    if (vD < 0) return std::nullopt;            // denominator ≡ 0
    // When asked only to pre-empt poles, defer an analytic ratio (denominator
    // non-vanishing at 0) to the Taylor path, which keeps its existing form.
    if (only_if_pole && vD == 0) return std::nullopt;
    if (vN < 0) return Expr{S::Zero()};         // numerator ≡ 0
    const long lead = vN - vD;                  // leading exponent of result
    const long L = static_cast<long>(n) - lead;  // terms for exponents lead..n−1
    if (L <= 0) return std::nullopt;

    // Power-series division (Ñ/D̃), Ñ = a[vN..], D̃ = b[vD..], D̃(0) = b[vD] ≠ 0.
    const Expr& b0 = b[static_cast<std::size_t>(vD)];
    std::vector<Expr> q(static_cast<std::size_t>(L), S::Zero());
    for (long i = 0; i < L; ++i) {
        Expr s = (vN + i < static_cast<long>(a.size()))
                     ? a[static_cast<std::size_t>(vN + i)]
                     : Expr{S::Zero()};
        for (long j = 1; j <= i; ++j) {
            const long bidx = vD + j;
            if (bidx < static_cast<long>(b.size())) {
                s = add(s, mul(S::NegativeOne(),
                               mul(b[static_cast<std::size_t>(bidx)],
                                   q[static_cast<std::size_t>(i - j)])));
            }
        }
        q[static_cast<std::size_t>(i)] =
            simplify(mul(s, pow(b0, S::NegativeOne())));
    }

    std::vector<Expr> terms;
    for (long i = 0; i < L; ++i) {
        const Expr& qi = q[static_cast<std::size_t>(i)];
        if (qi == S::Zero()) continue;
        terms.push_back(mul(qi, pow(var, integer(lead + i))));
    }
    if (terms.empty()) return Expr{S::Zero()};
    return add(std::move(terms));
}

// Series of a composite f(g(x)) by composition, for an outer single-argument
// function f and an inner g that is analytic at x₀ with a finite value g(x₀)=c.
// taylor_series differentiates f(g) directly and evaluates each derivative via a
// limit; for a g built from a removable form (e.g. sin(x)/x, with its 1/x factor)
// those derivative-limits get hard and return nan, so the Taylor path bails on
// log(sin x / x). Composition sidesteps that: expand g to a series, expand f about
// the constant c (where f is a single function the Taylor path handles), and
// substitute (t−c) → (g_series − c), which has positive valuation so only finitely
// many terms reach a given order.
[[nodiscard]] std::optional<Expr> try_composition_series(const Expr& expr,
                                                         const Expr& var,
                                                         const Expr& x0,
                                                         std::size_t n) {
    // The outer operation is either a single-argument function f, or a power
    // g^p with a var-free exponent (the function t ↦ t^p, e.g. √(tan x / x)).
    Expr g;
    if (expr->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*expr);
        if (fn.args().size() != 1) return std::nullopt;
        g = fn.args()[0];
    } else if (expr->type_id() == TypeId::Pow) {
        if (has(expr->args()[1], var)) return std::nullopt;  // exponent must be const
        g = expr->args()[0];
    } else {
        return std::nullopt;
    }
    if (!has(g, var)) return std::nullopt;

    // Inner series via the full engine, so a removable pole that needs the
    // Laurent-division path (e.g. tan(x)/x) still expands. The result must be
    // analytic at x₀: a polynomial in (x − x₀) with a finite value c = g(x₀).
    Expr gs = series(g, var, x0, n);
    // g must be analytic at x₀: c = g(x₀) finite. Evaluating the *expanded* series
    // distinguishes a removable form (tan x / x → 1) from a genuine pole (csc x,
    // whose Laurent series → ∞ here), where composition must not apply.
    Expr c = subs(gs, var, x0);
    if (is_non_finite(c)) return std::nullopt;
    try {
        // Outer series of f(t) about t = c — f is a single function, the case the
        // Taylor path handles. (For log(t) at c = 0 this returns nullopt, so a
        // genuine singularity such as log(x) stays unexpanded.)
        Expr t = symbol("__series_compose_t");
        Expr ft = subs(expr, g, t);
        if (has(ft, var)) return std::nullopt;  // g didn't factor out cleanly
        auto outer = taylor_series(ft, t, c, n);
        if (!outer) return std::nullopt;

        // Compose: (t − c) ← (g_series − c), then drop terms of order ≥ n.
        Expr u = expand(gs - c);
        Expr composed = expand(subs(*outer, t, add(c, u)));
        Expr dx = symbol("__series_compose_dx");
        Expr shifted = expand(subs(composed, var, add(x0, dx)));
        Poly p(shifted, dx);
        std::vector<Expr> terms;
        const long deg = static_cast<long>(p.degree());
        for (long i = 0; i <= deg && i < static_cast<long>(n); ++i) {
            const Expr& ci = p.coeffs()[static_cast<std::size_t>(i)];
            if (ci == S::Zero()) continue;
            terms.push_back(mul(ci, pow(var - x0, integer(i))));
        }
        if (terms.empty()) return Expr{S::Zero()};
        return add(std::move(terms));
    } catch (const std::exception&) {
        return std::nullopt;  // inner/composed not an analytic polynomial in var
    }
}

}  // namespace

Expr series(const Expr& expr, const Expr& var, const Expr& x0, std::size_t n) {
    if (n == 0) return S::Zero();
    // A ratio whose denominator vanishes at x0 (removable singularity or pole)
    // is expanded by Laurent division first: the Taylor path would otherwise
    // compute n derivative coefficients each via a hard 0/0 limit and can hang
    // (e.g. x/(exp(x)−1), x/sin(x), x²/(1−cos x)). Laurent expands numerator and
    // denominator separately, so no division-induced indeterminate ever arises.
    if (auto l = try_laurent_series(expr, var, x0, n, /*only_if_pole=*/true)) {
        return *l;
    }
    // Analytic functions expand as an ordinary Taylor polynomial.
    if (auto t = taylor_series(expr, var, x0, n)) return *t;
    // A composite f(g) whose direct Taylor derivatives hit hard limits — e.g.
    // log(sin x / x) — expands by composing the inner and outer series.
    if (auto c = try_composition_series(expr, var, x0, n)) return *c;
    // A pole (cot, csc, 1/sin, …) yields a Laurent series via power-series
    // division; genuine singularities (log x) stay unexpanded, matching SymPy.
    if (auto l = try_laurent_series(expr, var, x0, n)) return *l;
    return expr;
}

}  // namespace sympp
