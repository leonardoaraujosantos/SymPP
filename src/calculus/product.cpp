#include <sympp/calculus/product.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

namespace {

// Unevaluated product marker — an UndefinedFunction `Product(expr, var, lo, hi)`,
// mirroring summation's `Sum(...)`. Returned when no closed form is found, so a
// product is never silently dropped to its bare term.
[[nodiscard]] Expr product_marker(const Expr& expr, const Expr& var,
                                  const Expr& lo, const Expr& hi) {
    return function_symbol("Product")(std::vector<Expr>{expr, var, lo, hi});
}

[[nodiscard]] bool is_unevaluated(const Expr& e) {
    return e->str().find("Product(") != std::string::npos;
}

// ∏_{k=lo}^{hi} P(k) for a polynomial P in var with closed-form roots. Writing
// P = lead·∏ᵢ(k − rᵢ) and using ∏_{k=lo}^{hi}(k − r) = Γ(hi − r + 1)/Γ(lo − r),
// the whole product is lead^N · ∏ᵢ Γ(hi − rᵢ + 1)/Γ(lo − rᵢ), N = hi − lo + 1.
// Covers ∏ k = Γ(hi+1)/Γ(lo), ∏ (k+a), ∏ (k²−1) = ∏(k−1)(k+1), … Returns
// nullopt unless expr is a polynomial of degree ≥ 1 with every root in closed
// form. The result is exact (equivalent to SymPy's Gamma/RisingFactorial form).
[[nodiscard]] std::optional<Expr> prod_polynomial(const Expr& expr,
                                                  const Expr& var,
                                                  const Expr& lo,
                                                  const Expr& hi) {
    std::size_t deg = 0;
    std::vector<Expr> rts;
    Expr lead;
    try {
        Poly p(expand(expr), var);
        deg = p.degree();
        if (deg < 1) return std::nullopt;
        rts = p.roots();
        lead = p.leading_coeff();
    } catch (const std::exception&) {
        return std::nullopt;
    }
    if (rts.size() != deg) return std::nullopt;  // need every root in closed form

    Expr count = simplify(hi - lo + integer(1));
    Expr result = pow(lead, count);
    for (const auto& r : rts) {
        Expr ratio = mul(gamma(simplify(hi - r + integer(1))),
                         pow(gamma(simplify(lo - r)), integer(-1)));
        result = mul(result, ratio);
    }
    return simplify(result);
}

}  // namespace

Expr product(const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (!expr) return expr;

    // Single-term range (hi == lo): ∏_{k=a}^{a} f(k) = f(a).
    if (hi == lo) return simplify(subs(expr, var, lo));

    // Constant w.r.t. var: c^(hi − lo + 1).
    if (!has(expr, var)) {
        return simplify(pow(expr, simplify(hi - lo + integer(1))));
    }

    // Multiplicative split: ∏ (f·g) = (∏ f)·(∏ g); constant factors → c^N.
    if (expr->type_id() == TypeId::Mul) {
        Expr count = simplify(hi - lo + integer(1));
        std::vector<Expr> pieces;
        std::vector<Expr> const_factors;
        for (const auto& f : expr->args()) {
            if (!has(f, var)) {
                const_factors.push_back(f);
                continue;
            }
            Expr pf = product(f, var, lo, hi);
            if (is_unevaluated(pf)) return product_marker(expr, var, lo, hi);
            pieces.push_back(std::move(pf));
        }
        if (!const_factors.empty()) {
            pieces.push_back(pow(mul(std::move(const_factors)), count));
        }
        return simplify(mul(std::move(pieces)));
    }

    // Powers.
    if (expr->type_id() == TypeId::Pow) {
        const Expr& base = expr->args()[0];
        const Expr& e = expr->args()[1];
        // Constant exponent: ∏ base(k)^c = (∏ base(k))^c.
        if (!has(e, var)) {
            Expr pb = product(base, var, lo, hi);
            if (is_unevaluated(pb)) return product_marker(expr, var, lo, hi);
            return simplify(pow(pb, e));
        }
        // Geometric exponent: ∏ r^(c·k + d) = r^(c·Σk + d·N) for var-free base r.
        if (!has(base, var)) {
            Expr c = diff(e, var);
            if (!has(c, var)) {
                Expr d = simplify(e - mul(c, var));
                if (!has(d, var)) {
                    Expr count = simplify(hi - lo + integer(1));
                    Expr sumk = simplify(
                        mul(count, mul(lo + hi, pow(integer(2), integer(-1)))));
                    return simplify(
                        pow(base, simplify(mul(c, sumk) + mul(d, count))));
                }
            }
        }
        return product_marker(expr, var, lo, hi);
    }

    // Canonical ∏_{k=1}^{n} k = n! (cleaner than the Γ(n+1)/Γ(1) ratio).
    if (expr == var && lo == S::One()) return factorial(hi);

    // Polynomial in var → Gamma-ratio product over its linear factors.
    if (auto r = prod_polynomial(expr, var, lo, hi)) return *r;

    return product_marker(expr, var, lo, hi);
}

}  // namespace sympp
