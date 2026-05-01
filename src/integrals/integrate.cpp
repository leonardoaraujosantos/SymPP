#include <sympp/integrals/integrate.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>

namespace sympp {

namespace {

// Helper: does `e` depend on `var`? Uses has() from traversal.
[[nodiscard]] bool depends_on(const Expr& e, const Expr& var) noexcept {
    return has(e, var);
}

// Try the elementary table on a single term. Returns std::nullopt if the
// term is outside the table.
[[nodiscard]] std::optional<Expr> integrate_term(const Expr& term, const Expr& var) {
    // Constant w.r.t. var: ∫ c dx = c*x
    if (!depends_on(term, var)) {
        return mul(term, var);
    }

    // var itself: ∫ x dx = x^2/2
    if (term == var) {
        return mul(rational(1, 2), pow(var, integer(2)));
    }

    // x^n where n is a non-(-1) integer or symbolic constant: x^(n+1)/(n+1)
    if (term->type_id() == TypeId::Pow) {
        const auto& base = term->args()[0];
        const auto& exp = term->args()[1];
        if (base == var && !depends_on(exp, var)) {
            // Skip the n = -1 case → log(x).
            if (exp == S::NegativeOne()) {
                return log(var);
            }
            auto new_exp = add(exp, S::One());
            return mul(pow(var, new_exp), pow(new_exp, S::NegativeOne()));
        }
    }

    // 1/x: ∫ 1/x dx = log(x). 1/x is Pow(x, -1) — caught above.

    // sin(c*x), cos(c*x), exp(c*x) where the inner is c*var with c constant.
    if (term->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*term);
        if (fn.args().size() == 1) {
            const auto& inner = fn.args()[0];
            // Helper lambda: given inner = c*var (or var), extract c.
            auto extract_const_coeff = [&](const Expr& g) -> std::optional<Expr> {
                if (g == var) return S::One();
                if (g->type_id() == TypeId::Mul) {
                    std::vector<Expr> rest;
                    bool found_var = false;
                    for (const auto& a : g->args()) {
                        if (a == var) {
                            if (found_var) return std::nullopt;  // var^2 etc.
                            found_var = true;
                        } else if (depends_on(a, var)) {
                            return std::nullopt;
                        } else {
                            rest.push_back(a);
                        }
                    }
                    if (!found_var) return std::nullopt;
                    return mul(std::move(rest));
                }
                return std::nullopt;
            };

            switch (fn.function_id()) {
                case FunctionId::Sin: {
                    auto c = extract_const_coeff(inner);
                    if (!c.has_value()) break;
                    // ∫ sin(c*x) dx = -cos(c*x)/c
                    return mul({S::NegativeOne(), cos(inner), pow(*c, S::NegativeOne())});
                }
                case FunctionId::Cos: {
                    auto c = extract_const_coeff(inner);
                    if (!c.has_value()) break;
                    return mul(sin(inner), pow(*c, S::NegativeOne()));
                }
                case FunctionId::Exp: {
                    auto c = extract_const_coeff(inner);
                    if (!c.has_value()) break;
                    return mul(exp(inner), pow(*c, S::NegativeOne()));
                }
                default:
                    break;
            }
        }
    }

    // Mul: pull out factors that don't depend on var, then recurse on the rest.
    if (term->type_id() == TypeId::Mul) {
        std::vector<Expr> constant_factors;
        std::vector<Expr> var_factors;
        for (const auto& a : term->args()) {
            if (depends_on(a, var)) {
                var_factors.push_back(a);
            } else {
                constant_factors.push_back(a);
            }
        }
        if (var_factors.empty()) {
            // Entire term is constant — handled at the top of this function;
            // shouldn't reach here.
            return mul(mul(std::move(constant_factors)), var);
        }
        // Try integrating the remaining var-dependent product.
        Expr inner = mul(std::move(var_factors));
        if (auto sub = integrate_term(inner, var); sub.has_value()) {
            return mul(mul(std::move(constant_factors)), *sub);
        }
        return std::nullopt;
    }

    // Outside the table.
    return std::nullopt;
}

}  // namespace

Expr integrate(const Expr& expr, const Expr& var) {
    if (!expr || !var) return S::Zero();

    // Linearity: split Add into per-term integrals.
    if (expr->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(expr->args().size());
        for (const auto& a : expr->args()) {
            terms.push_back(integrate(a, var));
        }
        return add(std::move(terms));
    }

    if (auto r = integrate_term(expr, var); r.has_value()) {
        return *r;
    }

    // Outside the closed-form table — return an unevaluated marker. Use an
    // UndefinedFunction named "Integral" so callers can detect it.
    return function_symbol("Integral")(expr, var);
}

Expr integrate(const Expr& expr, const Expr& var,
              const Expr& lower, const Expr& upper) {
    auto antider = integrate(expr, var);
    return subs(antider, var, upper) - subs(antider, var, lower);
}

}  // namespace sympp
