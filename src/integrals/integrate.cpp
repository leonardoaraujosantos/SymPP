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
#include <sympp/core/expand.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/polys/poly.hpp>

namespace sympp {

namespace {

// Helper: does `e` depend on `var`? Uses has() from traversal.
[[nodiscard]] bool depends_on(const Expr& e, const Expr& var) noexcept {
    return has(e, var);
}

// Decompose `e` as (a, b) such that e == a*var + b, with a, b independent
// of var. Returns nullopt if e isn't affine in var.
[[nodiscard]] std::optional<std::pair<Expr, Expr>>
as_affine(const Expr& e, const Expr& var) {
    if (!e) return std::nullopt;
    if (!has(e, var)) return std::pair{S::Zero(), e};
    if (e == var) return std::pair{S::One(), S::Zero()};
    if (e->type_id() == TypeId::Mul) {
        // Split into constant and var-dependent factors. Any factor of
        // var must itself be affine, and exactly one such factor is
        // permitted (multiplying two affines in var produces a quadratic).
        std::vector<Expr> const_factors;
        std::vector<Expr> var_factors;
        for (const auto& a : e->args()) {
            if (has(a, var)) var_factors.push_back(a);
            else const_factors.push_back(a);
        }
        if (var_factors.size() != 1) return std::nullopt;
        auto sub = as_affine(var_factors[0], var);
        if (!sub) return std::nullopt;
        Expr scale = mul(const_factors);
        return std::pair{mul(scale, sub->first), mul(scale, sub->second)};
    }
    if (e->type_id() == TypeId::Add) {
        Expr a_acc = S::Zero();
        Expr b_acc = S::Zero();
        for (const auto& term : e->args()) {
            auto sub = as_affine(term, var);
            if (!sub) return std::nullopt;
            a_acc = add(a_acc, sub->first);
            b_acc = add(b_acc, sub->second);
        }
        return std::pair{a_acc, b_acc};
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> integrate_term(const Expr& term, const Expr& var);

// Forward-decl of the public entry — integrate_term needs to recurse
// through the full dispatch (Add linearity, trig, rational) when it
// pulls a constant out of a Mul whose remainder is an Add.
}  // namespace
Expr integrate(const Expr& expr, const Expr& var);
namespace {

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

    // (a*x + b)^n with n constant: ∫(ax+b)^n dx = (ax+b)^(n+1) / (a*(n+1))
    // and n == -1 → log(ax+b)/a.
    if (term->type_id() == TypeId::Pow) {
        const auto& base = term->args()[0];
        const auto& exp = term->args()[1];
        if (!depends_on(exp, var)) {
            auto aff = as_affine(base, var);
            if (aff && !(aff->first == S::Zero())) {
                if (exp == S::NegativeOne()) {
                    return log(base) / aff->first;
                }
                Expr n_plus_1 = add(exp, S::One());
                return pow(base, n_plus_1) / (aff->first * n_plus_1);
            }
        }
    }

    // 1/x: ∫ 1/x dx = log(x). 1/x is Pow(x, -1) — caught above.

    // sin(ax+b), cos(ax+b), exp(ax+b) where inner is affine in var.
    if (term->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*term);
        if (fn.args().size() == 1) {
            const auto& inner = fn.args()[0];
            auto aff = as_affine(inner, var);
            if (aff && !(aff->first == S::Zero())) {
                const Expr& a = aff->first;
                switch (fn.function_id()) {
                    case FunctionId::Sin:
                        // ∫sin(ax+b) dx = -cos(ax+b)/a
                        return mul(S::NegativeOne(), cos(inner)) / a;
                    case FunctionId::Cos:
                        return sin(inner) / a;
                    case FunctionId::Exp:
                        return exp(inner) / a;
                    default:
                        break;
                }
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
        // If we didn't pull anything constant out, recursing on
        // mul(var_factors) would just be the same term — bail to avoid
        // infinite recursion. Higher-level strategies (try_rational,
        // integration by parts, etc.) get a chance instead.
        if (constant_factors.empty()) {
            return std::nullopt;
        }
        // Try integrating the remaining var-dependent product. Use the
        // full integrate() dispatch (not just integrate_term) so an Add
        // inner gets linearity, and so trig / rational strategies still
        // apply.
        Expr inner = mul(std::move(var_factors));
        Expr sub = integrate(inner, var);
        // integrate() returns an Integral(_, _) marker on failure (an
        // UndefinedFunction named "Integral"); propagate as nullopt so
        // higher strategies can try other approaches.
        if (sub->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*sub);
            if (fn.name() == "Integral") return std::nullopt;
        }
        return mul(mul(std::move(constant_factors)), sub);
    }

    // Outside the table.
    return std::nullopt;
}

// Trig identity rewrites:
//   sin²(u) → (1 - cos(2u))/2
//   cos²(u) → (1 + cos(2u))/2
//   sin(p)cos(q) → (sin(p+q) + sin(p-q))/2
// After rewriting, recurse via integrate so the existing affine-argument
// rules pick up the resulting sin(linear) / cos(linear) terms.
[[nodiscard]] std::optional<Expr> try_trig_reduction(const Expr& expr, const Expr& var);

// Try to integrate `expr` as a rational function in var. Uses
// polynomial division for the improper part and apart() for the proper
// part. Returns nullopt when:
//   * expr has no top-level denominator structure, or
//   * num / den isn't pure-polynomial in var, or
//   * apart can't decompose the proper remainder (e.g. repeated roots
//     beyond what apart handles in this build).
[[nodiscard]] std::optional<Expr> try_rational(const Expr& expr, const Expr& var);

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
    if (auto r = try_trig_reduction(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_rational(expr, var); r.has_value()) {
        return *r;
    }

    // Outside the closed-form table — return an unevaluated marker. Use an
    // UndefinedFunction named "Integral" so callers can detect it.
    return function_symbol("Integral")(expr, var);
}

namespace {

std::optional<Expr> try_trig_reduction(const Expr& expr, const Expr& var) {
    // sin²(u) → (1 - cos(2u))/2 ; cos²(u) → (1 + cos(2u))/2
    if (expr->type_id() == TypeId::Pow
        && expr->args()[1] == integer(2)) {
        const Expr& base = expr->args()[0];
        if (base->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*base);
            if (fn.args().size() == 1 && depends_on(fn.args()[0], var)) {
                const Expr& u = fn.args()[0];
                if (fn.function_id() == FunctionId::Sin) {
                    Expr rewritten = (integer(1) - cos(integer(2) * u))
                                     / integer(2);
                    return integrate(rewritten, var);
                }
                if (fn.function_id() == FunctionId::Cos) {
                    Expr rewritten = (integer(1) + cos(integer(2) * u))
                                     / integer(2);
                    return integrate(rewritten, var);
                }
            }
        }
    }

    // sin(p)cos(q) → (sin(p+q) + sin(p-q))/2 — product-to-sum.
    if (expr->type_id() == TypeId::Mul) {
        Expr sin_arg, cos_arg;
        std::vector<Expr> rest;
        for (const auto& f : expr->args()) {
            if (!sin_arg && f->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*f);
                if (fn.function_id() == FunctionId::Sin
                    && fn.args().size() == 1) {
                    sin_arg = f->args()[0];
                    continue;
                }
            }
            if (!cos_arg && f->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*f);
                if (fn.function_id() == FunctionId::Cos
                    && fn.args().size() == 1) {
                    cos_arg = f->args()[0];
                    continue;
                }
            }
            rest.push_back(f);
        }
        if (sin_arg && cos_arg
            && depends_on(sin_arg, var) && depends_on(cos_arg, var)) {
            Expr rewritten = (sin(sin_arg + cos_arg)
                              + sin(sin_arg - cos_arg))
                             / integer(2);
            for (const auto& r : rest) rewritten = rewritten * r;
            return integrate(rewritten, var);
        }
    }
    return std::nullopt;
}

std::optional<Expr> try_rational(const Expr& expr, const Expr& var) {
    // Bring to single fraction.
    Expr t = together(expr);
    Expr num = S::One();
    Expr den;
    if (t->type_id() == TypeId::Pow
        && t->args()[1] == S::NegativeOne()) {
        den = t->args()[0];
    } else if (t->type_id() == TypeId::Mul) {
        std::vector<Expr> num_factors;
        for (const auto& f : t->args()) {
            if (f->type_id() == TypeId::Pow
                && f->args()[1] == S::NegativeOne()) {
                if (den) return std::nullopt;
                den = f->args()[0];
            } else {
                num_factors.push_back(f);
            }
        }
        if (!den) return std::nullopt;
        num = mul(num_factors);
    } else {
        return std::nullopt;
    }

    // Verify num and den are pure polynomials in var (constant
    // coefficients only — reject things like sin(x)/cos(x)).
    Poly num_p(expand(num), var);
    Poly den_p(expand(den), var);
    for (const auto& c : num_p.coeffs()) if (has(c, var)) return std::nullopt;
    for (const auto& c : den_p.coeffs()) if (has(c, var)) return std::nullopt;
    if (den_p.degree() == 0) return std::nullopt;

    // Polynomial division: num = q * den + r, deg(r) < deg(den).
    auto [q, r_poly] = num_p.divmod(den_p);

    if (r_poly.is_zero()) {
        // Pure polynomial result — let the linearity path handle it.
        return integrate(q.as_expr(), var);
    }

    // Decompose the proper remainder via apart().
    Expr proper = r_poly.as_expr() / den_p.as_expr();
    Expr apart_form = apart(proper, var);
    if (apart_form == proper) {
        // apart couldn't decompose further — fall through.
        return std::nullopt;
    }

    Expr poly_int = integrate(q.as_expr(), var);
    Expr apart_int = integrate(apart_form, var);
    return poly_int + apart_int;
}

}  // namespace

Expr integrate(const Expr& expr, const Expr& var,
              const Expr& lower, const Expr& upper) {
    auto antider = integrate(expr, var);
    return subs(antider, var, upper) - subs(antider, var, lower);
}

}  // namespace sympp
