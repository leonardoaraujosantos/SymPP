#include <sympp/integrals/transforms.hpp>

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
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/functions/combinatorial.hpp>

namespace sympp {

namespace {

// Try to extract a constant `a` such that `g == a*t`. Returns std::nullopt
// when g doesn't have that form.
[[nodiscard]] std::optional<Expr> extract_linear_coeff(const Expr& g, const Expr& t) {
    if (g == t) return S::One();
    if (g->type_id() == TypeId::Mul) {
        std::vector<Expr> rest;
        bool found_t = false;
        for (const auto& a : g->args()) {
            if (a == t) {
                if (found_t) return std::nullopt;
                found_t = true;
            } else if (has(a, t)) {
                return std::nullopt;
            } else {
                rest.push_back(a);
            }
        }
        if (!found_t) return std::nullopt;
        return mul(std::move(rest));
    }
    return std::nullopt;
}

// Table-lookup on a single term. Returns std::nullopt if no entry matches.
[[nodiscard]] std::optional<Expr> laplace_term(const Expr& f, const Expr& t, const Expr& s) {
    // Constant in t: ∫ c * exp(-s*t) dt from 0 to ∞ = c / s
    if (!has(f, t)) {
        return mul(f, pow(s, S::NegativeOne()));
    }

    // f == t: 1/s^2
    if (f == t) {
        return pow(s, integer(-2));
    }

    // f == t^n: n! / s^(n+1)
    if (f->type_id() == TypeId::Pow) {
        const auto& base = f->args()[0];
        const auto& exp = f->args()[1];
        if (base == t && exp->type_id() == TypeId::Integer) {
            const auto& n = static_cast<const Integer&>(*exp);
            if (!n.is_negative() && n.fits_long() && n.to_long() <= 100) {
                long ne = n.to_long();
                Expr fact = factorial(integer(ne));
                Expr sn1 = pow(s, integer(ne + 1));
                return mul(fact, pow(sn1, S::NegativeOne()));
            }
        }
    }

    // exp / sin / cos with linear inner.
    if (f->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.args().size() == 1) {
            const auto& inner = fn.args()[0];
            auto coeff = extract_linear_coeff(inner, t);
            if (coeff.has_value()) {
                switch (fn.function_id()) {
                    case FunctionId::Exp:
                        // L{exp(a*t)} = 1/(s - a)
                        return pow(s - *coeff, S::NegativeOne());
                    case FunctionId::Sin: {
                        // L{sin(a*t)} = a/(s^2 + a^2)
                        auto denom = pow(s, integer(2)) + pow(*coeff, integer(2));
                        return mul(*coeff, pow(denom, S::NegativeOne()));
                    }
                    case FunctionId::Cos: {
                        // L{cos(a*t)} = s/(s^2 + a^2)
                        auto denom = pow(s, integer(2)) + pow(*coeff, integer(2));
                        return mul(s, pow(denom, S::NegativeOne()));
                    }
                    default:
                        break;
                }
            }
        }
    }

    // Mul with a constant factor: pull it out and recurse.
    if (f->type_id() == TypeId::Mul) {
        std::vector<Expr> consts;
        std::vector<Expr> rest;
        for (const auto& a : f->args()) {
            if (has(a, t)) rest.push_back(a);
            else consts.push_back(a);
        }
        if (!consts.empty() && !rest.empty()) {
            Expr inner = mul(std::move(rest));
            if (auto sub = laplace_term(inner, t, s); sub.has_value()) {
                return mul(mul(std::move(consts)), *sub);
            }
        }
    }

    return std::nullopt;
}

}  // namespace

Expr laplace_transform(const Expr& f, const Expr& t, const Expr& s) {
    if (!f) return S::Zero();
    // Linearity over Add.
    if (f->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(f->args().size());
        for (const auto& a : f->args()) {
            terms.push_back(laplace_transform(a, t, s));
        }
        return add(std::move(terms));
    }
    if (auto r = laplace_term(f, t, s); r.has_value()) {
        return *r;
    }
    return function_symbol("LaplaceTransform")(f, t, s);
}

}  // namespace sympp
