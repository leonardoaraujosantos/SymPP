#include <sympp/calculus/diff.hpp>

#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/exponential.hpp>

namespace sympp {

namespace {

[[nodiscard]] Expr diff_pow(const Expr& base, const Expr& exp_, const Expr& var) {
    auto db = diff(base, var);
    auto de = diff(exp_, var);

    // Constant exponent (independent of var): d/dx (a^n) = n*a^(n-1)*da/dx
    if (de == S::Zero()) {
        // n * a^(n-1) * da/dx
        auto new_exp = exp_ - S::One();
        return exp_ * pow(base, new_exp) * db;
    }
    // Constant base, var-dep exp: d/dx (a^b) = a^b * log(a) * db/dx
    if (db == S::Zero()) {
        return pow(base, exp_) * log(base) * de;
    }
    // General: d/dx (a^b) = a^b * (b' * log(a) + b * a' / a)
    auto term1 = de * log(base);
    auto term2 = exp_ * db / base;
    return pow(base, exp_) * (term1 + term2);
}

[[nodiscard]] Expr diff_function(const Function& fn, const Expr& var) {
    const auto& fargs = fn.args();
    // Sum_i  partial_i_f(args) * d/dx args[i]
    std::vector<Expr> terms;
    terms.reserve(fargs.size());
    for (std::size_t i = 0; i < fargs.size(); ++i) {
        Expr partial = fn.diff_arg(i);
        Expr inner = diff(fargs[i], var);
        if (partial == S::Zero() || inner == S::Zero()) continue;
        terms.push_back(mul(partial, inner));
    }
    return add(std::move(terms));
}

}  // namespace

Expr diff(const Expr& e, const Expr& var) {
    if (!e || !var) return S::Zero();

    // Atomic cases.
    if (e == var) return S::One();
    if (e->type_id() == TypeId::Symbol) return S::Zero();
    if (is_number(e)) return S::Zero();
    if (e->type_id() == TypeId::NumberSymbol) return S::Zero();

    // Sum rule.
    if (e->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(e->args().size());
        for (const auto& a : e->args()) {
            terms.push_back(diff(a, var));
        }
        return add(std::move(terms));
    }

    // Product rule.
    if (e->type_id() == TypeId::Mul) {
        const auto& factors = e->args();
        std::vector<Expr> terms;
        terms.reserve(factors.size());
        for (std::size_t i = 0; i < factors.size(); ++i) {
            std::vector<Expr> rest;
            rest.reserve(factors.size());
            for (std::size_t j = 0; j < factors.size(); ++j) {
                if (j == i) {
                    rest.push_back(diff(factors[j], var));
                } else {
                    rest.push_back(factors[j]);
                }
            }
            terms.push_back(mul(std::move(rest)));
        }
        return add(std::move(terms));
    }

    // Power rule.
    if (e->type_id() == TypeId::Pow) {
        return diff_pow(e->args()[0], e->args()[1], var);
    }

    // Function chain rule.
    if (e->type_id() == TypeId::Function) {
        return diff_function(static_cast<const Function&>(*e), var);
    }

    // Boolean / Relational / Piecewise — derivatives don't apply.
    return S::Zero();
}

Expr diff(const Expr& e, const Expr& var, std::size_t order) {
    Expr cur = e;
    for (std::size_t i = 0; i < order; ++i) {
        cur = diff(cur, var);
    }
    return cur;
}

}  // namespace sympp
