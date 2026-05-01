#include <sympp/core/refine.hpp>

#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/assumption_key.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// (x^a)^b → x^(a*b) when valid given current assumptions on x.
[[nodiscard]] Expr refine_pow(const Expr& base, const Expr& exp) {
    if (base->type_id() != TypeId::Pow) {
        return pow(base, exp);
    }
    const auto& inner = static_cast<const Pow&>(*base);
    const auto& inner_base = inner.base();
    const auto& inner_exp = inner.exp();

    // Accept the rewrite if either:
    //   * the innermost base is positive (real-domain laws apply for any
    //     real exponent), or
    //   * both exponents are integers (integer power composition is always
    //     valid by convention).
    bool ok = false;
    if (is_positive(inner_base) == true) {
        ok = true;
    } else if (is_integer(inner_exp) == true && is_integer(exp) == true) {
        ok = true;
    }
    if (!ok) return pow(base, exp);

    return pow(inner_base, mul(inner_exp, exp));
}

[[nodiscard]] Expr rebuild_with_args(const Expr& original, std::vector<Expr> new_args) {
    switch (original->type_id()) {
        case TypeId::Add:
            return add(std::move(new_args));
        case TypeId::Mul:
            return mul(std::move(new_args));
        case TypeId::Pow:
            return refine_pow(new_args[0], new_args[1]);
        default:
            return original;
    }
}

}  // namespace

Expr refine(const Expr& e) {
    if (!e) return e;
    auto args = e->args();
    if (args.empty()) return e;

    std::vector<Expr> new_args;
    new_args.reserve(args.size());
    bool any_changed = false;
    for (const auto& a : args) {
        Expr r = refine(a);
        if (r.get() != a.get()) any_changed = true;
        new_args.push_back(std::move(r));
    }

    // Even if the children didn't change, this node might still simplify
    // (e.g. a Pow whose inner base just *is* positive at the leaf level).
    if (e->type_id() == TypeId::Pow) {
        return refine_pow(new_args[0], new_args[1]);
    }
    if (!any_changed) return e;
    return rebuild_with_args(e, std::move(new_args));
}

}  // namespace sympp
