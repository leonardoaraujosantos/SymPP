#include <sympp/core/refine.hpp>

#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/assumption_key.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/miscellaneous.hpp>

namespace sympp {

namespace {

[[nodiscard]] Expr refine_abs(const Expr& original, const Expr& arg);

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
    if (ok) return pow(inner_base, mul(inner_exp, exp));

    // √(g²) → |g| for a real (or imaginary) base whose sign is not yet known —
    // the principal square root of a square is the absolute value, not g. A
    // known sign then collapses |g| to ±g via refine_abs. Matches SymPy, which
    // leaves it as √(g²) without a real assumption.
    if (inner_exp == integer(2) && exp == rational(1, 2)
        && (is_real(inner_base) == true || is_imaginary(inner_base) == true)) {
        return refine_abs(abs(inner_base), inner_base);
    }
    return pow(base, exp);
}

// |g| → g when g is known nonnegative, −g when g is known nonpositive. The sign
// facts usually come from an active `assuming` scope (refine(Abs(x)) is x under
// Q.positive(x), −x under Q.negative(x)). For a product the absolute value
// distributes, |∏fᵢ| = ∏|fᵢ| (valid for any complex factors), so each factor's
// sign is applied independently — refine(|x·y|) is x·|y| under Q.positive(x), and
// x·y under Q.positive(x)∧Q.positive(y). The distributed form is returned only
// when a factor actually collapsed, so |x·y| with no facts stays |x·y|. Mirrors
// SymPy's refine_abs.
[[nodiscard]] Expr refine_abs(const Expr& original, const Expr& arg) {
    if (is_nonnegative(arg) == true) return arg;
    if (is_nonpositive(arg) == true) return mul(S::NegativeOne(), arg);
    if (arg->type_id() == TypeId::Mul) {
        std::vector<Expr> factors;
        factors.reserve(arg->args().size());
        bool collapsed = false;
        for (const auto& f : arg->args()) {
            if (is_nonnegative(f) == true) {
                factors.push_back(f);
                collapsed = true;
            } else if (is_nonpositive(f) == true) {
                factors.push_back(mul(S::NegativeOne(), f));
                collapsed = true;
            } else {
                factors.push_back(abs(f));  // sign unknown → keep |f|
            }
        }
        if (collapsed) return mul(std::move(factors));
    }
    return original;
}

[[nodiscard]] Expr rebuild_with_args(const Expr& original, std::vector<Expr> new_args) {
    switch (original->type_id()) {
        case TypeId::Add:
            return add(std::move(new_args));
        case TypeId::Mul:
            return mul(std::move(new_args));
        case TypeId::Pow:
            return refine_pow(new_args[0], new_args[1]);
        case TypeId::Function: {
            const auto& f = static_cast<const Function&>(*original);
            if (f.function_id() == FunctionId::Abs && new_args.size() == 1) {
                return refine_abs(f.rebuild(std::vector<Expr>{new_args[0]}),
                                  new_args[0]);
            }
            return f.rebuild(std::move(new_args));
        }
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
    // (e.g. a Pow whose inner base just *is* positive at the leaf level, or an
    // Abs whose argument has a known sign from an `assuming` scope).
    if (e->type_id() == TypeId::Pow) {
        return refine_pow(new_args[0], new_args[1]);
    }
    if (e->type_id() == TypeId::Function) {
        const auto& f = static_cast<const Function&>(*e);
        if (f.function_id() == FunctionId::Abs && new_args.size() == 1) {
            return refine_abs(e, new_args[0]);
        }
    }
    if (!any_changed) return e;
    return rebuild_with_args(e, std::move(new_args));
}

}  // namespace sympp
