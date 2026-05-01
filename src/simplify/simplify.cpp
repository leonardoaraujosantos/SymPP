#include <sympp/simplify/simplify.hpp>

#include <algorithm>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/boolean.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/expr_collections.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/piecewise.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/polys/poly.hpp>

namespace sympp {

namespace {

// Re-run the canonical-form factories on every node — many substitution
// chains leave nominally-canonical Expr that benefit from a fresh sweep.
[[nodiscard]] Expr re_canonicalize(const Expr& e) {
    auto args = e->args();
    if (args.empty()) return e;
    std::vector<Expr> new_args;
    new_args.reserve(args.size());
    for (const auto& a : args) new_args.push_back(re_canonicalize(a));
    switch (e->type_id()) {
        case TypeId::Add: return add(std::move(new_args));
        case TypeId::Mul: return mul(std::move(new_args));
        case TypeId::Pow: return pow(new_args[0], new_args[1]);
        default: return e;
    }
}

}  // namespace

Expr simplify(const Expr& e) {
    if (!e) return e;
    // Sweep 1: canonical form.
    Expr current = re_canonicalize(e);
    // Sweep 2: expand to flush nested products.
    Expr expanded = expand(current);
    // Sweep 3: canonical again to collect like terms after expansion.
    return re_canonicalize(expanded);
}

Expr collect(const Expr& e, const Expr& var) {
    // Build a Poly in var; reconstitute as a sum of (coeff * var^k).
    Poly p(expand(e), var);
    return p.as_expr();
}

// ----- powsimp ---------------------------------------------------------------

namespace {

// Reconstruct an expression of the same shape as `original` with new args.
// Mirrors traversal.cpp::rebuild_with_args; duplicated here to avoid making
// that helper public.
[[nodiscard]] Expr rebuild_node(const Expr& original,
                                std::vector<Expr> new_args) {
    switch (original->type_id()) {
        case TypeId::Add: return add(std::move(new_args));
        case TypeId::Mul: return mul(std::move(new_args));
        case TypeId::Pow: return pow(new_args[0], new_args[1]);
        case TypeId::Function: {
            const auto& f = static_cast<const Function&>(*original);
            return f.rebuild(std::move(new_args));
        }
        default: return original;
    }
}

// True iff combining bases under a shared non-integer exponent is safe.
// Conservative: requires every base to be definitely-positive.
[[nodiscard]] bool can_combine_under_exp(const std::vector<Expr>& bases,
                                          const Expr& exp) {
    if (exp->type_id() == TypeId::Integer) return true;
    for (const auto& b : bases) {
        auto p = is_positive(b);
        if (!p.has_value() || !*p) return false;
    }
    return true;
}

[[nodiscard]] Expr powsimp_node(const Expr& e) {
    if (e->type_id() != TypeId::Mul) return e;

    // Group factors by their exponent (as_base_exp).
    struct Group { Expr exp; std::vector<Expr> bases; };
    std::vector<Group> groups;
    auto find_or_create = [&](const Expr& exp) -> Group& {
        for (auto& g : groups) if (g.exp == exp) return g;
        groups.push_back(Group{exp, {}});
        return groups.back();
    };

    for (const auto& a : e->args()) {
        Expr base, exp;
        if (a->type_id() == TypeId::Pow) {
            base = a->args()[0];
            exp = a->args()[1];
        } else {
            base = a;
            exp = S::One();
        }
        find_or_create(exp).bases.push_back(base);
    }

    // Reassemble — merge groups with ≥ 2 bases when safe.
    std::vector<Expr> out;
    out.reserve(groups.size());
    for (auto& g : groups) {
        if (g.bases.size() == 1) {
            out.push_back(pow(g.bases[0], g.exp));
        } else if (can_combine_under_exp(g.bases, g.exp)) {
            out.push_back(pow(mul(g.bases), g.exp));
        } else {
            for (auto& b : g.bases) out.push_back(pow(b, g.exp));
        }
    }
    return mul(std::move(out));
}

}  // namespace

Expr powsimp(const Expr& e) {
    if (!e) return e;
    auto args = e->args();
    if (args.empty()) return e;
    std::vector<Expr> new_args;
    new_args.reserve(args.size());
    for (const auto& a : args) new_args.push_back(powsimp(a));
    Expr rebuilt = rebuild_node(e, std::move(new_args));
    return powsimp_node(rebuilt);
}

// ----- expand_power_exp / expand_power_base ----------------------------------

namespace {

[[nodiscard]] Expr expand_power_exp_node(const Expr& e) {
    if (e->type_id() != TypeId::Pow) return e;
    const Expr& base = e->args()[0];
    const Expr& exp = e->args()[1];
    if (exp->type_id() != TypeId::Add) return e;
    // x^(a + b + c) → x^a * x^b * x^c
    std::vector<Expr> factors;
    factors.reserve(exp->args().size());
    for (const auto& term : exp->args()) {
        factors.push_back(pow(base, term));
    }
    return mul(std::move(factors));
}

[[nodiscard]] Expr expand_power_base_node(const Expr& e) {
    if (e->type_id() != TypeId::Pow) return e;
    const Expr& base = e->args()[0];
    const Expr& exp = e->args()[1];
    if (base->type_id() != TypeId::Mul) return e;
    // Safe to distribute when exp is integer, or every factor of base is
    // positive.
    bool safe = (exp->type_id() == TypeId::Integer);
    if (!safe) {
        bool all_pos = true;
        for (const auto& f : base->args()) {
            auto p = is_positive(f);
            if (!p.has_value() || !*p) { all_pos = false; break; }
        }
        safe = all_pos;
    }
    if (!safe) return e;
    std::vector<Expr> factors;
    factors.reserve(base->args().size());
    for (const auto& f : base->args()) factors.push_back(pow(f, exp));
    return mul(std::move(factors));
}

template <typename NodeFn>
[[nodiscard]] Expr apply_recursive(const Expr& e, NodeFn fn) {
    if (!e) return e;
    auto args = e->args();
    if (args.empty()) return fn(e);
    std::vector<Expr> new_args;
    new_args.reserve(args.size());
    for (const auto& a : args) new_args.push_back(apply_recursive(a, fn));
    return fn(rebuild_node(e, std::move(new_args)));
}

}  // namespace

Expr expand_power_exp(const Expr& e) {
    return apply_recursive(e, expand_power_exp_node);
}

Expr expand_power_base(const Expr& e) {
    return apply_recursive(e, expand_power_base_node);
}

// ----- trigsimp --------------------------------------------------------------

namespace {

// If `e` is sin(arg)^2 or cos(arg)^2, return (is_sin, arg). Else nullopt.
[[nodiscard]] std::optional<std::pair<bool, Expr>>
detect_trig_pow2(const Expr& e) {
    if (e->type_id() != TypeId::Pow) return std::nullopt;
    if (!(e->args()[1] == integer(2))) return std::nullopt;
    const Expr& base = e->args()[0];
    if (base->type_id() != TypeId::Function) return std::nullopt;
    const auto& f = static_cast<const Function&>(*base);
    if (f.function_id() == FunctionId::Sin) return std::pair{true, base->args()[0]};
    if (f.function_id() == FunctionId::Cos) return std::pair{false, base->args()[0]};
    return std::nullopt;
}

// Decompose a term as (coef * sin(arg)^2) or (coef * cos(arg)^2). Returns
// nullopt when the term has neither, or has more than one such factor (an
// unusual structure we leave alone).
struct TrigSquareTerm {
    bool is_sin;
    Expr arg;
    Expr coef;
};
[[nodiscard]] std::optional<TrigSquareTerm>
as_trig_square_term(const Expr& term) {
    if (auto r = detect_trig_pow2(term); r) {
        return TrigSquareTerm{r->first, r->second, S::One()};
    }
    if (term->type_id() == TypeId::Mul) {
        std::optional<std::pair<bool, Expr>> trig;
        std::vector<Expr> coef_factors;
        for (const auto& f : term->args()) {
            if (auto r = detect_trig_pow2(f); r) {
                if (trig) return std::nullopt;  // multiple — leave alone
                trig = r;
            } else {
                coef_factors.push_back(f);
            }
        }
        if (trig) {
            return TrigSquareTerm{trig->first, trig->second, mul(coef_factors)};
        }
    }
    return std::nullopt;
}

[[nodiscard]] Expr trigsimp_add(const Expr& e) {
    if (e->type_id() != TypeId::Add) return e;

    struct CoefPair { Expr sin_coef = S::Zero(); Expr cos_coef = S::Zero(); };
    std::vector<std::pair<Expr, CoefPair>> by_arg;
    std::vector<Expr> non_trig;

    auto find_or_create = [&](const Expr& arg) -> CoefPair& {
        for (auto& [a, p] : by_arg) if (a == arg) return p;
        by_arg.push_back({arg, CoefPair{}});
        return by_arg.back().second;
    };

    for (const auto& term : e->args()) {
        auto t = as_trig_square_term(term);
        if (!t) { non_trig.push_back(term); continue; }
        auto& cp = find_or_create(t->arg);
        if (t->is_sin) cp.sin_coef = add(cp.sin_coef, t->coef);
        else cp.cos_coef = add(cp.cos_coef, t->coef);
    }

    std::vector<Expr> out = std::move(non_trig);
    for (auto& [arg, cp] : by_arg) {
        // Rewrite a*sin²(x) + b*cos²(x) as b + (a-b)*sin²(x). When a == b
        // this collapses to b — the Pythagorean identity. Otherwise it
        // produces an equivalent form with one fewer trig² (b absorbed).
        out.push_back(cp.cos_coef);
        Expr diff = cp.sin_coef - cp.cos_coef;
        if (!(diff == S::Zero())) {
            out.push_back(mul(diff, pow(sin(arg), integer(2))));
        }
    }
    if (out.empty()) return S::Zero();
    if (out.size() == 1) return out[0];
    return add(std::move(out));
}

}  // namespace

Expr trigsimp(const Expr& e) {
    return apply_recursive(e, trigsimp_add);
}

}  // namespace sympp
