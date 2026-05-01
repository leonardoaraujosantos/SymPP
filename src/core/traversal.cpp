#include <sympp/core/traversal.hpp>

#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/boolean.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/piecewise.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// Reconstruct an expression of the same shape as `original` but with new args.
// For atomics (no args), returns `original` unchanged.
//
// As more compound types are added (Function, Derivative, Matrix...), extend
// this dispatch. A virtual `Basic::rebuild` would also work; leaving it as a
// switch keeps Basic's interface narrow.
[[nodiscard]] Expr rebuild_with_args(const Expr& original, std::vector<Expr> new_args) {
    switch (original->type_id()) {
        case TypeId::Add:
            return add(std::move(new_args));
        case TypeId::Mul:
            return mul(std::move(new_args));
        case TypeId::Pow:
            return pow(new_args[0], new_args[1]);
        case TypeId::Function: {
            const auto& f = static_cast<const Function&>(*original);
            return f.rebuild(std::move(new_args));
        }
        case TypeId::Relational: {
            const auto& r = static_cast<const Relational&>(*original);
            switch (r.kind()) {
                case RelKind::Eq: return eq(new_args[0], new_args[1]);
                case RelKind::Ne: return ne(new_args[0], new_args[1]);
                case RelKind::Lt: return lt(new_args[0], new_args[1]);
                case RelKind::Le: return le(new_args[0], new_args[1]);
                case RelKind::Gt: return gt(new_args[0], new_args[1]);
                case RelKind::Ge: return ge(new_args[0], new_args[1]);
            }
            return original;
        }
        case TypeId::Piecewise: {
            std::vector<PiecewiseBranch> branches;
            for (std::size_t i = 0; i + 1 < new_args.size(); i += 2) {
                branches.push_back({new_args[i], new_args[i + 1]});
            }
            return piecewise(std::move(branches));
        }
        default:
            return original;  // atomic — no args to rewrite
    }
}

// Recursive xreplace driver. Stops descending when a node hits the map.
[[nodiscard]] Expr xreplace_impl(const Expr& e, const ExprMap<Expr>& mapping) {
    auto it = mapping.find(e);
    if (it != mapping.end()) {
        return it->second;
    }

    auto args = e->args();
    if (args.empty()) return e;  // atomic, not in map → unchanged

    bool any_changed = false;
    std::vector<Expr> new_args;
    new_args.reserve(args.size());
    for (const auto& a : args) {
        Expr replaced = xreplace_impl(a, mapping);
        if (replaced.get() != a.get()) any_changed = true;
        new_args.push_back(std::move(replaced));
    }
    if (!any_changed) return e;
    return rebuild_with_args(e, std::move(new_args));
}

}  // namespace

Expr xreplace(const Expr& expr, const ExprMap<Expr>& mapping) {
    if (mapping.empty()) return expr;
    return xreplace_impl(expr, mapping);
}

Expr subs(const Expr& expr, const Expr& old_value, const Expr& new_value) {
    ExprMap<Expr> m;
    m.emplace(old_value, new_value);
    return xreplace_impl(expr, m);
}

bool has(const Expr& expr, const Expr& target) {
    if (expr == target) return true;
    for (const auto& a : expr->args()) {
        if (has(a, target)) return true;
    }
    return false;
}

namespace {

void collect_symbols(const Expr& e, ExprSet& out) {
    if (e->type_id() == TypeId::Symbol) {
        out.insert(e);
        return;
    }
    for (const auto& a : e->args()) {
        collect_symbols(a, out);
    }
}

}  // namespace

ExprSet free_symbols(const Expr& expr) {
    ExprSet out;
    collect_symbols(expr, out);
    return out;
}

}  // namespace sympp
