#include <sympp/simplify/simplify.hpp>

#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/expr_collections.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
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

Expr together(const Expr& e) {
    // Phase 5 minimal: canonical form already brings most expressions over
    // a single denominator via Mul base collection. Return the canonicalized
    // expression. Full common-denominator combining lands in a follow-up.
    return simplify(e);
}

}  // namespace sympp
