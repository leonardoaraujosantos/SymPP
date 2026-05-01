#include <sympp/core/expand.hpp>

#include <cstddef>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// Distribute a single Add factor over the rest of a Mul.
// Precondition: `mul_args` is the canonical args of a Mul, and at least one
// of them is an Add. Recurses through `expand` so the result is fully
// distributed.
[[nodiscard]] Expr distribute_mul(std::vector<Expr> mul_args);

// Expand (base)^n when base is an Add and n is a non-negative Integer.
[[nodiscard]] Expr expand_int_pow_of_add(const Expr& base, long n);

[[nodiscard]] Expr distribute_mul(std::vector<Expr> mul_args) {
    // Find the first Add in the factor list.
    std::size_t add_idx = mul_args.size();
    for (std::size_t i = 0; i < mul_args.size(); ++i) {
        if (mul_args[i]->type_id() == TypeId::Add) {
            add_idx = i;
            break;
        }
    }
    if (add_idx == mul_args.size()) {
        // Nothing to distribute — re-canonicalize as Mul.
        return mul(std::move(mul_args));
    }

    // The Add factor we'll distribute over.
    Expr add_factor = std::move(mul_args[add_idx]);
    mul_args.erase(mul_args.begin() + static_cast<std::ptrdiff_t>(add_idx));

    // Build one Mul per Add term, expand each (other Adds in mul_args are
    // distributed by the recursive expand call).
    std::vector<Expr> sum_terms;
    sum_terms.reserve(add_factor->args().size());
    for (const auto& term : add_factor->args()) {
        std::vector<Expr> factors = mul_args;  // copy
        factors.push_back(term);
        sum_terms.push_back(expand(mul(std::move(factors))));
    }
    return add(std::move(sum_terms));
}

// Helper: split an Expr into its top-level Add terms (or itself if not an Add).
[[nodiscard]] std::vector<Expr> add_terms(const Expr& e) {
    if (e->type_id() == TypeId::Add) {
        return {e->args().begin(), e->args().end()};
    }
    return {e};
}

// Distribute (a) * (b) by computing every term-pair product. Both `a` and `b`
// may be Adds (or single terms). Each per-pair product goes through mul()
// (which collects bases) but never produces a new Pow-of-Add — neither
// factor is an Add — so we don't trigger expand_pow recursion.
[[nodiscard]] Expr distribute_two_adds(const Expr& a, const Expr& b) {
    auto av = add_terms(a);
    auto bv = add_terms(b);
    std::vector<Expr> terms;
    terms.reserve(av.size() * bv.size());
    for (const auto& at : av) {
        for (const auto& bt : bv) {
            terms.push_back(mul(at, bt));
        }
    }
    return add(std::move(terms));
}

[[nodiscard]] Expr expand_int_pow_of_add(const Expr& base, long n) {
    // Repeated multiplication, distributing manually so the Mul factory's
    // base-collection rule cannot rebuild a Pow(base, k) and re-trigger the
    // power-expansion path.
    if (n == 0) return S::One();
    if (n == 1) return base;
    Expr result = base;
    for (long i = 1; i < n; ++i) {
        result = distribute_two_adds(result, base);
    }
    return result;
}

[[nodiscard]] Expr expand_pow(const Expr& base, const Expr& exp) {
    // Only expand integer non-negative powers of an Add base.
    if (base->type_id() != TypeId::Add) {
        return pow(base, exp);
    }
    if (exp->type_id() != TypeId::Integer) {
        return pow(base, exp);
    }
    const auto& n = static_cast<const Integer&>(*exp);
    if (n.is_negative()) return pow(base, exp);
    if (!n.fits_long()) return pow(base, exp);
    return expand_int_pow_of_add(base, n.to_long());
}

}  // namespace

Expr expand(const Expr& e) {
    if (!e) return e;

    auto args = e->args();
    if (args.empty()) return e;  // atomic — nothing to expand

    // Expand each child first.
    std::vector<Expr> new_args;
    new_args.reserve(args.size());
    for (const auto& a : args) {
        new_args.push_back(expand(a));
    }

    // Apply expansion rule at this node.
    switch (e->type_id()) {
        case TypeId::Mul:
            return distribute_mul(std::move(new_args));
        case TypeId::Pow:
            return expand_pow(new_args[0], new_args[1]);
        case TypeId::Add: {
            // Sum doesn't distribute, but the sub-expansions may have changed
            // structure — re-canonicalize.
            return add(std::move(new_args));
        }
        default:
            // Other compound types (none yet in Phase 1) — return as-is.
            return e;
    }
}

}  // namespace sympp
