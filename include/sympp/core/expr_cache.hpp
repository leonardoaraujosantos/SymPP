#pragma once

// ExprCache — global hash-cons cache.
//
// Every `make<T>(args)` call routes through `ExprCache::intern(candidate)`.
// If a structurally-equal Expr is already alive somewhere in the program,
// `intern` returns that one and the freshly-constructed candidate is dropped.
// Otherwise the candidate is recorded (as a weak_ptr) and returned.
//
// Consequences:
//   * Two structurally-equal Exprs always have identical shared_ptr addresses.
//   * `expr_a == expr_b`  becomes a pointer compare in the common case.
//   * Memory footprint shrinks when expressions share substructure.
//
// The cache holds weak_ptrs so it does not extend Expr lifetimes. When the
// last external shared_ptr drops, the cached weak_ptr expires and is GC'd
// the next time someone hashes into the same bucket (or when cleanup() is
// called explicitly).
//
// Reference: SymPy uses cache decorators on construction
//            (sympy/core/cache.py) plus SymEngine's hash-cons implementation
//            (symengine/basic.h::make_rcp).

#include <cstddef>
#include <memory>
#include <mutex>
#include <unordered_map>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT ExprCache {
public:
    [[nodiscard]] static ExprCache& instance();

    // Intern `candidate`. Returns either:
    //   * the previously-cached Expr that is structurally equal, or
    //   * `candidate` itself (newly inserted into the cache).
    //
    // Thread-safe.
    [[nodiscard]] Expr intern(Expr candidate);

    // Walk the cache and drop expired weak_ptrs. Mostly useful in tests and
    // long-running daemons; ordinary use covers GC lazily inside intern().
    std::size_t cleanup();

    // Approximate count of cached entries (live + expired).
    [[nodiscard]] std::size_t size() const noexcept;

    // Disable interning — every call to `intern` returns the candidate
    // unchanged. Useful for benchmarking or debugging memory leaks.
    void set_enabled(bool enabled) noexcept;
    [[nodiscard]] bool enabled() const noexcept;

private:
    ExprCache() = default;

    mutable std::mutex mtx_;
    std::unordered_multimap<std::size_t, std::weak_ptr<const Basic>> buckets_;
    bool enabled_ = true;
};

}  // namespace sympp
