#pragma once

// Basic — abstract base for every node in the SymPP expression tree.
//
// Invariants:
//   * Immutable after construction. All "transforms" produce new nodes.
//   * Identity = structural equality + hash. Hash-cons cache (Phase 1) will
//     deduplicate structurally-equal subtrees.
//   * args() returns this node's direct children. Atomics return an empty span.

#include <cstddef>
#include <memory>
#include <span>
#include <string>

#include <sympp/core/api.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Basic : public std::enable_shared_from_this<Basic> {
public:
    Basic(const Basic&) = delete;
    Basic& operator=(const Basic&) = delete;
    Basic(Basic&&) = delete;
    Basic& operator=(Basic&&) = delete;

    virtual ~Basic() = default;

    [[nodiscard]] virtual TypeId type_id() const noexcept = 0;
    [[nodiscard]] virtual std::size_t hash() const noexcept = 0;
    [[nodiscard]] virtual std::span<const Expr> args() const noexcept = 0;

    // Structural equality. Default impl: type_id + hash + recursive args check.
    [[nodiscard]] virtual bool equals(const Basic& other) const noexcept;

    // String representation suitable for sympify-ing back through SymPy.
    // Phase 1 will introduce a dedicated Printer visitor; for now Basic owns it.
    [[nodiscard]] virtual std::string str() const = 0;

    [[nodiscard]] bool is_atomic() const noexcept { return args().empty(); }

protected:
    Basic() = default;
};

// Equality on Expr (shared_ptr) is structural, not pointer identity.
[[nodiscard]] inline bool operator==(const Expr& a, const Expr& b) noexcept {
    if (a.get() == b.get()) return true;
    if (!a || !b) return false;
    return a->equals(*b);
}

[[nodiscard]] inline bool operator!=(const Expr& a, const Expr& b) noexcept {
    return !(a == b);
}

// Same idea at the Basic level so dereferenced comparisons work too.
[[nodiscard]] inline bool operator==(const Basic& a, const Basic& b) noexcept {
    return a.equals(b);
}

[[nodiscard]] inline bool operator!=(const Basic& a, const Basic& b) noexcept {
    return !a.equals(b);
}

// Forward-declaration to break the include cycle between basic.hpp and
// expr_cache.hpp. The implementation lives in expr_cache.cpp.
class ExprCache;

namespace detail {
SYMPP_EXPORT Expr intern_through_cache(Expr candidate);
}

// Factory helper — every newly-constructed node is interned through the
// global ExprCache so structurally-equal subtrees share one shared_ptr.
//
// On a cache hit, `candidate` is deallocated immediately when this function
// returns; the returned shared_ptr aliases the cached object. The
// static_pointer_cast is safe because structural equality requires
// matching type_id, which in turn requires matching dynamic type.
template <typename T, typename... Args>
[[nodiscard]] std::shared_ptr<const T> make(Args&&... args) {
    Expr candidate = std::make_shared<const T>(std::forward<Args>(args)...);
    Expr interned = detail::intern_through_cache(std::move(candidate));
    return std::static_pointer_cast<const T>(std::move(interned));
}

}  // namespace sympp
