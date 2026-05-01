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

// Factory helper. In Phase 1 this routes through the hash-cons cache.
template <typename T, typename... Args>
[[nodiscard]] std::shared_ptr<const T> make(Args&&... args) {
    return std::make_shared<const T>(std::forward<Args>(args)...);
}

}  // namespace sympp
