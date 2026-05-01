#pragma once

// Order — Big-O asymptotic notation. O(expr, (var, point)) represents the
// asymptotic behavior of `expr` as `var` → `point`. Default point is 0.
//
// Minimal implementation:
//   * Class exists, hashes, equals, prints as O(expr) (variable inferred
//     from expr's free symbols when distinct).
//   * Auto-evaluation: O(0) → 0, O(constant) → O(1).
//   * Free function order(...) returns the canonical form.
//
// Full SymPy-style absorption rules (O(x^m) + O(x^n) → O(x^min(m,n)),
// constant * O(...) → O(...)) ship in a follow-up; this lays the
// foundation Series can use to denote a truncated remainder.
//
// Reference: sympy/series/order.py::Order

#include <array>
#include <cstddef>
#include <span>
#include <string>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Order final : public Basic {
public:
    // expr — the asymptotic generator (e.g. x**3 means O(x³)).
    // var, point — the limit specification; default point is 0.
    Order(Expr expr, Expr var, Expr point);

    [[nodiscard]] TypeId type_id() const noexcept override {
        return TypeId::Order;
    }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] std::span<const Expr> args() const noexcept override {
        return {args_.data(), args_.size()};
    }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;

    [[nodiscard]] const Expr& expr() const { return args_[0]; }
    [[nodiscard]] const Expr& var() const { return args_[1]; }
    [[nodiscard]] const Expr& point() const { return args_[2]; }

private:
    std::array<Expr, 3> args_;  // [expr, var, point]
    std::size_t hash_;
};

// Factory. Returns S::Zero() for O(0); otherwise an Order Expr. point
// defaults to S::Zero().
[[nodiscard]] SYMPP_EXPORT Expr order(const Expr& expr, const Expr& var);
[[nodiscard]] SYMPP_EXPORT Expr order(const Expr& expr, const Expr& var,
                                        const Expr& point);

}  // namespace sympp
