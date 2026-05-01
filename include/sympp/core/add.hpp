#pragma once

// Add — symbolic sum. Commutative, associative.
//
// Canonical form (Phase 1f.1):
//   * Flatten nested Adds                Add(Add(a,b), c) → Add(a, b, c)
//   * Combine all numeric args           Add(1, 2, x)     → Add(3, x)
//   * Drop zero terms                    Add(0, x)        → x
//   * Collapse single-arg                Add(x)           → x
//   * Stable arg ordering                (numbers last)
//
// Reference: sympy/core/add.py::Add and sympy/core/operations.py::AssocOp

#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Add final : public Basic {
public:
    explicit Add(std::vector<Expr> args);

    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Add; }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return args_; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;

private:
    std::vector<Expr> args_;
    std::size_t hash_;

    void compute_hash() noexcept;
};

// Factory — applies the canonical-form rules above.
[[nodiscard]] SYMPP_EXPORT Expr add(std::vector<Expr> args);
[[nodiscard]] SYMPP_EXPORT Expr add(const Expr& a, const Expr& b);

}  // namespace sympp
