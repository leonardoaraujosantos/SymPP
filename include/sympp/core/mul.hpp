#pragma once

// Mul — symbolic product. Commutative, associative.
//
// Canonical form (Phase 1f.1):
//   * Flatten nested Muls                Mul(Mul(a,b), c) → Mul(a, b, c)
//   * Combine all numeric args           Mul(2, 3, x)     → Mul(6, x)
//   * Zero short-circuit                 Mul(0, x, y)     → 0
//   * Drop one factors                   Mul(1, x)        → x
//   * Collapse single-arg                Mul(x)           → x
//   * Stable arg ordering                (numbers first)
//
// Reference: sympy/core/mul.py::Mul

#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Mul final : public Basic {
public:
    explicit Mul(std::vector<Expr> args);

    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Mul; }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return args_; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;

private:
    std::vector<Expr> args_;
    std::size_t hash_;

    void compute_hash() noexcept;
};

[[nodiscard]] SYMPP_EXPORT Expr mul(std::vector<Expr> args);
[[nodiscard]] SYMPP_EXPORT Expr mul(const Expr& a, const Expr& b);

}  // namespace sympp
