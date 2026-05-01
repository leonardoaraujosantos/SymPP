#pragma once

// Pow — power expression base^exp.
//
// Auto-eval rules (Phase 1f.1):
//   * x^0      → 1   (including 0^0 per SymPy convention)
//   * x^1      → x
//   * 1^x      → 1
//   * 0^x      → 0  when x is a positive Number
//   * Number^Integer-exp computed exactly via number_pow
//
// Cases that stay unevaluated:
//   * Symbolic base or exponent
//   * Non-Integer exponent on Integer/Rational base (no real-domain decision yet)
//   * Negative-zero exponent edge cases
//
// Reference: sympy/core/power.py::Pow

#include <cstddef>
#include <span>
#include <string>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Pow final : public Basic {
public:
    Pow(Expr base, Expr exp);

    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Pow; }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return args_; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;

    [[nodiscard]] const Expr& base() const noexcept { return args_[0]; }
    [[nodiscard]] const Expr& exp() const noexcept { return args_[1]; }

private:
    Expr args_[2];
    std::size_t hash_;
};

// Factory. Applies auto-eval rules; may return an Expr that is not a Pow.
[[nodiscard]] SYMPP_EXPORT Expr pow(const Expr& base, const Expr& exp);

}  // namespace sympp
