#pragma once

// Elementary trigonometric functions: sin, cos, tan.
//
// Auto-evaluation:
//   * sin(0) = 0, cos(0) = 1, tan(0) = 0
//   * sin(pi) = 0, cos(pi) = -1, tan(pi) = 0
//   * sin(pi/2) = 1, cos(pi/2) = 0
//   * sin(-x) = -sin(x), cos(-x) = cos(x), tan(-x) = -tan(x)
//   * Numeric arg → evalf via MPFR
//
// Inverse trigonometric (asin, acos, atan, atan2) and the reciprocal trio
// (cot, sec, csc) come in a follow-up sub-phase.
//
// Reference: sympy/functions/elementary/trigonometric.py

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Sin final : public Function {
public:
    explicit Sin(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Sin; }
    [[nodiscard]] std::string_view name() const noexcept override { return "sin"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Cos final : public Function {
public:
    explicit Cos(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Cos; }
    [[nodiscard]] std::string_view name() const noexcept override { return "cos"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Tan final : public Function {
public:
    explicit Tan(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Tan; }
    [[nodiscard]] std::string_view name() const noexcept override { return "tan"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// Factories — apply auto-eval rules. May return non-Sin/Cos/Tan results
// (numeric values, simpler expressions).
[[nodiscard]] SYMPP_EXPORT Expr sin(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr cos(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr tan(const Expr& arg);

}  // namespace sympp
