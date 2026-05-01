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

// ----- Inverse trigonometric -------------------------------------------------
//
// Auto-evaluation:
//   * asin(0) = 0, asin(1) = pi/2, asin(-1) = -pi/2
//   * acos(0) = pi/2, acos(1) = 0, acos(-1) = pi
//   * atan(0) = 0, atan(1) = pi/4, atan(-1) = -pi/4
//   * asin(-x) = -asin(x); atan(-x) = -atan(x)
//   * acos(-x) = pi - acos(x)
//   * atan2(0, 0) stays unevaluated; atan2(y, 0) for sign-known y → ±pi/2,
//     atan2(0, x) for sign-known x → 0 or pi.

class SYMPP_EXPORT Asin final : public Function {
public:
    explicit Asin(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Asin; }
    [[nodiscard]] std::string_view name() const noexcept override { return "asin"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Acos final : public Function {
public:
    explicit Acos(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Acos; }
    [[nodiscard]] std::string_view name() const noexcept override { return "acos"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Atan final : public Function {
public:
    explicit Atan(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Atan; }
    [[nodiscard]] std::string_view name() const noexcept override { return "atan"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Atan2 final : public Function {
public:
    Atan2(Expr y, Expr x);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Atan2; }
    [[nodiscard]] std::string_view name() const noexcept override { return "atan2"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

[[nodiscard]] SYMPP_EXPORT Expr asin(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr acos(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr atan(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr atan2(const Expr& y, const Expr& x);

}  // namespace sympp
