#pragma once

// Hyperbolic functions: sinh, cosh, tanh and inverses asinh, acosh, atanh.
//
// Auto-evaluation:
//   * sinh(0) = 0, cosh(0) = 1, tanh(0) = 0
//   * sinh(-x) = -sinh(x), tanh(-x) = -tanh(x), cosh(-x) = cosh(x)
//   * asinh(0) = 0, atanh(0) = 0, acosh(1) = 0
//   * asinh(-x) = -asinh(x), atanh(-x) = -atanh(x)
//   * Numeric Float arg → MPFR.
//
// Reference: sympy/functions/elementary/hyperbolic.py

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

class SYMPP_EXPORT Sinh final : public Function {
public:
    explicit Sinh(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Sinh; }
    [[nodiscard]] std::string_view name() const noexcept override { return "sinh"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT Cosh final : public Function {
public:
    explicit Cosh(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Cosh; }
    [[nodiscard]] std::string_view name() const noexcept override { return "cosh"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT Tanh final : public Function {
public:
    explicit Tanh(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Tanh; }
    [[nodiscard]] std::string_view name() const noexcept override { return "tanh"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT Asinh final : public Function {
public:
    explicit Asinh(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Asinh; }
    [[nodiscard]] std::string_view name() const noexcept override { return "asinh"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT Acosh final : public Function {
public:
    explicit Acosh(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Acosh; }
    [[nodiscard]] std::string_view name() const noexcept override { return "acosh"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT Atanh final : public Function {
public:
    explicit Atanh(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Atanh; }
    [[nodiscard]] std::string_view name() const noexcept override { return "atanh"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

[[nodiscard]] SYMPP_EXPORT Expr sinh(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr cosh(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr tanh(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr asinh(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr acosh(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr atanh(const Expr& arg);

}  // namespace sympp
