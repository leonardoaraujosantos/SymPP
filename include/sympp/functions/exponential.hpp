#pragma once

// exp(x) and log(x) — natural exponential and logarithm.
//
// Auto-evaluation:
//   * exp(0) = 1
//   * exp(1) = E
//   * exp(log(x)) = x   (when x is positive)
//   * log(1) = 0
//   * log(E) = 1
//   * log(exp(x)) = x   (when x is real)
//   * Numeric Float arg → evalf via mpfr_exp / mpfr_log.
//
// Multi-arg log (with base) is deferred — most calls are natural log.
//
// Reference: sympy/functions/elementary/exponential.py

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

class SYMPP_EXPORT Exp final : public Function {
public:
    explicit Exp(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Exp; }
    [[nodiscard]] std::string_view name() const noexcept override { return "exp"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT Log final : public Function {
public:
    explicit Log(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Log; }
    [[nodiscard]] std::string_view name() const noexcept override { return "log"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

[[nodiscard]] SYMPP_EXPORT Expr exp(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr log(const Expr& arg);

}  // namespace sympp
