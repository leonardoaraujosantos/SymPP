#pragma once

// sqrt, Abs, sign — workhorses for "real-number" reasoning.
//
// sqrt(x) returns Pow(x, 1/2). It is *not* a separate node type in SymPP
// (matches SymPy's modeling). The factory applies basic auto-eval:
//   * sqrt(0) = 0
//   * sqrt(1) = 1
//   * sqrt(perfect square integer) = integer
//
// Abs(x) is a Function. Auto-eval:
//   * Abs(real number) → numeric absolute value
//   * Abs(x) for x nonnegative → x
//   * Abs(x) for x nonpositive → -x
//   * Abs(-x) → Abs(x)
//
// sign(x) is a Function. Auto-eval:
//   * sign(0) = 0; sign(positive) = 1; sign(negative) = -1
//   * sign(-x) = -sign(x)
//
// Reference: sympy/functions/elementary/{miscellaneous,complexes}.py

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

[[nodiscard]] SYMPP_EXPORT Expr sqrt(const Expr& arg);

class SYMPP_EXPORT Abs final : public Function {
public:
    explicit Abs(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Abs; }
    [[nodiscard]] std::string_view name() const noexcept override { return "Abs"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Sign final : public Function {
public:
    explicit Sign(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Sign; }
    [[nodiscard]] std::string_view name() const noexcept override { return "sign"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

[[nodiscard]] SYMPP_EXPORT Expr abs(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr sign(const Expr& arg);

}  // namespace sympp
