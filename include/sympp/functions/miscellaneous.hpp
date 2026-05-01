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

// ----- Complex part extraction ----------------------------------------------
//
// Without a first-class Complex node yet, these stay symbolic on a generic
// argument. Auto-eval rules:
//   * re(real) = real;     im(real) = 0
//   * conjugate(real) = real
//   * conjugate(conjugate(x)) = x   (involution)
//   * arg_(positive) = 0;   arg_(negative) = pi
//
// Reference: sympy/functions/elementary/complexes.py

class SYMPP_EXPORT Re final : public Function {
public:
    explicit Re(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Re; }
    [[nodiscard]] std::string_view name() const noexcept override { return "re"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Im final : public Function {
public:
    explicit Im(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Im; }
    [[nodiscard]] std::string_view name() const noexcept override { return "im"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Conjugate final : public Function {
public:
    explicit Conjugate(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Conjugate;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "conjugate"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Arg final : public Function {
public:
    explicit Arg(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Arg; }
    [[nodiscard]] std::string_view name() const noexcept override { return "arg"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

[[nodiscard]] SYMPP_EXPORT Expr re(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr im(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr conjugate(const Expr& arg);
// Note: 'arg' clashes with the parameter name `const Expr& arg`, so the public
// factory is named `arg_`. The class type is `Arg`.
[[nodiscard]] SYMPP_EXPORT Expr arg_(const Expr& arg);

// ----- Min, Max --------------------------------------------------------------
//
// Variadic. All-numeric arguments collapse to the min/max value; mixed-numeric
// + symbolic retains the symbolic args (with the numeric extreme kept) — full
// bounds reasoning is deferred.
//
// Reference: sympy/functions/elementary/miscellaneous.py — Min, Max

class SYMPP_EXPORT MinFn final : public Function {
public:
    explicit MinFn(std::vector<Expr> args);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Min; }
    [[nodiscard]] std::string_view name() const noexcept override { return "Min"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT MaxFn final : public Function {
public:
    explicit MaxFn(std::vector<Expr> args);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Max; }
    [[nodiscard]] std::string_view name() const noexcept override { return "Max"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

[[nodiscard]] SYMPP_EXPORT Expr min(std::vector<Expr> args);
[[nodiscard]] SYMPP_EXPORT Expr min(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Expr max(std::vector<Expr> args);
[[nodiscard]] SYMPP_EXPORT Expr max(const Expr& a, const Expr& b);

}  // namespace sympp
