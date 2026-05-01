#pragma once

// Practical "special functions" needed for MATLAB toolbox examples.
//
// erf(x) — error function. erf(0)=0, erf(-x)=-erf(x). Numeric via mpfr_erf.
// erfc(x) = 1 - erf(x). erfc(0)=1.
// Heaviside(x) — unit step. Heaviside(0) = 1/2 (SymPy default H0=1/2),
//   Heaviside(positive)=1, Heaviside(negative)=0.
// DiracDelta(x) — distributional delta. DiracDelta(nonzero) = 0;
//   DiracDelta(0) stays unevaluated (it's not really a function value).
//
// Reference: sympy/functions/special/error_functions.py,
//            sympy/functions/special/delta_functions.py

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

class SYMPP_EXPORT Erf final : public Function {
public:
    explicit Erf(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Erf; }
    [[nodiscard]] std::string_view name() const noexcept override { return "erf"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Erfc final : public Function {
public:
    explicit Erfc(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Erfc; }
    [[nodiscard]] std::string_view name() const noexcept override { return "erfc"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT HeavisideFn final : public Function {
public:
    explicit HeavisideFn(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Heaviside;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "Heaviside"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT DiracDeltaFn final : public Function {
public:
    explicit DiracDeltaFn(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::DiracDelta;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "DiracDelta"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

[[nodiscard]] SYMPP_EXPORT Expr erf(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr erfc(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr heaviside(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr dirac_delta(const Expr& arg);

}  // namespace sympp
