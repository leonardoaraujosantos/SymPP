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
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT Erfc final : public Function {
public:
    explicit Erfc(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Erfc; }
    [[nodiscard]] std::string_view name() const noexcept override { return "erfc"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
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
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
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

// Riemann zeta. Exact values: zeta(0)=-1/2, zeta(1)=zoo (pole), even positive
// integers zeta(2n)=rₙ·π^(2n) (Basel etc.), negative integers zeta(-n) rational
// (with the trivial zeros zeta(-2k)=0). Odd positive and symbolic args stay.
class SYMPP_EXPORT Zeta final : public Function {
public:
    explicit Zeta(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Zeta; }
    [[nodiscard]] std::string_view name() const noexcept override { return "zeta"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// Principal-branch Lambert W (inverse of x·eˣ). Exact values: W(0)=0, W(e)=1,
// W(-1/e)=-1, W(oo)=oo; other arguments stay symbolic. Derivative
// W'(x)=W(x)/(x·(1+W(x))).
class SYMPP_EXPORT LambertWFn final : public Function {
public:
    explicit LambertWFn(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::LambertW; }
    [[nodiscard]] std::string_view name() const noexcept override { return "LambertW"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

// Exponential / sine / cosine integrals — the antiderivatives of eˣ/x,
// sin(x)/x, cos(x)/x. Ei'(x)=eˣ/x, Si'(x)=sin(x)/x, Ci'(x)=cos(x)/x. Si is odd
// with Si(0)=0, Si(∞)=π/2; Ci/Ei stay symbolic for most arguments.
class SYMPP_EXPORT Ei final : public Function {
public:
    explicit Ei(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Ei; }
    [[nodiscard]] std::string_view name() const noexcept override { return "Ei"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT Si final : public Function {
public:
    explicit Si(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Si; }
    [[nodiscard]] std::string_view name() const noexcept override { return "Si"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT Ci final : public Function {
public:
    explicit Ci(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Ci; }
    [[nodiscard]] std::string_view name() const noexcept override { return "Ci"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

[[nodiscard]] SYMPP_EXPORT Expr erf(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr erfc(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr heaviside(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr dirac_delta(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr zeta(const Expr& s);
[[nodiscard]] SYMPP_EXPORT Expr lambertw(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr expint_ei(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr sinint(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr cosint(const Expr& arg);

}  // namespace sympp
