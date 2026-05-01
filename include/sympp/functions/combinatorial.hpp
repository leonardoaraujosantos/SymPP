#pragma once

// factorial, binomial, gamma, loggamma.
//
// Auto-evaluation:
//   * factorial(n) for nonnegative Integer n -> mpz_fac_ui (closed-form)
//   * factorial(0) = factorial(1) = 1
//   * binomial(n, k) for nonneg Integer n,k with k<=n -> mpz_bin_uiui
//   * gamma(n+1) = factorial(n) for nonneg Integer n
//   * gamma(1) = 1, gamma(1/2) = sqrt(pi)
//   * Numeric Float -> mpfr_gamma / mpfr_lngamma
//
// Reference: sympy/functions/combinatorial/factorials.py
//            sympy/functions/special/gamma_functions.py

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

class SYMPP_EXPORT Factorial final : public Function {
public:
    explicit Factorial(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Factorial;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "factorial"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Binomial final : public Function {
public:
    Binomial(Expr n, Expr k);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Binomial;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "binomial"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT GammaFn final : public Function {
public:
    explicit GammaFn(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Gamma;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "gamma"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT LogGamma final : public Function {
public:
    explicit LogGamma(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::LogGamma;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "loggamma"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

[[nodiscard]] SYMPP_EXPORT Expr factorial(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr binomial(const Expr& n, const Expr& k);
[[nodiscard]] SYMPP_EXPORT Expr gamma(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr loggamma(const Expr& arg);

}  // namespace sympp
