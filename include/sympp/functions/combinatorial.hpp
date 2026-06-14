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
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT Beta final : public Function {
public:
    Beta(Expr a, Expr b);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Beta;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "beta"; }
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
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

// polygamma(n, x) = dⁿ⁺¹/dxⁿ⁺¹ log Γ(x); polygamma(0, x) is the digamma
// function ψ(x). Kept symbolic for symbolic arguments, as SymPy does.
class SYMPP_EXPORT PolyGammaFn final : public Function {
public:
    PolyGammaFn(Expr n, Expr x);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::PolyGamma;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "polygamma"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT Fibonacci final : public Function {
public:
    explicit Fibonacci(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Fibonacci;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "fibonacci"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// Euler's totient φ(n): the count of integers in [1, n] coprime to n. Evaluates
// for a positive integer; stays symbolic otherwise. Mirrors SymPy's totient.
class SYMPP_EXPORT Totient final : public Function {
public:
    explicit Totient(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Totient;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "totient";
    }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Catalan final : public Function {
public:
    explicit Catalan(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Catalan;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "catalan"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT RisingFactorial final : public Function {
public:
    RisingFactorial(Expr x, Expr n);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::RisingFactorial;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "RisingFactorial"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT FallingFactorial final : public Function {
public:
    FallingFactorial(Expr x, Expr n);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::FallingFactorial;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "FallingFactorial"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Subfactorial final : public Function {
public:
    explicit Subfactorial(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Subfactorial;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "subfactorial"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Gcd final : public Function {
public:
    Gcd(Expr a, Expr b);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Gcd;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "gcd"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Lcm final : public Function {
public:
    Lcm(Expr a, Expr b);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Lcm;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "lcm"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

[[nodiscard]] SYMPP_EXPORT Expr factorial(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr binomial(const Expr& n, const Expr& k);
[[nodiscard]] SYMPP_EXPORT Expr beta(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Expr fibonacci(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr totient(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr catalan(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr gcd(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Expr lcm(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Expr rising_factorial(const Expr& x, const Expr& n);
[[nodiscard]] SYMPP_EXPORT Expr falling_factorial(const Expr& x, const Expr& n);
[[nodiscard]] SYMPP_EXPORT Expr subfactorial(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr gamma(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr loggamma(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr polygamma(const Expr& n, const Expr& x);
[[nodiscard]] SYMPP_EXPORT Expr digamma(const Expr& x);

}  // namespace sympp
