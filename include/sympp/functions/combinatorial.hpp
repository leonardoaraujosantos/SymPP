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
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

// Incomplete gamma functions. lowergamma(s,x) = ∫₀ˣ tˢ⁻¹e⁻ᵗ dt, uppergamma(s,x)
// = ∫ₓ^∞ tˢ⁻¹e⁻ᵗ dt; together they sum to Γ(s). For a positive-integer first
// argument both reduce to a closed elementary form, as SymPy does; otherwise
// they stay symbolic.
class SYMPP_EXPORT LowerGamma final : public Function {
public:
    LowerGamma(Expr s, Expr x);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::LowerGamma;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "lowergamma";
    }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT UpperGamma final : public Function {
public:
    UpperGamma(Expr s, Expr x);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::UpperGamma;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "uppergamma";
    }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
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

// Lucas number L(n): L(0)=2, L(1)=1, L(n)=L(n-1)+L(n-2). Evaluates for a
// non-negative integer; stays symbolic otherwise. Mirrors SymPy's lucas.
class SYMPP_EXPORT Lucas final : public Function {
public:
    explicit Lucas(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Lucas;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "lucas"; }
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

// prime(n): the n-th prime (prime(1)=2, prime(5)=11). Evaluates for a positive
// integer index; stays symbolic otherwise. Mirrors SymPy's prime.
class SYMPP_EXPORT Prime final : public Function {
public:
    explicit Prime(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Prime;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "prime";
    }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// primepi(n): the count of primes ≤ n. Evaluates for an integer; stays symbolic
// otherwise. Mirrors SymPy's primepi.
class SYMPP_EXPORT PrimePi final : public Function {
public:
    explicit PrimePi(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::PrimePi;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "primepi";
    }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// mobius(n): the Möbius function μ(n) — 0 if n has a squared prime factor, else
// (−1)^(number of distinct primes). Evaluates for a positive integer. Mirrors
// SymPy's mobius.
class SYMPP_EXPORT Mobius final : public Function {
public:
    explicit Mobius(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Mobius;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "mobius";
    }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// divisor_count(n): the number of positive divisors of n (σ₀). Mirrors SymPy.
class SYMPP_EXPORT DivisorCount final : public Function {
public:
    explicit DivisorCount(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::DivisorCount;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "divisor_count";
    }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// divisor_sigma(n): the sum of the positive divisors of n (σ₁). Mirrors SymPy's
// single-argument divisor_sigma.
class SYMPP_EXPORT DivisorSigma final : public Function {
public:
    explicit DivisorSigma(Expr arg);
    DivisorSigma(Expr n, Expr k);  // σ_k(n) = Σ_{d|n} d^k
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::DivisorSigma;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "divisor_sigma";
    }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// harmonic(n): the n-th harmonic number Hₙ = Σ_{k=1}^n 1/k. Evaluates for a
// non-negative integer; stays symbolic otherwise. Mirrors SymPy's harmonic.
class SYMPP_EXPORT Harmonic final : public Function {
public:
    explicit Harmonic(Expr arg);
    Harmonic(Expr n, Expr m);  // generalized Hₙ⁽ᵐ⁾ = Σ_{k=1}^n k^(−m)
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Harmonic;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "harmonic";
    }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

// factorial2(n): the double factorial n!! = n(n−2)(n−4)…. Evaluates for an
// integer ≥ −1 (factorial2(0)=factorial2(−1)=1); stays symbolic otherwise.
// Mirrors SymPy's factorial2.
class SYMPP_EXPORT Factorial2 final : public Function {
public:
    explicit Factorial2(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Factorial2;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "factorial2";
    }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// bernoulli(n): the n-th Bernoulli number Bₙ (SymPy's convention B₁ = +1/2).
// Evaluates for a non-negative integer; stays symbolic otherwise.
class SYMPP_EXPORT Bernoulli final : public Function {
public:
    explicit Bernoulli(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Bernoulli;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "bernoulli";
    }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// euler(n): the n-th Euler number Eₙ (1, 0, −1, 0, 5, …; odd indices are 0).
// Evaluates for a non-negative integer; stays symbolic otherwise.
class SYMPP_EXPORT Euler final : public Function {
public:
    explicit Euler(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Euler;
    }
    [[nodiscard]] std::string_view name() const noexcept override {
        return "euler";
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
[[nodiscard]] SYMPP_EXPORT Expr lucas(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr totient(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr prime(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr primepi(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr mobius(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr divisor_count(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr divisor_sigma(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr divisor_sigma(const Expr& n, const Expr& k);
[[nodiscard]] SYMPP_EXPORT Expr harmonic(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr harmonic(const Expr& n, const Expr& m);
[[nodiscard]] SYMPP_EXPORT Expr factorial2(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr bernoulli(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr euler(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr catalan(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr gcd(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Expr lcm(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Expr rising_factorial(const Expr& x, const Expr& n);
[[nodiscard]] SYMPP_EXPORT Expr falling_factorial(const Expr& x, const Expr& n);
[[nodiscard]] SYMPP_EXPORT Expr subfactorial(const Expr& arg);
// Bell number Bₙ (number of set partitions of an n-set) and the tribonacci
// number Tₙ (T₀=0, T₁=T₂=1, Tₙ=Tₙ₋₁+Tₙ₋₂+Tₙ₋₃); evaluated for non-negative
// integer arguments, symbolic otherwise.
[[nodiscard]] SYMPP_EXPORT Expr bell(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr tribonacci(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr gamma(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr loggamma(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr polygamma(const Expr& n, const Expr& x);
[[nodiscard]] SYMPP_EXPORT Expr digamma(const Expr& x);
[[nodiscard]] SYMPP_EXPORT Expr lowergamma(const Expr& s, const Expr& x);
[[nodiscard]] SYMPP_EXPORT Expr uppergamma(const Expr& s, const Expr& x);

}  // namespace sympp
