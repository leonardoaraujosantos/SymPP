#pragma once

// Rational — exact rational number p/q with gcd(p,q) == 1 and q > 0.
// Backed by GMP mpq_class.
//
// Construction always normalizes:
//   * denominator > 0 (sign carried by numerator)
//   * gcd(num, denom) == 1
//
// The factory `rational(num, den)` collapses to an Integer when denom == 1
// after reduction (matches SymPy's behaviour).
//
// Reference: sympy/sympy/core/numbers.py::Rational

#include <cstddef>
#include <string>
#include <string_view>

#include <gmpxx.h>

#include <sympp/core/api.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Rational final : public Number {
public:
    explicit Rational(mpq_class v);
    Rational(long num, long den);

    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Rational; }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;

    [[nodiscard]] bool is_zero() const noexcept override { return value_ == 0; }
    [[nodiscard]] bool is_one() const noexcept override { return value_ == 1; }
    [[nodiscard]] bool is_negative() const noexcept override { return sgn(value_) < 0; }
    [[nodiscard]] bool is_positive() const noexcept override { return sgn(value_) > 0; }
    [[nodiscard]] bool is_integer() const noexcept override {
        return value_.get_den() == 1;
    }
    [[nodiscard]] bool is_rational() const noexcept override { return true; }
    [[nodiscard]] int sign() const noexcept override { return sgn(value_); }

    [[nodiscard]] const mpq_class& value() const noexcept { return value_; }
    [[nodiscard]] mpz_class numerator() const { return value_.get_num(); }
    [[nodiscard]] mpz_class denominator() const { return value_.get_den(); }

private:
    mpq_class value_;
    std::size_t hash_;

    void compute_hash() noexcept;
};

// Factory: returns Integer if denom reduces to 1, else Rational.
// Throws std::domain_error on Rational(*, 0). (Future: collapse to
// ComplexInfinity once Phase 1e singletons exist.)
[[nodiscard]] SYMPP_EXPORT Expr rational(long num, long den);
[[nodiscard]] SYMPP_EXPORT Expr rational(const mpz_class& num, const mpz_class& den);

// Parse "p/q" or "p" forms.
[[nodiscard]] SYMPP_EXPORT Expr rational(std::string_view s);

}  // namespace sympp
