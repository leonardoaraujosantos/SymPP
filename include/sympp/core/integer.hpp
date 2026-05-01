#pragma once

// Integer — arbitrary-precision integer. Backed by GMP mpz_class.
// Reference: sympy/sympy/core/numbers.py::Integer

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

#include <gmpxx.h>

#include <sympp/core/api.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Integer final : public Number {
public:
    explicit Integer(int v);
    explicit Integer(long v);
    explicit Integer(long long v);
    explicit Integer(unsigned long v);
    explicit Integer(unsigned long long v);
    explicit Integer(std::string_view decimal);
    explicit Integer(mpz_class v);

    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Integer; }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;

    [[nodiscard]] bool is_zero() const noexcept override { return value_ == 0; }
    [[nodiscard]] bool is_one() const noexcept override { return value_ == 1; }
    [[nodiscard]] bool is_negative() const noexcept override { return sgn(value_) < 0; }
    [[nodiscard]] bool is_positive() const noexcept override { return sgn(value_) > 0; }
    [[nodiscard]] bool is_integer() const noexcept override { return true; }
    [[nodiscard]] bool is_rational() const noexcept override { return true; }
    [[nodiscard]] int sign() const noexcept override { return sgn(value_); }

    [[nodiscard]] const mpz_class& value() const noexcept { return value_; }

    // Best-effort narrowing. Returns false if the value doesn't fit.
    [[nodiscard]] bool fits_long() const noexcept { return value_.fits_slong_p() != 0; }
    [[nodiscard]] long to_long() const noexcept { return value_.get_si(); }

private:
    mpz_class value_;
    std::size_t hash_;

    void compute_hash() noexcept;
};

// Convenience factories.
[[nodiscard]] inline Expr integer(int v) {
    return make<Integer>(v);
}
[[nodiscard]] inline Expr integer(long v) {
    return make<Integer>(v);
}
[[nodiscard]] inline Expr integer(long long v) {
    return make<Integer>(v);
}
[[nodiscard]] inline Expr integer(std::string_view decimal) {
    return make<Integer>(decimal);
}

}  // namespace sympp
