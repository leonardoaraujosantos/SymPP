#pragma once

// Float — variable-precision floating-point number. Backed by MPFR.
//
// Precision is held in *bits* internally (MPFR native) but the public API
// accepts/reports decimal places (dps) for SymPy/mpmath compatibility.
//
// Conversion (matches mpmath):
//     prec_bits = max(53, ceil(dps * 3.33))
//     dps       = max(15, floor((prec_bits - 2) * 0.30103))
//
// Reference: sympy/sympy/core/numbers.py::Float

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

#include <gmp.h>
#include <mpfr.h>

#include <sympp/core/api.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// Default = 53 bits ≈ 15 decimal digits (matches IEEE 754 double).
constexpr mpfr_prec_t kDefaultPrecBits = 53;
constexpr int kDefaultDps = 15;

// Bits needed to represent `dps` decimal digits, with mpmath's safety pad.
// Reference: mpmath/libmp/libmpf.py::dps_to_prec
//     return max(1, int(round((dps+1)*log2(10)))-1)
// We use the canonical 3.32193 ≈ log2(10).
[[nodiscard]] inline mpfr_prec_t dps_to_prec(int dps) noexcept {
    if (dps <= kDefaultDps) return kDefaultPrecBits;
    auto bits = static_cast<mpfr_prec_t>(
        static_cast<double>(dps + 1) * 3.32192809488736 + 0.5);
    return bits < kDefaultPrecBits ? kDefaultPrecBits : bits;
}

class SYMPP_EXPORT Float final : public Number {
public:
    explicit Float(double v, int dps = kDefaultDps);
    explicit Float(std::string_view decimal, int dps = kDefaultDps);
    explicit Float(long v, int dps = kDefaultDps);
    explicit Float(mpz_srcptr v, int dps = kDefaultDps);
    explicit Float(mpq_srcptr v, int dps = kDefaultDps);
    explicit Float(mpfr_srcptr v, int dps = kDefaultDps);
    // Re-precision an existing Float to a new dps. Round-half-to-even.
    Float(const Float& src, int new_dps);
    Float(const Float& other);
    Float(Float&& other) noexcept;
    Float& operator=(const Float& other) = delete;  // Numbers are immutable
    Float& operator=(Float&& other) = delete;
    ~Float() override;

    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Float; }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;

    [[nodiscard]] bool is_zero() const noexcept override { return mpfr_zero_p(value_) != 0; }
    [[nodiscard]] bool is_one() const noexcept override {
        return mpfr_cmp_si(value_, 1) == 0;
    }
    [[nodiscard]] bool is_negative() const noexcept override { return mpfr_signbit(value_) != 0 && !is_zero(); }
    [[nodiscard]] bool is_positive() const noexcept override { return mpfr_signbit(value_) == 0 && !is_zero(); }
    [[nodiscard]] bool is_integer() const noexcept override {
        return mpfr_integer_p(value_) != 0;
    }
    [[nodiscard]] bool is_rational() const noexcept override { return false; }
    [[nodiscard]] int sign() const noexcept override {
        if (is_zero()) return 0;
        return mpfr_signbit(value_) != 0 ? -1 : 1;
    }

    [[nodiscard]] mpfr_prec_t precision_bits() const noexcept { return mpfr_get_prec(value_); }
    // The decimal-digit precision the user *requested*. We store this directly
    // because the bits→dps round-trip is intentionally lossy in mpmath.
    [[nodiscard]] int precision_dps() const noexcept { return dps_; }

    // Direct accessor for arithmetic and conversions. Caller must NOT modify.
    [[nodiscard]] mpfr_srcptr value() const noexcept { return value_; }

private:
    mpfr_t value_;
    int dps_;
    std::size_t hash_;

    void compute_hash() noexcept;
};

// Factories.
[[nodiscard]] inline Expr float_value(double v, int dps = kDefaultDps) {
    return make<Float>(v, dps);
}
[[nodiscard]] inline Expr float_value(std::string_view decimal, int dps = kDefaultDps) {
    return make<Float>(decimal, dps);
}

// evalf — convert a numeric Expr to a Float at the requested precision.
// For Integer / Rational / Float, returns a Float. For non-numeric, returns
// the input unchanged. (Phase 1f extends this to compound expressions.)
[[nodiscard]] SYMPP_EXPORT Expr evalf(const Expr& e, int dps = kDefaultDps);

}  // namespace sympp
