#pragma once

// NumberSymbol — symbolic mathematical constant with arbitrary-precision
// numeric evaluation. Pi, E, EulerGamma, Catalan, ...
//
// Reference: sympy/core/numbers.py::NumberSymbol

#include <cstddef>
#include <string>

#include <mpfr.h>

#include <sympp/core/api.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// Distinguishes individual number-symbol singletons.
enum class NumberSymbolKind : std::uint8_t {
    Pi,
    E,
    EulerGamma,
    Catalan,
};

class SYMPP_EXPORT NumberSymbol : public Number {
public:
    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::NumberSymbol; }
    [[nodiscard]] virtual NumberSymbolKind kind() const noexcept = 0;

    // All standard NumberSymbols are real and positive (Pi, E, EulerGamma,
    // Catalan, GoldenRatio, ...). Override if a future symbol differs.
    [[nodiscard]] bool is_zero() const noexcept override { return false; }
    [[nodiscard]] bool is_one() const noexcept override { return false; }
    [[nodiscard]] bool is_negative() const noexcept override { return false; }
    [[nodiscard]] bool is_positive() const noexcept override { return true; }
    [[nodiscard]] bool is_integer() const noexcept override { return false; }
    [[nodiscard]] bool is_rational() const noexcept override { return false; }
    [[nodiscard]] int sign() const noexcept override { return 1; }

    [[nodiscard]] bool equals(const Basic& other) const noexcept override;

    // Write the constant's value into `out` at precision `prec` (bits).
    // Used by evalf() to materialize a Float at a requested precision.
    virtual void mpfr_evalf(mpfr_t out, mpfr_prec_t prec) const = 0;
};

// Concrete singletons. Each class has exactly one canonical instance returned
// by S::pi(), S::e(), etc.
class SYMPP_EXPORT PiSymbol final : public NumberSymbol {
public:
    PiSymbol();
    [[nodiscard]] NumberSymbolKind kind() const noexcept override { return NumberSymbolKind::Pi; }
    [[nodiscard]] std::size_t hash() const noexcept override;
    [[nodiscard]] std::string str() const override { return "pi"; }
    void mpfr_evalf(mpfr_t out, mpfr_prec_t prec) const override;
};

class SYMPP_EXPORT ESymbol final : public NumberSymbol {
public:
    ESymbol();
    [[nodiscard]] NumberSymbolKind kind() const noexcept override { return NumberSymbolKind::E; }
    [[nodiscard]] std::size_t hash() const noexcept override;
    [[nodiscard]] std::string str() const override { return "E"; }
    void mpfr_evalf(mpfr_t out, mpfr_prec_t prec) const override;
};

class SYMPP_EXPORT EulerGammaSymbol final : public NumberSymbol {
public:
    EulerGammaSymbol();
    [[nodiscard]] NumberSymbolKind kind() const noexcept override {
        return NumberSymbolKind::EulerGamma;
    }
    [[nodiscard]] std::size_t hash() const noexcept override;
    [[nodiscard]] std::string str() const override { return "EulerGamma"; }
    void mpfr_evalf(mpfr_t out, mpfr_prec_t prec) const override;
};

class SYMPP_EXPORT CatalanSymbol final : public NumberSymbol {
public:
    CatalanSymbol();
    [[nodiscard]] NumberSymbolKind kind() const noexcept override {
        return NumberSymbolKind::Catalan;
    }
    [[nodiscard]] std::size_t hash() const noexcept override;
    [[nodiscard]] std::string str() const override { return "Catalan"; }
    void mpfr_evalf(mpfr_t out, mpfr_prec_t prec) const override;
};

}  // namespace sympp
