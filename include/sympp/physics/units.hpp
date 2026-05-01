#pragma once

// Phase 12 — Units module.
//
// Layered design:
//
//   Dimension   — 7-tuple of integer exponents over the SI base
//                 dimensions (length L, mass M, time T, current I,
//                 temperature Θ, amount N, luminosity J). Arithmetic:
//                 add → multiply (exponents add); negate → divide;
//                 zero → dimensionless.
//
//   Unit        — (name, symbol, dimension, scale). `scale` is the
//                 multiplier you must apply to a quantity expressed
//                 in this unit to convert it into base SI. So
//                 km has scale 1000; minute has scale 60. A unit
//                 with dimension Dimensionless::zero() and scale 1
//                 is "no unit at all".
//
//   Quantity    — (Expr value, Unit). Multiplication / division
//                 combine units; addition / subtraction require
//                 dimensionally-compatible units (otherwise throw).
//
// Plus convert(), unit-system aware constants (c, G, h, ...), prefixes
// (yocto..yotta), CGS / US-customary unit definitions.
//
// Temperature is the one special case: conversion between °C, °F, K
// is *affine* (offsets, not just multipliers). Handled via dedicated
// to_kelvin() / from_kelvin() helpers — multiplicative convert() is
// only correct for temperature *differences*, not absolute values.
//
// Reference: sympy/physics/units/{quantities,dimensions,prefixes,
//            unitsystem,util}.py and sympy/physics/units/definitions/.

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp::physics {

// ---------------------------------------------------------------------------
// Dimension
// ---------------------------------------------------------------------------

// Indices into the Dimension exponent tuple, in this fixed order.
enum class BaseDim : std::uint8_t {
    Length = 0,       // L  meter
    Mass = 1,         // M  kilogram
    Time = 2,         // T  second
    Current = 3,      // I  ampere
    Temperature = 4,  // Θ  kelvin
    Amount = 5,       // N  mole
    Luminosity = 6,   // J  candela
};
constexpr std::size_t kBaseDimCount = 7;

class SYMPP_EXPORT Dimension {
public:
    using Exponents = std::array<int, kBaseDimCount>;

    Dimension() = default;
    explicit Dimension(Exponents exps) : exps_(exps) {}

    [[nodiscard]] static Dimension dimensionless() {
        return Dimension{};
    }
    [[nodiscard]] static Dimension base(BaseDim which) {
        Dimension d;
        d.exps_[static_cast<std::size_t>(which)] = 1;
        return d;
    }

    [[nodiscard]] const Exponents& exps() const noexcept { return exps_; }
    [[nodiscard]] int operator[](BaseDim d) const noexcept {
        return exps_[static_cast<std::size_t>(d)];
    }
    [[nodiscard]] bool is_dimensionless() const noexcept {
        for (auto e : exps_) if (e != 0) return false;
        return true;
    }
    [[nodiscard]] bool operator==(const Dimension& o) const noexcept {
        return exps_ == o.exps_;
    }
    [[nodiscard]] bool operator!=(const Dimension& o) const noexcept {
        return !(*this == o);
    }

    // Multiplication / division of dimensions: add / subtract exponents.
    [[nodiscard]] Dimension operator*(const Dimension& o) const noexcept {
        Dimension r;
        for (std::size_t i = 0; i < kBaseDimCount; ++i) {
            r.exps_[i] = exps_[i] + o.exps_[i];
        }
        return r;
    }
    [[nodiscard]] Dimension operator/(const Dimension& o) const noexcept {
        Dimension r;
        for (std::size_t i = 0; i < kBaseDimCount; ++i) {
            r.exps_[i] = exps_[i] - o.exps_[i];
        }
        return r;
    }
    [[nodiscard]] Dimension pow(int n) const noexcept {
        Dimension r;
        for (std::size_t i = 0; i < kBaseDimCount; ++i) {
            r.exps_[i] = exps_[i] * n;
        }
        return r;
    }

    [[nodiscard]] std::string str() const;

private:
    Exponents exps_ = {0, 0, 0, 0, 0, 0, 0};
};

// ---------------------------------------------------------------------------
// Unit
// ---------------------------------------------------------------------------

class SYMPP_EXPORT Unit {
public:
    Unit(std::string name, std::string symbol, Dimension dim, Expr scale);

    [[nodiscard]] const std::string& name() const noexcept { return name_; }
    [[nodiscard]] const std::string& symbol() const noexcept { return symbol_; }
    [[nodiscard]] const Dimension& dimension() const noexcept { return dim_; }
    [[nodiscard]] const Expr& scale() const noexcept { return scale_; }

    // Compose: combined scale = a.scale * b.scale; combined dim = a.dim + b.dim.
    [[nodiscard]] Unit operator*(const Unit& o) const;
    [[nodiscard]] Unit operator/(const Unit& o) const;
    [[nodiscard]] Unit pow(int n) const;

    [[nodiscard]] bool same_dim(const Unit& o) const noexcept {
        return dim_ == o.dim_;
    }

    [[nodiscard]] std::string str() const { return symbol_; }

private:
    std::string name_;
    std::string symbol_;
    Dimension dim_;
    Expr scale_;
};

// ---------------------------------------------------------------------------
// Quantity
// ---------------------------------------------------------------------------

class SYMPP_EXPORT Quantity {
public:
    Quantity(Expr value, Unit unit);

    [[nodiscard]] const Expr& value() const noexcept { return value_; }
    [[nodiscard]] const Unit& unit() const noexcept { return unit_; }

    // Numeric value at base SI: value * unit.scale.
    [[nodiscard]] Expr base_value() const;

    // Arithmetic. Add/sub require same dimension; otherwise throws.
    [[nodiscard]] Quantity operator+(const Quantity& o) const;
    [[nodiscard]] Quantity operator-(const Quantity& o) const;
    [[nodiscard]] Quantity operator*(const Quantity& o) const;
    [[nodiscard]] Quantity operator/(const Quantity& o) const;
    [[nodiscard]] Quantity operator*(const Expr& s) const;
    [[nodiscard]] Quantity operator/(const Expr& s) const;
    [[nodiscard]] Quantity pow(int n) const;

    [[nodiscard]] std::string str() const;

private:
    Expr value_;
    Unit unit_;
};

// ---------------------------------------------------------------------------
// convert
// ---------------------------------------------------------------------------

// Convert q to target unit. Throws if dimensions don't match.
[[nodiscard]] SYMPP_EXPORT Quantity convert(const Quantity& q,
                                              const Unit& target);

// Affine temperature conversions. Use these for absolute temperatures.
// (For temperature *differences*, multiplicative convert() is correct.)
[[nodiscard]] SYMPP_EXPORT Expr celsius_to_kelvin(const Expr& c);
[[nodiscard]] SYMPP_EXPORT Expr kelvin_to_celsius(const Expr& k);
[[nodiscard]] SYMPP_EXPORT Expr fahrenheit_to_kelvin(const Expr& f);
[[nodiscard]] SYMPP_EXPORT Expr kelvin_to_fahrenheit(const Expr& k);
[[nodiscard]] SYMPP_EXPORT Expr celsius_to_fahrenheit(const Expr& c);
[[nodiscard]] SYMPP_EXPORT Expr fahrenheit_to_celsius(const Expr& f);

// ---------------------------------------------------------------------------
// SI base, derived, prefixes
// ---------------------------------------------------------------------------

// Base SI units.
[[nodiscard]] SYMPP_EXPORT const Unit& meter();
[[nodiscard]] SYMPP_EXPORT const Unit& kilogram();
[[nodiscard]] SYMPP_EXPORT const Unit& second();
[[nodiscard]] SYMPP_EXPORT const Unit& ampere();
[[nodiscard]] SYMPP_EXPORT const Unit& kelvin();
[[nodiscard]] SYMPP_EXPORT const Unit& mole();
[[nodiscard]] SYMPP_EXPORT const Unit& candela();

// Derived SI.
[[nodiscard]] SYMPP_EXPORT const Unit& hertz();      // 1/s
[[nodiscard]] SYMPP_EXPORT const Unit& newton();     // kg·m/s²
[[nodiscard]] SYMPP_EXPORT const Unit& joule();      // N·m
[[nodiscard]] SYMPP_EXPORT const Unit& watt();       // J/s
[[nodiscard]] SYMPP_EXPORT const Unit& pascal();     // N/m²
[[nodiscard]] SYMPP_EXPORT const Unit& coulomb();    // A·s
[[nodiscard]] SYMPP_EXPORT const Unit& volt();       // W/A
[[nodiscard]] SYMPP_EXPORT const Unit& ohm();        // V/A
[[nodiscard]] SYMPP_EXPORT const Unit& farad();      // C/V
[[nodiscard]] SYMPP_EXPORT const Unit& tesla();      // V·s/m²
[[nodiscard]] SYMPP_EXPORT const Unit& weber();      // V·s
[[nodiscard]] SYMPP_EXPORT const Unit& henry();      // V·s/A
[[nodiscard]] SYMPP_EXPORT const Unit& lux();        // cd·sr/m²
[[nodiscard]] SYMPP_EXPORT const Unit& gram();       // 1/1000 kg

// Prefix factory: unit with scale × prefix multiplier.
//   prefix_unit(milli, gram) → milligram with scale 10^(-6) (vs SI base kg).
[[nodiscard]] SYMPP_EXPORT Unit prefix_unit(const Expr& multiplier,
                                              std::string_view prefix_symbol,
                                              const Unit& base);

// Standard prefix multipliers (constexpr literals for handy reuse).
[[nodiscard]] SYMPP_EXPORT Expr kilo();
[[nodiscard]] SYMPP_EXPORT Expr milli();
[[nodiscard]] SYMPP_EXPORT Expr micro();
[[nodiscard]] SYMPP_EXPORT Expr nano();
[[nodiscard]] SYMPP_EXPORT Expr pico();
[[nodiscard]] SYMPP_EXPORT Expr mega();
[[nodiscard]] SYMPP_EXPORT Expr giga();
[[nodiscard]] SYMPP_EXPORT Expr tera();

// ---------------------------------------------------------------------------
// CGS + US customary
// ---------------------------------------------------------------------------

// CGS.
[[nodiscard]] SYMPP_EXPORT const Unit& centimeter();
[[nodiscard]] SYMPP_EXPORT const Unit& cgs_gram();   // alias for gram in CGS context
[[nodiscard]] SYMPP_EXPORT const Unit& dyne();       // 10^-5 N
[[nodiscard]] SYMPP_EXPORT const Unit& erg();        // 10^-7 J

// US customary.
[[nodiscard]] SYMPP_EXPORT const Unit& foot();       // 0.3048 m
[[nodiscard]] SYMPP_EXPORT const Unit& inch();       // 0.0254 m
[[nodiscard]] SYMPP_EXPORT const Unit& mile();       // 1609.344 m
[[nodiscard]] SYMPP_EXPORT const Unit& yard();       // 0.9144 m
[[nodiscard]] SYMPP_EXPORT const Unit& pound();      // 0.45359237 kg (mass)
[[nodiscard]] SYMPP_EXPORT const Unit& ounce();      // 1/16 lb
[[nodiscard]] SYMPP_EXPORT const Unit& gallon_us();  // 3.785411784 L → m³
[[nodiscard]] SYMPP_EXPORT const Unit& btu();        // 1055.06 J

// Time helpers.
[[nodiscard]] SYMPP_EXPORT const Unit& minute();
[[nodiscard]] SYMPP_EXPORT const Unit& hour();
[[nodiscard]] SYMPP_EXPORT const Unit& day();
[[nodiscard]] SYMPP_EXPORT const Unit& year();

// ---------------------------------------------------------------------------
// Physical constants (numeric values in base SI)
// ---------------------------------------------------------------------------

[[nodiscard]] SYMPP_EXPORT Quantity speed_of_light();        // c
[[nodiscard]] SYMPP_EXPORT Quantity gravitational_constant(); // G
[[nodiscard]] SYMPP_EXPORT Quantity planck();                // h
[[nodiscard]] SYMPP_EXPORT Quantity hbar();                  // ℏ = h/(2π)
[[nodiscard]] SYMPP_EXPORT Quantity boltzmann();             // k_B
[[nodiscard]] SYMPP_EXPORT Quantity elementary_charge();     // e
[[nodiscard]] SYMPP_EXPORT Quantity avogadro();              // N_A
[[nodiscard]] SYMPP_EXPORT Quantity gas_constant();          // R
[[nodiscard]] SYMPP_EXPORT Quantity electron_mass();         // m_e
[[nodiscard]] SYMPP_EXPORT Quantity proton_mass();           // m_p
[[nodiscard]] SYMPP_EXPORT Quantity vacuum_permittivity();   // ε_0
[[nodiscard]] SYMPP_EXPORT Quantity vacuum_permeability();   // μ_0

}  // namespace sympp::physics
