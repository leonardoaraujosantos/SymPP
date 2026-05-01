#include <sympp/physics/units.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp::physics {

// ---------------------------------------------------------------------------
// Dimension::str
// ---------------------------------------------------------------------------

std::string Dimension::str() const {
    static const char* labels[kBaseDimCount] = {
        "L", "M", "T", "I", "Theta", "N", "J"
    };
    if (is_dimensionless()) return "1";
    std::string out;
    for (std::size_t i = 0; i < kBaseDimCount; ++i) {
        if (exps_[i] == 0) continue;
        if (!out.empty()) out += "*";
        out += labels[i];
        if (exps_[i] != 1) {
            out += "^" + std::to_string(exps_[i]);
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Unit
// ---------------------------------------------------------------------------

Unit::Unit(std::string name, std::string symbol, Dimension dim, Expr scale)
    : name_(std::move(name)), symbol_(std::move(symbol)),
      dim_(dim), scale_(std::move(scale)) {}

Unit Unit::operator*(const Unit& o) const {
    return Unit(name_ + "*" + o.name_, symbol_ + "*" + o.symbol_,
                dim_ * o.dim_, simplify(mul(scale_, o.scale_)));
}

Unit Unit::operator/(const Unit& o) const {
    return Unit(name_ + "/" + o.name_, symbol_ + "/" + o.symbol_,
                dim_ / o.dim_, simplify(scale_ / o.scale_));
}

Unit Unit::pow(int n) const {
    return Unit(name_ + "^" + std::to_string(n),
                symbol_ + "^" + std::to_string(n),
                dim_.pow(n), simplify(sympp::pow(scale_, integer(n))));
}

// ---------------------------------------------------------------------------
// Quantity
// ---------------------------------------------------------------------------

Quantity::Quantity(Expr value, Unit unit)
    : value_(std::move(value)), unit_(std::move(unit)) {}

Expr Quantity::base_value() const {
    return simplify(mul(value_, unit_.scale()));
}

Quantity Quantity::operator+(const Quantity& o) const {
    if (!unit_.same_dim(o.unit_)) {
        throw std::invalid_argument(
            "Quantity::+: incompatible dimensions");
    }
    // Convert o into our unit, then add values.
    Expr o_in_us = simplify(mul(o.value_, o.unit_.scale()) / unit_.scale());
    return Quantity(simplify(add(value_, o_in_us)), unit_);
}

Quantity Quantity::operator-(const Quantity& o) const {
    if (!unit_.same_dim(o.unit_)) {
        throw std::invalid_argument(
            "Quantity::-: incompatible dimensions");
    }
    Expr o_in_us = simplify(mul(o.value_, o.unit_.scale()) / unit_.scale());
    return Quantity(simplify(value_ - o_in_us), unit_);
}

Quantity Quantity::operator*(const Quantity& o) const {
    return Quantity(simplify(mul(value_, o.value_)), unit_ * o.unit_);
}

Quantity Quantity::operator/(const Quantity& o) const {
    return Quantity(simplify(value_ / o.value_), unit_ / o.unit_);
}

Quantity Quantity::operator*(const Expr& s) const {
    return Quantity(simplify(mul(value_, s)), unit_);
}

Quantity Quantity::operator/(const Expr& s) const {
    return Quantity(simplify(value_ / s), unit_);
}

Quantity Quantity::pow(int n) const {
    return Quantity(simplify(sympp::pow(value_, integer(n))), unit_.pow(n));
}

std::string Quantity::str() const {
    return value_->str() + " " + unit_.str();
}

// ---------------------------------------------------------------------------
// convert
// ---------------------------------------------------------------------------

Quantity convert(const Quantity& q, const Unit& target) {
    if (q.unit().dimension() != target.dimension()) {
        throw std::invalid_argument(
            "convert: source and target unit dimensions differ");
    }
    // q.value · q.unit.scale = target.value · target.scale
    // → target.value = q.value · q.unit.scale / target.scale
    Expr new_value = simplify(mul(q.value(), q.unit().scale())
                                / target.scale());
    return Quantity(new_value, target);
}

// ---------------------------------------------------------------------------
// Temperature affine conversions
// ---------------------------------------------------------------------------

Expr celsius_to_kelvin(const Expr& c) {
    return simplify(c + rational(54630, 200));   // 273.15 = 54630/200
}
Expr kelvin_to_celsius(const Expr& k) {
    return simplify(k - rational(54630, 200));
}
Expr fahrenheit_to_kelvin(const Expr& f) {
    // (F - 32) * 5/9 + 273.15
    return simplify((f - integer(32)) * rational(5, 9) + rational(54630, 200));
}
Expr kelvin_to_fahrenheit(const Expr& k) {
    return simplify((k - rational(54630, 200)) * rational(9, 5) + integer(32));
}
Expr celsius_to_fahrenheit(const Expr& c) {
    return simplify(c * rational(9, 5) + integer(32));
}
Expr fahrenheit_to_celsius(const Expr& f) {
    return simplify((f - integer(32)) * rational(5, 9));
}

// ---------------------------------------------------------------------------
// SI base units
// ---------------------------------------------------------------------------

const Unit& meter() {
    static const Unit u("meter", "m",
                          Dimension::base(BaseDim::Length), integer(1));
    return u;
}
const Unit& kilogram() {
    static const Unit u("kilogram", "kg",
                          Dimension::base(BaseDim::Mass), integer(1));
    return u;
}
const Unit& second() {
    static const Unit u("second", "s",
                          Dimension::base(BaseDim::Time), integer(1));
    return u;
}
const Unit& ampere() {
    static const Unit u("ampere", "A",
                          Dimension::base(BaseDim::Current), integer(1));
    return u;
}
const Unit& kelvin() {
    static const Unit u("kelvin", "K",
                          Dimension::base(BaseDim::Temperature), integer(1));
    return u;
}
const Unit& mole() {
    static const Unit u("mole", "mol",
                          Dimension::base(BaseDim::Amount), integer(1));
    return u;
}
const Unit& candela() {
    static const Unit u("candela", "cd",
                          Dimension::base(BaseDim::Luminosity), integer(1));
    return u;
}

// ---------------------------------------------------------------------------
// Derived SI units
// ---------------------------------------------------------------------------

const Unit& hertz() {
    static const Unit u("hertz", "Hz",
                          Dimension::base(BaseDim::Time).pow(-1), integer(1));
    return u;
}
const Unit& newton() {
    // N = kg·m/s²
    static const Unit u(
        "newton", "N",
        Dimension::base(BaseDim::Mass) * Dimension::base(BaseDim::Length)
            / Dimension::base(BaseDim::Time).pow(2),
        integer(1));
    return u;
}
const Unit& joule() {
    // J = N·m = kg·m²/s²
    static const Unit u(
        "joule", "J",
        Dimension::base(BaseDim::Mass)
            * Dimension::base(BaseDim::Length).pow(2)
            / Dimension::base(BaseDim::Time).pow(2),
        integer(1));
    return u;
}
const Unit& watt() {
    // W = J/s = kg·m²/s³
    static const Unit u(
        "watt", "W",
        Dimension::base(BaseDim::Mass)
            * Dimension::base(BaseDim::Length).pow(2)
            / Dimension::base(BaseDim::Time).pow(3),
        integer(1));
    return u;
}
const Unit& pascal() {
    // Pa = N/m² = kg/(m·s²)
    static const Unit u(
        "pascal", "Pa",
        Dimension::base(BaseDim::Mass)
            / Dimension::base(BaseDim::Length)
            / Dimension::base(BaseDim::Time).pow(2),
        integer(1));
    return u;
}
const Unit& coulomb() {
    // C = A·s
    static const Unit u(
        "coulomb", "C",
        Dimension::base(BaseDim::Current) * Dimension::base(BaseDim::Time),
        integer(1));
    return u;
}
const Unit& volt() {
    // V = W/A = kg·m²/(A·s³)
    static const Unit u(
        "volt", "V",
        Dimension::base(BaseDim::Mass)
            * Dimension::base(BaseDim::Length).pow(2)
            / (Dimension::base(BaseDim::Current)
               * Dimension::base(BaseDim::Time).pow(3)),
        integer(1));
    return u;
}
const Unit& ohm() {
    // Ω = V/A = kg·m²/(A²·s³)
    static const Unit u(
        "ohm", "Ohm",
        Dimension::base(BaseDim::Mass)
            * Dimension::base(BaseDim::Length).pow(2)
            / (Dimension::base(BaseDim::Current).pow(2)
               * Dimension::base(BaseDim::Time).pow(3)),
        integer(1));
    return u;
}
const Unit& farad() {
    // F = C/V = A²·s⁴/(kg·m²)
    static const Unit u(
        "farad", "F",
        Dimension::base(BaseDim::Current).pow(2)
            * Dimension::base(BaseDim::Time).pow(4)
            / (Dimension::base(BaseDim::Mass)
               * Dimension::base(BaseDim::Length).pow(2)),
        integer(1));
    return u;
}
const Unit& tesla() {
    // T = V·s/m² = kg/(A·s²)
    static const Unit u(
        "tesla", "T",
        Dimension::base(BaseDim::Mass)
            / (Dimension::base(BaseDim::Current)
               * Dimension::base(BaseDim::Time).pow(2)),
        integer(1));
    return u;
}
const Unit& weber() {
    // Wb = V·s = kg·m²/(A·s²)
    static const Unit u(
        "weber", "Wb",
        Dimension::base(BaseDim::Mass)
            * Dimension::base(BaseDim::Length).pow(2)
            / (Dimension::base(BaseDim::Current)
               * Dimension::base(BaseDim::Time).pow(2)),
        integer(1));
    return u;
}
const Unit& henry() {
    // H = Wb/A
    static const Unit u(
        "henry", "H",
        Dimension::base(BaseDim::Mass)
            * Dimension::base(BaseDim::Length).pow(2)
            / (Dimension::base(BaseDim::Current).pow(2)
               * Dimension::base(BaseDim::Time).pow(2)),
        integer(1));
    return u;
}
const Unit& lux() {
    // lx = cd/m²  (steradian is dimensionless)
    static const Unit u(
        "lux", "lx",
        Dimension::base(BaseDim::Luminosity)
            / Dimension::base(BaseDim::Length).pow(2),
        integer(1));
    return u;
}
const Unit& gram() {
    // g = 1/1000 kg
    static const Unit u(
        "gram", "g",
        Dimension::base(BaseDim::Mass), rational(1, 1000));
    return u;
}

// ---------------------------------------------------------------------------
// Prefixes
// ---------------------------------------------------------------------------

Unit prefix_unit(const Expr& multiplier, std::string_view prefix_symbol,
                   const Unit& base) {
    return Unit(std::string(prefix_symbol) + base.name(),
                std::string(prefix_symbol) + base.symbol(),
                base.dimension(), simplify(mul(multiplier, base.scale())));
}

Expr kilo()  { return integer(1000); }
Expr milli() { return rational(1, 1000); }
Expr micro() { return rational(1, 1'000'000); }
Expr nano()  { return rational(1, 1'000'000'000); }
Expr pico()  { return pow(integer(10), integer(-12)); }
Expr mega()  { return integer(1'000'000); }
Expr giga()  { return integer(1'000'000'000); }
Expr tera()  { return pow(integer(10), integer(12)); }

// ---------------------------------------------------------------------------
// CGS / US customary
// ---------------------------------------------------------------------------

const Unit& centimeter() {
    static const Unit u("centimeter", "cm",
                          Dimension::base(BaseDim::Length),
                          rational(1, 100));
    return u;
}
const Unit& cgs_gram() { return gram(); }
const Unit& dyne() {
    // 1 dyne = 1 g·cm/s² = 10^-5 N. Dimension same as newton.
    static const Unit u(
        "dyne", "dyn",
        Dimension::base(BaseDim::Mass) * Dimension::base(BaseDim::Length)
            / Dimension::base(BaseDim::Time).pow(2),
        rational(1, 100'000));
    return u;
}
const Unit& erg() {
    // 1 erg = 1 g·cm²/s² = 10^-7 J
    static const Unit u(
        "erg", "erg",
        Dimension::base(BaseDim::Mass)
            * Dimension::base(BaseDim::Length).pow(2)
            / Dimension::base(BaseDim::Time).pow(2),
        pow(integer(10), integer(-7)));
    return u;
}
const Unit& foot() {
    static const Unit u("foot", "ft",
                          Dimension::base(BaseDim::Length),
                          rational(3048, 10000));
    return u;
}
const Unit& inch() {
    static const Unit u("inch", "in",
                          Dimension::base(BaseDim::Length),
                          rational(254, 10000));
    return u;
}
const Unit& mile() {
    static const Unit u("mile", "mi",
                          Dimension::base(BaseDim::Length),
                          rational(1609344, 1000));
    return u;
}
const Unit& yard() {
    static const Unit u("yard", "yd",
                          Dimension::base(BaseDim::Length),
                          rational(9144, 10000));
    return u;
}
const Unit& pound() {
    // 1 lb = 0.45359237 kg = 45359237 / 100000000 kg
    static const Unit u("pound", "lb",
                          Dimension::base(BaseDim::Mass),
                          rational(45359237LL, 100'000'000LL));
    return u;
}
const Unit& ounce() {
    // 1 oz = 1/16 lb
    static const Unit u("ounce", "oz",
                          Dimension::base(BaseDim::Mass),
                          rational(45359237LL, 1'600'000'000LL));
    return u;
}
const Unit& gallon_us() {
    // 1 US gallon = 3.785411784 L = 0.003785411784 m³
    static const Unit u("gallon_us", "gal_US",
                          Dimension::base(BaseDim::Length).pow(3),
                          rational(3785411784LL, 1'000'000'000'000LL));
    return u;
}
const Unit& btu() {
    // 1 BTU (IT) = 1055.05585262 J — round to 105505585262 / 10^8
    static const Unit u(
        "btu", "BTU",
        Dimension::base(BaseDim::Mass)
            * Dimension::base(BaseDim::Length).pow(2)
            / Dimension::base(BaseDim::Time).pow(2),
        rational(105505585262LL, 100'000'000LL));
    return u;
}

// ---------------------------------------------------------------------------
// Time helpers
// ---------------------------------------------------------------------------

const Unit& minute() {
    static const Unit u("minute", "min",
                          Dimension::base(BaseDim::Time), integer(60));
    return u;
}
const Unit& hour() {
    static const Unit u("hour", "h",
                          Dimension::base(BaseDim::Time), integer(3600));
    return u;
}
const Unit& day() {
    static const Unit u("day", "d",
                          Dimension::base(BaseDim::Time), integer(86400));
    return u;
}
const Unit& year() {
    // Julian year = 365.25 days
    static const Unit u(
        "year", "yr",
        Dimension::base(BaseDim::Time),
        mul(rational(36525, 100), integer(86400)));
    return u;
}

// ---------------------------------------------------------------------------
// Physical constants — values in base SI
// ---------------------------------------------------------------------------

namespace {
[[nodiscard]] Expr float_const(const std::string& s) {
    // 30-digit Float for high-precision constants.
    return make<Float>(s, 30);
}
}  // namespace

Quantity speed_of_light() {
    // Exact: 299792458 m/s
    return Quantity(integer(299792458), meter() / second());
}

Quantity gravitational_constant() {
    // 6.67430e-11 m³/(kg·s²)
    Unit u = (meter().pow(3) / (kilogram() * second().pow(2)));
    return Quantity(float_const("6.67430e-11"), u);
}

Quantity planck() {
    // h = 6.62607015e-34 J·s (exact since 2019 redef)
    Unit u = joule() * second();
    return Quantity(float_const("6.62607015e-34"), u);
}

Quantity hbar() {
    // ℏ = h / (2π)
    Unit u = joule() * second();
    Expr h = float_const("6.62607015e-34");
    return Quantity(simplify(h / (integer(2) * S::Pi())), u);
}

Quantity boltzmann() {
    // k_B = 1.380649e-23 J/K (exact since 2019 redef)
    Unit u = joule() / kelvin();
    return Quantity(float_const("1.380649e-23"), u);
}

Quantity elementary_charge() {
    // e = 1.602176634e-19 C (exact since 2019 redef)
    return Quantity(float_const("1.602176634e-19"), coulomb());
}

Quantity avogadro() {
    // N_A = 6.02214076e+23 1/mol (exact since 2019 redef)
    return Quantity(float_const("6.02214076e+23"),
                     mole().pow(-1));
}

Quantity gas_constant() {
    // R = N_A · k_B = 8.31446261815324 J/(mol·K)
    Unit u = joule() / (mole() * kelvin());
    return Quantity(float_const("8.31446261815324"), u);
}

Quantity electron_mass() {
    // m_e = 9.1093837015e-31 kg
    return Quantity(float_const("9.1093837015e-31"), kilogram());
}

Quantity proton_mass() {
    // m_p = 1.67262192369e-27 kg
    return Quantity(float_const("1.67262192369e-27"), kilogram());
}

Quantity vacuum_permittivity() {
    // ε_0 = 8.8541878128e-12 F/m
    Unit u = farad() / meter();
    return Quantity(float_const("8.8541878128e-12"), u);
}

Quantity vacuum_permeability() {
    // μ_0 = 1.25663706212e-6 H/m
    Unit u = henry() / meter();
    return Quantity(float_const("1.25663706212e-6"), u);
}

}  // namespace sympp::physics
