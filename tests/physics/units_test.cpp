// Phase 12 — units module tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/physics/units.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using namespace sympp::physics;
using sympp::testing::Oracle;

// ----- Dimension ------------------------------------------------------------

TEST_CASE("Dimension: dimensionless is all zeros", "[12][units][dimension]") {
    REQUIRE(Dimension::dimensionless().is_dimensionless());
    REQUIRE(!Dimension::base(BaseDim::Length).is_dimensionless());
}

TEST_CASE("Dimension: multiplication adds exponents",
          "[12][units][dimension]") {
    auto L = Dimension::base(BaseDim::Length);
    auto T = Dimension::base(BaseDim::Time);
    auto velocity = L / T;
    REQUIRE(velocity[BaseDim::Length] == 1);
    REQUIRE(velocity[BaseDim::Time] == -1);
}

TEST_CASE("Dimension: pow scales exponents", "[12][units][dimension]") {
    auto L = Dimension::base(BaseDim::Length);
    auto vol = L.pow(3);
    REQUIRE(vol[BaseDim::Length] == 3);
}

TEST_CASE("Dimension: comparison", "[12][units][dimension]") {
    auto L = Dimension::base(BaseDim::Length);
    auto M = Dimension::base(BaseDim::Mass);
    REQUIRE(L == L);
    REQUIRE(L != M);
}

// ----- Unit -----------------------------------------------------------------

TEST_CASE("Unit: SI base have unit scale", "[12][units][si]") {
    REQUIRE(meter().scale() == integer(1));
    REQUIRE(kilogram().scale() == integer(1));
    REQUIRE(second().scale() == integer(1));
    REQUIRE(meter().dimension() == Dimension::base(BaseDim::Length));
}

TEST_CASE("Unit: derived dimensions match SI base products",
          "[12][units][si]") {
    // newton dim should equal kg·m/s²
    auto kg_m_s2 = kilogram().dimension() * meter().dimension()
                   / second().dimension().pow(2);
    REQUIRE(newton().dimension() == kg_m_s2);
    // joule = N·m
    REQUIRE(joule().dimension() == newton().dimension() * meter().dimension());
    // watt = J/s
    REQUIRE(watt().dimension() == joule().dimension() / second().dimension());
}

TEST_CASE("Unit: gram has scale 1/1000 vs kg", "[12][units][si]") {
    REQUIRE(gram().scale() == rational(1, 1000));
    REQUIRE(gram().dimension() == kilogram().dimension());
}

// ----- Quantity arithmetic --------------------------------------------------

TEST_CASE("Quantity: add same-dim", "[12][units][quantity]") {
    Quantity a(integer(5), meter());
    Quantity b(integer(3), meter());
    auto c = a + b;
    REQUIRE(c.value() == integer(8));
}

TEST_CASE("Quantity: add different-dim throws",
          "[12][units][quantity]") {
    Quantity a(integer(5), meter());
    Quantity b(integer(3), kilogram());
    REQUIRE_THROWS_AS(a + b, std::invalid_argument);
}

TEST_CASE("Quantity: cross-unit add converts to LHS unit",
          "[12][units][quantity]") {
    // 5 m + 50 cm = 5.5 m
    Quantity a(integer(5), meter());
    Quantity b(integer(50), centimeter());
    auto c = a + b;
    REQUIRE(c.value() == rational(11, 2));
}

TEST_CASE("Quantity: multiply combines units",
          "[12][units][quantity]") {
    Quantity a(integer(3), meter());
    Quantity b(integer(4), meter());
    auto area = a * b;
    REQUIRE(area.value() == integer(12));
    REQUIRE(area.unit().dimension() == meter().dimension().pow(2));
}

TEST_CASE("Quantity: pow", "[12][units][quantity]") {
    Quantity a(integer(2), meter());
    auto cube = a.pow(3);
    REQUIRE(cube.value() == integer(8));
    REQUIRE(cube.unit().dimension() == meter().dimension().pow(3));
}

// ----- Conversions ----------------------------------------------------------

TEST_CASE("convert: 1 km to m", "[12][units][convert]") {
    auto km = prefix_unit(kilo(), "k", meter());
    Quantity a(integer(1), km);
    auto in_m = convert(a, meter());
    REQUIRE(in_m.value() == integer(1000));
}

TEST_CASE("convert: 1 minute to seconds", "[12][units][convert]") {
    Quantity a(integer(1), minute());
    auto in_s = convert(a, second());
    REQUIRE(in_s.value() == integer(60));
}

TEST_CASE("convert: 1 mile to km via meters",
          "[12][units][convert][oracle]") {
    auto& oracle = Oracle::instance();
    auto km = prefix_unit(kilo(), "k", meter());
    Quantity a(integer(1), mile());
    auto in_km = convert(a, km);
    // 1 mile = 1.609344 km
    REQUIRE(oracle.equivalent(in_km.value()->str(), "1.609344"));
}

TEST_CASE("convert: 1 erg to joules", "[12][units][convert]") {
    Quantity a(integer(1), erg());
    auto in_J = convert(a, joule());
    // 1 erg = 10^-7 J
    REQUIRE(in_J.value() == pow(integer(10), integer(-7)));
}

TEST_CASE("convert: incompatible dimensions throws",
          "[12][units][convert]") {
    Quantity a(integer(5), meter());
    REQUIRE_THROWS_AS(convert(a, kilogram()), std::invalid_argument);
}

// ----- Temperature ---------------------------------------------------------

TEST_CASE("temperature: 0 °C → 273.15 K", "[12][units][temperature]") {
    auto k = celsius_to_kelvin(integer(0));
    // 273.15 = 54630/200 = 5463/20
    REQUIRE(k == rational(5463, 20));
}

TEST_CASE("temperature: 100 °C → 373.15 K",
          "[12][units][temperature]") {
    auto k = celsius_to_kelvin(integer(100));
    // 373.15 = 74630/200 = 7463/20
    REQUIRE(k == rational(7463, 20));
}

TEST_CASE("temperature: 32 °F → 0 °C", "[12][units][temperature]") {
    auto c = fahrenheit_to_celsius(integer(32));
    REQUIRE(c == integer(0));
}

TEST_CASE("temperature: 212 °F → 100 °C", "[12][units][temperature]") {
    auto c = fahrenheit_to_celsius(integer(212));
    REQUIRE(c == integer(100));
}

TEST_CASE("temperature: round trip K → °C → K",
          "[12][units][temperature]") {
    auto roundtrip = celsius_to_kelvin(kelvin_to_celsius(integer(300)));
    REQUIRE(roundtrip == integer(300));
}

// ----- Physical constants ---------------------------------------------------

TEST_CASE("constants: speed of light exact",
          "[12][units][constants]") {
    auto c = speed_of_light();
    REQUIRE(c.value() == integer(299792458));
    REQUIRE(c.unit().dimension() == meter().dimension() / second().dimension());
}

TEST_CASE("constants: planck has units of J·s",
          "[12][units][constants]") {
    auto h = planck();
    REQUIRE(h.unit().dimension()
            == (joule() * second()).dimension());
}

TEST_CASE("constants: avogadro per mole",
          "[12][units][constants]") {
    auto Na = avogadro();
    REQUIRE(Na.unit().dimension() == mole().dimension().pow(-1));
}

// ----- Real-world combo example --------------------------------------------

TEST_CASE("units: kinetic energy E = m·v²/2 dimensionally consistent",
          "[12][units][example]") {
    Quantity m(integer(2), kilogram());                 // 2 kg
    Quantity v(integer(3), meter() / second());          // 3 m/s
    auto E = m * v.pow(2) / Quantity(integer(2), kilogram() / kilogram());
    // E.value should be 2·9/2 = 9 kg·m²/s² = 9 J
    REQUIRE(E.unit().dimension() == joule().dimension());
}

TEST_CASE("units: Newton's law F = m·a dimensionally → newton",
          "[12][units][example]") {
    Quantity m(integer(10), kilogram());
    Quantity a(rational(98, 10), meter() / second().pow(2));  // 9.8 m/s²
    auto F = m * a;
    REQUIRE(F.unit().dimension() == newton().dimension());
    REQUIRE(F.value() == integer(98));
}
