// Singleton + NumberSymbol tests.
// Reference: sympy/core/tests/test_numbers.py — pi, E, EulerGamma, Catalan
//            sympy/core/singleton.py — S registry semantics

#include <string>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number_symbol.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Number singletons ------------------------------------------------------

TEST_CASE("S::Zero is Integer(0)", "[1][singletons]") {
    auto z = S::Zero();
    REQUIRE(z->type_id() == TypeId::Integer);
    REQUIRE(z->str() == "0");
    REQUIRE(static_cast<const Integer&>(*z).is_zero());
}

TEST_CASE("S::One is Integer(1)", "[1][singletons]") {
    auto o = S::One();
    REQUIRE(o->type_id() == TypeId::Integer);
    REQUIRE(o->str() == "1");
}

TEST_CASE("S::NegativeOne is Integer(-1)", "[1][singletons]") {
    auto n = S::NegativeOne();
    REQUIRE(n->type_id() == TypeId::Integer);
    REQUIRE(n->str() == "-1");
}

TEST_CASE("S::Half is Rational(1,2)", "[1][singletons]") {
    auto h = S::Half();
    REQUIRE(h->type_id() == TypeId::Rational);
    REQUIRE(h->str() == "1/2");
    REQUIRE(*h == *rational(1, 2));
}

TEST_CASE("S::NegativeHalf is Rational(-1,2)", "[1][singletons]") {
    auto h = S::NegativeHalf();
    REQUIRE(h->type_id() == TypeId::Rational);
    REQUIRE(h->str() == "-1/2");
}

TEST_CASE("Singletons return identical Expr across calls", "[1][singletons]") {
    // Address stability — the Meyer singleton always returns the same Expr.
    REQUIRE(S::Zero().get() == S::Zero().get());
    REQUIRE(S::One().get() == S::One().get());
    REQUIRE(S::Half().get() == S::Half().get());
    REQUIRE(S::Pi().get() == S::Pi().get());
}

// ----- NumberSymbols ----------------------------------------------------------

TEST_CASE("S::Pi has the expected attributes", "[1][singletons][pi]") {
    auto pi = S::Pi();
    REQUIRE(pi->type_id() == TypeId::NumberSymbol);
    REQUIRE(pi->str() == "pi");
    auto& p = static_cast<const NumberSymbol&>(*pi);
    REQUIRE(p.kind() == NumberSymbolKind::Pi);
    REQUIRE(p.is_positive());
    REQUIRE(p.is_real());
    REQUIRE_FALSE(p.is_rational());
    REQUIRE(p.sign() == 1);
}

TEST_CASE("S::E has the expected attributes", "[1][singletons][e]") {
    auto e = S::E();
    REQUIRE(e->type_id() == TypeId::NumberSymbol);
    REQUIRE(e->str() == "E");
    auto& es = static_cast<const NumberSymbol&>(*e);
    REQUIRE(es.kind() == NumberSymbolKind::E);
    REQUIRE(es.is_positive());
}

TEST_CASE("Different NumberSymbols compare unequal", "[1][singletons]") {
    REQUIRE_FALSE(S::Pi() == S::E());
    REQUIRE_FALSE(S::Pi() == S::EulerGamma());
    REQUIRE_FALSE(S::Pi() == S::Catalan());
    REQUIRE_FALSE(S::E() == S::EulerGamma());
}

TEST_CASE("Same NumberSymbol compares equal across instances",
          "[1][singletons]") {
    Expr a = make<PiSymbol>();
    Expr b = make<PiSymbol>();
    REQUIRE(a == b);
    REQUIRE(a->hash() == b->hash());
}

// ----- evalf on NumberSymbols -------------------------------------------------

TEST_CASE("evalf(Pi, 30) matches SymPy pi.evalf(30)", "[1][singletons][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto pi_50 = evalf(S::Pi(), 30);
    REQUIRE(pi_50->type_id() == TypeId::Float);
    auto sympy_value = oracle.evalf("pi", 30);
    INFO("SymPP pi: " << pi_50->str());
    INFO("SymPy pi: " << sympy_value);
    REQUIRE(oracle.equivalent(pi_50->str(), sympy_value));
}

TEST_CASE("evalf(Pi, 100) matches SymPy at high precision",
          "[1][singletons][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto pi_hi = evalf(S::Pi(), 100);
    auto sympy_value = oracle.evalf("pi", 100);
    INFO("SymPP: " << pi_hi->str());
    INFO("SymPy: " << sympy_value);
    REQUIRE(oracle.equivalent(pi_hi->str(), sympy_value));
}

TEST_CASE("evalf(E, 50) matches SymPy E.evalf(50)", "[1][singletons][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto e_val = evalf(S::E(), 50);
    auto sympy_value = oracle.evalf("E", 50);
    INFO("SymPP E: " << e_val->str());
    INFO("SymPy E: " << sympy_value);
    REQUIRE(oracle.equivalent(e_val->str(), sympy_value));
}

TEST_CASE("evalf(EulerGamma, 40) matches SymPy", "[1][singletons][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto g = evalf(S::EulerGamma(), 40);
    auto sympy_value = oracle.evalf("EulerGamma", 40);
    INFO("SymPP: " << g->str());
    INFO("SymPy: " << sympy_value);
    REQUIRE(oracle.equivalent(g->str(), sympy_value));
}

TEST_CASE("evalf(Catalan, 40) matches SymPy", "[1][singletons][evalf][oracle]") {
    auto& oracle = Oracle::instance();
    auto k = evalf(S::Catalan(), 40);
    auto sympy_value = oracle.evalf("Catalan", 40);
    INFO("SymPP: " << k->str());
    INFO("SymPy: " << sympy_value);
    REQUIRE(oracle.equivalent(k->str(), sympy_value));
}
