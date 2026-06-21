// Number-theory primitives — cross-checked against SymPy via the oracle's
// `ntheory` op.

#include <string>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/functions/ntheory.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
std::string fi_str(const Expr& n) {
    std::string s;
    bool first = true;
    for (const auto& [p, e] : factorint(n)) {
        if (!first) s += ";";
        first = false;
        s += p->str() + "^" + e->str();
    }
    return s;
}
std::string div_str(const Expr& n) {
    std::string s;
    bool first = true;
    for (const auto& d : divisors(n)) {
        if (!first) s += ",";
        first = false;
        s += d->str();
    }
    return s;
}
std::string oracle_nt(std::initializer_list<std::pair<const char*, std::string>> kv) {
    nlohmann::json req = {{"op", "ntheory"}};
    for (auto& [k, v] : kv) req[k] = v;
    auto resp = sympp::testing::Oracle::instance().send(req);
    REQUIRE(resp.ok);
    return resp.result_str();
}
}  // namespace

TEST_CASE("factorint matches SymPy", "[ntheory][oracle]") {
    for (long n : {360, 1024, 17, 1, -12, 999983, 2310}) {
        REQUIRE(fi_str(integer(n)) ==
                oracle_nt({{"fn", "factorint"}, {"n", std::to_string(n)}}));
    }
    // Large semiprime: exercises Pollard's rho past the trial-division bound.
    auto big = integer(104729) * integer(1299709);  // two primes
    REQUIRE(fi_str(big) == oracle_nt({{"fn", "factorint"}, {"n", big->str()}}));
}

TEST_CASE("divisors matches SymPy", "[ntheory][oracle]") {
    for (long n : {28, 12, 1, 100, 360, -36}) {
        REQUIRE(div_str(integer(n)) ==
                oracle_nt({{"fn", "divisors"}, {"n", std::to_string(n)}}));
    }
}

TEST_CASE("igcdex matches SymPy and the Bézout identity", "[ntheory][oracle]") {
    struct AB { long a, b; };
    for (auto [a, b] : {AB{12, 8}, AB{240, 46}, AB{7, 0}, AB{2, 3}, AB{-12, 18}}) {
        auto [x, y, g] = igcdex(integer(a), integer(b));
        // Bézout: a·x + b·y = g.
        auto lhs = integer(a) * x + integer(b) * y;
        REQUIRE(static_cast<const Integer&>(*lhs).value()
                == static_cast<const Integer&>(*g).value());
        // Same triple as SymPy.
        std::string got = x->str() + "," + y->str() + "," + g->str();
        REQUIRE(got == oracle_nt({{"fn", "igcdex"},
                                  {"a", std::to_string(a)},
                                  {"b", std::to_string(b)}}));
    }
}

TEST_CASE("jacobi_symbol matches SymPy", "[ntheory][oracle]") {
    struct AN { long a, n; };
    for (auto [a, n] : {AN{2, 15}, AN{7, 15}, AN{1001, 9907}, AN{0, 3}, AN{30, 7}}) {
        REQUIRE(jacobi_symbol(integer(a), integer(n))->str() ==
                oracle_nt({{"fn", "jacobi"},
                           {"a", std::to_string(a)},
                           {"n", std::to_string(n)}}));
    }
}

TEST_CASE("ntheory error handling", "[ntheory]") {
    REQUIRE_THROWS(factorint(integer(0)));
    REQUIRE_THROWS(jacobi_symbol(integer(3), integer(8)));   // even modulus
    REQUIRE_THROWS(jacobi_symbol(integer(3), integer(-5)));  // nonpositive
}
