// Probability distributions — moments / density / cdf cross-checked against
// sympy.stats.

#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/stats/stats.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
namespace st = sympp::stats;

namespace {
std::string stat(const std::string& dist, const std::string& fn,
                 const std::vector<std::string>& params, const std::string& x = "") {
    nlohmann::json req = {{"op", "stats"}, {"dist", dist}, {"fn", fn}, {"params", params}};
    if (!x.empty()) req["x"] = x;
    auto r = sympp::testing::Oracle::instance().send(req);
    REQUIRE(r.ok);
    return r.result_str();
}
bool equiv(const Expr& a, const std::string& sympy_str) {
    return sympp::testing::Oracle::instance().equivalent(a->str(), sympy_str);
}
}  // namespace

TEST_CASE("distribution moments match sympy.stats", "[stats][oracle]") {
    struct D { std::string name; st::Distribution dist; std::vector<std::string> ps; };
    std::vector<D> ds = {
        {"normal", st::normal(integer(0), integer(1)), {"0", "1"}},
        {"normal", st::normal(symbol("m"), symbol("s")), {"m", "s"}},
        {"uniform", st::uniform(integer(2), integer(8)), {"2", "8"}},
        {"exponential", st::exponential(integer(3)), {"3"}},
        {"bernoulli", st::bernoulli(rational(1, 3)), {"1/3"}},
        {"binomial", st::binomial(integer(10), rational(1, 4)), {"10", "1/4"}},
        {"poisson", st::poisson(integer(5)), {"5"}},
    };
    for (const auto& d : ds) {
        REQUIRE(equiv(st::mean(d.dist), stat(d.name, "E", d.ps)));
        REQUIRE(equiv(st::variance(d.dist), stat(d.name, "variance", d.ps)));
    }
}

TEST_CASE("continuous pdf/cdf match sympy.stats", "[stats][oracle]") {
    auto x = symbol("x");
    auto N = st::normal(integer(0), integer(1));
    REQUIRE(equiv(st::pdf(N, x), stat("normal", "density", {"0", "1"}, "x")));
    REQUIRE(equiv(st::cdf(N, x), stat("normal", "cdf", {"0", "1"}, "x")));

    // Exponential / Uniform densities are Piecewise in SymPy; compare at an
    // interior sample point (where the piecewise collapses to the value).
    auto E = st::exponential(integer(2));
    REQUIRE(equiv(st::pdf(E, rational(3, 2)), stat("exponential", "density", {"2"}, "3/2")));
    REQUIRE(equiv(st::cdf(E, integer(1)), stat("exponential", "cdf", {"2"}, "1")));

    auto U = st::uniform(integer(0), integer(4));
    REQUIRE(equiv(st::pdf(U, integer(2)), stat("uniform", "density", {"0", "4"}, "2")));
}

TEST_CASE("discrete pmf known values", "[stats][oracle]") {
    // Binomial(4, 1/2): P(2) = C(4,2)/16 = 3/8.
    auto B = st::binomial(integer(4), rational(1, 2));
    REQUIRE(equiv(st::pdf(B, integer(2)), "3/8"));
    // Bernoulli(1/3): P(1) = 1/3, P(0) = 2/3.
    auto Be = st::bernoulli(rational(1, 3));
    REQUIRE(equiv(st::pdf(Be, integer(1)), "1/3"));
    REQUIRE(equiv(st::pdf(Be, integer(0)), "2/3"));
}
