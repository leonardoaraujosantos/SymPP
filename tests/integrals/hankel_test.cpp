// Hankel transform — table pairs cross-checked against SymPy.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/integrals/transforms.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
auto& oracle() { return sympp::testing::Oracle::instance(); }
}  // namespace

TEST_CASE("hankel_transform: table pairs match SymPy", "[transforms][hankel][oracle]") {
    auto r = symbol("r"), k = symbol("k");
    auto a = symbol("a");
    auto H = [&](const Expr& f, const Expr& nu) {
        return hankel_transform(f, r, k, nu)->str();
    };
    auto& o = oracle();
    auto emar = exp(mul(mul(integer(-1), a), r));           // e^{−a·r}
    auto gauss = exp(mul(mul(integer(-1), a), pow(r, integer(2))));  // e^{−a·r²}

    // e^{−a·r}: ν=0 → a/(a²+k²)^{3/2}, ν=1 → k/(a²+k²)^{3/2}.
    REQUIRE(o.equivalent(H(emar, S::Zero()), "a/(a**2 + k**2)**(3/2)"));
    REQUIRE(o.equivalent(H(emar, integer(1)), "k/(a**2 + k**2)**(3/2)"));
    // e^{−a·r²}, ν=0 → e^{−k²/(4a)}/(2a).
    REQUIRE(o.equivalent(H(gauss, S::Zero()), "exp(-k**2/(4*a))/(2*a)"));
    // r·e^{−a·r²}, ν=1 → k·e^{−k²/(4a)}/(4a²).
    REQUIRE(o.equivalent(H(mul(r, gauss), integer(1)), "k*exp(-k**2/(4*a))/(4*a**2)"));
    // 1/r, ν=0 → 1/k.
    REQUIRE(o.equivalent(H(pow(r, integer(-1)), S::Zero()), "1/k"));

    // Linearity: H[2·e^{−r} + 1/r] (ν=0).
    auto lin = hankel_transform(add(mul(integer(2), exp(mul(integer(-1), r))),
                                    pow(r, integer(-1))), r, k, S::Zero());
    REQUIRE(o.equivalent(lin->str(), "2/(k**2 + 1)**(3/2) + 1/k"));

    // Unrecognized → unevaluated HankelTransform node.
    REQUIRE(hankel_transform(symbol("g"), r, k, S::Zero())->str().rfind("HankelTransform", 0)
            == 0);
}
