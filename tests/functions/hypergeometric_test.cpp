// Hypergeometric function tests — Hyper / MeijerG factories with
// auto-evaluation, plus hyperexpand closed-form rewrites.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hypergeometric.hpp>
#include <sympp/simplify/hyperexpand.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Auto-eval at construction --------------------------------------------

TEST_CASE("hyper((), (), z) auto-evaluates to exp(z)",
          "[5][hyper][autoeval][oracle]") {
    auto& oracle = Oracle::instance();
    auto z = symbol("z");
    auto h = hyper({}, {}, z);
    REQUIRE(oracle.equivalent(h->str(), "exp(z)"));
}

TEST_CASE("hyper((a,), (), z) auto-evaluates to (1 - z)^(-a)",
          "[5][hyper][autoeval][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a");
    auto z = symbol("z");
    auto h = hyper({a}, {}, z);
    REQUIRE(oracle.equivalent(h->str(), "(1 - z)**(-a)"));
}

TEST_CASE("hyper at z = 0 is 1",
          "[5][hyper][autoeval]") {
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto h = hyper({a, b}, {c}, S::Zero());
    REQUIRE(h == S::One());
}

TEST_CASE("hyper cancels matching numerator / denominator parameters",
          "[5][hyper][autoeval][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a");
    auto b = symbol("b");
    auto z = symbol("z");
    // ₂F₁(a, b; b; z) = ₁F₀(a; ; z) = (1 − z)^(−a).
    auto h = hyper({a, b}, {b}, z);
    REQUIRE(oracle.equivalent(h->str(), "(1 - z)**(-a)"));
}

TEST_CASE("hyper(a, b, c, z) 4-arg overload routes to ₂F₁",
          "[5][hyper][autoeval]") {
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto z = symbol("z");
    auto h = hyper(a, b, c, z);
    // No auto-eval rule fires; should be the ₂F₁ Hyper node.
    REQUIRE(h->type_id() == TypeId::Function);
}

// ----- hyperexpand: ₁F₁ and ₂F₁ closed forms --------------------------------

TEST_CASE("hyperexpand: ₁F₁(1; 2; z) = (exp(z) - 1)/z",
          "[5][hyperexpand][oracle]") {
    auto& oracle = Oracle::instance();
    auto z = symbol("z");
    auto h = hyper({integer(1)}, {integer(2)}, z);
    auto out = hyperexpand(h);
    REQUIRE(oracle.equivalent(out->str(), "(exp(z) - 1)/z"));
}

TEST_CASE("hyperexpand: ₂F₁(1, 1; 2; z) = -log(1 - z)/z",
          "[5][hyperexpand][oracle]") {
    auto& oracle = Oracle::instance();
    auto z = symbol("z");
    auto h = hyper({integer(1), integer(1)}, {integer(2)}, z);
    auto out = hyperexpand(h);
    REQUIRE(oracle.equivalent(out->str(), "-log(1 - z)/z"));
}

TEST_CASE("hyperexpand: leaves unrecognized hyper alone",
          "[5][hyperexpand]") {
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto z = symbol("z");
    auto h = hyper({a, b}, {c}, z);
    auto out = hyperexpand(h);
    REQUIRE(out == h);
}

TEST_CASE("hyperexpand: walks Mul / Add",
          "[5][hyperexpand][oracle]") {
    auto& oracle = Oracle::instance();
    auto z = symbol("z");
    auto h = hyper({integer(1), integer(1)}, {integer(2)}, z);
    auto e = z * h + integer(2);
    auto out = hyperexpand(e);
    REQUIRE(oracle.equivalent(out->str(), "z*(-log(1 - z)/z) + 2"));
}

TEST_CASE("simplify chains hyperexpand: ₁F₁(1; 2; z) collapses",
          "[5][hyperexpand][simplify][oracle]") {
    auto& oracle = Oracle::instance();
    auto z = symbol("z");
    auto h = hyper({integer(1)}, {integer(2)}, z);
    auto s = simplify(h);
    REQUIRE(oracle.equivalent(s->str(), "(exp(z) - 1)/z"));
}

// ----- Cross-validation: SymPP hyperexpand vs SymPy hyperexpand -------------
//
// These tests build a hypergeometric expression in SymPP, run *both*
// SymPP's and SymPy's hyperexpand on the same input, and verify the
// two outputs are algebraically equivalent. This is a stronger check
// than the hand-coded-expected-form tests above: the reference value
// is whatever SymPy's reference implementation produces, not what we
// hand-wrote.

namespace {

// Send `op: "hyperexpand"` to the oracle and return the string result.
[[nodiscard]] std::string sympy_hyperexpand(const std::string& s) {
    auto& oracle = Oracle::instance();
    auto resp = oracle.send({{"op", "hyperexpand"}, {"expr", s}});
    REQUIRE(resp.ok);
    return resp.raw.at("result").get<std::string>();
}

}  // namespace

TEST_CASE("cross-validate: SymPP hyperexpand matches SymPy on ₀F₀",
          "[5][hyperexpand][cross][oracle]") {
    auto& oracle = Oracle::instance();
    auto z = symbol("z");
    auto h = hyper({}, {}, z);
    auto sympp_out = hyperexpand(h);
    auto sympy_out = sympy_hyperexpand("hyper((), (), z)");
    REQUIRE(oracle.equivalent(sympp_out->str(), sympy_out));
}

TEST_CASE("cross-validate: SymPP hyperexpand matches SymPy on ₁F₀",
          "[5][hyperexpand][cross][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a");
    auto z = symbol("z");
    auto h = hyper({a}, {}, z);
    auto sympp_out = hyperexpand(h);
    auto sympy_out = sympy_hyperexpand("hyper((a,), (), z)");
    REQUIRE(oracle.equivalent(sympp_out->str(), sympy_out));
}

TEST_CASE("cross-validate: SymPP hyperexpand matches SymPy on ₁F₁(1; 2; z)",
          "[5][hyperexpand][cross][oracle]") {
    auto& oracle = Oracle::instance();
    auto z = symbol("z");
    auto h = hyper({integer(1)}, {integer(2)}, z);
    auto sympp_out = hyperexpand(h);
    auto sympy_out = sympy_hyperexpand("hyper((1,), (2,), z)");
    REQUIRE(oracle.equivalent(sympp_out->str(), sympy_out));
}

TEST_CASE("cross-validate: SymPP hyperexpand matches SymPy on ₂F₁(1, 1; 2; z)",
          "[5][hyperexpand][cross][oracle]") {
    auto& oracle = Oracle::instance();
    auto z = symbol("z");
    auto h = hyper({integer(1), integer(1)}, {integer(2)}, z);
    auto sympp_out = hyperexpand(h);
    auto sympy_out = sympy_hyperexpand("hyper((1, 1), (2,), z)");
    REQUIRE(oracle.equivalent(sympp_out->str(), sympy_out));
}

TEST_CASE("cross-validate: SymPP factory cancellation matches SymPy on ₂F₁(a, b; b; z)",
          "[5][hyperexpand][cross][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a");
    auto b = symbol("b");
    auto z = symbol("z");
    auto h = hyper({a, b}, {b}, z);                // SymPP cancels at construction
    auto sympp_out = hyperexpand(h);
    auto sympy_out = sympy_hyperexpand("hyper((a, b), (b,), z)");
    REQUIRE(oracle.equivalent(sympp_out->str(), sympy_out));
}

TEST_CASE("cross-validate: hyper(a, b; c; 0) → 1 matches SymPy",
          "[5][hyperexpand][cross][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto h = hyper({a, b}, {c}, S::Zero());
    auto sympp_out = hyperexpand(h);
    auto sympy_out = sympy_hyperexpand("hyper((a, b), (c,), 0)");
    REQUIRE(oracle.equivalent(sympp_out->str(), sympy_out));
}

TEST_CASE("cross-validate: SymPP leaves unrecognized ₂F₁ alone — SymPy too",
          "[5][hyperexpand][cross][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto z = symbol("z");
    auto h = hyper({a, b}, {c}, z);
    auto sympp_out = hyperexpand(h);
    auto sympy_out = sympy_hyperexpand("hyper((a, b), (c,), z)");
    // Both should leave it as a hyper(...). Equivalence still holds.
    REQUIRE(oracle.equivalent(sympp_out->str(), sympy_out));
}

// ----- MeijerG: factory + printability --------------------------------------

TEST_CASE("meijerg: factory builds an opaque node",
          "[5][meijerg]") {
    auto a = symbol("a");
    auto b = symbol("b");
    auto z = symbol("z");
    auto g = meijerg({}, {a}, {b}, {}, z);
    REQUIRE(g->type_id() == TypeId::Function);
    // The printable form should mention "meijerg".
    REQUIRE(g->str().find("meijerg") != std::string::npos);
}

TEST_CASE("meijerg: rebuild preserves structure",
          "[5][meijerg]") {
    auto a = symbol("a");
    auto b = symbol("b");
    auto z = symbol("z");
    auto g = meijerg({}, {a}, {b}, {}, z);
    // Substituting z → 2*z should produce a new meijerg node.
    auto g2 = subs(g, z, integer(2) * z);
    REQUIRE(g2->str().find("meijerg") != std::string::npos);
}
