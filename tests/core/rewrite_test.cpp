// rewrite(target) — converts between function families; each rewrite must be
// mathematically equal to its input (verified by SymPy via the oracle).

#include <nlohmann/json.hpp>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/rewrite.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/functions/hyperbolic.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
// A rewrite is correct iff (rewritten − original) is identically zero. SymPy's
// symbolic simplify can't always certify the complex-exponential identities, so
// verify numerically at sample points via the oracle's high-precision
// |value| < 10⁻²⁵ check.
bool value_preserved(const Expr& rewritten, const Expr& original) {
    auto& oracle = sympp::testing::Oracle::instance();
    Expr x = symbol("x");
    // Exact rational sample points keep everything symbolic (no premature
    // float-evaluation), so the oracle can certify |diff| < 10⁻²⁵ at 50 digits.
    for (auto [num, den] : {std::pair{7, 10}, std::pair{19, 10}, std::pair{-13, 10}}) {
        Expr diff = subs(rewritten + mul(integer(-1), original), x, rational(num, den));
        nlohmann::json req = {{"op", "evalf_is_zero"}, {"expr", diff->str()},
                              {"prec", 50}, {"tol", 25}};
        auto resp = oracle.send(req);
        if (!resp.ok || !resp.result_bool()) return false;
    }
    return true;
}
}  // namespace

TEST_CASE("rewrite to Exp preserves value (vs SymPy)", "[rewrite][oracle]") {
    auto x = symbol("x");
    for (const auto& e : {sin(x), cos(x), tan(x),
                          sinh(x), cosh(x), tanh(x),
                          sin(x) * cos(x), sin(x) + cos(x),
                          exp(sin(x))}) {
        REQUIRE(value_preserved(rewrite(e, RewriteTarget::Exp), e));
    }
}

TEST_CASE("rewrite to SinCos preserves value (vs SymPy)", "[rewrite][oracle]") {
    auto x = symbol("x");
    for (const auto& e : {tan(x), cot(x), sec(x), csc(x),
                          tan(x) + sec(x)}) {
        REQUIRE(value_preserved(rewrite(e, RewriteTarget::SinCos), e));
    }
}

TEST_CASE("rewrite leaves untargeted nodes intact", "[rewrite]") {
    auto x = symbol("x");
    // log is not a rewrite target; the expression is unchanged structurally.
    auto e = log(x) + integer(2) * x;
    REQUIRE(rewrite(e, RewriteTarget::Exp) == e);
    // but a trig argument inside is still rewritten (so the result differs).
    auto g = log(sin(x));
    REQUIRE_FALSE(rewrite(g, RewriteTarget::Exp) == g);
}
