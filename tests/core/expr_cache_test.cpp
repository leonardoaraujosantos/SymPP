// ExprCache (hash-consing) tests.
//
// The contract: any two structurally-equal Exprs constructed via make<T> /
// the public factory functions share a single shared_ptr instance. This
// strengthens equality from O(tree size) to O(1) on the common path.
//
// Reference: SymPy's @cacheit decorator on construction
//            (sympy/core/cache.py) and SymEngine's intern semantics.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/add.hpp>
#include <sympp/core/expr_cache.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>

using namespace sympp;

TEST_CASE("cache: two Symbol('x') share one shared_ptr", "[1g][cache]") {
    auto a = symbol("x");
    auto b = symbol("x");
    REQUIRE(a.get() == b.get());
}

TEST_CASE("cache: distinct names give distinct instances", "[1g][cache]") {
    auto a = symbol("x");
    auto b = symbol("y");
    REQUIRE(a.get() != b.get());
}

TEST_CASE("cache: Integer values intern by value", "[1g][cache]") {
    auto a = integer(42);
    auto b = integer(42);
    REQUIRE(a.get() == b.get());

    auto big1 = integer("123456789012345678901234567890");
    auto big2 = integer("123456789012345678901234567890");
    REQUIRE(big1.get() == big2.get());
}

TEST_CASE("cache: Rational values intern by reduced form", "[1g][cache]") {
    auto a = rational(1, 2);
    auto b = rational(2, 4);  // reduces to 1/2
    REQUIRE(a.get() == b.get());
}

TEST_CASE("cache: Add interns identical sums", "[1g][cache]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto a = x + y;
    auto b = x + y;
    REQUIRE(a.get() == b.get());

    // Different argument order produces the same canonical form, also shared.
    auto c = y + x;
    REQUIRE(a.get() == c.get());
}

TEST_CASE("cache: Mul interns identical products", "[1g][cache]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto a = x * y;
    auto b = y * x;
    REQUIRE(a.get() == b.get());
}

TEST_CASE("cache: Pow interns identical powers", "[1g][cache]") {
    auto x = symbol("x");
    auto a = pow(x, integer(2));
    auto b = pow(x, integer(2));
    REQUIRE(a.get() == b.get());
}

TEST_CASE("cache: deep tree shares all subexpressions", "[1g][cache]") {
    auto x = symbol("x");
    auto y = symbol("y");

    auto e1 = pow(x, integer(2)) + integer(3) * x + y;
    auto e2 = pow(x, integer(2)) + integer(3) * x + y;
    REQUIRE(e1.get() == e2.get());

    // Also: structural sub-equality preserved through expressions.
    auto p = pow(x, integer(2));
    auto fragment = e1->args()[0];  // first arg of Add — depending on order
    // (We don't know which arg is first without inspecting; just confirm at
    // least one sub-arg points at the cached Pow.)
    bool found = false;
    for (const auto& a : e1->args()) {
        if (a.get() == p.get()) {
            found = true;
            break;
        }
    }
    REQUIRE(found);
}

TEST_CASE("cache: singletons remain stable", "[1g][cache]") {
    REQUIRE(S::Zero().get() == S::Zero().get());
    REQUIRE(S::One().get() == S::One().get());
    REQUIRE(S::Half().get() == S::Half().get());
    REQUIRE(S::Pi().get() == S::Pi().get());

    // Constructing the same Number via the public factory hits the same
    // cache entry as the singleton.
    REQUIRE(integer(0).get() == S::Zero().get());
    REQUIRE(integer(1).get() == S::One().get());
    REQUIRE(integer(-1).get() == S::NegativeOne().get());
    REQUIRE(rational(1, 2).get() == S::Half().get());
}

TEST_CASE("cache: lazy GC reclaims expired entries", "[1g][cache]") {
    auto& cache = ExprCache::instance();
    auto baseline = cache.size();

    // Build ~50 unique short-lived Symbols and let them go out of scope.
    {
        for (int i = 0; i < 50; ++i) {
            auto s = symbol("__cache_test_sym_" + std::to_string(i));
            (void)s;
        }
    }
    // The cache currently holds 50 expired weak_ptrs. cleanup() drops them.
    auto removed = cache.cleanup();
    REQUIRE(removed >= 50);
    // After cleanup, size should be at or below the baseline.
    REQUIRE(cache.size() <= baseline);
}

TEST_CASE("cache: disable / re-enable", "[1g][cache]") {
    auto& cache = ExprCache::instance();
    REQUIRE(cache.enabled());

    cache.set_enabled(false);
    auto a = symbol("__cache_disabled_sym");
    auto b = symbol("__cache_disabled_sym");
    // With caching off, freshly-constructed Exprs are NOT pointer-equal,
    // but they're still structurally equal.
    REQUIRE(a.get() != b.get());
    REQUIRE(a == b);  // structural equality unchanged

    cache.set_enabled(true);
    auto c = symbol("__cache_disabled_sym");
    auto d = symbol("__cache_disabled_sym");
    REQUIRE(c.get() == d.get());  // back to interning
}
