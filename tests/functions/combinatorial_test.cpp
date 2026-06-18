// factorial / binomial / gamma / loggamma tests.
//
// References:
//   sympy/functions/combinatorial/tests/test_comb_factorials.py
//   sympy/functions/special/tests/test_gamma_functions.py

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/calculus/diff.hpp>
#include <cmath>

#include <sympp/functions/combinatorial.hpp>
#include <sympp/simplify/simplify.hpp>
#include <sympp/parsing/parser.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- factorial -------------------------------------------------------------

TEST_CASE("factorial: small integer values", "[3i][factorial]") {
    REQUIRE(factorial(integer(0)) == integer(1));
    REQUIRE(factorial(integer(1)) == integer(1));
    REQUIRE(factorial(integer(5)) == integer(120));
    REQUIRE(factorial(integer(10)) == integer(3628800));
}

TEST_CASE("factorial: large integer", "[3i][factorial][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = factorial(integer(20));
    REQUIRE(oracle.equivalent(e->str(), "factorial(20)"));
    // 20! = 2432902008176640000
    REQUIRE(e == integer("2432902008176640000"));
}

TEST_CASE("factorial: stays unevaluated on a symbol", "[3i][factorial]") {
    auto x = symbol("x");
    REQUIRE(factorial(x)->type_id() == TypeId::Function);
}

// FACT-NEGINT-1: (−n)! = Γ(−n+1) has a pole, so a negative-integer factorial is
// ComplexInfinity (matching SymPy). A non-integer negative arg keeps its node.
TEST_CASE("factorial: negative integers are zoo (FACT-NEGINT-1)",
          "[3i][factorial][regression]") {
    REQUIRE(factorial(integer(-1)) == S::ComplexInfinity());
    REQUIRE(factorial(integer(-2)) == S::ComplexInfinity());
    REQUIRE(factorial(integer(-10)) == S::ComplexInfinity());
    // 1/(−1)! = 0 (the reciprocal of the pole).
    REQUIRE(pow(factorial(integer(-1)), integer(-1)) == S::Zero());
    // A non-integer (no pole) stays a Factorial node.
    REQUIRE(factorial(rational(-1, 2))->type_id() == TypeId::Function);
}

// ----- binomial --------------------------------------------------------------

TEST_CASE("binomial: classic values", "[3i][binomial]") {
    REQUIRE(binomial(integer(5), integer(0)) == integer(1));
    REQUIRE(binomial(integer(5), integer(1)) == integer(5));
    REQUIRE(binomial(integer(5), integer(2)) == integer(10));
    REQUIRE(binomial(integer(5), integer(3)) == integer(10));
    REQUIRE(binomial(integer(5), integer(5)) == integer(1));
}

TEST_CASE("binomial: k > n returns 0", "[3i][binomial]") {
    REQUIRE(binomial(integer(3), integer(5)) == integer(0));
}

TEST_CASE("binomial: stays symbolic for symbolic args", "[3i][binomial]") {
    auto n = symbol("n");
    auto k = symbol("k");
    auto e = binomial(n, k);
    REQUIRE(e->type_id() == TypeId::Function);
    REQUIRE(e->str() == "binomial(n, k)");
}

TEST_CASE("binomial(n, 0) = 1", "[3i][binomial]") {
    auto n = symbol("n");
    REQUIRE(binomial(n, S::Zero()) == S::One());
}

// Regression (BINOM-1): binomial(n, 1) = n for any n (SymPy auto-evaluates it;
// SymPP previously kept it symbolic).
TEST_CASE("binomial(n, 1) = n", "[3i][binomial][regression]") {
    auto n = symbol("n");
    REQUIRE(binomial(n, S::One()) == n);
    REQUIRE(binomial(integer(5), S::One()) == integer(5));
    // binomial(n, 2) is not a special identity — must stay symbolic.
    REQUIRE(binomial(n, integer(2))->type_id() == TypeId::Function);
}

// ----- gamma -----------------------------------------------------------------

TEST_CASE("gamma: positive integer reduces to factorial", "[3i][gamma]") {
    // gamma(n) = (n-1)!
    REQUIRE(gamma(integer(1)) == integer(1));
    REQUIRE(gamma(integer(2)) == integer(1));
    REQUIRE(gamma(integer(3)) == integer(2));
    REQUIRE(gamma(integer(4)) == integer(6));
    REQUIRE(gamma(integer(5)) == integer(24));
}

// GAMMA-1: half-integer arguments stayed symbolic; SymPy gives closed forms
// with sqrt(pi). gamma(p/2) = C·sqrt(pi) via the duplication recurrence.
TEST_CASE("gamma: half-integer reduces to a multiple of sqrt(pi)",
          "[3i][gamma][oracle][regression]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(gamma(rational(1, 2))->str(), "sqrt(pi)"));
    REQUIRE(oracle.equivalent(gamma(rational(3, 2))->str(), "sqrt(pi)/2"));
    REQUIRE(oracle.equivalent(gamma(rational(5, 2))->str(), "3*sqrt(pi)/4"));
    REQUIRE(oracle.equivalent(gamma(rational(7, 2))->str(), "15*sqrt(pi)/8"));
}

TEST_CASE("gamma: negative half-integer reduces to a multiple of sqrt(pi)",
          "[3i][gamma][oracle][regression]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(gamma(rational(-1, 2))->str(), "-2*sqrt(pi)"));
    REQUIRE(oracle.equivalent(gamma(rational(-3, 2))->str(), "4*sqrt(pi)/3"));
}

TEST_CASE("gamma: numeric Float matches MPFR", "[3i][gamma]") {
    auto e = gamma(float_value(1.5));
    REQUIRE(e->type_id() == TypeId::Float);
}

TEST_CASE("gamma: high-precision Float matches SymPy",
          "[3i][gamma][oracle]") {
    auto& oracle = Oracle::instance();
    auto e = gamma(float_value("1.5", 30));
    auto sympy_value = oracle.evalf("gamma(Rational(3, 2))", 30);
    REQUIRE(oracle.equivalent(e->str(), sympy_value));
}

TEST_CASE("gamma: stays symbolic for symbolic arg", "[3i][gamma]") {
    auto x = symbol("x");
    REQUIRE(gamma(x)->str() == "gamma(x)");
}

// SPECVAL-1: gamma has a simple pole at every non-positive integer → zoo; and
// the polygamma special values at x = 1 (ψ⁽⁰⁾(1) = −γ,
// ψ⁽ⁿ⁾(1) = (−1)^(n+1) n! ζ(n+1)). Both match SymPy.
TEST_CASE("gamma: non-positive integers are poles → zoo (SPECVAL-1)",
          "[3i][gamma][regression]") {
    REQUIRE(gamma(integer(0)) == S::ComplexInfinity());
    REQUIRE(gamma(integer(-1)) == S::ComplexInfinity());
    REQUIRE(gamma(integer(-3)) == S::ComplexInfinity());
    // Positive integers / half-integers are unaffected.
    REQUIRE(gamma(integer(5)) == integer(24));
}

TEST_CASE("polygamma: special values at x = 1 (SPECVAL-1)",
          "[3i][polygamma][oracle][regression]") {
    auto& oracle = Oracle::instance();
    // ψ⁽⁰⁾(1) = −EulerGamma (also reachable via digamma(1)).
    REQUIRE(oracle.equivalent(polygamma(integer(0), integer(1))->str(),
                              "-EulerGamma"));
    REQUIRE(oracle.equivalent(digamma(integer(1))->str(), "-EulerGamma"));
    // ψ⁽¹⁾(1) = π²/6, ψ⁽²⁾(1) = −2ζ(3), ψ⁽³⁾(1) = π⁴/15.
    REQUIRE(oracle.equivalent(polygamma(integer(1), integer(1))->str(),
                              "pi**2/6"));
    REQUIRE(oracle.equivalent(polygamma(integer(2), integer(1))->str(),
                              "-2*zeta(3)"));
    REQUIRE(oracle.equivalent(polygamma(integer(3), integer(1))->str(),
                              "pi**4/15"));
    // A non-unit argument stays symbolic (the x = 1 rule must not over-fire).
    REQUIRE(polygamma(integer(1), symbol("x"))->str() == "polygamma(1, x)");
}

// SPECVAL-2: polygamma closed forms at any positive integer or half-integer
// argument, via the recurrence ψ⁽ⁿ⁾(y+1) = ψ⁽ⁿ⁾(y) + (−1)ⁿ·n!/y^(n+1) from the
// base values at 1 and 1/2. Previously only x = 1 was evaluated. Integer / digamma
// cases are oracle-checked against SymPy; the half-integer order-≥1 cases (which
// SymPy leaves unevaluated) are confirmed by the oracle's numeric .equals(0).
TEST_CASE("polygamma/digamma: integer and half-integer values (SPECVAL-2)",
          "[3i][polygamma][digamma][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto P = [](long n, const Expr& x) { return polygamma(integer(n), x); };
    auto half = [](long num) { return rational(num, 2); };

    // digamma (order 0) at integers: ψ(m) = −γ + H_{m−1}.
    REQUIRE(oracle.equivalent(digamma(integer(2))->str(), "1 - EulerGamma"));
    REQUIRE(oracle.equivalent(digamma(integer(3))->str(), "3/2 - EulerGamma"));
    REQUIRE(oracle.equivalent(digamma(integer(7))->str(),
                              "49/20 - EulerGamma"));
    // digamma at half-integers: ψ(1/2) = −γ − 2 log 2, ψ(m+1/2) adds 2·Σ 1/(2k−1).
    REQUIRE(oracle.equivalent(digamma(half(1))->str(),
                              "-2*log(2) - EulerGamma"));
    REQUIRE(oracle.equivalent(digamma(half(3))->str(),
                              "2 - 2*log(2) - EulerGamma"));
    REQUIRE(oracle.equivalent(digamma(half(5))->str(),
                              "8/3 - 2*log(2) - EulerGamma"));
    // Trigamma / higher order at integers.
    REQUIRE(oracle.equivalent(P(1, integer(2))->str(), "pi**2/6 - 1"));
    REQUIRE(oracle.equivalent(P(2, integer(2))->str(), "2 - 2*zeta(3)"));
    // Half-integer, order ≥ 1 (SymPy leaves these unevaluated; values exact).
    REQUIRE(oracle.equivalent(P(1, half(1))->str(), "pi**2/2"));
    REQUIRE(oracle.equivalent(P(1, half(3))->str(), "pi**2/2 - 4"));
    REQUIRE(oracle.equivalent(P(2, half(1))->str(), "-14*zeta(3)"));

    // Generic and pole arguments are untouched.
    REQUIRE(P(0, symbol("y"))->str() == "polygamma(0, y)");
    REQUIRE(digamma(integer(0)) == S::ComplexInfinity());
}

// POLYGAMMA-POLE-1: ψ⁽ⁿ⁾(x) = zoo at the nonpositive integers x ∈ {0, −1, −2, …}
// for any non-negative integer order n (the Γ pole); digamma inherits it via
// polygamma(0, ·). Matches SymPy.
TEST_CASE("polygamma/digamma: pole at nonpositive integers (POLYGAMMA-POLE-1)",
          "[3i][polygamma][digamma][regression]") {
    REQUIRE(polygamma(integer(0), integer(0)) == S::ComplexInfinity());
    REQUIRE(polygamma(integer(0), integer(-1)) == S::ComplexInfinity());
    REQUIRE(polygamma(integer(1), integer(0)) == S::ComplexInfinity());
    REQUIRE(polygamma(integer(2), integer(-3)) == S::ComplexInfinity());
    REQUIRE(digamma(integer(0)) == S::ComplexInfinity());
    REQUIRE(digamma(integer(-5)) == S::ComplexInfinity());
    // No over-reach: the pole rule must not fire for positive arguments. Positive
    // integers and half-integers now evaluate to finite closed forms (SPECVAL-2),
    // not zoo; a symbol stays symbolic.
    REQUIRE_FALSE(polygamma(integer(0), integer(2)) == S::ComplexInfinity());
    REQUIRE_FALSE(polygamma(integer(0), rational(1, 2))
                  == S::ComplexInfinity());
    REQUIRE(polygamma(integer(0), symbol("x"))->type_id() == TypeId::Function);
}

// ----- loggamma --------------------------------------------------------------

TEST_CASE("loggamma: classic values", "[3i][loggamma]") {
    REQUIRE(loggamma(integer(1)) == integer(0));   // log(0!) = 0
    REQUIRE(loggamma(integer(2)) == integer(0));   // log(1!) = 0
}

// LOGGAMMA-VALUES-1: loggamma(x) = log(Γ(x)) for x > 0 — log((n−1)!) at a positive
// integer, log(√π·…) at a positive half-integer; +∞ at the nonpositive-integer
// poles and at +∞. For x < 0 it stays loggamma(x) (branch cuts). Matches SymPy.
TEST_CASE("loggamma: positive-argument and pole values (LOGGAMMA-VALUES-1)",
          "[3i][loggamma][oracle][regression]") {
    auto& oracle = Oracle::instance();
    // Poles at the nonpositive integers and at +∞ → +∞.
    REQUIRE(loggamma(integer(0)) == S::Infinity());
    REQUIRE(loggamma(integer(-1)) == S::Infinity());
    REQUIRE(loggamma(integer(-5)) == S::Infinity());
    REQUIRE(loggamma(S::Infinity()) == S::Infinity());
    // Positive integers → log((n−1)!).
    REQUIRE(oracle.equivalent(loggamma(integer(3))->str(), "log(2)"));
    REQUIRE(oracle.equivalent(loggamma(integer(5))->str(), "log(24)"));
    // Positive half-integers → log(√π·…).
    REQUIRE(oracle.equivalent(loggamma(rational(1, 2))->str(), "log(sqrt(pi))"));
    REQUIRE(oracle.equivalent(loggamma(rational(3, 2))->str(),
                              "log(sqrt(pi)/2)"));
    // Negative (non-integer) and symbolic arguments stay symbolic.
    REQUIRE(loggamma(rational(-1, 2))->type_id() == TypeId::Function);
    REQUIRE(loggamma(rational(-3, 2))->type_id() == TypeId::Function);
    REQUIRE(loggamma(symbol("x"))->type_id() == TypeId::Function);
}

TEST_CASE("loggamma: numeric Float matches mpfr_lngamma", "[3i][loggamma]") {
    auto e = loggamma(float_value(5.0));
    REQUIRE(e->type_id() == TypeId::Float);
}

// Γ', (logΓ)' and ψ⁽ⁿ⁾' previously fell through to the default diff_arg = 0
// (DIGAMMA-1). They now produce polygamma, matching SymPy.
TEST_CASE("gamma/loggamma/polygamma derivatives (DIGAMMA-1)",
          "[3i][gamma][diff][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(diff(gamma(x), x)->str(),
                              "gamma(x)*polygamma(0, x)"));
    REQUIRE(oracle.equivalent(diff(loggamma(x), x)->str(), "polygamma(0, x)"));
    // The polygamma order increments under differentiation; the chain doesn't
    // collapse to 0.
    REQUIRE(oracle.equivalent(diff(polygamma(integer(0), x), x)->str(),
                              "polygamma(1, x)"));
    REQUIRE(oracle.equivalent(diff(polygamma(integer(1), x), x)->str(),
                              "polygamma(2, x)"));
    // Chain rule through a non-trivial argument.
    REQUIRE(oracle.equivalent(diff(gamma(pow(x, integer(2))), x)->str(),
                              "2*x*gamma(x**2)*polygamma(0, x**2)"));
    // digamma is sugar for polygamma(0, ·); parser round-trips polygamma.
    REQUIRE(digamma(x) == polygamma(integer(0), x));
    REQUIRE(parsing::parse("polygamma(0, x)") == polygamma(integer(0), x));
}

// ----- Substitution ----------------------------------------------------------

TEST_CASE("factorial: subs collapses", "[3i][subs]") {
    auto x = symbol("x");
    REQUIRE(subs(factorial(x), x, integer(6)) == integer(720));
}

TEST_CASE("binomial: subs collapses", "[3i][subs]") {
    auto n = symbol("n");
    REQUIRE(subs(binomial(n, integer(2)), n, integer(10)) == integer(45));
}

// ----- Oracle structural -----------------------------------------------------

TEST_CASE("combinatorial: structural output matches SymPy",
          "[3i][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(factorial(integer(5))->str(), "factorial(5)"));
    REQUIRE(oracle.equivalent(binomial(integer(10), integer(3))->str(),
                              "binomial(10, 3)"));
    REQUIRE(oracle.equivalent(gamma(integer(5))->str(), "gamma(5)"));
}

// ----- fibonacci / catalan (FIB-CAT) -----------------------------------------

TEST_CASE("fibonacci: integer values", "[3i][fibonacci]") {
    REQUIRE(fibonacci(integer(0)) == integer(0));
    REQUIRE(fibonacci(integer(1)) == integer(1));
    REQUIRE(fibonacci(integer(2)) == integer(1));
    REQUIRE(fibonacci(integer(10)) == integer(55));
    REQUIRE(fibonacci(integer(20)) == integer(6765));
    // Symbolic / negative arguments stay unevaluated.
    auto x = symbol("x");
    REQUIRE(fibonacci(x)->type_id() == TypeId::Function);
    REQUIRE(fibonacci(integer(-1))->type_id() == TypeId::Function);
}

// LUCAS-1: Lucas numbers L(0)=2, L(1)=1, L(n)=L(n-1)+L(n-2). Evaluates for a
// non-negative integer; symbolic and negative arguments stay symbolic. New
// function, matching SymPy's lucas.
TEST_CASE("lucas: integer values (LUCAS-1)", "[3i][lucas][regression]") {
    REQUIRE(lucas(integer(0)) == integer(2));
    REQUIRE(lucas(integer(1)) == integer(1));
    REQUIRE(lucas(integer(2)) == integer(3));
    REQUIRE(lucas(integer(5)) == integer(11));
    REQUIRE(lucas(integer(10)) == integer(123));
    REQUIRE(lucas(integer(20)) == integer(15127));
    // Symbolic / negative arguments stay unevaluated; parses round-trip.
    auto x = symbol("x");
    REQUIRE(lucas(x)->type_id() == TypeId::Function);
    REQUIRE(lucas(x)->str() == "lucas(x)");
    REQUIRE(parsing::parse("lucas(7)") == lucas(integer(7)));
    REQUIRE(lucas(integer(-3))->type_id() == TypeId::Function);
}

// TOTIENT-1: Euler's totient φ(n) evaluates for positive integers (φ(p)=p−1
// for prime p, φ(p^k·…) via the product formula), stays symbolic for symbols
// and non-positive integers. Matches SymPy's totient.
TEST_CASE("totient: Euler's phi (TOTIENT-1)", "[3i][totient][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(totient(integer(1)) == integer(1));
    REQUIRE(totient(integer(2)) == integer(1));
    REQUIRE(totient(integer(7)) == integer(6));    // prime: p − 1
    REQUIRE(totient(integer(12)) == integer(4));   // 2²·3
    REQUIRE(totient(integer(36)) == integer(12));  // 2²·3²
    REQUIRE(totient(integer(100)) == integer(40));
    REQUIRE(totient(integer(17)) == integer(16));
    // Cross-check a larger composite against SymPy.
    REQUIRE(oracle.equivalent(totient(integer(360360))->str(), "69120"));
    // Symbolic and non-positive arguments stay unevaluated.
    auto n = symbol("n");
    REQUIRE(totient(n)->type_id() == TypeId::Function);
    REQUIRE(totient(integer(0))->type_id() == TypeId::Function);
    REQUIRE(totient(integer(-5))->type_id() == TypeId::Function);
}

// PRIME-PRIMEPI-1: prime(n) is the n-th prime and primepi(n) counts primes ≤ n.
// Both evaluate for integers (primepi clamps to 0 below 2) and stay symbolic for
// symbols / non-positive prime() indices. Matches SymPy.
TEST_CASE("prime/primepi: nth prime and prime counting (PRIME-PRIMEPI-1)",
          "[3i][prime][primepi][oracle]") {
    auto& oracle = Oracle::instance();
    // prime(n): 1-indexed n-th prime.
    REQUIRE(prime(integer(1)) == integer(2));
    REQUIRE(prime(integer(5)) == integer(11));
    REQUIRE(prime(integer(10)) == integer(29));
    REQUIRE(prime(integer(100)) == integer(541));
    REQUIRE(oracle.equivalent(prime(integer(1000))->str(), "7919"));
    // primepi(n): count of primes ≤ n.
    REQUIRE(primepi(integer(1)) == integer(0));
    REQUIRE(primepi(integer(2)) == integer(1));
    REQUIRE(primepi(integer(10)) == integer(4));
    REQUIRE(primepi(integer(100)) == integer(25));
    REQUIRE(primepi(integer(-3)) == integer(0));
    REQUIRE(oracle.equivalent(primepi(integer(10000))->str(), "1229"));
    // prime(primepi(p)) == p for a prime p; primepi(prime(k)) == k.
    REQUIRE(prime(primepi(integer(13))) == integer(13));
    REQUIRE(primepi(prime(integer(7))) == integer(7));
    // Symbolic / out-of-domain arguments stay unevaluated.
    auto n = symbol("n");
    REQUIRE(prime(n)->type_id() == TypeId::Function);
    REQUIRE(primepi(n)->type_id() == TypeId::Function);
    REQUIRE(prime(integer(0))->type_id() == TypeId::Function);
}

// ARITH-FN-1: the multiplicative arithmetic functions computed from the prime
// factorization — mobius μ, divisor_count σ₀, divisor_sigma σ₁. All evaluate for
// positive integers and stay symbolic for symbols / non-positive arguments.
// Matches SymPy.
TEST_CASE("mobius/divisor_count/divisor_sigma (ARITH-FN-1)",
          "[3i][mobius][divisor][oracle]") {
    auto& oracle = Oracle::instance();
    // Möbius μ(n): 0 on a squared factor, else (−1)^#primes.
    REQUIRE(mobius(integer(1)) == integer(1));
    REQUIRE(mobius(integer(7)) == integer(-1));    // prime
    REQUIRE(mobius(integer(30)) == integer(-1));   // 2·3·5, three primes
    REQUIRE(mobius(integer(12)) == integer(0));    // 2²·3 — squared factor
    REQUIRE(mobius(integer(210)) == integer(1));   // 2·3·5·7, four primes
    // divisor_count σ₀(n) = ∏(eᵢ+1).
    REQUIRE(divisor_count(integer(1)) == integer(1));
    REQUIRE(divisor_count(integer(7)) == integer(2));
    REQUIRE(divisor_count(integer(12)) == integer(6));
    REQUIRE(divisor_count(integer(36)) == integer(9));
    // divisor_sigma σ₁(n) = sum of divisors.
    REQUIRE(divisor_sigma(integer(1)) == integer(1));
    REQUIRE(divisor_sigma(integer(6)) == integer(12));   // perfect: σ = 2n
    REQUIRE(divisor_sigma(integer(12)) == integer(28));
    REQUIRE(divisor_sigma(integer(28)) == integer(56));  // perfect
    // Cross-check a larger value against SymPy.
    REQUIRE(oracle.equivalent(divisor_sigma(integer(720))->str(), "2418"));
    // DIVSIGMA-GEN-1: generalized σ_k(n) = Σ_{d|n} d^k.
    REQUIRE(divisor_sigma(integer(6), integer(2)) == integer(50));   // 1+4+9+36
    REQUIRE(divisor_sigma(integer(10), integer(2)) == integer(130));
    REQUIRE(divisor_sigma(integer(28), integer(3)) == integer(25112));
    REQUIRE(divisor_sigma(integer(12), integer(0)) == integer(6));   // = divisor count
    REQUIRE(divisor_sigma(integer(6), integer(1))
            == divisor_sigma(integer(6)));                           // σ₁ = 1-arg
    // Symbolic / non-positive arguments stay unevaluated.
    auto n = symbol("n");
    REQUIRE(mobius(n)->type_id() == TypeId::Function);
    REQUIRE(divisor_count(integer(0))->type_id() == TypeId::Function);
    REQUIRE(divisor_sigma(integer(-4))->type_id() == TypeId::Function);
    REQUIRE(divisor_sigma(n, integer(2))->type_id() == TypeId::Function);
}

// HARMONIC-FACT2-1: harmonic(n) = Σ 1/k (a rational) and factorial2(n) = n!!
// (double factorial). Both evaluate for integers (factorial2 has the empty-
// product conventions factorial2(0)=factorial2(−1)=1) and stay symbolic for
// symbols / out-of-domain arguments. Matches SymPy.
TEST_CASE("harmonic/factorial2 (HARMONIC-FACT2-1)",
          "[3i][harmonic][factorial2][oracle]") {
    auto& oracle = Oracle::instance();
    // harmonic Hₙ.
    REQUIRE(harmonic(integer(0)) == integer(0));
    REQUIRE(harmonic(integer(1)) == integer(1));
    REQUIRE(oracle.equivalent(harmonic(integer(5))->str(), "137/60"));
    REQUIRE(oracle.equivalent(harmonic(integer(10))->str(), "7381/2520"));
    // HARMONIC-GEN-1: generalized harmonic Hₙ⁽ᵐ⁾ = Σ_{k=1}^n k^(−m).
    REQUIRE(oracle.equivalent(harmonic(integer(5), integer(2))->str(),
                              "5269/3600"));
    REQUIRE(oracle.equivalent(harmonic(integer(4), integer(3))->str(),
                              "2035/1728"));
    REQUIRE(harmonic(integer(5), integer(1)) == harmonic(integer(5)));  // m=1
    REQUIRE(harmonic(integer(3), integer(0)) == integer(3));            // Σ 1
    REQUIRE(harmonic(integer(3), integer(-1)) == integer(6));           // Σ k
    // Symbolic 2-arg stays unevaluated.
    REQUIRE(harmonic(symbol("n"), integer(2))->type_id() == TypeId::Function);
    // factorial2 n!!.
    REQUIRE(factorial2(integer(0)) == integer(1));
    REQUIRE(factorial2(integer(-1)) == integer(1));
    REQUIRE(factorial2(integer(5)) == integer(15));   // 5·3·1
    REQUIRE(factorial2(integer(6)) == integer(48));   // 6·4·2
    REQUIRE(factorial2(integer(7)) == integer(105));  // 7·5·3·1
    // Symbolic / out-of-domain arguments stay unevaluated.
    auto n = symbol("n");
    REQUIRE(harmonic(n)->type_id() == TypeId::Function);
    REQUIRE(harmonic(integer(-2))->type_id() == TypeId::Function);
    REQUIRE(factorial2(n)->type_id() == TypeId::Function);
    REQUIRE(factorial2(integer(-3))->type_id() == TypeId::Function);
}

// BERNOULLI-EULER-1: the Bernoulli numbers Bₙ (SymPy convention B₁ = +1/2) and
// Euler numbers Eₙ, computed from their binomial recurrences. Both evaluate for
// non-negative integers (odd Bₙ>1 and odd Eₙ are 0) and stay symbolic otherwise.
// Matches SymPy.
TEST_CASE("bernoulli/euler numbers (BERNOULLI-EULER-1)",
          "[3i][bernoulli][euler][oracle]") {
    auto& oracle = Oracle::instance();
    // Bernoulli: B₀=1, B₁=1/2, B₂=1/6, odd>1 vanish, B₄=−1/30, …
    REQUIRE(bernoulli(integer(0)) == integer(1));
    REQUIRE(oracle.equivalent(bernoulli(integer(1))->str(), "1/2"));
    REQUIRE(oracle.equivalent(bernoulli(integer(2))->str(), "1/6"));
    REQUIRE(bernoulli(integer(3)) == integer(0));
    REQUIRE(oracle.equivalent(bernoulli(integer(4))->str(), "-1/30"));
    REQUIRE(oracle.equivalent(bernoulli(integer(6))->str(), "1/42"));
    REQUIRE(oracle.equivalent(bernoulli(integer(12))->str(), "-691/2730"));
    // Euler: E₀=1, odd vanish, E₂=−1, E₄=5, E₆=−61, …
    REQUIRE(euler(integer(0)) == integer(1));
    REQUIRE(euler(integer(1)) == integer(0));
    REQUIRE(euler(integer(2)) == integer(-1));
    REQUIRE(euler(integer(4)) == integer(5));
    REQUIRE(euler(integer(6)) == integer(-61));
    REQUIRE(euler(integer(10)) == integer(-50521));
    // Symbolic / negative arguments stay unevaluated.
    auto n = symbol("n");
    REQUIRE(bernoulli(n)->type_id() == TypeId::Function);
    REQUIRE(euler(n)->type_id() == TypeId::Function);
    REQUIRE(bernoulli(integer(-1))->type_id() == TypeId::Function);
}

TEST_CASE("catalan: integer values", "[3i][catalan]") {
    REQUIRE(catalan(integer(0)) == integer(1));
    REQUIRE(catalan(integer(1)) == integer(1));
    REQUIRE(catalan(integer(2)) == integer(2));
    REQUIRE(catalan(integer(3)) == integer(5));
    REQUIRE(catalan(integer(5)) == integer(42));
    REQUIRE(catalan(integer(10)) == integer(16796));
    auto x = symbol("x");
    REQUIRE(catalan(x)->type_id() == TypeId::Function);
}

TEST_CASE("fibonacci/catalan: parse round-trip and subs",
          "[3i][fibonacci][catalan][parser]") {
    auto x = symbol("x");
    REQUIRE(parsing::parse("fibonacci(x)") == fibonacci(x));
    REQUIRE(parsing::parse("catalan(x)") == catalan(x));
    REQUIRE(fibonacci(x)->str() == "fibonacci(x)");
    REQUIRE(catalan(x)->str() == "catalan(x)");
    REQUIRE(subs(fibonacci(x), x, integer(12)) == integer(144));
    REQUIRE(subs(catalan(x), x, integer(4)) == integer(14));
}

TEST_CASE("fibonacci/catalan: values match SymPy", "[3i][fibonacci][catalan][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(oracle.equivalent(fibonacci(integer(15))->str(), "fibonacci(15)"));
    REQUIRE(oracle.equivalent(catalan(integer(7))->str(), "catalan(7)"));
}

// ----- gcd / lcm (GCD-LCM) ---------------------------------------------------

TEST_CASE("gcd: integer values", "[3i][gcd]") {
    REQUIRE(gcd(integer(12), integer(18)) == integer(6));
    REQUIRE(gcd(integer(17), integer(5)) == integer(1));
    REQUIRE(gcd(integer(-12), integer(8)) == integer(4));  // non-negative result
    REQUIRE(gcd(integer(0), integer(5)) == integer(5));
    REQUIRE(gcd(integer(0), integer(0)) == integer(0));
    // Symbolic (multivariate) args stay unevaluated.
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(gcd(x, y)->type_id() == TypeId::Function);
}

// GCD-POLY-1: univariate polynomial GCD. SymPy's convention is the primitive
// integer gcd (integer coefficients, content 1, positive leading) scaled by the
// gcd of the integer contents: gcd(x²−1, x−1) = x−1, gcd(2x²−2, 2x−2) = 2x−2,
// gcd(6x²+11x+3, 2x²−x−6) = 2x+3. Previously these stayed an unevaluated node.
TEST_CASE("gcd: univariate polynomial GCD (GCD-POLY-1)",
          "[3i][gcd][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto sq = [&](const Expr& e) { return pow(e, integer(2)); };
    auto chk = [&](const Expr& g, const std::string& want) {
        REQUIRE(oracle.equivalent(g->str(), want));
    };
    chk(gcd(sq(x) - integer(1), x - integer(1)), "x - 1");
    chk(gcd(sq(x) - integer(1), sq(x) + integer(2) * x + integer(1)), "x + 1");
    chk(gcd(pow(x, integer(3)) - integer(1), sq(x) - integer(1)), "x - 1");
    // Content is preserved.
    chk(gcd(integer(2) * sq(x) - integer(2), integer(2) * x - integer(2)),
        "2*x - 2");
    // Primitive integer form (not monic): gcd(6x²+11x+3, 2x²−x−6) = 2x+3.
    chk(gcd(integer(6) * sq(x) + integer(11) * x + integer(3),
            integer(2) * sq(x) - x - integer(6)),
        "2*x + 3");
    // Coprime polynomials → 1; constant vs polynomial → 1.
    chk(gcd(sq(x) + integer(1), x - integer(1)), "1");
    chk(gcd(sq(x) - integer(1), integer(2)), "1");
}

// LCM-POLY-1: univariate polynomial LCM via lcm(a,b) = a·b / gcd(a,b),
// reusing the polynomial gcd. lcm(x²−1, x−1) = x²−1, lcm(2x−2, 3x−3) = 6x−6.
TEST_CASE("lcm: univariate polynomial LCM (LCM-POLY-1)",
          "[3i][lcm][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto sq = [&](const Expr& e) { return pow(e, integer(2)); };
    auto chk = [&](const Expr& g, const std::string& want) {
        REQUIRE(oracle.equivalent(g->str(), want));
    };
    chk(lcm(sq(x) - integer(1), x - integer(1)), "x**2 - 1");
    chk(lcm(x - integer(1), x + integer(1)), "x**2 - 1");
    chk(lcm(sq(x) - integer(1), sq(x) + integer(2) * x + integer(1)),
        "x**3 + x**2 - x - 1");
    chk(lcm(integer(2) * x - integer(2), integer(3) * x - integer(3)),
        "6*x - 6");
    chk(lcm(x, sq(x)), "x**2");
    // lcm · gcd = a · b (up to content): consistency on a coprime pair.
    chk(lcm(sq(x) + integer(1), x - integer(1)),
        "x**3 - x**2 + x - 1");
}

TEST_CASE("lcm: integer values", "[3i][lcm]") {
    REQUIRE(lcm(integer(4), integer(6)) == integer(12));
    REQUIRE(lcm(integer(21), integer(6)) == integer(42));
    REQUIRE(lcm(integer(0), integer(5)) == integer(0));
    REQUIRE(lcm(integer(7), integer(5)) == integer(35));
    auto x = symbol("x");
    // lcm(x, n) = n·x over ℚ[x] (x and a constant are coprime), matching SymPy.
    REQUIRE(lcm(x, integer(4)) == integer(4) * x);
}

TEST_CASE("gcd/lcm: parse round-trip and subs", "[3i][gcd][lcm][parser]") {
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(parsing::parse("gcd(x, y)") == gcd(x, y));
    REQUIRE(parsing::parse("lcm(x, y)") == lcm(x, y));
    REQUIRE(gcd(x, y)->str() == "gcd(x, y)");
    REQUIRE(lcm(x, y)->str() == "lcm(x, y)");
    // gcd(x, n) = 1 over ℚ[x] (x and a constant are coprime), matching SymPy.
    REQUIRE(gcd(x, integer(18)) == integer(1));
    // subs reaches inside and re-evaluates: gcd(x, y)|_{y=x} = gcd(x, x) = x,
    // and lcm(x, y)|_{y=x} = lcm(x, x) = x.
    REQUIRE(subs(gcd(x, y), y, x) == x);
    REQUIRE(subs(lcm(x, y), y, x) == x);
}

// ----- rising/falling factorial & subfactorial (RFF-SUBF) --------------------

TEST_CASE("rising_factorial: products and symbolic expansion", "[3i][rising]") {
    auto x = symbol("x");
    REQUIRE(rising_factorial(integer(3), integer(2)) == integer(12));  // 3*4
    REQUIRE(rising_factorial(integer(2), integer(3)) == integer(24));  // 2*3*4
    REQUIRE(rising_factorial(x, integer(0)) == integer(1));
    REQUIRE(rising_factorial(x, integer(2)) == parsing::parse("x*(x + 1)"));
    // Symbolic n stays unevaluated.
    REQUIRE(rising_factorial(x, symbol("n"))->type_id() == TypeId::Function);
}

TEST_CASE("falling_factorial: products and symbolic expansion", "[3i][falling]") {
    auto x = symbol("x");
    REQUIRE(falling_factorial(integer(5), integer(2)) == integer(20));  // 5*4
    REQUIRE(falling_factorial(integer(7), integer(3)) == integer(210)); // 7*6*5
    REQUIRE(falling_factorial(x, integer(0)) == integer(1));
    REQUIRE(falling_factorial(x, integer(2)) == parsing::parse("x*(x - 1)"));
}

TEST_CASE("subfactorial: derangement counts", "[3i][subfactorial]") {
    REQUIRE(subfactorial(integer(0)) == integer(1));
    REQUIRE(subfactorial(integer(1)) == integer(0));
    REQUIRE(subfactorial(integer(2)) == integer(1));
    REQUIRE(subfactorial(integer(4)) == integer(9));
    REQUIRE(subfactorial(integer(5)) == integer(44));
    auto x = symbol("x");
    REQUIRE(subfactorial(x)->type_id() == TypeId::Function);
}

TEST_CASE("rising/falling/subfactorial: SymPy agreement & round-trip",
          "[3i][rising][falling][subfactorial][oracle][parser]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(rising_factorial(integer(3), integer(4))->str(),
                              "RisingFactorial(3, 4)"));
    REQUIRE(oracle.equivalent(falling_factorial(integer(8), integer(3))->str(),
                              "FallingFactorial(8, 3)"));
    REQUIRE(oracle.equivalent(subfactorial(integer(6))->str(), "subfactorial(6)"));
    REQUIRE(parsing::parse("subfactorial(x)") == subfactorial(x));
    REQUIRE(parsing::parse("RisingFactorial(x, 3)") == rising_factorial(x, integer(3)));
}

// ----- Beta function (BETA) --------------------------------------------------
// B(a,b)=Γ(a)Γ(b)/Γ(a+b); folds to a closed value when the gammas evaluate
// (positive integers, half-integers), else stays symbolic.
TEST_CASE("beta: evaluates via gamma ratio", "[3i][beta][oracle]") {
    auto& oracle = Oracle::instance();
    REQUIRE(beta(integer(2), integer(3)) == rational(1, 12));   // 1!*2!/4!
    REQUIRE(beta(integer(5), integer(2)) == rational(1, 30));
    REQUIRE(beta(integer(1), integer(1)) == integer(1));
    REQUIRE(beta(integer(3), integer(4)) == rational(1, 60));
    // Half-integers: B(1/2,1/2) = Γ(1/2)²/Γ(1) = π.
    REQUIRE(oracle.equivalent(beta(rational(1, 2), rational(1, 2))->str(), "pi"));
    // Symmetry and SymPy agreement.
    REQUIRE(beta(integer(2), integer(5)) == beta(integer(5), integer(2)));
}

TEST_CASE("beta: symbolic args stay unevaluated", "[3i][beta][parser]") {
    auto a = symbol("a");
    auto b = symbol("b");
    REQUIRE(beta(a, b)->type_id() == TypeId::Function);
    REQUIRE(parsing::parse("beta(a, b)") == beta(a, b));
    REQUIRE(beta(a, b)->str() == "beta(a, b)");
}

// FUNC-INCGAMMA-1: lowergamma/uppergamma as real function classes — the
// positive-integer first argument collapses to the closed elementary form, the
// derivatives are exact, and symbolic orders stay unevaluated and round-trip.
TEST_CASE("incomplete gamma: closed forms, derivatives, round-trip (FUNC-INCGAMMA-1)",
          "[3i][incompletegamma][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto s = symbol("s");

    // Closed forms for a positive-integer order (match SymPy).
    REQUIRE(oracle.equivalent(uppergamma(integer(1), x)->str(), "exp(-x)"));
    REQUIRE(oracle.equivalent(lowergamma(integer(1), x)->str(), "1 - exp(-x)"));
    REQUIRE(oracle.equivalent(uppergamma(integer(2), x)->str(), "(x + 1)*exp(-x)"));
    REQUIRE(oracle.equivalent(lowergamma(integer(2), x)->str(),
                              "1 - (x + 1)*exp(-x)"));
    REQUIRE(oracle.equivalent(uppergamma(integer(3), x)->str(),
                              "2*(x**2/2 + x + 1)*exp(-x)"));

    // γ(s,x) + Γ(s,x) = Γ(s); here Γ(3) = 2.
    REQUIRE(simplify(add(uppergamma(integer(3), x), lowergamma(integer(3), x)))
            == integer(2));

    // Special points.
    REQUIRE(uppergamma(integer(2), integer(0)) == integer(1));  // Γ(2,0)=Γ(2)
    REQUIRE(lowergamma(integer(2), integer(0)) == S::Zero());   // γ(s,0)=0
    REQUIRE(uppergamma(integer(3), S::Infinity()) == S::Zero());  // Γ(s,∞)=0

    // Derivatives: d/dx γ(s,x)= xˢ⁻¹e⁻ˣ, d/dx Γ(s,x)= −xˢ⁻¹e⁻ˣ.
    REQUIRE(oracle.equivalent(diff(lowergamma(s, x), x)->str(),
                              "x**(s-1)*exp(-x)"));
    REQUIRE(oracle.equivalent(diff(uppergamma(s, x), x)->str(),
                              "-x**(s-1)*exp(-x)"));

    // Symbolic order stays unevaluated and round-trips through the parser.
    REQUIRE(uppergamma(s, x)->type_id() == TypeId::Function);
    REQUIRE(lowergamma(s, x)->type_id() == TypeId::Function);
    REQUIRE(parsing::parse("uppergamma(s, x)") == uppergamma(s, x));
    REQUIRE(parsing::parse("lowergamma(s, x)") == lowergamma(s, x));
    REQUIRE(uppergamma(s, x)->str() == "uppergamma(s, x)");

    // Real for real order and argument.
    auto sr = symbol("sr", AssumptionMask{}.set_real(true));
    auto xr = symbol("xr", AssumptionMask{}.set_real(true));
    REQUIRE(is_real(uppergamma(sr, xr)) == true);
    REQUIRE(is_real(lowergamma(sr, xr)) == true);
}

// FUNC-INCGAMMA-HALF-1: half-integer orders reduce to erf/erfc closed forms via
// the recurrence V(s+1) = s·V(s) ± x^s·e^{-x}, up and down from the base
// Γ(1/2,x)=√π·erfc(√x), γ(1/2,x)=√π·erf(√x). Matches SymPy.
TEST_CASE("incomplete gamma: half-integer erf/erfc forms (FUNC-INCGAMMA-HALF-1)",
          "[3i][incompletegamma][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");

    // Base order 1/2.
    REQUIRE(oracle.equivalent(uppergamma(rational(1, 2), x)->str(),
                              "sqrt(pi)*erfc(sqrt(x))"));
    REQUIRE(oracle.equivalent(lowergamma(rational(1, 2), x)->str(),
                              "sqrt(pi)*erf(sqrt(x))"));
    // Climb: 3/2 and 5/2.
    REQUIRE(oracle.equivalent(uppergamma(rational(3, 2), x)->str(),
                              "sqrt(pi)*erfc(sqrt(x))/2 + sqrt(x)*exp(-x)"));
    REQUIRE(oracle.equivalent(lowergamma(rational(3, 2), x)->str(),
                              "sqrt(pi)*erf(sqrt(x))/2 - sqrt(x)*exp(-x)"));
    REQUIRE(oracle.equivalent(
        uppergamma(rational(5, 2), x)->str(),
        "3*sqrt(pi)*erfc(sqrt(x))/4 + x**(3/2)*exp(-x) + 3*sqrt(x)*exp(-x)/2"));
    // Descend: −1/2.
    REQUIRE(oracle.equivalent(uppergamma(rational(-1, 2), x)->str(),
                              "-2*sqrt(pi)*erfc(sqrt(x)) + 2*exp(-x)/sqrt(x)"));

    // γ + Γ = Γ(s); Γ(3/2) = √π/2. (SymPy's simplify won't fold erf+erfc=1, so
    // verify the identity numerically at a sample point instead.)
    Expr sum = add(uppergamma(rational(3, 2), x), lowergamma(rational(3, 2), x));
    Expr diff = subs(sum, x, integer(2)) - gamma(rational(3, 2));
    REQUIRE(std::fabs(std::stod(evalf(diff, 30)->str())) < 1e-15);

    // Γ(1/2, 0) = Γ(1/2) = √π.
    REQUIRE(oracle.equivalent(uppergamma(rational(1, 2), integer(0))->str(),
                              "sqrt(pi)"));
}
