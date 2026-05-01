#include <sympp/core/singletons.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/number_symbol.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/fwd.hpp>

namespace sympp::S {

// Meyer-singletons: each accessor caches an Expr in function-local storage so
// the address is stable for the program's lifetime.

const Expr& Zero() {
    static const Expr value = make<Integer>(0);
    return value;
}

const Expr& One() {
    static const Expr value = make<Integer>(1);
    return value;
}

const Expr& NegativeOne() {
    static const Expr value = make<Integer>(-1);
    return value;
}

const Expr& Half() {
    // Half is a Rational, not an Integer.
    static const Expr value = make<Rational>(1L, 2L);
    return value;
}

const Expr& NegativeHalf() {
    static const Expr value = make<Rational>(-1L, 2L);
    return value;
}

const Expr& Pi() {
    static const Expr value = make<PiSymbol>();
    return value;
}

const Expr& E() {
    static const Expr value = make<ESymbol>();
    return value;
}

const Expr& EulerGamma() {
    static const Expr value = make<EulerGammaSymbol>();
    return value;
}

const Expr& Catalan() {
    static const Expr value = make<CatalanSymbol>();
    return value;
}

}  // namespace sympp::S
