#pragma once

#include <cstdint>

namespace sympp {

// Stable enum tag identifying each concrete Basic subclass.
// Used for fast type-switch dispatch on hot paths and for stable canonical
// ordering of args in commutative operations.
//
// Enumerator values are assigned to roughly group related node families;
// the exact value is not stable across major versions but is stable within
// a release.
enum class TypeId : std::uint16_t {
    Unknown = 0,

    // -- Atoms --
    Symbol = 10,
    Dummy,
    Wild,

    // -- Numbers --
    Integer = 20,
    Rational,
    Float,
    Complex,
    NumberSymbol,  // Pi, E, EulerGamma, Catalan, GoldenRatio, ...
    Infinity,      // +oo
    NegativeInfinity,
    ComplexInfinity,  // zoo (1/0)
    NaN,
    ImaginaryUnit,  // I

    // -- Arithmetic --
    Add = 40,
    Mul,
    Pow,

    // -- Functions --
    Function = 60,
    UndefinedFunction,
    ElementaryFunction,
    SpecialFunction,

    // -- Calculus / unevaluated forms --
    Derivative = 80,
    Integral,
    Sum,
    Product,
    Limit,
    Order,
    RootOf,        // RootOf(poly, var, k) — k-th root of poly

    // -- Containers / control flow --
    Piecewise = 100,
    Tuple,

    // -- Boolean / set --
    Boolean = 120,
    Relational,
    Set,

    // -- Matrix --
    Matrix = 140,
    MatrixSymbol,
    MatrixExpr,
};

}  // namespace sympp
