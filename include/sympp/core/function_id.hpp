#pragma once

// FunctionId — stable tag identifying each named function. Two Function
// nodes with the same FunctionId and structurally-equal args are equal.
//
// New entries are appended to keep this stable across releases (and to
// keep the hash-cons cache valid across version bumps in the same process).

#include <cstdint>

namespace sympp {

enum class FunctionId : std::uint16_t {
    // -- User-declared (anonymous) -----------------------------------------
    Undefined = 0,

    // -- Trigonometric -----------------------------------------------------
    Sin = 100,
    Cos,
    Tan,
    Cot,
    Sec,
    Csc,
    Asin,
    Acos,
    Atan,
    Atan2,

    // -- Hyperbolic --------------------------------------------------------
    Sinh = 200,
    Cosh,
    Tanh,
    Coth,
    Sech,
    Csch,
    Asinh,
    Acosh,
    Atanh,

    // -- Exponential / log -------------------------------------------------
    Exp = 300,
    Log,
    LambertW,

    // -- Complex / sign ----------------------------------------------------
    Abs = 400,
    Sign,
    Re,
    Im,
    Conjugate,
    Arg,

    // -- Integer rounding --------------------------------------------------
    Floor = 500,
    Ceiling,
    Frac,
    Round,

    // -- Misc --------------------------------------------------------------
    Min = 600,
    Max,

    // -- Special functions (Phase 3+) --------------------------------------
    Gamma = 1000,
    LogGamma,
    DiGamma,
    Beta,
    Erf,
    Erfc,
    BesselJ,
    BesselY,
    BesselI,
    BesselK,
    Zeta,
    DiracDelta,
    Heaviside,
};

}  // namespace sympp
