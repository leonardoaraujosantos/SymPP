#pragma once

// Laplace, Fourier, Mellin, Z, and sine/cosine transforms — all
// table-based, with linearity over Add and constant-factor pull-out
// over Mul. Anything outside the table returns an UndefinedFunction
// marker named after the transform.
//
// Laplace pairs (variable s):
//     1        ->  1/s
//     t        ->  1/s^2
//     t^n      ->  n! / s^(n+1)                        (n nonneg integer)
//     exp(a*t) ->  1/(s - a)
//     sin(a*t) ->  a / (s^2 + a^2)
//     cos(a*t) ->  s / (s^2 + a^2)
//
// Inverse Laplace runs apart(F, s) and looks each summand up in the
// reverse table:
//     1/s             ->  1
//     1/s^n           ->  t^(n-1) / (n-1)!
//     1/(s - a)       ->  exp(a*t)
//     1/(s^2 + a^2)   ->  sin(a*t) / a
//     s/(s^2 + a^2)   ->  cos(a*t)
//
// Fourier (SymPy unitary "0, -2π" convention is verbose; here we use
// the engineering convention F(ω) = ∫ f(t) e^(-iωt) dt):
//     δ(t)            ->  1
//     1               ->  2π·δ(ω)
//     exp(i·a·t)      ->  2π·δ(ω - a)
//     exp(-a·t²)      ->  sqrt(π/a) · exp(-ω²/(4a))   (a > 0)
//
// Mellin (variable s):  M[f](s) = ∫_0^∞ f(t) t^(s-1) dt
//     exp(-t)         ->  Γ(s)
//     1/(1 + t)       ->  π/sin(πs)
//     t^a · g(t) ↔ shift on s.
//
// Z-transform (sequence variable n, transformed variable z):
//     1               ->  z/(z - 1)
//     n               ->  z/(z - 1)^2
//     a^n             ->  z/(z - a)
//     n·a^n           ->  a·z/(z - a)^2
//
// Sine/cosine transforms (half-line Fourier).
//     sine_transform     F_s(ω) = ∫_0^∞ f(t) sin(ωt) dt
//     cosine_transform   F_c(ω) = ∫_0^∞ f(t) cos(ωt) dt
//   Pairs (with a > 0):
//     exp(-a·t)  ─sin→  ω/(a² + ω²)
//     exp(-a·t)  ─cos→  a/(a² + ω²)
//
// Reference: sympy/integrals/{laplace,transforms}.py

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// ---- Laplace ---------------------------------------------------------------
[[nodiscard]] SYMPP_EXPORT Expr laplace_transform(const Expr& f,
                                                    const Expr& t,
                                                    const Expr& s);
[[nodiscard]] SYMPP_EXPORT Expr inverse_laplace_transform(const Expr& F,
                                                           const Expr& s,
                                                           const Expr& t);

// ---- Fourier ---------------------------------------------------------------
[[nodiscard]] SYMPP_EXPORT Expr fourier_transform(const Expr& f,
                                                    const Expr& t,
                                                    const Expr& w);
[[nodiscard]] SYMPP_EXPORT Expr inverse_fourier_transform(const Expr& F,
                                                            const Expr& w,
                                                            const Expr& t);

// ---- Mellin ----------------------------------------------------------------
[[nodiscard]] SYMPP_EXPORT Expr mellin_transform(const Expr& f,
                                                   const Expr& t,
                                                   const Expr& s);
[[nodiscard]] SYMPP_EXPORT Expr inverse_mellin_transform(const Expr& F,
                                                          const Expr& s,
                                                          const Expr& t);

// ---- Z ---------------------------------------------------------------------
[[nodiscard]] SYMPP_EXPORT Expr z_transform(const Expr& f,
                                              const Expr& n,
                                              const Expr& z);
[[nodiscard]] SYMPP_EXPORT Expr inverse_z_transform(const Expr& F,
                                                     const Expr& z,
                                                     const Expr& n);

// ---- Half-line sine / cosine ----------------------------------------------
[[nodiscard]] SYMPP_EXPORT Expr sine_transform(const Expr& f,
                                                 const Expr& t,
                                                 const Expr& w);
[[nodiscard]] SYMPP_EXPORT Expr cosine_transform(const Expr& f,
                                                   const Expr& t,
                                                   const Expr& w);

}  // namespace sympp
