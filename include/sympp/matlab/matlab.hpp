#pragma once

// MATLAB Symbolic Math Toolbox-named facade.
//
// Most of these are thin re-exports of `sympp::*` functions under
// MATLAB-friendly names so users porting MATLAB scripts can drop in
// the symbolic ops with minimal renaming. A few cases (variadic
// `syms`, `Int` because `int` is a C++ keyword, `taylor` mapping to
// `series`) need a real wrapper.
//
// Usage:
//
//     using namespace sympp::matlab;
//     auto [x, y] = syms("x", "y");
//     auto f = sin(x) * exp(x);
//     auto fp = diff(f, x);
//     auto F  = Int(f, x);
//     auto roots = solve(pow(x, integer(2)) - integer(4), x);
//     std::cout << ccode(fp) << "\n";
//
// MATLAB → SymPP renamings (only what differs from the SymPP API):
//
//     int        →  Int          (C++ keyword)
//     ilaplace   →  ilaplace     (matches; calls inverse_laplace_transform)
//     ifourier   →  ifourier     (matches; calls inverse_fourier_transform)
//     ztrans     →  ztrans       (calls z_transform)
//     iztrans    →  iztrans      (calls inverse_z_transform)
//     vpa        →  vpa          (calls evalf)
//     taylor     →  taylor       (calls series)
//     fortran    →  fortran      (calls printing::fcode)
//     matlabFunction → matlabFunction (calls printing::octave_function)
//     pretty     →  pretty       (calls printing::pretty)
//     ccode      →  ccode        (calls printing::ccode)
//     latex      →  latex        (calls printing::latex)
//
// Reference: MATLAB Symbolic Math Toolbox User's Guide, function list.

#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/calculus/limit.hpp>
#include <sympp/calculus/series.hpp>
#include <sympp/core/api.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/fwd.hpp>
#include <sympp/integrals/integrate.hpp>
#include <sympp/integrals/transforms.hpp>
#include <sympp/matlab/assumptions.hpp>
#include <sympp/matlab/ode.hpp>
#include <sympp/matlab/parsing.hpp>
#include <sympp/matlab/pde.hpp>
#include <sympp/matlab/solvers.hpp>
#include <sympp/ode/dsolve.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/printing/printing.hpp>
#include <sympp/simplify/simplify.hpp>
#include <sympp/solvers/solve.hpp>

namespace sympp::matlab {

// ---------------------------------------------------------------------------
// Symbol declaration
// ---------------------------------------------------------------------------

// syms("x") — single symbol
[[nodiscard]] inline Expr syms(std::string_view name) {
    return symbol(name);
}

// syms("x", "y", ...) — multiple symbols, structured-bindable.
//
//     auto [x, y, z] = syms("x", "y", "z");
template <typename... Names>
[[nodiscard]] inline auto syms(std::string_view first,
                                  std::string_view second,
                                  Names... rest) {
    return std::make_tuple(symbol(first), symbol(second), symbol(rest)...);
}

// sym(name) / sym(value) overloads live in matlab/parsing.hpp — a
// bare identifier becomes a Symbol, an expression string is parsed,
// and an integer literal becomes an Integer.

// ---------------------------------------------------------------------------
// Variable-precision arithmetic
// ---------------------------------------------------------------------------

[[nodiscard]] inline Expr vpa(const Expr& e, int dps = 32) {
    return evalf(e, dps);
}

// ---------------------------------------------------------------------------
// Calculus
// ---------------------------------------------------------------------------

// diff(f, x) — first derivative
[[nodiscard]] inline Expr diff(const Expr& f, const Expr& var) {
    return sympp::diff(f, var);
}

// diff(f, x, n) — n-th derivative
[[nodiscard]] inline Expr diff(const Expr& f, const Expr& var, std::size_t n) {
    Expr d = f;
    for (std::size_t i = 0; i < n; ++i) d = sympp::diff(d, var);
    return d;
}

// Int — integration. (`int` is a C++ keyword.)
[[nodiscard]] inline Expr Int(const Expr& f, const Expr& var) {
    return integrate(f, var);
}
[[nodiscard]] inline Expr Int(const Expr& f, const Expr& var,
                                const Expr& a, const Expr& b) {
    return integrate(f, var, a, b);
}

// limit(f, x, a)
[[nodiscard]] inline Expr limit(const Expr& f, const Expr& var,
                                  const Expr& target) {
    return sympp::limit(f, var, target);
}

// taylor(f, x, a, n) — Taylor series about x = a, truncated at order n.
[[nodiscard]] inline Expr taylor(const Expr& f, const Expr& var,
                                   const Expr& a, std::size_t n) {
    return series(f, var, a, n);
}

// ---------------------------------------------------------------------------
// Algebraic manipulation
// ---------------------------------------------------------------------------

[[nodiscard]] inline Expr simplify(const Expr& e) {
    return sympp::simplify(e);
}

[[nodiscard]] inline Expr subs(const Expr& e, const Expr& old_e,
                                  const Expr& new_e) {
    return sympp::subs(e, old_e, new_e);
}

[[nodiscard]] inline Expr expand(const Expr& e) {
    return sympp::expand(e);
}

[[nodiscard]] inline Expr factor(const Expr& e, const Expr& var) {
    return sympp::factor(e, var);
}

// ---------------------------------------------------------------------------
// Solvers
// ---------------------------------------------------------------------------

// solve(eq, var) — eq is treated as the LHS with implicit RHS = 0,
// matching MATLAB's solve(eq) when the user writes solve(x^2 - 4) etc.
// (MATLAB also accepts solve(lhs == rhs); for that, write
// solve(lhs - rhs).)
[[nodiscard]] inline std::vector<Expr> solve(const Expr& eq,
                                                const Expr& var) {
    return sympp::solve(eq, var);
}

// dsolve(...) overloads live in matlab/ode.hpp — first-order,
// second-order with auto-classification, IVP variants, and the
// linear-system overload dsolve(A, x).

// ---------------------------------------------------------------------------
// Transforms
// ---------------------------------------------------------------------------

[[nodiscard]] inline Expr laplace(const Expr& f, const Expr& t,
                                    const Expr& s) {
    return laplace_transform(f, t, s);
}
[[nodiscard]] inline Expr ilaplace(const Expr& F, const Expr& s,
                                     const Expr& t) {
    return inverse_laplace_transform(F, s, t);
}
[[nodiscard]] inline Expr fourier(const Expr& f, const Expr& t,
                                    const Expr& w) {
    return fourier_transform(f, t, w);
}
[[nodiscard]] inline Expr ifourier(const Expr& F, const Expr& w,
                                     const Expr& t) {
    return inverse_fourier_transform(F, w, t);
}
[[nodiscard]] inline Expr ztrans(const Expr& f, const Expr& n,
                                   const Expr& z) {
    return z_transform(f, n, z);
}
[[nodiscard]] inline Expr iztrans(const Expr& F, const Expr& z,
                                    const Expr& n) {
    return inverse_z_transform(F, z, n);
}

// ---------------------------------------------------------------------------
// Code generation / printing
// ---------------------------------------------------------------------------

[[nodiscard]] inline std::string pretty(const Expr& e) {
    return printing::pretty(e);
}
[[nodiscard]] inline std::string latex(const Expr& e) {
    return printing::latex(e);
}
[[nodiscard]] inline std::string ccode(const Expr& e) {
    return printing::ccode(e);
}
[[nodiscard]] inline std::string fortran(const Expr& e) {
    return printing::fcode(e);
}

// matlabFunction(f, vars) — emit an Octave/MATLAB function with the
// given name. MATLAB's matlabFunction also returns a function handle;
// we return the source string. Pair with file I/O on the consumer
// side if you want a .m file.
[[nodiscard]] inline std::string matlabFunction(
    const std::string& name, const Expr& body,
    const std::vector<Expr>& args) {
    return printing::octave_function(name, body, args);
}

}  // namespace sympp::matlab
