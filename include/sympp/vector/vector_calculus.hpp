#pragma once

// Vector calculus and (Riemannian) differential geometry over an explicit list
// of coordinate symbols.
//
//   auto x = symbol("x"), y = symbol("y"), z = symbol("z");
//   gradient(x*y*z, {x, y, z});                 // → {y*z, x*z, x*y}
//   divergence({x, y, z}, {x, y, z});           // → 3
//   curl({-y, x, 0}, {x, y, z});                // → {0, 0, 2}
//
// Differential geometry takes a metric (a symmetric Matrix g_{ij} in the given
// coordinates) and returns Christoffel symbols and curvature.
//
// Reference: sympy/vector/*, sympy/diffgeom/*.

#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>
#include <sympp/matrices/matrix.hpp>

namespace sympp::vector {

using VectorField = std::vector<Expr>;

// ----- Vector calculus -------------------------------------------------------
[[nodiscard]] SYMPP_EXPORT VectorField gradient(const Expr& f,
                                                const std::vector<Expr>& coords);
[[nodiscard]] SYMPP_EXPORT Expr divergence(const VectorField& field,
                                           const std::vector<Expr>& coords);
// 3D curl (requires exactly three coordinates and three components).
[[nodiscard]] SYMPP_EXPORT VectorField curl(const VectorField& field,
                                            const std::vector<Expr>& coords);
[[nodiscard]] SYMPP_EXPORT Expr laplacian(const Expr& f, const std::vector<Expr>& coords);

// ----- Differential geometry -------------------------------------------------
// Christoffel symbols of the second kind: result[k](i, j) = Γᵏᵢⱼ.
using ChristoffelSymbols = std::vector<Matrix>;
[[nodiscard]] SYMPP_EXPORT ChristoffelSymbols christoffel(const Matrix& metric,
                                                          const std::vector<Expr>& coords);
// Ricci curvature tensor Rᵢⱼ and the Ricci scalar R = gⁱʲ Rᵢⱼ.
[[nodiscard]] SYMPP_EXPORT Matrix ricci_tensor(const Matrix& metric,
                                               const std::vector<Expr>& coords);
[[nodiscard]] SYMPP_EXPORT Expr ricci_scalar(const Matrix& metric,
                                             const std::vector<Expr>& coords);

}  // namespace sympp::vector
