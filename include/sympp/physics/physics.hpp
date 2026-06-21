#pragma once

// A representative core of symbolic physics — NOT a full port of sympy.physics
// (mechanics/quantum/optics are each large modules). Covers the most reusable,
// self-contained pieces:
//
//   * quantum     — operator commutators, Pauli matrices, ladder/number ops
//   * optics      — ABCD ray-transfer matrices (free space, thin lens, mirror)
//   * mechanics   — canonical (conjugate) momentum from a Lagrangian
//
// For Lagrangian equations of motion use the existing euler_lagrange() in
// <calculus/euler_lagrange.hpp>.
//
// Reference: sympy/physics/{quantum,optics,mechanics}.

#include <cstddef>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>
#include <sympp/matrices/matrix.hpp>

namespace sympp::physics {

// ----- Quantum ---------------------------------------------------------------
[[nodiscard]] SYMPP_EXPORT Matrix commutator(const Matrix& a, const Matrix& b);      // AB − BA
[[nodiscard]] SYMPP_EXPORT Matrix anticommutator(const Matrix& a, const Matrix& b);  // AB + BA
[[nodiscard]] SYMPP_EXPORT Matrix pauli_x();
[[nodiscard]] SYMPP_EXPORT Matrix pauli_y();
[[nodiscard]] SYMPP_EXPORT Matrix pauli_z();
// Harmonic-oscillator ladder operators truncated to an n-level Fock basis.
[[nodiscard]] SYMPP_EXPORT Matrix annihilation_operator(std::size_t n);
[[nodiscard]] SYMPP_EXPORT Matrix creation_operator(std::size_t n);
[[nodiscard]] SYMPP_EXPORT Matrix number_operator(std::size_t n);

// ----- Geometric optics (ABCD ray-transfer matrices) -------------------------
[[nodiscard]] SYMPP_EXPORT Matrix abcd_free_space(const Expr& distance);
[[nodiscard]] SYMPP_EXPORT Matrix abcd_thin_lens(const Expr& focal_length);
[[nodiscard]] SYMPP_EXPORT Matrix abcd_curved_mirror(const Expr& radius);

// ----- Classical mechanics ---------------------------------------------------
// Canonical momentum conjugate to a velocity: p = ∂L/∂q̇.
[[nodiscard]] SYMPP_EXPORT Expr conjugate_momentum(const Expr& lagrangian,
                                                   const Expr& velocity);

}  // namespace sympp::physics
