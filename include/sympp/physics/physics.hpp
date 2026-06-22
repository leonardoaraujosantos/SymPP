#pragma once

// Symbolic physics — the reusable, self-contained pieces of sympy.physics.
// (The deep submodules — second quantization, full continuum mechanics,
// relativistic field theory — remain genuinely multi-week ports.)
//
//   * quantum     — commutators, Pauli matrices, ladder/number operators,
//                   arbitrary-spin angular-momentum operators (Jx,Jy,Jz,J±,J²),
//                   Wigner 3-j / 6-j and Clebsch–Gordan coupling coefficients
//   * atomic      — hydrogen energy levels & radial wavefunctions, 1D quantum
//                   harmonic-oscillator energies & wavefunctions
//   * optics      — ABCD ray-transfer matrices + Gaussian-beam q propagation
//   * mechanics   — canonical (conjugate) momentum and the Hamiltonian
//                   (Legendre transform) from a Lagrangian
//
// For Lagrangian equations of motion use the existing euler_lagrange() in
// <calculus/euler_lagrange.hpp>.
//
// Reference: sympy/physics/{quantum,hydrogen,qho_1d,wigner,optics,mechanics}.

#include <cstddef>
#include <vector>

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

// ----- Angular momentum / spin (in units of ℏ) -------------------------------
// Spin-j matrix representations in the standard basis |j,m⟩ ordered
// m = j, j−1, …, −j (SymPy's convention). j may be integer or half-integer.
[[nodiscard]] SYMPP_EXPORT Matrix spin_jz(const Expr& j);
[[nodiscard]] SYMPP_EXPORT Matrix spin_jplus(const Expr& j);   // raising  J₊
[[nodiscard]] SYMPP_EXPORT Matrix spin_jminus(const Expr& j);  // lowering J₋
[[nodiscard]] SYMPP_EXPORT Matrix spin_jx(const Expr& j);      // (J₊+J₋)/2
[[nodiscard]] SYMPP_EXPORT Matrix spin_jy(const Expr& j);      // (J₊−J₋)/2i
[[nodiscard]] SYMPP_EXPORT Matrix spin_j2(const Expr& j);      // J² = j(j+1)·I

// ----- Wigner / angular-momentum coupling coefficients -----------------------
// Wigner 3-j symbol; returns S::Zero() when selection rules are violated.
[[nodiscard]] SYMPP_EXPORT Expr wigner_3j(const Expr& j1, const Expr& j2, const Expr& j3,
                                          const Expr& m1, const Expr& m2, const Expr& m3);
// Wigner 6-j symbol {j1 j2 j3 ; j4 j5 j6}.
[[nodiscard]] SYMPP_EXPORT Expr wigner_6j(const Expr& j1, const Expr& j2, const Expr& j3,
                                          const Expr& j4, const Expr& j5, const Expr& j6);
// Clebsch–Gordan coefficient ⟨j1 m1 j2 m2 | j3 m3⟩.
[[nodiscard]] SYMPP_EXPORT Expr clebsch_gordan(const Expr& j1, const Expr& j2, const Expr& j3,
                                               const Expr& m1, const Expr& m2, const Expr& m3);
// Wigner 9-j symbol {j1 j2 j3 ; j4 j5 j6 ; j7 j8 j9} (sum over a 6-j product).
[[nodiscard]] SYMPP_EXPORT Expr wigner_9j(const Expr& j1, const Expr& j2, const Expr& j3,
                                          const Expr& j4, const Expr& j5, const Expr& j6,
                                          const Expr& j7, const Expr& j8, const Expr& j9);
// Racah W coefficient W(a b c d ; e f) = (−1)^(a+b+c+d)·{a b e ; d c f}.
[[nodiscard]] SYMPP_EXPORT Expr racah(const Expr& a, const Expr& b, const Expr& c, const Expr& d,
                                      const Expr& e, const Expr& f);
// Gaunt coefficient ∫ Yₗ₁ᵐ¹ Yₗ₂ᵐ² Yₗ₃ᵐ³ dΩ.
[[nodiscard]] SYMPP_EXPORT Expr gaunt(const Expr& l1, const Expr& l2, const Expr& l3,
                                      const Expr& m1, const Expr& m2, const Expr& m3);
// Dirac gamma matrix γ^μ (μ = 0,1,2,3) or the chirality matrix γ⁵ (mu = 5),
// in the Dirac representation (SymPy's sympy.physics.matrices.mgamma).
[[nodiscard]] SYMPP_EXPORT Matrix dirac_gamma(int mu);

// ----- Hydrogen atom (Hartree atomic units, ℏ=mₑ=e=1) ------------------------
// Bohr energy level E_{n} = −Z²/(2n²) (l-degenerate, non-relativistic).
[[nodiscard]] SYMPP_EXPORT Expr hydrogen_energy(const Expr& n, const Expr& Z);
// Normalized radial wavefunction R_{n,l}(r) in Bohr units.
[[nodiscard]] SYMPP_EXPORT Expr hydrogen_R_nl(int n, int l, const Expr& r, const Expr& Z);

// ----- 1D quantum harmonic oscillator (ℏ=1) ----------------------------------
[[nodiscard]] SYMPP_EXPORT Expr qho_energy(const Expr& n, const Expr& omega);  // (n+½)ω
// Stationary-state wavefunction ψ_n(x) in natural units (m=ω=ℏ=1).
[[nodiscard]] SYMPP_EXPORT Expr qho_wavefunction(int n, const Expr& x);

// ----- Quantum computing (qubit states & gates) ------------------------------
// Computational-basis kets |0⟩, |1⟩ as 2×1 column vectors.
[[nodiscard]] SYMPP_EXPORT Matrix ket0();
[[nodiscard]] SYMPP_EXPORT Matrix ket1();
// n-qubit computational basis state |bits[0] bits[1] …⟩ (each 0/1), big-endian.
[[nodiscard]] SYMPP_EXPORT Matrix qubit_state(const std::vector<int>& bits);
// Single-qubit gates. (Pauli X/Y/Z above double as the bit/phase-flip gates.)
[[nodiscard]] SYMPP_EXPORT Matrix gate_hadamard();   // H
[[nodiscard]] SYMPP_EXPORT Matrix gate_phase();      // S = diag(1, i)
[[nodiscard]] SYMPP_EXPORT Matrix gate_t();          // T = diag(1, e^{iπ/4})
// Two-qubit controlled-NOT (control = first/most-significant qubit).
[[nodiscard]] SYMPP_EXPORT Matrix gate_cnot();
// Inner product ⟨a|b⟩ = a†·b for two kets of equal dimension (scalar).
[[nodiscard]] SYMPP_EXPORT Expr braket(const Matrix& a, const Matrix& b);

// ----- Geometric optics (ABCD ray-transfer matrices) -------------------------
[[nodiscard]] SYMPP_EXPORT Matrix abcd_free_space(const Expr& distance);
[[nodiscard]] SYMPP_EXPORT Matrix abcd_thin_lens(const Expr& focal_length);
[[nodiscard]] SYMPP_EXPORT Matrix abcd_curved_mirror(const Expr& radius);
// Flat dielectric interface (refraction), n1 → n2.
[[nodiscard]] SYMPP_EXPORT Matrix abcd_flat_interface(const Expr& n1, const Expr& n2);
// Propagate a Gaussian-beam complex parameter q through an ABCD system:
//   q′ = (A·q + B)/(C·q + D).
[[nodiscard]] SYMPP_EXPORT Expr gaussian_beam_q(const Matrix& abcd, const Expr& q);

// ----- Classical mechanics ---------------------------------------------------
// Canonical momentum conjugate to a velocity: p = ∂L/∂q̇.
[[nodiscard]] SYMPP_EXPORT Expr conjugate_momentum(const Expr& lagrangian,
                                                   const Expr& velocity);
// Hamiltonian via Legendre transform of a one-coordinate Lagrangian:
//   H = p·q̇ − L,  with p = ∂L/∂q̇  (returned in terms of q̇).
[[nodiscard]] SYMPP_EXPORT Expr hamiltonian(const Expr& lagrangian, const Expr& velocity);

}  // namespace sympp::physics
