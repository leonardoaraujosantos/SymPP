#pragma once

// Symbolic physics — the reusable, self-contained pieces of sympy.physics.
// (The deep submodules — full continuum mechanics, relativistic field theory —
// remain genuinely multi-week ports. Second quantization is covered in part:
// finite-dim ladder/number matrices, Jordan–Wigner fermions, single-mode
// bosonic / fermionic Fock states, and multi-mode bosonic Fock states (an
// occupation vector with per-mode ladder operators) with
// apply_annihilation/creation/number below.)
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

// ----- Second quantization: bosonic Fock states ------------------------------
// A single-mode bosonic Fock state c·|n⟩: an occupation number n ≥ 0 paired with
// a symbolic scalar coefficient c. The vacuum / null result a|0⟩ = 0 is encoded
// as a state with a zero coefficient (occupation is then irrelevant).
struct FockState {
    std::size_t occupation = 0;  // n in |n⟩
    Expr coefficient;            // scalar prefactor c (c·|n⟩)

    // Whether this represents the zero vector (coefficient simplifies to 0).
    [[nodiscard]] SYMPP_EXPORT bool is_zero() const;
};

// Annihilation a (c·|n⟩) = c·√n·|n−1⟩;  a|0⟩ = 0 (zero coefficient).
[[nodiscard]] SYMPP_EXPORT FockState apply_annihilation(const FockState& s);
// Creation a† (c·|n⟩) = c·√(n+1)·|n+1⟩.
[[nodiscard]] SYMPP_EXPORT FockState apply_creation(const FockState& s);
// Number N (c·|n⟩) = c·n·|n⟩.
[[nodiscard]] SYMPP_EXPORT FockState apply_number(const FockState& s);
// Commutator [a, a†] (c·|n⟩) = (a·a† − a†·a)(c·|n⟩) = c·|n⟩ (canonical relation).
[[nodiscard]] SYMPP_EXPORT FockState apply_boson_commutator(const FockState& s);

// ----- Second quantization: multi-mode bosonic Fock states -------------------
// A multi-mode bosonic Fock state c·|n₀ n₁ … n_{M−1}⟩: an occupation vector over
// M independent modes paired with a symbolic scalar coefficient c. Each mode
// carries its own independent ladder operators aᵢ, aᵢ† satisfying
//   [aᵢ, aⱼ†] = δᵢⱼ,  [aᵢ, aⱼ] = 0,  [aᵢ†, aⱼ†] = 0.
// As for the single-mode case the zero vector is encoded by a zero coefficient.
struct MultiModeFockState {
    std::vector<int> occupation;  // (n₀, …, n_{M−1}), each nᵢ ≥ 0
    Expr coefficient;             // scalar prefactor c (c·|n₀ … n_{M−1}⟩)

    // Whether this represents the zero vector (coefficient simplifies to 0).
    [[nodiscard]] SYMPP_EXPORT bool is_zero() const;
};

// Annihilation aᵢ (c·|…nᵢ…⟩) = c·√nᵢ·|…nᵢ−1…⟩;  aᵢ with nᵢ=0 gives 0.
[[nodiscard]] SYMPP_EXPORT MultiModeFockState apply_annihilation(const MultiModeFockState& s,
                                                                 std::size_t mode);
// Creation aᵢ† (c·|…nᵢ…⟩) = c·√(nᵢ+1)·|…nᵢ+1…⟩.
[[nodiscard]] SYMPP_EXPORT MultiModeFockState apply_creation(const MultiModeFockState& s,
                                                             std::size_t mode);
// Number Nᵢ (c·|…nᵢ…⟩) = c·nᵢ·|…nᵢ…⟩.
[[nodiscard]] SYMPP_EXPORT MultiModeFockState apply_number(const MultiModeFockState& s,
                                                           std::size_t mode);
// Commutator [aᵢ, aⱼ†] applied to a state: (aᵢ·aⱼ† − aⱼ†·aᵢ)(c·|n⟩) = δᵢⱼ·c·|n⟩.
[[nodiscard]] SYMPP_EXPORT MultiModeFockState apply_boson_commutator(const MultiModeFockState& s,
                                                                     std::size_t i, std::size_t j);
// Commutator [aᵢ, aⱼ] applied to a state: (aᵢ·aⱼ − aⱼ·aᵢ)(c·|n⟩) = 0.
[[nodiscard]] SYMPP_EXPORT MultiModeFockState apply_annihilation_commutator(
    const MultiModeFockState& s, std::size_t i, std::size_t j);

// ----- Second quantization: fermionic Fock states ----------------------------
// A single-mode fermionic Fock state c·|n⟩ with Pauli exclusion: the occupation
// number n is restricted to {0, 1}. As for the bosonic case, the zero vector is
// encoded by a zero coefficient (occupation then irrelevant).
struct FermionState {
    std::size_t occupation = 0;  // n in |n⟩, must be 0 or 1
    Expr coefficient;            // scalar prefactor c (c·|n⟩)

    // Whether this represents the zero vector (coefficient simplifies to 0).
    [[nodiscard]] SYMPP_EXPORT bool is_zero() const;
};

// Annihilation a (c·|1⟩) = c·|0⟩;  a|0⟩ = 0 (zero coefficient).
[[nodiscard]] SYMPP_EXPORT FermionState apply_fermion_annihilation(const FermionState& s);
// Creation a† (c·|0⟩) = c·|1⟩;  a†|1⟩ = 0 (Pauli exclusion).
[[nodiscard]] SYMPP_EXPORT FermionState apply_fermion_creation(const FermionState& s);
// Number N (c·|n⟩) = c·n·|n⟩ (eigenvalue n ∈ {0, 1}).
[[nodiscard]] SYMPP_EXPORT FermionState apply_fermion_number(const FermionState& s);
// Anticommutator {a, a†} (c·|n⟩) = (a·a† + a†·a)(c·|n⟩) = c·|n⟩ (CAR).
[[nodiscard]] SYMPP_EXPORT FermionState apply_fermion_anticommutator(const FermionState& s);

// ----- Second quantization: normal ordering of operator words ----------------
// An operator word is a product of single-particle ladder operators written
// left-to-right. Each factor is a (mode, dagger?) pair: dagger=true is a
// creation operator a_mode†, dagger=false an annihilation operator a_mode.
//   e.g. {{0,false},{0,true}} represents  a₀ a₀†.
struct LadderOp {
    std::size_t mode = 0;  // single-particle mode index
    bool dagger = false;   // true: creation a†, false: annihilation a

    [[nodiscard]] bool operator==(const LadderOp& o) const {
        return mode == o.mode && dagger == o.dagger;
    }
    [[nodiscard]] bool operator!=(const LadderOp& o) const { return !(*this == o); }
};
using OperatorWord = std::vector<LadderOp>;

// A single normal-ordered term: a scalar coefficient times an operator word
// whose creation operators all stand to the left of its annihilation operators.
struct OperatorTerm {
    Expr coefficient;   // scalar prefactor
    OperatorWord word;  // normal-ordered factors (a†…a† a…a)
};

// The particle statistics used by normal_order.
enum class Statistics { Boson, Fermion };

// Normal-order a product of ladder operators, returning the equivalent linear
// combination of normal-ordered words (creation operators moved to the left).
// Uses the canonical (anti)commutation relations and accumulates contractions:
//   bosons:   a_i a_j† = δ_ij + a_j† a_i             ([a_i,a_j†]=δ_ij)
//   fermions: a_i a_j† = δ_ij − a_j† a_i             ({a_i,a_j†}=δ_ij)
// For fermions a word with two equal adjacent operators vanishes after ordering
// (a_i† a_i† = 0, a_i a_i = 0). Zero-coefficient terms are dropped; the empty
// result {} denotes the operator 0, a term with an empty word the identity.
[[nodiscard]] SYMPP_EXPORT std::vector<OperatorTerm> normal_order(const OperatorWord& word,
                                                                  Statistics stats);

// ----- Second quantization: fermionic ladder operators -----------------------
// Jordan–Wigner representation of fermionic mode `p` (0-based) among `n` modes
// on the 2ⁿ-dimensional Fock space: c_p = Z⊗…⊗Z⊗σ⁻⊗I⊗…⊗I (p copies of Z).
// These satisfy the canonical anticommutation relations exactly (finite-dim):
//   {c_p, c_q†} = δ_pq·I,  {c_p, c_q} = 0,  c_p² = 0.
[[nodiscard]] SYMPP_EXPORT Matrix fermion_annihilation(std::size_t p, std::size_t n);
[[nodiscard]] SYMPP_EXPORT Matrix fermion_creation(std::size_t p, std::size_t n);
// Occupation-number operator nₚ = c_p†·c_p (eigenvalues 0/1).
[[nodiscard]] SYMPP_EXPORT Matrix fermion_number(std::size_t p, std::size_t n);

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

// Lagrangian mechanics: equations of motion for a system with generalized
// coordinates q[i], velocities qdot[i], accelerations qddot[i] and time t. The
// Lagrangian L is written in those symbols; for each coordinate the returned
// expression is the Euler–Lagrange LHS  ∂L/∂qᵢ − D_t(∂L/∂q̇ᵢ)  (set = 0 for the
// EOM), where D_t is the total time derivative
//   D_t f = ∂f/∂t + Σⱼ (∂f/∂qⱼ)·q̇ⱼ + Σⱼ (∂f/∂q̇ⱼ)·q̈ⱼ.
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> lagrange_equations(
    const Expr& lagrangian, const std::vector<Expr>& q, const std::vector<Expr>& qdot,
    const std::vector<Expr>& qddot, const Expr& t);

}  // namespace sympp::physics
