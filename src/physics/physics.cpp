#include <sympp/physics/physics.hpp>

#include <stdexcept>
#include <vector>

#include <gmpxx.h>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/orthopolys.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp::physics {

namespace {

// Read an Integer/Rational Expr as an exact rational. Throws on anything else.
[[nodiscard]] mpq_class to_mpq(const Expr& e) {
    if (e->type_id() == TypeId::Integer) {
        return mpq_class(static_cast<const Integer&>(*e).value());
    }
    if (e->type_id() == TypeId::Rational) {
        return static_cast<const Rational&>(*e).value();
    }
    throw std::invalid_argument("physics: expected a numeric (integer/rational) argument");
}

[[nodiscard]] Expr from_mpq(const mpq_class& q) {
    return rational(q.get_num(), q.get_den());
}

[[nodiscard]] bool is_nonneg_integer(const mpq_class& q) {
    return q.get_den() == 1 && sgn(q) >= 0;
}

// Factorial of a non-negative integer-valued rational.
[[nodiscard]] mpz_class fac(const mpq_class& q) {
    mpz_class r;
    mpz_fac_ui(r.get_mpz_t(), q.get_num().get_ui());
    return r;
}

// (-1)^k for integer-valued k.
[[nodiscard]] int phase(const mpq_class& k) {
    return (k.get_num() % 2 == 0) ? 1 : -1;
}

// Triangle coefficient Δ(a,b,c) = (a+b−c)!(a−b+c)!(−a+b+c)! / (a+b+c+1)!  (rational).
[[nodiscard]] mpq_class triangle_delta(const mpq_class& a, const mpq_class& b, const mpq_class& c) {
    mpq_class num = mpq_class(fac(a + b - c)) * fac(a - b + c) * fac(-a + b + c);
    mpq_class res = num / mpq_class(fac(a + b + c + 1));
    res.canonicalize();
    return res;
}

// A triad (a,b,c) satisfies the angular-momentum triangle rule.
[[nodiscard]] bool triangle_ok(const mpq_class& a, const mpq_class& b, const mpq_class& c) {
    return is_nonneg_integer(a + b - c) && is_nonneg_integer(a - b + c) &&
           is_nonneg_integer(-a + b + c);
}

}  // namespace

Matrix commutator(const Matrix& a, const Matrix& b) { return a * b - b * a; }
Matrix anticommutator(const Matrix& a, const Matrix& b) { return a * b + b * a; }

Matrix pauli_x() {
    return Matrix{{integer(0), integer(1)}, {integer(1), integer(0)}};
}
Matrix pauli_y() {
    Expr i = S::I();
    return Matrix{{integer(0), mul(integer(-1), i)}, {i, integer(0)}};
}
Matrix pauli_z() {
    return Matrix{{integer(1), integer(0)}, {integer(0), integer(-1)}};
}

Matrix annihilation_operator(std::size_t n) {
    Matrix a = Matrix::zeros(n, n);
    for (std::size_t k = 0; k + 1 < n; ++k) {
        a.set(k, k + 1, sqrt(integer(static_cast<long>(k + 1))));  // a|k+1> = √(k+1)|k>
    }
    return a;
}
Matrix creation_operator(std::size_t n) { return annihilation_operator(n).transpose(); }
Matrix number_operator(std::size_t n) {
    Matrix N = Matrix::zeros(n, n);
    for (std::size_t k = 0; k < n; ++k) N.set(k, k, integer(static_cast<long>(k)));
    return N;
}

// ----- Second quantization: bosonic Fock states ------------------------------

bool FockState::is_zero() const {
    if (!coefficient) return true;  // default-constructed (null) coefficient
    return simplify(coefficient) == S::Zero();
}

FockState apply_annihilation(const FockState& s) {
    Expr c = s.coefficient ? s.coefficient : S::One();
    if (s.occupation == 0) {
        return FockState{0, S::Zero()};  // a|0⟩ = 0
    }
    Expr factor = sqrt(integer(static_cast<long>(s.occupation)));  // √n
    return FockState{s.occupation - 1, simplify(c * factor)};
}

FockState apply_creation(const FockState& s) {
    Expr c = s.coefficient ? s.coefficient : S::One();
    Expr factor = sqrt(integer(static_cast<long>(s.occupation + 1)));  // √(n+1)
    return FockState{s.occupation + 1, simplify(c * factor)};
}

FockState apply_number(const FockState& s) {
    Expr c = s.coefficient ? s.coefficient : S::One();
    Expr factor = integer(static_cast<long>(s.occupation));  // n
    return FockState{s.occupation, simplify(c * factor)};
}

FockState apply_boson_commutator(const FockState& s) {
    // [a, a†] = a·a† − a†·a, applied to a state.
    FockState aad = apply_annihilation(apply_creation(s));  // a·a†|n⟩ = (n+1)|n⟩
    FockState ada = apply_creation(apply_annihilation(s));  // a†·a|n⟩ = n|n⟩
    // Both terms (when nonzero) share occupation n; combine coefficients.
    Expr coeff = simplify(aad.coefficient - ada.coefficient);
    return FockState{s.occupation, coeff};
}

// ----- Quantum computing -----------------------------------------------------

Matrix ket0() { return Matrix{{integer(1)}, {integer(0)}}; }
Matrix ket1() { return Matrix{{integer(0)}, {integer(1)}}; }

Matrix qubit_state(const std::vector<int>& bits) {
    if (bits.empty()) throw std::invalid_argument("qubit_state: need at least one qubit");
    Matrix state = (bits[0] == 0) ? ket0() : ket1();
    for (std::size_t i = 1; i < bits.size(); ++i) {
        state = kronecker_product(state, bits[i] == 0 ? ket0() : ket1());
    }
    return state;
}

Matrix gate_hadamard() {
    Expr s = rational(1, 2);  // 1/√2 = (1/2)^{1/2}
    Expr inv_sqrt2 = pow(integer(2), rational(-1, 2));
    (void)s;
    return Matrix{{inv_sqrt2, inv_sqrt2}, {inv_sqrt2, mul(integer(-1), inv_sqrt2)}};
}

Matrix gate_phase() {  // S = diag(1, i)
    return Matrix{{integer(1), integer(0)}, {integer(0), S::I()}};
}

Matrix gate_t() {  // T = diag(1, e^{iπ/4})
    Expr ph = exp(mul(S::I(), mul(rational(1, 4), S::Pi())));
    return Matrix{{integer(1), integer(0)}, {integer(0), ph}};
}

Matrix gate_cnot() {
    Expr o = integer(0), l = integer(1);
    // Basis |00>,|01>,|10>,|11>: flip the second qubit when the first is 1.
    return Matrix{{l, o, o, o}, {o, l, o, o}, {o, o, o, l}, {o, o, l, o}};
}

Expr braket(const Matrix& a, const Matrix& b) {
    Matrix inner = a.conjugate_transpose() * b;  // (1×n)·(n×1) → 1×1
    if (inner.rows() != 1 || inner.cols() != 1) {
        throw std::invalid_argument("braket: arguments must be kets of equal dimension");
    }
    return inner.at(0, 0);
}

Matrix abcd_free_space(const Expr& d) {
    return Matrix{{integer(1), d}, {integer(0), integer(1)}};
}
Matrix abcd_thin_lens(const Expr& f) {
    return Matrix{{integer(1), integer(0)},
                  {mul(integer(-1), pow(f, integer(-1))), integer(1)}};
}
Matrix abcd_curved_mirror(const Expr& r) {
    // Concave mirror, radius R: lower-left = −2/R.
    return Matrix{{integer(1), integer(0)},
                  {mul(integer(-2), pow(r, integer(-1))), integer(1)}};
}
Matrix abcd_flat_interface(const Expr& n1, const Expr& n2) {
    // Refraction at a planar dielectric boundary: D = n1/n2.
    return Matrix{{integer(1), integer(0)},
                  {integer(0), mul(n1, pow(n2, integer(-1)))}};
}

Expr gaussian_beam_q(const Matrix& abcd, const Expr& q) {
    // q′ = (A·q + B)/(C·q + D).
    Expr a = abcd.at(0, 0), b = abcd.at(0, 1), c = abcd.at(1, 0), d = abcd.at(1, 1);
    return mul(add(mul(a, q), b), pow(add(mul(c, q), d), integer(-1)));
}

Expr conjugate_momentum(const Expr& lagrangian, const Expr& velocity) {
    return diff(lagrangian, velocity);
}

Expr hamiltonian(const Expr& lagrangian, const Expr& velocity) {
    Expr p = conjugate_momentum(lagrangian, velocity);  // p = ∂L/∂q̇
    return mul(p, velocity) - lagrangian;               // H = p·q̇ − L
}

std::vector<Expr> lagrange_equations(const Expr& L, const std::vector<Expr>& q,
                                     const std::vector<Expr>& qdot,
                                     const std::vector<Expr>& qddot, const Expr& t) {
    const std::size_t dof = q.size();
    if (qdot.size() != dof || qddot.size() != dof) {
        throw std::invalid_argument("lagrange_equations: coordinate/velocity/accel "
                                    "lists must have equal length");
    }
    std::vector<Expr> eom;
    eom.reserve(dof);
    for (std::size_t i = 0; i < dof; ++i) {
        Expr p = diff(L, qdot[i]);  // generalized momentum pᵢ = ∂L/∂q̇ᵢ
        // Total time derivative D_t pᵢ via the chain rule over all coordinates.
        Expr dt_p = diff(p, t);
        for (std::size_t j = 0; j < dof; ++j) {
            dt_p = add(dt_p, mul(diff(p, q[j]), qdot[j]));
            dt_p = add(dt_p, mul(diff(p, qdot[j]), qddot[j]));
        }
        eom.push_back(simplify(diff(L, q[i]) - dt_p));  // ∂L/∂qᵢ − D_t pᵢ
    }
    return eom;
}

// ----- Second quantization: fermionic ladder operators -----------------------

Matrix fermion_annihilation(std::size_t p, std::size_t n) {
    if (p >= n) throw std::invalid_argument("fermion_annihilation: mode index out of range");
    Matrix Z = pauli_z();
    Matrix sigma_minus{{integer(0), integer(1)}, {integer(0), integer(0)}};  // σ⁻
    Matrix I2 = Matrix::identity(2);
    // Build Z^{⊗p} ⊗ σ⁻ ⊗ I^{⊗(n−p−1)} left-to-right.
    Matrix result{{integer(1)}};  // 1×1 identity seed for the Kronecker fold
    for (std::size_t i = 0; i < n; ++i) {
        const Matrix& factor = (i < p) ? Z : (i == p) ? sigma_minus : I2;
        result = kronecker_product(result, factor);
    }
    return result;
}

Matrix fermion_creation(std::size_t p, std::size_t n) {
    return fermion_annihilation(p, n).conjugate_transpose();
}

Matrix fermion_number(std::size_t p, std::size_t n) {
    return fermion_creation(p, n) * fermion_annihilation(p, n);
}

// ----- Angular momentum / spin -----------------------------------------------

namespace {
// Common setup: validate j, return dimension and the list of m-values (desc).
[[nodiscard]] std::vector<mpq_class> spin_m_values(const Expr& j, std::size_t& dim) {
    mpq_class jq = to_mpq(j);
    mpq_class two_j_plus_1 = 2 * jq + 1;
    if (!is_nonneg_integer(two_j_plus_1) || sgn(jq) < 0) {
        throw std::invalid_argument("spin: j must be a non-negative (half-)integer");
    }
    dim = two_j_plus_1.get_num().get_ui();
    std::vector<mpq_class> ms;
    ms.reserve(dim);
    for (std::size_t k = 0; k < dim; ++k) ms.push_back(jq - mpq_class(static_cast<long>(k)));
    return ms;
}
}  // namespace

Matrix spin_jz(const Expr& j) {
    std::size_t dim = 0;
    auto ms = spin_m_values(j, dim);
    Matrix M = Matrix::zeros(dim, dim);
    for (std::size_t k = 0; k < dim; ++k) M.set(k, k, from_mpq(ms[k]));
    return M;
}

Matrix spin_jplus(const Expr& j) {
    std::size_t dim = 0;
    auto ms = spin_m_values(j, dim);
    mpq_class jq = to_mpq(j);
    Matrix M = Matrix::zeros(dim, dim);
    // ⟨m+1|J₊|m⟩ = √(j(j+1) − m(m+1)); in descending basis m=ms[k] sits at row k−1.
    for (std::size_t k = 1; k < dim; ++k) {
        mpq_class m = ms[k];
        mpq_class val = jq * (jq + 1) - m * (m + 1);
        M.set(k - 1, k, sqrt(from_mpq(val)));
    }
    return M;
}

Matrix spin_jminus(const Expr& j) { return spin_jplus(j).transpose(); }

Matrix spin_jx(const Expr& j) {
    return (spin_jplus(j) + spin_jminus(j)).scalar_mul(rational(1, 2));
}

Matrix spin_jy(const Expr& j) {
    // Jy = (J₊ − J₋)/(2i) = −i/2·(J₊ − J₋).
    Expr coeff = mul(rational(-1, 2), S::I());
    return (spin_jplus(j) - spin_jminus(j)).scalar_mul(coeff);
}

Matrix spin_j2(const Expr& j) {
    std::size_t dim = 0;
    (void)spin_m_values(j, dim);
    mpq_class jq = to_mpq(j);
    return Matrix::identity(dim).scalar_mul(from_mpq(jq * (jq + 1)));
}

// ----- Wigner / Clebsch–Gordan -----------------------------------------------

Expr wigner_3j(const Expr& j1e, const Expr& j2e, const Expr& j3e, const Expr& m1e,
               const Expr& m2e, const Expr& m3e) {
    mpq_class j1 = to_mpq(j1e), j2 = to_mpq(j2e), j3 = to_mpq(j3e);
    mpq_class m1 = to_mpq(m1e), m2 = to_mpq(m2e), m3 = to_mpq(m3e);

    // Selection rules.
    if (m1 + m2 + m3 != 0) return S::Zero();
    if (!triangle_ok(j1, j2, j3)) return S::Zero();
    if (!is_nonneg_integer(j1 + m1) || !is_nonneg_integer(j1 - m1) ||
        !is_nonneg_integer(j2 + m2) || !is_nonneg_integer(j2 - m2) ||
        !is_nonneg_integer(j3 + m3) || !is_nonneg_integer(j3 - m3)) {
        return S::Zero();
    }

    mpq_class delta = triangle_delta(j1, j2, j3);
    mpq_class prod = mpq_class(fac(j1 + m1)) * fac(j1 - m1) * fac(j2 + m2) * fac(j2 - m2) *
                     fac(j3 + m3) * fac(j3 - m3);
    mpq_class radicand = delta * prod;
    radicand.canonicalize();

    // Summation over t with all six factorial arguments non-negative integers.
    mpq_class sum = 0;
    mpq_class jsum = j1 + j2 + j3;
    long tmax = jsum.get_num().get_si();
    for (long t = 0; t <= tmax; ++t) {
        mpq_class tt(t);
        mpq_class a1 = tt, a2 = j3 - j2 + m1 + tt, a3 = j3 - j1 - m2 + tt;
        mpq_class a4 = j1 + j2 - j3 - tt, a5 = j1 - m1 - tt, a6 = j2 + m2 - tt;
        if (!is_nonneg_integer(a1) || !is_nonneg_integer(a2) || !is_nonneg_integer(a3) ||
            !is_nonneg_integer(a4) || !is_nonneg_integer(a5) || !is_nonneg_integer(a6)) {
            continue;
        }
        mpq_class denom = mpq_class(fac(a1)) * fac(a2) * fac(a3) * fac(a4) * fac(a5) * fac(a6);
        sum += mpq_class(phase(tt)) / denom;
    }
    sum.canonicalize();

    int ph = phase(j1 - j2 - m3);
    Expr result = mul(mul(integer(ph), from_mpq(sum)), sqrt(from_mpq(radicand)));
    return result;
}

Expr clebsch_gordan(const Expr& j1e, const Expr& j2e, const Expr& j3e, const Expr& m1e,
                    const Expr& m2e, const Expr& m3e) {
    mpq_class j1 = to_mpq(j1e), j2 = to_mpq(j2e), j3 = to_mpq(j3e), m3 = to_mpq(m3e);
    // CG = (−1)^(j1−j2+m3)·√(2j3+1)·( j1 j2  j3 ; m1 m2 −m3 ).
    Expr w3 = wigner_3j(j1e, j2e, j3e, m1e, m2e, from_mpq(-m3));
    int ph = phase(j1 - j2 + m3);
    Expr factor = mul(integer(ph), sqrt(from_mpq(2 * j3 + 1)));
    return mul(factor, w3);
}

Expr wigner_6j(const Expr& a, const Expr& b, const Expr& c, const Expr& d, const Expr& e,
               const Expr& f) {
    mpq_class j1 = to_mpq(a), j2 = to_mpq(b), j3 = to_mpq(c);
    mpq_class j4 = to_mpq(d), j5 = to_mpq(e), j6 = to_mpq(f);

    // All four triads must satisfy the triangle rule.
    if (!triangle_ok(j1, j2, j3) || !triangle_ok(j1, j5, j6) || !triangle_ok(j4, j2, j6) ||
        !triangle_ok(j4, j5, j3)) {
        return S::Zero();
    }

    mpq_class radicand = triangle_delta(j1, j2, j3) * triangle_delta(j1, j5, j6) *
                         triangle_delta(j4, j2, j6) * triangle_delta(j4, j5, j3);
    radicand.canonicalize();

    mpq_class s1 = j1 + j2 + j3, s2 = j1 + j5 + j6, s3 = j4 + j2 + j6, s4 = j4 + j5 + j3;
    mpq_class p1 = j1 + j2 + j4 + j5, p2 = j2 + j3 + j5 + j6, p3 = j1 + j3 + j4 + j6;

    mpq_class lo = std::max(std::max(s1, s2), std::max(s3, s4));
    mpq_class hi = std::min(std::min(p1, p2), p3);
    long tmin = lo.get_num().get_si();
    long tmax = hi.get_num().get_si();

    mpq_class sum = 0;
    for (long t = tmin; t <= tmax; ++t) {
        mpq_class tt(t);
        mpq_class b1 = tt - s1, b2 = tt - s2, b3 = tt - s3, b4 = tt - s4;
        mpq_class b5 = p1 - tt, b6 = p2 - tt, b7 = p3 - tt;
        if (!is_nonneg_integer(b1) || !is_nonneg_integer(b2) || !is_nonneg_integer(b3) ||
            !is_nonneg_integer(b4) || !is_nonneg_integer(b5) || !is_nonneg_integer(b6) ||
            !is_nonneg_integer(b7)) {
            continue;
        }
        mpq_class denom = mpq_class(fac(b1)) * fac(b2) * fac(b3) * fac(b4) * fac(b5) * fac(b6) *
                          fac(b7);
        sum += mpq_class(phase(tt)) * mpq_class(fac(tt + 1)) / denom;
    }
    sum.canonicalize();

    return mul(from_mpq(sum), sqrt(from_mpq(radicand)));
}

Expr wigner_9j(const Expr& j1e, const Expr& j2e, const Expr& j3e, const Expr& j4e,
               const Expr& j5e, const Expr& j6e, const Expr& j7e, const Expr& j8e,
               const Expr& j9e) {
    mpq_class j1 = to_mpq(j1e), j2 = to_mpq(j2e), j3 = to_mpq(j3e);
    mpq_class j4 = to_mpq(j4e), j5 = to_mpq(j5e), j6 = to_mpq(j6e);
    mpq_class j7 = to_mpq(j7e), j8 = to_mpq(j8e), j9 = to_mpq(j9e);

    // {j1 j2 j3; j4 j5 j6; j7 j8 j9} =
    //   Σ_k (−1)^{2k}(2k+1)·{j1 j4 j7; j8 j9 k}·{j2 j5 j8; j4 k j6}·{j3 j6 j9; k j1 j2}.
    // k runs in unit steps over the overlap of the three coupling ranges; the
    // 6-j factors vanish outside their triangles, so a safe superset suffices.
    mpq_class kmin = std::max(std::max(abs(j1 - j9), abs(j2 - j6)), abs(j4 - j8));
    mpq_class kmax = std::min(std::min(j1 + j9, j2 + j6), j4 + j8);

    Expr sum = integer(0);
    for (mpq_class k = kmin; k <= kmax; k += 1) {
        Expr ek = from_mpq(k);
        Expr w1 = wigner_6j(j1e, j4e, j7e, j8e, j9e, ek);
        Expr w2 = wigner_6j(j2e, j5e, j8e, j4e, ek, j6e);
        Expr w3 = wigner_6j(j3e, j6e, j9e, ek, j1e, j2e);
        int sign = phase(2 * k);  // (−1)^{2k}
        Expr weight = mul(integer(sign), from_mpq(2 * k + 1));
        sum = add(sum, mul(weight, mul(w1, mul(w2, w3))));
    }
    return sum;
}

Expr racah(const Expr& a, const Expr& b, const Expr& c, const Expr& d, const Expr& e,
           const Expr& f) {
    int sign = phase(to_mpq(a) + to_mpq(b) + to_mpq(c) + to_mpq(d));
    return mul(integer(sign), wigner_6j(a, b, e, d, c, f));
}

Expr gaunt(const Expr& l1, const Expr& l2, const Expr& l3, const Expr& m1, const Expr& m2,
           const Expr& m3) {
    mpq_class a = to_mpq(l1), b = to_mpq(l2), cc = to_mpq(l3);
    // √[(2l1+1)(2l2+1)(2l3+1)/(4π)]·( l1 l2 l3 ; 0 0 0 )·( l1 l2 l3 ; m1 m2 m3 ).
    mpq_class deg = (2 * a + 1) * (2 * b + 1) * (2 * cc + 1);
    Expr prefactor = sqrt(mul(from_mpq(deg), pow(mul(integer(4), S::Pi()), integer(-1))));
    Expr w0 = wigner_3j(l1, l2, l3, integer(0), integer(0), integer(0));
    Expr wm = wigner_3j(l1, l2, l3, m1, m2, m3);
    return mul(prefactor, mul(w0, wm));
}

// ----- Dirac matrices --------------------------------------------------------

Matrix dirac_gamma(int mu) {
    Expr O = integer(0), P = integer(1), N = integer(-1), i = S::I(), ni = mul(integer(-1), i);
    switch (mu) {
        case 0:
            return Matrix{{P, O, O, O}, {O, P, O, O}, {O, O, N, O}, {O, O, O, N}};
        case 1:
            return Matrix{{O, O, O, P}, {O, O, P, O}, {O, N, O, O}, {N, O, O, O}};
        case 2:
            return Matrix{{O, O, O, ni}, {O, O, i, O}, {O, i, O, O}, {ni, O, O, O}};
        case 3:
            return Matrix{{O, O, P, O}, {O, O, O, N}, {N, O, O, O}, {O, P, O, O}};
        case 5:  // γ⁵ = i·γ⁰γ¹γ²γ³
            return Matrix{{O, O, P, O}, {O, O, O, P}, {P, O, O, O}, {O, P, O, O}};
        default:
            throw std::invalid_argument("dirac_gamma: mu must be 0,1,2,3 or 5");
    }
}

// ----- Hydrogen atom ---------------------------------------------------------

Expr hydrogen_energy(const Expr& n, const Expr& Z) {
    // E_n = −Z²/(2n²).
    return mul(rational(-1, 2), mul(pow(Z, integer(2)), pow(n, integer(-2))));
}

Expr hydrogen_R_nl(int n, int l, const Expr& r, const Expr& Z) {
    // R_nl = √[(2Z/n)³·(n−l−1)!/(2n·(n+l)!)]·e^{−Zr/n}·(2Zr/n)ˡ·L_{n−l−1}^{2l+1}(2Zr/n)
    Expr two_Z_over_n = mul(integer(2), mul(Z, pow(integer(n), integer(-1))));  // 2Z/n
    Expr rho = mul(two_Z_over_n, r);                                            // 2Zr/n

    Expr norm_sq = mul(pow(two_Z_over_n, integer(3)),
                       mul(factorial(integer(n - l - 1)),
                           pow(mul(integer(2 * n), factorial(integer(n + l))), integer(-1))));
    Expr norm = sqrt(norm_sq);

    // Associated Laguerre L_k^α(x) = Σ_{i=0}^{k} (−1)ⁱ·C(k+α, k−i)·xⁱ/i!.
    int k = n - l - 1;
    int alpha = 2 * l + 1;
    Expr laguerre_sum = integer(0);
    for (int i = 0; i <= k; ++i) {
        Expr term = mul(binomial(integer(k + alpha), integer(k - i)),
                        mul(pow(rho, integer(i)), pow(factorial(integer(i)), integer(-1))));
        if (i % 2 == 1) term = mul(integer(-1), term);
        laguerre_sum = add(laguerre_sum, term);
    }

    Expr decay = exp(mul(integer(-1), mul(Z, mul(r, pow(integer(n), integer(-1))))));  // e^{−Zr/n}
    return mul(norm, mul(decay, mul(pow(rho, integer(l)), laguerre_sum)));
}

// ----- 1D quantum harmonic oscillator ----------------------------------------

Expr qho_energy(const Expr& n, const Expr& omega) {
    return mul(add(n, rational(1, 2)), omega);  // (n+½)ω
}

Expr qho_wavefunction(int n, const Expr& x) {
    // ψ_n(x) = (2ⁿ n!)^{−1/2}·π^{−1/4}·e^{−x²/2}·H_n(x)   (m=ω=ℏ=1).
    Expr norm = pow(mul(pow(integer(2), integer(n)), factorial(integer(n))), rational(-1, 2));
    Expr pi_factor = pow(S::Pi(), rational(-1, 4));
    Expr gaussian = exp(mul(rational(-1, 2), pow(x, integer(2))));
    Expr H = hermite(integer(n), x);
    return mul(mul(norm, pi_factor), mul(gaussian, H));
}

}  // namespace sympp::physics
