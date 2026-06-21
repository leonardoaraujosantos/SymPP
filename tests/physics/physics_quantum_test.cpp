// Extended physics — angular-momentum/spin operators, Wigner 3-j/6-j and
// Clebsch–Gordan coefficients, hydrogen & QHO states, Hamiltonian, Gaussian
// optics. Closed forms are cross-checked against SymPy's sympy.physics.

#include <string>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/physics/physics.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
namespace ph = sympp::physics;

namespace {
auto& oracle() { return sympp::testing::Oracle::instance(); }

bool mat_equiv(const Matrix& a, const Matrix& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    for (std::size_t i = 0; i < a.rows(); ++i) {
        for (std::size_t j = 0; j < a.cols(); ++j) {
            if (!(simplify(a.at(i, j) - b.at(i, j)) == S::Zero())) return false;
        }
    }
    return true;
}

// Run a physics oracle op and compare its result to `mine` for symbolic equality.
void check_physics(const Expr& mine, nlohmann::json req) {
    auto r = oracle().send(req);
    REQUIRE(r.ok);
    REQUIRE(oracle().equivalent(mine->str(), r.result_str()));
}
}  // namespace

TEST_CASE("spin: angular-momentum operator algebra", "[physics][quantum][spin]") {
    // Spin-½: Jz = diag(½,−½), and the su(2) commutator [Jx,Jy] = i·Jz.
    auto half = rational(1, 2);
    auto Jz = ph::spin_jz(half);
    REQUIRE(Jz.at(0, 0) == rational(1, 2));
    REQUIRE(Jz.at(1, 1) == rational(-1, 2));

    auto Jx = ph::spin_jx(half), Jy = ph::spin_jy(half);
    REQUIRE(mat_equiv(ph::commutator(Jx, Jy), Jz.scalar_mul(S::I())));

    // Spin-1: J² = j(j+1)·I = 2·I, and [Jx,Jy] = i·Jz holds too.
    auto one = integer(1);
    REQUIRE(mat_equiv(ph::spin_j2(one), Matrix::identity(3).scalar_mul(integer(2))));
    REQUIRE(mat_equiv(ph::commutator(ph::spin_jx(one), ph::spin_jy(one)),
                      ph::spin_jz(one).scalar_mul(S::I())));

    // J₊|m=−1⟩ = √2|0⟩ for spin-1 (descending basis: row 1, col 2).
    REQUIRE(ph::spin_jplus(one).at(1, 2) == sqrt(integer(2)));
}

TEST_CASE("wigner: 3-j symbols match SymPy", "[physics][quantum][wigner]") {
    auto t = [](int n, int d) { return d == 1 ? integer(n) : rational(n, d); };
    // Integer arguments.
    check_physics(ph::wigner_3j(integer(2), integer(6), integer(4), integer(0), integer(0),
                                integer(0)),
                  {{"op", "physics"}, {"fn", "wigner_3j"}, {"j1", "2"}, {"j2", "6"}, {"j3", "4"},
                   {"m1", "0"}, {"m2", "0"}, {"m3", "0"}});
    // Half-integer arguments.
    check_physics(ph::wigner_3j(t(1, 2), t(1, 2), integer(1), t(1, 2), t(-1, 2), integer(0)),
                  {{"op", "physics"}, {"fn", "wigner_3j"}, {"j1", "1/2"}, {"j2", "1/2"},
                   {"j3", "1"}, {"m1", "1/2"}, {"m2", "-1/2"}, {"m3", "0"}});
    check_physics(ph::wigner_3j(integer(1), integer(1), integer(2), integer(1), integer(-1),
                                integer(0)),
                  {{"op", "physics"}, {"fn", "wigner_3j"}, {"j1", "1"}, {"j2", "1"}, {"j3", "2"},
                   {"m1", "1"}, {"m2", "-1"}, {"m3", "0"}});
    // Selection-rule violation → 0.
    REQUIRE(ph::wigner_3j(integer(1), integer(1), integer(1), integer(0), integer(0),
                          integer(0)) == S::Zero());
}

TEST_CASE("wigner: 6-j symbols match SymPy", "[physics][quantum][wigner]") {
    check_physics(ph::wigner_6j(integer(1), integer(1), integer(1), integer(1), integer(1),
                                integer(1)),
                  {{"op", "physics"}, {"fn", "wigner_6j"}, {"j1", "1"}, {"j2", "1"}, {"j3", "1"},
                   {"j4", "1"}, {"j5", "1"}, {"j6", "1"}});
    check_physics(ph::wigner_6j(integer(3), integer(3), integer(3), integer(3), integer(3),
                                integer(3)),
                  {{"op", "physics"}, {"fn", "wigner_6j"}, {"j1", "3"}, {"j2", "3"}, {"j3", "3"},
                   {"j4", "3"}, {"j5", "3"}, {"j6", "3"}});
}

TEST_CASE("wigner: Clebsch–Gordan match SymPy", "[physics][quantum][wigner]") {
    auto h = rational(1, 2);
    check_physics(ph::clebsch_gordan(h, h, integer(1), h, h, integer(1)),
                  {{"op", "physics"}, {"fn", "clebsch_gordan"}, {"j1", "1/2"}, {"j2", "1/2"},
                   {"j3", "1"}, {"m1", "1/2"}, {"m2", "1/2"}, {"m3", "1"}});
    check_physics(ph::clebsch_gordan(h, h, integer(1), h, rational(-1, 2), integer(0)),
                  {{"op", "physics"}, {"fn", "clebsch_gordan"}, {"j1", "1/2"}, {"j2", "1/2"},
                   {"j3", "1"}, {"m1", "1/2"}, {"m2", "-1/2"}, {"m3", "0"}});
}

TEST_CASE("hydrogen: energies and radial wavefunctions", "[physics][atomic]") {
    check_physics(ph::hydrogen_energy(integer(2), integer(1)),
                  {{"op", "physics"}, {"fn", "hydrogen_E"}, {"n", "2"}, {"Z", "1"}});
    check_physics(ph::hydrogen_energy(integer(3), integer(2)),
                  {{"op", "physics"}, {"fn", "hydrogen_E"}, {"n", "3"}, {"Z", "2"}});

    auto r = symbol("r");
    check_physics(ph::hydrogen_R_nl(1, 0, r, integer(1)),
                  {{"op", "physics"}, {"fn", "hydrogen_R"}, {"n", "1"}, {"l", "0"}, {"Z", "1"}});
    check_physics(ph::hydrogen_R_nl(2, 1, r, integer(1)),
                  {{"op", "physics"}, {"fn", "hydrogen_R"}, {"n", "2"}, {"l", "1"}, {"Z", "1"}});
    check_physics(ph::hydrogen_R_nl(2, 0, r, integer(1)),
                  {{"op", "physics"}, {"fn", "hydrogen_R"}, {"n", "2"}, {"l", "0"}, {"Z", "1"}});
}

TEST_CASE("qho: energies and wavefunctions", "[physics][atomic]") {
    auto w = symbol("omega");
    check_physics(ph::qho_energy(integer(3), w),
                  {{"op", "physics"}, {"fn", "qho_E"}, {"n", "3"}});

    auto x = symbol("x");
    check_physics(ph::qho_wavefunction(0, x),
                  {{"op", "physics"}, {"fn", "qho_psi"}, {"n", "0"}});
    check_physics(ph::qho_wavefunction(2, x),
                  {{"op", "physics"}, {"fn", "qho_psi"}, {"n", "2"}});
}

TEST_CASE("mechanics: Hamiltonian via Legendre transform", "[physics][mechanics]") {
    auto m = symbol("m"), v = symbol("v");
    // L = ½mv² → H = p·v − L = mv² − ½mv² = ½mv².
    Expr L = mul(rational(1, 2), mul(m, pow(v, integer(2))));
    REQUIRE(oracle().equivalent(ph::hamiltonian(L, v)->str(), L->str()));
}

TEST_CASE("optics: refraction and Gaussian-beam propagation", "[physics][optics]") {
    auto n1 = symbol("n1"), n2 = symbol("n2");
    // A flat interface preserves the position, scaling the angle by n1/n2.
    REQUIRE(simplify(ph::abcd_flat_interface(n1, n2).det() - mul(n1, pow(n2, integer(-1)))) ==
            S::Zero());

    // Free-space propagation by d sends q → q + d.
    auto q = symbol("q"), d = symbol("d");
    Expr out = ph::gaussian_beam_q(ph::abcd_free_space(d), q);
    REQUIRE(oracle().equivalent(out->str(), (q + d)->str()));
}
