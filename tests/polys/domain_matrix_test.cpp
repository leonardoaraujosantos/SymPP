// DomainMatrix — fraction-free det / rank / rref cross-checked against Matrix.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/polys/domain_matrix.hpp>
#include <sympp/simplify/simplify.hpp>

using namespace sympp;

namespace {
bool det_matches(const Matrix& m) {
    auto dm = DomainMatrix::from_matrix(m);
    return simplify(dm.det() - m.det()) == S::Zero();
}
}  // namespace

TEST_CASE("DomainMatrix determinant matches Matrix (fraction-free)", "[domainmatrix]") {
    REQUIRE(det_matches(Matrix{{integer(1), integer(2)}, {integer(3), integer(4)}}));
    REQUIRE(det_matches(Matrix{{integer(2), integer(0), integer(1)},
                               {integer(3), integer(1), integer(4)},
                               {integer(1), integer(5), integer(9)}}));
    // 4×4 — Bareiss keeps it integer-only.
    REQUIRE(det_matches(Matrix{{integer(1), integer(2), integer(3), integer(4)},
                               {integer(5), integer(6), integer(0), integer(8)},
                               {integer(9), integer(1), integer(2), integer(3)},
                               {integer(4), integer(5), integer(6), integer(0)}}));
    // Rational entries (domain ℚ).
    REQUIRE(det_matches(Matrix{{rational(1, 2), integer(1)}, {integer(3), rational(2, 3)}}));
    // Singular matrix → 0.
    auto dm = DomainMatrix::from_matrix(
        Matrix{{integer(1), integer(2)}, {integer(2), integer(4)}});
    REQUIRE(dm.det() == integer(0));
}

TEST_CASE("DomainMatrix rank and domain inference", "[domainmatrix]") {
    Matrix full = {{integer(1), integer(0)}, {integer(0), integer(2)}};
    Matrix rank1 = {{integer(1), integer(2)}, {integer(2), integer(4)}};
    REQUIRE(DomainMatrix::from_matrix(full).rank() == full.rank());
    REQUIRE(DomainMatrix::from_matrix(rank1).rank() == rank1.rank());
    REQUIRE(DomainMatrix::from_matrix(rank1).rank() == 1);

    REQUIRE(DomainMatrix::from_matrix(full).domain() == DomainMatrix::Domain::ZZ);
    Matrix q = {{rational(1, 2), integer(1)}, {integer(0), integer(1)}};
    REQUIRE(DomainMatrix::from_matrix(q).domain() == DomainMatrix::Domain::QQ);
}

TEST_CASE("DomainMatrix multiply and round-trip", "[domainmatrix]") {
    Matrix a = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    Matrix b = {{integer(0), integer(1)}, {integer(1), integer(0)}};
    auto prod = (DomainMatrix::from_matrix(a) * DomainMatrix::from_matrix(b)).to_matrix();
    Matrix expected = a * b;
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            REQUIRE(simplify(prod.at(i, j) - expected.at(i, j)) == S::Zero());
        }
    }
    // round-trip preserves entries
    REQUIRE(simplify(DomainMatrix::from_matrix(a).to_matrix().at(1, 1) - integer(4)) == S::Zero());
}
