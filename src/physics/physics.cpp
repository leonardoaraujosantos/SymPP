#include <sympp/physics/physics.hpp>

#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/functions/miscellaneous.hpp>

namespace sympp::physics {

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

Expr conjugate_momentum(const Expr& lagrangian, const Expr& velocity) {
    return diff(lagrangian, velocity);
}

}  // namespace sympp::physics
