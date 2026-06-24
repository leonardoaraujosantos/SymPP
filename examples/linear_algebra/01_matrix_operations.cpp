// Engineering Linear Algebra — 01: Matrix operations.
//
// Matrices represent linear maps, networks and state-space models. Covers the
// core algebra: products, transpose, determinant, trace, inverse and powers,
// plus the inverse-times-matrix identity A·A^-1 = I.

#include <sympp/sympp.hpp>

#include <iostream>

namespace {
using namespace sympp;
void print(const std::string& name, const Matrix& M) {
    std::cout << name << " =\n";
    for (std::size_t i = 0; i < M.rows(); ++i) {
        std::cout << "  [ ";
        for (std::size_t j = 0; j < M.cols(); ++j)
            std::cout << M.at(i, j)->str() << (j + 1 < M.cols() ? "  " : "");
        std::cout << " ]\n";
    }
}
}  // namespace

int main() {
    using namespace sympp;

    Matrix A = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    Matrix B = {{integer(0), integer(1)}, {integer(1), integer(0)}};
    print("A", A);
    print("B", B);

    std::cout << "\n=== Sum & product ===\n";
    print("A + B", A + B);
    print("A * B", A * B);
    print("B * A (note: A*B != B*A)", B * A);

    std::cout << "\n=== Transpose ===\n";
    print("A^T", A.transpose());

    std::cout << "\n=== Determinant, trace ===\n";
    std::cout << "det(A)   = " << A.det()->str() << "\n";
    std::cout << "trace(A) = " << A.trace()->str() << "\n";

    std::cout << "\n=== Inverse ===\n";
    Matrix Ainv = A.inverse();
    print("A^-1", Ainv);
    // Verify A * A^-1 = I.
    print("A * A^-1 (= I)", A * Ainv);

    std::cout << "\n=== Matrix powers (discrete-time state update) ===\n";
    // A 90-degree rotation matrix: R^4 = I.
    Matrix R = {{integer(0), integer(-1)}, {integer(1), integer(0)}};
    print("R", R);
    print("R^2", R * R);
    print("R^4 (= I)", R * R * R * R);
    return 0;
}
