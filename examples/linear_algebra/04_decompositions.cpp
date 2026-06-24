// Engineering Linear Algebra — 04: Matrix decompositions.
//
// Factorisations turn one hard problem into a chain of easy ones: LU for fast
// repeated solves, QR for least squares and orthogonalisation, SVD for rank,
// conditioning and data compression. Each is verified by reconstructing the
// original matrix.

#include <sympp/sympp.hpp>

#include <iostream>

namespace {
using namespace sympp;
void print(const std::string& name, const Matrix& M) {
    std::cout << name << " =\n";
    for (std::size_t i = 0; i < M.rows(); ++i) {
        std::cout << "  [ ";
        for (std::size_t j = 0; j < M.cols(); ++j)
            std::cout << simplify(M.at(i, j))->str()
                      << (j + 1 < M.cols() ? "  " : "");
        std::cout << " ]\n";
    }
}
}  // namespace

int main() {
    using namespace sympp;

    Matrix A = {{integer(4), integer(3)}, {integer(6), integer(3)}};
    print("A", A);

    std::cout << "\n=== LU decomposition (A = L U) ===\n";
    auto [L, U] = A.lu();
    print("L (lower)", L);
    print("U (upper)", U);
    print("L * U (reconstructed)", L * U);

    std::cout << "\n=== QR decomposition (A = Q R, Q orthogonal) ===\n";
    auto [Q, Rqr] = A.qr();
    print("Q", Q);
    print("R", Rqr);
    print("Q * R (reconstructed)", Q * Rqr);
    // Q^T Q = I (orthonormal columns).
    print("Q^T Q (= I)", Q.transpose() * Q);

    std::cout << "\n=== SVD (A = U S V^T) ===\n";
    auto [Usvd, Sig, V] = A.svd();
    print("singular values (diag of S)", Sig);
    print("U S V^T (reconstructed)", Usvd * Sig * V.transpose());
    return 0;
}
