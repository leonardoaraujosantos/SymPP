// Engineering Linear Algebra — 03: Eigenvalues & vibration modes.
//
// Eigen-analysis gives natural frequencies and mode shapes of vibrating
// structures, principal axes of stress/inertia, and the stability of dynamic
// systems. Covers eigenvalues, eigenvectors, the characteristic polynomial and
// a 2-DOF spring–mass modal analysis.

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
    auto lam = symbol("lambda");

    std::cout << "=== Eigenvalues & eigenvectors ===\n";
    // Symmetric matrix -> real eigenvalues, orthogonal modes.
    Matrix A = {{integer(2), integer(1)}, {integer(1), integer(2)}};
    print("A", A);
    std::cout << "characteristic poly = " << A.charpoly(lam)->str() << "\n";
    std::cout << "eigenvalues:";
    for (const auto& e : A.eigenvals()) std::cout << " " << e->str();
    std::cout << "\n";
    for (const auto& [val, basis] : A.eigenvects()) {
        std::cout << "  lambda = " << val->str() << ", eigenvector = ("
                  << basis[0].at(0, 0)->str() << ", "
                  << basis[0].at(1, 0)->str() << ")\n";
    }

    std::cout << "\n=== 2-DOF spring-mass modal analysis ===\n";
    // Two unit masses, three identical springs (k=1): the classic
    // dynamical matrix K = [[2,-1],[-1,2]]. Eigenvalues are omega^2.
    Matrix K = {{integer(2), integer(-1)}, {integer(-1), integer(2)}};
    print("stiffness/mass matrix K", K);
    std::cout << "modal omega^2 (eigenvalues):";
    for (const auto& e : K.eigenvals()) std::cout << " " << e->str();
    std::cout << "\n";
    std::cout << "mode shapes:\n";
    for (const auto& [val, basis] : K.eigenvects()) {
        std::cout << "  omega^2 = " << val->str() << ":  ("
                  << basis[0].at(0, 0)->str() << ", "
                  << basis[0].at(1, 0)->str() << ")"
                  << (val->str() == "1" ? "  in-phase\n" : "  out-of-phase\n");
    }

    std::cout << "\n=== Diagonalizability ===\n";
    std::cout << "A diagonalizable? "
              << (A.is_diagonalizable() ? "yes" : "no") << "\n";
    return 0;
}
