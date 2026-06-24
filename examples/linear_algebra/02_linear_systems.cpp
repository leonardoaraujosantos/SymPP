// Engineering Linear Algebra — 02: Linear systems.
//
// Solving A·x = b is the heart of statics, circuit nodal analysis and finite
// elements. Covers a determined system, rank / row-reduction, the null space
// (homogeneous solutions) and a consistency check for an over/under-determined
// system.

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

    std::cout << "=== Determined system  A x = b ===\n";
    // Node equations:  3x + 2y +  z = 1
    //                  2x + 3y +  z = 2
    //                   x +  y + 4z = 3
    Matrix A = {{integer(3), integer(2), integer(1)},
                {integer(2), integer(3), integer(1)},
                {integer(1), integer(1), integer(4)}};
    Matrix b = {{integer(1)}, {integer(2)}, {integer(3)}};
    print("A", A);
    std::cout << "det(A) = " << A.det()->str() << "  (nonzero -> unique soln)\n";
    Matrix x = linsolve(A, b);
    std::cout << "solution x = (" << x.at(0, 0)->str() << ", " << x.at(1, 0)->str()
              << ", " << x.at(2, 0)->str() << ")\n";

    std::cout << "\n=== Rank & reduced row echelon form ===\n";
    std::cout << "rank(A) = " << A.rank() << "\n";
    print("rref(A)", A.rref().first);

    std::cout << "\n=== Null space (homogeneous solutions) ===\n";
    // A rank-deficient matrix: rows 2 and 3 are multiples of row 1.
    Matrix S = {{integer(1), integer(2), integer(3)},
                {integer(2), integer(4), integer(6)},
                {integer(3), integer(6), integer(9)}};
    print("S (rank-deficient)", S);
    std::cout << "rank(S) = " << S.rank() << "\n";
    auto ns = S.nullspace();
    std::cout << "null space dimension = " << ns.size() << "\n";
    for (std::size_t i = 0; i < ns.size(); ++i)
        print("  null vector " + std::to_string(i + 1), ns[i]);
    return 0;
}
