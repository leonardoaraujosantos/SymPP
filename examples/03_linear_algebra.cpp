// 03 — Linear algebra: det, inverse, eigendecomposition, jacobian.

#include <sympp/sympp.hpp>
#include <iostream>

int main() {
    using namespace sympp;

    Matrix A = {{integer(2), integer(1)}, {integer(1), integer(3)}};

    std::cout << "A = " << A.str() << "\n";
    std::cout << "det(A)   = " << A.det()->str()   << "\n";
    std::cout << "trace(A) = " << A.trace()->str() << "\n";

    // Eigendecomposition
    auto evs = A.eigenvals();
    std::cout << "eigenvalues: ";
    for (auto& e : evs) std::cout << e->str() << " ";
    std::cout << "\n";

    // LU and Cholesky decompositions
    auto [L, U] = A.lu();
    std::cout << "LU decomposition done (L lower-tri, U upper-tri)\n";

    // Multivariate Jacobian
    auto x = symbol("x"), y = symbol("y");
    auto J = jacobian(
        {pow(x, integer(2)) + y, integer(2) * x * y},
        {x, y});
    std::cout << "jacobian shape: " << J.rows() << "x" << J.cols() << "\n";
    std::cout << "  J[0,0] = " << J.at(0, 0)->str() << "\n";
    std::cout << "  J[1,0] = " << J.at(1, 0)->str() << "\n";
}
