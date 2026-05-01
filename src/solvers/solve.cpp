#include <sympp/solvers/solve.hpp>

#include <stdexcept>
#include <vector>

#include <sympp/core/operators.hpp>
#include <sympp/polys/poly.hpp>

namespace sympp {

std::vector<Expr> solve(const Expr& expr, const Expr& var) {
    Poly p(expr, var);
    return p.roots();
}

std::vector<Expr> solve(const Expr& lhs, const Expr& rhs, const Expr& var) {
    return solve(lhs - rhs, var);
}

Matrix linsolve(const Matrix& A, const Matrix& b) {
    if (!A.is_square()) {
        throw std::invalid_argument("linsolve: A must be square");
    }
    if (b.rows() != A.rows()) {
        throw std::invalid_argument("linsolve: b dimension mismatch");
    }
    return A.inverse() * b;
}

}  // namespace sympp
