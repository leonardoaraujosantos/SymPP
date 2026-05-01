#include <sympp/calculus/pade.hpp>

#include <stdexcept>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

Expr pade(const Expr& expr, const Expr& var, std::size_t m, std::size_t n) {
    const std::size_t order = m + n;

    // Extract Taylor coefficients c_0, ..., c_{m+n} of expr at var = 0.
    // c_k = (1/k!) * d^k expr / d var^k  evaluated at 0.
    std::vector<Expr> c(order + 1);
    Expr current = expr;
    Expr factorial = S::One();
    for (std::size_t k = 0; k <= order; ++k) {
        c[k] = simplify(subs(current, var, S::Zero()) / factorial);
        if (k < order) {
            current = diff(current, var);
            factorial = simplify(factorial * integer(static_cast<long>(k + 1)));
        }
    }

    // Solve for Q(x) = 1 + b_1 x + ... + b_n x^n via the n×n system
    //   Σ_{k=1}^n c_{j-k} b_k = -c_j   for j = m+1, m+2, ..., m+n.
    if (n == 0) {
        // Pure polynomial case: P(x) = c_0 + c_1 x + ... + c_m x^m, Q = 1.
        Expr P = S::Zero();
        for (std::size_t j = 0; j <= m; ++j) {
            P = P + c[j] * pow(var, integer(static_cast<long>(j)));
        }
        return simplify(P);
    }

    Matrix M(n, n);
    Matrix rhs(n, 1);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 1; k <= n; ++k) {
            std::size_t idx = m + 1 + i - k;
            M.set(i, k - 1, c[idx]);
        }
        rhs.set(i, 0, mul(S::NegativeOne(), c[m + 1 + i]));
    }

    Matrix b = Matrix::zeros(n, 1);
    try {
        Matrix Minv = M.inverse();
        b = Minv * rhs;
    } catch (const std::exception&) {
        // Singular system — Padé table degeneracy. Return input unchanged.
        return expr;
    }

    // Compute a_j for j = 0 .. m:
    //   a_j = c_j + Σ_{k=1}^min(j, n) c_{j-k} b_k.
    std::vector<Expr> a(m + 1);
    for (std::size_t j = 0; j <= m; ++j) {
        a[j] = c[j];
        for (std::size_t k = 1; k <= n && k <= j; ++k) {
            a[j] = a[j] + c[j - k] * b.at(k - 1, 0);
        }
    }

    // Assemble P(x) and Q(x).
    Expr P = S::Zero();
    for (std::size_t j = 0; j <= m; ++j) {
        if (a[j] == S::Zero()) continue;
        P = P + a[j] * pow(var, integer(static_cast<long>(j)));
    }
    Expr Q = S::One();
    for (std::size_t k = 1; k <= n; ++k) {
        Q = Q + b.at(k - 1, 0) * pow(var, integer(static_cast<long>(k)));
    }
    return simplify(P / Q);
}

}  // namespace sympp
