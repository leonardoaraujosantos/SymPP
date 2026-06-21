#include <sympp/functions/orthopolys.hpp>

#include <functional>
#include <stdexcept>

#include <sympp/core/expand.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// Validate and extract a nonnegative machine-sized degree.
[[nodiscard]] long degree(const Expr& n, const char* fn) {
    if (!n || n->type_id() != TypeId::Integer) {
        throw std::invalid_argument(std::string(fn) +
                                    ": degree must be a nonnegative integer");
    }
    const auto& z = static_cast<const Integer&>(*n);
    if (z.value() < 0 || !z.fits_long()) {
        throw std::invalid_argument(std::string(fn) +
                                    ": degree must be a nonnegative machine integer");
    }
    return z.to_long();
}

// Drive a three-term recurrence p_{k+1} = rec(k, p_k, p_{k-1}) from the two
// seed values, expanding at each step to keep the canonical polynomial form.
[[nodiscard]] Expr recur(long n, const Expr& p0, const Expr& p1,
                         const std::function<Expr(long, const Expr&, const Expr&)>& rec) {
    if (n == 0) return p0;
    if (n == 1) return p1;
    Expr prev = p0, cur = p1;
    for (long k = 1; k < n; ++k) {
        Expr next = expand(rec(k, cur, prev));
        prev = cur;
        cur = next;
    }
    return cur;
}

}  // namespace

Expr legendre(const Expr& n, const Expr& x) {
    long d = degree(n, "legendre");
    // (k+1)·P_{k+1} = (2k+1)·x·P_k − k·P_{k-1}
    return recur(d, integer(1), x, [&](long k, const Expr& pk, const Expr& pm) {
        return mul(rational(1, k + 1),
                   mul(mul(integer(2 * k + 1), x), pk) - mul(integer(k), pm));
    });
}

Expr chebyshevt(const Expr& n, const Expr& x) {
    long d = degree(n, "chebyshevt");
    // T_{k+1} = 2·x·T_k − T_{k-1}
    return recur(d, integer(1), x, [&](long, const Expr& pk, const Expr& pm) {
        return mul(mul(integer(2), x), pk) - pm;
    });
}

Expr chebyshevu(const Expr& n, const Expr& x) {
    long d = degree(n, "chebyshevu");
    // U_{k+1} = 2·x·U_k − U_{k-1}, seeds U0 = 1, U1 = 2x
    return recur(d, integer(1), mul(integer(2), x),
                 [&](long, const Expr& pk, const Expr& pm) {
                     return mul(mul(integer(2), x), pk) - pm;
                 });
}

Expr hermite(const Expr& n, const Expr& x) {
    long d = degree(n, "hermite");
    // H_{k+1} = 2·x·H_k − 2k·H_{k-1}, seeds H0 = 1, H1 = 2x
    return recur(d, integer(1), mul(integer(2), x),
                 [&](long k, const Expr& pk, const Expr& pm) {
                     return mul(mul(integer(2), x), pk) - mul(integer(2 * k), pm);
                 });
}

Expr laguerre(const Expr& n, const Expr& x) {
    long d = degree(n, "laguerre");
    // (k+1)·L_{k+1} = (2k+1−x)·L_k − k·L_{k-1}, seeds L0 = 1, L1 = 1 − x
    return recur(d, integer(1), integer(1) - x,
                 [&](long k, const Expr& pk, const Expr& pm) {
                     return mul(rational(1, k + 1),
                                mul(integer(2 * k + 1) - x, pk) - mul(integer(k), pm));
                 });
}

}  // namespace sympp
