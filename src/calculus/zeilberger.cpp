#include "zeilberger.hpp"

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include <sympp/calculus/product.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp::detail {

namespace {

// Number of integer sample points and the largest recurrence order / coeff
// degree we look for. Generous enough for the standard binomial sums, bounded
// so discovery is cheap and always terminates.
constexpr int kSamples = 18;
constexpr int kMaxOrder = 2;
constexpr int kMaxDeg = 3;

[[nodiscard]] std::optional<mpq_class> expr_to_q(const Expr& e) {
    Expr s = simplify(e);
    if (s->type_id() == TypeId::Integer) {
        return mpq_class(static_cast<const Integer&>(*s).value());
    }
    if (s->type_id() == TypeId::Rational) {
        return static_cast<const Rational&>(*s).value();
    }
    return std::nullopt;
}

[[nodiscard]] Expr q_to_expr(const mpq_class& q) {
    return rational(q.get_num(), q.get_den());
}

// Evaluate S(n0) = Σ_{k=lo0}^{hi0} F(n0, k) as an exact rational, or nullopt if
// a term is not numeric / the range is not finite.
[[nodiscard]] std::optional<mpq_class> eval_S(const Expr& F, const Expr& k,
                                              const Expr& lo, const Expr& hi,
                                              const Expr& n, long n0) {
    auto lo_q = expr_to_q(subs(lo, n, integer(n0)));
    auto hi_q = expr_to_q(subs(hi, n, integer(n0)));
    if (!lo_q || !hi_q || lo_q->get_den() != 1 || hi_q->get_den() != 1) {
        return std::nullopt;
    }
    long k_lo = lo_q->get_num().get_si();
    long k_hi = hi_q->get_num().get_si();
    if (k_hi - k_lo > 4000) return std::nullopt;  // guard
    Expr Fn = subs(F, n, integer(n0));
    mpq_class acc(0);
    for (long kk = k_lo; kk <= k_hi; ++kk) {
        auto t = expr_to_q(subs(Fn, k, integer(kk)));
        if (!t) return std::nullopt;
        acc += *t;
    }
    return acc;
}

// One basis vector of the right null space of the exact-rational matrix M
// (rows × cols), or nullopt when M has full column rank (only the zero vector).
[[nodiscard]] std::optional<std::vector<mpq_class>> nullspace_vector(
    std::vector<std::vector<mpq_class>> M, std::size_t cols) {
    const std::size_t rows = M.size();
    std::vector<long> pivot_col(rows, -1);
    std::vector<bool> is_pivot(cols, false);
    std::size_t r = 0;
    for (std::size_t c = 0; c < cols && r < rows; ++c) {
        std::size_t sel = r;
        bool found = false;
        for (std::size_t i = r; i < rows; ++i) {
            if (M[i][c] != 0) { sel = i; found = true; break; }
        }
        if (!found) continue;
        std::swap(M[r], M[sel]);
        mpq_class piv = M[r][c];
        for (std::size_t j = c; j < cols; ++j) M[r][j] /= piv;
        for (std::size_t i = 0; i < rows; ++i) {
            if (i != r && M[i][c] != 0) {
                mpq_class f = M[i][c];
                for (std::size_t j = c; j < cols; ++j) M[i][j] -= f * M[r][j];
            }
        }
        pivot_col[r] = static_cast<long>(c);
        is_pivot[c] = true;
        ++r;
    }
    // First free column → a null-space basis vector.
    long free_c = -1;
    for (std::size_t c = 0; c < cols; ++c) {
        if (!is_pivot[c]) { free_c = static_cast<long>(c); break; }
    }
    if (free_c < 0) return std::nullopt;  // full column rank
    std::vector<mpq_class> v(cols, mpq_class(0));
    v[static_cast<std::size_t>(free_c)] = 1;
    for (std::size_t i = 0; i < r; ++i) {
        long pc = pivot_col[i];
        if (pc >= 0) v[static_cast<std::size_t>(pc)] = -M[i][static_cast<std::size_t>(free_c)];
    }
    return v;
}

// Build the polynomial a_j(n) = Σ_d coeff[j*(D+1)+d] · n^d as an Expr in `var`.
[[nodiscard]] Expr build_poly(const std::vector<mpq_class>& coeff, int j, int D,
                              const Expr& var) {
    std::vector<Expr> terms;
    for (int d = 0; d <= D; ++d) {
        const mpq_class& c = coeff[static_cast<std::size_t>(j * (D + 1) + d)];
        if (c == 0) continue;
        terms.push_back(mul(q_to_expr(c), pow(var, integer(d))));
    }
    return terms.empty() ? Expr{S::Zero()} : add(std::move(terms));
}

}  // namespace

std::optional<Expr> zeilberger_summation(const Expr& F, const Expr& k,
                                         const Expr& lo, const Expr& hi) {
    if (!F || !has(F, k)) return std::nullopt;
    if (hi->type_id() == TypeId::Infinity) return std::nullopt;  // finite hi only

    // Identify the single parameter n: free symbols of F and hi, minus k.
    ExprSet syms = free_symbols(F);
    for (const auto& s : free_symbols(hi)) syms.insert(s);
    Expr n;
    for (const auto& s : syms) {
        if (s == k) continue;
        if (n) return std::nullopt;  // more than one extra parameter
        n = s;
    }
    if (!n) return std::nullopt;  // no parameter → nothing to recurse on

    // Sample S(n0). Start at 1; n0 = 0 is often an exceptional initial value.
    std::vector<std::optional<mpq_class>> S(static_cast<std::size_t>(kSamples) + 1);
    int n_lo = 1;
    for (int n0 = n_lo; n0 <= kSamples; ++n0) {
        S[static_cast<std::size_t>(n0)] = eval_S(F, k, lo, hi, n, n0);
    }

    auto m = symbol("__zeil_m");

    // Search recurrence order J and coefficient degree D, smallest first.
    for (int J = 1; J <= kMaxOrder; ++J) {
        for (int D = 0; D <= kMaxDeg; ++D) {
            const std::size_t cols =
                static_cast<std::size_t>((J + 1) * (D + 1));
            // Rows: every base n0 with S(n0 .. n0+J) all defined.
            std::vector<std::vector<mpq_class>> M;
            for (int n0 = n_lo; n0 + J <= kSamples; ++n0) {
                bool ok = true;
                for (int j = 0; j <= J; ++j) {
                    if (!S[static_cast<std::size_t>(n0 + j)]) { ok = false; break; }
                }
                if (!ok) continue;
                std::vector<mpq_class> row(cols);
                for (int j = 0; j <= J; ++j) {
                    mpq_class npow(1);
                    const mpq_class& sv = *S[static_cast<std::size_t>(n0 + j)];
                    for (int d = 0; d <= D; ++d) {
                        row[static_cast<std::size_t>(j * (D + 1) + d)] = npow * sv;
                        npow *= n0;
                    }
                }
                M.push_back(std::move(row));
            }
            // Need a comfortably overdetermined system so a found null vector is
            // already verified on many points.
            if (M.size() < cols + 3) continue;

            auto v = nullspace_vector(M, cols);
            if (!v) continue;

            // Only first-order recurrences are solved to a closed form here.
            if (J != 1) return std::nullopt;  // higher order → leave unevaluated

            Expr a0 = build_poly(*v, 0, D, m);
            Expr a1 = build_poly(*v, 1, D, m);
            if (a1 == S::Zero()) continue;

            // ρ(m) = S(m+1)/S(m) = −a0(m)/a1(m). Choose the base n_b: the first
            // sampled n with S(n_b) ≠ 0.
            int n_b = -1;
            for (int n0 = n_lo; n0 <= kSamples; ++n0) {
                if (S[static_cast<std::size_t>(n0)]
                    && *S[static_cast<std::size_t>(n0)] != 0) { n_b = n0; break; }
            }
            if (n_b < 0) continue;

            // Keep ρ as a single reduced fraction p(m)/q(m) — cancel() lowers
            // it to lowest terms without expanding into an Add (which product()
            // could not factor).
            Expr rho = cancel(mul(mul(integer(-1), a0), pow(a1, integer(-1))), m);
            Expr prod = product(rho, m, integer(n_b), add(n, integer(-1)));
            if (prod->str().find("Product(") != std::string::npos) continue;
            Expr closed = simplify(
                mul(q_to_expr(*S[static_cast<std::size_t>(n_b)]), prod));

            // Verify on held-out points (every sampled n0 > n_b).
            int checked = 0;
            bool good = true;
            for (int n0 = n_b + 1; n0 <= kSamples && good; ++n0) {
                if (!S[static_cast<std::size_t>(n0)]) continue;
                auto cv = expr_to_q(subs(closed, n, integer(n0)));
                if (!cv || *cv != *S[static_cast<std::size_t>(n0)]) good = false;
                else ++checked;
            }
            if (good && checked >= 4) return closed;
            return std::nullopt;  // recurrence found but closed form unverified
        }
    }
    return std::nullopt;
}

}  // namespace sympp::detail
