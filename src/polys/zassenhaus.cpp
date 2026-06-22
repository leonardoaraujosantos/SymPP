#include <sympp/polys/zassenhaus.hpp>

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <vector>

#include <gmpxx.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/polys/poly.hpp>

namespace sympp {

namespace {

using Z = std::vector<mpz_class>;  // dense polynomial, coeff[0] = constant term

void trim(Z& a) {
    while (a.size() > 1 && a.back() == 0) a.pop_back();
}
[[nodiscard]] long degree(const Z& a) { return static_cast<long>(a.size()) - 1; }
[[nodiscard]] bool is_zero(const Z& a) { return a.size() == 1 && a[0] == 0; }

// --- arithmetic modulo m -----------------------------------------------------
[[nodiscard]] mpz_class mod(const mpz_class& x, const mpz_class& m) {
    mpz_class r = x % m;
    if (r < 0) r += m;
    return r;
}
[[nodiscard]] Z modp(Z a, const mpz_class& m) {
    for (auto& c : a) c = mod(c, m);
    trim(a);
    return a;
}
[[nodiscard]] mpz_class inv_mod(const mpz_class& a, const mpz_class& m) {
    mpz_class r;
    if (mpz_invert(r.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t()) == 0) {
        throw std::runtime_error("zassenhaus: non-invertible coefficient");
    }
    return r;
}
[[nodiscard]] Z mul_mod(const Z& a, const Z& b, const mpz_class& m) {
    if (is_zero(a) || is_zero(b)) return {mpz_class(0)};
    Z r(a.size() + b.size() - 1, mpz_class(0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j) r[i + j] += a[i] * b[j];
    return modp(r, m);
}
[[nodiscard]] Z sub_mod(const Z& a, const Z& b, const mpz_class& m) {
    Z r(std::max(a.size(), b.size()), mpz_class(0));
    for (std::size_t i = 0; i < a.size(); ++i) r[i] += a[i];
    for (std::size_t i = 0; i < b.size(); ++i) r[i] -= b[i];
    return modp(r, m);
}
[[nodiscard]] Z rem_mod(Z a, const Z& b, const mpz_class& m) {
    a = modp(std::move(a), m);
    Z bb = modp(b, m);
    if (is_zero(bb)) return a;
    mpz_class linv = inv_mod(bb.back(), m);
    while (!is_zero(a) && degree(a) >= degree(bb)) {
        long shift = degree(a) - degree(bb);
        mpz_class factor = mod(a.back() * linv, m);
        for (std::size_t i = 0; i < bb.size(); ++i) {
            std::size_t idx = static_cast<std::size_t>(shift) + i;
            a[idx] = mod(a[idx] - factor * bb[i], m);
        }
        trim(a);
    }
    return a;
}
[[nodiscard]] Z div_mod(const Z& a, const Z& b, const mpz_class& m) {  // exact-ish quotient
    Z rem = modp(a, m), bb = modp(b, m);
    Z quo(static_cast<std::size_t>(std::max<long>(0, degree(a) - degree(bb)) + 1), mpz_class(0));
    mpz_class linv = inv_mod(bb.back(), m);
    while (!is_zero(rem) && degree(rem) >= degree(bb)) {
        long shift = degree(rem) - degree(bb);
        mpz_class co = mod(rem.back() * linv, m);
        quo[static_cast<std::size_t>(shift)] = co;
        for (std::size_t i = 0; i < bb.size(); ++i) {
            std::size_t idx = static_cast<std::size_t>(shift) + i;
            rem[idx] = mod(rem[idx] - co * bb[i], m);
        }
        trim(rem);
    }
    trim(quo);
    return quo;
}
[[nodiscard]] Z gcd_mod(Z a, Z b, const mpz_class& p) {
    a = modp(std::move(a), p);
    b = modp(std::move(b), p);
    while (!is_zero(b)) {
        Z r = rem_mod(a, b, p);
        a = std::move(b);
        b = std::move(r);
    }
    if (!is_zero(a)) {
        mpz_class linv = inv_mod(a.back(), p);
        for (auto& c : a) c = mod(c * linv, p);
    }
    return a;
}
[[nodiscard]] Z monic_mod(Z a, const mpz_class& p) {
    a = modp(std::move(a), p);
    mpz_class linv = inv_mod(a.back(), p);
    for (auto& c : a) c = mod(c * linv, p);
    return a;
}
[[nodiscard]] Z deriv(const Z& a) {
    if (a.size() <= 1) return {mpz_class(0)};
    Z d(a.size() - 1);
    for (std::size_t i = 1; i < a.size(); ++i) d[i - 1] = a[i] * static_cast<long>(i);
    return d;
}

// --- integer (non-modular) polynomial arithmetic, used by Hensel lifting ------
[[nodiscard]] Z mul_z(const Z& a, const Z& b) {
    if (is_zero(a) || is_zero(b)) return {mpz_class(0)};
    Z r(a.size() + b.size() - 1, mpz_class(0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j) r[i + j] += a[i] * b[j];
    trim(r);
    return r;
}
[[nodiscard]] Z sub_z(const Z& a, const Z& b) {
    Z r(std::max(a.size(), b.size()), mpz_class(0));
    for (std::size_t i = 0; i < a.size(); ++i) r[i] += a[i];
    for (std::size_t i = 0; i < b.size(); ++i) r[i] -= b[i];
    trim(r);
    return r;
}
// Reduce coefficients to the symmetric residue system modulo m: (−m/2, m/2].
[[nodiscard]] Z trunc_sym(Z a, const mpz_class& m) {
    for (auto& c : a) {
        c = mod(c, m);
        if (2 * c > m) c -= m;
    }
    trim(a);
    return a;
}

// Extended Euclid in F_p[x]: returns g = gcd and s, t with s·a + t·b ≡ g (mod p).
struct GcdEx {
    Z g, s, t;
};
[[nodiscard]] GcdEx gcdex_mod(Z a, Z b, const mpz_class& p) {
    Z old_r = modp(std::move(a), p), r = modp(std::move(b), p);
    Z old_s = {mpz_class(1)}, s = {mpz_class(0)};
    Z old_t = {mpz_class(0)}, t = {mpz_class(1)};
    while (!is_zero(r)) {
        Z q = div_mod(old_r, r, p);
        Z nr = sub_mod(old_r, mul_mod(q, r, p), p);
        old_r = std::move(r);
        r = std::move(nr);
        Z ns = sub_mod(old_s, mul_mod(q, s, p), p);
        old_s = std::move(s);
        s = std::move(ns);
        Z nt = sub_mod(old_t, mul_mod(q, t, p), p);
        old_t = std::move(t);
        t = std::move(nt);
    }
    return {old_r, old_s, old_t};
}

// Inverse of `a` in F_p[x]/(u) (assumes gcd(a,u)=1 mod p).
[[nodiscard]] Z inverse_mod_poly(const Z& a, const Z& u, const mpz_class& p) {
    GcdEx g = gcdex_mod(a, u, p);          // g.g is a nonzero constant
    mpz_class cinv = inv_mod(g.g[0], p);   // make it 1
    Z inv = g.s;
    for (auto& c : inv) c = mod(c * cinv, p);
    return rem_mod(inv, u, p);
}

// Linear multifactor Hensel lifting of a MONIC target.
// Given F monic with F ≡ ∏ uᵢ (mod p), the uᵢ monic and pairwise coprime mod p,
// return monic Uᵢ ≡ uᵢ (mod p) with ∏ Uᵢ ≡ F (mod pˡ).
[[nodiscard]] std::vector<Z> hensel_lift_monic(const Z& F, const std::vector<Z>& u,
                                               const mpz_class& p, long l) {
    std::size_t r = u.size();
    // Cofactor inverses cᵢ = (∏_{k≠i} u_k)⁻¹ mod uᵢ, computed once over F_p.
    std::vector<Z> c(r);
    for (std::size_t i = 0; i < r; ++i) {
        Z w = {mpz_class(1)};
        for (std::size_t k = 0; k < r; ++k)
            if (k != i) w = rem_mod(mul_mod(w, u[k], p), u[i], p);
        c[i] = inverse_mod_poly(w, u[i], p);
    }
    std::vector<Z> U(u.begin(), u.end());  // current lift, integer coefficients
    mpz_class pj = p;                       // pʲ at step j
    for (long step = 1; step < l; ++step) {
        Z prod = {mpz_class(1)};
        for (const auto& Ui : U) prod = mul_z(prod, Ui);
        Z err = sub_z(F, prod);             // divisible by pʲ
        for (auto& e : err) e /= pj;        // E_full / pʲ
        Z E = modp(err, p);
        if (is_zero(E)) {                   // already exact to higher precision
            pj *= p;
            continue;
        }
        for (std::size_t i = 0; i < r; ++i) {
            Z delta = rem_mod(mul_mod(E, c[i], p), u[i], p);  // deg < deg uᵢ
            for (std::size_t k = 0; k < delta.size(); ++k) {
                if (k >= U[i].size()) U[i].resize(k + 1, mpz_class(0));
                U[i][k] += pj * delta[k];
            }
        }
        pj *= p;
    }
    mpz_class M = pj;  // = pˡ
    for (auto& Ui : U) Ui = trunc_sym(std::move(Ui), M);
    return U;
}

// --- Berlekamp factorization of a monic squarefree poly mod prime p ----------
[[nodiscard]] Z x_pow_p(const Z& f, const mpz_class& p) {
    Z base = {mpz_class(0), mpz_class(1)}, result = {mpz_class(1)};
    mpz_class e = p;
    while (e > 0) {
        if (mpz_odd_p(e.get_mpz_t())) result = rem_mod(mul_mod(result, base, p), f, p);
        e /= 2;
        if (e > 0) base = rem_mod(mul_mod(base, base, p), f, p);
    }
    return result;
}

[[nodiscard]] std::vector<Z> berlekamp(const Z& f, const mpz_class& p) {
    long n = degree(f);
    if (n <= 1) return {f};
    Z xp = x_pow_p(f, p);
    auto N = static_cast<std::size_t>(n);
    // Q[i][j] = coefficient j of (x^p)^i mod f.
    std::vector<std::vector<mpz_class>> Q(N, std::vector<mpz_class>(N, 0));
    Z cur = {mpz_class(1)};
    for (std::size_t i = 0; i < N; ++i) {
        for (long j = 0; j <= degree(cur) && static_cast<std::size_t>(j) < N; ++j)
            Q[i][static_cast<std::size_t>(j)] = cur[static_cast<std::size_t>(j)];
        cur = rem_mod(mul_mod(cur, xp, p), f, p);
    }
    // Solve (Qᵀ − I) v = 0 over F_p.
    std::vector<std::vector<mpz_class>> M(N, std::vector<mpz_class>(N, 0));
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            M[i][j] = mod(Q[j][i] - (i == j ? 1 : 0), p);

    std::vector<long> pivot_col(N, -1);
    std::vector<bool> is_pivot(N, false);
    std::size_t row = 0;
    for (std::size_t col = 0; col < N && row < N; ++col) {
        std::size_t sel = N;
        for (std::size_t r = row; r < N; ++r) if (M[r][col] != 0) { sel = r; break; }
        if (sel == N) continue;
        std::swap(M[sel], M[row]);
        mpz_class inv = inv_mod(M[row][col], p);
        for (std::size_t j = 0; j < N; ++j) M[row][j] = mod(M[row][j] * inv, p);
        for (std::size_t r = 0; r < N; ++r) {
            if (r == row || M[r][col] == 0) continue;
            mpz_class fac = M[r][col];
            for (std::size_t j = 0; j < N; ++j) M[r][j] = mod(M[r][j] - fac * M[row][j], p);
        }
        pivot_col[row] = static_cast<long>(col);
        is_pivot[col] = true;
        ++row;
    }
    std::vector<Z> basis;
    for (std::size_t col = 0; col < N; ++col) {
        if (is_pivot[col]) continue;
        Z v(N, mpz_class(0));
        v[col] = 1;
        for (std::size_t r = 0; r < row; ++r) {
            long pc = pivot_col[r];
            if (pc >= 0) v[static_cast<std::size_t>(pc)] = mod(-M[r][col], p);
        }
        trim(v);
        basis.push_back(v);
    }

    std::size_t r_count = basis.size();
    std::vector<Z> factors = {monic_mod(f, p)};
    for (std::size_t bi = 0; bi < basis.size() && factors.size() < r_count; ++bi) {
        const Z& v = basis[bi];
        if (degree(v) < 1) continue;
        std::vector<Z> next;
        for (const Z& g : factors) {
            if (degree(g) <= 1) { next.push_back(g); continue; }
            // Partition g by the value of v on each root: collect gcd(g, v−s)
            // over all s ∈ F_p (this fully splits g by this basis element).
            Z rest = g;
            for (mpz_class s = 0; s < p && degree(rest) > 0; ++s) {
                Z d = gcd_mod(rest, sub_mod(v, {s}, p), p);
                if (degree(d) > 0) {
                    next.push_back(monic_mod(d, p));
                    rest = div_mod(rest, d, p);
                }
            }
            if (degree(rest) > 0) next.push_back(monic_mod(rest, p));
        }
        factors = std::move(next);
    }
    return factors;
}

// --- integer helpers ---------------------------------------------------------
[[nodiscard]] mpz_class content_of(const Z& a) {
    mpz_class g = 0;
    for (const auto& c : a) mpz_gcd(g.get_mpz_t(), g.get_mpz_t(), c.get_mpz_t());
    return g == 0 ? mpz_class(1) : g;
}
void make_primitive(Z& a) {
    mpz_class g = content_of(a);
    for (auto& c : a) c /= g;
    if (a.back() < 0) for (auto& c : a) c = -c;
}
[[nodiscard]] bool divides_z(Z a, const Z& b, Z& quo) {
    trim(a);
    Z bb = b;
    trim(bb);
    if (is_zero(bb)) return false;
    quo.assign(static_cast<std::size_t>(std::max<long>(0, degree(a) - degree(bb)) + 1), mpz_class(0));
    while (!is_zero(a) && degree(a) >= degree(bb)) {
        if (a.back() % bb.back() != 0) return false;
        mpz_class co = a.back() / bb.back();
        long sh = degree(a) - degree(bb);
        quo[static_cast<std::size_t>(sh)] = co;
        for (std::size_t i = 0; i < bb.size(); ++i) a[static_cast<std::size_t>(sh) + i] -= co * bb[i];
        trim(a);
    }
    return is_zero(a);
}
// Primitive integer polynomial from rational coefficients (clears denominators).
[[nodiscard]] Z rational_coeffs_to_primitive(const std::vector<Expr>& coeffs) {
    std::vector<mpq_class> q;
    for (const auto& c : coeffs) {
        if (c && c->type_id() == TypeId::Integer)
            q.emplace_back(static_cast<const Integer&>(*c).value());
        else if (c && c->type_id() == TypeId::Rational)
            q.push_back(static_cast<const Rational&>(*c).value());
        else
            throw std::invalid_argument("factor_zassenhaus: rational polynomial required");
    }
    mpz_class lcm = 1;
    for (const auto& v : q) mpz_lcm(lcm.get_mpz_t(), lcm.get_mpz_t(), v.get_den().get_mpz_t());
    Z out;
    for (const auto& v : q) out.push_back(v.get_num() * (lcm / v.get_den()));
    if (out.empty()) out.push_back(mpz_class(0));
    trim(out);
    make_primitive(out);
    return out;
}
[[nodiscard]] Expr to_expr(const Z& a, const Expr& var) {
    Expr e = S::Zero();
    for (long i = 0; i <= degree(a); ++i) {
        if (a[static_cast<std::size_t>(i)] == 0) continue;
        Expr term = make<Integer>(a[static_cast<std::size_t>(i)]);
        if (i > 0) term = mul(term, pow(var, integer(i)));
        e = add(e, term);
    }
    return e;
}

// Factor a squarefree primitive integer polynomial into irreducibles over ℤ:
// Berlekamp modulo a *small* prime, Hensel-lift the factorization to a prime
// power above the Landau–Mignotte bound, then recombine over symmetric residues.
[[nodiscard]] std::vector<Z> zz_factor_sqf(const Z& f) {
    long n = degree(f);
    if (n <= 1) return {f};
    mpz_class lc = f.back();

    mpz_class norm_sq = 0;
    for (const auto& c : f) norm_sq += c * c;
    mpz_class norm = 0;
    mpz_sqrt(norm.get_mpz_t(), norm_sq.get_mpz_t());
    norm += 1;
    mpz_class binom;
    mpz_bin_uiui(binom.get_mpz_t(), static_cast<unsigned long>(n - 1),
                 static_cast<unsigned long>((n - 1) / 2));
    mpz_class bound = binom * (norm + abs(lc));      // |coeff| of any factor
    mpz_class target = 2 * abs(lc) * bound + 1;       // symmetric-residue range

    // Smallest prime p with lc ≢ 0 and f squarefree mod p — keeps Berlekamp cheap.
    mpz_class p = 1;
    Z fp;
    for (;;) {
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
        if (mod(lc, p) == 0) continue;
        fp = monic_mod(f, p);
        if (degree(gcd_mod(fp, deriv(fp), p)) == 0) break;  // squarefree mod p
    }
    std::vector<Z> modfac = berlekamp(fp, p);
    if (modfac.size() == 1) return {f};

    // Lift the mod-p factorization to mod M = pˡ with pˡ ≥ target.
    long l = 1;
    mpz_class M = p;
    while (M < target) { M *= p; ++l; }
    if (l > 1) {
        mpz_class lc_inv = inv_mod(mod(lc, M), M);
        Z target_monic = trunc_sym(mul_z({lc_inv}, f), M);  // ≡ lc⁻¹·f, monic mod M
        modfac = hensel_lift_monic(target_monic, modfac, p, l);
    }

    auto symm = [&](Z a) {
        for (auto& c : a) { c = mod(c, M); if (2 * c > M) c -= M; }
        trim(a);
        return a;
    };
    std::vector<Z> result;
    std::vector<bool> used(modfac.size(), false);
    Z remaining = f;
    std::size_t r = 1;
    auto avail = [&] {
        return modfac.size() - static_cast<std::size_t>(std::count(used.begin(), used.end(), true));
    };
    while (r <= avail() / 2) {
        std::vector<std::size_t> idx;
        bool found = false;
        std::function<bool(std::size_t, std::size_t)> rec = [&](std::size_t start, std::size_t need) {
            if (need == 0) {
                Z cand = {remaining.back()};  // lc(remaining) · ∏ chosen
                for (std::size_t i : idx) cand = mul_mod(cand, modfac[i], M);
                cand = symm(cand);
                make_primitive(cand);
                Z quo;
                if (degree(cand) >= 1 && divides_z(remaining, cand, quo)) {
                    result.push_back(cand);
                    remaining = quo;
                    for (std::size_t i : idx) used[i] = true;
                    return true;
                }
                return false;
            }
            for (std::size_t i = start; i < modfac.size(); ++i) {
                if (used[i]) continue;
                idx.push_back(i);
                if (rec(i + 1, need - 1)) return true;
                idx.pop_back();
            }
            return false;
        };
        found = rec(0, r);
        if (found) { r = 1; continue; }
        ++r;
    }
    if (degree(remaining) >= 1) {
        make_primitive(remaining);
        result.push_back(remaining);
    }
    if (result.empty()) result.push_back(f);
    return result;
}

}  // namespace

std::vector<Expr> factor_zassenhaus(const Expr& p, const Expr& var) {
    Poly f(expand(p), var);
    if (f.degree() < 1) return {p};
    SqfList sqf = sqf_list(f);
    std::vector<Expr> out;
    for (const auto& [g_poly, mult] : sqf.factors) {
        Z g = rational_coeffs_to_primitive(g_poly.coeffs());
        if (degree(g) < 1) continue;
        for (const Z& irr : zz_factor_sqf(g)) {
            for (std::size_t t = 0; t < mult; ++t) out.push_back(to_expr(irr, var));
        }
    }
    if (out.empty()) out.push_back(p);
    return out;
}

}  // namespace sympp
