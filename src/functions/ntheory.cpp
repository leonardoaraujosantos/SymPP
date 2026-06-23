#include <sympp/functions/ntheory.hpp>

#include <algorithm>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

[[nodiscard]] const mpz_class& as_int(const Expr& e, const char* fn) {
    if (!e || e->type_id() != TypeId::Integer) {
        throw std::invalid_argument(std::string(fn) + ": integer argument required");
    }
    return static_cast<const Integer&>(*e).value();
}

[[nodiscard]] bool is_prime_mpz(const mpz_class& n) {
    return mpz_probab_prime_p(n.get_mpz_t(), 25) != 0;
}

// Pollard's rho (Floyd cycle) — returns a nontrivial factor of an odd composite
// n. Retries with a different polynomial constant on failure.
[[nodiscard]] mpz_class pollard_rho(const mpz_class& n) {
    if (n % 2 == 0) return mpz_class(2);
    for (mpz_class c = 1; ; ++c) {
        mpz_class x = 2, y = 2, d = 1;
        auto step = [&](const mpz_class& v) { return (v * v + c) % n; };
        while (d == 1) {
            x = step(x);
            y = step(step(y));
            mpz_class diff = x - y;
            if (diff < 0) diff = -diff;
            mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());
        }
        if (d != n) return d;  // success; else cycle hit n — retry with new c
    }
}

// Factor |n| (n != 0, |n| != 1) into a multiset of primes.
[[nodiscard]] std::vector<mpz_class> prime_factors(mpz_class m) {
    if (m < 0) m = -m;
    std::vector<mpz_class> primes;
    // Strip small factors by trial division first (fast, and keeps rho on odd
    // composites only).
    for (mpz_class p = 2; p * p <= m && p < 1000000; ++p) {
        while (m % p == 0) {
            primes.push_back(p);
            m /= p;
        }
        if (p > 2) ++p;  // after 2, only test odd candidates
    }
    // Factor the remaining cofactor with rho.
    std::vector<mpz_class> stack;
    if (m > 1) stack.push_back(m);
    while (!stack.empty()) {
        mpz_class x = stack.back();
        stack.pop_back();
        if (x == 1) continue;
        if (is_prime_mpz(x)) {
            primes.push_back(x);
            continue;
        }
        mpz_class d = pollard_rho(x);
        stack.push_back(d);
        stack.push_back(x / d);
    }
    std::sort(primes.begin(), primes.end());
    return primes;
}

// Distinct prime factors of m (> 0).
[[nodiscard]] std::vector<mpz_class> distinct_primes(const mpz_class& m) {
    auto all = prime_factors(m);
    all.erase(std::unique(all.begin(), all.end()), all.end());
    return all;
}

// Euler's totient φ(m) for m > 0.
[[nodiscard]] mpz_class euler_phi(const mpz_class& m) {
    if (m == 1) return mpz_class(1);
    mpz_class phi = m;
    for (const mpz_class& p : distinct_primes(m)) phi = phi / p * (p - 1);
    return phi;
}

// Non-negative remainder x mod m (m > 0).
[[nodiscard]] mpz_class mod(const mpz_class& x, const mpz_class& m) {
    mpz_class r = x % m;
    if (r < 0) r += m;
    return r;
}

[[nodiscard]] mpz_class mod_pow(const mpz_class& base, const mpz_class& exp,
                                const mpz_class& mod) {
    mpz_class r;
    mpz_powm(r.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), mod.get_mpz_t());
    return r;
}

[[nodiscard]] mpz_class gcd_mpz(const mpz_class& a, const mpz_class& b) {
    mpz_class g;
    mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return g;
}

// Multiplicative order of a modulo n, given gcd(a, n) = 1.
[[nodiscard]] mpz_class order_mod(const mpz_class& a, const mpz_class& n) {
    mpz_class o = euler_phi(n);
    for (const mpz_class& q : distinct_primes(o)) {
        while (o % q == 0 && mod_pow(a, o / q, n) == 1) o /= q;
    }
    return o;
}

// Does (ℤ/nℤ)* admit a primitive root? (n = 2, 4, pᵏ, or 2·pᵏ, p odd prime.)
[[nodiscard]] bool has_primitive_root(const mpz_class& n) {
    if (n == 2 || n == 4) return true;
    mpz_class m = n;
    if (m % 2 == 0) m /= 2;
    if (m % 2 == 0) return false;             // divisible by 4 ⇒ no
    auto ps = distinct_primes(m);
    return m > 1 && ps.size() == 1 && ps[0] != 2;  // m is an odd prime power
}

}  // namespace

std::vector<std::pair<Expr, Expr>> factorint(const Expr& n) {
    const mpz_class& z = as_int(n, "factorint");
    if (z == 0) throw std::invalid_argument("factorint: n must be nonzero");
    std::vector<std::pair<Expr, Expr>> out;
    // Match SymPy: a negative n carries a (-1)^1 sign factor.
    if (z < 0) out.emplace_back(make<Integer>(mpz_class(-1)), make<Integer>(1));
    if (z == 1 || z == -1) return out;
    auto primes = prime_factors(z);
    for (std::size_t i = 0; i < primes.size();) {
        std::size_t j = i;
        while (j < primes.size() && primes[j] == primes[i]) ++j;
        out.emplace_back(make<Integer>(primes[i]),
                         make<Integer>(mpz_class(static_cast<long>(j - i))));
        i = j;
    }
    return out;
}

std::vector<Expr> divisors(const Expr& n) {
    const mpz_class& z = as_int(n, "divisors");
    std::vector<Expr> out;
    mpz_class m = z < 0 ? mpz_class(-z) : z;
    if (m == 0) return out;
    std::vector<mpz_class> divs{1};
    if (m > 1) {
        auto primes = prime_factors(m);
        std::size_t i = 0;
        while (i < primes.size()) {
            std::size_t j = i;
            while (j < primes.size() && primes[j] == primes[i]) ++j;
            const mpz_class& p = primes[i];
            std::size_t e = j - i;
            std::vector<mpz_class> next;
            for (const mpz_class& d : divs) {
                mpz_class pk = 1;
                for (std::size_t k = 0; k <= e; ++k) {
                    next.push_back(d * pk);
                    pk *= p;
                }
            }
            divs.swap(next);
            i = j;
        }
    }
    std::sort(divs.begin(), divs.end());
    for (const mpz_class& d : divs) out.push_back(make<Integer>(d));
    return out;
}

std::tuple<Expr, Expr, Expr> igcdex(const Expr& a, const Expr& b) {
    const mpz_class& za = as_int(a, "igcdex");
    const mpz_class& zb = as_int(b, "igcdex");
    mpz_class g, s, t;
    mpz_gcdext(g.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(),
               za.get_mpz_t(), zb.get_mpz_t());
    return {make<Integer>(std::move(s)), make<Integer>(std::move(t)),
            make<Integer>(std::move(g))};
}

Expr jacobi_symbol(const Expr& a, const Expr& n) {
    const mpz_class& za = as_int(a, "jacobi_symbol");
    const mpz_class& zn = as_int(n, "jacobi_symbol");
    if (zn <= 0 || zn % 2 == 0) {
        throw std::invalid_argument("jacobi_symbol: n must be a positive odd integer");
    }
    return make<Integer>(mpz_class(mpz_jacobi(za.get_mpz_t(), zn.get_mpz_t())));
}

std::vector<Expr> continued_fraction(const Expr& x) {
    mpz_class p, q;
    if (x && x->type_id() == TypeId::Integer) {
        p = static_cast<const Integer&>(*x).value();
        q = 1;
    } else if (x && x->type_id() == TypeId::Rational) {
        const mpq_class& r = static_cast<const Rational&>(*x).value();
        p = r.get_num();
        q = r.get_den();  // canonical: q > 0
    } else {
        throw std::invalid_argument("continued_fraction: rational argument required");
    }
    std::vector<Expr> out;
    while (q != 0) {
        mpz_class a, rem;
        mpz_fdiv_qr(a.get_mpz_t(), rem.get_mpz_t(), p.get_mpz_t(), q.get_mpz_t());
        out.push_back(make<Integer>(a));
        p = q;
        q = rem;
    }
    return out;
}

Expr n_order(const Expr& a, const Expr& n) {
    const mpz_class& za = as_int(a, "n_order");
    const mpz_class& zn = as_int(n, "n_order");
    if (zn <= 0) throw std::invalid_argument("n_order: n must be positive");
    mpz_class am = ((za % zn) + zn) % zn;
    if (gcd_mpz(am, zn) != 1) {
        throw std::invalid_argument("n_order: a and n must be coprime");
    }
    return make<Integer>(order_mod(am, zn));
}

std::optional<Expr> primitive_root(const Expr& n) {
    const mpz_class& zn = as_int(n, "primitive_root");
    if (zn < 2) throw std::invalid_argument("primitive_root: n must be > 1");
    if (!has_primitive_root(zn)) return std::nullopt;
    mpz_class phi = euler_phi(zn);
    for (mpz_class g = 1; g < zn; ++g) {
        if (gcd_mpz(g, zn) != 1) continue;
        if (order_mod(g, zn) == phi) return make<Integer>(g);
    }
    return std::nullopt;  // unreachable when has_primitive_root() is true
}

std::optional<Expr> sqrt_mod(const Expr& a, const Expr& p) {
    const mpz_class& za = as_int(a, "sqrt_mod");
    const mpz_class& zp = as_int(p, "sqrt_mod");
    if (zp < 2 || !is_prime_mpz(zp)) {
        throw std::invalid_argument("sqrt_mod: p must be prime");
    }
    mpz_class am = ((za % zp) + zp) % zp;
    if (am == 0) return make<Integer>(0);
    if (zp == 2) return make<Integer>(am);  // am == 1
    if (mod_pow(am, (zp - 1) / 2, zp) != 1) return std::nullopt;  // non-residue
    if (zp % 4 == 3) return make<Integer>(mod_pow(am, (zp + 1) / 4, zp));
    // Tonelli–Shanks for p ≡ 1 (mod 4).
    mpz_class q = zp - 1, s = 0;
    while (q % 2 == 0) { q /= 2; ++s; }
    mpz_class z = 2;
    while (mod_pow(z, (zp - 1) / 2, zp) != zp - 1) ++z;  // a quadratic non-residue
    mpz_class m = s, c = mod_pow(z, q, zp), t = mod_pow(am, q, zp);
    mpz_class r = mod_pow(am, (q + 1) / 2, zp);
    while (t != 1) {
        mpz_class i = 0, t2 = t;
        while (t2 != 1) { t2 = (t2 * t2) % zp; ++i; }
        mpz_class b = c;
        for (mpz_class j = 0; j < m - i - 1; ++j) b = (b * b) % zp;
        m = i;
        c = (b * b) % zp;
        t = (t * c) % zp;
        r = (r * b) % zp;
    }
    return make<Integer>(r);
}

// ----- Chinese Remainder Theorem ---------------------------------------------

std::optional<std::pair<Expr, Expr>> crt(const std::vector<Expr>& moduli,
                                         const std::vector<Expr>& residues) {
    if (moduli.size() != residues.size()) {
        throw std::invalid_argument("crt: moduli and residues must have equal length");
    }
    if (moduli.empty()) return std::make_pair(integer(0), integer(1));

    // Fold the congruences pairwise: combine (x ≡ r mod M) with (x ≡ rᵢ mod mᵢ).
    mpz_class M = 1, r = 0;
    for (std::size_t i = 0; i < moduli.size(); ++i) {
        mpz_class mi = as_int(moduli[i], "crt");
        if (mi <= 0) throw std::invalid_argument("crt: moduli must be positive");
        mpz_class ri = ((as_int(residues[i], "crt") % mi) + mi) % mi;

        mpz_class g, p, q;  // g = gcd(M, mi) = p·M + q·mi
        mpz_gcdext(g.get_mpz_t(), p.get_mpz_t(), q.get_mpz_t(), M.get_mpz_t(), mi.get_mpz_t());
        mpz_class diff = ri - r;
        if (diff % g != 0) return std::nullopt;  // inconsistent

        mpz_class lcm = M / g * mi;
        // x = r + M·p·(diff/g)  (mod lcm).
        mpz_class x = r + M * p % lcm * (diff / g) % lcm;
        r = ((x % lcm) + lcm) % lcm;
        M = lcm;
    }
    return std::make_pair(make<Integer>(r), make<Integer>(M));
}

// ----- Discrete logarithm (baby-step / giant-step) ---------------------------

std::optional<Expr> discrete_log(const Expr& n, const Expr& a, const Expr& b) {
    mpz_class zn = as_int(n, "discrete_log");
    if (zn <= 1) throw std::invalid_argument("discrete_log: modulus must be > 1");
    mpz_class za = ((as_int(a, "discrete_log") % zn) + zn) % zn;
    mpz_class zb = ((as_int(b, "discrete_log") % zn) + zn) % zn;

    if (za == 1) return make<Integer>(mpz_class(0));  // b⁰ ≡ 1

    mpz_class m;  // m = ceil(√n)
    mpz_sqrt(m.get_mpz_t(), zn.get_mpz_t());
    m += 1;

    // Baby steps: first occurrence of bʲ for j in [0, m).
    std::map<mpz_class, mpz_class> table;
    mpz_class e = 1;
    for (mpz_class j = 0; j < m; ++j) {
        if (table.find(e) == table.end()) table.emplace(e, j);
        e = (e * zb) % zn;
    }

    // Giant steps: factor = b^(−m) mod n; scan a·factorⁱ for i in [0, m).
    mpz_class bm;
    mpz_powm(bm.get_mpz_t(), zb.get_mpz_t(), m.get_mpz_t(), zn.get_mpz_t());
    mpz_class factor;
    if (mpz_invert(factor.get_mpz_t(), bm.get_mpz_t(), zn.get_mpz_t()) == 0) {
        // b not invertible mod n: fall back to a direct bounded scan.
        mpz_class val = 1;
        for (mpz_class x = 0; x < zn; ++x) {
            if (val == za) return make<Integer>(x);
            val = (val * zb) % zn;
        }
        return std::nullopt;
    }
    mpz_class gamma = za;
    for (mpz_class i = 0; i < m; ++i) {
        auto it = table.find(gamma);
        if (it != table.end()) return make<Integer>(i * m + it->second);
        gamma = (gamma * factor) % zn;
    }
    return std::nullopt;
}

// ----- Linear Diophantine ----------------------------------------------------

std::optional<DiophantineLinear> diop_linear(const Expr& a, const Expr& b, const Expr& c) {
    mpz_class za = as_int(a, "diop_linear");
    mpz_class zb = as_int(b, "diop_linear");
    mpz_class zc = as_int(c, "diop_linear");

    if (za == 0 && zb == 0) {
        if (zc != 0) return std::nullopt;  // 0 = c ≠ 0
        return DiophantineLinear{integer(0), integer(0), integer(0), integer(0)};
    }

    // g = gcd(a, b) = x·a + y·b.
    mpz_class g, x, y;
    mpz_gcdext(g.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), za.get_mpz_t(), zb.get_mpz_t());
    if (zc % g != 0) return std::nullopt;  // no integer solution

    mpz_class k = zc / g;
    mpz_class x0 = x * k, y0 = y * k;        // particular solution
    mpz_class dx = zb / g, dy = -(za / g);   // homogeneous step
    return DiophantineLinear{make<Integer>(x0), make<Integer>(y0), make<Integer>(dx),
                             make<Integer>(dy)};
}

// ----- Quadratic residues, Legendre symbol, Carmichael λ ---------------------

Expr legendre_symbol(const Expr& a, const Expr& p) {
    mpz_class zp = as_int(p, "legendre_symbol");
    if (zp <= 2 || zp % 2 == 0 || !is_prime_mpz(zp)) {
        throw std::invalid_argument("legendre_symbol: p must be an odd prime");
    }
    mpz_class za = mod(as_int(a, "legendre_symbol"), zp);
    if (za == 0) return integer(0);
    mpz_class r = mod_pow(za, (zp - 1) / 2, zp);  // Euler's criterion
    return r == 1 ? integer(1) : integer(-1);
}

bool is_quadratic_residue(const Expr& a, const Expr& n) {
    mpz_class zn = as_int(n, "is_quadratic_residue");
    if (zn <= 0) throw std::invalid_argument("is_quadratic_residue: n must be positive");
    mpz_class za = mod(as_int(a, "is_quadratic_residue"), zn);
    for (mpz_class x = 0; x <= zn / 2; ++x) {
        if (mod(x * x, zn) == za) return true;
    }
    return false;
}

std::vector<Expr> quadratic_residues(const Expr& n) {
    mpz_class zn = as_int(n, "quadratic_residues");
    if (zn <= 0) throw std::invalid_argument("quadratic_residues: n must be positive");
    std::set<mpz_class> qr;
    for (mpz_class x = 0; x <= zn / 2; ++x) qr.insert(mod(x * x, zn));
    std::vector<Expr> out;
    for (const auto& v : qr) out.push_back(make<Integer>(v));  // std::set is ordered
    return out;
}

Expr reduced_totient(const Expr& n) {
    mpz_class zn = as_int(n, "reduced_totient");
    if (zn <= 0) throw std::invalid_argument("reduced_totient: n must be positive");
    if (zn == 1) return integer(1);
    auto lcm = [](mpz_class a, mpz_class b) {
        if (a == 0) return b;
        mpz_class g;
        mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
        return mpz_class(a / g * b);
    };
    mpz_class result = 1;
    // λ over each prime power in the factorization.
    for (const auto& [pe, ee] : factorint(n)) {
        mpz_class pp = as_int(pe, "reduced_totient");
        long e = as_int(ee, "reduced_totient").get_si();
        mpz_class lam;
        if (pp == 2) {
            lam = (e == 1) ? mpz_class(1)
                  : (e == 2)
                      ? mpz_class(2)
                      : mpz_class(mpz_class(1) << static_cast<mp_bitcnt_t>(e - 2));
        } else {
            mpz_class pe_pow;
            mpz_pow_ui(pe_pow.get_mpz_t(), pp.get_mpz_t(), static_cast<unsigned long>(e - 1));
            lam = pe_pow * (pp - 1);  // p^{e−1}(p−1)
        }
        result = lcm(result, lam);
    }
    return make<Integer>(result);
}

// ----- Pell equation x² − D·y² = 1 -------------------------------------------

std::optional<std::pair<Expr, Expr>> diop_pell(const Expr& D) {
    mpz_class d = as_int(D, "diop_pell");
    if (d <= 0) return std::nullopt;
    mpz_class a0;
    mpz_sqrt(a0.get_mpz_t(), d.get_mpz_t());
    if (a0 * a0 == d) return std::nullopt;  // perfect square → no nontrivial solution

    // Convergents h/k of the continued fraction of √D; the first with
    // h² − D·k² = 1 is the fundamental solution.
    mpz_class m = 0, q = 1, a = a0;
    mpz_class hp = a0, hpp = 1, kp = 1, kpp = 0;
    while (hp * hp - d * kp * kp != 1) {
        m = q * a - m;
        q = (d - m * m) / q;
        a = (a0 + m) / q;
        mpz_class h = a * hp + hpp;
        hpp = hp;
        hp = h;
        mpz_class k = a * kp + kpp;
        kpp = kp;
        kp = k;
    }
    return std::make_pair(make<Integer>(hp), make<Integer>(kp));
}

// ----- Sums of squares -------------------------------------------------------

std::optional<std::pair<Expr, Expr>> sum_of_two_squares(const Expr& n) {
    mpz_class zn = as_int(n, "sum_of_two_squares");
    if (zn < 0) return std::nullopt;
    // a ≤ b ⇒ 2a² ≤ n; test whether n − a² is a perfect square.
    for (mpz_class a = 0; 2 * a * a <= zn; ++a) {
        mpz_class r = zn - a * a;
        mpz_class b;
        mpz_sqrt(b.get_mpz_t(), r.get_mpz_t());
        if (b * b == r) return std::make_pair(make<Integer>(a), make<Integer>(b));
    }
    return std::nullopt;
}

std::optional<std::vector<Expr>> sum_of_three_squares(const Expr& n) {
    mpz_class zn = as_int(n, "sum_of_three_squares");
    if (zn < 0) return std::nullopt;
    // Legendre: impossible iff n = 4ᵃ·(8b+7).
    {
        mpz_class t = zn;
        while (t > 0 && t % 4 == 0) t /= 4;
        if (t % 8 == 7) return std::nullopt;
    }
    // Brute search with 0 ≤ a ≤ b ≤ c (guaranteed to succeed by the test above).
    for (mpz_class a = 0; 3 * a * a <= zn; ++a) {
        for (mpz_class b = a; a * a + 2 * b * b <= zn; ++b) {
            mpz_class r = zn - a * a - b * b;
            mpz_class c;
            mpz_sqrt(c.get_mpz_t(), r.get_mpz_t());
            if (c >= b && c * c == r) {
                return std::vector<Expr>{make<Integer>(a), make<Integer>(b),
                                         make<Integer>(c)};
            }
        }
    }
    return std::nullopt;  // unreachable when the Legendre test passed
}

std::optional<std::vector<Expr>> sum_of_four_squares(const Expr& n) {
    mpz_class zn = as_int(n, "sum_of_four_squares");
    if (zn < 0) return std::nullopt;
    // Greedy descent: pick a, b (a ≥ b) then a two-square rep of the remainder.
    mpz_class amax;
    mpz_sqrt(amax.get_mpz_t(), zn.get_mpz_t());
    for (mpz_class a = amax; a >= 0; --a) {
        mpz_class r1 = zn - a * a;
        mpz_class bmax;
        mpz_sqrt(bmax.get_mpz_t(), r1.get_mpz_t());
        if (bmax > a) bmax = a;  // keep b ≤ a
        for (mpz_class b = bmax; b >= 0; --b) {
            mpz_class r2 = r1 - b * b;
            // Two-square rep of r2 with both ≤ b (so the quadruple is sorted).
            for (mpz_class c = b; 2 * c * c >= r2 && c >= 0; --c) {
                mpz_class r3 = r2 - c * c;
                if (r3 < 0) continue;
                mpz_class dd;
                mpz_sqrt(dd.get_mpz_t(), r3.get_mpz_t());
                if (dd * dd == r3 && dd <= c) {
                    return std::vector<Expr>{make<Integer>(a), make<Integer>(b),
                                             make<Integer>(c), make<Integer>(dd)};
                }
            }
        }
    }
    return std::nullopt;
}

}  // namespace sympp
