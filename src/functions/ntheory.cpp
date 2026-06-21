#include <sympp/functions/ntheory.hpp>

#include <algorithm>
#include <stdexcept>
#include <vector>

#include <gmpxx.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
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

}  // namespace sympp
