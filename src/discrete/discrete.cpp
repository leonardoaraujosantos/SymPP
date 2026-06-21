#include <sympp/discrete/discrete.hpp>

#include <stdexcept>
#include <vector>

#include <gmpxx.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp::discrete {

namespace {

[[nodiscard]] std::size_t next_pow2(std::size_t n) {
    std::size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}
[[nodiscard]] std::vector<Expr> pad_pow2(std::vector<Expr> v) {
    std::size_t n = next_pow2(v.size());
    while (v.size() < n) v.push_back(S::Zero());
    return v;
}

// Discrete Fourier transform via exact roots of unity built from cos/sin of
// rational multiples of π (which SymPP evaluates to exact radicals). `sign` is
// −1 for the forward transform, +1 for the inverse.
[[nodiscard]] std::vector<Expr> dft(std::vector<Expr> x, int sign) {
    x = pad_pow2(std::move(x));
    long n = static_cast<long>(x.size());
    std::vector<Expr> out(static_cast<std::size_t>(n), S::Zero());
    for (long k = 0; k < n; ++k) {
        Expr acc = S::Zero();
        for (long j = 0; j < n; ++j) {
            Expr ang = mul(rational(2 * j * k, n), S::Pi());          // 2π jk/n
            Expr w = add(cos(ang), mul(integer(sign), mul(S::I(), sin(ang))));
            acc = add(acc, mul(x[static_cast<std::size_t>(j)], w));
        }
        out[static_cast<std::size_t>(k)] = simplify(acc);
    }
    return out;
}

[[nodiscard]] const mpz_class& as_int(const Expr& e, const char* fn) {
    if (!e || e->type_id() != TypeId::Integer) {
        throw std::invalid_argument(std::string(fn) + ": integer entries required");
    }
    return static_cast<const Integer&>(*e).value();
}
[[nodiscard]] mpz_class powm(const mpz_class& b, const mpz_class& e, const mpz_class& m) {
    mpz_class r;
    mpz_powm(r.get_mpz_t(), b.get_mpz_t(), e.get_mpz_t(), m.get_mpz_t());
    return r;
}
// Smallest primitive root modulo prime p.
[[nodiscard]] mpz_class primitive_root_mod(const mpz_class& p) {
    mpz_class phi = p - 1, m = phi;
    std::vector<mpz_class> primes;
    for (mpz_class d = 2; d * d <= m; ++d) {
        if (m % d == 0) { primes.push_back(d); while (m % d == 0) m /= d; }
    }
    if (m > 1) primes.push_back(m);
    for (mpz_class g = 2; g < p; ++g) {
        bool ok = true;
        for (const auto& q : primes) {
            if (powm(g, phi / q, p) == 1) { ok = false; break; }
        }
        if (ok) return g;
    }
    throw std::invalid_argument("ntt: no primitive root (p not prime?)");
}

[[nodiscard]] std::vector<Expr> ntt_impl(std::vector<Expr> seq, const mpz_class& p,
                                         bool inverse) {
    seq = pad_pow2(std::move(seq));
    long n = static_cast<long>(seq.size());
    if ((p - 1) % n != 0) {
        throw std::invalid_argument("ntt: transform length must divide p - 1");
    }
    std::vector<mpz_class> a;
    for (const auto& e : seq) a.push_back(((as_int(e, "ntt") % p) + p) % p);

    mpz_class g = primitive_root_mod(p);
    mpz_class w = powm(g, (p - 1) / n, p);  // primitive n-th root of unity mod p
    if (inverse) mpz_invert(w.get_mpz_t(), w.get_mpz_t(), p.get_mpz_t());

    std::vector<Expr> out(static_cast<std::size_t>(n));
    mpz_class ninv;
    mpz_invert(ninv.get_mpz_t(), mpz_class(n).get_mpz_t(), p.get_mpz_t());
    for (long k = 0; k < n; ++k) {
        mpz_class acc = 0;
        for (long j = 0; j < n; ++j) {
            acc = (acc + a[static_cast<std::size_t>(j)] * powm(w, mpz_class(j * k), p)) % p;
        }
        if (inverse) acc = (acc * ninv) % p;
        out[static_cast<std::size_t>(k)] = make<Integer>(acc);
    }
    return out;
}

}  // namespace

std::vector<Expr> fft(const std::vector<Expr>& seq) { return dft(seq, +1); }

std::vector<Expr> ifft(const std::vector<Expr>& seq) {
    auto r = dft(seq, -1);
    Expr inv = pow(integer(static_cast<long>(r.size())), integer(-1));
    for (auto& e : r) e = simplify(mul(e, inv));
    return r;
}

std::vector<Expr> ntt(const std::vector<Expr>& seq, const Expr& prime) {
    return ntt_impl(seq, as_int(prime, "ntt"), /*inverse=*/false);
}
std::vector<Expr> intt(const std::vector<Expr>& seq, const Expr& prime) {
    return ntt_impl(seq, as_int(prime, "intt"), /*inverse=*/true);
}

std::vector<Expr> convolution(const std::vector<Expr>& a, const std::vector<Expr>& b) {
    if (a.empty() || b.empty()) return {};
    std::vector<Expr> c(a.size() + b.size() - 1, S::Zero());
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            c[i + j] = add(c[i + j], mul(a[i], b[j]));
        }
    }
    for (auto& e : c) e = simplify(e);
    return c;
}

std::vector<Expr> mobius_transform(const std::vector<Expr>& seq) {
    std::vector<Expr> a = pad_pow2(seq);
    std::size_t n = a.size();
    for (std::size_t bit = 1; bit < n; bit <<= 1) {
        for (std::size_t mask = 0; mask < n; ++mask) {
            if (mask & bit) a[mask] = add(a[mask], a[mask ^ bit]);
        }
    }
    for (auto& e : a) e = simplify(e);
    return a;
}

std::vector<Expr> inverse_mobius_transform(const std::vector<Expr>& seq) {
    std::vector<Expr> a = pad_pow2(seq);
    std::size_t n = a.size();
    for (std::size_t bit = 1; bit < n; bit <<= 1) {
        for (std::size_t mask = 0; mask < n; ++mask) {
            if (mask & bit) a[mask] = add(a[mask], mul(integer(-1), a[mask ^ bit]));
        }
    }
    for (auto& e : a) e = simplify(e);
    return a;
}

}  // namespace sympp::discrete
