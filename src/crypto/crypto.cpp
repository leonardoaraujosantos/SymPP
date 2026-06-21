#include <sympp/crypto/crypto.hpp>

#include <stdexcept>

#include <gmpxx.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp::crypto {

namespace {
[[nodiscard]] const mpz_class& as_int(const Expr& e, const char* fn) {
    if (!e || e->type_id() != TypeId::Integer) {
        throw std::invalid_argument(std::string(fn) + ": integer argument required");
    }
    return static_cast<const Integer&>(*e).value();
}
[[nodiscard]] mpz_class powm(const mpz_class& b, const mpz_class& e, const mpz_class& m) {
    mpz_class r;
    mpz_powm(r.get_mpz_t(), b.get_mpz_t(), e.get_mpz_t(), m.get_mpz_t());
    return r;
}
}  // namespace

Expr mod_pow(const Expr& base, const Expr& exp, const Expr& m) {
    return make<Integer>(powm(as_int(base, "mod_pow"), as_int(exp, "mod_pow"),
                              as_int(m, "mod_pow")));
}

Expr mod_inverse(const Expr& a, const Expr& m) {
    mpz_class r;
    if (mpz_invert(r.get_mpz_t(), as_int(a, "mod_inverse").get_mpz_t(),
                   as_int(m, "mod_inverse").get_mpz_t()) == 0) {
        throw std::invalid_argument("mod_inverse: argument not invertible modulo m");
    }
    return make<Integer>(std::move(r));
}

RSAKey rsa_key(const Expr& p, const Expr& q, const Expr& e) {
    const mpz_class& zp = as_int(p, "rsa_key");
    const mpz_class& zq = as_int(q, "rsa_key");
    const mpz_class& ze = as_int(e, "rsa_key");
    mpz_class n = zp * zq;
    mpz_class phi = (zp - 1) * (zq - 1);
    mpz_class d;
    if (mpz_invert(d.get_mpz_t(), ze.get_mpz_t(), phi.get_mpz_t()) == 0) {
        throw std::invalid_argument("rsa_key: e is not coprime with φ(n)");
    }
    return {make<Integer>(std::move(n)), e, make<Integer>(std::move(d))};
}

Expr rsa_encrypt(const Expr& msg, const Expr& n, const Expr& e) { return mod_pow(msg, e, n); }
Expr rsa_decrypt(const Expr& cipher, const Expr& n, const Expr& d) {
    return mod_pow(cipher, d, n);
}

Expr dh_public(const Expr& g, const Expr& secret, const Expr& p) {
    return mod_pow(g, secret, p);
}
Expr dh_shared_secret(const Expr& other_public, const Expr& secret, const Expr& p) {
    return mod_pow(other_public, secret, p);
}

ElGamalCipher elgamal_encrypt(const Expr& m, const Expr& p, const Expr& g, const Expr& y,
                              const Expr& k) {
    const mpz_class& zp = as_int(p, "elgamal_encrypt");
    mpz_class c1 = powm(as_int(g, "elgamal_encrypt"), as_int(k, "elgamal_encrypt"), zp);
    mpz_class s = powm(as_int(y, "elgamal_encrypt"), as_int(k, "elgamal_encrypt"), zp);
    mpz_class c2 = (as_int(m, "elgamal_encrypt") * s) % zp;
    return {make<Integer>(std::move(c1)), make<Integer>(std::move(c2))};
}

Expr elgamal_decrypt(const ElGamalCipher& c, const Expr& x, const Expr& p) {
    const mpz_class& zp = as_int(p, "elgamal_decrypt");
    mpz_class s = powm(as_int(c.c1, "elgamal_decrypt"), as_int(x, "elgamal_decrypt"), zp);
    mpz_class sinv;
    if (mpz_invert(sinv.get_mpz_t(), s.get_mpz_t(), zp.get_mpz_t()) == 0) {
        throw std::invalid_argument("elgamal_decrypt: shared secret not invertible");
    }
    mpz_class m = (as_int(c.c2, "elgamal_decrypt") * sinv) % zp;
    return make<Integer>(std::move(m));
}

}  // namespace sympp::crypto
