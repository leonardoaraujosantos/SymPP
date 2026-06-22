#include <sympp/crypto/crypto.hpp>

#include <stdexcept>

#include <gmpxx.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/singletons.hpp>
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

// ----- Elliptic-curve cryptography -------------------------------------------

namespace {
// Reduce x into [0, m) for a positive modulus m.
[[nodiscard]] mpz_class mod_pos(const mpz_class& x, const mpz_class& m) {
    mpz_class r = x % m;
    if (r < 0) r += m;
    return r;
}
[[nodiscard]] mpz_class inv_mod_z(const mpz_class& a, const mpz_class& m, const char* fn) {
    mpz_class r;
    if (mpz_invert(r.get_mpz_t(), mod_pos(a, m).get_mpz_t(), m.get_mpz_t()) == 0) {
        throw std::invalid_argument(std::string(fn) + ": value not invertible");
    }
    return r;
}
}  // namespace

ECPoint ec_infinity() { return ECPoint{S::Zero(), S::Zero(), true}; }

bool ec_equal(const ECPoint& P, const ECPoint& Q) {
    if (P.at_infinity || Q.at_infinity) return P.at_infinity == Q.at_infinity;
    return P.x == Q.x && P.y == Q.y;
}

bool ec_is_on_curve(const ECCurve& c, const ECPoint& P) {
    if (P.at_infinity) return true;
    const mpz_class& p = as_int(c.p, "ec_is_on_curve");
    mpz_class x = as_int(P.x, "ec_is_on_curve");
    mpz_class y = as_int(P.y, "ec_is_on_curve");
    mpz_class lhs = mod_pos(y * y, p);
    mpz_class rhs = mod_pos(x * x % p * x + as_int(c.a, "ec_is_on_curve") * x +
                                as_int(c.b, "ec_is_on_curve"),
                            p);
    return lhs == rhs;
}

ECPoint ec_neg(const ECCurve& c, const ECPoint& P) {
    if (P.at_infinity) return P;
    const mpz_class& p = as_int(c.p, "ec_neg");
    return ECPoint{P.x, make<Integer>(mod_pos(-as_int(P.y, "ec_neg"), p)), false};
}

namespace {
// Combine two points given the chord/tangent slope λ (mod p).
[[nodiscard]] ECPoint ec_from_slope(const mpz_class& lambda, const mpz_class& x1,
                                    const mpz_class& y1, const mpz_class& x2,
                                    const mpz_class& p) {
    mpz_class x3 = mod_pos(lambda * lambda - x1 - x2, p);
    mpz_class y3 = mod_pos(lambda * (x1 - x3) - y1, p);
    return ECPoint{make<Integer>(std::move(x3)), make<Integer>(std::move(y3)), false};
}
}  // namespace

ECPoint ec_double(const ECCurve& c, const ECPoint& P) {
    if (P.at_infinity) return P;
    const mpz_class& p = as_int(c.p, "ec_double");
    mpz_class x = as_int(P.x, "ec_double");
    mpz_class y = as_int(P.y, "ec_double");
    if (mod_pos(y, p) == 0) return ec_infinity();  // vertical tangent
    // λ = (3x² + a) / (2y).
    mpz_class num = mod_pos(3 * x * x + as_int(c.a, "ec_double"), p);
    mpz_class lambda = mod_pos(num * inv_mod_z(2 * y, p, "ec_double"), p);
    return ec_from_slope(lambda, x, y, x, p);
}

ECPoint ec_add(const ECCurve& c, const ECPoint& P, const ECPoint& Q) {
    if (P.at_infinity) return Q;
    if (Q.at_infinity) return P;
    const mpz_class& p = as_int(c.p, "ec_add");
    mpz_class x1 = as_int(P.x, "ec_add"), y1 = as_int(P.y, "ec_add");
    mpz_class x2 = as_int(Q.x, "ec_add"), y2 = as_int(Q.y, "ec_add");
    if (mod_pos(x1 - x2, p) == 0) {
        // Same x: either doubling (equal points) or mutual inverses → identity.
        if (mod_pos(y1 + y2, p) == 0) return ec_infinity();
        return ec_double(c, P);
    }
    mpz_class lambda = mod_pos((y2 - y1) * inv_mod_z(x2 - x1, p, "ec_add"), p);
    return ec_from_slope(lambda, x1, y1, x2, p);
}

ECPoint ec_scalar_mul(const ECCurve& c, const Expr& k, const ECPoint& P) {
    mpz_class scalar = as_int(k, "ec_scalar_mul");
    ECPoint base = P;
    if (scalar < 0) {  // (−k)·P = k·(−P)
        scalar = -scalar;
        base = ec_neg(c, P);
    }
    ECPoint result = ec_infinity();
    while (scalar > 0) {
        if (mpz_odd_p(scalar.get_mpz_t())) result = ec_add(c, result, base);
        base = ec_double(c, base);
        scalar /= 2;
    }
    return result;
}

ECPoint ecdh_public(const ECCurve& c, const Expr& secret, const ECPoint& G) {
    return ec_scalar_mul(c, secret, G);
}
ECPoint ecdh_shared(const ECCurve& c, const Expr& secret, const ECPoint& other_public) {
    return ec_scalar_mul(c, secret, other_public);
}

ECDSASignature ecdsa_sign(const ECCurve& c, const ECPoint& G, const Expr& n, const Expr& d,
                          const Expr& z, const Expr& k) {
    const mpz_class& zn = as_int(n, "ecdsa_sign");
    ECPoint kg = ec_scalar_mul(c, k, G);
    if (kg.at_infinity) throw std::invalid_argument("ecdsa_sign: k·G is the identity");
    mpz_class r = mod_pos(as_int(kg.x, "ecdsa_sign"), zn);
    if (r == 0) throw std::invalid_argument("ecdsa_sign: r == 0, choose another k");
    mpz_class kinv = inv_mod_z(as_int(k, "ecdsa_sign"), zn, "ecdsa_sign");
    mpz_class s = mod_pos(kinv * (as_int(z, "ecdsa_sign") +
                                  r * as_int(d, "ecdsa_sign")),
                          zn);
    if (s == 0) throw std::invalid_argument("ecdsa_sign: s == 0, choose another k");
    return {make<Integer>(std::move(r)), make<Integer>(std::move(s))};
}

bool ecdsa_verify(const ECCurve& c, const ECPoint& G, const Expr& n, const ECPoint& Q,
                  const Expr& z, const ECDSASignature& sig) {
    const mpz_class& zn = as_int(n, "ecdsa_verify");
    mpz_class r = as_int(sig.r, "ecdsa_verify");
    mpz_class s = as_int(sig.s, "ecdsa_verify");
    if (r <= 0 || r >= zn || s <= 0 || s >= zn) return false;
    mpz_class w = inv_mod_z(s, zn, "ecdsa_verify");
    mpz_class u1 = mod_pos(as_int(z, "ecdsa_verify") * w, zn);
    mpz_class u2 = mod_pos(r * w, zn);
    ECPoint pt = ec_add(c, ec_scalar_mul(c, make<Integer>(u1), G),
                        ec_scalar_mul(c, make<Integer>(u2), Q));
    if (pt.at_infinity) return false;
    return mod_pos(as_int(pt.x, "ecdsa_verify"), zn) == r;
}

}  // namespace sympp::crypto
