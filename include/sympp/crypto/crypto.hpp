#pragma once

// Public-key cryptography primitives over arbitrary-precision integers (GMP).
//
//   auto k = crypto::rsa_key(integer(61), integer(53), integer(17));
//   auto c = crypto::rsa_encrypt(integer(65), k.n, k.e);   // → 2790
//   crypto::rsa_decrypt(c, k.n, k.d);                       // → 65
//
// Reference: sympy/crypto/crypto.py (rsa_*_key, encipher_rsa, decipher_rsa),
//            Diffie–Hellman and ElGamal.

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp::crypto {

// Modular exponentiation base^exp mod m, and modular inverse a⁻¹ mod m (throws
// when a is not invertible mod m). All arguments must be integers.
[[nodiscard]] SYMPP_EXPORT Expr mod_pow(const Expr& base, const Expr& exp, const Expr& m);
[[nodiscard]] SYMPP_EXPORT Expr mod_inverse(const Expr& a, const Expr& m);

// ----- RSA -------------------------------------------------------------------
struct RSAKey {
    Expr n;  // modulus p·q
    Expr e;  // public exponent
    Expr d;  // private exponent = e⁻¹ mod φ(n)
};
// Build an RSA key from primes p, q and public exponent e.
[[nodiscard]] SYMPP_EXPORT RSAKey rsa_key(const Expr& p, const Expr& q, const Expr& e);
[[nodiscard]] SYMPP_EXPORT Expr rsa_encrypt(const Expr& msg, const Expr& n, const Expr& e);
[[nodiscard]] SYMPP_EXPORT Expr rsa_decrypt(const Expr& cipher, const Expr& n, const Expr& d);

// ----- Diffie–Hellman --------------------------------------------------------
// Public value gˢ mod p, and the shared secret (other_public)ˢ mod p.
[[nodiscard]] SYMPP_EXPORT Expr dh_public(const Expr& g, const Expr& secret, const Expr& p);
[[nodiscard]] SYMPP_EXPORT Expr dh_shared_secret(const Expr& other_public,
                                                 const Expr& secret, const Expr& p);

// ----- ElGamal ---------------------------------------------------------------
struct ElGamalCipher {
    Expr c1, c2;
};
// Encrypt m with public key y = gˣ mod p using ephemeral key k.
[[nodiscard]] SYMPP_EXPORT ElGamalCipher elgamal_encrypt(const Expr& m, const Expr& p,
                                                         const Expr& g, const Expr& y,
                                                         const Expr& k);
[[nodiscard]] SYMPP_EXPORT Expr elgamal_decrypt(const ElGamalCipher& c, const Expr& x,
                                                const Expr& p);

// ----- Elliptic-curve cryptography -------------------------------------------
// Short Weierstrass curve y² = x³ + a·x + b over the prime field 𝔽ₚ. Points form
// an abelian group under the chord-and-tangent law with the point at infinity as
// identity. Reference: sympy.crypto / standard ECC (SEC 1).
struct ECCurve {
    Expr a, b, p;  // y² = x³ + a·x + b (mod p)
};
struct ECPoint {
    Expr x, y;            // affine coordinates (ignored when at_infinity)
    bool at_infinity = false;
};

// Group identity (point at infinity).
[[nodiscard]] SYMPP_EXPORT ECPoint ec_infinity();
// Whether P lies on the curve (the identity always does).
[[nodiscard]] SYMPP_EXPORT bool ec_is_on_curve(const ECCurve& c, const ECPoint& P);
[[nodiscard]] SYMPP_EXPORT ECPoint ec_neg(const ECCurve& c, const ECPoint& P);
[[nodiscard]] SYMPP_EXPORT ECPoint ec_add(const ECCurve& c, const ECPoint& P, const ECPoint& Q);
[[nodiscard]] SYMPP_EXPORT ECPoint ec_double(const ECCurve& c, const ECPoint& P);
// Scalar multiple k·P via double-and-add (k may be any integer; negatives flip).
[[nodiscard]] SYMPP_EXPORT ECPoint ec_scalar_mul(const ECCurve& c, const Expr& k,
                                                 const ECPoint& P);
[[nodiscard]] SYMPP_EXPORT bool ec_equal(const ECPoint& P, const ECPoint& Q);

// ECDH: a party's public point secret·G, and the shared secret secret·other_public.
[[nodiscard]] SYMPP_EXPORT ECPoint ecdh_public(const ECCurve& c, const Expr& secret,
                                               const ECPoint& G);
[[nodiscard]] SYMPP_EXPORT ECPoint ecdh_shared(const ECCurve& c, const Expr& secret,
                                               const ECPoint& other_public);

// ECDSA over a base point G of prime order n. `z` is the (already reduced) message
// hash, `k` the per-signature ephemeral nonce, `d` the private key.
struct ECDSASignature {
    Expr r, s;
};
[[nodiscard]] SYMPP_EXPORT ECDSASignature ecdsa_sign(const ECCurve& c, const ECPoint& G,
                                                     const Expr& n, const Expr& d,
                                                     const Expr& z, const Expr& k);
[[nodiscard]] SYMPP_EXPORT bool ecdsa_verify(const ECCurve& c, const ECPoint& G, const Expr& n,
                                             const ECPoint& Q, const Expr& z,
                                             const ECDSASignature& sig);

}  // namespace sympp::crypto
