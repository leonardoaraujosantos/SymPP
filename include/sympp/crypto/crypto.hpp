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

}  // namespace sympp::crypto
