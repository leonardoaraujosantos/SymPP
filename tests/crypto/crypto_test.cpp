// Public-key cryptography — round-trips plus RSA cross-checked against
// sympy.crypto.

#include <stdexcept>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/crypto/crypto.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
namespace cr = sympp::crypto;

namespace {
std::string crypto_oracle(nlohmann::json req) {
    req["op"] = "crypto";
    auto r = sympp::testing::Oracle::instance().send(req);
    REQUIRE(r.ok);
    return r.result_str();
}
}  // namespace

TEST_CASE("modular arithmetic primitives", "[crypto]") {
    REQUIRE(cr::mod_pow(integer(7), integer(13), integer(11)) == integer(2));  // 7^13 mod 11
    REQUIRE(cr::mod_inverse(integer(3), integer(11)) == integer(4));           // 3*4=12≡1
    REQUIRE_THROWS_AS(cr::mod_inverse(integer(2), integer(4)), std::invalid_argument);
}

TEST_CASE("RSA matches sympy.crypto and round-trips", "[crypto][oracle]") {
    auto k = cr::rsa_key(integer(61), integer(53), integer(17));
    REQUIRE(k.n == integer(3233));
    REQUIRE(k.d->str() == crypto_oracle({{"fn", "rsa_d"}, {"p", "61"}, {"q", "53"}, {"e", "17"}}));

    auto c = cr::rsa_encrypt(integer(65), k.n, k.e);
    REQUIRE(c->str() == crypto_oracle({{"fn", "rsa_encrypt"}, {"m", "65"},
                                       {"p", "61"}, {"q", "53"}, {"e", "17"}}));
    REQUIRE(cr::rsa_decrypt(c, k.n, k.d) == integer(65));  // round-trip
}

TEST_CASE("Diffie-Hellman shared secret agrees", "[crypto]") {
    Expr g = integer(5), p = integer(23), a = integer(6), b = integer(15);
    Expr A = cr::dh_public(g, a, p);  // 8
    Expr B = cr::dh_public(g, b, p);  // 19
    Expr sa = cr::dh_shared_secret(B, a, p);
    Expr sb = cr::dh_shared_secret(A, b, p);
    REQUIRE(sa == sb);
    REQUIRE(sa == integer(2));
}

TEST_CASE("ElGamal round-trip", "[crypto]") {
    Expr p = integer(23), g = integer(5), x = integer(6);
    Expr y = cr::dh_public(g, x, p);              // public key
    auto c = cr::elgamal_encrypt(integer(10), p, g, y, integer(3));
    REQUIRE(cr::elgamal_decrypt(c, x, p) == integer(10));
}

// ----- Elliptic-curve cryptography -------------------------------------------
// Textbook curve y² = x³ + 2x + 2 over 𝔽₁₇ with base point G = (5,1) of order 19.
namespace {
cr::ECCurve test_curve() {
    return {integer(2), integer(2), integer(17)};
}
cr::ECPoint P(int x, int y) { return {integer(x), integer(y), false}; }
}  // namespace

TEST_CASE("ECC: points lie on the curve and group law holds", "[crypto][ecc]") {
    auto c = test_curve();
    auto G = P(5, 1);
    REQUIRE(cr::ec_is_on_curve(c, G));
    REQUIRE(cr::ec_is_on_curve(c, cr::ec_infinity()));
    REQUIRE_FALSE(cr::ec_is_on_curve(c, P(5, 2)));

    // Identity laws: P + O = P, P + (−P) = O.
    REQUIRE(cr::ec_equal(cr::ec_add(c, G, cr::ec_infinity()), G));
    REQUIRE(cr::ec_add(c, G, cr::ec_neg(c, G)).at_infinity);

    // 2G via doubling and via addition agree, and match the reference (6,3).
    REQUIRE(cr::ec_equal(cr::ec_double(c, G), P(6, 3)));
    REQUIRE(cr::ec_equal(cr::ec_add(c, G, G), P(6, 3)));
}

TEST_CASE("ECC: scalar multiples match reference and wrap at the order",
          "[crypto][ecc]") {
    auto c = test_curve();
    auto G = P(5, 1);
    // k·G for k = 1..7 (independently computed reference values).
    REQUIRE(cr::ec_equal(cr::ec_scalar_mul(c, integer(1), G), P(5, 1)));
    REQUIRE(cr::ec_equal(cr::ec_scalar_mul(c, integer(2), G), P(6, 3)));
    REQUIRE(cr::ec_equal(cr::ec_scalar_mul(c, integer(3), G), P(10, 6)));
    REQUIRE(cr::ec_equal(cr::ec_scalar_mul(c, integer(4), G), P(3, 1)));
    REQUIRE(cr::ec_equal(cr::ec_scalar_mul(c, integer(5), G), P(9, 16)));
    REQUIRE(cr::ec_equal(cr::ec_scalar_mul(c, integer(7), G), P(0, 6)));

    // Order 19: 19·G = O, and (n+k)·G = k·G.
    REQUIRE(cr::ec_scalar_mul(c, integer(19), G).at_infinity);
    REQUIRE(cr::ec_equal(cr::ec_scalar_mul(c, integer(22), G),
                         cr::ec_scalar_mul(c, integer(3), G)));
    // Negative scalar: (−3)·G = −(3·G).
    REQUIRE(cr::ec_equal(cr::ec_scalar_mul(c, integer(-3), G),
                         cr::ec_neg(c, P(10, 6))));
    // All multiples stay on the curve.
    for (int k = 1; k <= 18; ++k)
        REQUIRE(cr::ec_is_on_curve(c, cr::ec_scalar_mul(c, integer(k), G)));
}

TEST_CASE("ECDH: both parties derive the same shared point", "[crypto][ecc]") {
    auto c = test_curve();
    auto G = P(5, 1);
    Expr da = integer(7), db = integer(11);          // private keys
    auto QA = cr::ecdh_public(c, da, G);             // public keys
    auto QB = cr::ecdh_public(c, db, G);
    auto sA = cr::ecdh_shared(c, da, QB);
    auto sB = cr::ecdh_shared(c, db, QA);
    REQUIRE(cr::ec_equal(sA, sB));                    // shared secret agrees
    REQUIRE(cr::ec_is_on_curve(c, sA));
    REQUIRE(cr::ec_equal(sA, cr::ec_scalar_mul(c, integer(77), G)));  // (7·11)·G
}

TEST_CASE("ECDSA: valid signatures verify, tampered ones do not",
          "[crypto][ecc]") {
    auto c = test_curve();
    auto G = P(5, 1);
    Expr n = integer(19);          // order of G
    Expr d = integer(7);           // private key
    auto Q = cr::ec_scalar_mul(c, d, G);  // public key
    Expr z = integer(10);          // message hash (already reduced mod n)
    Expr k = integer(5);           // ephemeral nonce

    auto sig = cr::ecdsa_sign(c, G, n, d, z, k);
    REQUIRE(cr::ecdsa_verify(c, G, n, Q, z, sig));

    // Wrong message → reject.
    REQUIRE_FALSE(cr::ecdsa_verify(c, G, n, Q, integer(11), sig));
    // Tampered signature (in-range but wrong) → reject.
    REQUIRE_FALSE(cr::ecdsa_verify(c, G, n, Q, z, cr::ECDSASignature{sig.r, integer(1)}));
    REQUIRE_FALSE(cr::ecdsa_verify(c, G, n, Q, z, cr::ECDSASignature{integer(1), sig.s}));
    // Wrong public key → reject.
    auto Qbad = cr::ec_scalar_mul(c, integer(8), G);
    REQUIRE_FALSE(cr::ecdsa_verify(c, G, n, Qbad, z, sig));
}
