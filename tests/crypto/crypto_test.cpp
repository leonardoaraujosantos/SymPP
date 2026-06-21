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
