#include <sympp/core/rational.hpp>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include <gmpxx.h>

#include <sympp/core/integer.hpp>

namespace sympp {

namespace {

constexpr std::size_t kRationalHashSeed = 0xB4B4'4B4B'D2D2'2D2DULL;

[[nodiscard]] std::size_t hash_mpz_into(std::size_t h, const mpz_t z) noexcept {
    int s = mpz_sgn(z);
    h ^= static_cast<std::size_t>(s + 1);
    std::size_t n = mpz_size(z);
    for (std::size_t i = 0; i < n; ++i) {
        h ^= static_cast<std::size_t>(mpz_getlimbn(z, static_cast<mp_size_t>(i)));
        h *= 0x100'0000'01B3ULL;
    }
    return h;
}

[[nodiscard]] std::size_t hash_mpq(const mpq_class& q) noexcept {
    std::size_t h = kRationalHashSeed;
    h = hash_mpz_into(h, mpq_numref(q.get_mpq_t()));
    h = hash_mpz_into(h, mpq_denref(q.get_mpq_t()));
    return h;
}

}  // namespace

Rational::Rational(mpq_class v) : value_(std::move(v)) {
    value_.canonicalize();
    compute_hash();
}

Rational::Rational(long num, long den) : value_(num, den) {
    if (den == 0) {
        throw std::domain_error("Rational: division by zero");
    }
    value_.canonicalize();
    compute_hash();
}

void Rational::compute_hash() noexcept {
    hash_ = hash_mpq(value_);
}

bool Rational::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Rational) return false;
    return value_ == static_cast<const Rational&>(other).value_;
}

std::string Rational::str() const {
    return value_.get_str(10);
}

std::optional<bool> Rational::ask(AssumptionKey k) const noexcept {
    int s = sgn(value_);
    bool is_int = (value_.get_den() == 1);
    switch (k) {
        case AssumptionKey::Real: return true;
        case AssumptionKey::Rational: return true;
        case AssumptionKey::Integer: return is_int;
        case AssumptionKey::Finite: return true;
        case AssumptionKey::Positive: return s > 0;
        case AssumptionKey::Negative: return s < 0;
        case AssumptionKey::Zero: return s == 0;
        case AssumptionKey::Nonzero: return s != 0;
        case AssumptionKey::Nonnegative: return s >= 0;
        case AssumptionKey::Nonpositive: return s <= 0;
    }
    return std::nullopt;
}

// ---- Factories ------------------------------------------------------------

Expr rational(long num, long den) {
    if (den == 0) {
        throw std::domain_error("Rational: division by zero");
    }
    mpq_class q(num, den);
    q.canonicalize();
    if (q.get_den() == 1) {
        return make<Integer>(q.get_num());
    }
    return make<Rational>(std::move(q));
}

Expr rational(const mpz_class& num, const mpz_class& den) {
    if (den == 0) {
        throw std::domain_error("Rational: division by zero");
    }
    mpq_class q;
    mpq_set_num(q.get_mpq_t(), num.get_mpz_t());
    mpq_set_den(q.get_mpq_t(), den.get_mpz_t());
    q.canonicalize();
    if (q.get_den() == 1) {
        return make<Integer>(q.get_num());
    }
    return make<Rational>(std::move(q));
}

Expr rational(std::string_view s) {
    auto slash = s.find('/');
    if (slash == std::string_view::npos) {
        // No denominator — treat as integer.
        return integer(s);
    }
    auto num_str = std::string{s.substr(0, slash)};
    auto den_str = std::string{s.substr(slash + 1)};

    mpz_class num;
    mpz_class den;
    if (num.set_str(num_str, 10) != 0) {
        throw std::invalid_argument("Rational: invalid numerator: " + num_str);
    }
    if (den.set_str(den_str, 10) != 0) {
        throw std::invalid_argument("Rational: invalid denominator: " + den_str);
    }
    return rational(num, den);
}

}  // namespace sympp
