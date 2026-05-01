#include <sympp/core/integer.hpp>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace sympp {

namespace {

// Mix the type tag in so Integer(7) and a future Float(7.0) can never share
// a hash bucket inadvertently.
constexpr std::size_t kIntegerHashSeed = 0xA5A5'5A5A'C3C3'3C3CULL;

// FNV-1a-style mix over the limbs of an mpz_t. Stable across runs and
// platform-independent (limb width is normalized via mpz_get_ui).
[[nodiscard]] std::size_t hash_mpz(const mpz_class& v) noexcept {
    std::size_t h = kIntegerHashSeed;
    int s = mpz_sgn(v.get_mpz_t());
    h ^= static_cast<std::size_t>(s + 1);  // -1, 0, +1 → 0, 1, 2
    std::size_t n_limbs = mpz_size(v.get_mpz_t());
    for (std::size_t i = 0; i < n_limbs; ++i) {
        auto limb = mpz_getlimbn(v.get_mpz_t(), static_cast<mp_size_t>(i));
        h ^= static_cast<std::size_t>(limb);
        h *= 0x100'0000'01B3ULL;  // FNV prime (64-bit)
    }
    return h;
}

}  // namespace

Integer::Integer(int v) : value_(v) {
    compute_hash();
}

Integer::Integer(long v) : value_(v) {
    compute_hash();
}

Integer::Integer(long long v) {
    // mpz_class doesn't have a long-long ctor on all platforms; route via string
    // when the value doesn't fit in long.
    if (v >= LONG_MIN && v <= LONG_MAX) {
        value_ = static_cast<long>(v);
    } else {
        value_.set_str(std::to_string(v), 10);
    }
    compute_hash();
}

Integer::Integer(unsigned long v) : value_(v) {
    compute_hash();
}

Integer::Integer(unsigned long long v) {
    if (v <= ULONG_MAX) {
        value_ = static_cast<unsigned long>(v);
    } else {
        value_.set_str(std::to_string(v), 10);
    }
    compute_hash();
}

Integer::Integer(std::string_view decimal) {
    // Reject anything that looks like a float — SymPy raises on Integer("10.5").
    // Reference: sympy/core/tests/test_numbers.py::test_Integer_new
    //     raises(ValueError, lambda: Integer("10.5"))
    for (char c : decimal) {
        if (c == '.' || c == 'e' || c == 'E') {
            throw std::invalid_argument("Integer: non-integer literal: " + std::string{decimal});
        }
    }
    if (value_.set_str(std::string{decimal}, 10) != 0) {
        throw std::invalid_argument("Integer: invalid decimal literal: " + std::string{decimal});
    }
    compute_hash();
}

Integer::Integer(mpz_class v) : value_(std::move(v)) {
    compute_hash();
}

void Integer::compute_hash() noexcept {
    hash_ = hash_mpz(value_);
}

bool Integer::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Integer) return false;
    return value_ == static_cast<const Integer&>(other).value_;
}

std::string Integer::str() const {
    return value_.get_str(10);
}

}  // namespace sympp
