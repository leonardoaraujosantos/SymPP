#include <sympp/core/float.hpp>

#include <array>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>

#include <mpfr.h>

#include <sympp/core/integer.hpp>
#include <sympp/polys/rootof.hpp>
#include <sympp/core/number_symbol.hpp>
#include <sympp/core/rational.hpp>

namespace sympp {

namespace {

constexpr std::size_t kFloatHashSeed = 0xDEAD'BEEF'CAFE'BABEULL;

[[nodiscard]] std::size_t hash_mpfr(mpfr_srcptr x) noexcept {
    if (mpfr_zero_p(x) != 0) {
        return kFloatHashSeed;
    }
    if (mpfr_nan_p(x) != 0) {
        return kFloatHashSeed ^ 1ULL;
    }
    if (mpfr_inf_p(x) != 0) {
        return kFloatHashSeed ^ (mpfr_signbit(x) != 0 ? 2ULL : 3ULL);
    }
    // Mix sign + exponent + mantissa bits
    std::size_t h = kFloatHashSeed;
    h ^= static_cast<std::size_t>(mpfr_signbit(x) != 0 ? 1 : 0);
    auto exp = static_cast<long>(mpfr_get_exp(x));
    h ^= static_cast<std::size_t>(exp);
    h *= 0x100'0000'01B3ULL;

    // Pull the mantissa via mpfr_get_str into a deterministic representation
    mpfr_exp_t e_unused = 0;
    char* s = mpfr_get_str(nullptr, &e_unused, 16, 0, x, MPFR_RNDN);
    if (s != nullptr) {
        for (const char* p = s; *p != '\0'; ++p) {
            h ^= static_cast<std::size_t>(static_cast<unsigned char>(*p));
            h *= 0x100'0000'01B3ULL;
        }
        mpfr_free_str(s);
    }
    return h;
}

// Format as SymPy/mpmath does: `dps` significant digits, padded with trailing
// zeros if necessary. Fixed notation when the resulting decimal magnitude
// keeps the layout reasonable (mpmath threshold), otherwise scientific.
//
// Uses mpfr_get_str to extract digits + exponent, then composes the output
// manually so we control trailing-zero padding (which %g would strip).
[[nodiscard]] std::string format_mpfr(mpfr_srcptr x, int dps) {
    if (mpfr_nan_p(x) != 0) return "nan";
    if (mpfr_inf_p(x) != 0) return mpfr_signbit(x) != 0 ? "-inf" : "inf";
    if (mpfr_zero_p(x) != 0) return "0";

    mpfr_exp_t mpfr_exp = 0;
    char* raw = mpfr_get_str(nullptr, &mpfr_exp, 10,
                             static_cast<std::size_t>(dps), x, MPFR_RNDN);
    if (raw == nullptr) return "0";

    std::string digits{raw};
    mpfr_free_str(raw);

    bool negative = false;
    if (!digits.empty() && digits[0] == '-') {
        negative = true;
        digits.erase(0, 1);
    }

    // mpfr_get_str returns a string where the implicit decimal point sits
    // *before* digits[0], scaled by mpfr_exp. So the value = 0.<digits> * 10^exp.
    long exp = static_cast<long>(mpfr_exp);

    // Decide fixed vs scientific. mpmath threshold: fixed if exp ∈ [-4, dps+5].
    bool use_scientific = (exp <= -4) || (exp > dps + 4);

    std::string body;
    if (use_scientific) {
        // 1.<rest>e<exp-1>  — first digit before the point.
        body.reserve(digits.size() + 8);
        body.push_back(digits[0]);
        body.push_back('.');
        if (digits.size() > 1) {
            body.append(digits, 1, std::string::npos);
        } else {
            body.push_back('0');
        }
        body.push_back('e');
        long sci_exp = exp - 1;
        if (sci_exp >= 0) body.push_back('+');
        body.append(std::to_string(sci_exp));
    } else if (exp <= 0) {
        // 0.000<digits>  — leading zeros before the digits.
        body.reserve(digits.size() + 4 + static_cast<std::size_t>(-exp));
        body.append("0.");
        for (long i = 0; i < -exp; ++i) body.push_back('0');
        body.append(digits);
    } else if (static_cast<std::size_t>(exp) >= digits.size()) {
        // <digits>000.0  — integer with trailing zeros, then ".0".
        body.reserve(static_cast<std::size_t>(exp) + 2);
        body.append(digits);
        for (std::size_t i = digits.size(); i < static_cast<std::size_t>(exp); ++i) {
            body.push_back('0');
        }
        body.append(".0");
    } else {
        // <int>.<frac>
        body.reserve(digits.size() + 1);
        body.append(digits, 0, static_cast<std::size_t>(exp));
        body.push_back('.');
        body.append(digits, static_cast<std::size_t>(exp), std::string::npos);
    }

    if (negative) body.insert(body.begin(), '-');
    return body;
}

}  // namespace

Float::Float(double v, int dps) : dps_(dps) {
    mpfr_init2(value_, dps_to_prec(dps));
    mpfr_set_d(value_, v, MPFR_RNDN);
    compute_hash();
}

Float::Float(long v, int dps) : dps_(dps) {
    mpfr_init2(value_, dps_to_prec(dps));
    mpfr_set_si(value_, v, MPFR_RNDN);
    compute_hash();
}

Float::Float(std::string_view decimal, int dps) : dps_(dps) {
    mpfr_init2(value_, dps_to_prec(dps));
    std::string s{decimal};
    if (mpfr_set_str(value_, s.c_str(), 10, MPFR_RNDN) != 0) {
        mpfr_clear(value_);
        throw std::invalid_argument("Float: invalid decimal literal: " + s);
    }
    compute_hash();
}

Float::Float(mpz_srcptr v, int dps) : dps_(dps) {
    mpfr_init2(value_, dps_to_prec(dps));
    mpfr_set_z(value_, v, MPFR_RNDN);
    compute_hash();
}

Float::Float(mpq_srcptr v, int dps) : dps_(dps) {
    mpfr_init2(value_, dps_to_prec(dps));
    mpfr_set_q(value_, v, MPFR_RNDN);
    compute_hash();
}

Float::Float(mpfr_srcptr v, int dps) : dps_(dps) {
    mpfr_init2(value_, dps_to_prec(dps));
    mpfr_set(value_, v, MPFR_RNDN);
    compute_hash();
}

Float::Float(const Float& src, int new_dps) : dps_(new_dps) {
    mpfr_init2(value_, dps_to_prec(new_dps));
    mpfr_set(value_, src.value_, MPFR_RNDN);
    compute_hash();
}

Float::Float(const Float& other) : dps_(other.dps_) {
    mpfr_init2(value_, mpfr_get_prec(other.value_));
    mpfr_set(value_, other.value_, MPFR_RNDN);
    hash_ = other.hash_;
}

Float::Float(Float&& other) noexcept : dps_(other.dps_) {
    // mpfr_t is an array — swap the trivial struct contents.
    std::memcpy(&value_, &other.value_, sizeof(mpfr_t));
    // Leave `other` in a valid, freshly-init'd state so its destructor's
    // mpfr_clear is a no-op-equivalent.
    mpfr_init2(other.value_, kDefaultPrecBits);
    mpfr_set_zero(other.value_, 1);
    hash_ = other.hash_;
}

Float::~Float() {
    mpfr_clear(value_);
}

void Float::compute_hash() noexcept {
    hash_ = hash_mpfr(value_);
}

bool Float::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Float) return false;
    const auto& f = static_cast<const Float&>(other);
    if (mpfr_get_prec(value_) != mpfr_get_prec(f.value_)) return false;
    return mpfr_equal_p(value_, f.value_) != 0;
}

std::string Float::str() const {
    return format_mpfr(value_, dps_);
}

std::optional<bool> Float::ask(AssumptionKey k) const noexcept {
    if (mpfr_nan_p(value_) != 0) return std::nullopt;
    bool finite = mpfr_inf_p(value_) == 0;
    bool zero = mpfr_zero_p(value_) != 0;
    bool negative_signbit = mpfr_signbit(value_) != 0;
    int s = zero ? 0 : (negative_signbit ? -1 : 1);
    switch (k) {
        case AssumptionKey::Real: return true;
        // mpmath/SymPy: Floats are *not* known to be rational/integer in the
        // assumption sense even when the value would fit. Match that.
        case AssumptionKey::Rational: return std::nullopt;
        case AssumptionKey::Integer: return std::nullopt;
        case AssumptionKey::Finite: return finite;
        case AssumptionKey::Positive: return s > 0;
        case AssumptionKey::Negative: return s < 0;
        case AssumptionKey::Zero: return s == 0;
        case AssumptionKey::Nonzero: return s != 0;
        case AssumptionKey::Nonnegative: return s >= 0;
        case AssumptionKey::Nonpositive: return s <= 0;
    }
    return std::nullopt;
}

// ---- evalf ---------------------------------------------------------------

Expr evalf(const Expr& e, int dps) {
    if (!e) return e;
    switch (e->type_id()) {
        case TypeId::Float: {
            const auto& f = static_cast<const Float&>(*e);
            return make<Float>(f, dps);
        }
        case TypeId::Integer: {
            const auto& z = static_cast<const Integer&>(*e);
            return make<Float>(z.value().get_mpz_t(), dps);
        }
        case TypeId::Rational: {
            const auto& q = static_cast<const Rational&>(*e);
            return make<Float>(q.value().get_mpq_t(), dps);
        }
        case TypeId::NumberSymbol: {
            const auto& sym = static_cast<const NumberSymbol&>(*e);
            mpfr_t tmp;
            mpfr_init2(tmp, dps_to_prec(dps));
            sym.mpfr_evalf(tmp, dps_to_prec(dps));
            auto out = make<Float>(static_cast<mpfr_srcptr>(tmp), dps);
            mpfr_clear(tmp);
            return out;
        }
        case TypeId::RootOf: {
            const auto& r = static_cast<const RootOf&>(*e);
            if (auto v = r.try_evalf(dps)) return *v;
            return e;
        }
        default:
            // Non-numeric — identity. Phase 1f will recurse into Add/Mul/Pow.
            return e;
    }
}

}  // namespace sympp
