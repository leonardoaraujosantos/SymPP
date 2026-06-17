#include <sympp/core/infinity.hpp>

#include <gmpxx.h>

#include <sympp/core/assumption_key.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

// Distinct stable per-run hash constants so the atoms never collide with each
// other or with numeric hashes.
std::size_t Infinity::hash() const noexcept { return 0x00'1F'17'00'00'00'0001ULL; }
std::size_t NegativeInfinity::hash() const noexcept {
    return 0x00'1F'17'00'00'00'0002ULL;
}
std::size_t ComplexInfinity::hash() const noexcept {
    return 0x00'1F'17'00'00'00'0003ULL;
}
std::size_t NaN::hash() const noexcept { return 0x00'1F'17'00'00'00'0004ULL; }

std::optional<bool> Infinity::ask(AssumptionKey k) const noexcept {
    // +oo is positive and non-zero, but neither finite nor (ordinary) real.
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return false;  // not a finite complex number

        case AssumptionKey::Positive:
        case AssumptionKey::Nonnegative:
        case AssumptionKey::Nonzero:
            return true;
        case AssumptionKey::Finite:
        case AssumptionKey::Real:
        case AssumptionKey::Integer:
        case AssumptionKey::Rational:
        case AssumptionKey::Negative:
        case AssumptionKey::Nonpositive:
        case AssumptionKey::Zero:
        case AssumptionKey::Even:
        case AssumptionKey::Odd:
        case AssumptionKey::Prime:
        case AssumptionKey::Composite:
        case AssumptionKey::Irrational:  // infinities are not real numbers
        case AssumptionKey::Algebraic:   // nor finite complex numbers
        case AssumptionKey::Transcendental:
            return false;
    }
    return std::nullopt;
}

std::optional<bool> NegativeInfinity::ask(AssumptionKey k) const noexcept {
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return false;  // not a finite complex number

        case AssumptionKey::Negative:
        case AssumptionKey::Nonpositive:
        case AssumptionKey::Nonzero:
            return true;
        case AssumptionKey::Finite:
        case AssumptionKey::Real:
        case AssumptionKey::Integer:
        case AssumptionKey::Rational:
        case AssumptionKey::Positive:
        case AssumptionKey::Nonnegative:
        case AssumptionKey::Zero:
        case AssumptionKey::Even:
        case AssumptionKey::Odd:
        case AssumptionKey::Prime:
        case AssumptionKey::Composite:
        case AssumptionKey::Irrational:  // infinities are not real numbers
        case AssumptionKey::Algebraic:   // nor finite complex numbers
        case AssumptionKey::Transcendental:
            return false;
    }
    return std::nullopt;
}

std::optional<bool> ComplexInfinity::ask(AssumptionKey k) const noexcept {
    // zoo has no determined sign and is not real/finite; it is non-zero.
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return false;  // not a finite complex number

        case AssumptionKey::Nonzero:
            return true;
        case AssumptionKey::Finite:
        case AssumptionKey::Real:
        case AssumptionKey::Integer:
        case AssumptionKey::Rational:
        case AssumptionKey::Positive:
        case AssumptionKey::Negative:
        case AssumptionKey::Nonnegative:
        case AssumptionKey::Nonpositive:
        case AssumptionKey::Zero:
        case AssumptionKey::Even:
        case AssumptionKey::Odd:
        case AssumptionKey::Prime:
        case AssumptionKey::Composite:
        case AssumptionKey::Irrational:  // infinities are not real numbers
        case AssumptionKey::Algebraic:   // nor finite complex numbers
        case AssumptionKey::Transcendental:
            return false;
    }
    return std::nullopt;
}

std::optional<bool> NaN::ask(AssumptionKey /*k*/) const noexcept {
    // nan is indeterminate under every predicate.
    return std::nullopt;
}

bool is_infinity(const Expr& e) noexcept {
    if (!e) return false;
    auto t = e->type_id();
    return t == TypeId::Infinity || t == TypeId::NegativeInfinity
           || t == TypeId::ComplexInfinity;
}

std::optional<Expr> pow_with_infinity(const Expr& base, const Expr& exp) {
    if (base->type_id() == TypeId::NaN || exp->type_id() == TypeId::NaN) {
        return S::NaN();
    }
    if (!is_infinity(base) && !is_infinity(exp)) return std::nullopt;

    const auto exp_pos = is_positive(exp);
    const auto exp_neg = is_negative(exp);

    // base = +oo
    if (base->type_id() == TypeId::Infinity) {
        if (exp == S::Zero()) return S::One();
        if (exp_pos == true) return S::Infinity();  // oo^positive = oo
        if (exp_neg == true) return S::Zero();       // oo^negative = 0
        return std::nullopt;
    }
    // base = -oo
    if (base->type_id() == TypeId::NegativeInfinity) {
        if (exp == S::Zero()) return S::One();
        if (exp_neg == true) return S::Zero();
        if (exp->type_id() == TypeId::Infinity) return S::NaN();  // (-oo)^oo
        if (exp->type_id() == TypeId::Integer) {
            const auto& z = static_cast<const Integer&>(*exp);
            if (z.is_positive()) {
                const bool even = mpz_even_p(z.value().get_mpz_t()) != 0;
                return even ? S::Infinity() : S::NegativeInfinity();
            }
        }
        return std::nullopt;  // non-integer positive exponent → complex; leave
    }
    // base = zoo
    if (base->type_id() == TypeId::ComplexInfinity) {
        if (exp == S::Zero()) return S::One();
        if (exp_pos == true) return S::ComplexInfinity();
        if (exp_neg == true) return S::Zero();
        return std::nullopt;
    }

    // From here, base is finite and exp is an infinity.
    if (is_number(base) && static_cast<const Number&>(*base).is_zero()) {
        if (exp->type_id() == TypeId::Infinity) return S::Zero();          // 0^oo
        if (exp->type_id() == TypeId::NegativeInfinity) {
            return S::ComplexInfinity();                                   // 0^-oo
        }
        return std::nullopt;
    }
    if (is_number(base)) {
        const auto& bn = static_cast<const Number&>(*base);
        if (bn.sign() < 0) {
            if (base == S::NegativeOne()) return S::NaN();  // (-1)^(±oo)
            return std::nullopt;  // other negative bases → complex; leave
        }
        // Positive base: compare its magnitude to 1 via sign(base - 1).
        auto bm1 = number_add(bn, static_cast<const Number&>(*S::NegativeOne()));
        if (!bm1) return std::nullopt;
        const int cmp = static_cast<const Number&>(**bm1).sign();
        if (exp->type_id() == TypeId::Infinity) {
            if (cmp > 0) return S::Infinity();  // base>1 → oo
            if (cmp < 0) return S::Zero();       // base<1 → 0
            return S::NaN();                     // base==1 → nan
        }
        if (exp->type_id() == TypeId::NegativeInfinity) {
            if (cmp > 0) return S::Zero();
            if (cmp < 0) return S::Infinity();
            return S::NaN();
        }
        return std::nullopt;
    }
    return std::nullopt;  // symbolic finite base, infinite exponent → leave
}

}  // namespace sympp
