#include <sympp/core/number_symbol.hpp>

#include <cstddef>
#include <cstdint>

#include <mpfr.h>

namespace sympp {

bool NumberSymbol::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::NumberSymbol) return false;
    return kind() == static_cast<const NumberSymbol&>(other).kind();
}

namespace {

// Stable per-kind hashes — chosen at random, kept constant across versions.
constexpr std::size_t kHashPi = 0x3243'F6A8'8854'46A2ULL;          // π first digits
constexpr std::size_t kHashE = 0xB7E1'5162'8AED'2A6BULL;           // e first digits
constexpr std::size_t kHashEulerGamma = 0x9325'1A4F'A4D2'8A9CULL;  // γ first digits
constexpr std::size_t kHashCatalan = 0x9159'6553'5BBC'BB87ULL;     // K first digits

}  // namespace

PiSymbol::PiSymbol() = default;
std::size_t PiSymbol::hash() const noexcept { return kHashPi; }
void PiSymbol::mpfr_evalf(mpfr_t out, mpfr_prec_t prec) const {
    mpfr_set_prec(out, prec);
    mpfr_const_pi(out, MPFR_RNDN);
}

ESymbol::ESymbol() = default;
std::size_t ESymbol::hash() const noexcept { return kHashE; }
void ESymbol::mpfr_evalf(mpfr_t out, mpfr_prec_t prec) const {
    mpfr_set_prec(out, prec);
    // MPFR doesn't ship a const_e; compute as exp(1).
    mpfr_set_si(out, 1, MPFR_RNDN);
    mpfr_exp(out, out, MPFR_RNDN);
}

EulerGammaSymbol::EulerGammaSymbol() = default;
std::size_t EulerGammaSymbol::hash() const noexcept { return kHashEulerGamma; }
void EulerGammaSymbol::mpfr_evalf(mpfr_t out, mpfr_prec_t prec) const {
    mpfr_set_prec(out, prec);
    mpfr_const_euler(out, MPFR_RNDN);
}

CatalanSymbol::CatalanSymbol() = default;
std::size_t CatalanSymbol::hash() const noexcept { return kHashCatalan; }
void CatalanSymbol::mpfr_evalf(mpfr_t out, mpfr_prec_t prec) const {
    mpfr_set_prec(out, prec);
    mpfr_const_catalan(out, MPFR_RNDN);
}

}  // namespace sympp
