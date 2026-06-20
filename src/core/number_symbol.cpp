#include <sympp/core/number_symbol.hpp>

#include <cstddef>
#include <cstdint>

#include <mpfr.h>

namespace sympp {

bool NumberSymbol::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::NumberSymbol) return false;
    return kind() == static_cast<const NumberSymbol&>(other).kind();
}

std::optional<bool> NumberSymbol::ask(AssumptionKey k) const noexcept {
    // Pi, E, EulerGamma, Catalan — all positive, real, finite, nonzero. Pi and
    // E are proven irrational (hence ¬rational/¬integer/¬parity/¬prime), so the
    // generic layer derives Irrational=true for them. The rationality of
    // EulerGamma and Catalan is an open problem, so every rationality-dependent
    // predicate is Unknown for those two — matching SymPy (γ.is_rational is None).
    const bool proven_irrational = kind() == NumberSymbolKind::Pi
                                   || kind() == NumberSymbolKind::E;
    const auto rationality = proven_irrational ? std::optional<bool>{false}
                                               : std::nullopt;
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
        case AssumptionKey::Commutative:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real: return true;
        case AssumptionKey::Finite: return true;
        case AssumptionKey::Positive: return true;
        case AssumptionKey::Negative: return false;
        case AssumptionKey::Zero: return false;
        case AssumptionKey::Nonzero: return true;
        case AssumptionKey::Nonnegative: return true;
        case AssumptionKey::Nonpositive: return false;

        // Rationality-dependent — known only for the proven-irrational constants.
        case AssumptionKey::Rational: return rationality;
        case AssumptionKey::Integer: return rationality;
        case AssumptionKey::Even: return rationality;
        case AssumptionKey::Odd: return rationality;
        case AssumptionKey::Prime: return rationality;
        case AssumptionKey::Composite: return rationality;
        // Irrational follows from real ∧ ¬rational in the generic layer; leaving
        // it Unknown here lets that derivation answer (true for Pi/E, Unknown
        // for EulerGamma/Catalan).
        case AssumptionKey::Irrational: return std::nullopt;
        // Pi and e are proven transcendental (hence not algebraic). EulerGamma
        // and Catalan are not known to be either — Unknown.
        case AssumptionKey::Algebraic:
            return proven_irrational ? std::optional<bool>{false} : std::nullopt;
        case AssumptionKey::Transcendental:
            return proven_irrational ? std::optional<bool>{true} : std::nullopt;
        // All of π, e, γ, Catalan are finite real constants.
        case AssumptionKey::ExtendedReal: return true;
        case AssumptionKey::Infinite: return false;
    }
    return std::nullopt;
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
