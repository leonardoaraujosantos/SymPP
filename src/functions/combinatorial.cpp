#include <sympp/functions/combinatorial.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

[[nodiscard]] Expr float_unary_op(int (*op)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t),
                                  const Expr& arg) {
    const auto& f = static_cast<const Float&>(*arg);
    int dps = f.precision_dps();
    mpfr_t r;
    mpfr_init2(r, dps_to_prec(dps));
    op(r, f.value(), MPFR_RNDN);
    auto out = make<Float>(static_cast<mpfr_srcptr>(r), dps);
    mpfr_clear(r);
    return out;
}

}  // namespace

// ============================================================================
// Factorial
// ============================================================================

Factorial::Factorial(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Factorial);
}
Expr Factorial::rebuild(std::vector<Expr> new_args) const {
    return factorial(new_args[0]);
}
std::optional<bool> Factorial::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
        case AssumptionKey::Integer:
        case AssumptionKey::Rational:
            if (is_integer(a) == true && is_nonnegative(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Positive:
            if (is_integer(a) == true && is_nonnegative(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr factorial(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (z.is_negative()) {
            // factorial of a negative integer: SymPy returns ComplexInfinity;
            // we keep it unevaluated for now since ComplexInfinity isn't
            // wired into the singletons yet.
            return make<Factorial>(arg);
        }
        if (!z.fits_long()) {
            // Astronomically huge — keep symbolic.
            return make<Factorial>(arg);
        }
        long n = z.to_long();
        if (n > 100'000) {
            return make<Factorial>(arg);  // safety bound
        }
        mpz_class r;
        mpz_fac_ui(r.get_mpz_t(), static_cast<unsigned long>(n));
        return make<Integer>(std::move(r));
    }
    return make<Factorial>(arg);
}

// ============================================================================
// Binomial
// ============================================================================

Binomial::Binomial(Expr n, Expr k)
    : Function(std::vector<Expr>{std::move(n), std::move(k)}) {
    compute_hash(FunctionId::Binomial);
}
Expr Binomial::rebuild(std::vector<Expr> new_args) const {
    return binomial(new_args[0], new_args[1]);
}
std::optional<bool> Binomial::ask(AssumptionKey k) const noexcept {
    const auto& n = args_[0];
    const auto& kk = args_[1];
    switch (k) {
        case AssumptionKey::Integer:
            if (is_integer(n) == true && is_integer(kk) == true) return true;
            return std::nullopt;
        case AssumptionKey::Nonnegative:
            if (is_nonnegative(n) == true && is_nonnegative(kk) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr binomial(const Expr& n, const Expr& k) {
    if (n->type_id() == TypeId::Integer && k->type_id() == TypeId::Integer) {
        const auto& zn = static_cast<const Integer&>(*n);
        const auto& zk = static_cast<const Integer&>(*k);
        if (!zn.is_negative() && !zk.is_negative()
            && zn.fits_long() && zk.fits_long()) {
            long nv = zn.to_long();
            long kv = zk.to_long();
            if (kv > nv) return S::Zero();
            mpz_class r;
            mpz_bin_uiui(r.get_mpz_t(),
                         static_cast<unsigned long>(nv),
                         static_cast<unsigned long>(kv));
            return make<Integer>(std::move(r));
        }
    }
    // binomial(n, 0) = 1
    if (k == S::Zero()) return S::One();
    // binomial(n, n) = 1 when both same and known nonneg integer
    if (n == k && is_integer(n) == true && is_nonnegative(n) == true) {
        return S::One();
    }
    return make<Binomial>(n, k);
}

// ============================================================================
// Gamma
// ============================================================================

GammaFn::GammaFn(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Gamma);
}
Expr GammaFn::rebuild(std::vector<Expr> new_args) const {
    return gamma(new_args[0]);
}
std::optional<bool> GammaFn::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Positive:
            // Gamma(x) > 0 for x > 0.
            if (is_positive(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Real:
            // Real for real positive arg; complications at non-positive
            // integers — defer.
            if (is_positive(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr gamma(const Expr& arg) {
    // gamma(n+1) = n! for nonneg Integer n means gamma(positive integer) = (n-1)!
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (z.is_positive() && z.fits_long()) {
            long n = z.to_long();
            if (n <= 100'000) {
                // gamma(n) = (n-1)!
                mpz_class r;
                mpz_fac_ui(r.get_mpz_t(), static_cast<unsigned long>(n - 1));
                return make<Integer>(std::move(r));
            }
        }
        // gamma(0) and gamma(negative integer) = pole; keep symbolic.
    }
    if (arg->type_id() == TypeId::Float) {
        return float_unary_op(mpfr_gamma, arg);
    }
    return make<GammaFn>(arg);
}

// ============================================================================
// LogGamma
// ============================================================================

LogGamma::LogGamma(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::LogGamma);
}
Expr LogGamma::rebuild(std::vector<Expr> new_args) const {
    return loggamma(new_args[0]);
}
std::optional<bool> LogGamma::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
            if (is_positive(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr loggamma(const Expr& arg) {
    if (arg == S::One()) return S::Zero();      // log(0!) = 0
    if (arg == integer(2)) return S::Zero();    // log(1!) = 0
    if (arg->type_id() == TypeId::Float) {
        return float_unary_op(mpfr_lngamma, arg);
    }
    return make<LogGamma>(arg);
}

}  // namespace sympp
