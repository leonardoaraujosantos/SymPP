#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>  // log
#include <sympp/functions/special.hpp>  // zeta

#include <optional>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/polys/poly.hpp>

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

// Exact rational coefficient C such that gamma(z) = C·sqrt(pi) for a
// half-integer z (odd numerator over 2). Reduces to the base gamma(1/2) via
// gamma(z) = (z-1)·gamma(z-1) when z > 1/2 and its inverse gamma(z) =
// gamma(z+1)/z when z < 1/2.
[[nodiscard]] mpq_class half_integer_gamma_coeff(mpq_class z) {
    const mpq_class half(1, 2);
    mpq_class coeff(1);
    while (z > half) {
        mpq_class z_minus_1 = z - 1;
        coeff *= z_minus_1;
        z = z_minus_1;
    }
    while (z < half) {
        coeff /= z;
        z += 1;
    }
    return coeff;
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
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

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
    // (+∞)! = +∞.
    if (arg->type_id() == TypeId::Infinity) return S::Infinity();
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (z.is_negative()) {
            // (−n)! = Γ(−n+1) has a pole for every positive integer n, so a
            // negative-integer factorial is ComplexInfinity (matching SymPy).
            return S::ComplexInfinity();
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
// Fibonacci
// ============================================================================

Fibonacci::Fibonacci(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Fibonacci);
}
Expr Fibonacci::rebuild(std::vector<Expr> new_args) const {
    return fibonacci(new_args[0]);
}
std::optional<bool> Fibonacci::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Integer:
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
            if (is_integer(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Nonnegative:
            if (is_integer(a) == true && is_nonnegative(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr fibonacci(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        // Negafibonacci F(-n) = (-1)^(n+1)·F(n) deferred; keep symbolic.
        if (z.is_negative()) return make<Fibonacci>(arg);
        if (!z.fits_long()) return make<Fibonacci>(arg);
        long n = z.to_long();
        if (n > 1'000'000) return make<Fibonacci>(arg);  // safety bound
        mpz_class r;
        mpz_fib_ui(r.get_mpz_t(), static_cast<unsigned long>(n));
        return make<Integer>(std::move(r));
    }
    return make<Fibonacci>(arg);
}

Lucas::Lucas(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Lucas);
}
Expr Lucas::rebuild(std::vector<Expr> new_args) const {
    return lucas(new_args[0]);
}
std::optional<bool> Lucas::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Integer:
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
            if (is_integer(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Positive:
            // L(n) ≥ 1 for every non-negative integer n (L(0)=2, L(1)=1, …).
            if (is_integer(a) == true && is_nonnegative(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr lucas(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        // Negative-index Lucas L(-n) = (-1)^n·L(n) deferred; keep symbolic.
        if (z.is_negative()) return make<Lucas>(arg);
        if (!z.fits_long()) return make<Lucas>(arg);
        long n = z.to_long();
        if (n > 1'000'000) return make<Lucas>(arg);  // safety bound
        mpz_class r;
        mpz_lucnum_ui(r.get_mpz_t(), static_cast<unsigned long>(n));
        return make<Integer>(std::move(r));
    }
    return make<Lucas>(arg);
}

// ============================================================================
// Totient (Euler's φ)
// ============================================================================

Totient::Totient(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Totient);
}
Expr Totient::rebuild(std::vector<Expr> new_args) const {
    return totient(new_args[0]);
}
std::optional<bool> Totient::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Integer:
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
            if (is_integer(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Positive:
            if (is_integer(a) == true && is_positive(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr totient(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        // φ is defined for n ≥ 1; non-positive n stays symbolic.
        if (!z.is_positive()) return make<Totient>(arg);
        mpz_class n = z.value();
        // φ(n) = n · ∏_{p|n} (1 − 1/p), computed by trial-dividing out each
        // distinct prime: when p divides the remaining m, do result -= result/p.
        mpz_class result = n;
        mpz_class m = n;
        for (mpz_class p = 2; p * p <= m; ++p) {
            if (mpz_divisible_p(m.get_mpz_t(), p.get_mpz_t())) {
                while (mpz_divisible_p(m.get_mpz_t(), p.get_mpz_t())) m /= p;
                result -= result / p;
            }
        }
        if (m > 1) result -= result / m;  // a remaining prime factor > √n
        return make<Integer>(std::move(result));
    }
    return make<Totient>(arg);
}

// ============================================================================
// Prime / PrimePi
// ============================================================================

namespace {
// Shared integer-arg / assumption boilerplate for the prime functions.
[[nodiscard]] std::optional<bool> prime_fn_ask(const Expr& a,
                                               AssumptionKey k) noexcept {
    switch (k) {
        case AssumptionKey::Integer:
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
            if (is_integer(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Positive:
            if (is_integer(a) == true && is_positive(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}
}  // namespace

Prime::Prime(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Prime);
}
Expr Prime::rebuild(std::vector<Expr> new_args) const {
    return prime(new_args[0]);
}
std::optional<bool> Prime::ask(AssumptionKey k) const noexcept {
    return prime_fn_ask(args_[0], k);
}

Expr prime(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        // prime(n) is the n-th prime, defined for n ≥ 1.
        if (!z.is_positive() || !z.fits_long()) return make<Prime>(arg);
        long n = z.to_long();
        if (n > 1'000'000) return make<Prime>(arg);  // safety bound
        mpz_class p = 1;
        for (long i = 0; i < n; ++i) {
            mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
        }
        return make<Integer>(std::move(p));
    }
    return make<Prime>(arg);
}

PrimePi::PrimePi(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::PrimePi);
}
Expr PrimePi::rebuild(std::vector<Expr> new_args) const {
    return primepi(new_args[0]);
}
std::optional<bool> PrimePi::ask(AssumptionKey k) const noexcept {
    return prime_fn_ask(args_[0], k);
}

Expr primepi(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        // π(n) counts primes ≤ n; π(n) = 0 for n < 2.
        if (z.is_negative() || !z.fits_long()) {
            if (z.is_negative()) return integer(0);
            return make<PrimePi>(arg);
        }
        long n = z.to_long();
        if (n < 2) return integer(0);
        if (n > 100'000'000) return make<PrimePi>(arg);  // safety bound
        long count = 0;
        mpz_class p = 1;
        const mpz_class limit = z.value();
        for (;;) {
            mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
            if (p > limit) break;
            ++count;
        }
        return integer(count);
    }
    return make<PrimePi>(arg);
}

// ============================================================================
// Multiplicative arithmetic functions (mobius / divisor_count / divisor_sigma)
// ============================================================================

namespace {
// Trial-division factorization of n > 0 into (prime, exponent) pairs.
[[nodiscard]] std::vector<std::pair<mpz_class, long>> factorize(mpz_class n) {
    std::vector<std::pair<mpz_class, long>> factors;
    for (mpz_class p = 2; p * p <= n; ++p) {
        if (mpz_divisible_p(n.get_mpz_t(), p.get_mpz_t())) {
            long e = 0;
            while (mpz_divisible_p(n.get_mpz_t(), p.get_mpz_t())) {
                n /= p;
                ++e;
            }
            factors.emplace_back(p, e);
        }
    }
    if (n > 1) factors.emplace_back(n, 1L);  // remaining prime factor > √n
    return factors;
}
}  // namespace

Mobius::Mobius(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Mobius);
}
Expr Mobius::rebuild(std::vector<Expr> new_args) const {
    return mobius(new_args[0]);
}
std::optional<bool> Mobius::ask(AssumptionKey k) const noexcept {
    return prime_fn_ask(args_[0], k);
}

Expr mobius(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (!z.is_positive()) return make<Mobius>(arg);
        if (z.value() == 1) return integer(1);  // μ(1) = 1
        auto factors = factorize(z.value());
        for (const auto& [p, e] : factors) {
            (void)p;
            if (e > 1) return integer(0);  // a squared factor ⇒ μ = 0
        }
        return integer(factors.size() % 2 == 0 ? 1 : -1);
    }
    return make<Mobius>(arg);
}

DivisorCount::DivisorCount(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::DivisorCount);
}
Expr DivisorCount::rebuild(std::vector<Expr> new_args) const {
    return divisor_count(new_args[0]);
}
std::optional<bool> DivisorCount::ask(AssumptionKey k) const noexcept {
    return prime_fn_ask(args_[0], k);
}

Expr divisor_count(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (!z.is_positive()) return make<DivisorCount>(arg);
        // σ₀(n) = ∏ (eᵢ + 1).
        mpz_class result = 1;
        for (const auto& [p, e] : factorize(z.value())) {
            (void)p;
            result *= (e + 1);
        }
        return make<Integer>(std::move(result));
    }
    return make<DivisorCount>(arg);
}

DivisorSigma::DivisorSigma(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::DivisorSigma);
}
DivisorSigma::DivisorSigma(Expr n, Expr k)
    : Function(std::vector<Expr>{std::move(n), std::move(k)}) {
    compute_hash(FunctionId::DivisorSigma);
}
Expr DivisorSigma::rebuild(std::vector<Expr> new_args) const {
    if (new_args.size() == 2) return divisor_sigma(new_args[0], new_args[1]);
    return divisor_sigma(new_args[0]);
}
std::optional<bool> DivisorSigma::ask(AssumptionKey k) const noexcept {
    return prime_fn_ask(args_[0], k);
}

Expr divisor_sigma(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (!z.is_positive()) return make<DivisorSigma>(arg);
        // σ₁(n) = ∏ (p^(eᵢ+1) − 1)/(p − 1).
        mpz_class result = 1;
        for (const auto& [p, e] : factorize(z.value())) {
            mpz_class num;
            mpz_pow_ui(num.get_mpz_t(), p.get_mpz_t(),
                       static_cast<unsigned long>(e + 1));
            num -= 1;
            result *= num / (p - 1);
        }
        return make<Integer>(std::move(result));
    }
    return make<DivisorSigma>(arg);
}

// σ_k(n) = Σ_{d|n} d^k = ∏_p (p^(k(eᵢ+1)) − 1)/(p^k − 1) for k ≥ 1; σ_0(n) =
// ∏ (eᵢ+1) is the divisor count. Evaluated for a positive integer n and a
// non-negative integer order k; symbolic otherwise.
Expr divisor_sigma(const Expr& n, const Expr& k) {
    if (n->type_id() == TypeId::Integer && k->type_id() == TypeId::Integer) {
        const auto& zn = static_cast<const Integer&>(*n);
        const auto& zk = static_cast<const Integer&>(*k);
        if (zn.is_positive() && !zk.is_negative() && zk.fits_long()) {
            const unsigned long ek = static_cast<unsigned long>(zk.to_long());
            if (ek == 1) return divisor_sigma(n);  // σ₁ — reuse the 1-arg path
            mpz_class result = 1;
            for (const auto& [p, e] : factorize(zn.value())) {
                if (ek == 0) {
                    result *= (e + 1);  // divisor count
                    continue;
                }
                mpz_class pk;
                mpz_pow_ui(pk.get_mpz_t(), p.get_mpz_t(), ek);  // p^k
                mpz_class num;
                mpz_pow_ui(num.get_mpz_t(), p.get_mpz_t(),
                           ek * static_cast<unsigned long>(e + 1));  // p^(k(e+1))
                num -= 1;
                result *= num / (pk - 1);
            }
            return make<Integer>(std::move(result));
        }
    }
    return make<DivisorSigma>(n, k);
}

// ============================================================================
// Harmonic / Factorial2
// ============================================================================

Harmonic::Harmonic(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Harmonic);
}
Harmonic::Harmonic(Expr n, Expr m)
    : Function(std::vector<Expr>{std::move(n), std::move(m)}) {
    compute_hash(FunctionId::Harmonic);
}
Expr Harmonic::rebuild(std::vector<Expr> new_args) const {
    if (new_args.size() == 2) return harmonic(new_args[0], new_args[1]);
    return harmonic(new_args[0]);
}
std::optional<bool> Harmonic::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
            if (is_integer(a) == true && is_nonnegative(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr harmonic(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (z.is_negative() || !z.fits_long()) return make<Harmonic>(arg);
        long n = z.to_long();
        if (n > 1'000'000) return make<Harmonic>(arg);  // safety bound
        // Hₙ = Σ_{k=1}^n 1/k, accumulated exactly as a rational.
        mpq_class sum(0);
        for (long k = 1; k <= n; ++k) sum += mpq_class(1, k);
        sum.canonicalize();
        return rational(sum.get_num(), sum.get_den());
    }
    return make<Harmonic>(arg);
}

// Generalized harmonic number Hₙ⁽ᵐ⁾ = Σ_{k=1}^n k^(−m). For m = 1 it is the
// ordinary Hₙ; m ≤ 0 gives a power sum. Evaluated exactly as a rational when n
// is a non-negative integer and m a (bounded) integer; symbolic otherwise.
Expr harmonic(const Expr& n, const Expr& m) {
    if (n->type_id() == TypeId::Integer && m->type_id() == TypeId::Integer) {
        const auto& zn = static_cast<const Integer&>(*n);
        const auto& zm = static_cast<const Integer&>(*m);
        if (!zn.is_negative() && zn.fits_long() && zm.fits_long()) {
            long N = zn.to_long();
            long M = zm.to_long();
            if (N <= 100'000 && M >= -1000 && M <= 1000) {
                mpq_class sum(0);
                for (long k = 1; k <= N; ++k) {
                    mpz_class kp;
                    mpz_class base(k);
                    mpz_pow_ui(kp.get_mpz_t(), base.get_mpz_t(),
                               static_cast<unsigned long>(M >= 0 ? M : -M));
                    // 1/kᴹ for M ≥ 0, else k^(−M) = k^|M|.
                    sum += (M >= 0) ? mpq_class(mpz_class(1), kp) : mpq_class(kp);
                }
                sum.canonicalize();
                return rational(sum.get_num(), sum.get_den());
            }
        }
    }
    return make<Harmonic>(n, m);
}

Factorial2::Factorial2(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Factorial2);
}
Expr Factorial2::rebuild(std::vector<Expr> new_args) const {
    return factorial2(new_args[0]);
}
std::optional<bool> Factorial2::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Integer:
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
            if (is_integer(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Positive:
            if (is_integer(a) == true && is_nonnegative(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr factorial2(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (!z.fits_long()) return make<Factorial2>(arg);
        long n = z.to_long();
        if (n == 0 || n == -1) return integer(1);  // empty product conventions
        if (n < -1) return make<Factorial2>(arg);   // negatives stay symbolic
        if (n > 1'000'000) return make<Factorial2>(arg);  // safety bound
        // n!! = n(n−2)(n−4)… down to 2 (n even) or 1 (n odd).
        mpz_class result(1);
        for (long k = n; k >= 1; k -= 2) result *= k;
        return make<Integer>(std::move(result));
    }
    return make<Factorial2>(arg);
}

// ============================================================================
// Bernoulli / Euler numbers
// ============================================================================

namespace {
[[nodiscard]] mpz_class binom(unsigned long n, unsigned long k) {
    mpz_class c;
    mpz_bin_uiui(c.get_mpz_t(), n, k);
    return c;
}
}  // namespace

Bernoulli::Bernoulli(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Bernoulli);
}
Expr Bernoulli::rebuild(std::vector<Expr> new_args) const {
    return bernoulli(new_args[0]);
}
std::optional<bool> Bernoulli::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
            if (is_integer(a) == true && is_nonnegative(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr bernoulli(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (z.is_negative() || !z.fits_long()) return make<Bernoulli>(arg);
        long n = z.to_long();
        if (n > 5000) return make<Bernoulli>(arg);  // safety bound (O(n²))
        if (n & 1L) return (n == 1) ? rational(1, 2) : integer(0);  // odd: B₁=½
        // Recurrence (B₁ = −1/2 internally): B_m = −1/(m+1) Σ_{k<m} C(m+1,k) B_k.
        std::vector<mpq_class> b(static_cast<std::size_t>(n) + 1);
        b[0] = 1;
        for (long m = 1; m <= n; ++m) {
            mpq_class s(0);
            for (long k = 0; k < m; ++k) {
                s += mpq_class(binom(static_cast<unsigned long>(m + 1),
                                     static_cast<unsigned long>(k)))
                     * b[static_cast<std::size_t>(k)];
            }
            b[static_cast<std::size_t>(m)] = -s / mpq_class(m + 1);
        }
        mpq_class r = b[static_cast<std::size_t>(n)];
        r.canonicalize();
        return rational(r.get_num(), r.get_den());
    }
    return make<Bernoulli>(arg);
}

Euler::Euler(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Euler);
}
Expr Euler::rebuild(std::vector<Expr> new_args) const {
    return euler(new_args[0]);
}
std::optional<bool> Euler::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Integer:
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
            if (is_integer(a) == true && is_nonnegative(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr euler(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (z.is_negative() || !z.fits_long()) return make<Euler>(arg);
        long n = z.to_long();
        if (n & 1L) return integer(0);          // odd Euler numbers vanish
        if (n > 5000) return make<Euler>(arg);  // safety bound (O(n²))
        // E_{2m} = −Σ_{k=0}^{m-1} C(2m, 2k) E_{2k}, indexing by half (e[j]=E_{2j}).
        const long half = n / 2;
        std::vector<mpz_class> e(static_cast<std::size_t>(half) + 1);
        e[0] = 1;
        for (long m = 1; m <= half; ++m) {
            mpz_class s(0);
            for (long k = 0; k < m; ++k) {
                s += binom(static_cast<unsigned long>(2 * m),
                           static_cast<unsigned long>(2 * k))
                     * e[static_cast<std::size_t>(k)];
            }
            e[static_cast<std::size_t>(m)] = -s;
        }
        return make<Integer>(e[static_cast<std::size_t>(half)]);
    }
    return make<Euler>(arg);
}

// ============================================================================
// Catalan
// ============================================================================

Catalan::Catalan(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Catalan);
}
Expr Catalan::rebuild(std::vector<Expr> new_args) const {
    return catalan(new_args[0]);
}
std::optional<bool> Catalan::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Integer:
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
        case AssumptionKey::Positive:
            if (is_integer(a) == true && is_nonnegative(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr catalan(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (z.is_negative()) return make<Catalan>(arg);
        if (!z.fits_long()) return make<Catalan>(arg);
        long n = z.to_long();
        if (n > 100'000) return make<Catalan>(arg);  // safety bound
        // Catalan(n) = binomial(2n, n) / (n + 1).
        mpz_class bin;
        mpz_bin_uiui(bin.get_mpz_t(), static_cast<unsigned long>(2 * n),
                     static_cast<unsigned long>(n));
        mpz_class r = bin / (n + 1);
        return make<Integer>(std::move(r));
    }
    return make<Catalan>(arg);
}

// ============================================================================
// Gcd / Lcm
// ============================================================================

Gcd::Gcd(Expr a, Expr b)
    : Function(std::vector<Expr>{std::move(a), std::move(b)}) {
    compute_hash(FunctionId::Gcd);
}
Expr Gcd::rebuild(std::vector<Expr> new_args) const {
    return gcd(new_args[0], new_args[1]);
}
std::optional<bool> Gcd::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    const auto& b = args_[1];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Integer:
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
        case AssumptionKey::Nonnegative:
            if (is_integer(a) == true && is_integer(b) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

namespace {
// Integer content of a polynomial: gcd of its coefficients. Returns nullopt if
// any coefficient is not an integer (so the caller skips content scaling).
[[nodiscard]] std::optional<mpz_class> integer_content(const Poly& p) {
    mpz_class g = 0;
    for (const auto& c : p.coeffs()) {
        if (c->type_id() != TypeId::Integer) return std::nullopt;
        mpz_gcd(g.get_mpz_t(), g.get_mpz_t(),
                static_cast<const Integer&>(*c).value().get_mpz_t());
    }
    return g;
}

// Convert the monic gcd `g` to SymPy's convention: the primitive integer
// polynomial (integer coefficients, content 1, positive leading) scaled by
// `content` — e.g. monic x+3/2 with content 2 → 2·(2x+3) is wrong; with
// content 1 → 2x+3. Clears denominators, divides by the integer content, then
// multiplies by the requested content factor.
[[nodiscard]] Expr gcd_to_primitive(const Poly& g, const mpz_class& content,
                                    const Expr& var) {
    const auto& cs = g.coeffs();
    mpz_class lcm_den = 1;
    for (const auto& c : cs) {
        if (c->type_id() == TypeId::Rational) {
            mpz_lcm(lcm_den.get_mpz_t(), lcm_den.get_mpz_t(),
                    static_cast<const Rational&>(*c).denominator().get_mpz_t());
        }
    }
    std::vector<mpz_class> ints(cs.size());
    mpz_class g_int = 0;
    for (std::size_t i = 0; i < cs.size(); ++i) {
        mpz_class m;
        if (cs[i]->type_id() == TypeId::Integer) {
            m = static_cast<const Integer&>(*cs[i]).value() * lcm_den;
        } else {
            const auto& q = static_cast<const Rational&>(*cs[i]);
            m = q.numerator() * lcm_den / q.denominator();
        }
        ints[i] = m;
        mpz_gcd(g_int.get_mpz_t(), g_int.get_mpz_t(), m.get_mpz_t());
    }
    if (g_int == 0) g_int = 1;
    std::vector<Expr> terms;
    for (std::size_t i = 0; i < ints.size(); ++i) {
        mpz_class coeff = content * (ints[i] / g_int);
        if (coeff == 0) continue;
        terms.push_back(mul(make<Integer>(std::move(coeff)),
                            pow(var, integer(static_cast<long>(i)))));
    }
    if (terms.empty()) return S::Zero();
    return add(std::move(terms));
}
}  // namespace

Expr gcd(const Expr& a, const Expr& b) {
    if (a->type_id() == TypeId::Integer && b->type_id() == TypeId::Integer) {
        // mpz_gcd yields the non-negative gcd (handles signs and zero:
        // gcd(0,0)=0, gcd(±n,0)=|n|), matching SymPy.
        mpz_class r;
        mpz_gcd(r.get_mpz_t(), static_cast<const Integer&>(*a).value().get_mpz_t(),
                static_cast<const Integer&>(*b).value().get_mpz_t());
        return make<Integer>(std::move(r));
    }
    // Univariate polynomial GCD. SymPy's convention is
    //   gcd = gcd(int-content a, int-content b) · monic-gcd(primitive parts),
    // e.g. gcd(x²−1, x−1) = x−1 and gcd(2x²−2, 2x−2) = 2x−2.
    {
        ExprSet vars;
        for (const auto& s : free_symbols(a)) vars.insert(s);
        for (const auto& s : free_symbols(b)) vars.insert(s);
        if (vars.size() == 1) {
            const Expr& var = *vars.begin();
            try {
                Poly pa(expand(a), var);
                Poly pb(expand(b), var);
                Poly g = gcd(pa, pb);  // monic gcd over ℚ
                // Content factor: gcd of the integer contents (1 if either side
                // has non-integer coefficients).
                mpz_class content = 1;
                if (auto ca = integer_content(pa)) {
                    if (auto cb = integer_content(pb)) {
                        mpz_gcd(content.get_mpz_t(), ca->get_mpz_t(),
                                cb->get_mpz_t());
                    }
                }
                return gcd_to_primitive(g, content, var);
            } catch (const std::exception&) {
                // not polynomial in `var` — fall through to the symbolic node
            }
        }
    }
    return make<Gcd>(a, b);
}

Lcm::Lcm(Expr a, Expr b)
    : Function(std::vector<Expr>{std::move(a), std::move(b)}) {
    compute_hash(FunctionId::Lcm);
}
Expr Lcm::rebuild(std::vector<Expr> new_args) const {
    return lcm(new_args[0], new_args[1]);
}
std::optional<bool> Lcm::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    const auto& b = args_[1];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Integer:
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
        case AssumptionKey::Nonnegative:
            if (is_integer(a) == true && is_integer(b) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr lcm(const Expr& a, const Expr& b) {
    if (a->type_id() == TypeId::Integer && b->type_id() == TypeId::Integer) {
        // mpz_lcm is non-negative; lcm(0,n)=0, matching SymPy.
        mpz_class r;
        mpz_lcm(r.get_mpz_t(), static_cast<const Integer&>(*a).value().get_mpz_t(),
                static_cast<const Integer&>(*b).value().get_mpz_t());
        return make<Integer>(std::move(r));
    }
    // Univariate polynomial LCM via lcm(a, b) = a·b / gcd(a, b), reusing the
    // SymPy-convention polynomial gcd. The exact division restores the right
    // content, so lcm(x²−1, x−1) = x²−1 and lcm(2x−2, 3x−3) = 6x−6.
    {
        ExprSet vars;
        for (const auto& s : free_symbols(a)) vars.insert(s);
        for (const auto& s : free_symbols(b)) vars.insert(s);
        if (vars.size() == 1) {
            const Expr& var = *vars.begin();
            try {
                Poly pa(expand(a), var);
                Poly pb(expand(b), var);
                Poly pg(expand(gcd(a, b)), var);
                Poly prod = pa * pb;
                auto [q, r] = prod.divmod(pg);
                if (!(r.as_expr() == S::Zero())) {
                    return make<Lcm>(a, b);  // not an exact division (shouldn't happen)
                }
                return q.as_expr();
            } catch (const std::exception&) {
                // not polynomial in `var` — fall through to the symbolic node
            }
        }
    }
    return make<Lcm>(a, b);
}

// ============================================================================
// RisingFactorial / FallingFactorial / Subfactorial
// ============================================================================

namespace {

// Build the rising (sign=+1) or falling (sign=-1) factorial product
// x·(x±1)·…·(x±(k-1)) for a non-negative integer k. Returns nullopt when n is
// not a usable non-negative integer (kept symbolic) or k exceeds the bound.
[[nodiscard]] std::optional<Expr> factorial_product(const Expr& x, const Expr& n,
                                                    int sign) {
    if (n->type_id() != TypeId::Integer) return std::nullopt;
    const auto& zn = static_cast<const Integer&>(*n);
    if (zn.is_negative() || !zn.fits_long()) return std::nullopt;
    long k = zn.to_long();
    if (k > 100'000) return std::nullopt;  // safety bound
    if (k == 0) return S::One();
    Expr prod = x;
    for (long i = 1; i < k; ++i) {
        prod = mul(prod, add(x, integer(sign * i)));
    }
    return prod;
}

}  // namespace

RisingFactorial::RisingFactorial(Expr x, Expr n)
    : Function(std::vector<Expr>{std::move(x), std::move(n)}) {
    compute_hash(FunctionId::RisingFactorial);
}
Expr RisingFactorial::rebuild(std::vector<Expr> new_args) const {
    return rising_factorial(new_args[0], new_args[1]);
}
std::optional<bool> RisingFactorial::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Integer && is_integer(args_[0]) == true
        && is_integer(args_[1]) == true) {
        return true;
    }
    return std::nullopt;
}

Expr rising_factorial(const Expr& x, const Expr& n) {
    // rf(x,n) = x·(x+1)·…·(x+n-1) for a non-negative integer n (n=0 → 1).
    if (auto p = factorial_product(x, n, +1)) return *p;
    return make<RisingFactorial>(x, n);
}

FallingFactorial::FallingFactorial(Expr x, Expr n)
    : Function(std::vector<Expr>{std::move(x), std::move(n)}) {
    compute_hash(FunctionId::FallingFactorial);
}
Expr FallingFactorial::rebuild(std::vector<Expr> new_args) const {
    return falling_factorial(new_args[0], new_args[1]);
}
std::optional<bool> FallingFactorial::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Integer && is_integer(args_[0]) == true
        && is_integer(args_[1]) == true) {
        return true;
    }
    return std::nullopt;
}

Expr falling_factorial(const Expr& x, const Expr& n) {
    // ff(x,n) = x·(x-1)·…·(x-n+1) for a non-negative integer n (n=0 → 1).
    if (auto p = factorial_product(x, n, -1)) return *p;
    return make<FallingFactorial>(x, n);
}

Subfactorial::Subfactorial(Expr arg)
    : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Subfactorial);
}
Expr Subfactorial::rebuild(std::vector<Expr> new_args) const {
    return subfactorial(new_args[0]);
}
std::optional<bool> Subfactorial::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Integer:
        case AssumptionKey::Nonnegative:
            if (is_integer(a) == true && is_nonnegative(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr subfactorial(const Expr& arg) {
    // !n = number of derangements, via the recurrence
    // !0=1, !1=0, !k=(k-1)(!(k-1)+!(k-2)).
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (z.is_negative() || !z.fits_long()) return make<Subfactorial>(arg);
        long n = z.to_long();
        if (n > 100'000) return make<Subfactorial>(arg);  // safety bound
        if (n == 0) return S::One();
        if (n == 1) return S::Zero();
        mpz_class prev2(1), prev1(0), cur;  // !0, !1
        for (long k = 2; k <= n; ++k) {
            cur = mpz_class(k - 1) * (prev1 + prev2);
            prev2 = prev1;
            prev1 = cur;
        }
        return make<Integer>(std::move(cur));
    }
    return make<Subfactorial>(arg);
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
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

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
    // binomial(n, 1) = n  (valid for any n)
    if (k == S::One()) return n;
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
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

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
// Γ'(x) = Γ(x)·ψ(x) = Γ(x)·polygamma(0, x). diff()'s chain rule supplies arg'.
Expr GammaFn::diff_arg(std::size_t /*i*/) const {
    return mul(make<GammaFn>(args_[0]), polygamma(S::Zero(), args_[0]));
}

Expr gamma(const Expr& arg) {
    // gamma(+∞) = +∞ (Γ grows without bound). gamma(−∞) is indeterminate
    // (poles accumulate), so only the positive infinity folds.
    if (arg->type_id() == TypeId::Infinity) return S::Infinity();
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
        // gamma has a simple pole at every non-positive integer → zoo.
        if (!z.is_positive()) return S::ComplexInfinity();
    }
    // Half-integer: gamma(p/2) = C·sqrt(pi) with C exact rational.
    // Bound the recurrence so a huge numerator can never spin.
    if (arg->type_id() == TypeId::Rational) {
        const auto& q = static_cast<const Rational&>(*arg);
        if (q.denominator() == 2) {
            const mpz_class& num = q.numerator();  // odd
            if (num >= -100'001 && num <= 100'001) {
                mpq_class coeff =
                    half_integer_gamma_coeff(mpq_class(num, 2));
                coeff.canonicalize();
                return mul(rational(coeff.get_num(), coeff.get_den()),
                           sqrt(S::Pi()));
            }
        }
    }
    if (arg->type_id() == TypeId::Float) {
        return float_unary_op(mpfr_gamma, arg);
    }
    return make<GammaFn>(arg);
}

// ============================================================================
// Beta
// ============================================================================

namespace {
// True if `e` is or contains an unevaluated gamma(...) node — used to decide
// whether the Β(a,b)=Γ(a)Γ(b)/Γ(a+b) ratio fully reduced to a closed value.
[[nodiscard]] bool has_gamma(const Expr& e) {
    if (e->type_id() == TypeId::Function
        && static_cast<const Function&>(*e).function_id() == FunctionId::Gamma) {
        return true;
    }
    for (const auto& a : e->args()) {
        if (has_gamma(a)) return true;
    }
    return false;
}
}  // namespace

Beta::Beta(Expr a, Expr b)
    : Function(std::vector<Expr>{std::move(a), std::move(b)}) {
    compute_hash(FunctionId::Beta);
}
Expr Beta::rebuild(std::vector<Expr> new_args) const {
    return beta(new_args[0], new_args[1]);
}
std::optional<bool> Beta::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Positive && is_positive(args_[0]) == true
        && is_positive(args_[1]) == true) {
        return true;
    }
    return std::nullopt;
}

Expr beta(const Expr& a, const Expr& b) {
    // Β(a,b) = Γ(a)·Γ(b)/Γ(a+b). Fold to the gamma ratio only when all three
    // gammas evaluate to a closed form (no gamma node left) — e.g. positive
    // integers or half-integers; otherwise keep Β(a,b) symbolic, as SymPy does.
    Expr ga = gamma(a);
    Expr gb = gamma(b);
    Expr gab = gamma(add(a, b));
    if (!has_gamma(ga) && !has_gamma(gb) && !has_gamma(gab)) {
        return mul(mul(ga, gb), pow(gab, S::NegativeOne()));
    }
    return make<Beta>(a, b);
}

// ============================================================================
// Incomplete gamma: lowergamma(s, x), uppergamma(s, x)
// ============================================================================

namespace {
// Closed form of the upper incomplete gamma for a positive integer first
// argument n ≥ 1: Γ(n, x) = (n−1)!·e⁻ˣ·Σ_{k=0}^{n−1} xᵏ/k!. Used by both
// factories; the caller guarantees n ≥ 1 and handles x = 0 separately (the
// series term x⁰ would otherwise hit 0⁰ for a literal zero argument).
[[nodiscard]] Expr upper_gamma_closed(long n, const Expr& x) {
    std::vector<Expr> terms;
    terms.reserve(static_cast<std::size_t>(n));
    mpz_class kfac = 1;  // k!
    for (long k = 0; k < n; ++k) {
        if (k > 0) kfac *= k;
        terms.push_back(mul(pow(x, integer(k)), rational(mpz_class(1), kfac)));
    }
    mpz_class nm1fac;
    mpz_fac_ui(nm1fac.get_mpz_t(), static_cast<unsigned long>(n - 1));
    Expr series = mul(exp(mul(S::NegativeOne(), x)), add(std::move(terms)));
    return mul(make<Integer>(std::move(nm1fac)), std::move(series));
}

// True when s is a positive integer that fits a long and is small enough that
// the (n−1)-term series stays cheap.
[[nodiscard]] std::optional<long> small_positive_int(const Expr& s) {
    if (s->type_id() != TypeId::Integer) return std::nullopt;
    const auto& z = static_cast<const Integer&>(*s);
    if (!z.is_positive() || !z.fits_long()) return std::nullopt;
    long n = z.to_long();
    if (n > 1000) return std::nullopt;
    return n;
}
}  // namespace

LowerGamma::LowerGamma(Expr s, Expr x)
    : Function(std::vector<Expr>{std::move(s), std::move(x)}) {
    compute_hash(FunctionId::LowerGamma);
}
Expr LowerGamma::rebuild(std::vector<Expr> new_args) const {
    return lowergamma(new_args[0], new_args[1]);
}
std::optional<bool> LowerGamma::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_real(args_[0]) == true
        && is_real(args_[1]) == true) {
        return true;
    }
    return std::nullopt;
}
// ∂/∂x γ(s, x) = xˢ⁻¹·e⁻ˣ. The ∂/∂s direction is non-elementary (Meijer-G); as
// with polygamma's order argument, return 0 so diff()'s chain rule (× s′ = 0 for
// a constant s) leaves the usual case correct.
Expr LowerGamma::diff_arg(std::size_t i) const {
    if (i == 1) {
        return mul(pow(args_[1], add(args_[0], S::NegativeOne())),
                   exp(mul(S::NegativeOne(), args_[1])));
    }
    return S::Zero();
}

Expr lowergamma(const Expr& s, const Expr& x) {
    // γ(s, 0) = 0 for any order (SymPy folds this unconditionally).
    if (x == S::Zero()) return S::Zero();
    if (auto n = small_positive_int(s)) {
        // γ(s, x) = Γ(s) − Γ(s, x). Built directly (not via the uppergamma
        // factory) so the x = +∞ form folds to nan like SymPy, rather than Γ(s).
        return add(gamma(s), mul(S::NegativeOne(), upper_gamma_closed(*n, x)));
    }
    return make<LowerGamma>(s, x);
}

UpperGamma::UpperGamma(Expr s, Expr x)
    : Function(std::vector<Expr>{std::move(s), std::move(x)}) {
    compute_hash(FunctionId::UpperGamma);
}
Expr UpperGamma::rebuild(std::vector<Expr> new_args) const {
    return uppergamma(new_args[0], new_args[1]);
}
std::optional<bool> UpperGamma::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_real(args_[0]) == true
        && is_real(args_[1]) == true) {
        return true;
    }
    return std::nullopt;
}
// ∂/∂x Γ(s, x) = −xˢ⁻¹·e⁻ˣ; ∂/∂s is non-elementary → 0 (see LowerGamma).
Expr UpperGamma::diff_arg(std::size_t i) const {
    if (i == 1) {
        return mul(S::NegativeOne(),
                   mul(pow(args_[1], add(args_[0], S::NegativeOne())),
                       exp(mul(S::NegativeOne(), args_[1]))));
    }
    return S::Zero();
}

Expr uppergamma(const Expr& s, const Expr& x) {
    // Γ(s, +∞) = 0 (the tail integral vanishes).
    if (x->type_id() == TypeId::Infinity) return S::Zero();
    if (auto n = small_positive_int(s)) {
        // Γ(s, 0) = Γ(s) = (n−1)!.
        if (x == S::Zero()) return gamma(s);
        return upper_gamma_closed(*n, x);
    }
    return make<UpperGamma>(s, x);
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
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
            if (is_positive(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

// (log Γ)'(x) = ψ(x) = polygamma(0, x).
Expr LogGamma::diff_arg(std::size_t /*i*/) const {
    return polygamma(S::Zero(), args_[0]);
}

Expr loggamma(const Expr& arg) {
    if (arg->type_id() == TypeId::Infinity) return S::Infinity();  // loggamma(oo)=oo
    if (arg->type_id() == TypeId::Float) {
        return float_unary_op(mpfr_lngamma, arg);
    }
    // Pole: loggamma at a nonpositive integer is +∞ (the Γ pole, log of |Γ|→∞).
    if (arg->type_id() == TypeId::Integer
        && !static_cast<const Integer&>(*arg).is_positive()) {
        return S::Infinity();
    }
    // loggamma(x) = log(Γ(x)) for x > 0, where Γ(x) is a positive closed form:
    // log((n−1)!) for a positive integer (loggamma(3)=log 2), log(√π·…) for a
    // positive half-integer (loggamma(½)=log√π). Restricted to positive x — for
    // x < 0, loggamma ≠ log∘Γ (branch cuts), so it is left as loggamma(x), matching
    // SymPy. A symbolic positive p keeps loggamma(p) (Γ(p) does not reduce).
    if (is_positive(arg) == true) {
        Expr g = gamma(arg);
        if (!has_gamma(g)) return log(g);
    }
    return make<LogGamma>(arg);
}

// ============================================================================
// PolyGamma / DiGamma
// ============================================================================

PolyGammaFn::PolyGammaFn(Expr n, Expr x)
    : Function(std::vector<Expr>{std::move(n), std::move(x)}) {
    compute_hash(FunctionId::PolyGamma);
}
Expr PolyGammaFn::rebuild(std::vector<Expr> new_args) const {
    return polygamma(new_args[0], new_args[1]);
}
std::optional<bool> PolyGammaFn::ask(AssumptionKey k) const noexcept {
    // ψ⁽ⁿ⁾(x) is real for real x > 0 and n a non-negative integer.
    if (k == AssumptionKey::Real && is_integer(args_[0]) == true
        && is_nonnegative(args_[0]) == true && is_positive(args_[1]) == true) {
        return true;
    }
    return std::nullopt;
}
// ∂/∂x polygamma(n, x) = polygamma(n+1, x). The ∂/∂n direction is not
// elementary; diff()'s chain rule multiplies by n' (= 0 when n is constant),
// so returning 0 there is harmless for the usual constant-order case.
Expr PolyGammaFn::diff_arg(std::size_t i) const {
    if (i == 1) return polygamma(add(args_[0], S::One()), args_[1]);
    return S::Zero();
}

Expr polygamma(const Expr& n, const Expr& x) {
    // Pole at the nonpositive integers: ψ⁽ⁿ⁾(x) → zoo for x ∈ {0, −1, −2, …} and a
    // non-negative integer order n (the underlying Γ pole). digamma(0)/digamma(−k)
    // inherit this via polygamma(0, ·). Matches SymPy.
    if (x->type_id() == TypeId::Integer
        && !static_cast<const Integer&>(*x).is_positive()
        && n->type_id() == TypeId::Integer
        && !static_cast<const Integer&>(*n).is_negative()) {
        return S::ComplexInfinity();
    }
    // Closed forms at a positive integer or positive half-integer argument, for a
    // non-negative integer order n. Start from the base values at b ∈ {1, 1/2}
    //   ψ⁽⁰⁾(1)   = −γ,                 ψ⁽ⁿ⁾(1)   = (−1)^(n+1)·n!·ζ(n+1),
    //   ψ⁽⁰⁾(1/2) = −γ − 2·log 2,       ψ⁽ⁿ⁾(1/2) = (−1)^(n+1)·n!·(2^(n+1)−1)·ζ(n+1)
    // and walk up with the recurrence ψ⁽ⁿ⁾(y+1) = ψ⁽ⁿ⁾(y) + (−1)ⁿ·n!/y^(n+1):
    //   ψ⁽ⁿ⁾(b+j) = ψ⁽ⁿ⁾(b) + Σ_{i=0}^{j−1} (−1)ⁿ·n!/(b+i)^(n+1).
    // So ψ(m) = −γ + H_{m−1}, ψ(m+1/2) = −γ − 2 log 2 + 2·Σ 1/(2k−1), etc.
    if (n->type_id() == TypeId::Integer
        && !static_cast<const Integer&>(*n).is_negative()) {
        // Resolve x to a base b and a non-negative integer step count j (x = b+j).
        Expr base;
        long j = -1;
        if (x->type_id() == TypeId::Integer) {
            const auto& xi = static_cast<const Integer&>(*x);
            if (xi.is_positive() && xi.fits_long()) {
                base = S::One();
                j = xi.to_long() - 1;
            }
        } else if (x->type_id() == TypeId::Rational) {
            const auto& xr = static_cast<const Rational&>(*x);
            if (xr.denominator() == 2 && xr.numerator() > 0) {
                base = rational(1, 2);
                mpz_class jz = (xr.numerator() - 1) / 2;
                j = jz.get_si();
            }
        }
        if (j >= 0) {
            const bool half = !(base == S::One());
            Expr np1 = add(n, S::One());
            Expr base_val;
            if (n == S::Zero()) {
                base_val = half ? add(mul(S::NegativeOne(), S::EulerGamma()),
                                      mul(integer(-2), log(integer(2))))
                                : mul(S::NegativeOne(), S::EulerGamma());
            } else {
                // (−1)^(n+1) folds to ±1 via the parity rule in the pow factory.
                Expr v = mul(pow(S::NegativeOne(), np1),
                             mul(factorial(n), zeta(np1)));
                if (half) {
                    v = mul(v, add(pow(integer(2), np1), S::NegativeOne()));
                }
                base_val = v;
            }
            std::vector<Expr> terms;
            terms.push_back(base_val);
            for (long i = 0; i < j; ++i) {
                Expr y = add(base, integer(i));
                terms.push_back(mul(pow(S::NegativeOne(), n),
                                    mul(factorial(n),
                                        pow(y, mul(S::NegativeOne(), np1)))));
            }
            return add(std::move(terms));
        }
    }
    return make<PolyGammaFn>(n, x);
}

// ψ(x) = polygamma(0, x). SymPy emits polygamma(0, x) (e.g. from diff loggamma),
// so digamma is sugar for that canonical form.
Expr digamma(const Expr& x) {
    return polygamma(S::Zero(), x);
}

}  // namespace sympp
