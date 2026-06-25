#include <sympp/polys/poly.hpp>

#include <algorithm>
#include <array>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/imaginary_unit.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/infinity.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/matrices/matrix.hpp>

#include <gmpxx.h>

namespace sympp {

namespace {

// Walk an Expr looking for `var^k` factors and return the degree (long) and
// the residual coefficient. Handles a bare `var` (degree 1, coeff 1),
// `var^n` for non-negative Integer n, and Mul of these with other factors
// considered as the coefficient.
// Returns std::nullopt when the term isn't a clean polynomial monomial in `var`.
struct DegCoeff {
    long degree;
    Expr coeff;
};

[[nodiscard]] std::optional<DegCoeff> term_degree_coeff(const Expr& term, const Expr& var) {
    if (term == var) return DegCoeff{1, S::One()};
    if (term->type_id() == TypeId::Pow) {
        const auto& base = term->args()[0];
        const auto& exp = term->args()[1];
        if (base == var && exp->type_id() == TypeId::Integer) {
            const auto& z = static_cast<const Integer&>(*exp);
            if (z.is_negative() || !z.fits_long()) return std::nullopt;
            return DegCoeff{z.to_long(), S::One()};
        }
    }
    if (term->type_id() == TypeId::Mul) {
        long deg = 0;
        std::vector<Expr> coeff_factors;
        for (const auto& a : term->args()) {
            auto sub = term_degree_coeff(a, var);
            if (!sub.has_value()) {
                // The factor doesn't involve var as a clean monomial — treat
                // it as a coefficient if it doesn't contain var, else fail.
                // Quick check: if `var` doesn't appear in `a`, accept it.
                // (We don't have free_symbols-with-comparison here; accept as
                // coefficient unconditionally — caller should pre-flatten.)
                coeff_factors.push_back(a);
                continue;
            }
            deg += sub->degree;
            if (!(sub->coeff == S::One())) {
                coeff_factors.push_back(sub->coeff);
            }
        }
        Expr coeff = mul(std::move(coeff_factors));
        return DegCoeff{deg, coeff};
    }
    // Term has no `var` — treat as constant coefficient at degree 0.
    return DegCoeff{0, term};
}

}  // namespace

Poly::Poly(std::vector<Expr> coeffs, Expr var)
    : coeffs_(std::move(coeffs)), var_(std::move(var)) {
    if (coeffs_.empty()) coeffs_.push_back(S::Zero());
    trim_leading_zeros();
}

Poly::Poly(const Expr& expr, const Expr& var) : var_(var) {
    // Collect terms.
    std::vector<Expr> terms;
    if (expr->type_id() == TypeId::Add) {
        for (const auto& a : expr->args()) terms.push_back(a);
    } else {
        terms.push_back(expr);
    }
    for (const auto& t : terms) {
        auto dc = term_degree_coeff(t, var);
        if (!dc.has_value()) {
            // Couldn't parse cleanly; treat the whole term as a deg-0 coeff.
            if (coeffs_.empty()) coeffs_.push_back(S::Zero());
            coeffs_[0] = add(coeffs_[0], t);
            continue;
        }
        if (static_cast<std::size_t>(dc->degree) >= coeffs_.size()) {
            coeffs_.resize(static_cast<std::size_t>(dc->degree) + 1, S::Zero());
        }
        coeffs_[static_cast<std::size_t>(dc->degree)] =
            add(coeffs_[static_cast<std::size_t>(dc->degree)], dc->coeff);
    }
    if (coeffs_.empty()) coeffs_.push_back(S::Zero());
    trim_leading_zeros();
}

void Poly::trim_leading_zeros() {
    while (coeffs_.size() > 1 && coeffs_.back() == S::Zero()) {
        coeffs_.pop_back();
    }
}

std::size_t Poly::degree() const noexcept {
    return coeffs_.empty() ? 0 : coeffs_.size() - 1;
}

const Expr& Poly::leading_coeff() const {
    if (coeffs_.empty()) {
        throw std::out_of_range("Poly: no coefficients");
    }
    return coeffs_.back();
}

const Expr& Poly::constant_term() const {
    if (coeffs_.empty()) {
        throw std::out_of_range("Poly: no coefficients");
    }
    return coeffs_.front();
}

Poly Poly::operator+(const Poly& other) const {
    std::size_t n = std::max(coeffs_.size(), other.coeffs_.size());
    std::vector<Expr> out(n, S::Zero());
    for (std::size_t i = 0; i < n; ++i) {
        Expr a = i < coeffs_.size() ? coeffs_[i] : S::Zero();
        Expr b = i < other.coeffs_.size() ? other.coeffs_[i] : S::Zero();
        out[i] = add(a, b);
    }
    return Poly{std::move(out), var_};
}

Poly Poly::operator-(const Poly& other) const {
    std::size_t n = std::max(coeffs_.size(), other.coeffs_.size());
    std::vector<Expr> out(n, S::Zero());
    for (std::size_t i = 0; i < n; ++i) {
        Expr a = i < coeffs_.size() ? coeffs_[i] : S::Zero();
        Expr b = i < other.coeffs_.size() ? other.coeffs_[i] : S::Zero();
        out[i] = a - b;
    }
    return Poly{std::move(out), var_};
}

Poly Poly::operator*(const Poly& other) const {
    std::size_t n = coeffs_.size() + other.coeffs_.size() - 1;
    std::vector<Expr> out(n, S::Zero());
    for (std::size_t i = 0; i < coeffs_.size(); ++i) {
        for (std::size_t j = 0; j < other.coeffs_.size(); ++j) {
            out[i + j] = add(out[i + j], mul(coeffs_[i], other.coeffs_[j]));
        }
    }
    return Poly{std::move(out), var_};
}

bool Poly::is_zero() const {
    return coeffs_.size() == 1 && coeffs_[0] == S::Zero();
}

Poly Poly::monic() const {
    if (is_zero()) return *this;
    const Expr& lc = coeffs_.back();
    if (lc == S::One()) return *this;
    std::vector<Expr> out;
    out.reserve(coeffs_.size());
    for (const auto& c : coeffs_) {
        out.push_back(c / lc);
    }
    return Poly{std::move(out), var_};
}

Poly Poly::diff() const {
    if (coeffs_.size() <= 1) {
        return Poly{std::vector<Expr>{S::Zero()}, var_};
    }
    std::vector<Expr> out;
    out.reserve(coeffs_.size() - 1);
    for (std::size_t i = 1; i < coeffs_.size(); ++i) {
        out.push_back(mul(integer(static_cast<long>(i)), coeffs_[i]));
    }
    return Poly{std::move(out), var_};
}

std::pair<Poly, Poly> Poly::divmod(const Poly& other) const {
    if (other.is_zero()) {
        throw std::invalid_argument("Poly::divmod: division by zero polynomial");
    }
    if (!(var_ == other.var_)) {
        throw std::invalid_argument("Poly::divmod: variable mismatch");
    }
    Poly r = *this;
    if (r.degree() < other.degree() || r.is_zero()) {
        return {Poly{std::vector<Expr>{S::Zero()}, var_}, std::move(r)};
    }
    const std::size_t qd = r.degree() - other.degree();
    std::vector<Expr> q(qd + 1, S::Zero());
    const Expr& lc_other = other.leading_coeff();
    while (!r.is_zero() && r.degree() >= other.degree()) {
        const std::size_t prev_degree = r.degree();
        const std::size_t shift = prev_degree - other.degree();
        Expr factor = expand(r.leading_coeff() / lc_other);
        q[shift] = factor;
        // Subtract factor·x^shift·other, expanding each coefficient so that
        // symbolic cancellations fold to a structural zero. Without expand,
        // (b+b²) − (b+b²) stays as an unmerged Add (the bare Add flattens but
        // the −1·Add subtrahend does not), the leading term never zeroes, and
        // the loop never terminates.
        std::vector<Expr> rc = r.coeffs_;
        for (std::size_t i = 0; i < other.coeffs_.size(); ++i) {
            const std::size_t idx = i + shift;
            rc[idx] = expand(rc[idx] - mul(factor, other.coeffs_[i]));
        }
        r = Poly{std::move(rc), var_};  // constructor trims structural zeros
        // Backstop: a valid reduction always drops the degree. If a coefficient
        // could not be cancelled (e.g. a non-polynomial coefficient), stop with
        // the current remainder rather than spinning forever.
        if (!r.is_zero() && r.degree() >= prev_degree) break;
    }
    return {Poly{std::move(q), var_}, std::move(r)};
}

Poly Poly::operator/(const Poly& other) const { return divmod(other).first; }
Poly Poly::operator%(const Poly& other) const { return divmod(other).second; }

Expr Poly::eval(const Expr& at) const {
    // Horner's scheme.
    if (coeffs_.empty()) return S::Zero();
    Expr acc = coeffs_.back();
    for (std::size_t i = coeffs_.size() - 1; i > 0; --i) {
        acc = add(mul(acc, at), coeffs_[i - 1]);
    }
    return acc;
}

Expr Poly::as_expr() const {
    std::vector<Expr> terms;
    terms.reserve(coeffs_.size());
    for (std::size_t i = 0; i < coeffs_.size(); ++i) {
        if (coeffs_[i] == S::Zero()) continue;
        if (i == 0) {
            terms.push_back(coeffs_[i]);
        } else if (i == 1) {
            terms.push_back(mul(coeffs_[i], var_));
        } else {
            terms.push_back(mul(coeffs_[i], pow(var_, integer(static_cast<long>(i)))));
        }
    }
    return add(std::move(terms));
}

std::string Poly::str() const {
    return "Poly(" + as_expr()->str() + ", " + var_->str() + ")";
}

Poly gcd(const Poly& a, const Poly& b) {
    if (!(a.var() == b.var())) {
        throw std::invalid_argument("gcd(Poly, Poly): variable mismatch");
    }
    if (a.is_zero() && b.is_zero()) {
        return Poly{std::vector<Expr>{S::Zero()}, a.var()};
    }
    if (a.is_zero()) return b.monic();
    if (b.is_zero()) return a.monic();
    Poly p = a;
    Poly q = b;
    while (!q.is_zero()) {
        Poly r = p % q;
        p = std::move(q);
        q = std::move(r);
    }
    return p.monic();
}

namespace {
[[nodiscard]] bool poly_is_one(const Poly& p) {
    return p.coeffs().size() == 1 && p.coeffs()[0] == S::One();
}

// Extract a Number coefficient as mpq_class. Returns nullopt if the
// coefficient isn't a numeric Integer/Rational.
[[nodiscard]] std::optional<mpq_class> coeff_as_mpq(const Expr& e) {
    if (!e) return std::nullopt;
    if (e->type_id() == TypeId::Integer) {
        return mpq_class(static_cast<const Integer&>(*e).value());
    }
    if (e->type_id() == TypeId::Rational) {
        return static_cast<const Rational&>(*e).value();
    }
    return std::nullopt;
}

// Build an Expr from an mpq_class (canonicalized). Returns Integer when
// denominator is 1.
[[nodiscard]] Expr mpq_to_expr(const mpq_class& q) {
    mpq_class r = q;
    r.canonicalize();
    if (r.get_den() == 1) return make<Integer>(r.get_num());
    return make<Rational>(r);
}

// All positive divisors of |n|, ascending. Empty if n == 0.
[[nodiscard]] std::vector<mpz_class> positive_divisors(const mpz_class& n) {
    std::vector<mpz_class> ds;
    if (n == 0) return ds;
    mpz_class abs_n = abs(n);
    for (mpz_class i = 1; i * i <= abs_n; ++i) {
        if (abs_n % i == 0) {
            ds.push_back(i);
            mpz_class other = abs_n / i;
            if (other != i) ds.push_back(other);
        }
    }
    std::sort(ds.begin(), ds.end());
    return ds;
}

// Try to interpret `f` as having all-rational coefficients; return them as
// mpq_class. Returns nullopt if any coefficient is non-numeric.
[[nodiscard]] std::optional<std::vector<mpq_class>>
poly_coeffs_as_mpq(const Poly& f) {
    std::vector<mpq_class> out;
    out.reserve(f.coeffs().size());
    for (const auto& c : f.coeffs()) {
        auto q = coeff_as_mpq(c);
        if (!q) return std::nullopt;
        out.push_back(*q);
    }
    return out;
}

// Multiply through by LCM of denominators so that integer coefficients
// remain. Returns the integer-coefficient list (highest degree last) and
// preserves the sign of the leading coefficient.
[[nodiscard]] std::vector<mpz_class>
clear_denominators(const std::vector<mpq_class>& qs) {
    mpz_class lcm = 1;
    for (const auto& q : qs) {
        mpz_class d = q.get_den();
        mpz_class g;
        mpz_gcd(g.get_mpz_t(), lcm.get_mpz_t(), d.get_mpz_t());
        lcm = (lcm / g) * d;
    }
    std::vector<mpz_class> out;
    out.reserve(qs.size());
    for (const auto& q : qs) {
        mpq_class r = q * mpq_class(lcm);
        out.push_back(r.get_num());
    }
    return out;
}

}  // namespace

std::vector<Expr> rational_roots(const Poly& f) {
    auto qs = poly_coeffs_as_mpq(f);
    if (!qs) return {};
    if (qs->size() == 1) return {};  // constant
    auto zs = clear_denominators(*qs);
    // Strip out powers of x (zero roots).
    std::vector<Expr> roots;
    while (!zs.empty() && zs.front() == 0) {
        roots.push_back(S::Zero());
        zs.erase(zs.begin());
    }
    if (zs.size() < 2) return roots;  // constant remains, no more roots

    // Now zs[0] != 0. Find rational roots p/q.
    // p divides |zs[0]|, q divides |zs.back()|.
    auto p_divs = positive_divisors(zs.front());
    auto q_divs = positive_divisors(zs.back());

    // Helper: evaluate polynomial at p/q (rational), check if zero.
    auto eval_at = [&](const std::vector<mpz_class>& coeffs,
                       const mpq_class& at) -> mpq_class {
        mpq_class acc = 0;
        mpq_class pow_at = 1;
        for (const auto& c : coeffs) {
            mpq_class term = mpq_class(c) * pow_at;
            acc += term;
            pow_at *= at;
        }
        return acc;
    };

    // Synthetic division: given coeffs (lowest-degree first) and a
    // candidate root r = p/q, divide poly(x) by (x - r). Returns the
    // reduced coefficient list (one shorter) — assumes r really is a
    // root of the integer-cleared poly multiplied by appropriate factor.
    // Operates on mpq.
    auto divide_out_root = [](std::vector<mpq_class> coeffs,
                              const mpq_class& r) -> std::vector<mpq_class> {
        // coeffs is lowest-degree first: a0 + a1 x + a2 x² + ... + an x^n
        // dividing by (x - r): synthetic division on the high-to-low form.
        if (coeffs.size() < 2) return {};
        std::vector<mpq_class> out(coeffs.size() - 1);
        // Synthetic division highest-degree first.
        out.back() = coeffs.back();
        for (std::size_t i = coeffs.size() - 1; i-- > 1;) {
            out[i - 1] = coeffs[i] + r * out[i];
        }
        return out;
    };

    // Test all candidates p/q with both signs against the current
    // coefficient list; deflate when found.
    std::vector<mpq_class> current(qs->begin() + static_cast<std::ptrdiff_t>(roots.size()),
                                   qs->end());
    bool changed = true;
    while (changed) {
        changed = false;
        // Re-clear denominators after each deflation for divisor enumeration.
        if (current.size() < 2) break;
        auto z = clear_denominators(current);
        auto pds = positive_divisors(z.front());
        auto qds = positive_divisors(z.back());
        for (const auto& p : pds) {
            for (const auto& q : qds) {
                for (int sgn : {1, -1}) {
                    mpq_class cand(sgn * p, q);
                    cand.canonicalize();
                    mpq_class val = eval_at(z, cand);
                    if (val == 0) {
                        roots.push_back(mpq_to_expr(cand));
                        current = divide_out_root(current, cand);
                        changed = true;
                        goto next_pass;
                    }
                }
            }
        }
        next_pass:;
    }

    return roots;
}

namespace {

// Solve a depressed cubic t^3 + p*t + q = 0 in radicals using Cardano.
// Returns three (possibly complex) roots in t. Handles p == 0 and q == 0
// specially.
[[nodiscard]] std::vector<Expr> depressed_cubic_roots(const Expr& p, const Expr& q) {
    if (p == S::Zero()) {
        // t^3 = -q  ⇒  t = ∛(-q) ω_k
        Expr base = sqrt(S::Zero());  // placeholder
        Expr cbrt_neg_q = pow(mul(S::NegativeOne(), q), rational(1, 3));
        Expr omega = rational(-1, 2) + S::I() * sqrt(integer(3)) / integer(2);
        Expr omega2 = rational(-1, 2) - S::I() * sqrt(integer(3)) / integer(2);
        return {cbrt_neg_q, omega * cbrt_neg_q, omega2 * cbrt_neg_q};
    }
    if (q == S::Zero()) {
        // t(t² + p) = 0  ⇒  t = 0, ±√(-p)
        Expr s = sqrt(mul(S::NegativeOne(), p));
        return {S::Zero(), s, mul(S::NegativeOne(), s)};
    }
    // General case: u^3 = -q/2 + sqrt((q/2)² + (p/3)³)
    Expr q_half = q / integer(2);
    Expr disc = pow(q_half, integer(2)) + pow(p / integer(3), integer(3));
    Expr u = pow(mul(S::NegativeOne(), q_half) + sqrt(disc), rational(1, 3));
    // v = -p / (3u)  (uses uv = -p/3 constraint)
    Expr v = mul(S::NegativeOne(), p) / (integer(3) * u);
    Expr omega = rational(-1, 2) + S::I() * sqrt(integer(3)) / integer(2);
    Expr omega2 = rational(-1, 2) - S::I() * sqrt(integer(3)) / integer(2);
    return {u + v, omega * u + omega2 * v, omega2 * u + omega * v};
}

// Solve a*x^3 + b*x^2 + c*x + d = 0 by depressing then Cardano.
[[nodiscard]] std::vector<Expr> cubic_roots(const Expr& a, const Expr& b,
                                            const Expr& c, const Expr& d) {
    // Depress: x = t - b/(3a). Then t^3 + p*t + q = 0 where
    //   p = (3ac - b²) / (3a²)
    //   q = (2b³ - 9abc + 27a²d) / (27a³)
    Expr p = (integer(3) * a * c - b * b) / (integer(3) * a * a);
    Expr q = (integer(2) * b * b * b - integer(9) * a * b * c
              + integer(27) * a * a * d) / (integer(27) * a * a * a);
    auto ts = depressed_cubic_roots(p, q);
    Expr shift = b / (integer(3) * a);
    std::vector<Expr> xs;
    xs.reserve(3);
    for (auto& t : ts) xs.push_back(t - shift);
    return xs;
}

// Solve a depressed quartic t^4 + p*t² + q*t + r = 0 via Ferrari's resolvent.
[[nodiscard]] std::vector<Expr> depressed_quartic_roots(const Expr& p,
                                                        const Expr& q,
                                                        const Expr& r) {
    if (q == S::Zero()) {
        // Biquadratic: t^4 + pt² + r = 0 — solve as quadratic in u = t².
        // u = (-p ± √(p² - 4r)) / 2; then t = ±√u.
        Expr d = pow(p, integer(2)) - integer(4) * r;
        Expr u1 = (mul(S::NegativeOne(), p) + sqrt(d)) / integer(2);
        Expr u2 = (mul(S::NegativeOne(), p) - sqrt(d)) / integer(2);
        Expr s1 = sqrt(u1);
        Expr s2 = sqrt(u2);
        return {s1, mul(S::NegativeOne(), s1), s2, mul(S::NegativeOne(), s2)};
    }
    // Resolvent cubic: 8y³ + 8p*y² + (2p² - 8r)*y - q² = 0.
    // Find any root y.
    Expr a3 = integer(8);
    Expr b3 = integer(8) * p;
    Expr c3 = integer(2) * pow(p, integer(2)) - integer(8) * r;
    Expr d3 = mul(S::NegativeOne(), pow(q, integer(2)));
    auto ys = cubic_roots(a3, b3, c3, d3);
    Expr y = ys.front();  // any root works for Ferrari

    // With y chosen so that 2y > 0 (or treat as complex), the depressed
    // quartic factors as
    //   (t² + p/2 + y)² = 2y t² - q t + (y² + p y + p²/4 - r)
    // The RHS is the perfect square (√(2y) t - q/(2√(2y)))² when y is the
    // resolvent root, so
    //   t² + p/2 + y = ±(√(2y) t - q/(2√(2y)))
    // giving two quadratics; each contributes two roots.
    Expr sqrt_2y = sqrt(integer(2) * y);
    Expr alpha = sqrt_2y;                           // coeff of t in RHS
    Expr beta = q / (integer(2) * sqrt_2y);         // constant in RHS (with sign)
    // Quadratic A: t² - α t + (p/2 + y + β) = 0
    // Quadratic B: t² + α t + (p/2 + y - β) = 0
    Expr half_p_plus_y = p / integer(2) + y;

    auto quad_roots = [](const Expr& A, const Expr& B, const Expr& C)
        -> std::pair<Expr, Expr> {
        Expr disc = pow(B, integer(2)) - integer(4) * A * C;
        Expr root = sqrt(disc);
        Expr two_a = integer(2) * A;
        return {(mul(S::NegativeOne(), B) + root) / two_a,
                (mul(S::NegativeOne(), B) - root) / two_a};
    };
    auto [r1, r2] = quad_roots(S::One(), mul(S::NegativeOne(), alpha),
                               half_p_plus_y + beta);
    auto [r3, r4] = quad_roots(S::One(), alpha, half_p_plus_y - beta);
    return {r1, r2, r3, r4};
}

// ---------------------------------------------------------------------------
// Kronecker factorization helpers
// ---------------------------------------------------------------------------

// Evaluate an integer-coefficient poly (lowest-degree first) at integer x.
[[nodiscard]] mpz_class eval_int_poly(const std::vector<mpz_class>& coeffs,
                                      const mpz_class& at) {
    mpz_class acc = 0;
    mpz_class power = 1;
    for (const auto& c : coeffs) {
        acc += c * power;
        power *= at;
    }
    return acc;
}

// Lagrange interpolation. Given samples (x_i, y_i), return rational
// coefficients (lowest-degree first) of the interpolating polynomial of
// degree at most n - 1 where n = samples.size().
[[nodiscard]] std::optional<std::vector<mpq_class>>
lagrange_interpolate(const std::vector<std::pair<mpz_class, mpz_class>>& samples) {
    const std::size_t n = samples.size();
    if (n == 0) return std::nullopt;
    std::vector<mpq_class> result(n, mpq_class(0));
    for (std::size_t i = 0; i < n; ++i) {
        // Denominator: Π_{j≠i} (x_i - x_j)
        mpz_class denom = 1;
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            denom *= (samples[i].first - samples[j].first);
        }
        if (denom == 0) return std::nullopt;
        // Numerator polynomial: Π_{j≠i} (x - x_j), coefficients lowest-first.
        std::vector<mpq_class> num{mpq_class(1)};
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            std::vector<mpq_class> next(num.size() + 1, mpq_class(0));
            for (std::size_t a = 0; a < num.size(); ++a) {
                next[a] -= num[a] * mpq_class(samples[j].first);
                next[a + 1] += num[a];
            }
            num = std::move(next);
        }
        mpq_class scale = mpq_class(samples[i].second) / mpq_class(denom);
        scale.canonicalize();
        for (std::size_t a = 0; a < num.size(); ++a) {
            result[a] += num[a] * scale;
            result[a].canonicalize();
        }
    }
    return result;
}

// Convert mpq coefficients to a Poly with Integer coefficients. Returns
// nullopt if any coefficient is non-integer.
[[nodiscard]] std::optional<Poly>
mpq_coeffs_to_z_poly(const std::vector<mpq_class>& coeffs, const Expr& var) {
    std::vector<Expr> es;
    es.reserve(coeffs.size());
    for (const auto& q : coeffs) {
        mpq_class c = q;
        c.canonicalize();
        if (c.get_den() != 1) return std::nullopt;
        es.push_back(make<Integer>(c.get_num()));
    }
    return Poly{std::move(es), var};
}

// Try to find a non-trivial integer-coefficient factor of f via Kronecker.
// f is assumed to have rational coefficients. Returns nullopt if no factor
// of degree in [1, deg(f)/2] exists with integer coefficients.
[[nodiscard]] std::optional<Poly> kronecker_find_factor(const Poly& f) {
    if (f.degree() < 2) return std::nullopt;
    auto qcoeffs = poly_coeffs_as_mpq(f);
    if (!qcoeffs) return std::nullopt;
    auto zcoeffs = clear_denominators(*qcoeffs);

    const long n = static_cast<long>(f.degree());
    const long max_k = n / 2;

    for (long k = 1; k <= max_k; ++k) {
        // Pick k+1 distinct integer points where f does not vanish.
        std::vector<mpz_class> xs;
        std::vector<mpz_class> ys;
        long candidate = 0;
        long sign = 1;
        std::size_t guard = 0;
        while (static_cast<long>(xs.size()) < k + 1) {
            mpz_class x_pt(candidate * sign);
            mpz_class y = eval_int_poly(zcoeffs, x_pt);
            if (y != 0) {
                xs.push_back(x_pt);
                ys.push_back(y);
            }
            // Spiral: 0, 1, -1, 2, -2, 3, -3, ...
            if (sign == 1 && candidate != 0) sign = -1;
            else { sign = 1; ++candidate; }
            if (++guard > 200) return std::nullopt;
        }

        // Build divisor-with-sign lists per sample.
        std::vector<std::vector<mpz_class>> divs(static_cast<std::size_t>(k + 1));
        bool any_huge = false;
        for (long i = 0; i <= k; ++i) {
            auto pos = positive_divisors(ys[static_cast<std::size_t>(i)]);
            // Cap divisor count to keep enumeration tractable.
            if (pos.size() > 32) { any_huge = true; break; }
            std::vector<mpz_class>& slot = divs[static_cast<std::size_t>(i)];
            slot.reserve(2 * pos.size());
            for (auto& d : pos) slot.push_back(d);
            for (auto& d : pos) slot.push_back(-d);
        }
        if (any_huge) continue;

        // Enumerate all index combinations.
        std::vector<std::size_t> indices(static_cast<std::size_t>(k + 1), 0);
        while (true) {
            std::vector<std::pair<mpz_class, mpz_class>> samples(
                static_cast<std::size_t>(k + 1));
            for (long i = 0; i <= k; ++i) {
                samples[static_cast<std::size_t>(i)] = {
                    xs[static_cast<std::size_t>(i)],
                    divs[static_cast<std::size_t>(i)][indices[static_cast<std::size_t>(i)]]};
            }
            if (auto interp = lagrange_interpolate(samples)) {
                if (auto candidate_poly = mpq_coeffs_to_z_poly(*interp, f.var())) {
                    if (candidate_poly->degree() == static_cast<std::size_t>(k)
                        && !candidate_poly->is_zero()) {
                        // Trial-divide f (over ℚ) by the candidate.
                        auto [q, r] = f.divmod(*candidate_poly);
                        if (r.is_zero()) {
                            return *candidate_poly;
                        }
                    }
                }
            }
            // Increment index combination.
            long c = k;
            while (c >= 0) {
                ++indices[static_cast<std::size_t>(c)];
                if (indices[static_cast<std::size_t>(c)]
                    < divs[static_cast<std::size_t>(c)].size()) break;
                indices[static_cast<std::size_t>(c)] = 0;
                --c;
            }
            if (c < 0) break;
        }
    }
    return std::nullopt;
}

// Split f into a non-trivial factor g (over ℚ) and its cofactor, returning the
// union of their roots; nullopt if f is irreducible over ℚ. Lets a quartic such
// as x⁴+x²+1 = (x²+x+1)(x²−x+1) yield clean quadratic-formula roots instead of
// Ferrari's nested radicals — matching SymPy, which factors before solving.
[[nodiscard]] std::optional<std::vector<Expr>> roots_by_factoring(const Poly& f) {
    auto g = kronecker_find_factor(f);
    if (!g) return std::nullopt;
    auto [cofactor, rem] = f.divmod(*g);
    if (!rem.is_zero()) return std::nullopt;  // defensive: g must divide f
    std::vector<Expr> all = g->roots();
    for (auto& r : cofactor.roots()) all.push_back(std::move(r));
    return all;
}

// Returns (primitive_integer_poly, scalar) such that scalar * primitive == g
// (as Exprs after as_expr), with primitive having integer coefficients whose
// gcd is 1 (and positive leading coefficient when possible).
[[nodiscard]] std::pair<Poly, mpq_class> primitive_part(const Poly& g) {
    auto qs = poly_coeffs_as_mpq(g);
    if (!qs) return {g, mpq_class(1)};
    if (qs->empty()) return {g, mpq_class(1)};
    // Find LCM of denominators.
    mpz_class lcm = 1;
    for (const auto& q : *qs) {
        mpz_class d = q.get_den();
        mpz_class gd;
        mpz_gcd(gd.get_mpz_t(), lcm.get_mpz_t(), d.get_mpz_t());
        lcm = (lcm / gd) * d;
    }
    // Multiply through to integers.
    std::vector<mpz_class> zs;
    zs.reserve(qs->size());
    for (const auto& q : *qs) {
        mpq_class r = q * mpq_class(lcm);
        zs.push_back(r.get_num());
    }
    // GCD of integer coefficients.
    mpz_class g_int = 0;
    for (auto& z : zs) {
        mpz_class abs_z = abs(z);
        mpz_class gn;
        mpz_gcd(gn.get_mpz_t(), g_int.get_mpz_t(), abs_z.get_mpz_t());
        g_int = gn;
    }
    if (g_int == 0) return {g, mpq_class(1)};
    // Make leading coefficient positive.
    int sign = 1;
    if (zs.back() < 0) sign = -1;
    std::vector<Expr> prim_coeffs;
    prim_coeffs.reserve(zs.size());
    for (auto& z : zs) {
        mpz_class divided = z / g_int;
        if (sign == -1) divided = -divided;
        prim_coeffs.push_back(make<Integer>(std::move(divided)));
    }
    Poly prim{std::move(prim_coeffs), g.var()};
    // g_orig = (g_int / lcm) * sign * prim   (because zs = lcm * coeffs(g))
    // and prim = (sign * zs) / g_int = (sign * lcm * coeffs(g)) / g_int.
    // So coeffs(g) = (g_int * sign / lcm) * coeffs(prim).
    mpq_class scalar(sign * g_int, lcm);
    scalar.canonicalize();
    return {prim, scalar};
}

// Solve a*x^4 + b*x^3 + c*x^2 + d*x + e = 0 via depression + Ferrari.
[[nodiscard]] std::vector<Expr> quartic_roots(const Expr& a, const Expr& b,
                                              const Expr& c, const Expr& d,
                                              const Expr& e) {
    // Depress: x = t - b/(4a). Then t^4 + p*t² + q*t + r = 0 where
    //   p = c/a - 3b²/(8a²)
    //   q = d/a - bc/(2a²) + b³/(8a³)
    //   r = e/a - bd/(4a²) + b²c/(16a³) - 3b⁴/(256a⁴)
    Expr a2 = a * a;
    Expr a3 = a2 * a;
    Expr a4 = a2 * a2;
    Expr b2 = b * b;
    Expr b3 = b2 * b;
    Expr b4 = b2 * b2;
    Expr p = c / a - integer(3) * b2 / (integer(8) * a2);
    Expr q = d / a - b * c / (integer(2) * a2) + b3 / (integer(8) * a3);
    Expr r = e / a - b * d / (integer(4) * a2)
             + b2 * c / (integer(16) * a3)
             - integer(3) * b4 / (integer(256) * a4);
    auto ts = depressed_quartic_roots(p, q, r);
    Expr shift = b / (integer(4) * a);
    std::vector<Expr> xs;
    xs.reserve(4);
    for (auto& t : ts) xs.push_back(t - shift);
    return xs;
}

}  // namespace


FactorList factor_list(const Poly& f) {
    auto sqf = sqf_list(f);
    mpq_class running_content;
    {
        // sqf.content is an Expr; convert to mpq if numeric, else leave as 1
        // and just keep the Expr (rare for textbook inputs).
        auto q = coeff_as_mpq(sqf.content);
        running_content = q ? *q : mpq_class(1);
    }
    std::vector<std::pair<Poly, std::size_t>> out;

    // A primitive-part scalar relates a monic root/Kronecker factor to its
    // integer-content form; since that factor is stored with multiplicity `mult`,
    // the scalar enters the content `mult` times — qᵐ, not q. (Without this,
    // (2n+1)² = 4n²+4n+1 mis-factored as 2·(2n+1)², which is 2× too large.)
    auto pow_mq = [](mpq_class b, std::size_t m) {
        mpq_class r(1);
        for (std::size_t i = 0; i < m; ++i) r *= b;
        return r;
    };

    for (auto& [g_monic, mult] : sqf.factors) {
        // g_monic is square-free monic over ℚ. Pull out rational roots first.
        Poly remaining = g_monic;
        auto rats = rational_roots(remaining);
        for (const auto& r : rats) {
            // Build (x - r) — may have rational coefficients; primitivize.
            Poly lin(std::vector<Expr>{mul(S::NegativeOne(), r), S::One()},
                     remaining.var());
            auto [prim, scalar] = primitive_part(lin);
            running_content *= pow_mq(scalar, mult);
            // Multiplicity: rational_roots returns repeats already, so each
            // call here corresponds to one (x - r) factor; track separately.
            // We'll consolidate at the end.
            out.emplace_back(prim, mult);
            // Deflate `remaining` by (x - r) (use the rational factor, not prim).
            remaining = remaining / lin;
        }
        // Now `remaining` has no rational roots. Apply Kronecker for any
        // higher-degree factors.
        while (remaining.degree() >= 2) {
            auto k_factor = kronecker_find_factor(remaining);
            if (!k_factor) {
                // remaining is irreducible.
                auto [prim, scalar] = primitive_part(remaining);
                running_content *= pow_mq(scalar, mult);
                out.emplace_back(prim, mult);
                break;
            }
            auto [prim, scalar] = primitive_part(*k_factor);
            running_content *= pow_mq(scalar, mult);
            out.emplace_back(prim, mult);
            remaining = remaining / *k_factor;
        }
        // Trailing constant from deflation (when leading coefficient ≠ 1)
        // gets absorbed into content.
        if (remaining.degree() == 0 && !remaining.is_zero()) {
            if (auto qc = coeff_as_mpq(remaining.leading_coeff()); qc) {
                running_content *= pow_mq(*qc, mult);
            }
        }
    }

    // Consolidate duplicate (factor, mult) entries (rational_roots may repeat
    // the same root for higher multiplicity — combine same-Poly entries).
    std::vector<std::pair<Poly, std::size_t>> consolidated;
    for (auto& [p, m] : out) {
        bool merged = false;
        for (auto& [q, m2] : consolidated) {
            if (p.as_expr() == q.as_expr()) {
                m2 += m;
                merged = true;
                break;
            }
        }
        if (!merged) consolidated.emplace_back(p, m);
    }
    running_content.canonicalize();
    return FactorList{mpq_to_expr(running_content), std::move(consolidated)};
}

namespace {

// Homogeneous bivariate factorization: a polynomial all of whose monomials share
// the same total degree (x²−y², x²+2xy+y², x³−y³, …) factors by dehomogenizing
// — set the second variable to 1, factor the resulting univariate ℚ-polynomial
// with the existing machinery, then re-homogenize each factor back to its own
// degree. The product is verified to expand to the input, so a non-homogeneous
// or otherwise-unfactorable polynomial is rejected (returns nullopt) rather than
// risking a wrong answer. Covers the common textbook multivariate cases without
// the full Wang/multivariate-GCD port.
[[nodiscard]] std::optional<Expr> factor_homogeneous_bivariate(const Expr& expr,
                                                               const Expr& var) {
    auto syms = free_symbols(expr);
    if (syms.size() != 2 || syms.find(var) == syms.end()) return std::nullopt;
    Expr other;
    for (const auto& s : syms) {
        if (!(s == var)) { other = s; break; }
    }
    if (!other) return std::nullopt;

    // Total (homogeneous) degree n: max over var-powers k of k + deg_other(c_k).
    Poly px(expand(expr), var);
    long n = 0;
    const auto& cs = px.coeffs();
    for (std::size_t k = 0; k < cs.size(); ++k) {
        if (cs[k] == S::Zero()) continue;
        long dk = static_cast<long>(Poly(expand(cs[k]), other).degree());
        n = std::max(n, static_cast<long>(k) + dk);
    }

    // Dehomogenize and factor over ℚ.
    Poly q(expand(subs(expr, other, S::One())), var);
    if (!poly_coeffs_as_mpq(q)) return std::nullopt;
    const long d = static_cast<long>(q.degree());
    FactorList fl = factor_list(q);

    // Re-homogenize a univariate-in-var factor to its own degree:
    //   Σ aₖ·varᵏ  ↦  Σ aₖ·varᵏ·other^(deg−k).
    auto homogenize = [&](const Poly& f) -> Expr {
        const auto& fc = f.coeffs();
        const long fd = static_cast<long>(f.degree());
        std::vector<Expr> terms;
        for (std::size_t k = 0; k < fc.size(); ++k) {
            if (fc[k] == S::Zero()) continue;
            Expr t = fc[k];
            if (k > 0) t = mul(t, pow(var, integer(static_cast<long>(k))));
            long oe = fd - static_cast<long>(k);
            if (oe > 0) t = mul(t, pow(other, integer(oe)));
            terms.push_back(std::move(t));
        }
        if (terms.empty()) return S::Zero();
        return terms.size() == 1 ? terms[0] : add(std::move(terms));
    };

    std::vector<Expr> parts;
    if (!(fl.content == S::One())) parts.push_back(fl.content);
    if (n - d > 0) parts.push_back(pow(other, integer(n - d)));  // pure-other factor
    for (const auto& [fp, m] : fl.factors) {
        Expr hf = homogenize(fp);
        parts.push_back(m == 1 ? hf : pow(hf, integer(static_cast<long>(m))));
    }
    if (parts.empty()) return std::nullopt;
    Expr result = parts.size() == 1 ? parts[0] : mul(parts);

    // Self-verify: only accept a factorization that reproduces the input.
    if (expand(result) == expand(expr)) return result;
    return std::nullopt;
}

// ---------------------------------------------------------------------------
// Cyclotomic decomposition of x^N − 1.
//
// x^N − 1 = ∏_{d | N} Φ_d(x), where Φ_d is the d-th cyclotomic polynomial and
// each Φ_d is irreducible over ℚ. This sidesteps the super-linear Kronecker
// path for the pure binomial x^N − 1 (factor(x^60−1) becomes instant).
//
// Φ_d is built from the standard recurrence
//   Φ_d(x) = (x^d − 1) / ∏_{e | d, e < d} Φ_e(x)        (exact division),
// memoized across calls.
// ---------------------------------------------------------------------------

// Build the Poly x^n − 1 directly (no expansion needed).
[[nodiscard]] Poly poly_x_pow_minus_one(long n, const Expr& var) {
    std::vector<Expr> c(static_cast<std::size_t>(n) + 1, S::Zero());
    c[0] = S::NegativeOne();
    c[static_cast<std::size_t>(n)] = S::One();
    return Poly(std::move(c), var);
}

// Memoized d-th cyclotomic polynomial as a Poly in `var`. The cache key is d;
// the stored var is irrelevant to the coefficients (all ±integers), so we just
// reconstruct a Poly in the requested var from cached coefficients.
[[nodiscard]] const Poly& cyclotomic_poly(long d, const Expr& var) {
    static std::map<long, std::vector<Expr>> cache;  // d -> coeffs (increasing degree)
    static std::map<long, Poly> typed;               // d -> Poly in current var
    // Reuse a typed cache only while the variable matches; otherwise rebuild.
    auto it = typed.find(d);
    if (it != typed.end() && it->second.var() == var) return it->second;

    auto build = [&](long dd) -> Poly {
        if (auto c = cache.find(dd); c != cache.end())
            return Poly(c->second, var);
        // Φ_d = (x^d − 1) / ∏_{e|d, e<d} Φ_e.
        Poly num = poly_x_pow_minus_one(dd, var);
        Poly denom(std::vector<Expr>{S::One()}, var);  // constant 1
        for (long e = 1; e < dd; ++e) {
            if (dd % e == 0) denom = denom * cyclotomic_poly(e, var);
        }
        Poly phi = num / denom;
        cache[dd] = phi.coeffs();
        return phi;
    };

    Poly phi = build(d);
    typed.insert_or_assign(d, phi);
    return typed.at(d);
}

// If `f` is exactly x^N − 1 (monic, constant term −1, all interior coeffs 0,
// N ≥ 1), return ∏_{d|N} Φ_d(x). Otherwise nullopt. Only applied for
// ℚ-coefficient polynomials by the caller.
[[nodiscard]] std::optional<Expr> factor_cyclotomic_binomial(const Poly& f) {
    const std::size_t n = f.degree();
    if (n < 1) return std::nullopt;
    const auto& cs = f.coeffs();
    if (cs.size() != n + 1) return std::nullopt;
    // Leading coeff 1, constant −1, everything else 0.
    if (!(cs[n] == S::One())) return std::nullopt;
    if (!(cs[0] == S::NegativeOne())) return std::nullopt;
    for (std::size_t i = 1; i < n; ++i)
        if (!(cs[i] == S::Zero())) return std::nullopt;

    const long N = static_cast<long>(n);
    std::vector<Expr> terms;
    for (long d = 1; d <= N; ++d) {
        if (N % d == 0) terms.push_back(cyclotomic_poly(d, f.var()).as_expr());
    }
    if (terms.empty()) return std::nullopt;
    if (terms.size() == 1) return terms[0];
    return mul(std::move(terms));
}

}  // namespace

// Bivariate Wang-style factorization (content extraction; quadratic-in-var via
// a perfect-square discriminant; cubic-and-higher via monomial-root deflation),
// defined after as_numer_denom below.
[[nodiscard]] std::optional<Expr> factor_multivariate(const Expr& expr, const Expr& var);

Expr factor(const Expr& expr, const Expr& var) {
    Poly f(expr, var);
    // factor_list (square-free + rational-root + Kronecker) is defined only
    // over ℚ-coefficient polynomials. A genuinely multivariate polynomial has
    // symbolic coefficients; try the homogeneous-bivariate path (verified) and
    // otherwise return the input unfactored. Full multivariate (Wang)
    // factorization is tracked separately.
    if (!poly_coeffs_as_mpq(f)) {
        if (auto h = factor_homogeneous_bivariate(expr, var)) return *h;
        if (auto w = factor_multivariate(expr, var)) return *w;
        return expr;
    }
    // Fast path: the pure binomial x^N − 1 factors into cyclotomic polynomials
    // ∏_{d|N} Φ_d(x), avoiding the super-linear Kronecker search for high N.
    if (auto cyc = factor_cyclotomic_binomial(f)) return *cyc;
    auto fl = factor_list(f);
    std::vector<Expr> terms;
    if (!(fl.content == S::One())) terms.push_back(fl.content);
    for (auto& [p, m] : fl.factors) {
        Expr ge = p.as_expr();
        if (m == 1) terms.push_back(ge);
        else terms.push_back(pow(ge, integer(static_cast<long>(m))));
    }
    if (terms.empty()) return S::One();
    if (terms.size() == 1) return terms[0];
    return mul(std::move(terms));
}

namespace {

struct NumerDenom { Expr numer; Expr denom; };

// Decompose a denominator into base^power factors (powers accumulated by
// structural base equality). Pow(base, +int) and Mul are peeled; anything else
// (a Symbol, an Add like (x+1), a Function, a non-integer power) is an opaque
// base of power 1. Used to combine a sum of fractions over the LCM of the
// denominators instead of their product.
using FactorPowers = std::vector<std::pair<Expr, long>>;

void accumulate_denom_factors(const Expr& d, long mult, FactorPowers& fp) {
    if (d == S::One()) return;
    if (d->type_id() == TypeId::Pow) {
        const Expr& base = d->args()[0];
        const Expr& ex = d->args()[1];
        if (ex->type_id() == TypeId::Integer) {
            const auto& z = static_cast<const Integer&>(*ex);
            if (z.is_positive() && z.fits_long()) {
                accumulate_denom_factors(base, mult * z.to_long(), fp);
                return;
            }
        }
    }
    if (d->type_id() == TypeId::Mul) {
        for (const auto& f : d->args()) accumulate_denom_factors(f, mult, fp);
        return;
    }
    for (auto& [b, p] : fp) {
        if (b == d) { p += mult; return; }
    }
    fp.push_back({d, mult});
}

// Build Π base^power from a factor-power list (skipping zero powers).
[[nodiscard]] Expr build_from_factors(const FactorPowers& fp) {
    std::vector<Expr> factors;
    for (const auto& [b, p] : fp) {
        if (p == 0) continue;
        factors.push_back(p == 1 ? b : pow(b, integer(p)));
    }
    if (factors.empty()) return S::One();
    if (factors.size() == 1) return factors[0];
    return mul(std::move(factors));
}

// Walk an Expr and split into (numer, denom) such that numer/denom == e.
// No reduction beyond what the canonical factories provide — caller should
// run cancel() afterwards to get a reduced fraction.
[[nodiscard]] NumerDenom as_numer_denom(const Expr& e) {
    if (!e) return {S::One(), S::One()};
    auto t = e->type_id();
    if (t == TypeId::Rational) {
        const auto& r = static_cast<const Rational&>(*e);
        return {make<Integer>(r.numerator()), make<Integer>(r.denominator())};
    }
    if (t == TypeId::Pow) {
        const auto& base = e->args()[0];
        const auto& exp = e->args()[1];
        // Negative integer or rational exponent moves the term to the
        // denominator.
        bool neg = false;
        Expr abs_exp;
        if (exp->type_id() == TypeId::Integer) {
            const auto& z = static_cast<const Integer&>(*exp);
            if (z.is_negative()) {
                neg = true;
                abs_exp = make<Integer>(mpz_class(-z.value()));
            }
        } else if (exp->type_id() == TypeId::Rational) {
            const auto& q = static_cast<const Rational&>(*exp);
            if (q.is_negative()) {
                neg = true;
                mpq_class neg_q = -q.value();
                abs_exp = make<Rational>(neg_q);
            }
        }
        if (neg) return {S::One(), pow(base, abs_exp)};
        return {e, S::One()};
    }
    if (t == TypeId::Mul) {
        std::vector<Expr> nums, dens;
        for (const auto& a : e->args()) {
            auto nd = as_numer_denom(a);
            nums.push_back(nd.numer);
            dens.push_back(nd.denom);
        }
        return {mul(std::move(nums)), mul(std::move(dens))};
    }
    if (t == TypeId::Add) {
        std::vector<NumerDenom> parts;
        parts.reserve(e->args().size());
        for (const auto& a : e->args()) parts.push_back(as_numer_denom(a));

        // Common denominator = LCM of the part denominators (max power per
        // base), not their product — so a shared factor isn't duplicated
        // (a/b + c/b → (a+c)/b, not (a·b + b·c)/b²).
        std::vector<FactorPowers> part_fp;
        part_fp.reserve(parts.size());
        FactorPowers lcm;
        for (const auto& p : parts) {
            FactorPowers fp;
            accumulate_denom_factors(p.denom, 1, fp);
            for (const auto& [b, pw] : fp) {
                bool found = false;
                for (auto& [lb, lp] : lcm) {
                    if (lb == b) { if (pw > lp) lp = pw; found = true; break; }
                }
                if (!found) lcm.push_back({b, pw});
            }
            part_fp.push_back(std::move(fp));
        }
        Expr common = build_from_factors(lcm);

        // Numerator term i = numer_i · (LCM / denom_i), the cofactor being the
        // per-base power deficit of denom_i relative to the LCM.
        std::vector<Expr> result_nums;
        result_nums.reserve(parts.size());
        for (std::size_t i = 0; i < parts.size(); ++i) {
            FactorPowers cof;
            for (const auto& [b, lp] : lcm) {
                long dp = 0;
                for (const auto& [db, de] : part_fp[i]) {
                    if (db == b) { dp = de; break; }
                }
                if (lp - dp > 0) cof.push_back({b, lp - dp});
            }
            Expr cofactor = build_from_factors(cof);
            result_nums.push_back(cofactor == S::One()
                                      ? parts[i].numer
                                      : mul(parts[i].numer, cofactor));
        }
        return {add(std::move(result_nums)), common};
    }
    return {e, S::One()};
}

}  // namespace

// Perfect-square root of a univariate polynomial d(y) over ℚ: returns s with
// s² == d, or nullopt. Greedy top-down coefficient solve, then verified.
[[nodiscard]] static std::optional<Expr> perfect_square_poly(const Expr& d, const Expr& y) {
    if (expand(d) == S::Zero()) return S::Zero();  // disc = 0 → repeated root
    Poly dp(expand(d), y);
    long n = static_cast<long>(dp.degree());
    if (n % 2 != 0) return std::nullopt;
    long m = n / 2;
    const auto& dc = dp.coeffs();
    auto at = [&](long i) -> std::optional<mpq_class> {
        if (i < 0 || i >= static_cast<long>(dc.size())) return mpq_class(0);
        return coeff_as_mpq(dc[static_cast<std::size_t>(i)]);
    };
    auto lead = at(n);
    if (!lead || sgn(*lead) <= 0) return std::nullopt;
    // leading coefficient must be a perfect rational square
    mpz_class sn, sd;
    if (!mpz_perfect_square_p(lead->get_num().get_mpz_t()) ||
        !mpz_perfect_square_p(lead->get_den().get_mpz_t())) {
        return std::nullopt;
    }
    mpz_sqrt(sn.get_mpz_t(), lead->get_num().get_mpz_t());
    mpz_sqrt(sd.get_mpz_t(), lead->get_den().get_mpz_t());
    std::vector<mpq_class> s(static_cast<std::size_t>(m) + 1, mpq_class(0));
    s[static_cast<std::size_t>(m)] = mpq_class(sn, sd);
    for (long k = m - 1; k >= 0; --k) {
        mpq_class known = 0;
        for (long i = k + 1; i < m; ++i) {
            long j = m + k - i;
            if (j >= 0 && j <= m) known += s[static_cast<std::size_t>(i)] * s[static_cast<std::size_t>(j)];
        }
        auto dval = at(m + k);
        if (!dval) return std::nullopt;
        mpq_class sk = (*dval - known) / (2 * s[static_cast<std::size_t>(m)]);
        sk.canonicalize();
        s[static_cast<std::size_t>(k)] = sk;
    }
    Expr se = S::Zero();
    for (long k = 0; k <= m; ++k) {
        if (s[static_cast<std::size_t>(k)] == 0) continue;
        Expr c = mpq_to_expr(s[static_cast<std::size_t>(k)]);
        se = add(se, k == 0 ? c : mul(c, pow(y, integer(k))));
    }
    if (!(expand(add(mul(se, se), mul(integer(-1), expand(d)))) == S::Zero())) return std::nullopt;
    return se;
}

// Rational content (gcd of numerators / lcm of denominators) of a polynomial's
// numeric coefficients, always positive.
[[nodiscard]] static mpq_class numeric_content(const Expr& e) {
    Expr ex = expand(e);
    std::vector<Expr> terms;
    if (ex->type_id() == TypeId::Add) {
        for (const auto& t : ex->args()) terms.push_back(t);
    } else {
        terms.push_back(ex);
    }
    mpz_class g = 0, l = 1;
    for (const auto& t : terms) {
        mpq_class c = 1;
        if (t->type_id() == TypeId::Mul) {
            for (const auto& f : t->args()) {
                if (auto q = coeff_as_mpq(f)) c *= *q;
            }
        } else if (auto q = coeff_as_mpq(t)) {
            c = *q;
        }
        mpz_gcd(g.get_mpz_t(), g.get_mpz_t(), c.get_num().get_mpz_t());
        mpz_lcm(l.get_mpz_t(), l.get_mpz_t(), c.get_den().get_mpz_t());
    }
    if (g == 0) return mpq_class(1);
    mpq_class r(g, l);
    r.canonicalize();
    return r;
}

// Search for a closed-form root r(y) of a bivariate primitive polynomial pp in
// `var` (coefficients polynomial in `y`), among the monomial/rational family
//   r = ±(p/q)·yᵏ ,  p | content(constant term),  q | content(leading coeff),
//                     k ∈ [0, deg].
// Each candidate is accepted only when pp(r) expands to exactly zero, so a wrong
// guess is never returned. This catches the textbook bivariate cubics whose
// roots are simple monomials (x³−y³ has root y; x³+y·x²−x−y has roots ±1, −y).
// Search for a root of `pp` (a polynomial in `var` with coefficients polynomial
// in the variables `vars`) that is a LINEAR FORM r = a₀ + Σ cᵢ·vᵢ with small
// integer coefficients aᵢ, cᵢ ∈ {−1, 0, 1}. This catches symmetric multivariate
// factors such as x³ + y³ + z³ − 3xyz, whose root x = −(y+z) is a linear form in
// {y, z} rather than a monomial. The search is bounded to coefficients in
// {−1, 0, 1} over up to two variables (3³·3 ≈ 81 candidates), enough for the
// symmetric cubic and similar textbook cases. A candidate is accepted only when
// pp(r) expands to exactly zero, so a wrong guess is never returned.
[[nodiscard]] static std::optional<Expr> find_linear_form_root(
    const Poly& pp, const std::vector<Expr>& vars) {
    if (pp.coeffs().size() < 2) return std::nullopt;
    if (vars.empty() || vars.size() > 2) return std::nullopt;
    const std::array<int, 3> choices{0, -1, 1};
    // Enumerate constant term a₀ and one coefficient per variable.
    for (int a0 : choices) {
        std::vector<int> c(vars.size(), 0);
        // Iterate the |choices|^|vars| grid of coefficient tuples.
        long total = 1;
        for (std::size_t i = 0; i < vars.size(); ++i) total *= 3;
        for (long idx = 0; idx < total; ++idx) {
            long t = idx;
            bool all_zero = (a0 == 0);
            for (std::size_t i = 0; i < vars.size(); ++i) {
                c[i] = choices[static_cast<std::size_t>(t % 3)];
                t /= 3;
                if (c[i] != 0) all_zero = false;
            }
            if (all_zero) continue;  // r ≡ 0 is the trivial monomial case
            std::vector<Expr> terms;
            if (a0 != 0) terms.push_back(integer(a0));
            for (std::size_t i = 0; i < vars.size(); ++i) {
                if (c[i] == 0) continue;
                terms.push_back(c[i] == 1 ? vars[i]
                                          : mul(integer(c[i]), vars[i]));
            }
            Expr r = terms.size() == 1 ? terms[0] : add(std::move(terms));
            if (expand(pp.eval(r)) == S::Zero()) return r;
        }
    }
    return std::nullopt;
}

[[nodiscard]] static std::optional<Expr> find_bivariate_root(const Poly& pp,
                                                             const Expr& y) {
    const auto& cs = pp.coeffs();
    if (cs.size() < 2) return std::nullopt;
    long deg = static_cast<long>(pp.degree());

    mpq_class cc = numeric_content(cs.front());        // constant-term content
    mpq_class lc = numeric_content(cs.back());         // leading-coeff content
    // Numerator candidates from |num(cc)·den(lc)|, denominator from |num(lc)·den(cc)|.
    mpz_class num_src = cc.get_num() * lc.get_den();
    mpz_class den_src = lc.get_num() * cc.get_den();
    if (num_src == 0) num_src = 1;
    if (den_src == 0) den_src = 1;
    std::vector<mpz_class> nums = positive_divisors(num_src);
    std::vector<mpz_class> dens = positive_divisors(den_src);
    if (nums.empty()) nums.push_back(1);
    if (dens.empty()) dens.push_back(1);

    for (long k = 0; k <= deg; ++k) {
        Expr yk = (k == 0) ? S::One() : pow(y, integer(k));
        for (const auto& p : nums) {
            for (const auto& q : dens) {
                mpq_class mag(p, q);
                mag.canonicalize();
                for (int sign : {+1, -1}) {
                    mpq_class v = sign > 0 ? mag : -mag;
                    Expr r = (k == 0) ? mpq_to_expr(v) : mul(mpq_to_expr(v), yk);
                    if (expand(pp.eval(r)) == S::Zero()) return r;
                }
            }
        }
    }
    return std::nullopt;
}

// ----- multivariate content extraction (general N-variable) ------------------

// Numeric (rational) GCD of two numbers, as an Expr; 0 acts as identity.
[[nodiscard]] static Expr numeric_gcd_expr(const mpq_class& a, const mpq_class& b) {
    if (a == 0) return mpq_to_expr(abs(b));
    if (b == 0) return mpq_to_expr(abs(a));
    // gcd of two rationals = gcd(num)/lcm(den).
    mpz_class gn, ld;
    mpz_gcd(gn.get_mpz_t(), a.get_num().get_mpz_t(), b.get_num().get_mpz_t());
    mpz_lcm(ld.get_mpz_t(), a.get_den().get_mpz_t(), b.get_den().get_mpz_t());
    mpq_class r(gn, ld);
    r.canonicalize();
    return mpq_to_expr(r);
}

// Exact division of expression `num` by `den` over the polynomial ring, treating
// every free symbol as a polynomial variable. Returns the quotient only when the
// division is exact (verified by re-expansion); otherwise nullopt.
[[nodiscard]] static std::optional<Expr> mv_exact_div(const Expr& num,
                                                      const Expr& den) {
    Expr n = expand(num);
    Expr d = expand(den);
    if (d == S::Zero()) return std::nullopt;
    if (n == S::Zero()) return S::Zero();
    if (auto qd = coeff_as_mpq(d)) {
        // Divide every term by the numeric divisor.
        Expr q = expand(mul(n, mpq_to_expr(mpq_class(qd->get_den(), qd->get_num()))));
        return q;
    }
    // Pick a variable present in the divisor.
    Expr v;
    for (const auto& s : free_symbols(d)) { v = s; break; }
    if (!v) return std::nullopt;
    try {
        Poly pn(n, v), pd(d, v);
        auto [quo, rem] = pn.divmod(pd);
        if (!rem.is_zero()) return std::nullopt;
        Expr q = expand(quo.as_expr());
        if (expand(add(mul(q, d), mul(integer(-1), n))) == S::Zero()) return q;
        return std::nullopt;
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

// GCD of two expressions viewed as multivariate polynomials over ℚ. Recursion is
// on the variable list: content (GCD of var-coefficients, computed recursively)
// times the primitive-part GCD (Euclidean in the chosen variable, with each
// remainder made primitive again to keep coefficients polynomial). The result is
// returned without a fixed sign/normalization beyond being content-positive; the
// caller verifies any factorization that uses it. Bounded recursion depth guards
// against pathological inputs.
[[nodiscard]] static Expr mv_gcd(const Expr& a_in, const Expr& b_in, int depth);

// Content (GCD of all coefficients) of `p` viewed in variable `v`.
[[nodiscard]] static Expr mv_content(const Poly& p, int depth) {
    Expr g = S::Zero();
    for (const auto& c : p.coeffs()) {
        if (expand(c) == S::Zero()) continue;
        g = mv_gcd(g, c, depth + 1);
    }
    if (g == S::Zero()) return S::One();
    return g;
}

[[nodiscard]] static Expr mv_gcd(const Expr& a_in, const Expr& b_in, int depth) {
    Expr a = expand(a_in), b = expand(b_in);
    if (a == S::Zero()) return b == S::Zero() ? S::One() : b;
    if (b == S::Zero()) return a;
    if (depth > 8) return S::One();
    auto qa = coeff_as_mpq(a);
    auto qb = coeff_as_mpq(b);
    if (qa && qb) return numeric_gcd_expr(*qa, *qb);

    // Choose a variable common to both; if none is shared the polynomial GCD is
    // the GCD of their numeric contents (no shared variable can divide both).
    ExprSet sa, sb;
    for (const auto& s : free_symbols(a)) sa.insert(s);
    for (const auto& s : free_symbols(b)) sb.insert(s);
    Expr v;
    for (const auto& s : sa) {
        if (sb.find(s) != sb.end()) { v = s; break; }
    }
    if (!v) return numeric_gcd_expr(numeric_content(a), numeric_content(b));

    Poly pa(a, v), pb(b, v);
    // Contents and primitive parts in v.
    Expr ca = mv_content(pa, depth), cb = mv_content(pb, depth);
    Expr cont = mv_gcd(ca, cb, depth + 1);
    auto ppa_opt = mv_exact_div(a, ca);
    auto ppb_opt = mv_exact_div(b, cb);
    if (!ppa_opt || !ppb_opt) return cont;  // fall back to the content GCD
    Poly ra(*ppa_opt, v), rb(*ppb_opt, v);
    // Euclidean loop on primitive parts; re-primitivize each remainder so the
    // symbolic coefficients stay polynomial rather than blowing up as fractions.
    for (int guard = 0; guard < 64 && !rb.is_zero(); ++guard) {
        Poly rem = ra % rb;
        if (rem.is_zero()) { ra = rb; rb = rem; break; }
        Expr rc = mv_content(rem, depth);
        if (!(rc == S::One())) {
            if (auto pr = mv_exact_div(rem.as_expr(), rc)) rem = Poly(*pr, v);
        }
        ra = rb;
        rb = rem;
    }
    Expr pp_gcd = ra.degree() == 0 ? S::One() : expand(ra.as_expr());
    // Strip pp_gcd's own numeric content so it is primitive.
    mpq_class k = numeric_content(pp_gcd);
    if (!(k == 1) && k != 0) {
        if (auto s = mv_exact_div(pp_gcd, mpq_to_expr(k))) pp_gcd = *s;
    }
    return expand(mul(cont, pp_gcd));
}

// Exact polynomial square root over ℚ in any number of variables: returns s with
// s² == e (s sign-normalized non-negative), or nullopt. Recurses on one variable
// at a time (a perfect square is a perfect square in each variable), then
// verifies s² == e by re-expansion so a wrong guess is never returned.
[[nodiscard]] static std::optional<Expr> mv_sqrt(const Expr& e_in) {
    Expr e = expand(e_in);
    if (e == S::Zero()) return S::Zero();
    if (auto q = coeff_as_mpq(e)) {
        if (sgn(*q) < 0) return std::nullopt;
        mpz_class sn, sd;
        if (!mpz_perfect_square_p(q->get_num().get_mpz_t())
            || !mpz_perfect_square_p(q->get_den().get_mpz_t()))
            return std::nullopt;
        mpz_sqrt(sn.get_mpz_t(), q->get_num().get_mpz_t());
        mpz_sqrt(sd.get_mpz_t(), q->get_den().get_mpz_t());
        return mpq_to_expr(mpq_class(sn, sd));
    }
    Expr v;
    for (const auto& s : free_symbols(e)) { v = s; break; }
    if (!v) return std::nullopt;
    Poly p(e, v);
    long n = static_cast<long>(p.degree());
    if (n % 2 != 0) return std::nullopt;
    long m = n / 2;
    const auto& pc = p.coeffs();
    auto at = [&](long i) -> Expr {
        if (i < 0 || i >= static_cast<long>(pc.size())) return S::Zero();
        return pc[static_cast<std::size_t>(i)];
    };
    // Leading coefficient must itself be a perfect square (recursively).
    auto sm_opt = mv_sqrt(at(n));
    if (!sm_opt) return std::nullopt;
    std::vector<Expr> s(static_cast<std::size_t>(m) + 1, S::Zero());
    s[static_cast<std::size_t>(m)] = *sm_opt;
    Expr two_sm = expand(mul(integer(2), *sm_opt));
    for (long k = m - 1; k >= 0; --k) {
        Expr known = S::Zero();
        for (long i = k + 1; i < m; ++i) {
            long j = m + k - i;
            if (j >= 0 && j <= m)
                known = add(known, mul(s[static_cast<std::size_t>(i)],
                                       s[static_cast<std::size_t>(j)]));
        }
        Expr numerator = expand(add(at(m + k), mul(integer(-1), known)));
        auto sk = mv_exact_div(numerator, two_sm);
        if (!sk) return std::nullopt;
        s[static_cast<std::size_t>(k)] = *sk;
    }
    Expr se = S::Zero();
    for (long k = 0; k <= m; ++k) {
        if (s[static_cast<std::size_t>(k)] == S::Zero()) continue;
        Expr ck = s[static_cast<std::size_t>(k)];
        se = add(se, k == 0 ? ck : mul(ck, pow(v, integer(k))));
    }
    if (expand(add(mul(se, se), mul(integer(-1), e))) == S::Zero()) return se;
    return std::nullopt;
}

// General multivariate factorization for >2 other variables, covering the two
// patterns reachable without a full Wang port:
//   (1) content extraction: P (degree d in var) shares a non-trivial polynomial
//       factor g(other vars) across all its var-coefficients; pull g out and
//       factor the primitive part recursively. Handles linear/bilinear forms
//       like a·b−a·c+b·d−c·d = (b−c)(a+d).
//   (2) generalized difference of squares: P = A·var² + C with A and −C perfect
//       squares (as polynomials in the other vars) → (√A·var−√(−C))(√A·var+√(−C)).
//       Handles x²y²−z²w² = (xy−wz)(xy+wz).
// Every result is verified to expand back to the input and to be strictly more
// factored; otherwise nullopt.
[[nodiscard]] static std::optional<Expr> factor_multivariate_general(
    const Expr& expr, const Expr& var) {
    Expr ex = expand(expr);
    Poly f(ex, var);
    if (f.degree() < 1) return std::nullopt;
    const auto& cs = f.coeffs();

    // (1) Content extraction across the var-coefficients.
    Expr cont = mv_content(f, 0);
    if (!(cont == S::One()) && !coeff_as_mpq(cont)) {
        if (auto pp = mv_exact_div(ex, cont)) {
            Expr cont_f = factor(cont, var);  // recurse on the content (no var inside)
            Expr pp_f = *pp;
            // Recurse on the primitive part (may itself factor, e.g. via case 2).
            if (auto inner = factor_multivariate_general(pp_f, var)) pp_f = *inner;
            Expr result = mul(cont_f, pp_f);
            if (expand(add(result, mul(integer(-1), ex))) == S::Zero()
                && !(result == expr) && !(result == ex)) {
                return result;
            }
        }
    }

    // (2) Generalized difference of squares: A·var² + C (no linear term).
    if (f.degree() == 2 && expand(cs[1]) == S::Zero()) {
        Expr A = cs[2], C = cs[0];
        Expr negC = expand(mul(integer(-1), C));
        auto sa = mv_sqrt(A);
        auto sc = mv_sqrt(negC);
        if (sa && sc && !(*sa == S::Zero()) && !(*sc == S::Zero())) {
            Expr l1 = expand(add(mul(*sa, var), mul(integer(-1), *sc)));
            Expr l2 = expand(add(mul(*sa, var), *sc));
            Expr result = mul(l1, l2);
            if (expand(add(result, mul(integer(-1), ex))) == S::Zero()
                && !(result == expr) && !(result == ex)) {
                return result;
            }
        }
    }

    return std::nullopt;
}

std::optional<Expr> factor_multivariate(const Expr& expr, const Expr& var) {
    // Collect the other free symbols. One other variable y enables the full
    // bivariate path (content extraction + quadratic disc + monomial deflation);
    // two enable the trivariate linear-form deflation (e.g. the symmetric cubic
    // x³+y³+z³−3xyz, root x = −(y+z)). Beyond two we bail out.
    std::vector<Expr> vars;
    for (const auto& s : free_symbols(expr)) {
        if (!(s == var)) vars.push_back(s);
    }
    if (vars.empty()) return std::nullopt;
    if (vars.size() > 2) return factor_multivariate_general(expr, var);
    const bool bivariate = (vars.size() == 1);
    Expr y = bivariate ? vars[0] : Expr{};

    Expr ex = expand(expr);  // normalize so (x+1)²-style inputs parse as polynomials
    Poly f(ex, var);
    if (f.degree() < 1) return std::nullopt;
    const auto& cs = f.coeffs();

    // Content = gcd of the var-coefficients as polynomials in y; primitive part
    // = the original with that content divided out of each coefficient. Content
    // extraction needs a single other variable (univariate GCD); with two we
    // treat the whole polynomial as already primitive.
    Expr cont_e = S::One();
    std::vector<Expr> ppc(cs.begin(), cs.end());
    if (bivariate) {
        Poly cont(cs[0], y);
        for (std::size_t i = 1; i < cs.size(); ++i) cont = gcd(cont, Poly(cs[i], y));
        cont_e = cont.as_expr();
        ppc.clear();
        ppc.reserve(cs.size());
        for (const auto& c : cs) ppc.push_back((Poly(c, y) / cont).as_expr());
    }
    Poly pp(ppc, var);
    Expr pp_e = pp.as_expr();

    // Factor the primitive part — currently the quadratic-in-var case.
    Expr pp_factored = pp_e;
    if (f.degree() == 2 && bivariate) {
        Expr a = ppc[2], b = ppc[1], c = ppc[0];
        Expr disc = expand(add(mul(b, b), mul(integer(-4), mul(a, c))));
        if (auto sopt = perfect_square_poly(disc, y)) {
            const Expr& srt = *sopt;
            // roots (−b±s)/(2a); clear each to an integer-polynomial linear factor.
            Expr lead = a;  // accumulates a / (den₁·den₂)
            std::vector<Expr> lin;
            for (int sign : {+1, -1}) {
                Expr num = add(mul(integer(-1), b), mul(integer(sign), srt));  // −b ± s
                Expr root = cancel(mul(num, pow(mul(integer(2), a), integer(-1))));  // (−b±s)/(2a)
                NumerDenom nd = as_numer_denom(root);
                Expr factor_lin = expand(add(mul(nd.denom, var), mul(integer(-1), nd.numer)));
                mpq_class k = numeric_content(factor_lin);  // pull numeric content into lead
                if (!(k == 1)) {
                    factor_lin = expand(mul(factor_lin, mpq_to_expr(mpq_class(k.get_den(), k.get_num()))));
                    lead = mul(lead, mpq_to_expr(k));
                }
                lin.push_back(factor_lin);
                lead = mul(lead, pow(nd.denom, integer(-1)));
            }
            Expr cand = cancel(mul(lead, mul(lin[0], lin[1])));
            if (expand(add(cand, mul(integer(-1), pp_e))) == S::Zero()) {
                Expr lead_s = cancel(lead);
                std::vector<Expr> parts;
                if (!(lead_s == S::One())) parts.push_back(lead_s);
                parts.push_back(lin[0]);
                parts.push_back(lin[1]);
                pp_factored = mul(parts);
            }
        }
    } else if (f.degree() >= 3 || (f.degree() == 2 && !bivariate)) {
        // Cubic and higher in `var`: find one closed-form root r(y), deflate by
        // (var − r), and recurse on the (lower-degree) quotient. The linear
        // factor is denominator-cleared to (den·var − num) so it stays a
        // polynomial; the leftover numeric scale is folded into the quotient.
        // For a single other variable the root is sought among monomials ±(p/q)yᵏ;
        // for two it is sought among small linear forms a₀+Σcᵢvᵢ. The degree-2
        // multivariate case (e.g. (x+y)²−z², two other free vars) also takes this
        // linear-form path since the bivariate-only quadratic-disc branch above
        // does not apply.
        std::optional<Expr> ropt = bivariate ? find_bivariate_root(pp, y)
                                             : find_linear_form_root(pp, vars);
        if (ropt) {
            NumerDenom nd = as_numer_denom(*ropt);
            Expr lin = expand(add(mul(nd.denom, var), mul(integer(-1), nd.numer)));
            mpq_class k = numeric_content(lin);
            if (!(k == 1)) {
                lin = expand(mul(lin, mpq_to_expr(mpq_class(k.get_den(), k.get_num()))));
            }
            Poly lin_p(lin, var);
            auto [quo_p, rem_p] = pp.divmod(lin_p);
            if (rem_p.is_zero()) {
                Expr quo_e = quo_p.as_expr();
                // Recurse on the quotient: a cubic deflates to a quadratic that
                // the disc path (or a further deflation) may split.
                Expr quo_factored = quo_e;
                if (auto qf = factor_multivariate(quo_e, var)) quo_factored = *qf;
                Expr cand = expand(mul(lin, quo_factored));
                if (expand(add(cand, mul(integer(-1), pp_e))) == S::Zero()) {
                    std::vector<Expr> parts{lin};
                    if (quo_factored->type_id() == TypeId::Mul) {
                        for (const auto& a : quo_factored->args()) parts.push_back(a);
                    } else {
                        parts.push_back(quo_factored);
                    }
                    pp_factored = mul(parts);
                }
            }
        }
    }

    std::vector<Expr> out;
    if (!(cont_e == S::One())) out.push_back(factor(cont_e, y));  // recurse on the content
    out.push_back(pp_factored);
    Expr result = mul(out);
    // Accept only a verified, genuinely-more-factored result.
    if (!(expand(add(result, mul(integer(-1), ex))) == S::Zero())) return std::nullopt;
    if (result == expr || result == ex) return std::nullopt;
    return result;
}

namespace {
// Recursive together: combine each Add/Mul/Pow child into a single fraction
// first, so a nested compound fraction collapses bottom-up before the current
// level is recombined. This is what lets a reciprocal-of-a-sum-of-fractions
// reduce — together(1/(1+1/x)) recurses into the base (1+1/x → (x+1)/x), and the
// resulting ((x+1)/x)⁻¹ flattens to x/(x+1). Function arguments are left intact
// (together stays shallow inside sin/exp/…), matching SymPy's default.
[[nodiscard]] Expr together_recursive(const Expr& e) {
    if (!e) return e;
    Expr rebuilt = e;
    switch (e->type_id()) {
        case TypeId::Add: {
            std::vector<Expr> a;
            a.reserve(e->args().size());
            for (const auto& c : e->args()) a.push_back(together_recursive(c));
            rebuilt = add(std::move(a));
            break;
        }
        case TypeId::Mul: {
            std::vector<Expr> a;
            a.reserve(e->args().size());
            for (const auto& c : e->args()) a.push_back(together_recursive(c));
            rebuilt = mul(std::move(a));
            break;
        }
        case TypeId::Pow: {
            const Expr& exp = e->args()[1];
            Expr base = together_recursive(e->args()[0]);
            auto bnd = as_numer_denom(base);
            if (bnd.denom == S::One()) {
                rebuilt = pow(base, exp);
            } else {
                // (bn/bd)^exp = bn^exp · bd^(−exp); a negative exp thus flips the
                // fraction so 1/((x+1)/x) → x/(x+1).
                rebuilt = mul(pow(bnd.numer, exp),
                              pow(bnd.denom, mul(integer(-1), exp)));
            }
            break;
        }
        default:
            break;
    }
    auto nd = as_numer_denom(rebuilt);
    if (nd.denom == S::One()) return nd.numer;
    return nd.numer / nd.denom;
}
}  // namespace

Expr together(const Expr& expr) {
    return together_recursive(expr);
}

Expr cancel(const Expr& expr, const Expr& var) {
    auto nd = as_numer_denom(expr);
    Expr num_expanded = expand(nd.numer);
    Expr den_expanded = expand(nd.denom);
    Poly num_p(num_expanded, var);
    Poly den_p(den_expanded, var);
    if (den_p.is_zero()) return expr;  // can't cancel by zero
    auto g = gcd(num_p, den_p);
    Poly reduced_num = num_p / g;
    Poly reduced_den = den_p / g;
    if (reduced_den.degree() == 0
        && reduced_den.coeffs().front() == S::One()) {
        return reduced_num.as_expr();
    }
    return reduced_num.as_expr() / reduced_den.as_expr();
}

namespace {
// The single free symbol of the given expressions, or nullopt if there is not
// exactly one (used to infer the polynomial variable for the parser-facing
// poly operations).
[[nodiscard]] std::optional<Expr> inferred_var(const Expr& a, const Expr& b) {
    ExprSet vars;
    for (const auto& s : free_symbols(a)) vars.insert(s);
    if (b) {
        for (const auto& s : free_symbols(b)) vars.insert(s);
    }
    if (vars.size() != 1) return std::nullopt;
    return *vars.begin();
}
}  // namespace

Expr cancel(const Expr& expr) {
    auto nd = as_numer_denom(expr);
    auto var = inferred_var(nd.numer, nd.denom);
    if (!var) return expr;  // no/multiple vars: nothing univariate to cancel
    return cancel(expr, *var);
}

Expr degree(const Expr& expr) {
    auto var = inferred_var(expr, Expr{});
    if (!var) {
        // A constant: deg(0) = −∞, deg(c≠0) = 0 (matching SymPy).
        if (free_symbols(expr).empty()) {
            return expr == S::Zero() ? S::NegativeInfinity()
                                     : Expr{integer(0)};
        }
        return function_symbol("degree")(std::vector<Expr>{expr});
    }
    try {
        Poly p(expand(expr), *var);
        return integer(static_cast<long>(p.degree()));
    } catch (const std::exception&) {
        return function_symbol("degree")(std::vector<Expr>{expr});
    }
}

Expr quo(const Expr& a, const Expr& b) {
    auto var = inferred_var(a, b);
    if (var) {
        try {
            Poly pa(expand(a), *var);
            Poly pb(expand(b), *var);
            return pa.divmod(pb).first.as_expr();
        } catch (const std::exception&) {
        }
    }
    return function_symbol("quo")(std::vector<Expr>{a, b});
}

Expr rem(const Expr& a, const Expr& b) {
    auto var = inferred_var(a, b);
    if (var) {
        try {
            Poly pa(expand(a), *var);
            Poly pb(expand(b), *var);
            return pa.divmod(pb).second.as_expr();
        } catch (const std::exception&) {
        }
    }
    return function_symbol("rem")(std::vector<Expr>{a, b});
}

Expr resultant(const Expr& p, const Expr& q) {
    auto var = inferred_var(p, q);
    if (var) {
        try {
            return resultant(p, q, *var);
        } catch (const std::exception&) {
        }
    }
    return function_symbol("resultant")(std::vector<Expr>{p, q});
}

Expr discriminant(const Expr& p) {
    auto var = inferred_var(p, Expr{});
    if (var) {
        try {
            return discriminant(p, *var);
        } catch (const std::exception&) {
        }
    }
    return function_symbol("discriminant")(std::vector<Expr>{p});
}

namespace {

// Partial-fraction decomposition of num/den over the irreducible factorization
// of den, including repeated factors: a factor fᵢ of multiplicity mᵢ contributes
// terms Pᵢⱼ/fᵢʲ for j = 1..mᵢ, each with deg(Pᵢⱼ) < deg(fᵢ). Found by undetermined
// coefficients: num = Σᵢ Σⱼ Pᵢⱼ·(den/fᵢʲ) is an N×N rational linear system
// (N = deg den), solved via rref. Returns nullopt for non-rational coefficients,
// a single multiplicity-1 irreducible factor (nothing to split), or an
// impractical degree.
[[nodiscard]] std::optional<Expr> partial_fractions_squarefree(
    const Poly& num_p, const Poly& den_p, const Expr& var) {
    constexpr std::size_t kMaxDegree = 10;
    const std::size_t N = den_p.degree();
    if (N < 2 || N > kMaxDegree) return std::nullopt;
    if (!poly_coeffs_as_mpq(den_p)) return std::nullopt;

    FactorList fl = factor_list(den_p);
    // Nothing to separate: a lone irreducible factor of multiplicity 1.
    if (fl.factors.size() == 1 && fl.factors[0].second == 1) {
        return std::nullopt;
    }

    const Poly one_poly(std::vector<Expr>{S::One()}, var);

    // One unknown per coefficient of each Pᵢⱼ (degree deg(fᵢ)−1) for every power
    // j = 1..mᵢ. The basis polynomial for coefficient p of (factor i, power j) is
    // xᵖ · (den/fᵢʲ).
    struct Unknown { std::size_t factor; std::size_t power; std::size_t coeff; };
    std::vector<Unknown> unknowns;
    std::vector<std::vector<Expr>> basis;  // per unknown: coeff vector, low→high
    for (std::size_t i = 0; i < fl.factors.size(); ++i) {
        const Poly& f = fl.factors[i].first;
        const std::size_t m = fl.factors[i].second;
        // base = den / fᵢ^mᵢ (the part coprime to fᵢ); den/fᵢʲ = base · fᵢ^(mᵢ−j).
        Poly fpow = one_poly;
        for (std::size_t k = 0; k < m; ++k) fpow = fpow * f;
        auto [base, rem] = den_p.divmod(fpow);
        if (!rem.is_zero()) return std::nullopt;  // fᵢ^mᵢ must divide den exactly
        Poly den_over = base;  // starts at j = m (den/fᵢ^m); multiply by f as j↓
        for (std::size_t j = m; j >= 1; --j) {
            const auto& cc = den_over.coeffs();
            for (std::size_t p = 0; p < f.degree(); ++p) {
                unknowns.push_back({i, j, p});
                std::vector<Expr> col(N, S::Zero());
                for (std::size_t k = 0; k < cc.size() && p + k < N; ++k) {
                    col[p + k] = cc[k];
                }
                basis.push_back(std::move(col));
            }
            if (j > 1) den_over = den_over * f;  // den/fᵢ^(j−1)
        }
    }
    if (unknowns.size() != N) return std::nullopt;

    // Solve the augmented system [M | num] by rref; M column = basis vector.
    Matrix aug = Matrix::zeros(N, N + 1);
    for (std::size_t col = 0; col < N; ++col) {
        for (std::size_t row = 0; row < N; ++row) {
            aug.set(row, col, basis[col][row]);
        }
    }
    const auto& nc = num_p.coeffs();
    for (std::size_t row = 0; row < N; ++row) {
        aug.set(row, N, row < nc.size() ? nc[row] : S::Zero());
    }
    auto [R, pivots] = aug.rref();
    if (pivots.size() != N) return std::nullopt;  // singular — bail
    for (std::size_t k = 0; k < N; ++k) {
        if (pivots[k] != k) return std::nullopt;  // not a clean identity block
    }

    // Reassemble each Pᵢⱼ from the solution column and build Σᵢ Σⱼ Pᵢⱼ/fᵢʲ. Key
    // the coefficient buckets by (factor i, power j).
    std::map<std::pair<std::size_t, std::size_t>, std::vector<Expr>> p_coeffs;
    for (std::size_t i = 0; i < fl.factors.size(); ++i) {
        const std::size_t m = fl.factors[i].second;
        for (std::size_t j = 1; j <= m; ++j) {
            p_coeffs[{i, j}].assign(fl.factors[i].first.degree(), S::Zero());
        }
    }
    for (std::size_t k = 0; k < N; ++k) {
        p_coeffs[{unknowns[k].factor, unknowns[k].power}][unknowns[k].coeff] =
            R.at(k, N);
    }
    std::vector<Expr> terms;
    for (auto& [key, coeffs] : p_coeffs) {
        Poly numer(std::move(coeffs), var);
        Expr ne = numer.as_expr();
        if (ne == S::Zero()) continue;
        Expr f_expr = fl.factors[key.first].first.as_expr();
        terms.push_back(ne / pow(f_expr, integer(static_cast<long>(key.second))));
    }
    if (terms.empty()) return S::Zero();
    // A single term means nothing actually split — the input was already a
    // partial fraction (e.g. 1/(x²+1)² or (x+1)/(x²+1)²). Report "no
    // decomposition" so callers don't re-process an unchanged fraction (which
    // would loop in the integration pipeline).
    if (terms.size() == 1) return std::nullopt;
    return add(std::move(terms));
}

}  // namespace

Expr apart(const Expr& expr, const Expr& var) {
    Expr canceled = cancel(expr, var);
    auto nd = as_numer_denom(canceled);
    Poly num_p(expand(nd.numer), var);
    Poly den_p(expand(nd.denom), var);
    if (den_p.is_zero() || den_p.degree() == 0) return canceled;

    // Polynomial part.
    Expr poly_part = S::Zero();
    if (num_p.degree() >= den_p.degree()) {
        auto [q, r] = num_p.divmod(den_p);
        poly_part = q.as_expr();
        num_p = r;
    }
    if (num_p.is_zero()) return poly_part;

    // Distinct linear factor case via Heaviside cover-up.
    auto roots = rational_roots(den_p);
    // Detect repeats.
    std::vector<Expr> distinct;
    bool any_repeat = false;
    for (const auto& r : roots) {
        bool seen = false;
        for (const auto& dr : distinct) {
            if (dr == r) { seen = true; break; }
        }
        if (seen) { any_repeat = true; break; }
        distinct.push_back(r);
    }
    if (any_repeat || distinct.size() != den_p.degree()) {
        // An irreducible quadratic (or higher) factor is present: try full
        // partial fractions over the squarefree factorization.
        if (auto pf = partial_fractions_squarefree(num_p, den_p, var)) {
            return poly_part == S::Zero() ? *pf : poly_part + *pf;
        }
        // Repeated roots, or otherwise outside scope — return unreduced.
        if (poly_part == S::Zero()) return num_p.as_expr() / den_p.as_expr();
        return poly_part + num_p.as_expr() / den_p.as_expr();
    }

    Expr leading = den_p.leading_coeff();
    std::vector<Expr> terms;
    if (!(poly_part == S::Zero())) terms.push_back(poly_part);
    for (std::size_t i = 0; i < distinct.size(); ++i) {
        const Expr& r_i = distinct[i];
        Expr a_num = num_p.eval(r_i);
        Expr a_den = leading;
        for (std::size_t j = 0; j < distinct.size(); ++j) {
            if (i == j) continue;
            a_den = mul(a_den, distinct[i] - distinct[j]);
        }
        Expr a_i = a_num / a_den;
        terms.push_back(a_i / (var - r_i));
    }
    if (terms.empty()) return S::Zero();
    if (terms.size() == 1) return terms[0];
    return add(std::move(terms));
}

Expr horner(const Expr& expr, const Expr& var) {
    Poly p(expr, var);
    return p.eval(var);  // Poly::eval already builds the nested Horner form.
}

Expr resultant(const Poly& p, const Poly& q) {
    if (!(p.var() == q.var())) {
        throw std::invalid_argument("resultant: variable mismatch");
    }
    if (p.is_zero() || q.is_zero()) return S::Zero();
    const std::size_t m = p.degree();
    const std::size_t n = q.degree();
    if (m == 0 && n == 0) {
        // Both constants; resultant convention: 1 (no common roots).
        return S::One();
    }
    const std::size_t size = m + n;
    Matrix syl = Matrix::zeros(size, size);
    // p contributes n rows. Row i covers columns [i, i+m].
    // p.coeffs() is lowest-degree first; Sylvester wants highest-degree first.
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= m; ++j) {
            syl.set(i, i + j, p.coeffs()[m - j]);
        }
    }
    // q contributes m rows. Row n + i covers columns [i, i+n].
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j <= n; ++j) {
            syl.set(n + i, i + j, q.coeffs()[n - j]);
        }
    }
    return syl.det();
}

Expr resultant(const Expr& p, const Expr& q, const Expr& var) {
    return resultant(Poly(p, var), Poly(q, var));
}

Expr discriminant(const Poly& p) {
    if (p.degree() < 1) return S::Zero();
    Expr res = resultant(p, p.diff());
    const std::size_t n = p.degree();
    const Expr& lc = p.leading_coeff();
    // (-1)^(n(n-1)/2) sign factor
    Expr sign = ((n * (n - 1) / 2) % 2 == 0) ? S::One() : S::NegativeOne();
    return sign * res / lc;
}

Expr discriminant(const Expr& p, const Expr& var) {
    return discriminant(Poly(p, var));
}

SqfList sqf_list(const Poly& f) {
    if (f.is_zero()) return SqfList{S::Zero(), {}};
    // Pull out the leading coefficient as content; work on the monic primitive.
    Expr content = f.leading_coeff();
    Poly p = f.monic();
    if (p.degree() == 0) return SqfList{content, {}};

    Poly fprime = p.diff();
    Poly c = gcd(p, fprime);  // monic
    Poly w = p / c;
    Poly y = fprime / c;

    std::vector<std::pair<Poly, std::size_t>> factors;
    std::size_t i = 1;
    while (!poly_is_one(w)) {
        Poly z = y - w.diff();
        Poly g = gcd(w, z);
        if (!poly_is_one(g)) factors.emplace_back(g, i);
        w = w / g;
        y = z / g;
        ++i;
    }
    return SqfList{content, std::move(factors)};
}

std::vector<Expr> Poly::roots() const {
    if (degree() == 0) return {};
    if (degree() == 1) {
        // a*x + b = 0 -> x = -b/a
        const auto& a = coeffs_[1];
        const auto& b = coeffs_[0];
        return {mul(S::NegativeOne(), b) / a};
    }
    if (degree() == 2) {
        // a*x^2 + b*x + c = 0 -> x = (-b ± sqrt(b^2 - 4ac)) / (2a)
        const auto& a = coeffs_[2];
        const auto& b = coeffs_[1];
        const auto& c = coeffs_[0];
        auto disc = b * b - integer(4) * a * c;
        auto sqrt_disc = sqrt(disc);
        auto two_a = integer(2) * a;
        auto neg_b = mul(S::NegativeOne(), b);
        return {(neg_b + sqrt_disc) / two_a, (neg_b - sqrt_disc) / two_a};
    }
    if (degree() == 3) {
        // Try rational roots first to keep output clean when possible.
        auto rat = rational_roots(*this);
        if (!rat.empty()) {
            // Deflate by each rational root; recurse on the cofactor.
            Poly remaining = *this;
            for (const auto& r : rat) {
                Poly factor(std::vector<Expr>{mul(S::NegativeOne(), r), S::One()},
                            var_);
                remaining = remaining / factor;
            }
            std::vector<Expr> all = rat;
            for (auto& r : remaining.roots()) all.push_back(r);
            return all;
        }
        return cubic_roots(coeffs_[3], coeffs_[2], coeffs_[1], coeffs_[0]);
    }
    if (degree() == 4) {
        auto rat = rational_roots(*this);
        if (!rat.empty()) {
            Poly remaining = *this;
            for (const auto& r : rat) {
                Poly factor(std::vector<Expr>{mul(S::NegativeOne(), r), S::One()},
                            var_);
                remaining = remaining / factor;
            }
            std::vector<Expr> all = rat;
            for (auto& r : remaining.roots()) all.push_back(r);
            return all;
        }
        // No rational root, but the quartic may still factor into two quadratics
        // over ℚ — prefer that to Ferrari's nested radicals.
        if (auto split = roots_by_factoring(*this)) return *split;
        return quartic_roots(coeffs_[4], coeffs_[3], coeffs_[2], coeffs_[1],
                             coeffs_[0]);
    }
    // Degree ≥ 5: by Abel-Ruffini there's no general radical formula.
    // Pull out rational roots and recurse on the cofactor — the cofactor
    // may itself be degree ≤ 4 and fully solvable.
    auto rat = rational_roots(*this);
    if (rat.empty()) {
        // No rational root: try a non-trivial factorization over ℚ before
        // giving up. A degree-≥5 polynomial with no radical formula may still
        // split into solvable (≤4) factors, e.g. (x²+x+1)(x³+2).
        if (auto split = roots_by_factoring(*this)) return *split;
        return {};
    }
    Poly remaining = *this;
    for (const auto& r : rat) {
        Poly factor(std::vector<Expr>{mul(S::NegativeOne(), r), S::One()},
                    var_);
        remaining = remaining / factor;
    }
    std::vector<Expr> all = rat;
    for (auto& r : remaining.roots()) all.push_back(r);
    return all;
}

}  // namespace sympp
