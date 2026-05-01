#include <sympp/polys/poly.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/imaginary_unit.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/miscellaneous.hpp>

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
        const std::size_t shift = r.degree() - other.degree();
        Expr factor = r.leading_coeff() / lc_other;
        q[shift] = factor;
        std::vector<Expr> sub(shift + other.coeffs_.size(), S::Zero());
        for (std::size_t i = 0; i < other.coeffs_.size(); ++i) {
            sub[i + shift] = mul(factor, other.coeffs_[i]);
        }
        Poly subp{std::move(sub), var_};
        r = r - subp;
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
        return quartic_roots(coeffs_[4], coeffs_[3], coeffs_[2], coeffs_[1],
                             coeffs_[0]);
    }
    // Degree ≥ 5: by Abel-Ruffini there's no general radical formula.
    // Pull out rational roots and recurse on the cofactor — the cofactor
    // may itself be degree ≤ 4 and fully solvable.
    auto rat = rational_roots(*this);
    if (rat.empty()) return {};
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
