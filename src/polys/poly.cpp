#include <sympp/polys/poly.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/miscellaneous.hpp>

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
    // Degree >= 3: closed-form roots are non-trivial; defer.
    return {};
}

}  // namespace sympp
