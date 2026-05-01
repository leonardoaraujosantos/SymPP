#include <sympp/polys/rootof.hpp>

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/polys/poly.hpp>

namespace sympp {

namespace {

// Hash combining helper (matching style used elsewhere in core).
[[nodiscard]] std::size_t hash_combine(std::size_t a, std::size_t b) noexcept {
    return a ^ (b + 0x9e3779b9 + (a << 6) + (a >> 2));
}

// Pull rational coefficients out of `poly_expr` interpreted as a polynomial
// in `var`. Returns nullopt if any coefficient isn't a numeric Integer or
// Rational.
[[nodiscard]] std::optional<std::vector<mpq_class>>
extract_rational_coeffs(const Expr& poly_expr, const Expr& var) {
    Poly p(poly_expr, var);
    std::vector<mpq_class> out;
    out.reserve(p.coeffs().size());
    for (const auto& c : p.coeffs()) {
        if (!c) return std::nullopt;
        if (c->type_id() == TypeId::Integer) {
            out.push_back(mpq_class(static_cast<const Integer&>(*c).value()));
        } else if (c->type_id() == TypeId::Rational) {
            out.push_back(static_cast<const Rational&>(*c).value());
        } else {
            return std::nullopt;
        }
    }
    return out;
}

}  // namespace

RootOf::RootOf(Expr poly_expr, Expr var, std::size_t index)
    : args_{std::move(poly_expr), std::move(var)}, index_(index) {
    hash_ = hash_combine(args_[0]->hash(), args_[1]->hash());
    hash_ = hash_combine(hash_, index_);
    hash_ = hash_combine(hash_, 0x1007'0F'1D'EA'FA47ULL);
}

bool RootOf::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::RootOf) return false;
    const auto& o = static_cast<const RootOf&>(other);
    return index_ == o.index_ && args_[0] == o.args_[0]
           && args_[1] == o.args_[1];
}

std::string RootOf::str() const {
    // SymPy printable form: CRootOf(<poly>, <index>) — variable inferred.
    return "CRootOf(" + args_[0]->str() + ", " + std::to_string(index_) + ")";
}

std::optional<bool> RootOf::ask(AssumptionKey k) const noexcept {
    // Without numerically resolving, we can't say much. Real roots are real,
    // but we don't know which roots are real here. Be conservative.
    switch (k) {
        case AssumptionKey::Finite:
            return true;
        default:
            return std::nullopt;
    }
}

std::optional<Expr> RootOf::try_evalf(int dps) const {
    auto qcoeffs = extract_rational_coeffs(args_[0], args_[1]);
    if (!qcoeffs || qcoeffs->size() < 2) return std::nullopt;

    const mpfr_prec_t prec = dps_to_prec(dps);
    const std::size_t n = qcoeffs->size();

    // Cauchy bound on root magnitude: 1 + max|a_i|/|a_n|.
    double bound = 1.0;
    {
        double max_abs = 0.0;
        for (std::size_t i = 0; i + 1 < n; ++i) {
            double v = std::fabs((*qcoeffs)[i].get_d());
            if (v > max_abs) max_abs = v;
        }
        double leading = std::fabs(qcoeffs->back().get_d());
        if (leading > 0) bound = 1.0 + max_abs / leading;
        if (bound < 2.0) bound = 2.0;  // tiny safety margin
    }

    // Step 1: locate sign changes via dense double-precision sampling.
    auto eval_double = [&](double x) {
        double y = 0.0, power = 1.0;
        for (auto& c : *qcoeffs) {
            y += c.get_d() * power;
            power *= x;
        }
        return y;
    };
    const int samples = static_cast<int>(200 * n + 200);
    std::vector<std::pair<double, double>> brackets;
    double prev_x = -bound;
    double prev_y = eval_double(prev_x);
    for (int i = 1; i <= samples; ++i) {
        double x = -bound + (2.0 * bound) * static_cast<double>(i)
                   / static_cast<double>(samples);
        double y = eval_double(x);
        if ((prev_y < 0 && y > 0) || (prev_y > 0 && y < 0)) {
            brackets.emplace_back(prev_x, x);
        } else if (y == 0) {
            // Exact zero; record a tiny bracket around x.
            brackets.emplace_back(x - 1e-9, x + 1e-9);
        }
        prev_x = x; prev_y = y;
    }

    if (index_ >= brackets.size()) return std::nullopt;
    auto [a0, b0] = brackets[index_];

    // Step 2: bisect in MPFR for high precision.
    mpfr_t a, b, mid, fa, fmid, power, term;
    mpfr_init2(a, prec); mpfr_init2(b, prec); mpfr_init2(mid, prec);
    mpfr_init2(fa, prec); mpfr_init2(fmid, prec);
    mpfr_init2(power, prec); mpfr_init2(term, prec);
    mpfr_set_d(a, a0, MPFR_RNDN);
    mpfr_set_d(b, b0, MPFR_RNDN);

    auto eval_mpfr = [&](mpfr_t out, mpfr_srcptr at) {
        mpfr_set_d(out, 0.0, MPFR_RNDN);
        mpfr_set_d(power, 1.0, MPFR_RNDN);
        for (auto& c : *qcoeffs) {
            mpfr_set_q(term, c.get_mpq_t(), MPFR_RNDN);
            mpfr_mul(term, term, power, MPFR_RNDN);
            mpfr_add(out, out, term, MPFR_RNDN);
            mpfr_mul(power, power, at, MPFR_RNDN);
        }
    };

    eval_mpfr(fa, a);
    const int iters = 4 * dps + 100;
    for (int it = 0; it < iters; ++it) {
        mpfr_add(mid, a, b, MPFR_RNDN);
        mpfr_div_2ui(mid, mid, 1, MPFR_RNDN);
        eval_mpfr(fmid, mid);
        if (mpfr_zero_p(fmid)) break;
        if (mpfr_signbit(fmid) == mpfr_signbit(fa)) {
            mpfr_set(a, mid, MPFR_RNDN);
            mpfr_set(fa, fmid, MPFR_RNDN);
        } else {
            mpfr_set(b, mid, MPFR_RNDN);
        }
    }
    Expr result = make<Float>(static_cast<mpfr_srcptr>(mid), dps);
    mpfr_clear(a); mpfr_clear(b); mpfr_clear(mid);
    mpfr_clear(fa); mpfr_clear(fmid);
    mpfr_clear(power); mpfr_clear(term);
    return result;
}

Expr root_of(const Expr& poly_expr, const Expr& var, std::size_t index) {
    return make<RootOf>(poly_expr, var, index);
}

}  // namespace sympp
