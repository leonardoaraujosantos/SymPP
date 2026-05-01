#include <sympp/integrals/quadrature.hpp>

#include <stdexcept>
#include <string>
#include <utility>

#include <mpfr.h>

#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// Materialize a numeric expression to a double for sample-point evaluation.
// We use mpfr internally to avoid double's precision issues, but the
// driver below works in MPFR throughout.
void to_mpfr(const Expr& e, mpfr_t out, int dps) {
    Expr v = evalf(e, dps);
    if (v->type_id() != TypeId::Float) {
        // Last resort: try a string round-trip via the str() form.
        // For non-numeric inputs the caller will catch the error.
        try {
            mpfr_set_str(out, v->str().c_str(), 10, MPFR_RNDN);
            return;
        } catch (...) {
            mpfr_set_d(out, 0.0, MPFR_RNDN);
            return;
        }
    }
    // Float-typed Expr — round-trip via str() to get the value into MPFR.
    mpfr_set_str(out, v->str().c_str(), 10, MPFR_RNDN);
}

// Sample f(at) at high precision.
void sample_at(const Expr& f, const Expr& var, mpfr_srcptr at_val,
               int dps, mpfr_t out) {
    // Build an Expr for `at_val` and substitute.
    char buf[256];
    mpfr_snprintf(buf, sizeof(buf), "%.*Re", dps + 5, at_val);
    auto at_expr = make<Float>(std::string(buf), dps + 5);
    Expr substituted = subs(f, var, at_expr);
    to_mpfr(substituted, out, dps + 5);
}

// Simpson's 1/3 rule on [a, b] using already-known fa, fm, fb.
void simpson(mpfr_srcptr a, mpfr_srcptr b,
             mpfr_srcptr fa, mpfr_srcptr fm, mpfr_srcptr fb,
             mpfr_t out, mpfr_prec_t prec) {
    mpfr_t scratch, h;
    mpfr_init2(scratch, prec);
    mpfr_init2(h, prec);
    mpfr_sub(h, b, a, MPFR_RNDN);            // h = b - a
    // out = h/6 * (fa + 4*fm + fb)
    mpfr_mul_ui(scratch, fm, 4, MPFR_RNDN);  // scratch = 4*fm
    mpfr_add(scratch, scratch, fa, MPFR_RNDN);
    mpfr_add(scratch, scratch, fb, MPFR_RNDN);
    mpfr_mul(out, h, scratch, MPFR_RNDN);
    mpfr_div_ui(out, out, 6, MPFR_RNDN);
    mpfr_clear(scratch);
    mpfr_clear(h);
}

// Recursive adaptive Simpson. Returns the result via `out`.
void adaptive_simpson(const Expr& f, const Expr& var, int dps,
                      mpfr_prec_t prec,
                      mpfr_srcptr a, mpfr_srcptr b,
                      mpfr_srcptr fa, mpfr_srcptr fm, mpfr_srcptr fb,
                      mpfr_srcptr whole,
                      mpfr_srcptr eps, int depth, mpfr_t out) {
    mpfr_t mid, lm, rm, flm, frm, left, right, sum, diff, m, abs_diff;
    mpfr_inits2(prec, mid, lm, rm, flm, frm, left, right, sum, diff, m,
                abs_diff, (mpfr_ptr)nullptr);

    mpfr_add(m, a, b, MPFR_RNDN);
    mpfr_div_2ui(m, m, 1, MPFR_RNDN);                      // m = (a+b)/2

    mpfr_add(lm, a, m, MPFR_RNDN);
    mpfr_div_2ui(lm, lm, 1, MPFR_RNDN);                    // lm = (a+m)/2
    mpfr_add(rm, m, b, MPFR_RNDN);
    mpfr_div_2ui(rm, rm, 1, MPFR_RNDN);                    // rm = (m+b)/2

    sample_at(f, var, lm, dps, flm);
    sample_at(f, var, rm, dps, frm);

    simpson(a, m, fa, flm, fm, left, prec);
    simpson(m, b, fm, frm, fb, right, prec);
    mpfr_add(sum, left, right, MPFR_RNDN);

    mpfr_sub(diff, sum, whole, MPFR_RNDN);
    mpfr_abs(abs_diff, diff, MPFR_RNDN);

    // 15*eps is the standard adaptive-Simpson tolerance (3rd-order term).
    mpfr_t fifteen_eps;
    mpfr_init2(fifteen_eps, prec);
    mpfr_mul_ui(fifteen_eps, eps, 15, MPFR_RNDN);

    if (depth <= 0 || mpfr_cmp(abs_diff, fifteen_eps) <= 0) {
        // Accept: sum + (sum - whole)/15 (Richardson correction).
        mpfr_div_ui(out, diff, 15, MPFR_RNDN);
        mpfr_add(out, sum, out, MPFR_RNDN);
    } else {
        mpfr_t l_out, r_out, half_eps;
        mpfr_inits2(prec, l_out, r_out, half_eps, (mpfr_ptr)nullptr);
        mpfr_div_2ui(half_eps, eps, 1, MPFR_RNDN);
        adaptive_simpson(f, var, dps, prec, a, m, fa, flm, fm, left,
                         half_eps, depth - 1, l_out);
        adaptive_simpson(f, var, dps, prec, m, b, fm, frm, fb, right,
                         half_eps, depth - 1, r_out);
        mpfr_add(out, l_out, r_out, MPFR_RNDN);
        mpfr_clears(l_out, r_out, half_eps, (mpfr_ptr)nullptr);
    }

    mpfr_clear(fifteen_eps);
    mpfr_clears(mid, lm, rm, flm, frm, left, right, sum, diff, m, abs_diff,
                (mpfr_ptr)nullptr);
}

}  // namespace

Expr vpaintegral(const Expr& f, const Expr& var, const Expr& lower,
                 const Expr& upper, int dps) {
    if (dps < 6) dps = 6;
    const mpfr_prec_t prec = dps_to_prec(dps + 10);

    mpfr_t a, b, fa, fm, fb, whole, m, eps, out;
    mpfr_inits2(prec, a, b, fa, fm, fb, whole, m, eps, out, (mpfr_ptr)nullptr);

    to_mpfr(lower, a, dps + 10);
    to_mpfr(upper, b, dps + 10);
    sample_at(f, var, a, dps, fa);
    sample_at(f, var, b, dps, fb);
    mpfr_add(m, a, b, MPFR_RNDN);
    mpfr_div_2ui(m, m, 1, MPFR_RNDN);
    sample_at(f, var, m, dps, fm);
    simpson(a, b, fa, fm, fb, whole, prec);

    // tolerance: 10^(-dps).
    mpfr_set_d(eps, 10.0, MPFR_RNDN);
    mpfr_pow_si(eps, eps, -dps, MPFR_RNDN);

    adaptive_simpson(f, var, dps, prec, a, b, fa, fm, fb, whole, eps,
                     30, out);

    Expr result = make<Float>(static_cast<mpfr_srcptr>(out), dps);
    mpfr_clears(a, b, fa, fm, fb, whole, m, eps, out, (mpfr_ptr)nullptr);
    return result;
}

}  // namespace sympp
