#include <sympp/functions/miscellaneous.hpp>

#include <algorithm>
#include <functional>
#include <optional>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/expr_collections.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/imaginary_unit.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/trigonometric.hpp>

namespace sympp {

namespace {

// True if the expression tree contains an unevaluated Re or Im node — used to
// decide whether a complex value's real/imaginary split fully resolved.
[[nodiscard]] bool contains_re_im(const Expr& e) {
    if (e->type_id() == TypeId::Function) {
        auto id = static_cast<const Function&>(*e).function_id();
        if (id == FunctionId::Re || id == FunctionId::Im) return true;
    }
    for (const auto& a : e->args()) {
        if (contains_re_im(a)) return true;
    }
    return false;
}

// Detect a leading minus sign on a Mul: -1 * rest. Returns the stripped tail
// if so, std::nullopt otherwise. Number negatives stay as-is for the abs path
// (we handle them separately).
[[nodiscard]] std::optional<Expr> strip_neg_factor(const Expr& e) {
    if (e->type_id() == TypeId::Mul) {
        const auto& args = e->args();
        if (!args.empty() && args[0]->type_id() == TypeId::Integer) {
            const auto& n = static_cast<const Integer&>(*args[0]);
            if (n.value() == -1) {
                if (args.size() == 2) return args[1];
                std::vector<Expr> rest(args.begin() + 1, args.end());
                return mul(std::move(rest));
            }
        }
    }
    return std::nullopt;
}

// Negate a Number — return its absolute value.
[[nodiscard]] Expr abs_number(const Number& n) {
    switch (n.type_id()) {
        case TypeId::Integer: {
            mpz_class v = static_cast<const Integer&>(n).value();
            if (sgn(v) < 0) v = -v;
            return make<Integer>(std::move(v));
        }
        case TypeId::Rational: {
            mpq_class q = static_cast<const Rational&>(n).value();
            if (sgn(q) < 0) q = -q;
            q.canonicalize();
            if (q.get_den() == 1) return make<Integer>(q.get_num());
            return make<Rational>(std::move(q));
        }
        case TypeId::Float: {
            const auto& f = static_cast<const Float&>(n);
            int dps = f.precision_dps();
            mpfr_t r;
            mpfr_init2(r, dps_to_prec(dps));
            mpfr_abs(r, f.value(), MPFR_RNDN);
            auto out = make<Float>(static_cast<mpfr_srcptr>(r), dps);
            mpfr_clear(r);
            return out;
        }
        default:
            // Should not happen in practice — Number subclasses are exhaustive.
            return Expr{};
    }
}

// A rational (Integer or Rational) literal.
[[nodiscard]] bool is_rational_literal(const Expr& e) {
    return e->type_id() == TypeId::Integer || e->type_id() == TypeId::Rational;
}

// If `m` is a pure-imaginary term with a rational coefficient (b·I, or I
// itself), return b; std::nullopt otherwise.
[[nodiscard]] std::optional<Expr> rational_imag_coeff(const Expr& m) {
    if (m == S::I()) return S::One();
    if (m->type_id() != TypeId::Mul) return std::nullopt;
    bool has_i = false;
    Expr coeff = S::One();
    for (const auto& f : m->args()) {
        if (f == S::I()) {
            if (has_i) return std::nullopt;  // I² etc.
            has_i = true;
        } else if (is_rational_literal(f)) {
            coeff = mul(coeff, f);
        } else {
            return std::nullopt;  // symbolic / irrational factor
        }
    }
    if (!has_i) return std::nullopt;
    return coeff;
}

// If `e` is a complex number a + b·I with rational real and imaginary parts,
// return (a, b); std::nullopt if any term is symbolic or irrational.
[[nodiscard]] std::optional<std::pair<Expr, Expr>> rational_complex(
    const Expr& e) {
    if (e == S::I()) return std::make_pair(S::Zero(), S::One());
    if (auto b = rational_imag_coeff(e)) {
        return std::make_pair(S::Zero(), *b);
    }
    if (e->type_id() == TypeId::Add) {
        Expr re = S::Zero();
        Expr im = S::Zero();
        for (const auto& t : e->args()) {
            if (is_rational_literal(t)) {
                re = add(re, t);
            } else if (auto b = rational_imag_coeff(t)) {
                im = add(im, *b);
            } else {
                return std::nullopt;
            }
        }
        return std::make_pair(re, im);
    }
    return std::nullopt;
}

}  // namespace

// ----- Abs -------------------------------------------------------------------

Abs::Abs(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Abs);
}

Expr Abs::rebuild(std::vector<Expr> new_args) const {
    return abs(new_args[0]);
}

// d/darg |arg| = sign(arg). diff() applies the chain rule, multiplying this
// by darg/dvar. (For real arg; matches SymPy's Abs(x).diff(x) == sign(x).)
Expr Abs::diff_arg(std::size_t /*i*/) const {
    return sign(args_[0]);
}

std::optional<bool> Abs::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
            return true;  // |x| is always real (for any complex x, |x| ∈ ℝ_≥0)
        case AssumptionKey::Negative:
            return false;
        case AssumptionKey::Nonnegative:
            return true;
        case AssumptionKey::Positive:
            if (is_nonzero(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Zero:
            if (is_zero(a) == true) return true;
            if (is_nonzero(a) == true) return false;
            return std::nullopt;
        case AssumptionKey::Integer:
            if (is_integer(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Rational:
            if (is_rational(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Finite:
            if (is_finite(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr abs(const Expr& arg) {
    // Numeric: take the absolute value directly.
    if (is_number(arg)) {
        return abs_number(static_cast<const Number&>(*arg));
    }

    // NumberSymbol — Pi, E, EulerGamma, Catalan are all positive.
    if (arg->type_id() == TypeId::NumberSymbol) {
        return arg;
    }

    // |exp(z)| = exp(re(z)). With re() already evaluating the imaginary part, this
    // gives the unit modulus |exp(I·x)| = 1 for real x (re(I·x) = 0), the general
    // |exp(x)| = exp(re(x)), and |exp(I·x)| = exp(−im(x)) for a complex x. SymPy
    // does the same.
    if (arg->type_id() == TypeId::Function
        && static_cast<const Function&>(*arg).function_id() == FunctionId::Exp) {
        return exp(re(arg->args()[0]));
    }

    // Complex number a + b·I with rational parts → sqrt(a² + b²) (the modulus).
    // Covers Abs(I) = 1, Abs(b·I) = |b|, Abs(3 + 4·I) = 5, etc.
    if (auto z = rational_complex(arg); z.has_value() && z->second != S::Zero()) {
        const Expr& re = z->first;
        const Expr& im = z->second;
        return sqrt(add(mul(re, re), mul(im, im)));
    }

    // Abs(x) for nonnegative x → x.
    if (is_nonnegative(arg) == true) return arg;
    // Abs(x) for nonpositive x → -x.
    if (is_nonpositive(arg) == true) return mul(S::NegativeOne(), arg);

    // |·| is multiplicative: Abs(∏fᵢ) = ∏|fᵢ|. Pull out every factor whose
    // modulus is known — numeric coefficients (|c|), positive factors (= f), and
    // negative factors (= −f) — leaving the rest under a single Abs. Covers
    // Abs(c·rest)=|c|·Abs(rest), Abs(-x)=Abs(x), and Abs(p·x)=p·Abs(x) for a
    // positive symbol p.
    if (arg->type_id() == TypeId::Mul) {
        std::vector<Expr> pulled;
        std::vector<Expr> kept;
        for (const auto& f : arg->args()) {
            // I has modulus 1, so |I·x| = |x|.
            if (is_number(f) || is_positive(f) == true
                || is_negative(f) == true || f == S::I()) {
                pulled.push_back(abs(f));
            } else {
                kept.push_back(f);
            }
        }
        if (!pulled.empty() && !kept.empty()) {
            return mul(mul(std::move(pulled)),
                       make<Abs>(mul(std::move(kept))));
        }
    }

    // Symbolic complex modulus: |a + b·I| = √(a² + b²) when re/im resolve to
    // expressions free of unevaluated Re/Im (e.g. real-symbol parts) and the
    // imaginary part is nonzero. A generic Abs(z) keeps its Re(z)/Im(z) split
    // and so is left unevaluated, matching SymPy.
    {
        Expr a = re(arg);
        Expr b = im(arg);
        if (!(b == S::Zero()) && !contains_re_im(a) && !contains_re_im(b)) {
            return sqrt(add(mul(a, a), mul(b, b)));
        }
    }

    return make<Abs>(arg);
}

// ----- Sign ------------------------------------------------------------------

Sign::Sign(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Sign);
}

Expr Sign::rebuild(std::vector<Expr> new_args) const {
    return sign(new_args[0]);
}

std::optional<bool> Sign::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Integer:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Finite:
            return true;
        default:
            return std::nullopt;
    }
}

Expr sign(const Expr& arg) {
    // Numeric: use the sign predicate.
    if (is_number(arg)) {
        int s = static_cast<const Number&>(*arg).sign();
        if (s > 0) return S::One();
        if (s < 0) return S::NegativeOne();
        return S::Zero();
    }

    // NumberSymbol → +1 (Pi, E, etc.).
    if (arg->type_id() == TypeId::NumberSymbol) {
        return S::One();
    }

    // Idempotent: sign(sign(z)) = sign(z). sign(z) is already in {−1, 0, 1} for
    // real z (and a unit complex for complex z), whose sign is itself.
    if (arg->type_id() == TypeId::Function
        && static_cast<const Function&>(*arg).function_id()
               == FunctionId::Sign) {
        return arg;
    }

    if (is_zero(arg) == true) return S::Zero();
    if (is_positive(arg) == true) return S::One();
    if (is_negative(arg) == true) return S::NegativeOne();

    // sign(-x) = -sign(x)
    if (auto pos = strip_neg_factor(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Sign>(*pos));
    }

    return make<Sign>(arg);
}

// ----- sqrt — uses Pow under the hood -----------------------------------------

// ============================================================================
// Complex part extraction
// ============================================================================

Re::Re(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Re);
}
Expr Re::rebuild(std::vector<Expr> new_args) const { return re(new_args[0]); }
std::optional<bool> Re::ask(AssumptionKey k) const noexcept {
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
            return true;  // re(z) is always real for any complex z
        case AssumptionKey::Finite:
            if (is_finite(args_[0]) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Im::Im(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Im);
}
Expr Im::rebuild(std::vector<Expr> new_args) const { return im(new_args[0]); }
std::optional<bool> Im::ask(AssumptionKey k) const noexcept {
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
            return true;  // im(z) ∈ ℝ
        case AssumptionKey::Finite:
            if (is_finite(args_[0]) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Conjugate::Conjugate(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Conjugate);
}
Expr Conjugate::rebuild(std::vector<Expr> new_args) const {
    return conjugate(new_args[0]);
}
std::optional<bool> Conjugate::ask(AssumptionKey k) const noexcept {
    // Conjugation preserves real-ness, integer-ness, etc. on real arguments;
    // mostly mirrors arg's properties.
    return args_[0]->ask(k);
}

Arg::Arg(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Arg);
}
Expr Arg::rebuild(std::vector<Expr> new_args) const { return arg_(new_args[0]); }
std::optional<bool> Arg::ask(AssumptionKey k) const noexcept {
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
            return true;
        default:
            return std::nullopt;
    }
}

namespace {

// Decompose a Mul as c · Iᵏ · w, with c the product of the real factors,
// k = (number of I factors) mod 4, and w the remaining (non-real, non-I) part.
// Returns nullopt unless something was actually peeled off (k > 0 or c ≠ 1), so
// the re()/im() recursion below terminates.
struct MulComplexParts { Expr c; int k; Expr w; };
[[nodiscard]] std::optional<MulComplexParts> decompose_mul_complex(
    const Expr& arg) {
    if (arg->type_id() != TypeId::Mul) return std::nullopt;
    Expr c = S::One();
    int k = 0;
    std::vector<Expr> rest;
    for (const auto& f : arg->args()) {
        if (f == S::I()) {
            k = (k + 1) % 4;
        } else if (is_real(f) == true) {
            c = mul(c, f);
        } else {
            rest.push_back(f);
        }
    }
    if (k == 0 && c == S::One()) return std::nullopt;  // nothing extracted
    Expr w = rest.empty() ? Expr{S::One()}
                          : (rest.size() == 1 ? rest[0] : mul(std::move(rest)));
    return MulComplexParts{std::move(c), k, std::move(w)};
}

}  // namespace

// Rationalize a complex denominator: 1/(a+bI) = (a−bI)/(a²+b²). Rewrites every
// Pow(d, −m) whose base d carries the imaginary unit and whose |d|² = d·conj(d)
// is provably real (the numeric-complex case), so the expanded result is a real
// + imaginary split. Symbolic denominators (where |d|² stays complex) are left
// untouched.
[[nodiscard]] Expr rationalize_complex(const Expr& e) {
    ExprMap<Expr> repl;
    std::function<void(const Expr&)> scan = [&](const Expr& x) {
        if (x->type_id() == TypeId::Pow
            && x->args()[1]->type_id() == TypeId::Integer) {
            const Expr& base = x->args()[0];
            const auto& z = static_cast<const Integer&>(*x->args()[1]);
            if (z.is_negative() && z.fits_long() && has(base, S::I())
                && is_real(base) != true) {
                Expr cj = conjugate(base);
                Expr denom = expand(mul(base, cj));
                if (is_real(denom) == true) {
                    const long m = -z.to_long();
                    repl.emplace(x, mul(pow(cj, integer(m)),
                                        pow(denom, integer(-m))));
                }
            }
        }
        for (const auto& a : x->args()) scan(a);
    };
    scan(e);
    if (repl.empty()) return e;
    return xreplace(e, repl);
}

Expr re(const Expr& arg) {
    // Multiply out any complex denominator, then re-enter on the a+bI form.
    if (Expr n = expand(rationalize_complex(arg)); !(n == arg)) return re(n);
    // Numeric: real -> identity.
    if (is_real(arg) == true) return arg;
    // Numeric complex a + b·I → a.
    if (auto z = rational_complex(arg); z.has_value()) return z->first;
    // Linearity over a sum: re(Σ aᵢ) = Σ re(aᵢ).
    if (arg->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(arg->args().size());
        for (const auto& a : arg->args()) terms.push_back(re(a));
        return add(std::move(terms));
    }
    // c · Iᵏ · w: re pulls out the real c and rotates by Iᵏ —
    //   re(c·w) = c·re(w),  re(c·I·w) = −c·im(w),  (k mod 4 cycles the sign).
    if (auto d = decompose_mul_complex(arg)) {
        switch (d->k) {
            case 0: return mul(d->c, re(d->w));
            case 1: return mul(mul(S::NegativeOne(), d->c), im(d->w));
            case 2: return mul(mul(S::NegativeOne(), d->c), re(d->w));
            case 3: return mul(d->c, im(d->w));
        }
    }
    // Linearity: re(-x) = -re(x)
    if (auto p = strip_neg_factor(arg); p.has_value()) {
        return mul(S::NegativeOne(), make<Re>(*p));
    }
    return make<Re>(arg);
}

Expr im(const Expr& arg) {
    if (Expr n = expand(rationalize_complex(arg)); !(n == arg)) return im(n);
    if (is_real(arg) == true) return S::Zero();
    // Numeric complex a + b·I → b.
    if (auto z = rational_complex(arg); z.has_value()) return z->second;
    // Linearity over a sum: im(Σ aᵢ) = Σ im(aᵢ).
    if (arg->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(arg->args().size());
        for (const auto& a : arg->args()) terms.push_back(im(a));
        return add(std::move(terms));
    }
    // c · Iᵏ · w: im(c·w) = c·im(w),  im(c·I·w) = c·re(w),  (k mod 4 cycles).
    if (auto d = decompose_mul_complex(arg)) {
        switch (d->k) {
            case 0: return mul(d->c, im(d->w));
            case 1: return mul(d->c, re(d->w));
            case 2: return mul(mul(S::NegativeOne(), d->c), im(d->w));
            case 3: return mul(mul(S::NegativeOne(), d->c), re(d->w));
        }
    }
    if (auto p = strip_neg_factor(arg); p.has_value()) {
        return mul(S::NegativeOne(), make<Im>(*p));
    }
    return make<Im>(arg);
}

Expr conjugate(const Expr& arg) {
    if (is_real(arg) == true) return arg;
    // Numeric complex a + b·I → a - b·I.
    if (auto z = rational_complex(arg); z.has_value()) {
        return add(z->first, mul(mul(S::NegativeOne(), z->second), S::I()));
    }
    if (arg->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*arg);
        // conjugate(conjugate(x)) = x.
        if (fn.function_id() == FunctionId::Conjugate) {
            return arg->args()[0];
        }
        // conjugate(f(g)) = f(conjugate(g)) for an entire function with real
        // Taylor coefficients (exp / sin / cos / tan / sinh / cosh / tanh) —
        // e.g. conjugate(exp(I·x)) = exp(−I·x). log is excluded (branch cut).
        if (fn.args().size() == 1) {
            Expr cg = conjugate(fn.args()[0]);
            switch (fn.function_id()) {
                case FunctionId::Exp:  return exp(cg);
                case FunctionId::Sin:  return sin(cg);
                case FunctionId::Cos:  return cos(cg);
                case FunctionId::Tan:  return tan(cg);
                case FunctionId::Sinh: return sinh(cg);
                case FunctionId::Cosh: return cosh(cg);
                case FunctionId::Tanh: return tanh(cg);
                default: break;
            }
        }
    }
    // conjugate is a ring homomorphism: it distributes over products and sums
    // (conjugate(a·b) = ā·b̄, conjugate(a+b) = ā+b̄) and over an integer power.
    // Recursion reduces each part (conjugate(I) = −I, conjugate(real) = real),
    // so conjugate(I·x) = −I·conjugate(x).
    if (arg->type_id() == TypeId::Mul) {
        std::vector<Expr> factors;
        factors.reserve(arg->args().size());
        for (const auto& f : arg->args()) factors.push_back(conjugate(f));
        return mul(std::move(factors));
    }
    if (arg->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(arg->args().size());
        for (const auto& a : arg->args()) terms.push_back(conjugate(a));
        return add(std::move(terms));
    }
    if (arg->type_id() == TypeId::Pow
        && arg->args()[1]->type_id() == TypeId::Integer) {
        return pow(conjugate(arg->args()[0]), arg->args()[1]);
    }
    return make<Conjugate>(arg);
}

Expr arg_(const Expr& arg) {
    // arg(0) is undefined — the origin has no argument. SymPy returns nan.
    if (arg == S::Zero()) return S::NaN();
    if (is_positive(arg) == true) return S::Zero();
    if (is_negative(arg) == true) return S::Pi();
    // arg(z) = atan2(im z, re z), applied when there is a (resolved) imaginary
    // part: arg(I) = atan2(1, 0) = π/2, arg(1+I) = atan2(1, 1) = π/4. A real
    // value (im = 0) of unknown sign stays unevaluated.
    Expr a = re(arg);
    Expr b = im(arg);
    if (!(b == S::Zero()) && !contains_re_im(a) && !contains_re_im(b)) {
        return atan2(b, a);
    }
    return make<Arg>(arg);
}

// ============================================================================
// Min, Max
// ============================================================================

MinFn::MinFn(std::vector<Expr> args)
    : Function(std::move(args)) {
    compute_hash(FunctionId::Min);
}
Expr MinFn::rebuild(std::vector<Expr> new_args) const { return min(std::move(new_args)); }
std::optional<bool> MinFn::ask(AssumptionKey k) const noexcept {
    // Min preserves real / integer / rational under all-true closure.
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
        case AssumptionKey::Integer:
        case AssumptionKey::Rational:
        case AssumptionKey::Finite:
            for (const auto& a : args_) {
                if (a->ask(k) != true) return std::nullopt;
            }
            return true;
        default:
            return std::nullopt;
    }
}

MaxFn::MaxFn(std::vector<Expr> args)
    : Function(std::move(args)) {
    compute_hash(FunctionId::Max);
}
Expr MaxFn::rebuild(std::vector<Expr> new_args) const { return max(std::move(new_args)); }
std::optional<bool> MaxFn::ask(AssumptionKey k) const noexcept {
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
        case AssumptionKey::Integer:
        case AssumptionKey::Rational:
        case AssumptionKey::Finite:
            for (const auto& a : args_) {
                if (a->ask(k) != true) return std::nullopt;
            }
            return true;
        default:
            return std::nullopt;
    }
}

namespace {

// Compare two Number Exprs. Returns -1, 0, +1. Caller verifies both are Numbers.
[[nodiscard]] int compare_numbers(const Number& a, const Number& b) noexcept {
    // Promote to mpq for the exact tower; Float comparisons compute via mpfr.
    if (a.type_id() == TypeId::Float || b.type_id() == TypeId::Float) {
        // At least one is Float — promote the other to Float at the higher dps.
        int dps = kDefaultDps;
        if (a.type_id() == TypeId::Float) {
            dps = std::max(dps, static_cast<const Float&>(a).precision_dps());
        }
        if (b.type_id() == TypeId::Float) {
            dps = std::max(dps, static_cast<const Float&>(b).precision_dps());
        }
        mpfr_t ax, bx;
        mpfr_init2(ax, dps_to_prec(dps));
        mpfr_init2(bx, dps_to_prec(dps));

        if (a.type_id() == TypeId::Float) {
            mpfr_set(ax, static_cast<const Float&>(a).value(), MPFR_RNDN);
        } else if (a.type_id() == TypeId::Integer) {
            mpfr_set_z(ax, static_cast<const Integer&>(a).value().get_mpz_t(), MPFR_RNDN);
        } else {
            mpfr_set_q(ax, static_cast<const Rational&>(a).value().get_mpq_t(), MPFR_RNDN);
        }
        if (b.type_id() == TypeId::Float) {
            mpfr_set(bx, static_cast<const Float&>(b).value(), MPFR_RNDN);
        } else if (b.type_id() == TypeId::Integer) {
            mpfr_set_z(bx, static_cast<const Integer&>(b).value().get_mpz_t(), MPFR_RNDN);
        } else {
            mpfr_set_q(bx, static_cast<const Rational&>(b).value().get_mpq_t(), MPFR_RNDN);
        }
        int c = mpfr_cmp(ax, bx);
        mpfr_clear(ax);
        mpfr_clear(bx);
        return c < 0 ? -1 : (c > 0 ? 1 : 0);
    }
    // Both exact: promote to mpq.
    mpq_class qa, qb;
    if (a.type_id() == TypeId::Integer) qa = static_cast<const Integer&>(a).value();
    else qa = static_cast<const Rational&>(a).value();
    if (b.type_id() == TypeId::Integer) qb = static_cast<const Integer&>(b).value();
    else qb = static_cast<const Rational&>(b).value();
    int c = cmp(qa, qb);
    return c < 0 ? -1 : (c > 0 ? 1 : 0);
}

// Reduce a Min/Max args list by combining all numeric args into one extreme,
// keeping symbolic args as-is.
template <bool IsMin>
[[nodiscard]] Expr min_max_impl(std::vector<Expr> args) {
    // Flatten a nested same-kind Min/Max: Max(x, Max(y, z)) → Max(x, y, z).
    {
        const FunctionId same = IsMin ? FunctionId::Min : FunctionId::Max;
        std::vector<Expr> flat;
        flat.reserve(args.size());
        for (auto& a : args) {
            if (a->type_id() == TypeId::Function
                && static_cast<const Function&>(*a).function_id() == same) {
                for (const auto& inner : a->args()) flat.push_back(inner);
            } else {
                flat.push_back(std::move(a));
            }
        }
        args = std::move(flat);
    }

    std::vector<Expr> non_numeric;
    Expr extreme;            // running min or max numeric
    bool dropped_identity = false;  // saw the identity ±∞ (kept as fallback)
    for (auto& a : args) {
        // Infinity identities/absorbers: Max(…, +∞) = +∞ and −∞ is its identity;
        // Min mirrors. The identity is dropped but remembered so Max(−∞, −∞) = −∞.
        if (a->type_id() == TypeId::Infinity) {
            if constexpr (!IsMin) return S::Infinity();
            dropped_identity = true;
            continue;
        }
        if (a->type_id() == TypeId::NegativeInfinity) {
            if constexpr (IsMin) return S::NegativeInfinity();
            dropped_identity = true;
            continue;
        }
        if (is_number(a)) {
            if (!extreme) {
                extreme = a;
            } else {
                int c = compare_numbers(static_cast<const Number&>(*extreme),
                                       static_cast<const Number&>(*a));
                if constexpr (IsMin) {
                    if (c > 0) extreme = a;  // a < current min
                } else {
                    if (c < 0) extreme = a;  // a > current max
                }
            }
        } else {
            non_numeric.push_back(std::move(a));
        }
    }

    // Drop duplicate symbolic args — Min(x, x) = x.
    std::sort(non_numeric.begin(), non_numeric.end(),
              [](const Expr& l, const Expr& r) noexcept {
                  return l.get() < r.get();
              });
    non_numeric.erase(std::unique(non_numeric.begin(), non_numeric.end()),
                      non_numeric.end());

    std::vector<Expr> out;
    out.reserve(non_numeric.size() + 1);
    if (extreme) out.push_back(std::move(extreme));
    for (auto& a : non_numeric) out.push_back(std::move(a));

    if (out.empty()) {
        // Everything collapsed to the dropped identity ±∞ (e.g. Max(−∞, −∞)).
        if (dropped_identity) {
            return IsMin ? S::Infinity() : S::NegativeInfinity();
        }
        // Otherwise an empty Min/Max — by convention undefined; SymPy raises.
        // Build the empty form to surface the misuse; rarely hit in practice.
    }
    if (out.size() == 1) return std::move(out[0]);

    if constexpr (IsMin) return make<MinFn>(std::move(out));
    else return make<MaxFn>(std::move(out));
}

}  // namespace

Expr min(std::vector<Expr> args) { return min_max_impl<true>(std::move(args)); }
Expr min(const Expr& a, const Expr& b) {
    return min_max_impl<true>(std::vector<Expr>{a, b});
}
Expr max(std::vector<Expr> args) { return min_max_impl<false>(std::move(args)); }
Expr max(const Expr& a, const Expr& b) {
    return min_max_impl<false>(std::vector<Expr>{a, b});
}

// ============================================================================
// sqrt — uses Pow under the hood
// ============================================================================

Expr sqrt(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    if (arg == S::One()) return S::One();

    // Perfect-square integers fold to integers when they're nonnegative.
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (!z.is_negative()) {
            mpz_class s;
            mpz_class r;
            mpz_sqrtrem(s.get_mpz_t(), r.get_mpz_t(), z.value().get_mpz_t());
            if (sgn(r) == 0) {
                return make<Integer>(std::move(s));
            }
        }
    }

    return pow(arg, S::Half());
}

}  // namespace sympp
