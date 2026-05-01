#include <sympp/integrals/transforms.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/imaginary_unit.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/special.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

namespace {

// Try to extract a constant `a` such that `g == a*t`. Returns std::nullopt
// when g doesn't have that form.
[[nodiscard]] std::optional<Expr> extract_linear_coeff(const Expr& g, const Expr& t) {
    if (g == t) return S::One();
    if (g->type_id() == TypeId::Mul) {
        std::vector<Expr> rest;
        bool found_t = false;
        for (const auto& a : g->args()) {
            if (a == t) {
                if (found_t) return std::nullopt;
                found_t = true;
            } else if (has(a, t)) {
                return std::nullopt;
            } else {
                rest.push_back(a);
            }
        }
        if (!found_t) return std::nullopt;
        return mul(std::move(rest));
    }
    return std::nullopt;
}

// Table-lookup on a single term. Returns std::nullopt if no entry matches.
[[nodiscard]] std::optional<Expr> laplace_term(const Expr& f, const Expr& t, const Expr& s) {
    // Constant in t: ∫ c * exp(-s*t) dt from 0 to ∞ = c / s
    if (!has(f, t)) {
        return mul(f, pow(s, S::NegativeOne()));
    }

    // f == t: 1/s^2
    if (f == t) {
        return pow(s, integer(-2));
    }

    // f == t^n: n! / s^(n+1)
    if (f->type_id() == TypeId::Pow) {
        const auto& base = f->args()[0];
        const auto& exp = f->args()[1];
        if (base == t && exp->type_id() == TypeId::Integer) {
            const auto& n = static_cast<const Integer&>(*exp);
            if (!n.is_negative() && n.fits_long() && n.to_long() <= 100) {
                long ne = n.to_long();
                Expr fact = factorial(integer(ne));
                Expr sn1 = pow(s, integer(ne + 1));
                return mul(fact, pow(sn1, S::NegativeOne()));
            }
        }
    }

    // exp / sin / cos with linear inner.
    if (f->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.args().size() == 1) {
            const auto& inner = fn.args()[0];
            auto coeff = extract_linear_coeff(inner, t);
            if (coeff.has_value()) {
                switch (fn.function_id()) {
                    case FunctionId::Exp:
                        // L{exp(a*t)} = 1/(s - a)
                        return pow(s - *coeff, S::NegativeOne());
                    case FunctionId::Sin: {
                        // L{sin(a*t)} = a/(s^2 + a^2)
                        auto denom = pow(s, integer(2)) + pow(*coeff, integer(2));
                        return mul(*coeff, pow(denom, S::NegativeOne()));
                    }
                    case FunctionId::Cos: {
                        // L{cos(a*t)} = s/(s^2 + a^2)
                        auto denom = pow(s, integer(2)) + pow(*coeff, integer(2));
                        return mul(s, pow(denom, S::NegativeOne()));
                    }
                    default:
                        break;
                }
            }
        }
    }

    // Mul with a constant factor: pull it out and recurse.
    if (f->type_id() == TypeId::Mul) {
        std::vector<Expr> consts;
        std::vector<Expr> rest;
        for (const auto& a : f->args()) {
            if (has(a, t)) rest.push_back(a);
            else consts.push_back(a);
        }
        if (!consts.empty() && !rest.empty()) {
            Expr inner = mul(std::move(rest));
            if (auto sub = laplace_term(inner, t, s); sub.has_value()) {
                return mul(mul(std::move(consts)), *sub);
            }
        }
    }

    return std::nullopt;
}

}  // namespace

Expr laplace_transform(const Expr& f, const Expr& t, const Expr& s) {
    if (!f) return S::Zero();
    // Linearity over Add.
    if (f->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(f->args().size());
        for (const auto& a : f->args()) {
            terms.push_back(laplace_transform(a, t, s));
        }
        return add(std::move(terms));
    }
    if (auto r = laplace_term(f, t, s); r.has_value()) {
        return *r;
    }
    return function_symbol("LaplaceTransform")(f, t, s);
}

// ---------------------------------------------------------------------------
// Inverse Laplace
// ---------------------------------------------------------------------------

namespace {

// Try to detect a Pow(s - a, k) factor in `den`. Returns (a, k) when found;
// requires a to be free of `s` and k to be a positive integer.
[[nodiscard]] std::optional<std::pair<Expr, long>>
match_linear_pow(const Expr& den, const Expr& s) {
    Expr base;
    long exponent = 1;
    if (den->type_id() == TypeId::Pow) {
        const auto& exp = den->args()[1];
        if (exp->type_id() != TypeId::Integer) return std::nullopt;
        const auto& z = static_cast<const Integer&>(*exp);
        if (z.is_negative() || !z.fits_long()) return std::nullopt;
        exponent = z.to_long();
        base = den->args()[0];
    } else {
        base = den;
    }
    // base should be (s - a) or (s + a) — i.e., affine with leading 1 in s.
    // Decompose into (coeff, const).
    if (base == s) {
        return std::pair<Expr, long>{S::Zero(), exponent};
    }
    // Try base = s + b: Add of s + something.
    if (base->type_id() == TypeId::Add) {
        std::vector<Expr> rest;
        bool found_s = false;
        for (const auto& term : base->args()) {
            if (term == s) {
                if (found_s) return std::nullopt;
                found_s = true;
            } else if (has(term, s)) {
                return std::nullopt;
            } else {
                rest.push_back(term);
            }
        }
        if (!found_s) return std::nullopt;
        Expr b = add(std::move(rest));
        // base = s + b, so a = -b in (s - a).
        return std::pair<Expr, long>{mul(S::NegativeOne(), b), exponent};
    }
    return std::nullopt;
}

// Single-summand inverse Laplace lookup. Returns nullopt when the
// summand isn't recognized.
[[nodiscard]] std::optional<Expr> inverse_laplace_term(
    const Expr& F, const Expr& s, const Expr& t) {
    // Constant in s: that's nominally the transform of c·δ(t).
    if (!has(F, s)) {
        return mul(F, dirac_delta(t));
    }
    // F == 1/s^n form: pull numerator (constant w.r.t. s) and denominator.
    Expr num_const = S::One();
    Expr den;
    if (F->type_id() == TypeId::Mul) {
        std::vector<Expr> num_factors;
        for (const auto& f : F->args()) {
            if (f->type_id() == TypeId::Pow
                && f->args()[1]->type_id() == TypeId::Integer
                && static_cast<const Integer&>(*f->args()[1]).is_negative()) {
                if (den) return std::nullopt;
                long e = static_cast<const Integer&>(*f->args()[1]).to_long();
                den = pow(f->args()[0], integer(-e));
                continue;
            }
            num_factors.push_back(f);
        }
        if (!den) return std::nullopt;
        num_const = mul(num_factors);
    } else if (F->type_id() == TypeId::Pow
               && F->args()[1]->type_id() == TypeId::Integer
               && static_cast<const Integer&>(*F->args()[1]).is_negative()) {
        long e = static_cast<const Integer&>(*F->args()[1]).to_long();
        den = pow(F->args()[0], integer(-e));
        num_const = S::One();
    } else {
        return std::nullopt;
    }
    if (has(num_const, s)) return std::nullopt;

    // Case 1: den = (s - a)^n.
    if (auto m = match_linear_pow(den, s); m) {
        const Expr& a = m->first;
        long n = m->second;
        if (n == 1) {
            // 1/(s - a) → exp(a*t)
            return mul(num_const, exp(mul(a, t)));
        }
        // 1/(s - a)^n → t^(n-1) * exp(a*t) / (n-1)!
        Expr tn = pow(t, integer(n - 1));
        Expr fact = factorial(integer(n - 1));
        return mul({num_const, tn, exp(mul(a, t)),
                    pow(fact, S::NegativeOne())});
    }

    // Case 2: den = s^2 + a^2 (irreducible quadratic over ℝ).
    // Or numerator = s on top of the same denominator.
    // Detect den = s² + c with c free of s and positive.
    if (den->type_id() == TypeId::Add) {
        Expr s_squared_coeff;
        Expr constant_part;
        bool ok = true;
        for (const auto& term : den->args()) {
            if (term->type_id() == TypeId::Pow
                && term->args()[0] == s
                && term->args()[1] == integer(2)) {
                if (s_squared_coeff) { ok = false; break; }
                s_squared_coeff = S::One();
                continue;
            }
            // Coefficient * s^2.
            if (term->type_id() == TypeId::Mul) {
                bool has_ssq = false;
                std::vector<Expr> coef;
                for (const auto& f : term->args()) {
                    if (!has_ssq && f->type_id() == TypeId::Pow
                        && f->args()[0] == s && f->args()[1] == integer(2)) {
                        has_ssq = true;
                    } else {
                        coef.push_back(f);
                    }
                }
                if (has_ssq) {
                    if (s_squared_coeff) { ok = false; break; }
                    s_squared_coeff = mul(coef);
                    continue;
                }
            }
            if (has(term, s)) { ok = false; break; }
            constant_part = constant_part ? add(constant_part, term) : term;
        }
        if (ok && s_squared_coeff && constant_part
            && s_squared_coeff == S::One()) {
            // den == s^2 + constant_part. constant_part should be a^2.
            Expr a_sq = constant_part;
            // If the numerator is num_const (constant in s):
            // num_const / (s^2 + a^2) → num_const * sin(a*t) / a
            //   where a = sqrt(a_sq).
            Expr a = sqrt(a_sq);
            return mul({num_const, sin(mul(a, t)), pow(a, S::NegativeOne())});
        }
    }
    return std::nullopt;
}

// Detect a Mul factor of the form s * (something_inverse) so we can
// match the s / (s² + a²) → cos(a·t) pair.
[[nodiscard]] std::optional<Expr>
inverse_laplace_s_over_quad(const Expr& F, const Expr& s, const Expr& t) {
    if (F->type_id() != TypeId::Mul) return std::nullopt;
    bool has_s_factor = false;
    Expr den_factor_inverse;  // a Pow(_, -1)
    std::vector<Expr> rest;
    for (const auto& f : F->args()) {
        if (!has_s_factor && f == s) { has_s_factor = true; continue; }
        if (!den_factor_inverse
            && f->type_id() == TypeId::Pow
            && f->args()[1] == S::NegativeOne()) {
            den_factor_inverse = f->args()[0];
            continue;
        }
        rest.push_back(f);
    }
    if (!has_s_factor || !den_factor_inverse) return std::nullopt;
    // den should be s² + a²: Add of s² and a constant.
    if (den_factor_inverse->type_id() != TypeId::Add) return std::nullopt;
    Expr s_sq, constant_part;
    for (const auto& term : den_factor_inverse->args()) {
        if (term->type_id() == TypeId::Pow
            && term->args()[0] == s && term->args()[1] == integer(2)) {
            if (s_sq) return std::nullopt;
            s_sq = term;
            continue;
        }
        if (has(term, s)) return std::nullopt;
        constant_part = constant_part ? add(constant_part, term) : term;
    }
    if (!s_sq || !constant_part) return std::nullopt;
    Expr a = sqrt(constant_part);
    Expr coef = mul(rest);
    return mul(coef, cos(mul(a, t)));
}

}  // namespace

Expr inverse_laplace_transform(const Expr& F, const Expr& s, const Expr& t) {
    if (!F) return S::Zero();
    // Linearity.
    if (F->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(F->args().size());
        for (const auto& a : F->args()) {
            terms.push_back(inverse_laplace_transform(a, s, t));
        }
        return add(std::move(terms));
    }
    // s/(s²+a²) — handle before generic numerator-only case.
    if (auto r = inverse_laplace_s_over_quad(F, s, t); r.has_value()) {
        return *r;
    }
    if (auto r = inverse_laplace_term(F, s, t); r.has_value()) {
        return *r;
    }
    // Try apart() to split into table-recognizable pieces.
    Expr decomposed = apart(F, s);
    if (!(decomposed == F)) {
        return inverse_laplace_transform(decomposed, s, t);
    }
    return function_symbol("InverseLaplaceTransform")(F, s, t);
}

// ---------------------------------------------------------------------------
// Fourier
// ---------------------------------------------------------------------------

namespace {

[[nodiscard]] std::optional<Expr> fourier_term(
    const Expr& f, const Expr& t, const Expr& w) {
    // δ(t) → 1
    if (f->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.function_id() == FunctionId::DiracDelta && fn.args().size() == 1) {
            if (fn.args()[0] == t) return S::One();
        }
    }
    // Constant → 2π · δ(ω)
    if (!has(f, t)) {
        return mul({integer(2), S::Pi(), f, dirac_delta(w)});
    }
    // exp(I·a·t) → 2π·δ(ω - a)
    if (f->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1) {
            const Expr& inner = fn.args()[0];
            // inner = c·t. Pull out the coefficient; it should contain I.
            auto coeff = extract_linear_coeff(inner, t);
            if (coeff.has_value()) {
                // exp(c·t): treat c = I·a, return 2π·δ(ω - a) where a = c/I = -I·c.
                Expr a = mul(mul(S::NegativeOne(), S::I()), *coeff);
                return mul({integer(2), S::Pi(), dirac_delta(w - a)});
            }
        }
    }
    // exp(-a·t²) → sqrt(π/a) · exp(-ω²/(4a))
    if (f->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1) {
            const Expr& inner = fn.args()[0];
            // inner should be -a·t² with a positive (we don't enforce sign).
            // Match -a · t².
            if (inner->type_id() == TypeId::Mul) {
                std::vector<Expr> coef;
                bool has_t2 = false;
                for (const auto& factor : inner->args()) {
                    if (!has_t2 && factor->type_id() == TypeId::Pow
                        && factor->args()[0] == t
                        && factor->args()[1] == integer(2)) {
                        has_t2 = true;
                    } else if (has(factor, t)) {
                        coef.clear();
                        break;
                    } else {
                        coef.push_back(factor);
                    }
                }
                if (has_t2 && !coef.empty()) {
                    Expr neg_a = mul(coef);  // expect this to be -a
                    Expr a = mul(S::NegativeOne(), neg_a);
                    Expr norm = sqrt(S::Pi() / a);
                    Expr exp_arg = mul(S::NegativeOne(),
                                        pow(w, integer(2)) / (integer(4) * a));
                    return mul(norm, exp(exp_arg));
                }
            }
        }
    }
    // Mul: pull constants out.
    if (f->type_id() == TypeId::Mul) {
        std::vector<Expr> consts, rest;
        for (const auto& a : f->args()) {
            (has(a, t) ? rest : consts).push_back(a);
        }
        if (!consts.empty() && !rest.empty()) {
            Expr inner = mul(rest);
            if (auto sub = fourier_term(inner, t, w); sub.has_value()) {
                return mul(mul(consts), *sub);
            }
        }
    }
    return std::nullopt;
}

}  // namespace

Expr fourier_transform(const Expr& f, const Expr& t, const Expr& w) {
    if (!f) return S::Zero();
    if (f->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(f->args().size());
        for (const auto& a : f->args()) {
            terms.push_back(fourier_transform(a, t, w));
        }
        return add(std::move(terms));
    }
    if (auto r = fourier_term(f, t, w); r.has_value()) return *r;
    return function_symbol("FourierTransform")(f, t, w);
}

// Inverse Fourier: F(ω) → (1/2π) · Fourier(F)(t) with sign flip on ω.
// The cleanest implementation reuses fourier_transform with t and ω swapped
// (and the result divided by 2π) — that holds for the conventions chosen.
Expr inverse_fourier_transform(const Expr& F, const Expr& w, const Expr& t) {
    // Use the symmetry: (1/2π) · Fourier(F)(-t).
    Expr negt = mul(S::NegativeOne(), t);
    Expr forward = fourier_transform(F, w, negt);
    return mul(forward, pow(integer(2) * S::Pi(), S::NegativeOne()));
}

// ---------------------------------------------------------------------------
// Mellin
// ---------------------------------------------------------------------------

namespace {

[[nodiscard]] std::optional<Expr> mellin_term(
    const Expr& f, const Expr& t, const Expr& s) {
    // exp(-t) → Γ(s)
    if (f->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1) {
            const Expr& inner = fn.args()[0];
            // inner should be -c·t, c positive. Extract c.
            if (inner->type_id() == TypeId::Mul) {
                std::vector<Expr> coef;
                bool has_t = false;
                for (const auto& factor : inner->args()) {
                    if (!has_t && factor == t) {
                        has_t = true;
                    } else if (has(factor, t)) {
                        coef.clear();
                        break;
                    } else {
                        coef.push_back(factor);
                    }
                }
                if (has_t && !coef.empty()) {
                    Expr neg_c = mul(coef);  // expect -c
                    Expr c = mul(S::NegativeOne(), neg_c);
                    // M[exp(-c·t)] = Γ(s) / c^s
                    return mul(gamma(s), pow(c, mul(S::NegativeOne(), s)));
                }
            }
            if (inner == mul(S::NegativeOne(), t)) {
                return gamma(s);
            }
        }
    }
    // 1/(1+t) → π/sin(πs)
    if (f->type_id() == TypeId::Pow
        && f->args()[1] == S::NegativeOne()) {
        Expr base = f->args()[0];
        if (base == add(integer(1), t) || base == add(t, integer(1))) {
            return S::Pi() / sin(mul(S::Pi(), s));
        }
    }
    // Mul: pull constants out.
    if (f->type_id() == TypeId::Mul) {
        std::vector<Expr> consts, rest;
        for (const auto& a : f->args()) {
            (has(a, t) ? rest : consts).push_back(a);
        }
        if (!consts.empty() && !rest.empty()) {
            Expr inner = mul(rest);
            if (auto sub = mellin_term(inner, t, s); sub.has_value()) {
                return mul(mul(consts), *sub);
            }
        }
    }
    return std::nullopt;
}

}  // namespace

Expr mellin_transform(const Expr& f, const Expr& t, const Expr& s) {
    if (!f) return S::Zero();
    if (f->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(f->args().size());
        for (const auto& a : f->args()) {
            terms.push_back(mellin_transform(a, t, s));
        }
        return add(std::move(terms));
    }
    if (auto r = mellin_term(f, t, s); r.has_value()) return *r;
    return function_symbol("MellinTransform")(f, t, s);
}

Expr inverse_mellin_transform(const Expr& F, const Expr& s, const Expr& t) {
    // Reverse the table for the two pairs we know.
    // Γ(s) → exp(-t).
    if (F->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*F);
        if (fn.function_id() == FunctionId::Gamma && fn.args().size() == 1
            && fn.args()[0] == s) {
            return exp(mul(S::NegativeOne(), t));
        }
    }
    return function_symbol("InverseMellinTransform")(F, s, t);
}

// ---------------------------------------------------------------------------
// Z-transform
// ---------------------------------------------------------------------------

namespace {

// Try to match a^n where a is constant (free of n).
[[nodiscard]] std::optional<Expr> match_geometric(const Expr& f, const Expr& n) {
    if (f->type_id() != TypeId::Pow) return std::nullopt;
    const Expr& base = f->args()[0];
    const Expr& exp = f->args()[1];
    if (exp != n) return std::nullopt;
    if (has(base, n)) return std::nullopt;
    return base;
}

[[nodiscard]] std::optional<Expr> z_term(
    const Expr& f, const Expr& n, const Expr& z) {
    // 1 → z/(z - 1)
    if (!has(f, n)) {
        return mul(f, mul(z, pow(z - integer(1), S::NegativeOne())));
    }
    // n → z/(z-1)^2
    if (f == n) {
        return mul(z, pow(pow(z - integer(1), integer(2)), S::NegativeOne()));
    }
    // a^n → z/(z - a)
    if (auto a = match_geometric(f, n); a) {
        return mul(z, pow(z - *a, S::NegativeOne()));
    }
    // n·a^n → a·z/(z - a)^2
    if (f->type_id() == TypeId::Mul) {
        bool has_n_factor = false;
        Expr a_pow_n;
        std::vector<Expr> rest;
        for (const auto& factor : f->args()) {
            if (!has_n_factor && factor == n) {
                has_n_factor = true;
                continue;
            }
            if (!a_pow_n) {
                if (auto a = match_geometric(factor, n); a) {
                    a_pow_n = factor;
                    continue;
                }
            }
            rest.push_back(factor);
        }
        if (has_n_factor && a_pow_n) {
            const Expr& a = a_pow_n->args()[0];
            Expr result = mul(a, z) / pow(z - a, integer(2));
            if (!rest.empty()) result = mul(mul(rest), result);
            return result;
        }
        // Pure constant-factor pull-out.
        std::vector<Expr> consts, var_factors;
        for (const auto& a : f->args()) {
            (has(a, n) ? var_factors : consts).push_back(a);
        }
        if (!consts.empty() && !var_factors.empty()) {
            Expr inner = mul(var_factors);
            if (auto sub = z_term(inner, n, z); sub.has_value()) {
                return mul(mul(consts), *sub);
            }
        }
    }
    return std::nullopt;
}

}  // namespace

Expr z_transform(const Expr& f, const Expr& n, const Expr& z) {
    if (!f) return S::Zero();
    if (f->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(f->args().size());
        for (const auto& a : f->args()) {
            terms.push_back(z_transform(a, n, z));
        }
        return add(std::move(terms));
    }
    if (auto r = z_term(f, n, z); r.has_value()) return *r;
    return function_symbol("ZTransform")(f, n, z);
}

Expr inverse_z_transform(const Expr& F, const Expr& z, const Expr& n) {
    // Detect z/(z-a) → a^n and z/(z-a)^2 → n·a^(n-1) by inspection.
    if (F->type_id() == TypeId::Mul) {
        bool has_z_num = false;
        Expr inner_inv;
        long power = 0;
        std::vector<Expr> rest;
        for (const auto& factor : F->args()) {
            if (!has_z_num && factor == z) { has_z_num = true; continue; }
            if (!inner_inv && factor->type_id() == TypeId::Pow) {
                const Expr& exp = factor->args()[1];
                if (exp->type_id() == TypeId::Integer) {
                    long e = static_cast<const Integer&>(*exp).to_long();
                    if (e < 0) {
                        inner_inv = factor->args()[0];
                        power = -e;
                        continue;
                    }
                }
            }
            rest.push_back(factor);
        }
        if (has_z_num && inner_inv) {
            // Determine a from inner_inv, possibly a Pow base.
            Expr base = inner_inv;
            long base_pow = 1;
            if (base->type_id() == TypeId::Pow
                && base->args()[1]->type_id() == TypeId::Integer) {
                long e = static_cast<const Integer&>(*base->args()[1]).to_long();
                if (e > 0) {
                    base_pow = e;
                    base = base->args()[0];
                }
            }
            // base should be (z - a). Extract a.
            if (base->type_id() == TypeId::Add) {
                Expr a_neg;
                bool has_z_inside = false;
                for (const auto& term : base->args()) {
                    if (term == z) { has_z_inside = true; continue; }
                    if (has(term, z)) { has_z_inside = false; break; }
                    a_neg = a_neg ? add(a_neg, term) : term;
                }
                if (has_z_inside) {
                    Expr a = mul(S::NegativeOne(), a_neg ? a_neg : S::Zero());
                    long total = power == 1 ? base_pow : (power * base_pow);
                    if (total == 1) {
                        // z/(z-a) → a^n
                        Expr result = pow(a, n);
                        if (!rest.empty()) result = mul(mul(rest), result);
                        return result;
                    }
                    if (total == 2) {
                        // z/(z-a)^2 → n·a^(n-1)
                        Expr result = mul(n, pow(a, n - integer(1)));
                        if (!rest.empty()) result = mul(mul(rest), result);
                        return result;
                    }
                }
            }
        }
    }
    return function_symbol("InverseZTransform")(F, z, n);
}

// ---------------------------------------------------------------------------
// Sine and cosine transforms (half-line Fourier)
// ---------------------------------------------------------------------------

namespace {

// Match exp(-a·t) with a constant in t. Returns a or nullopt.
[[nodiscard]] std::optional<Expr> match_exp_neg(
    const Expr& f, const Expr& t) {
    if (f->type_id() != TypeId::Function) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*f);
    if (fn.function_id() != FunctionId::Exp || fn.args().size() != 1) {
        return std::nullopt;
    }
    const Expr& inner = fn.args()[0];
    // inner should be -a·t. Extract -a.
    auto coeff = extract_linear_coeff(inner, t);
    if (!coeff) return std::nullopt;
    return mul(S::NegativeOne(), *coeff);
}

}  // namespace

Expr sine_transform(const Expr& f, const Expr& t, const Expr& w) {
    if (f->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(f->args().size());
        for (const auto& a : f->args()) {
            terms.push_back(sine_transform(a, t, w));
        }
        return add(std::move(terms));
    }
    if (f->type_id() == TypeId::Mul) {
        std::vector<Expr> consts, rest;
        for (const auto& a : f->args()) {
            (has(a, t) ? rest : consts).push_back(a);
        }
        if (!consts.empty() && !rest.empty()) {
            Expr inner_t = mul(rest);
            return mul(mul(consts), sine_transform(inner_t, t, w));
        }
    }
    // exp(-a·t) → ω/(a² + ω²)
    if (auto a_opt = match_exp_neg(f, t); a_opt) {
        const Expr& a = *a_opt;
        return w / (pow(a, integer(2)) + pow(w, integer(2)));
    }
    return function_symbol("SineTransform")(f, t, w);
}

Expr cosine_transform(const Expr& f, const Expr& t, const Expr& w) {
    if (f->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(f->args().size());
        for (const auto& a : f->args()) {
            terms.push_back(cosine_transform(a, t, w));
        }
        return add(std::move(terms));
    }
    if (f->type_id() == TypeId::Mul) {
        std::vector<Expr> consts, rest;
        for (const auto& a : f->args()) {
            (has(a, t) ? rest : consts).push_back(a);
        }
        if (!consts.empty() && !rest.empty()) {
            Expr inner_t = mul(rest);
            return mul(mul(consts), cosine_transform(inner_t, t, w));
        }
    }
    // exp(-a·t) → a/(a² + ω²)
    if (auto a_opt = match_exp_neg(f, t); a_opt) {
        const Expr& a = *a_opt;
        return a / (pow(a, integer(2)) + pow(w, integer(2)));
    }
    return function_symbol("CosineTransform")(f, t, w);
}

}  // namespace sympp
