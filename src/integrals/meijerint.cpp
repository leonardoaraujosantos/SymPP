#include <sympp/integrals/meijerint.hpp>

#include <cstddef>
#include <vector>

#include <sympp/core/basic.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/combinatorial.hpp>  // gamma
#include <sympp/functions/hypergeometric.hpp>  // meijerg
#include <sympp/simplify/simplify.hpp>

namespace sympp {

namespace {

[[nodiscard]] bool is_meijerg(const Expr& e) {
    return e->type_id() == TypeId::Function &&
           static_cast<const Function&>(*e).function_id() == FunctionId::MeijerG;
}

[[nodiscard]] long parse_count(const Expr& e) {
    if (e->type_id() != TypeId::Integer) return -1;
    return static_cast<const Integer&>(*e).value().get_si();
}

// Parsed Meijer-G: the four parameter groups, plus the scale η and validity.
struct MGData {
    std::vector<Expr> an, ap_rest, bm, bq_rest;
    Expr eta;       // argument = η·var
    bool ok = false;
};

[[nodiscard]] MGData parse_meijerg(const Expr& G, const Expr& var) {
    MGData d;
    if (!is_meijerg(G)) return d;
    const auto& a = G->args();
    long n = parse_count(a[0]), pr = parse_count(a[1]);
    long m = parse_count(a[2]), qr = parse_count(a[3]);
    if (n < 0 || pr < 0 || m < 0 || qr < 0) return d;

    auto slice = [&](std::size_t start, long count) {
        return std::vector<Expr>(a.begin() + static_cast<std::ptrdiff_t>(start),
                                 a.begin() + static_cast<std::ptrdiff_t>(start) + count);
    };
    std::size_t off = 4;
    d.an = slice(off, n);            off += static_cast<std::size_t>(n);
    d.ap_rest = slice(off, pr);      off += static_cast<std::size_t>(pr);
    d.bm = slice(off, m);            off += static_cast<std::size_t>(m);
    d.bq_rest = slice(off, qr);
    Expr z = a.back();

    // Argument must be η·var with η free of var.
    Expr eta = simplify(mul(z, pow(var, integer(-1))));
    if (has(eta, var)) return d;
    d.eta = eta;
    d.ok = true;
    return d;
}

}  // namespace

std::optional<Expr> meijerg_mellin_transform(const Expr& G, const Expr& var, const Expr& s) {
    MGData d = parse_meijerg(G, var);
    if (!d.ok) return std::nullopt;
    Expr one = S::One();
    Expr num = one, den = one;
    for (const auto& b : d.bm) num = mul(num, gamma(b + s));            // Πᵐ Γ(b_j+s)
    for (const auto& aa : d.an) num = mul(num, gamma(one - aa - s));    // Πⁿ Γ(1−a_j−s)
    for (const auto& b : d.bq_rest) den = mul(den, gamma(one - b - s)); // Π_{m+1}^q Γ(1−b_j−s)
    for (const auto& aa : d.ap_rest) den = mul(den, gamma(aa + s));     // Π_{n+1}^p Γ(a_j+s)
    Expr ratio = mul(num, pow(den, integer(-1)));
    Expr scale = pow(d.eta, mul(integer(-1), s));                       // η^{−s}
    return simplify(mul(scale, ratio));
}

std::optional<Expr> meijerg_integrate_0_inf(const Expr& G, const Expr& var) {
    // Convergence at the +∞ end needs the exponential-decay scale η > 0. If the
    // scale's sign is unknown (e.g. ∫₀^∞ e^{−a·x} with generic a), abstain so the
    // caller can keep the integral unevaluated rather than return a value that is
    // only valid for η > 0. (The 0-end divergence is caught by the Gamma-pole
    // guard below.)
    MGData d = parse_meijerg(G, var);
    if (!d.ok || is_positive(d.eta) != std::optional<bool>{true}) return std::nullopt;

    auto mt = meijerg_mellin_transform(G, var, S::One());  // transform at s = 1
    if (!mt) return std::nullopt;
    Expr v = simplify(*mt);
    // Reject non-convergent results (a Gamma pole leaves a complex-infinity / nan).
    if (v->type_id() == TypeId::ComplexInfinity || v->type_id() == TypeId::NaN ||
        v->type_id() == TypeId::Infinity || v->type_id() == TypeId::NegativeInfinity) {
        return std::nullopt;
    }
    return v;
}

namespace {
[[nodiscard]] bool is_exp(const Expr& e) {
    return e->type_id() == TypeId::Function &&
           static_cast<const Function&>(*e).function_id() == FunctionId::Exp;
}
// Match var^a with `a` free of var; on success set a_out and return true.
[[nodiscard]] bool match_power_of(const Expr& e, const Expr& var, Expr& a_out) {
    if (e == var) {
        a_out = S::One();
        return true;
    }
    if (e->type_id() == TypeId::Pow && e->args()[0] == var && !has(e->args()[1], var)) {
        a_out = e->args()[1];
        return true;
    }
    return false;
}
}  // namespace

std::optional<Expr> to_meijerg(const Expr& f, const Expr& var) {
    if (!has(f, var)) return std::nullopt;

    // e^{g}  →  G^{1,0}_{0,1}([],[],[0],[], −g).
    if (is_exp(f)) {
        Expr g = f->args()[0];
        return meijerg({}, {}, {integer(0)}, {}, mul(integer(-1), g));
    }

    // 1/(1+var)  →  G^{1,1}_{1,1}([0],[],[0],[], var).
    if (f->type_id() == TypeId::Pow && f->args()[1] == integer(-1)) {
        Expr base = simplify(f->args()[0] - var);  // base == 1 + var ?
        if (base == S::One()) {
            return meijerg({integer(0)}, {}, {integer(0)}, {}, var);
        }
    }

    // var^a · e^{−var}  →  G^{1,0}_{0,1}([],[],[a],[], var).
    if (f->type_id() == TypeId::Mul) {
        Expr a_exp;
        bool have_pow = false, have_exp = false;
        for (const auto& factor : f->args()) {
            Expr a;
            if (!have_pow && match_power_of(factor, var, a)) {
                a_exp = a;
                have_pow = true;
            } else if (!have_exp && is_exp(factor) &&
                       simplify(factor->args()[0] + var) == S::Zero()) {
                have_exp = true;  // exp(−var)
            } else {
                have_pow = have_exp = false;  // an unexpected factor → bail
                break;
            }
        }
        if (have_pow && have_exp) {
            return meijerg({}, {}, {a_exp}, {}, var);
        }
    }

    return std::nullopt;
}

std::optional<Expr> meijerg_integrate_0_inf_of(const Expr& f, const Expr& var) {
    auto g = to_meijerg(f, var);
    if (!g) return std::nullopt;
    return meijerg_integrate_0_inf(*g, var);
}

}  // namespace sympp
