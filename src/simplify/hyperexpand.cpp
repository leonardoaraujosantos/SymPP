#include <sympp/simplify/hyperexpand.hpp>

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/hypergeometric.hpp>
#include <sympp/functions/bessel.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/special.hpp>
#include <sympp/functions/trigonometric.hpp>

namespace sympp {

namespace {

[[nodiscard]] bool is_hyper(const Expr& e) {
    if (e->type_id() != TypeId::Function) return false;
    const auto& fn = static_cast<const Function&>(*e);
    return fn.function_id() == FunctionId::Hyper;
}

[[nodiscard]] bool is_meijerg(const Expr& e) {
    if (e->type_id() != TypeId::Function) return false;
    return static_cast<const Function&>(*e).function_id() == FunctionId::MeijerG;
}

[[nodiscard]] long parse_count(const Expr& e) {
    if (e->type_id() != TypeId::Integer) return -1;
    return static_cast<const Integer&>(*e).value().get_si();
}

// Reduce a Meijer-G node to a sum of hypergeometric functions via Slater's
// theorem, in the generic lower-parameter case (no two of b₁…b_m differ by an
// integer). Returns std::nullopt — leaving the node opaque — for the confluent
// case or a malformed/empty-lower-set node.
//
//   G^{m,n}_{p,q}(z) = Σ_{k=1}^{m} A_k·z^{b_k}·pF_{q−1}(1+b_k−a ; 1+b_k−b' ; σz)
//   A_k = Π_{j≠k}Γ(b_j−b_k)·Π_{j≤n}Γ(1+b_k−a_j)
//         / [ Π_{j>m}Γ(1+b_k−b_j)·Π_{j>n}Γ(a_j−b_k) ],   σ = (−1)^{p−m−n}
[[nodiscard]] std::optional<Expr> meijerg_to_hyper(const Expr& e) {
    if (!is_meijerg(e)) return std::nullopt;
    const auto& a = e->args();
    long n = parse_count(a[0]);    // |an|   (numerator, first group)
    long pr = parse_count(a[1]);   // |ap_rest|
    long m = parse_count(a[2]);    // |bm|   (denominator poles)
    long qr = parse_count(a[3]);   // |bq_rest|
    if (n < 0 || pr < 0 || m < 1 || qr < 0) return std::nullopt;  // m≥1 (else empty sum)
    long p = n + pr, q = m + qr;

    auto slice = [&](std::size_t start, long count) {
        return std::vector<Expr>(a.begin() + static_cast<std::ptrdiff_t>(start),
                                 a.begin() + static_cast<std::ptrdiff_t>(start) + count);
    };
    std::size_t off = 4;
    std::vector<Expr> an = slice(off, n);          off += static_cast<std::size_t>(n);
    std::vector<Expr> ap_rest = slice(off, pr);    off += static_cast<std::size_t>(pr);
    std::vector<Expr> bm = slice(off, m);          off += static_cast<std::size_t>(m);
    std::vector<Expr> bq_rest = slice(off, qr);
    Expr z = a.back();

    std::vector<Expr> A = an;                       // a₁…aₚ
    A.insert(A.end(), ap_rest.begin(), ap_rest.end());
    std::vector<Expr> B = bm;                       // b₁…b_q
    B.insert(B.end(), bq_rest.begin(), bq_rest.end());

    // Confluent G^{2,0}_{0,2}(z|;;b₁,b₂) with b₁−b₂ ∈ ℤ — the generic Slater form
    // breaks (repeated poles), but the closed form is a modified Bessel K:
    //   G^{2,0}_{0,2}(z) = 2·z^{(b₁+b₂)/2}·K_{b₁−b₂}(2√z).
    if (n == 0 && p == 0 && m == 2 && q == 2 &&
        (bm[0] - bm[1])->type_id() == TypeId::Integer) {
        Expr order = bm[0] - bm[1];
        Expr expo = mul(rational(1, 2), add(bm[0], bm[1]));
        return mul(integer(2), mul(pow(z, expo), besselk(order, mul(integer(2), sqrt(z)))));
    }

    // Generic guard: no two lower poles b₁…b_m may differ by an integer.
    for (std::size_t i = 0; i < bm.size(); ++i)
        for (std::size_t j = i + 1; j < bm.size(); ++j)
            if ((bm[i] - bm[j])->type_id() == TypeId::Integer) return std::nullopt;

    Expr one = S::One();
    bool neg = (((p - m - n) % 2) != 0);           // σ = (−1)^{p−m−n}
    Expr w = neg ? mul(integer(-1), z) : z;

    Expr result = S::Zero();
    for (long k = 0; k < m; ++k) {
        const Expr& bk = B[static_cast<std::size_t>(k)];
        Expr num = one, den = one;
        for (long j = 0; j < m; ++j)
            if (j != k) num = mul(num, gamma(B[static_cast<std::size_t>(j)] - bk));
        for (long j = 0; j < n; ++j)
            num = mul(num, gamma(one + bk - A[static_cast<std::size_t>(j)]));
        for (long j = m; j < q; ++j)
            den = mul(den, gamma(one + bk - B[static_cast<std::size_t>(j)]));
        for (long j = n; j < p; ++j)
            den = mul(den, gamma(A[static_cast<std::size_t>(j)] - bk));
        Expr coeff = mul(num, pow(den, integer(-1)));

        std::vector<Expr> top, bottom;
        for (long j = 0; j < p; ++j) top.push_back(one + bk - A[static_cast<std::size_t>(j)]);
        for (long j = 0; j < q; ++j)
            if (j != k) bottom.push_back(one + bk - B[static_cast<std::size_t>(j)]);

        result = add(result, mul(coeff, mul(pow(z, bk), hyper(top, bottom, w))));
    }
    return result;
}

[[nodiscard]] bool both_match(const std::vector<Expr>& v,
                                  const std::vector<Expr>& target) {
    if (v.size() != target.size()) return false;
    // Order-independent match — Pochhammer products commute, so SymPy
    // also normalizes the parameter sets before checking.
    std::vector<bool> used(target.size(), false);
    for (const auto& a : v) {
        bool matched = false;
        for (std::size_t j = 0; j < target.size(); ++j) {
            if (used[j]) continue;
            if (a == target[j]) {
                used[j] = true;
                matched = true;
                break;
            }
        }
        if (!matched) return false;
    }
    return true;
}

// Rewrite a hyper(ap, bq, z) node using the recognized-identity table.
// Returns std::nullopt if no rule fires; the caller leaves the node alone.
[[nodiscard]] std::optional<Expr> rewrite_hyper_node(const Expr& e) {
    if (!is_hyper(e)) return std::nullopt;
    const auto& h = static_cast<const Hyper&>(*e);
    auto a = h.ap();
    auto b = h.bq();
    auto z = h.z();
    auto p = a.size();
    auto q = b.size();

    Expr half = rational(1, 2);
    Expr threehalf = rational(3, 2);
    Expr one = S::One();
    Expr two = integer(2);

    Expr sqz = sqrt(z);  // √z, used by the radical-argument closed forms below

    // ₀F₁ family (p == 0, q == 1): ₀F₁(; b; z) = Σ zᵏ/((b)ₖ k!).
    if (p == 0 && q == 1) {
        // ₀F₁(; 1/2; z) = cosh(2√z).
        if (b[0] == half) return cosh(mul(two, sqz));
        // ₀F₁(; 3/2; z) = sinh(2√z)/(2√z).
        if (b[0] == threehalf) return sinh(mul(two, sqz)) / mul(two, sqz);
    }

    // ₁F₁(a; b; z) family (p == 1, q == 1).
    if (p == 1 && q == 1) {
        // ₁F₁(1; 2; z) = (eᶻ − 1) / z.
        if (a[0] == one && b[0] == two) {
            return (exp(z) - one) / z;
        }
        // ₁F₁(1; 3/2; z) = √π·eᶻ·erf(√z) / (2√z).
        if (a[0] == one && b[0] == threehalf) {
            return mul(sqrt(S::Pi()), mul(exp(z), erf(sqz))) / mul(two, sqz);
        }
    }

    // ₂F₁ family (p == 2, q == 1).
    if (p == 2 && q == 1) {
        // ₂F₁(1, 1; 2; z) = −log(1 − z)/z.
        if (both_match(a, {one, one}) && b[0] == two) {
            return -log(one - z) / z;
        }
        // ₂F₁(1/2, 1; 3/2; z) = atanh(√z)/√z   (→ arctan(y)/y when z = −y²).
        if (both_match(a, {half, one}) && b[0] == threehalf) {
            return atanh(sqz) / sqz;
        }
        // ₂F₁(1/2, 1/2; 3/2; z) = asin(√z)/√z  (→ arcsinh(y)/y when z = −y²).
        if (both_match(a, {half, half}) && b[0] == threehalf) {
            return asin(sqz) / sqz;
        }
    }

    return std::nullopt;
}

[[nodiscard]] Expr apply_recursive(const Expr& e) {
    if (!e) return e;
    // Meijer-G: reduce to a hypergeometric sum first, then expand that.
    if (auto mg = meijerg_to_hyper(e)) return apply_recursive(*mg);
    auto args = e->args();
    if (args.empty()) {
        if (auto r = rewrite_hyper_node(e); r) return *r;
        return e;
    }
    std::vector<Expr> new_args;
    new_args.reserve(args.size());
    for (const auto& a : args) new_args.push_back(apply_recursive(a));
    Expr rebuilt;
    switch (e->type_id()) {
        case TypeId::Add: rebuilt = add(std::move(new_args)); break;
        case TypeId::Mul: rebuilt = mul(std::move(new_args)); break;
        case TypeId::Pow: rebuilt = pow(new_args[0], new_args[1]); break;
        case TypeId::Function:
            // Rebuild the function with rewritten args via its rebuild()
            // hook (which goes back through the auto-eval factories).
            rebuilt = static_cast<const Function&>(*e).rebuild(
                std::move(new_args));
            break;
        default:
            rebuilt = e;
            break;
    }
    if (auto r = rewrite_hyper_node(rebuilt); r) return *r;
    return rebuilt;
}

}  // namespace

Expr hyperexpand(const Expr& e) { return apply_recursive(e); }

}  // namespace sympp
