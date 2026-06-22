#include <sympp/integrals/meijerint.hpp>

#include <cstddef>
#include <vector>

#include <sympp/core/basic.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/combinatorial.hpp>  // gamma
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

}  // namespace sympp
