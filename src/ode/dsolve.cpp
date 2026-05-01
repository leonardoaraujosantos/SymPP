#include <sympp/ode/dsolve.hpp>

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/integrals/integrate.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>
#include <sympp/solvers/solve.hpp>

namespace sympp {

namespace {

// Generate a fresh free constant __C{i}.
[[nodiscard]] Expr fresh_constant(int& counter) {
    return symbol("__C" + std::to_string(counter++));
}

// Detect Integral(_, _) failure marker from integrate().
[[nodiscard]] bool is_integral_marker(const Expr& e) {
    if (!e || e->type_id() != TypeId::Function) return false;
    const auto& fn = static_cast<const Function&>(*e);
    return fn.name() == "Integral";
}

// Try to write `eq` as yp - rhs(y, x) = 0 and return rhs. Recognizes
// the canonical Add(yp, -rhs) and Add(-yp, rhs) forms.
[[nodiscard]] std::optional<Expr> isolate_yp(
    const Expr& eq, const Expr& yp) {
    // After simplify, eq should be a polynomial in yp. Build Poly in yp.
    Poly p(simplify(eq), yp);
    if (p.degree() != 1) return std::nullopt;
    // Linear: p[1] * yp + p[0] = 0 → yp = -p[0] / p[1].
    return simplify(mul(S::NegativeOne(), p.coeffs()[0]) / p.coeffs()[1]);
}

// Given rhs = expr in (y, x), check if it factors as f(x) * g(y).
// Returns (f, g) if separable, nullopt otherwise.
[[nodiscard]] std::optional<std::pair<Expr, Expr>>
try_separate(const Expr& rhs, const Expr& y, const Expr& x) {
    // Cases:
    //   * No y: rhs = f(x), trivially separable with g(y) = 1.
    //   * No x: rhs = g(y), trivially separable with f(x) = 1.
    //   * Mul: split factors by which variable they contain.
    if (!has(rhs, y)) return std::pair{rhs, S::One()};
    if (!has(rhs, x)) return std::pair{S::One(), rhs};
    if (rhs->type_id() == TypeId::Mul) {
        std::vector<Expr> fx, gy;
        for (const auto& f : rhs->args()) {
            if (has(f, y) && has(f, x)) return std::nullopt;
            if (has(f, y)) gy.push_back(f);
            else fx.push_back(f);
        }
        return std::pair{mul(fx), mul(gy)};
    }
    return std::nullopt;
}

}  // namespace

// ---------------------------------------------------------------------------
// Separable: y' = f(x) * g(y)  →  ∫ 1/g(y) dy = ∫ f(x) dx + C
// ---------------------------------------------------------------------------

Expr dsolve_separable(const Expr& eq, const Expr& y, const Expr& yp,
                      const Expr& x) {
    auto rhs = isolate_yp(eq, yp);
    if (!rhs) return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
    auto sep = try_separate(*rhs, y, x);
    if (!sep) return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
    Expr f = sep->first;
    Expr g = sep->second;
    Expr lhs_int = integrate(pow(g, S::NegativeOne()), y);
    Expr rhs_int = integrate(f, x);
    if (is_integral_marker(lhs_int) || is_integral_marker(rhs_int)) {
        return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
    }
    int cnt = 0;
    Expr C = fresh_constant(cnt);
    // Explicit-form inversion patterns. lhs_int(y) = rhs_int(x) + C; try
    // to invert lhs_int and return y = lhs_int⁻¹(rhs_int + C).
    //   * lhs_int == log(y)        → y = exp(rhs_int + C)  ·  absorb C: K·exp(rhs_int)
    //   * lhs_int == y             → y = rhs_int + C
    //   * lhs_int == y^k / k       → y = (k(rhs_int + C))^(1/k)
    if (lhs_int->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*lhs_int);
        if (fn.function_id() == FunctionId::Log && fn.args().size() == 1
            && fn.args()[0] == y) {
            return simplify(C * exp(rhs_int));
        }
    }
    if (lhs_int == y) {
        return simplify(rhs_int + C);
    }
    // Try the polynomial-solve path as a fallback.
    Expr implicit = lhs_int - rhs_int - C;
    auto explicit_sols = solve(implicit, y);
    if (!explicit_sols.empty()) {
        return simplify(explicit_sols[0]);
    }
    return implicit;  // F(y) = G(x) + C
}

// ---------------------------------------------------------------------------
// Linear first-order: y' + p(x) y = q(x)
//   Integrating factor μ = exp(∫p dx). Then y = (∫μ·q dx + C) / μ.
// ---------------------------------------------------------------------------

Expr dsolve_linear_first_order(const Expr& eq, const Expr& y, const Expr& yp,
                                const Expr& x) {
    auto rhs = isolate_yp(eq, yp);
    if (!rhs) return function_symbol("Dsolve")(eq, y, x);
    // Express *rhs as -p(x) * y + q(x). Build a polynomial in y.
    Poly poly_in_y(simplify(*rhs), y);
    if (poly_in_y.degree() > 1) return function_symbol("Dsolve")(eq, y, x);
    Expr q = poly_in_y.coeffs()[0];                    // constant w.r.t. y
    Expr neg_p = poly_in_y.degree() == 1
                     ? poly_in_y.coeffs()[1] : S::Zero();  // coeff of y
    if (has(q, y) || has(neg_p, y)) {
        return function_symbol("Dsolve")(eq, y, x);
    }
    Expr p = mul(S::NegativeOne(), neg_p);             // y' = -p y + q  →  y' + p y = q
    Expr int_p = integrate(p, x);
    if (is_integral_marker(int_p)) {
        return function_symbol("Dsolve")(eq, y, x);
    }
    Expr mu = exp(int_p);
    Expr int_mu_q = integrate(mu * q, x);
    if (is_integral_marker(int_mu_q)) {
        return function_symbol("Dsolve")(eq, y, x);
    }
    int cnt = 0;
    Expr C = fresh_constant(cnt);
    Expr sol = (int_mu_q + C) / mu;
    return simplify(sol);
}

// ---------------------------------------------------------------------------
// Bernoulli: y' + p(x) y = q(x) y^n  with n ≠ 0, 1.
//   Substitute v = y^(1-n). Then v' + (1-n) p(x) v = (1-n) q(x).
//   Solve as linear in v, then back-substitute y = v^(1/(1-n)).
// ---------------------------------------------------------------------------

Expr dsolve_bernoulli(const Expr& eq, const Expr& y, const Expr& yp,
                       const Expr& x) {
    auto rhs = isolate_yp(eq, yp);
    if (!rhs) return function_symbol("Dsolve")(eq, y, x);
    // We want rhs = -p(x) y + q(x) y^n. Try to detect by inspecting
    // each Add term for its power of y.
    Expr e = simplify(*rhs);
    if (e->type_id() != TypeId::Add) return function_symbol("Dsolve")(eq, y, x);

    Expr lin_coef = S::Zero();          // coefficient of y^1
    Expr nonlin_coef;                   // coefficient of y^n for some n
    long n_exponent = 0;
    bool ambiguous = false;

    auto get_y_power = [&](const Expr& term, long& out_exp,
                            Expr& out_coef) -> bool {
        Expr coef = S::One();
        long exp_y = 0;
        std::vector<Expr> rest;
        if (term->type_id() == TypeId::Mul) {
            for (const auto& f : term->args()) {
                if (f == y) {
                    if (exp_y != 0) return false;  // y twice in same term
                    exp_y = 1;
                    continue;
                }
                if (f->type_id() == TypeId::Pow && f->args()[0] == y
                    && f->args()[1]->type_id() == TypeId::Integer) {
                    if (exp_y != 0) return false;
                    exp_y = static_cast<const Integer&>(*f->args()[1]).to_long();
                    continue;
                }
                if (has(f, y)) return false;
                rest.push_back(f);
            }
        } else if (term == y) {
            exp_y = 1;
        } else if (term->type_id() == TypeId::Pow && term->args()[0] == y
                   && term->args()[1]->type_id() == TypeId::Integer) {
            exp_y = static_cast<const Integer&>(*term->args()[1]).to_long();
        } else if (!has(term, y)) {
            return false;  // not a multiple of y^k
        } else {
            return false;
        }
        out_exp = exp_y;
        out_coef = mul(rest);
        return true;
    };

    for (const auto& term : e->args()) {
        long exp_y;
        Expr coef;
        if (!get_y_power(term, exp_y, coef)) {
            ambiguous = true;
            break;
        }
        if (exp_y == 1) {
            lin_coef = add(lin_coef, coef);
        } else {
            if (n_exponent != 0 && n_exponent != exp_y) { ambiguous = true; break; }
            n_exponent = exp_y;
            nonlin_coef = nonlin_coef ? add(nonlin_coef, coef) : coef;
        }
    }
    if (ambiguous || n_exponent == 0 || n_exponent == 1
        || !nonlin_coef) {
        return function_symbol("Dsolve")(eq, y, x);
    }
    // y' = lin_coef * y + nonlin_coef * y^n
    // → y' + (-lin_coef) y = nonlin_coef * y^n
    // Substitution v = y^(1-n): v' = (1-n) y^(-n) y'.
    // Resulting linear ODE in v:
    //   v' + (1-n)(-lin_coef) v = (1-n) nonlin_coef
    //   v' + ((n-1) lin_coef) v = (1-n) nonlin_coef
    Expr p = mul(integer(n_exponent - 1), lin_coef);
    Expr q = mul(integer(1 - n_exponent), nonlin_coef);
    // Solve linear ODE for v:
    Expr int_p = integrate(p, x);
    if (is_integral_marker(int_p)) return function_symbol("Dsolve")(eq, y, x);
    Expr mu = exp(int_p);
    Expr int_mu_q = integrate(mu * q, x);
    if (is_integral_marker(int_mu_q)) return function_symbol("Dsolve")(eq, y, x);
    int cnt = 0;
    Expr C = fresh_constant(cnt);
    Expr v_sol = (int_mu_q + C) / mu;
    // y = v^(1/(1-n)).
    Expr inv_exp = pow(integer(1 - n_exponent), S::NegativeOne());
    return simplify(pow(v_sol, inv_exp));
}

// ---------------------------------------------------------------------------
// Exact ODE: M(x, y) dx + N(x, y) dy = 0 with ∂M/∂y = ∂N/∂x.
//   Potential function ψ(x, y) such that ψ_x = M, ψ_y = N. Solution
//   is ψ(x, y) = C.
// ---------------------------------------------------------------------------

Expr dsolve_exact(const Expr& M, const Expr& N, const Expr& y, const Expr& x) {
    Expr My = diff(M, y);
    Expr Nx = diff(N, x);
    if (!(simplify(My - Nx) == S::Zero())) {
        return function_symbol("Dsolve")(std::vector<Expr>{M, N, y, x});
    }
    // ψ = ∫M dx + h(y), with h'(y) = N - ∂(∫M dx)/∂y.
    Expr int_M = integrate(M, x);
    if (is_integral_marker(int_M)) {
        return function_symbol("Dsolve")(std::vector<Expr>{M, N, y, x});
    }
    Expr partial_y = diff(int_M, y);
    Expr h_prime = simplify(N - partial_y);
    Expr h = integrate(h_prime, y);
    if (is_integral_marker(h)) {
        return function_symbol("Dsolve")(std::vector<Expr>{M, N, y, x});
    }
    int cnt = 0;
    Expr C = fresh_constant(cnt);
    return simplify(int_M + h - C);   // ψ(x, y) - C = 0
}

// ---------------------------------------------------------------------------
// First-order classifier — try strategies in order.
// ---------------------------------------------------------------------------

Expr dsolve_first_order(const Expr& eq, const Expr& y, const Expr& yp,
                         const Expr& x) {
    auto try_one = [&](Expr (*fn)(const Expr&, const Expr&, const Expr&,
                                    const Expr&)) -> std::optional<Expr> {
        Expr r = fn(eq, y, yp, x);
        if (r->type_id() == TypeId::Function) {
            const auto& f = static_cast<const Function&>(*r);
            if (f.name() == "Dsolve") return std::nullopt;
        }
        return r;
    };
    if (auto r = try_one(dsolve_linear_first_order); r) return *r;
    if (auto r = try_one(dsolve_separable); r) return *r;
    if (auto r = try_one(dsolve_bernoulli); r) return *r;
    return function_symbol("Dsolve")(eq, y, x);
}

// ---------------------------------------------------------------------------
// Constant-coefficient linear ODE.
//   coeffs[n] y^(n) + ... + coeffs[0] y = 0
//   → characteristic poly Σ coeffs[k] r^k. Each root r of multiplicity m
//     contributes (C0 + C1 x + ... + C_{m-1} x^{m-1}) e^(rx).
// ---------------------------------------------------------------------------

Expr dsolve_constant_coeff(const std::vector<Expr>& coeffs, const Expr& x) {
    if (coeffs.size() < 2) {
        throw std::invalid_argument("dsolve_constant_coeff: order ≥ 1 required");
    }
    auto r_sym = symbol("__r");
    Poly char_poly(std::vector<Expr>(coeffs.begin(), coeffs.end()), r_sym);
    auto roots = char_poly.roots();
    if (roots.empty()) {
        throw std::runtime_error(
            "dsolve_constant_coeff: characteristic poly has no closed-form roots");
    }
    // Group roots by multiplicity.
    std::vector<std::pair<Expr, std::size_t>> grouped;
    for (const auto& r : roots) {
        bool merged = false;
        for (auto& [val, mult] : grouped) {
            if (val == r) { ++mult; merged = true; break; }
        }
        if (!merged) grouped.emplace_back(r, 1);
    }
    int cnt = 0;
    Expr y = S::Zero();
    for (const auto& [root, mult] : grouped) {
        Expr poly_x = S::Zero();
        for (std::size_t k = 0; k < mult; ++k) {
            Expr C = fresh_constant(cnt);
            poly_x = poly_x + C * pow(x, integer(static_cast<long>(k)));
        }
        y = y + poly_x * exp(mul(root, x));
    }
    return y;
}

// ---------------------------------------------------------------------------
// Cauchy-Euler. Substitute x = e^t. The kth derivative of y w.r.t. x
// becomes a polynomial in (d/dt) operators; the resulting equation in t
// is constant-coefficient. Specifically, with D = d/dt:
//   x^k y^(k) = D(D-1)(D-2)...(D-k+1) y(t)
// So coefficients in t become polynomials in D. We expand the resulting
// polynomial and pass it to dsolve_constant_coeff.
// ---------------------------------------------------------------------------

Expr dsolve_cauchy_euler(const std::vector<Expr>& coeffs, const Expr& x) {
    if (coeffs.size() < 2) {
        throw std::invalid_argument("dsolve_cauchy_euler: order ≥ 1 required");
    }
    // Build the polynomial in D corresponding to a_k x^k y^(k) → a_k D(D-1)...(D-k+1) y.
    auto D = symbol("__D");
    // We'll accumulate the constant coefficients of the resulting poly in D.
    // Easier: form the symbolic polynomial p(D) = Σ coeffs[k] · Π_{j=0..k-1}(D - j).
    Expr p_in_D = S::Zero();
    for (std::size_t k = 0; k < coeffs.size(); ++k) {
        Expr term = coeffs[k];
        for (long j = 0; j < static_cast<long>(k); ++j) {
            term = term * (D - integer(j));
        }
        p_in_D = p_in_D + term;
    }
    // Build a Poly in D.
    Poly poly_D(expand(p_in_D), D);
    std::vector<Expr> char_coeffs(poly_D.coeffs().begin(),
                                    poly_D.coeffs().end());
    auto t = symbol("__t");
    Expr y_t = dsolve_constant_coeff(char_coeffs, t);
    // Back-substitute t = log(x).
    return simplify(subs(y_t, t, log(x)));
}

// ---------------------------------------------------------------------------
// Linear ODE system: y' = A · y for n×n constant matrix A.
//   Diagonalize A = P D P⁻¹. Solution: y(t) = P · diag(exp(λᵢ t)) · P⁻¹ · y₀
//   General: y(t) = Σᵢ Cᵢ vᵢ exp(λᵢ t) where vᵢ are eigenvectors.
// ---------------------------------------------------------------------------

Matrix dsolve_system(const Matrix& A, const Expr& x) {
    if (!A.is_square()) {
        throw std::invalid_argument("dsolve_system: A must be square");
    }
    auto evs = A.eigenvects();
    int cnt = 0;
    const std::size_t n = A.rows();
    std::vector<Expr> sol_components(n, S::Zero());
    for (const auto& ev_pair : evs) {
        const Expr& lambda = ev_pair.first;
        for (const auto& vec : ev_pair.second) {
            Expr C = fresh_constant(cnt);
            Expr factor = mul(C, exp(mul(lambda, x)));
            for (std::size_t i = 0; i < n; ++i) {
                sol_components[i] = add(sol_components[i],
                                         mul(factor, vec.at(i, 0)));
            }
        }
    }
    Matrix out(n, 1);
    for (std::size_t i = 0; i < n; ++i) out.set(i, 0, simplify(sol_components[i]));
    return out;
}

// ---------------------------------------------------------------------------
// IVP / checkodesol
// ---------------------------------------------------------------------------

Expr apply_ivp(const Expr& general_solution, const Expr& x,
                const std::vector<std::pair<Expr, Expr>>& conditions) {
    // Collect free constants present in the solution: any Symbol named
    // "__C0", "__C1", ... (a simple convention from our other generators).
    // We'll detect them via a syntactic walk.
    Expr current = general_solution;
    int idx = 0;
    for (const auto& [x_val, y_val] : conditions) {
        // Substitute x = x_val into the solution; equate to y_val; solve
        // for the next free constant.
        std::string name = "__C" + std::to_string(idx);
        auto C = symbol(name);
        Expr at = simplify(subs(current, x, x_val));
        // Solve at = y_val for C.
        auto sols = solve(at - y_val, C);
        if (sols.empty()) break;
        current = simplify(subs(current, C, sols[0]));
        ++idx;
    }
    return current;
}

Expr checkodesol(const Expr& eq, const Expr& sol, const Expr& y,
                  const Expr& yp, const Expr& x) {
    Expr sol_p = diff(sol, x);
    Expr substituted = subs(subs(eq, yp, sol_p), y, sol);
    return simplify(substituted);
}

// ---------------------------------------------------------------------------
// First-order linear PDE: a · u_x + b · u_y = c (constants).
//   General solution: u(x, y) = (c/a) · x + F(b·x - a·y) for any F.
//   We emit F as an UndefinedFunction "F" applied to the characteristic.
// ---------------------------------------------------------------------------

Expr pdsolve_first_order_linear(const Expr& a, const Expr& b, const Expr& c,
                                 const Expr& x, const Expr& y) {
    if (has(a, x) || has(a, y) || has(b, x) || has(b, y)
        || has(c, x) || has(c, y)) {
        // Variable-coefficient PDE — beyond this minimal scope.
        return function_symbol("Pdsolve")(std::vector<Expr>{a, b, c, x, y});
    }
    if (a == S::Zero()) {
        // b · u_y = c  →  u = (c/b) y + F(x)
        Expr F = function_symbol("F")(x);
        return c / b * y + F;
    }
    Expr characteristic = mul(b, x) - mul(a, y);
    Expr F = function_symbol("F")(characteristic);
    return mul(c / a, x) + F;
}

}  // namespace sympp
