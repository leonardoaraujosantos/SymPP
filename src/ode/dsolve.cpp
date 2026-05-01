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
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
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

// Try to express the equation `lhs_int(y) = rhs_in_x + C` as an explicit
// y = h(x). Patterns supported:
//   * lhs_int == log(y)        →  y = exp(rhs)
//   * lhs_int == y             →  y = rhs
//   * lhs_int == coef · y^k    →  y = (rhs/coef)^(1/k)  (for k ≠ 0)
// Returns nullopt when no pattern matches.
[[nodiscard]] std::optional<Expr> invert_for_y(
    const Expr& lhs_int, const Expr& rhs_in_x, const Expr& y) {
    // log(y)
    if (lhs_int->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*lhs_int);
        if (fn.function_id() == FunctionId::Log && fn.args().size() == 1
            && fn.args()[0] == y) {
            return simplify(exp(rhs_in_x));
        }
    }
    // y itself
    if (lhs_int == y) return simplify(rhs_in_x);
    // y^k or Mul(coef, y^k)
    auto match_y_power = [&y](const Expr& e) -> std::optional<std::pair<Expr, Expr>> {
        // returns (coef, k) such that e == coef · y^k
        if (e->type_id() == TypeId::Pow && e->args()[0] == y) {
            return std::pair{S::One(), e->args()[1]};
        }
        if (e == y) return std::pair{S::One(), S::One()};
        if (e->type_id() == TypeId::Mul) {
            Expr coef = S::One();
            std::optional<Expr> k;
            for (const auto& f : e->args()) {
                if (f == y) {
                    if (k) return std::nullopt;
                    k = S::One();
                } else if (f->type_id() == TypeId::Pow && f->args()[0] == y) {
                    if (k) return std::nullopt;
                    k = f->args()[1];
                } else if (has(f, y)) {
                    return std::nullopt;
                } else {
                    coef = coef * f;
                }
            }
            if (k) return std::pair{coef, *k};
        }
        return std::nullopt;
    };
    if (auto m = match_y_power(lhs_int); m) {
        const Expr& coef = m->first;
        const Expr& k = m->second;
        if (k == S::Zero()) return std::nullopt;
        // y^k = rhs / coef → y = (rhs / coef)^(1/k)
        Expr inner = rhs_in_x / coef;
        Expr inv_k = pow(k, S::NegativeOne());
        return simplify(pow(inner, inv_k));
    }
    return std::nullopt;
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
    // Expand first so Mul(Pow(x, -1), Add(... yp ...)) flattens into
    // separate yp-bearing terms, letting Poly see the structure.
    Poly p(expand(simplify(eq)), yp);
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
    // Try the structural inversion helper first — it covers log(y), y,
    // and coef·y^k patterns. For log(y) specifically we absorb the
    // constant of integration: y = C · exp(rhs_int).
    if (lhs_int->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*lhs_int);
        if (fn.function_id() == FunctionId::Log && fn.args().size() == 1
            && fn.args()[0] == y) {
            return simplify(C * exp(rhs_int));
        }
    }
    if (auto inv = invert_for_y(lhs_int, rhs_int + C, y); inv) {
        return *inv;
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
    if (auto r = try_one(dsolve_lie_autonomous); r) return *r;
    if (auto r = try_one(dsolve_homogeneous); r) return *r;
    if (auto r = try_one(dsolve_bernoulli); r) return *r;
    return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
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
// Riccati: y' = P + Qy + Ry²
//   Substitute y = -u'/(R u). Then the linear ODE for u is:
//     u'' - (R'/R + Q) u' + R P u = 0
//   which is constant-coefficient when R, Q, R'/R + Q, R P are constants.
// ---------------------------------------------------------------------------

Expr dsolve_riccati(const Expr& P, const Expr& Q, const Expr& R, const Expr& x) {
    // Build the linear 2nd-order ODE coefficients: a₀ y + a₁ y' + a₂ y'' = 0.
    // u'' - (R'/R + Q) u' + R P u = 0 → a₂ = 1, a₁ = -(R'/R + Q), a₀ = R*P.
    Expr R_prime = simplify(diff(R, x));
    Expr coef_u_prime = simplify(mul(S::NegativeOne(), R_prime / R + Q));
    Expr coef_u = simplify(mul(R, P));
    // dsolve_constant_coeff requires constant coefficients. For the
    // common case where P, Q, R are constants in x, this works.
    if (has(coef_u_prime, x) || has(coef_u, x)) {
        return function_symbol("Dsolve")(std::vector<Expr>{P, Q, R, x});
    }
    Expr u = dsolve_constant_coeff({coef_u, coef_u_prime, S::One()}, x);
    // y = -u' / (R u)
    Expr u_prime = diff(u, x);
    return simplify(mul(S::NegativeOne(), u_prime) / (R * u));
}

// ---------------------------------------------------------------------------
// Homogeneous y' = f(y/x). Detect: substitute y → v·x in rhs; if the
// result is independent of x then rhs == f(v). Then v + x·v' = f(v) is
// separable as x·dv/(f(v) - v) = dx/x.
// ---------------------------------------------------------------------------

Expr dsolve_homogeneous(const Expr& eq, const Expr& y, const Expr& yp,
                         const Expr& x) {
    auto rhs = isolate_yp(eq, yp);
    if (!rhs) return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
    auto v = symbol("__hom_v");
    Expr rhs_in_v = simplify(subs(*rhs, y, mul(v, x)));
    if (has(rhs_in_v, x)) {
        return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
    }
    // Now rhs_in_v == f(v). Separable in (v, x): integrate
    //   1 / (f(v) - v) dv = ∫ 1/x dx + C  =  log(x) + C
    Expr denom = simplify(rhs_in_v - v);
    if (denom == S::Zero()) {
        // f(v) - v = 0 ⇒ v' = 0 ⇒ v = C ⇒ y = C·x.
        int cnt0 = 0;
        Expr C = fresh_constant(cnt0);
        return mul(C, x);
    }
    Expr lhs_int = integrate(pow(denom, S::NegativeOne()), v);
    if (is_integral_marker(lhs_int)) {
        return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
    }
    int cnt = 0;
    Expr C = fresh_constant(cnt);
    // Implicit: lhs_int(v) - log(x) - C = 0, with v = y/x.
    Expr implicit_in_v = simplify(lhs_int - log(x) - C);
    Expr implicit = subs(implicit_in_v, v, y / x);
    // Try explicit y when invertible.
    auto sols = solve(implicit, y);
    if (!sols.empty()) return simplify(sols[0]);
    return implicit;
}

// ---------------------------------------------------------------------------
// Lie symmetry — autonomous ODE y' = f(y), x-translation invariant.
// Reduces to separable: dy/f(y) = dx + C.
// ---------------------------------------------------------------------------

Expr dsolve_lie_autonomous(const Expr& eq, const Expr& y, const Expr& yp,
                            const Expr& x) {
    auto rhs = isolate_yp(eq, yp);
    if (!rhs) return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
    if (has(*rhs, x)) {
        return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
    }
    // y' = f(y) ⇒ ∫ 1/f(y) dy = x + C.
    Expr lhs_int = integrate(pow(*rhs, S::NegativeOne()), y);
    if (is_integral_marker(lhs_int)) {
        return function_symbol("Dsolve")(std::vector<Expr>{eq, y, x});
    }
    int cnt = 0;
    Expr C = fresh_constant(cnt);
    // log(y) absorbs C as a multiplicative constant: y = C · exp(x).
    if (lhs_int->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*lhs_int);
        if (fn.function_id() == FunctionId::Log && fn.args().size() == 1
            && fn.args()[0] == y) {
            return simplify(C * exp(x));
        }
    }
    if (auto inv = invert_for_y(lhs_int, x + C, y); inv) return *inv;
    Expr implicit = simplify(lhs_int - x - C);
    auto sols = solve(implicit, y);
    if (!sols.empty()) return simplify(sols[0]);
    return implicit;
}

// ---------------------------------------------------------------------------
// Hypergeometric ODE recognition:
//   x(1-x) y'' + (c - (a+b+1)x) y' - a·b·y = 0
// The user passes [coef_y, coef_yp, coef_ypp] as Exprs in x. We match
// the standard form by extracting polynomial coefficients in x.
// ---------------------------------------------------------------------------

Expr hyper(const Expr& a, const Expr& b, const Expr& c, const Expr& z) {
    return function_symbol("hyper")(std::vector<Expr>{a, b, c, z});
}

Expr dsolve_hypergeometric(const std::vector<Expr>& coeffs_y_ypp,
                            const Expr& x) {
    if (coeffs_y_ypp.size() != 3) {
        throw std::invalid_argument(
            "dsolve_hypergeometric: need [coef_y, coef_yp, coef_ypp]");
    }
    // coef_ypp should be x(1-x) = x - x²
    Poly p_ypp(expand(coeffs_y_ypp[2]), x);
    if (p_ypp.degree() != 2) {
        return function_symbol("Dsolve")(std::vector<Expr>{x});
    }
    // Expect coefficients [0, 1, -1] (i.e. x - x²) — check ratio.
    // Allow leading coefficient = 1 for now (caller normalizes).
    if (!(p_ypp.coeffs()[0] == S::Zero())
        || !(p_ypp.coeffs()[1] == S::One())
        || !(p_ypp.coeffs()[2] == S::NegativeOne())) {
        return function_symbol("Dsolve")(std::vector<Expr>{x});
    }
    // coef_yp = c - (a+b+1)*x
    Poly p_yp(expand(coeffs_y_ypp[1]), x);
    if (p_yp.degree() != 1) {
        return function_symbol("Dsolve")(std::vector<Expr>{x});
    }
    Expr c_param = p_yp.coeffs()[0];
    Expr neg_apb1 = p_yp.coeffs()[1];                      // -(a+b+1)
    Expr a_plus_b_plus_1 = mul(S::NegativeOne(), neg_apb1);
    Expr a_plus_b = a_plus_b_plus_1 - integer(1);
    // coef_y = -a·b
    Poly p_y(expand(coeffs_y_ypp[0]), x);
    if (p_y.degree() != 0) {
        return function_symbol("Dsolve")(std::vector<Expr>{x});
    }
    Expr neg_ab = p_y.coeffs()[0];
    Expr ab = mul(S::NegativeOne(), neg_ab);
    // a, b are roots of t² - (a+b)t + ab = 0.
    auto a_var = symbol("__hyper_t");
    Poly char_poly({ab, mul(S::NegativeOne(), a_plus_b), S::One()}, a_var);
    auto roots = char_poly.roots();
    if (roots.size() != 2) {
        return function_symbol("Dsolve")(std::vector<Expr>{x});
    }
    return hyper(roots[0], roots[1], c_param, x);
}

// ---------------------------------------------------------------------------
// Variable-coefficient first-order linear PDE: a(x,y)·u_x + b(x,y)·u_y = c(x,y)
// Method of characteristics: dy/dx = b/a along characteristic curves.
// We only handle the homogeneous case c == 0 here (otherwise we'd need
// to integrate du/dx = c/a along the characteristic).
// ---------------------------------------------------------------------------

Expr pdsolve_first_order_variable(const Expr& a, const Expr& b, const Expr& c,
                                    const Expr& y, const Expr& yp,
                                    const Expr& x) {
    if (!(c == S::Zero())) {
        return function_symbol("Pdsolve")(std::vector<Expr>{a, b, c, x, y});
    }
    // dy/dx = b/a → solve as ODE; first-order classifier handles common cases.
    Expr characteristic_eq = yp - b / a;
    Expr y_of_x = dsolve_first_order(characteristic_eq, y, yp, x);
    // The implicit form y - <general sol in x> = 0 defines a one-parameter
    // family; the constant of integration parameterizes characteristics.
    // u(x, y) = F(C) where C is the constant we'd need to solve for in the
    // characteristic. For the simple case where y_of_x = exp(...) * __C0,
    // we extract the characteristic invariant.
    // Coarse approach: u(x, y) = F(<the implicit relation>).
    auto C0 = symbol("__C0");
    auto sols = solve(y - y_of_x, C0);
    if (sols.empty()) {
        return function_symbol("Pdsolve")(std::vector<Expr>{a, b, c, x, y});
    }
    Expr characteristic = sols[0];
    Expr F = function_symbol("F")(characteristic);
    return F;
}

// ---------------------------------------------------------------------------
// Heat equation u_t = k·u_xx: separation of variables.
//   u(x, t) = X(x) T(t), then T'/T = k X''/X = -λ.
//   T(t) = exp(-k λ t); X(x) = A sin(√λ x) + B cos(√λ x).
// ---------------------------------------------------------------------------

Expr pdsolve_heat(const Expr& k, const Expr& lambda, const Expr& x,
                   const Expr& t) {
    int cnt = 0;
    Expr A = fresh_constant(cnt);
    Expr B = fresh_constant(cnt);
    Expr T = exp(mul({S::NegativeOne(), k, lambda, t}));
    Expr sqrtl = sqrt(lambda);
    Expr X = A * sin(sqrtl * x) + B * cos(sqrtl * x);
    return T * X;
}

// ---------------------------------------------------------------------------
// Wave equation u_tt = c²·u_xx: d'Alembert.
//   u(x, t) = F(x - c·t) + G(x + c·t).
// ---------------------------------------------------------------------------

Expr pdsolve_wave(const Expr& c, const Expr& x, const Expr& t) {
    Expr F = function_symbol("F")(x - mul(c, t));
    Expr G = function_symbol("G")(x + mul(c, t));
    return F + G;
}

// ---------------------------------------------------------------------------
// DAE Jacobian + structural index
// ---------------------------------------------------------------------------

std::pair<Matrix, Matrix> dae_jacobians(
    const std::vector<Expr>& F, const std::vector<Expr>& y,
    const std::vector<Expr>& yp) {
    // M = ∂F/∂yp, K = ∂F/∂y. Both are |F|×|y| matrices.
    if (y.size() != yp.size()) {
        throw std::invalid_argument(
            "dae_jacobians: y and yp must have the same length");
    }
    Matrix M = jacobian(F, yp);
    Matrix K = jacobian(F, y);
    return {M, K};
}

std::size_t dae_structural_index(
    const std::vector<Expr>& F, const std::vector<Expr>& y,
    const std::vector<Expr>& yp) {
    auto [M, K] = dae_jacobians(F, y, yp);
    Expr d = simplify(M.det());
    if (!(d == S::Zero())) return 0;  // index-1 (effectively ODE)
    // Coarse: if M is singular, we'd need to differentiate constraint
    // equations to recover a solvable system. Full Pantelides analysis
    // is deep-deferred; report 1 for "needs differentiation" when at
    // least one row of M is identically zero (algebraic constraint).
    std::size_t algebraic_rows = 0;
    for (std::size_t i = 0; i < M.rows(); ++i) {
        bool all_zero = true;
        for (std::size_t j = 0; j < M.cols(); ++j) {
            if (!(M.at(i, j) == S::Zero())) { all_zero = false; break; }
        }
        if (all_zero) ++algebraic_rows;
    }
    return algebraic_rows > 0 ? 1 : 2;  // ≥2 → full Pantelides needed
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
