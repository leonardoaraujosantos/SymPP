#include <sympp/calculus/extrema.hpp>

#include <cmath>
#include <utility>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/simplify/simplify.hpp>
#include <sympp/solvers/solve.hpp>

namespace sympp {

namespace {

// Try to render a value as a finite double for ordering. Returns nullopt
// when the value isn't reducible to a real Float.
[[nodiscard]] std::optional<double> as_double(const Expr& e) {
    if (e->type_id() == TypeId::Integer) {
        return static_cast<const Integer&>(*e).value().get_d();
    }
    if (e->type_id() == TypeId::Rational) {
        return static_cast<const Rational&>(*e).value().get_d();
    }
    Expr v = evalf(e, 30);
    if (v->type_id() == TypeId::Float) {
        try { return std::stod(v->str()); } catch (...) { return std::nullopt; }
    }
    return std::nullopt;
}

}  // namespace

std::vector<Expr> stationary_points(const Expr& f, const Expr& var) {
    return solve(diff(f, var), var);
}

Expr minimum(const Expr& f, const Expr& var) {
    auto pts = stationary_points(f, var);
    if (pts.empty()) return f;
    Expr best = simplify(subs(f, var, pts[0]));
    auto best_d = as_double(best);
    for (std::size_t i = 1; i < pts.size(); ++i) {
        Expr v = simplify(subs(f, var, pts[i]));
        auto d = as_double(v);
        if (best_d && d && *d < *best_d) {
            best = v;
            best_d = d;
        } else if (!best_d && d) {
            // Prefer a numerically-comparable value over an opaque one.
            best = v;
            best_d = d;
        }
    }
    return best;
}

Expr maximum(const Expr& f, const Expr& var) {
    auto pts = stationary_points(f, var);
    if (pts.empty()) return f;
    Expr best = simplify(subs(f, var, pts[0]));
    auto best_d = as_double(best);
    for (std::size_t i = 1; i < pts.size(); ++i) {
        Expr v = simplify(subs(f, var, pts[i]));
        auto d = as_double(v);
        if (best_d && d && *d > *best_d) {
            best = v;
            best_d = d;
        } else if (!best_d && d) {
            best = v;
            best_d = d;
        }
    }
    return best;
}

std::optional<bool> is_increasing(const Expr& f, const Expr& var) {
    Expr fp = simplify(diff(f, var));
    if (auto p = is_positive(fp); p && *p) return true;
    if (auto n = is_negative(fp); n && *n) return false;
    return std::nullopt;
}

std::optional<bool> is_decreasing(const Expr& f, const Expr& var) {
    Expr fp = simplify(diff(f, var));
    if (auto n = is_negative(fp); n && *n) return true;
    if (auto p = is_positive(fp); p && *p) return false;
    return std::nullopt;
}

}  // namespace sympp
