#include <sympp/calculus/asymptotes.hpp>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/solvers/solve.hpp>

namespace sympp {

std::vector<Expr> vertical_asymptotes(const Expr& f, const Expr& var) {
    // 1. cancel — strip shared factors so removable singularities are gone.
    Expr c = cancel(f, var);
    // 2. together — present as num * den^(-1) (or just Pow(den, -1) when
    //    num == 1).
    Expr t = together(c);
    Expr denom;
    if (t->type_id() == TypeId::Pow
        && t->args()[1] == S::NegativeOne()) {
        denom = t->args()[0];
    } else if (t->type_id() == TypeId::Mul) {
        for (const auto& factor : t->args()) {
            if (factor->type_id() == TypeId::Pow
                && factor->args()[1] == S::NegativeOne()) {
                denom = factor->args()[0];
                break;
            }
        }
    }
    if (!denom) return {};
    return solve(denom, var);
}

std::vector<Expr> inflection_points(const Expr& f, const Expr& var) {
    Expr fpp = diff(diff(f, var), var);
    return solve(fpp, var);
}

}  // namespace sympp
