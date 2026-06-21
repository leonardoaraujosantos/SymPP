#include <sympp/stats/stats.hpp>

#include <stdexcept>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/special.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp::stats {

using Kind = Distribution::Kind;

Distribution normal(const Expr& mean, const Expr& std_dev) {
    return {Kind::Normal, {mean, std_dev}};
}
Distribution uniform(const Expr& a, const Expr& b) { return {Kind::Uniform, {a, b}}; }
Distribution exponential(const Expr& rate) { return {Kind::Exponential, {rate}}; }
Distribution bernoulli(const Expr& p) { return {Kind::Bernoulli, {p}}; }
Distribution binomial(const Expr& n, const Expr& p) { return {Kind::Binomial, {n, p}}; }
Distribution poisson(const Expr& lambda) { return {Kind::Poisson, {lambda}}; }

namespace {
[[nodiscard]] Expr sq(const Expr& e) { return pow(e, integer(2)); }
[[nodiscard]] Expr one_minus(const Expr& p) { return integer(1) - p; }
}  // namespace

Expr mean(const Distribution& d) {
    const auto& p = d.params();
    switch (d.kind()) {
        case Kind::Normal: return p[0];
        case Kind::Uniform: return simplify((p[0] + p[1]) * rational(1, 2));
        case Kind::Exponential: return pow(p[0], integer(-1));
        case Kind::Bernoulli: return p[0];
        case Kind::Binomial: return p[0] * p[1];
        case Kind::Poisson: return p[0];
    }
    throw std::logic_error("mean: unknown distribution");
}

Expr variance(const Distribution& d) {
    const auto& p = d.params();
    switch (d.kind()) {
        case Kind::Normal: return sq(p[1]);
        case Kind::Uniform: return simplify(sq(p[1] - p[0]) * rational(1, 12));
        case Kind::Exponential: return pow(p[0], integer(-2));
        case Kind::Bernoulli: return p[0] * one_minus(p[0]);
        case Kind::Binomial: return p[0] * p[1] * one_minus(p[1]);
        case Kind::Poisson: return p[0];
    }
    throw std::logic_error("variance: unknown distribution");
}

Expr pdf(const Distribution& d, const Expr& x) {
    const auto& p = d.params();
    switch (d.kind()) {
        case Kind::Normal: {
            Expr mu = p[0], s = p[1];
            Expr norm = pow(s * sqrt(integer(2) * S::Pi()), integer(-1));
            Expr e = exp(mul(integer(-1), sq(x - mu)) * pow(integer(2) * sq(s), integer(-1)));
            return simplify(norm * e);
        }
        case Kind::Uniform:
            return simplify(pow(p[1] - p[0], integer(-1)));  // density on the support
        case Kind::Exponential:
            return simplify(p[0] * exp(mul(integer(-1), p[0]) * x));
        case Kind::Bernoulli:  // p^x (1-p)^(1-x)
            return pow(p[0], x) * pow(one_minus(p[0]), one_minus(x));
        case Kind::Binomial:  // C(n,x) p^x (1-p)^(n-x)
            return sympp::binomial(p[0], x) * pow(p[1], x) * pow(one_minus(p[1]), p[0] - x);
        case Kind::Poisson:  // λ^x e^{-λ} / x!
            return p[0] == S::Zero()
                       ? S::Zero()
                       : pow(p[0], x) * exp(mul(integer(-1), p[0])) *
                             pow(factorial(x), integer(-1));
    }
    throw std::logic_error("pdf: unknown distribution");
}

Expr cdf(const Distribution& d, const Expr& x) {
    const auto& p = d.params();
    switch (d.kind()) {
        case Kind::Normal: {
            Expr z = (x - p[0]) * pow(p[1] * sqrt(integer(2)), integer(-1));
            return simplify((integer(1) + erf(z)) * rational(1, 2));
        }
        case Kind::Uniform:
            return simplify((x - p[0]) * pow(p[1] - p[0], integer(-1)));
        case Kind::Exponential:
            return simplify(integer(1) - exp(mul(integer(-1), p[0]) * x));
        case Kind::Bernoulli:
        case Kind::Binomial:
        case Kind::Poisson:
            throw std::runtime_error("cdf: not provided for discrete distributions");
    }
    throw std::logic_error("cdf: unknown distribution");
}

}  // namespace sympp::stats
