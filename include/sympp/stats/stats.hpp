#pragma once

// Probability distributions with symbolic moments and densities — a compact
// analogue of sympy.stats.
//
//   auto N = normal(integer(0), integer(1));
//   mean(N);          // → 0
//   variance(N);      // → 1
//   pdf(N, x);        // → sqrt(2)*exp(-x**2/2)/(2*sqrt(pi))
//
// Continuous (Normal, Uniform, Exponential) expose pdf and cdf; discrete
// (Bernoulli, Binomial, Poisson) expose the probability mass via pdf and have
// no cdf (throws). mean/variance are available for all.
//
// Reference: sympy/stats/{crv_types,drv_types}.py.

#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp::stats {

class SYMPP_EXPORT Distribution {
public:
    enum class Kind { Normal, Uniform, Exponential, Bernoulli, Binomial, Poisson };
    Distribution(Kind k, std::vector<Expr> params) : kind_(k), params_(std::move(params)) {}
    [[nodiscard]] Kind kind() const noexcept { return kind_; }
    [[nodiscard]] const std::vector<Expr>& params() const noexcept { return params_; }
    [[nodiscard]] bool is_discrete() const noexcept {
        return kind_ == Kind::Bernoulli || kind_ == Kind::Binomial || kind_ == Kind::Poisson;
    }

private:
    Kind kind_;
    std::vector<Expr> params_;
};

// Constructors.
[[nodiscard]] SYMPP_EXPORT Distribution normal(const Expr& mean, const Expr& std_dev);
[[nodiscard]] SYMPP_EXPORT Distribution uniform(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Distribution exponential(const Expr& rate);
[[nodiscard]] SYMPP_EXPORT Distribution bernoulli(const Expr& p);
[[nodiscard]] SYMPP_EXPORT Distribution binomial(const Expr& n, const Expr& p);
[[nodiscard]] SYMPP_EXPORT Distribution poisson(const Expr& lambda);

// Moments.
[[nodiscard]] SYMPP_EXPORT Expr mean(const Distribution& d);
[[nodiscard]] SYMPP_EXPORT Expr variance(const Distribution& d);

// Density (continuous) / probability mass (discrete) at x.
[[nodiscard]] SYMPP_EXPORT Expr pdf(const Distribution& d, const Expr& x);
// Cumulative distribution at x (continuous only; throws for discrete).
[[nodiscard]] SYMPP_EXPORT Expr cdf(const Distribution& d, const Expr& x);

}  // namespace sympp::stats
