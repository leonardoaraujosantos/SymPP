#include <sympp/core/rewrite.hpp>

#include <optional>
#include <vector>

#include <sympp/core/basic.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>

namespace sympp {

namespace {

[[nodiscard]] Expr recip(const Expr& e) { return pow(e, integer(-1)); }

// Exponential building blocks for argument a.
struct ExpForms {
    Expr sin_, cos_, sinh_, cosh_;
};
[[nodiscard]] ExpForms exp_forms(const Expr& a) {
    Expr eIa = exp(mul(S::I(), a));
    Expr emIa = exp(mul(integer(-1), mul(S::I(), a)));
    Expr ea = exp(a);
    Expr ema = exp(mul(integer(-1), a));
    return {
        mul(eIa - emIa, recip(mul(integer(2), S::I()))),  // sin  = (e^{ia}-e^{-ia})/(2i)
        mul(eIa + emIa, rational(1, 2)),                  // cos  = (e^{ia}+e^{-ia})/2
        mul(ea - ema, rational(1, 2)),                    // sinh = (e^a-e^{-a})/2
        mul(ea + ema, rational(1, 2)),                    // cosh = (e^a+e^{-a})/2
    };
}

// Transform a single unary trig/hyperbolic function to exponentials; nullopt if
// `id` is not one we rewrite.
[[nodiscard]] std::optional<Expr> to_exp(FunctionId id, const Expr& a) {
    ExpForms f = exp_forms(a);
    switch (id) {
        case FunctionId::Sin: return f.sin_;
        case FunctionId::Cos: return f.cos_;
        case FunctionId::Tan: return mul(f.sin_, recip(f.cos_));
        case FunctionId::Cot: return mul(f.cos_, recip(f.sin_));
        case FunctionId::Sec: return recip(f.cos_);
        case FunctionId::Csc: return recip(f.sin_);
        case FunctionId::Sinh: return f.sinh_;
        case FunctionId::Cosh: return f.cosh_;
        case FunctionId::Tanh: return mul(f.sinh_, recip(f.cosh_));
        case FunctionId::Coth: return mul(f.cosh_, recip(f.sinh_));
        case FunctionId::Sech: return recip(f.cosh_);
        case FunctionId::Csch: return recip(f.sinh_);
        default: return std::nullopt;
    }
}

// Transform tan/cot/sec/csc into sin & cos; nullopt otherwise.
[[nodiscard]] std::optional<Expr> to_sincos(FunctionId id, const Expr& a) {
    switch (id) {
        case FunctionId::Tan: return mul(sin(a), recip(cos(a)));
        case FunctionId::Cot: return mul(cos(a), recip(sin(a)));
        case FunctionId::Sec: return recip(cos(a));
        case FunctionId::Csc: return recip(sin(a));
        default: return std::nullopt;
    }
}

// factorial / binomial / (rising|falling) factorial → Γ:
//   n!                    = Γ(n+1)
//   C(n,k)                = Γ(n+1) / (Γ(k+1)·Γ(n−k+1))
//   RisingFactorial(x,n)  = Γ(x+n) / Γ(x)
//   FallingFactorial(x,n) = Γ(x+1) / Γ(x−n+1)
[[nodiscard]] std::optional<Expr> to_gamma(FunctionId id,
                                           const std::vector<Expr>& a) {
    if (id == FunctionId::Factorial && a.size() == 1) {
        return gamma(a[0] + integer(1));
    }
    if (id == FunctionId::Binomial && a.size() == 2) {
        return mul(gamma(a[0] + integer(1)),
                   recip(mul(gamma(a[1] + integer(1)),
                             gamma(a[0] - a[1] + integer(1)))));
    }
    if (id == FunctionId::RisingFactorial && a.size() == 2) {
        return mul(gamma(a[0] + a[1]), recip(gamma(a[0])));
    }
    if (id == FunctionId::FallingFactorial && a.size() == 2) {
        return mul(gamma(a[0] + integer(1)),
                   recip(gamma(a[0] - a[1] + integer(1))));
    }
    return std::nullopt;
}

// Γ / binomial / (rising|falling) factorial → factorial:
//   Γ(z)                  = (z−1)!
//   C(n,k)                = n! / (k!·(n−k)!)
//   RisingFactorial(x,n)  = (x+n−1)! / (x−1)!
//   FallingFactorial(x,n) = x! / (x−n)!
[[nodiscard]] std::optional<Expr> to_factorial(FunctionId id,
                                               const std::vector<Expr>& a) {
    if (id == FunctionId::Gamma && a.size() == 1) {
        return factorial(a[0] - integer(1));
    }
    if (id == FunctionId::Binomial && a.size() == 2) {
        return mul(factorial(a[0]),
                   recip(mul(factorial(a[1]), factorial(a[0] - a[1]))));
    }
    if (id == FunctionId::RisingFactorial && a.size() == 2) {
        return mul(factorial(a[0] + a[1] - integer(1)),
                   recip(factorial(a[0] - integer(1))));
    }
    if (id == FunctionId::FallingFactorial && a.size() == 2) {
        return mul(factorial(a[0]), recip(factorial(a[0] - a[1])));
    }
    return std::nullopt;
}

[[nodiscard]] Expr rw(const Expr& e, RewriteTarget target) {
    if (!e) return e;
    switch (e->type_id()) {
        case TypeId::Add: {
            std::vector<Expr> parts;
            for (const auto& a : e->args()) parts.push_back(rw(a, target));
            return add(std::move(parts));
        }
        case TypeId::Mul: {
            std::vector<Expr> parts;
            for (const auto& a : e->args()) parts.push_back(rw(a, target));
            return mul(std::move(parts));
        }
        case TypeId::Pow:
            return pow(rw(e->args()[0], target), rw(e->args()[1], target));
        case TypeId::Function: {
            const auto& fn = static_cast<const Function&>(*e);
            std::vector<Expr> new_args;
            for (const auto& a : fn.args()) new_args.push_back(rw(a, target));
            const FunctionId id = fn.function_id();
            std::optional<Expr> t;
            switch (target) {
                case RewriteTarget::Exp:
                    if (new_args.size() == 1) t = to_exp(id, new_args[0]);
                    break;
                case RewriteTarget::SinCos:
                    if (new_args.size() == 1) t = to_sincos(id, new_args[0]);
                    break;
                case RewriteTarget::Gamma:
                    t = to_gamma(id, new_args);
                    break;
                case RewriteTarget::Factorial:
                    t = to_factorial(id, new_args);
                    break;
            }
            if (t.has_value()) return *t;
            return fn.rebuild(std::move(new_args));
        }
        default:
            return e;
    }
}

}  // namespace

Expr rewrite(const Expr& e, RewriteTarget target) { return rw(e, target); }

}  // namespace sympp
