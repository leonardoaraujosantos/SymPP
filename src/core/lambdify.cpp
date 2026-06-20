#include <sympp/core/lambdify.hpp>

#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <mpfr.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number_symbol.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// A compiled node: maps the variable vector (as a raw pointer) to a value.
using Node = std::function<double(const double*)>;
using Index = std::unordered_map<std::string, std::size_t>;

[[noreturn]] void unsupported(const std::string& what) {
    throw std::runtime_error("lambdify: unsupported in real double evaluation: " + what);
}

[[nodiscard]] double number_symbol_value(NumberSymbolKind k) {
    switch (k) {
        case NumberSymbolKind::Pi: return M_PI;
        case NumberSymbolKind::E: return M_E;
        case NumberSymbolKind::EulerGamma: return 0.57721566490153286060651209;
        case NumberSymbolKind::Catalan: return 0.91596559417721901505460351;
    }
    unsupported("number symbol");
}

[[nodiscard]] Node compile(const Expr& e, const Index& idx);

// Fold a child list with an associative binary op and an identity element.
template <class Op>
[[nodiscard]] Node fold(const Expr& e, const Index& idx, double identity, Op op) {
    std::vector<Node> parts;
    parts.reserve(e->args().size());
    for (const auto& a : e->args()) parts.push_back(compile(a, idx));
    return [parts = std::move(parts), identity, op](const double* v) {
        double acc = identity;
        for (const auto& p : parts) acc = op(acc, p(v));
        return acc;
    };
}

[[nodiscard]] Node compile_unary(const Expr& arg, const Index& idx, double (*fn)(double)) {
    auto a = compile(arg, idx);
    return [a, fn](const double* v) { return fn(a(v)); };
}

[[nodiscard]] Node compile_function(const Expr& e, const Index& idx) {
    const auto& fn = static_cast<const Function&>(*e);
    const auto& args = fn.args();
    auto unary = [&](double (*f)(double)) { return compile_unary(args[0], idx, f); };

    switch (fn.function_id()) {
        case FunctionId::Sin: return unary(std::sin);
        case FunctionId::Cos: return unary(std::cos);
        case FunctionId::Tan: return unary(std::tan);
        case FunctionId::Asin: return unary(std::asin);
        case FunctionId::Acos: return unary(std::acos);
        case FunctionId::Atan: return unary(std::atan);
        case FunctionId::Sinh: return unary(std::sinh);
        case FunctionId::Cosh: return unary(std::cosh);
        case FunctionId::Tanh: return unary(std::tanh);
        case FunctionId::Asinh: return unary(std::asinh);
        case FunctionId::Acosh: return unary(std::acosh);
        case FunctionId::Atanh: return unary(std::atanh);
        case FunctionId::Exp: return unary(std::exp);
        case FunctionId::Abs: return unary(std::fabs);
        case FunctionId::Floor: return unary(std::floor);
        case FunctionId::Ceiling: return unary(std::ceil);
        case FunctionId::Gamma: return unary(std::tgamma);
        case FunctionId::LogGamma: return unary(std::lgamma);
        case FunctionId::Erf: return unary(std::erf);
        case FunctionId::Erfc: return unary(std::erfc);
        // Reciprocal trig / hyperbolic.
        case FunctionId::Cot: { auto a = compile(args[0], idx);
            return [a](const double* v) { return 1.0 / std::tan(a(v)); }; }
        case FunctionId::Sec: { auto a = compile(args[0], idx);
            return [a](const double* v) { return 1.0 / std::cos(a(v)); }; }
        case FunctionId::Csc: { auto a = compile(args[0], idx);
            return [a](const double* v) { return 1.0 / std::sin(a(v)); }; }
        case FunctionId::Coth: { auto a = compile(args[0], idx);
            return [a](const double* v) { return 1.0 / std::tanh(a(v)); }; }
        case FunctionId::Sech: { auto a = compile(args[0], idx);
            return [a](const double* v) { return 1.0 / std::cosh(a(v)); }; }
        case FunctionId::Csch: { auto a = compile(args[0], idx);
            return [a](const double* v) { return 1.0 / std::sinh(a(v)); }; }
        // Complex projections on a real argument.
        case FunctionId::Re:
        case FunctionId::Conjugate: return compile(args[0], idx);
        case FunctionId::Im: { (void)args; return [](const double*) { return 0.0; }; }
        case FunctionId::Sign: { auto a = compile(args[0], idx);
            return [a](const double* v) { double x = a(v); return (x > 0) - (x < 0); }; }
        // log(x) or log(x, base).
        case FunctionId::Log: {
            auto a = compile(args[0], idx);
            if (args.size() == 1) return [a](const double* v) { return std::log(a(v)); };
            auto base = compile(args[1], idx);
            return [a, base](const double* v) { return std::log(a(v)) / std::log(base(v)); };
        }
        case FunctionId::Atan2: {
            auto y = compile(args[0], idx);
            auto x = compile(args[1], idx);
            return [y, x](const double* v) { return std::atan2(y(v), x(v)); };
        }
        case FunctionId::Min: return fold(e, idx, HUGE_VAL,
            [](double a, double b) { return std::fmin(a, b); });
        case FunctionId::Max: return fold(e, idx, -HUGE_VAL,
            [](double a, double b) { return std::fmax(a, b); });
        case FunctionId::Factorial: { auto a = compile(args[0], idx);
            return [a](const double* v) { return std::tgamma(a(v) + 1.0); }; }
        default:
            unsupported("function '" + std::string(fn.name()) + "'");
    }
}

[[nodiscard]] Node compile(const Expr& e, const Index& idx) {
    if (!e) unsupported("null expression");
    switch (e->type_id()) {
        case TypeId::Integer: {
            double c = static_cast<const Integer&>(*e).value().get_d();
            return [c](const double*) { return c; };
        }
        case TypeId::Rational: {
            double c = static_cast<const Rational&>(*e).value().get_d();
            return [c](const double*) { return c; };
        }
        case TypeId::Float: {
            double c = mpfr_get_d(static_cast<const Float&>(*e).value(), MPFR_RNDN);
            return [c](const double*) { return c; };
        }
        case TypeId::NumberSymbol: {
            double c = number_symbol_value(static_cast<const NumberSymbol&>(*e).kind());
            return [c](const double*) { return c; };
        }
        case TypeId::Symbol: {
            const auto& name = static_cast<const Symbol&>(*e).name();
            auto it = idx.find(name);
            if (it == idx.end()) unsupported("free symbol '" + name + "'");
            std::size_t i = it->second;
            return [i](const double* v) { return v[i]; };
        }
        case TypeId::Add:
            return fold(e, idx, 0.0, [](double a, double b) { return a + b; });
        case TypeId::Mul:
            return fold(e, idx, 1.0, [](double a, double b) { return a * b; });
        case TypeId::Pow: {
            auto base = compile(e->args()[0], idx);
            auto exp = compile(e->args()[1], idx);
            return [base, exp](const double* v) { return std::pow(base(v), exp(v)); };
        }
        case TypeId::Function:
            return compile_function(e, idx);
        case TypeId::Infinity: return [](const double*) { return HUGE_VAL; };
        case TypeId::NegativeInfinity: return [](const double*) { return -HUGE_VAL; };
        case TypeId::NaN: return [](const double*) { return std::nan(""); };
        case TypeId::ImaginaryUnit:
            unsupported("imaginary unit I");
        default:
            unsupported("node type");
    }
}

}  // namespace

CompiledExpr lambdify(const std::vector<Expr>& vars, const Expr& body) {
    Index idx;
    for (std::size_t i = 0; i < vars.size(); ++i) {
        if (!vars[i] || vars[i]->type_id() != TypeId::Symbol) {
            throw std::runtime_error("lambdify: variables must be Symbols");
        }
        idx.emplace(static_cast<const Symbol&>(*vars[i]).name(), i);
    }
    Node program = compile(body, idx);
    std::size_t arity = vars.size();
    return [program = std::move(program), arity](const std::vector<double>& args) {
        if (args.size() != arity) {
            throw std::runtime_error("lambdify: wrong number of arguments");
        }
        return program(args.data());
    };
}

CompiledExpr lambdify(const Expr& var, const Expr& body) {
    return lambdify(std::vector<Expr>{var}, body);
}

}  // namespace sympp
