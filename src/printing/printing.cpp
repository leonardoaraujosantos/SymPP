#include <sympp/printing/printing.hpp>

#include <algorithm>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number_symbol.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp::printing {

// ---------------------------------------------------------------------------
// CodePrinter — visitor base for textual code generation.
//
// Each subclass overrides the small handful of symbols that differ across
// languages: function name lookup, power operator, comment syntax, number
// formatting. Operator precedence and the recursive walk are shared.
// ---------------------------------------------------------------------------

namespace {

// Operator precedence — lower number = lower precedence (binds less
// tightly). Loosely follows C precedence for printing comparability.
enum class Prec : int {
    None = 0,
    Add = 10,
    Mul = 20,
    Pow = 30,
    Atom = 100,
};

class CodePrinter {
public:
    virtual ~CodePrinter() = default;

    [[nodiscard]] std::string print(const Expr& e) {
        return walk(e, Prec::None);
    }

protected:
    // Look up the language-specific name for a function id. Empty string
    // means "fall back to the symbolic name".
    [[nodiscard]] virtual std::string func_name(FunctionId id) const = 0;

    // Power: ccode wraps with pow(), Fortran/Latex/Octave use **.
    [[nodiscard]] virtual std::string format_pow(const std::string& base,
                                                   const std::string& exp) const = 0;

    // Number-symbol literals (Pi, E).
    [[nodiscard]] virtual std::string format_pi() const = 0;
    [[nodiscard]] virtual std::string format_e() const = 0;

    // Format an integer value.
    [[nodiscard]] virtual std::string format_integer(const std::string& s) const {
        return s;
    }

    // Format a rational p/q.
    [[nodiscard]] virtual std::string format_rational(const std::string& num,
                                                       const std::string& den) const {
        return num + "/" + den;
    }

    // Format an Imaginary unit literal.
    [[nodiscard]] virtual std::string format_imaginary_unit() const = 0;

    // Format the * operator. Most languages use "*", Octave uses ".*".
    [[nodiscard]] virtual std::string mul_op() const { return "*"; }

    // Format the / operator. Most languages use "/", Octave uses "./".
    [[nodiscard]] virtual std::string div_op() const { return "/"; }

    // Format a symbol name (LaTeX may convert greek letters).
    [[nodiscard]] virtual std::string format_symbol(const std::string& name) const {
        return name;
    }

    // Override to format Add specially (LaTeX wraps minuses, etc.).
    [[nodiscard]] virtual std::string format_add(
        const std::vector<std::string>& terms) const {
        std::string out;
        for (std::size_t i = 0; i < terms.size(); ++i) {
            if (i == 0) out = terms[i];
            else {
                if (terms[i].size() > 0 && terms[i][0] == '-') {
                    out += " - " + terms[i].substr(1);
                } else {
                    out += " + " + terms[i];
                }
            }
        }
        return out;
    }

    // The recursive printer.
    [[nodiscard]] virtual std::string walk(const Expr& e, Prec parent) {
        if (!e) return "0";
        switch (e->type_id()) {
            case TypeId::Integer: {
                const auto& z = static_cast<const Integer&>(*e);
                return format_integer(z.value().get_str());
            }
            case TypeId::Rational: {
                const auto& q = static_cast<const Rational&>(*e);
                std::string n = q.numerator().get_str();
                std::string d = q.denominator().get_str();
                std::string r = format_rational(n, d);
                if (parent > Prec::Mul) r = "(" + r + ")";
                return r;
            }
            case TypeId::Float: {
                const auto& f = static_cast<const Float&>(*e);
                return f.str();  // already in decimal scientific form
            }
            case TypeId::NumberSymbol: {
                const auto& sym = static_cast<const NumberSymbol&>(*e);
                switch (sym.kind()) {
                    case NumberSymbolKind::Pi: return format_pi();
                    case NumberSymbolKind::E:  return format_e();
                    default: return e->str();
                }
            }
            case TypeId::ImaginaryUnit:
                return format_imaginary_unit();
            case TypeId::Symbol:
                return format_symbol(e->str());
            case TypeId::Add: {
                std::vector<std::string> terms;
                terms.reserve(e->args().size());
                for (const auto& a : e->args()) terms.push_back(walk(a, Prec::Add));
                std::string out = format_add(terms);
                if (parent > Prec::Add) out = "(" + out + ")";
                return out;
            }
            case TypeId::Mul: {
                // Split into numerator and denominator factors.
                std::vector<std::string> nums, dens;
                std::string sign;  // accumulated sign
                for (const auto& a : e->args()) {
                    if (a->type_id() == TypeId::Integer) {
                        const auto& z = static_cast<const Integer&>(*a);
                        if (z.value() == -1 && e->args().size() > 1) {
                            sign = "-";
                            continue;
                        }
                    }
                    if (a->type_id() == TypeId::Pow
                        && a->args()[1]->type_id() == TypeId::Integer) {
                        const auto& exp_int =
                            static_cast<const Integer&>(*a->args()[1]);
                        if (exp_int.value() < 0) {
                            // 1/base^|exp|
                            mpz_class abs_exp = -exp_int.value();
                            std::string base_s = walk(a->args()[0], Prec::Pow);
                            if (abs_exp == 1) {
                                dens.push_back(base_s);
                            } else {
                                dens.push_back(format_pow(base_s,
                                                            abs_exp.get_str()));
                            }
                            continue;
                        }
                    }
                    nums.push_back(walk(a, Prec::Mul));
                }
                std::string num_str;
                if (nums.empty()) num_str = "1";
                else {
                    for (std::size_t i = 0; i < nums.size(); ++i) {
                        if (i > 0) num_str += mul_op();
                        num_str += nums[i];
                    }
                }
                std::string out;
                if (!dens.empty()) {
                    std::string den_str;
                    for (std::size_t i = 0; i < dens.size(); ++i) {
                        if (i > 0) den_str += mul_op();
                        den_str += dens[i];
                    }
                    if (dens.size() > 1
                        || (dens.size() == 1
                             && dens[0].find_first_of("+-*/") != std::string::npos)) {
                        den_str = "(" + den_str + ")";
                    }
                    if (nums.size() > 1) num_str = "(" + num_str + ")";
                    out = num_str + div_op() + den_str;
                } else {
                    out = num_str;
                }
                if (!sign.empty()) out = sign + out;
                if (parent > Prec::Mul) out = "(" + out + ")";
                return out;
            }
            case TypeId::Pow: {
                std::string base_s = walk(e->args()[0], Prec::Pow);
                std::string exp_s = walk(e->args()[1], Prec::Pow);
                std::string out = format_pow(base_s, exp_s);
                if (parent > Prec::Pow) out = "(" + out + ")";
                return out;
            }
            case TypeId::Function: {
                const auto& fn = static_cast<const Function&>(*e);
                std::string name = func_name(fn.function_id());
                if (name.empty()) name = std::string(fn.name());
                std::string out = name + "(";
                for (std::size_t i = 0; i < e->args().size(); ++i) {
                    if (i > 0) out += ", ";
                    out += walk(e->args()[i], Prec::None);
                }
                out += ")";
                return out;
            }
            default:
                return e->str();
        }
    }
};

// ---- C printer ----
class CPrinter final : public CodePrinter {
protected:
    std::string func_name(FunctionId id) const override {
        static const std::map<FunctionId, std::string> table = {
            {FunctionId::Sin, "sin"}, {FunctionId::Cos, "cos"},
            {FunctionId::Tan, "tan"}, {FunctionId::Asin, "asin"},
            {FunctionId::Acos, "acos"}, {FunctionId::Atan, "atan"},
            {FunctionId::Atan2, "atan2"},
            {FunctionId::Sinh, "sinh"}, {FunctionId::Cosh, "cosh"},
            {FunctionId::Tanh, "tanh"},
            {FunctionId::Exp, "exp"}, {FunctionId::Log, "log"},
            {FunctionId::Abs, "fabs"}, {FunctionId::Sign, "copysign"},
            {FunctionId::Floor, "floor"}, {FunctionId::Ceiling, "ceil"},
        };
        auto it = table.find(id);
        return it != table.end() ? it->second : "";
    }
    std::string format_pow(const std::string& base,
                            const std::string& exp) const override {
        return "pow(" + base + ", " + exp + ")";
    }
    std::string format_pi() const override { return "M_PI"; }
    std::string format_e() const override { return "M_E"; }
    std::string format_imaginary_unit() const override { return "I"; }
    std::string format_rational(const std::string& num,
                                 const std::string& den) const override {
        // Force floating-point division by appending .0 to one operand.
        return num + ".0/" + den + ".0";
    }
};

// ---- C++ printer ----
class CXXPrinter final : public CodePrinter {
protected:
    std::string func_name(FunctionId id) const override {
        static const std::map<FunctionId, std::string> table = {
            {FunctionId::Sin, "std::sin"}, {FunctionId::Cos, "std::cos"},
            {FunctionId::Tan, "std::tan"}, {FunctionId::Asin, "std::asin"},
            {FunctionId::Acos, "std::acos"}, {FunctionId::Atan, "std::atan"},
            {FunctionId::Atan2, "std::atan2"},
            {FunctionId::Sinh, "std::sinh"}, {FunctionId::Cosh, "std::cosh"},
            {FunctionId::Tanh, "std::tanh"},
            {FunctionId::Exp, "std::exp"}, {FunctionId::Log, "std::log"},
            {FunctionId::Abs, "std::abs"}, {FunctionId::Sign, "std::copysign"},
            {FunctionId::Floor, "std::floor"}, {FunctionId::Ceiling, "std::ceil"},
        };
        auto it = table.find(id);
        return it != table.end() ? it->second : "";
    }
    std::string format_pow(const std::string& base,
                            const std::string& exp) const override {
        return "std::pow(" + base + ", " + exp + ")";
    }
    std::string format_pi() const override {
        return "std::numbers::pi_v<double>";
    }
    std::string format_e() const override {
        return "std::numbers::e_v<double>";
    }
    std::string format_imaginary_unit() const override {
        return "std::complex<double>(0, 1)";
    }
    std::string format_rational(const std::string& num,
                                 const std::string& den) const override {
        return num + ".0/" + den + ".0";
    }
};

// ---- Fortran printer ----
class FPrinter final : public CodePrinter {
protected:
    std::string func_name(FunctionId id) const override {
        static const std::map<FunctionId, std::string> table = {
            {FunctionId::Sin, "sin"}, {FunctionId::Cos, "cos"},
            {FunctionId::Tan, "tan"}, {FunctionId::Asin, "asin"},
            {FunctionId::Acos, "acos"}, {FunctionId::Atan, "atan"},
            {FunctionId::Atan2, "atan2"},
            {FunctionId::Sinh, "sinh"}, {FunctionId::Cosh, "cosh"},
            {FunctionId::Tanh, "tanh"},
            {FunctionId::Exp, "exp"}, {FunctionId::Log, "log"},
            {FunctionId::Abs, "abs"}, {FunctionId::Sign, "sign"},
            {FunctionId::Floor, "floor"}, {FunctionId::Ceiling, "ceiling"},
        };
        auto it = table.find(id);
        return it != table.end() ? it->second : "";
    }
    std::string format_pow(const std::string& base,
                            const std::string& exp) const override {
        return base + "**" + exp;
    }
    std::string format_pi() const override { return "3.14159265358979d0"; }
    std::string format_e() const override { return "2.71828182845905d0"; }
    std::string format_imaginary_unit() const override { return "(0,1)"; }
    std::string format_integer(const std::string& s) const override {
        // Append "d0" to keep double-precision context in Fortran 90+.
        return s + ".0d0";
    }
    std::string format_rational(const std::string& num,
                                 const std::string& den) const override {
        return num + ".0d0/" + den + ".0d0";
    }
};

// ---- Octave / MATLAB printer ----
class OctavePrinter final : public CodePrinter {
protected:
    std::string func_name(FunctionId id) const override {
        static const std::map<FunctionId, std::string> table = {
            {FunctionId::Sin, "sin"}, {FunctionId::Cos, "cos"},
            {FunctionId::Tan, "tan"}, {FunctionId::Asin, "asin"},
            {FunctionId::Acos, "acos"}, {FunctionId::Atan, "atan"},
            {FunctionId::Atan2, "atan2"},
            {FunctionId::Sinh, "sinh"}, {FunctionId::Cosh, "cosh"},
            {FunctionId::Tanh, "tanh"},
            {FunctionId::Exp, "exp"}, {FunctionId::Log, "log"},
            {FunctionId::Abs, "abs"}, {FunctionId::Sign, "sign"},
            {FunctionId::Floor, "floor"}, {FunctionId::Ceiling, "ceil"},
        };
        auto it = table.find(id);
        return it != table.end() ? it->second : "";
    }
    std::string format_pow(const std::string& base,
                            const std::string& exp) const override {
        return base + ".^" + exp;
    }
    std::string format_pi() const override { return "pi"; }
    std::string format_e() const override { return "exp(1)"; }
    std::string format_imaginary_unit() const override { return "1i"; }
    std::string mul_op() const override { return ".*"; }
    std::string div_op() const override { return "./"; }
};

// ---- LaTeX printer ----
class LatexPrinter final : public CodePrinter {
protected:
    std::string func_name(FunctionId id) const override {
        static const std::map<FunctionId, std::string> table = {
            {FunctionId::Sin, "\\sin"}, {FunctionId::Cos, "\\cos"},
            {FunctionId::Tan, "\\tan"}, {FunctionId::Asin, "\\arcsin"},
            {FunctionId::Acos, "\\arccos"}, {FunctionId::Atan, "\\arctan"},
            {FunctionId::Sinh, "\\sinh"}, {FunctionId::Cosh, "\\cosh"},
            {FunctionId::Tanh, "\\tanh"},
            {FunctionId::Exp, "\\exp"}, {FunctionId::Log, "\\log"},
            {FunctionId::Floor, "\\lfloor"}, {FunctionId::Ceiling, "\\lceil"},
        };
        auto it = table.find(id);
        return it != table.end() ? it->second : "";
    }
    std::string format_pow(const std::string& base,
                            const std::string& exp) const override {
        // sqrt-special-case via 1/2 detection done at higher level (Pow case
        // in walk is generic). Emit base^{exp}.
        return base + "^{" + exp + "}";
    }
    std::string format_pi() const override { return "\\pi"; }
    std::string format_e() const override { return "e"; }
    std::string format_imaginary_unit() const override { return "i"; }
    std::string format_rational(const std::string& num,
                                 const std::string& den) const override {
        return "\\frac{" + num + "}{" + den + "}";
    }
    std::string mul_op() const override { return " \\cdot "; }
    std::string format_symbol(const std::string& name) const override {
        // Map common Greek letter names to TeX commands.
        static const std::map<std::string, std::string> greek = {
            {"alpha", "\\alpha"}, {"beta", "\\beta"}, {"gamma", "\\gamma"},
            {"delta", "\\delta"}, {"epsilon", "\\epsilon"},
            {"zeta", "\\zeta"}, {"eta", "\\eta"}, {"theta", "\\theta"},
            {"iota", "\\iota"}, {"kappa", "\\kappa"}, {"lambda", "\\lambda"},
            {"mu", "\\mu"}, {"nu", "\\nu"}, {"xi", "\\xi"}, {"pi", "\\pi"},
            {"rho", "\\rho"}, {"sigma", "\\sigma"}, {"tau", "\\tau"},
            {"phi", "\\phi"}, {"chi", "\\chi"}, {"psi", "\\psi"},
            {"omega", "\\omega"},
        };
        auto it = greek.find(name);
        return it != greek.end() ? it->second : name;
    }

    // Override Pow to emit \sqrt for ^(1/2).
    std::string walk(const Expr& e, Prec parent) override {
        if (e && e->type_id() == TypeId::Pow) {
            const auto& exp_e = e->args()[1];
            if (exp_e->type_id() == TypeId::Rational) {
                const auto& q = static_cast<const Rational&>(*exp_e);
                if (q.numerator() == 1 && q.denominator() == 2) {
                    return "\\sqrt{" + walk(e->args()[0], Prec::None) + "}";
                }
            }
        }
        return CodePrinter::walk(e, parent);
    }
};

// ---- ASCII pretty printer (single-line for now) ----
class PrettyPrinter final : public CodePrinter {
protected:
    std::string func_name(FunctionId id) const override {
        static const std::map<FunctionId, std::string> table = {
            {FunctionId::Sin, "sin"}, {FunctionId::Cos, "cos"},
            {FunctionId::Tan, "tan"}, {FunctionId::Asin, "asin"},
            {FunctionId::Acos, "acos"}, {FunctionId::Atan, "atan"},
            {FunctionId::Sinh, "sinh"}, {FunctionId::Cosh, "cosh"},
            {FunctionId::Tanh, "tanh"},
            {FunctionId::Exp, "exp"}, {FunctionId::Log, "log"},
            {FunctionId::Abs, "abs"},
            {FunctionId::Floor, "floor"}, {FunctionId::Ceiling, "ceiling"},
        };
        auto it = table.find(id);
        return it != table.end() ? it->second : "";
    }
    std::string format_pow(const std::string& base,
                            const std::string& exp) const override {
        return base + "**" + exp;
    }
    std::string format_pi() const override { return "pi"; }
    std::string format_e() const override { return "E"; }
    std::string format_imaginary_unit() const override { return "I"; }
};

}  // namespace

// --- Public entry points ---

std::string ccode(const Expr& e) { CPrinter p; return p.print(e); }
std::string cxxcode(const Expr& e) { CXXPrinter p; return p.print(e); }
std::string fcode(const Expr& e) { FPrinter p; return p.print(e); }
std::string latex(const Expr& e) { LatexPrinter p; return p.print(e); }
std::string octave_code(const Expr& e) { OctavePrinter p; return p.print(e); }
std::string pretty(const Expr& e) { PrettyPrinter p; return p.print(e); }

// --- Function emission ---

namespace {

std::string symbol_name(const Expr& s) {
    return s->str();
}

}  // namespace

std::string c_function(const std::string& name, const Expr& body,
                        const std::vector<Expr>& args) {
    std::ostringstream os;
    os << "double " << name << "(";
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i > 0) os << ", ";
        os << "double " << symbol_name(args[i]);
    }
    os << ") {\n    return " << ccode(body) << ";\n}";
    return os.str();
}

std::string cxx_function(const std::string& name, const Expr& body,
                          const std::vector<Expr>& args) {
    std::ostringstream os;
    os << "double " << name << "(";
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i > 0) os << ", ";
        os << "double " << symbol_name(args[i]);
    }
    os << ") {\n    return " << cxxcode(body) << ";\n}";
    return os.str();
}

std::string fortran_function(const std::string& name, const Expr& body,
                               const std::vector<Expr>& args) {
    std::ostringstream os;
    os << "function " << name << "(";
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i > 0) os << ", ";
        os << symbol_name(args[i]);
    }
    os << ")\n    implicit none\n    real(8) :: " << name;
    for (const auto& a : args) {
        os << "\n    real(8), intent(in) :: " << symbol_name(a);
    }
    os << "\n    " << name << " = " << fcode(body) << "\nend function "
       << name;
    return os.str();
}

std::string octave_function(const std::string& name, const Expr& body,
                              const std::vector<Expr>& args) {
    std::ostringstream os;
    os << "function y = " << name << "(";
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i > 0) os << ", ";
        os << symbol_name(args[i]);
    }
    os << ")\n    y = " << octave_code(body) << ";\nend";
    return os.str();
}

}  // namespace sympp::printing
