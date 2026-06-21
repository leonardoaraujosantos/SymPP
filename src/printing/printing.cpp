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
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
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
            {FunctionId::Cot, "\\cot"}, {FunctionId::Sec, "\\sec"},
            {FunctionId::Csc, "\\csc"},
            {FunctionId::Acot, "\\operatorname{acot}"},
            {FunctionId::Asec, "\\operatorname{asec}"},
            {FunctionId::Acsc, "\\operatorname{acsc}"},
            {FunctionId::Sinh, "\\sinh"}, {FunctionId::Cosh, "\\cosh"},
            {FunctionId::Tanh, "\\tanh"}, {FunctionId::Coth, "\\coth"},
            {FunctionId::Sech, "\\operatorname{sech}"},
            {FunctionId::Csch, "\\operatorname{csch}"},
            {FunctionId::Acoth, "\\operatorname{acoth}"},
            {FunctionId::Asech, "\\operatorname{asech}"},
            {FunctionId::Acsch, "\\operatorname{acsch}"},
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


}  // namespace

// ---------------------------------------------------------------------------
// 2D pretty printer — multi-line ASCII layout (stacked fractions, superscript
// powers, radical signs), matching SymPy's `pretty(use_unicode=False)` on the
// common constructs.
// ---------------------------------------------------------------------------
namespace {

struct Block {
    std::vector<std::string> lines;
    std::size_t baseline = 0;
    [[nodiscard]] std::size_t width() const { return lines.empty() ? 0 : lines[0].size(); }
    [[nodiscard]] std::size_t height() const { return lines.size(); }
};

[[nodiscard]] std::string sp(std::size_t n) { return std::string(n, ' '); }
[[nodiscard]] Block atom2d(const std::string& s) { return {{s}, 0}; }
[[nodiscard]] std::string center(const std::string& s, std::size_t w) {
    if (s.size() >= w) return s;
    std::size_t total = w - s.size(), left = total / 2;
    return sp(left) + s + sp(total - left);
}

[[nodiscard]] Block hcat(const Block& a, const Block& b) {
    std::size_t above = std::max(a.baseline, b.baseline);
    std::size_t below = std::max(a.height() - 1 - a.baseline, b.height() - 1 - b.baseline);
    std::size_t h = above + below + 1, wa = a.width(), wb = b.width();
    std::size_t topA = above - a.baseline, topB = above - b.baseline;
    Block r;
    r.baseline = above;
    r.lines.resize(h);
    for (std::size_t i = 0; i < h; ++i) {
        std::string la = (i >= topA && i < topA + a.height()) ? a.lines[i - topA] : sp(wa);
        std::string lb = (i >= topB && i < topB + b.height()) ? b.lines[i - topB] : sp(wb);
        if (la.size() < wa) la += sp(wa - la.size());
        if (lb.size() < wb) lb += sp(wb - lb.size());
        r.lines[i] = la + lb;
    }
    return r;
}

[[nodiscard]] Block superscript(const Block& base, const Block& exp) {
    Block r;
    r.baseline = exp.height() + base.baseline;
    r.lines.resize(exp.height() + base.height());
    for (std::size_t i = 0; i < exp.height(); ++i) {
        std::string l = exp.lines[i];
        if (l.size() < exp.width()) l += sp(exp.width() - l.size());
        r.lines[i] = sp(base.width()) + l;
    }
    for (std::size_t i = 0; i < base.height(); ++i) {
        std::string l = base.lines[i];
        if (l.size() < base.width()) l += sp(base.width() - l.size());
        r.lines[exp.height() + i] = l + sp(exp.width());
    }
    return r;
}

[[nodiscard]] Block frac2d(const Block& n, const Block& d) {
    std::size_t w = std::max(n.width(), d.width());
    Block r;
    r.baseline = n.height();
    for (const auto& l : n.lines) r.lines.push_back(center(l, w));
    r.lines.push_back(std::string(w, '-'));
    for (const auto& l : d.lines) r.lines.push_back(center(l, w));
    return r;
}

// Radical for a single-line radicand (the common case): \/ over an overline.
[[nodiscard]] Block sqrt2d(const Block& rad) {
    if (rad.height() == 1) {
        std::size_t w = rad.width();
        Block r;
        r.baseline = 1;
        r.lines.push_back(sp(2) + std::string(w + 2, '_'));
        r.lines.push_back("\\/ " + rad.lines[0] + " ");
        return r;
    }
    // Taller radicands: a vertical bar on the left + overline (approximate).
    Block r;
    std::size_t w = rad.width();
    r.lines.push_back(sp(2) + std::string(w + 1, '_'));
    for (std::size_t i = 0; i < rad.height(); ++i) {
        r.lines.push_back((i + 1 == rad.height() ? "\\/" : " |") + (" " + rad.lines[i]));
    }
    r.baseline = rad.baseline + 1;
    return r;
}

[[nodiscard]] bool is_compound(const Expr& e) {
    auto t = e->type_id();
    return t == TypeId::Add || t == TypeId::Mul || t == TypeId::Pow;
}

[[nodiscard]] Block render2d(const Expr& e);

// Wrap a block in parentheses (baseline row only — adequate for short args).
[[nodiscard]] Block parens(const Block& b) {
    return hcat(hcat(atom2d("("), b), atom2d(")"));
}

// Numeric sign test for a leading factor / term.
[[nodiscard]] bool is_negative_number(const Expr& e) {
    if (e->type_id() == TypeId::Integer) {
        return static_cast<const Integer&>(*e).value() < 0;
    }
    if (e->type_id() == TypeId::Rational) {
        return static_cast<const Rational&>(*e).value() < 0;
    }
    return false;
}

// Negate a numeric exponent.
[[nodiscard]] Expr negate_num(const Expr& e) {
    if (e->type_id() == TypeId::Integer) {
        return make<Integer>(mpz_class(-static_cast<const Integer&>(*e).value()));
    }
    const auto& q = static_cast<const Rational&>(*e).value();
    return rational(-q.get_num(), q.get_den());
}

[[nodiscard]] Block render_pow(const Expr& base, const Expr& exp) {
    // x**(1/2) → radical.
    if (exp->type_id() == TypeId::Rational) {
        const auto& q = static_cast<const Rational&>(*exp).value();
        if (q.get_num() == 1 && q.get_den() == 2) return sqrt2d(render2d(base));
    }
    // negative exponent → 1 / x**|exp|.
    if (is_negative_number(exp)) {
        Expr posexp = negate_num(exp);
        Expr denom = (posexp == integer(1)) ? base : pow(base, posexp);
        return frac2d(atom2d("1"), render2d(denom));
    }
    Block b = render2d(base);
    if (is_compound(base)) b = parens(b);
    return superscript(b, render2d(exp));
}

// Split a Mul into (sign, numerator factors, denominator factors).
struct MulParts {
    bool negative = false;
    std::vector<Expr> num, den;
};
[[nodiscard]] MulParts split_mul(const Expr& e) {
    MulParts p;
    for (const auto& f : e->args()) {
        if (f->type_id() == TypeId::Pow && is_negative_number(f->args()[1])) {
            Expr pe = negate_num(f->args()[1]);
            p.den.push_back(pe == integer(1) ? f->args()[0] : pow(f->args()[0], pe));
        } else if (f->type_id() == TypeId::Rational) {
            const auto& q = static_cast<const Rational&>(*f).value();
            mpz_class n = q.get_num();
            if (n < 0) { p.negative = !p.negative; n = -n; }
            if (n != 1) p.num.push_back(make<Integer>(n));
            p.den.push_back(make<Integer>(q.get_den()));
        } else if (f == integer(-1)) {
            p.negative = !p.negative;
        } else if (is_negative_number(f)) {
            p.negative = !p.negative;
            p.num.push_back(negate_num(f));
        } else {
            p.num.push_back(f);
        }
    }
    return p;
}

[[nodiscard]] Block product_block(const std::vector<Expr>& factors) {
    if (factors.empty()) return atom2d("1");
    // A single factor needs no '*' and (as a fraction numerator/denominator) no
    // parentheses — the bar already groups it.
    if (factors.size() == 1) return render2d(factors[0]);
    Block r = render2d(factors[0]);
    if (factors[0]->type_id() == TypeId::Add) r = parens(r);
    for (std::size_t i = 1; i < factors.size(); ++i) {
        Block f = render2d(factors[i]);
        if (factors[i]->type_id() == TypeId::Add) f = parens(f);
        r = hcat(hcat(r, atom2d("*")), f);
    }
    return r;
}

[[nodiscard]] Block mul_body(const MulParts& p) {
    return p.den.empty() ? product_block(p.num)
                         : frac2d(product_block(p.num), product_block(p.den));
}

// Sign of a term, decided from the expression (never from rendered '-', which
// would clash with a fraction bar).
[[nodiscard]] bool is_neg_term(const Expr& e) {
    if (e->type_id() == TypeId::Integer) return static_cast<const Integer&>(*e).value() < 0;
    if (e->type_id() == TypeId::Rational) return static_cast<const Rational&>(*e).value() < 0;
    if (e->type_id() == TypeId::Mul) return split_mul(e).negative;
    return false;
}

// Render the magnitude of a (possibly negative) term.
[[nodiscard]] Block render_abs(const Expr& e) {
    if (e->type_id() == TypeId::Integer || e->type_id() == TypeId::Rational) {
        return render2d(negate_num(e));
    }
    if (e->type_id() == TypeId::Mul) return mul_body(split_mul(e));
    return render2d(e);
}

Block render2d(const Expr& e) {
    switch (e->type_id()) {
        case TypeId::Integer:
        case TypeId::Symbol:
        case TypeId::NumberSymbol:
        case TypeId::ImaginaryUnit:
        case TypeId::Float:
        case TypeId::Boolean:
            return atom2d(e->str());
        case TypeId::Rational: {
            const auto& q = static_cast<const Rational&>(*e).value();
            if (q.get_den() == 1) return atom2d(e->str());
            bool neg = q.get_num() < 0;
            mpz_class n = neg ? mpz_class(-q.get_num()) : q.get_num();
            Block f = frac2d(atom2d(n.get_str()), atom2d(q.get_den().get_str()));
            return neg ? hcat(atom2d("-"), f) : f;
        }
        case TypeId::Pow:
            return render_pow(e->args()[0], e->args()[1]);
        case TypeId::Mul: {
            MulParts p = split_mul(e);
            Block body = mul_body(p);
            return p.negative ? hcat(atom2d("-"), body) : body;
        }
        case TypeId::Add: {
            const auto& args = e->args();
            Block r = render2d(args[0]);
            for (std::size_t i = 1; i < args.size(); ++i) {
                if (is_neg_term(args[i])) {
                    r = hcat(hcat(r, atom2d(" - ")), render_abs(args[i]));
                } else {
                    r = hcat(hcat(r, atom2d(" + ")), render2d(args[i]));
                }
            }
            return r;
        }
        case TypeId::Function: {
            const auto& fn = static_cast<const Function&>(*e);
            if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1) {
                return superscript(atom2d("e"), render2d(fn.args()[0]));
            }
            Block args_block;
            const auto& as = fn.args();
            for (std::size_t i = 0; i < as.size(); ++i) {
                Block a = render2d(as[i]);
                args_block = (i == 0) ? a : hcat(hcat(args_block, atom2d(", ")), a);
            }
            return hcat(hcat(atom2d(std::string(fn.name()) + "("), args_block), atom2d(")"));
        }
        default:
            return atom2d(e->str());
    }
}

[[nodiscard]] std::string render_pretty(const Expr& e) {
    Block b = render2d(e);
    std::string out;
    for (std::size_t i = 0; i < b.lines.size(); ++i) {
        if (i) out += "\n";
        out += b.lines[i];
    }
    return out;
}

}  // namespace

// --- Public entry points ---

std::string ccode(const Expr& e) { CPrinter p; return p.print(e); }
std::string cxxcode(const Expr& e) { CXXPrinter p; return p.print(e); }
std::string fcode(const Expr& e) { FPrinter p; return p.print(e); }
std::string latex(const Expr& e) { LatexPrinter p; return p.print(e); }
std::string octave_code(const Expr& e) { OctavePrinter p; return p.print(e); }
std::string pretty(const Expr& e) { return render_pretty(e); }

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
