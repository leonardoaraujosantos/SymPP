#include <sympp/parsing/latex_parser.hpp>

#include <cctype>
#include <stdexcept>
#include <string>
#include <vector>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>

namespace sympp {

namespace {

enum class Tok { Num, Ident, Cmd, Plus, Minus, Times, Slash, Caret,
                 LBrace, RBrace, LParen, RParen, End };
struct Token { Tok type; std::string text; };

// --- Tokenizer ---------------------------------------------------------------
std::vector<Token> tokenize(const std::string& s) {
    std::vector<Token> out;
    std::size_t i = 0, n = s.size();
    while (i < n) {
        char c = s[i];
        if (std::isspace(static_cast<unsigned char>(c))) { ++i; continue; }
        if (c == '\\') {
            ++i;
            if (i < n && !std::isalpha(static_cast<unsigned char>(s[i]))) {
                // spacing command (\, \! \; …) or \\ — skip the single char.
                ++i;
                continue;
            }
            std::string name;
            while (i < n && std::isalpha(static_cast<unsigned char>(s[i]))) name += s[i++];
            if (name == "cdot" || name == "times") out.push_back({Tok::Times, name});
            else if (name == "left" || name == "right") { /* skip */ }
            else out.push_back({Tok::Cmd, name});
            continue;
        }
        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
            std::string num;
            while (i < n && (std::isdigit(static_cast<unsigned char>(s[i])) || s[i] == '.'))
                num += s[i++];
            out.push_back({Tok::Num, num});
            continue;
        }
        if (std::isalpha(static_cast<unsigned char>(c))) {
            out.push_back({Tok::Ident, std::string(1, c)});  // each letter is a symbol
            ++i;
            continue;
        }
        switch (c) {
            case '+': out.push_back({Tok::Plus, "+"}); break;
            case '-': out.push_back({Tok::Minus, "-"}); break;
            case '*': out.push_back({Tok::Times, "*"}); break;
            case '/': out.push_back({Tok::Slash, "/"}); break;
            case '^': out.push_back({Tok::Caret, "^"}); break;
            case '{': out.push_back({Tok::LBrace, "{"}); break;
            case '}': out.push_back({Tok::RBrace, "}"}); break;
            case '(': out.push_back({Tok::LParen, "("}); break;
            case ')': out.push_back({Tok::RParen, ")"}); break;
            case '_': /* subscript: fold into the identifier name below */ break;
            default: throw std::runtime_error(std::string("parse_latex: bad char '") + c + "'");
        }
        ++i;
    }
    out.push_back({Tok::End, ""});
    return out;
}

// --- Recursive-descent parser ------------------------------------------------
class Parser {
public:
    explicit Parser(std::vector<Token> toks) : t_(std::move(toks)) {}
    Expr parse() {
        Expr e = expr();
        if (peek().type != Tok::End) throw std::runtime_error("parse_latex: trailing input");
        return e;
    }

private:
    [[nodiscard]] const Token& peek() const { return t_[pos_]; }
    Token next() { return t_[pos_++]; }
    bool accept(Tok k) { if (peek().type == k) { ++pos_; return true; } return false; }
    void expect(Tok k, const char* what) {
        if (!accept(k)) throw std::runtime_error(std::string("parse_latex: expected ") + what);
    }

    [[nodiscard]] bool starts_factor() const {
        switch (peek().type) {
            case Tok::Num: case Tok::Ident: case Tok::Cmd:
            case Tok::LParen: case Tok::LBrace: case Tok::Minus:
                return true;
            default: return false;
        }
    }

    Expr expr() {
        Expr e = term();
        while (peek().type == Tok::Plus || peek().type == Tok::Minus) {
            bool minus = next().type == Tok::Minus;
            Expr r = term();
            e = add(e, minus ? mul(integer(-1), r) : r);
        }
        return e;
    }

    Expr term() {
        Expr e = factor();
        for (;;) {
            if (accept(Tok::Times)) {
                e = mul(e, factor());
            } else if (accept(Tok::Slash)) {
                e = mul(e, pow(factor(), integer(-1)));
            } else if (peek().type != Tok::Plus && peek().type != Tok::Minus &&
                       starts_factor()) {
                e = mul(e, factor());  // implicit multiplication
            } else {
                break;
            }
        }
        return e;
    }

    Expr factor() {
        if (accept(Tok::Minus)) return mul(integer(-1), factor());
        Expr base = atom();
        if (accept(Tok::Caret)) base = pow(base, exponent());
        return base;
    }

    Expr exponent() {
        if (peek().type == Tok::LBrace) return braced();
        // single token exponent (e.g. x^2)
        return factor();
    }

    Expr braced() {
        expect(Tok::LBrace, "{");
        Expr e = expr();
        expect(Tok::RBrace, "}");
        return e;
    }

    Expr atom() {
        const Token& tk = peek();
        switch (tk.type) {
            case Tok::Num: {
                std::string s = next().text;
                if (s.find('.') == std::string::npos) return integer(std::stoll(s));
                throw std::runtime_error("parse_latex: decimals not supported");
            }
            case Tok::Ident: return symbol(next().text);
            case Tok::LParen: { next(); Expr e = expr(); expect(Tok::RParen, ")"); return e; }
            case Tok::LBrace: return braced();
            case Tok::Cmd: return command();
            default: throw std::runtime_error("parse_latex: unexpected token");
        }
    }

    // Argument of a function command: \sin(x), \sin{x} or \sin x.
    Expr func_arg() {
        if (peek().type == Tok::LParen) { next(); Expr e = expr(); expect(Tok::RParen, ")"); return e; }
        if (peek().type == Tok::LBrace) return braced();
        return atom();
    }

    Expr command() {
        std::string name = next().text;
        if (name == "pi") return S::Pi();
        if (name == "frac") {
            Expr num = braced();  // sequenced: numerator first, then denominator
            Expr den = braced();
            return mul(num, pow(den, integer(-1)));
        }
        if (name == "sqrt") {
            if (accept(Tok::LParen)) {  // \sqrt[...] uses '[' but we only do {}
                throw std::runtime_error("parse_latex: \\sqrt[ not supported");
            }
            return pow(braced(), pow(integer(2), integer(-1)));
        }
        if (name == "sin") return sin(func_arg());
        if (name == "cos") return cos(func_arg());
        if (name == "tan") return tan(func_arg());
        if (name == "cot") return cot(func_arg());
        if (name == "sec") return sec(func_arg());
        if (name == "csc") return csc(func_arg());
        if (name == "arcsin") return asin(func_arg());
        if (name == "arccos") return acos(func_arg());
        if (name == "arctan") return atan(func_arg());
        if (name == "sinh") return sinh(func_arg());
        if (name == "cosh") return cosh(func_arg());
        if (name == "tanh") return tanh(func_arg());
        if (name == "coth") return coth(func_arg());
        if (name == "exp") return exp(func_arg());
        if (name == "log" || name == "ln") return log(func_arg());
        // Otherwise treat as a symbol named after the command (Greek letters etc.).
        return symbol(name);
    }

    std::vector<Token> t_;
    std::size_t pos_ = 0;
};

}  // namespace

Expr parse_latex(const std::string& latex) {
    return Parser(tokenize(latex)).parse();
}

}  // namespace sympp
