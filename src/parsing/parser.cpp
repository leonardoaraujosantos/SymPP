#include <sympp/parsing/parser.hpp>

#include <cctype>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/integers.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/special.hpp>
#include <sympp/functions/trigonometric.hpp>

namespace sympp::parsing {

namespace {

// ---------------------------------------------------------------------------
// Lexer
// ---------------------------------------------------------------------------

enum class TokKind {
    End,
    Number,        // 42, 3.14, 1e-3
    Identifier,    // x, alpha, sin
    Plus, Minus,
    Star, Slash,
    StarStar,      // **
    Caret,         // ^ (alias for **)
    LParen, RParen,
    Comma,
};

struct Token {
    TokKind kind;
    std::string text;
    std::size_t pos;
};

class Lexer {
public:
    explicit Lexer(std::string_view source) : src_(source) {}

    [[nodiscard]] Token next() {
        skip_ws();
        if (i_ >= src_.size()) return {TokKind::End, "", i_};
        std::size_t start = i_;
        char c = src_[i_];
        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
            return number(start);
        }
        if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
            return identifier(start);
        }
        ++i_;
        switch (c) {
            case '+': return {TokKind::Plus, "+", start};
            case '-': return {TokKind::Minus, "-", start};
            case '*':
                if (i_ < src_.size() && src_[i_] == '*') {
                    ++i_;
                    return {TokKind::StarStar, "**", start};
                }
                return {TokKind::Star, "*", start};
            case '/': return {TokKind::Slash, "/", start};
            case '^': return {TokKind::Caret, "^", start};
            case '(': return {TokKind::LParen, "(", start};
            case ')': return {TokKind::RParen, ")", start};
            case ',': return {TokKind::Comma, ",", start};
        }
        throw ParseError(std::string("unexpected character '") + c + "'",
                         start);
    }

private:
    void skip_ws() {
        while (i_ < src_.size()
               && std::isspace(static_cast<unsigned char>(src_[i_]))) {
            ++i_;
        }
    }

    [[nodiscard]] Token number(std::size_t start) {
        bool seen_dot = false;
        bool seen_exp = false;
        while (i_ < src_.size()) {
            char c = src_[i_];
            if (std::isdigit(static_cast<unsigned char>(c))) {
                ++i_;
            } else if (c == '.' && !seen_dot && !seen_exp) {
                seen_dot = true;
                ++i_;
            } else if ((c == 'e' || c == 'E') && !seen_exp) {
                seen_exp = true;
                ++i_;
                if (i_ < src_.size() && (src_[i_] == '+' || src_[i_] == '-')) {
                    ++i_;
                }
            } else {
                break;
            }
        }
        return {TokKind::Number, std::string(src_.substr(start, i_ - start)),
                start};
    }

    [[nodiscard]] Token identifier(std::size_t start) {
        while (i_ < src_.size()
               && (std::isalnum(static_cast<unsigned char>(src_[i_]))
                   || src_[i_] == '_')) {
            ++i_;
        }
        return {TokKind::Identifier,
                std::string(src_.substr(start, i_ - start)), start};
    }

    std::string_view src_;
    std::size_t i_ = 0;
};

// ---------------------------------------------------------------------------
// Function registry: identifier → factory
// ---------------------------------------------------------------------------

using OneArgFn = Expr (*)(const Expr&);
using TwoArgFn = Expr (*)(const Expr&, const Expr&);

[[nodiscard]] const std::map<std::string, OneArgFn>& one_arg_funcs() {
    static const std::map<std::string, OneArgFn> table = {
        {"sin", &sin}, {"cos", &cos}, {"tan", &tan},
        {"asin", &asin}, {"acos", &acos}, {"atan", &atan},
        {"sinh", &sinh}, {"cosh", &cosh}, {"tanh", &tanh},
        {"asinh", &asinh}, {"acosh", &acosh}, {"atanh", &atanh},
        {"exp", &exp}, {"log", &log},
        {"sqrt", &sqrt}, {"abs", &abs}, {"sign", &sign},
        {"floor", &floor}, {"ceiling", &ceiling},
        {"factorial", &factorial},
        {"gamma", &gamma}, {"loggamma", &loggamma},
        {"erf", &erf}, {"erfc", &erfc},
        {"heaviside", &heaviside}, {"dirac_delta", &dirac_delta},
        {"re", &re}, {"im", &im}, {"conjugate", &conjugate},
    };
    return table;
}

[[nodiscard]] const std::map<std::string, TwoArgFn>& two_arg_funcs() {
    static const std::map<std::string, TwoArgFn> table = {
        {"atan2", &atan2},
        {"binomial", &binomial},
    };
    return table;
}

[[nodiscard]] std::optional<Expr> recognized_constant(const std::string& name) {
    if (name == "pi" || name == "Pi") return S::Pi();
    if (name == "E") return S::E();
    if (name == "I") return S::I();
    if (name == "True") return S::True();
    if (name == "False") return S::False();
    if (name == "EulerGamma") return S::EulerGamma();
    if (name == "Catalan") return S::Catalan();
    return std::nullopt;
}

// ---------------------------------------------------------------------------
// Number parsing — pick Integer / Rational / Float by syntactic shape.
// ---------------------------------------------------------------------------

[[nodiscard]] Expr parse_number(const std::string& s, std::size_t pos) {
    bool has_dot = s.find('.') != std::string::npos;
    bool has_exp = s.find_first_of("eE") != std::string::npos;
    if (!has_dot && !has_exp) {
        return integer(s);
    }
    return make<Float>(s, 15);
    (void)pos;
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

class Parser {
public:
    explicit Parser(std::string_view source)
        : lex_(source), src_(source) {
        advance();
    }

    [[nodiscard]] Expr parse_top() {
        Expr e = parse_add();
        if (cur_.kind != TokKind::End) {
            throw ParseError("unexpected trailing token: '" + cur_.text + "'",
                              cur_.pos);
        }
        return e;
    }

private:
    void advance() { cur_ = lex_.next(); }

    [[nodiscard]] Expr parse_add() {
        Expr lhs = parse_mul();
        while (cur_.kind == TokKind::Plus || cur_.kind == TokKind::Minus) {
            bool minus = cur_.kind == TokKind::Minus;
            advance();
            Expr rhs = parse_mul();
            lhs = minus ? lhs - rhs : lhs + rhs;
        }
        return lhs;
    }

    [[nodiscard]] Expr parse_mul() {
        Expr lhs = parse_unary();
        while (cur_.kind == TokKind::Star || cur_.kind == TokKind::Slash) {
            bool divide = cur_.kind == TokKind::Slash;
            advance();
            Expr rhs = parse_unary();
            lhs = divide ? lhs / rhs : lhs * rhs;
        }
        return lhs;
    }

    // Power binds tighter than unary minus (Python/SymPy convention:
    // -x**2 == -(x**2), not (-x)**2). So unary applies *after* pow_e.
    [[nodiscard]] Expr parse_unary() {
        if (cur_.kind == TokKind::Plus) {
            advance();
            return parse_unary();
        }
        if (cur_.kind == TokKind::Minus) {
            advance();
            return mul(S::NegativeOne(), parse_unary());
        }
        return parse_pow();
    }

    [[nodiscard]] Expr parse_pow() {
        Expr base = parse_atom();
        if (cur_.kind == TokKind::StarStar || cur_.kind == TokKind::Caret) {
            advance();
            // Right-associative: parse_unary on the rhs so
            // x**-2 parses as x ** (-2).
            Expr exp = parse_unary();
            return pow(base, exp);
        }
        return base;
    }

    [[nodiscard]] Expr parse_atom() {
        if (cur_.kind == TokKind::Number) {
            Expr v = parse_number(cur_.text, cur_.pos);
            advance();
            return v;
        }
        if (cur_.kind == TokKind::LParen) {
            advance();
            Expr e = parse_add();
            expect(TokKind::RParen, "expected ')'");
            return e;
        }
        if (cur_.kind == TokKind::Identifier) {
            std::string name = cur_.text;
            std::size_t name_pos = cur_.pos;
            advance();
            if (cur_.kind == TokKind::LParen) {
                advance();
                std::vector<Expr> args;
                if (cur_.kind != TokKind::RParen) {
                    args.push_back(parse_add());
                    while (cur_.kind == TokKind::Comma) {
                        advance();
                        args.push_back(parse_add());
                    }
                }
                expect(TokKind::RParen, "expected ')'");
                return apply_function(name, args, name_pos);
            }
            // Plain identifier — constant or symbol.
            if (auto c = recognized_constant(name); c) return *c;
            return symbol(name);
        }
        throw ParseError("unexpected token: '" + cur_.text + "'", cur_.pos);
    }

    [[nodiscard]] Expr apply_function(const std::string& name,
                                       const std::vector<Expr>& args,
                                       std::size_t pos) {
        if (args.size() == 1) {
            const auto& t = one_arg_funcs();
            auto it = t.find(name);
            if (it != t.end()) return it->second(args[0]);
        }
        if (args.size() == 2) {
            const auto& t = two_arg_funcs();
            auto it = t.find(name);
            if (it != t.end()) return it->second(args[0], args[1]);
        }
        // Unknown: treat as user-defined function.
        (void)pos;
        return function_symbol(name)(args);
    }

    void expect(TokKind k, const std::string& msg) {
        if (cur_.kind != k) throw ParseError(msg, cur_.pos);
        advance();
    }

    Lexer lex_;
    std::string_view src_;
    Token cur_{TokKind::End, "", 0};
};

}  // namespace

Expr parse(std::string_view source) {
    Parser p(source);
    return p.parse_top();
}

}  // namespace sympp::parsing
