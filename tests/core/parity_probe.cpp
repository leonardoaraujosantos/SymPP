// Assumption parity probe (hidden test, tag [.parity]).
//
// Prints SymPP's full predicate table for a fixed battery of expressions, one
// line per (label, predicate): "<label>\t<predicate>\t<T|F|?>". The companion
// script tools/parity_sympy.py builds the *identical* battery in SymPy and emits
// the same format; `diff` of the two outputs is the divergence report.
//
// Run: ./sympp_tests "[.parity]"  (writes to stdout)

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/miscellaneous.hpp>

using namespace sympp;

namespace {

struct Pred {
    const char* name;
    AssumptionKey key;
};

const std::vector<Pred>& predicates() {
    static const std::vector<Pred> p = {
        {"real", AssumptionKey::Real},
        {"rational", AssumptionKey::Rational},
        {"integer", AssumptionKey::Integer},
        {"positive", AssumptionKey::Positive},
        {"negative", AssumptionKey::Negative},
        {"zero", AssumptionKey::Zero},
        {"nonzero", AssumptionKey::Nonzero},
        {"nonnegative", AssumptionKey::Nonnegative},
        {"nonpositive", AssumptionKey::Nonpositive},
        {"finite", AssumptionKey::Finite},
        {"even", AssumptionKey::Even},
        {"odd", AssumptionKey::Odd},
        {"complex", AssumptionKey::Complex},
        {"imaginary", AssumptionKey::Imaginary},
        {"prime", AssumptionKey::Prime},
        {"composite", AssumptionKey::Composite},
        {"irrational", AssumptionKey::Irrational},
        {"algebraic", AssumptionKey::Algebraic},
        {"transcendental", AssumptionKey::Transcendental},
        {"extended_real", AssumptionKey::ExtendedReal},
        {"infinite", AssumptionKey::Infinite},
        {"extended_positive", AssumptionKey::ExtendedPositive},
        {"extended_negative", AssumptionKey::ExtendedNegative},
        {"extended_nonnegative", AssumptionKey::ExtendedNonnegative},
        {"extended_nonpositive", AssumptionKey::ExtendedNonpositive},
        {"hermitian", AssumptionKey::Hermitian},
        {"antihermitian", AssumptionKey::Antihermitian},
        {"commutative", AssumptionKey::Commutative},
    };
    return p;
}

char tri(std::optional<bool> v) { return !v.has_value() ? '?' : (*v ? 'T' : 'F'); }

std::vector<std::pair<std::string, Expr>> battery() {
    auto S_ = [](AssumptionMask m) { return m; };
    std::vector<std::pair<std::string, Expr>> b;
    auto add = [&](const std::string& label, const Expr& e) { b.emplace_back(label, e); };

    // Symbols carrying a single declared fact.
    add("sym_generic", symbol("g"));
    add("sym_real", symbol("s", S_(AssumptionMask{}.set_real(true))));
    add("sym_positive", symbol("s", S_(AssumptionMask{}.set_positive(true))));
    add("sym_negative", symbol("s", S_(AssumptionMask{}.set_negative(true))));
    add("sym_nonnegative", symbol("s", S_(AssumptionMask{}.set_nonnegative(true))));
    add("sym_nonpositive", symbol("s", S_(AssumptionMask{}.set_nonpositive(true))));
    add("sym_nonzero", symbol("s", S_(AssumptionMask{}.set_nonzero(true))));
    add("sym_integer", symbol("s", S_(AssumptionMask{}.set_integer(true))));
    add("sym_rational", symbol("s", S_(AssumptionMask{}.set_rational(true))));
    add("sym_even", symbol("s", S_(AssumptionMask{}.set_even(true))));
    add("sym_odd", symbol("s", S_(AssumptionMask{}.set_odd(true))));
    add("sym_prime", symbol("s", S_(AssumptionMask{}.set_prime(true))));
    add("sym_composite", symbol("s", S_(AssumptionMask{}.set_composite(true))));
    add("sym_irrational", symbol("s", S_(AssumptionMask{}.set_irrational(true))));
    add("sym_imaginary", symbol("s", S_(AssumptionMask{}.set_imaginary(true))));
    add("sym_complex", symbol("s", S_(AssumptionMask{}.set_complex(true))));
    add("sym_finite", symbol("s", S_(AssumptionMask{}.set_finite(true))));
    add("sym_extended_real", symbol("s", S_(AssumptionMask{}.set_extended_real(true))));
    add("sym_hermitian", symbol("s", S_(AssumptionMask{}.set_hermitian(true))));
    add("sym_antihermitian", symbol("s", S_(AssumptionMask{}.set_antihermitian(true))));
    add("sym_ext_positive", symbol("s", S_(AssumptionMask{}.set_extended_positive(true))));
    add("sym_ext_negative", symbol("s", S_(AssumptionMask{}.set_extended_negative(true))));

    // Symbols carrying a pair of facts.
    add("sym_real_nonzero",
        symbol("s", S_(AssumptionMask{}.set_real(true).set_nonzero(true))));
    add("sym_integer_positive",
        symbol("s", S_(AssumptionMask{}.set_integer(true).set_positive(true))));
    add("sym_integer_nonnegative",
        symbol("s", S_(AssumptionMask{}.set_integer(true).set_nonnegative(true))));
    add("sym_ext_positive_finite",
        symbol("s", S_(AssumptionMask{}.set_extended_positive(true).set_finite(true))));

    // Numbers.
    add("int_0", integer(0));
    add("int_1", integer(1));
    add("int_2", integer(2));
    add("int_3", integer(3));
    add("int_4", integer(4));
    add("int_neg1", integer(-1));
    add("int_neg2", integer(-2));
    add("rat_half", rational(1, 2));
    add("rat_neg3_2", rational(-3, 2));
    add("oo", S::Infinity());
    add("neg_oo", S::NegativeInfinity());
    add("zoo", S::ComplexInfinity());
    add("I", S::I());
    add("pi", S::Pi());
    add("E", S::E());

    // Compounds.
    auto g = symbol("g");
    auto p = symbol("p", S_(AssumptionMask{}.set_positive(true)));
    auto q = symbol("q", S_(AssumptionMask{}.set_positive(true)));
    auto nn = symbol("n", S_(AssumptionMask{}.set_integer(true)));
    auto r = symbol("r", S_(AssumptionMask{}.set_real(true)));
    add("abs_generic", abs(g));
    add("neg_positive", -p);
    add("p_squared", pow(p, integer(2)));
    add("g_squared", pow(g, integer(2)));
    add("g_squared_plus1", pow(g, integer(2)) + integer(1));
    add("p_plus_q", p + q);
    add("p_times_q", p * q);
    add("two_n_plus_1", integer(2) * nn + integer(1));
    add("i_times_real", S::I() * r);
    return b;
}

}  // namespace

TEST_CASE("assumption parity probe", "[.parity]") {
    for (const auto& [label, e] : battery()) {
        for (const auto& pr : predicates()) {
            std::cout << label << '\t' << pr.name << '\t' << tri(ask(e, pr.key)) << '\n';
        }
    }
    std::cout.flush();
    SUCCEED();
}
