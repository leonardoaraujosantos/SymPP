// Assumption parity regression guard.
//
// Compares SymPP's full predicate table over a fixed battery of expressions
// against a committed SymPy golden table (tools/parity_expected.tsv, produced by
// tools/parity_sympy.py from the IDENTICAL battery). Every (label, predicate)
// must match SymPy, except a small whitelist of *intentional* divergences:
//
//   * SymPP is strictly more decisive and still correct (e.g. Abs(x) is
//     real / nonnegative / complex / hermitian for a generic x);
//   * SymPP's `nonzero` means "≠ 0" rather than SymPy's "real ∧ ≠ 0";
//   * SymPP treats ±∞ as (strictly) positive / negative — required by the
//     limit engine — where SymPy ties strict sign to finiteness.
//
// A NEW divergence (anything not whitelisted) fails the test. To refresh the
// golden after an intended SymPy-facing change: tools/parity_sympy.py > the tsv.

#include <fstream>
#include <map>
#include <set>
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
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>

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
    add("sqrt2", pow(integer(2), rational(1, 2)));
    add("abs_real", abs(r));
    add("sym_transcendental", symbol("t", S_(AssumptionMask{}.set_transcendental(true))));
    add("two_pi", integer(2) * S::Pi());
    add("pi_plus_one", S::Pi() + integer(1));
    auto nsym = symbol("m", S_(AssumptionMask{}.set_negative(true)));
    add("neg_times_neg", nsym * symbol("m2", S_(AssumptionMask{}.set_negative(true))));
    add("pos_times_neg", p * nsym);
    add("pow_neg_2", pow(nsym, integer(2)));
    add("pow_real_2", pow(r, integer(2)));
    add("p_over_q", p * pow(q, integer(-1)));
    add("exp_real", exp(r));
    add("log_positive", log(p));
    add("sin_real", sin(r));
    return b;
}

// Intentional divergences: (label, predicate) pairs where SymPP deliberately
// differs from SymPy. See the file header for the three reasons.
const std::set<std::pair<std::string, std::string>>& whitelist() {
    static const std::set<std::pair<std::string, std::string>> w = {
        // Abs(x): SymPP is strictly more decisive (and correct) for generic x.
        {"abs_generic", "real"}, {"abs_generic", "nonnegative"},
        {"abs_generic", "complex"}, {"abs_generic", "hermitian"},
        // `nonzero` = "≠ 0" in SymPP (not SymPy's real ∧ ≠ 0).
        {"I", "nonzero"}, {"zoo", "nonzero"}, {"sym_imaginary", "nonzero"},
        {"oo", "nonzero"}, {"neg_oo", "nonzero"},
        {"sym_ext_positive", "nonzero"}, {"sym_ext_negative", "nonzero"},
        {"sym_transcendental", "nonzero"},
        {"sym_nonzero", "real"}, {"sym_nonzero", "complex"},
        {"sym_nonzero", "finite"}, {"sym_nonzero", "extended_real"},
        {"sym_nonzero", "hermitian"}, {"sym_nonzero", "imaginary"},
        {"sym_nonzero", "infinite"},
        // ±∞ are (strictly) signed in SymPP — load-bearing for the limit engine.
        {"oo", "positive"}, {"oo", "nonnegative"},
        {"neg_oo", "negative"}, {"neg_oo", "nonpositive"},
    };
    return w;
}

std::map<std::pair<std::string, std::string>, char> load_golden() {
    std::map<std::pair<std::string, std::string>, char> g;
    std::ifstream in(SYMPP_TEST_PARITY_DATA);
    REQUIRE(in.is_open());
    std::string label, pred, val;
    while (in >> label >> pred >> val) {
        g[{label, pred}] = val.empty() ? '?' : val[0];
    }
    return g;
}

}  // namespace

TEST_CASE("assumption parity with SymPy (no new divergences)", "[parity]") {
    const auto golden = load_golden();
    REQUIRE(golden.size() > 1000);  // sanity: the table loaded

    std::vector<std::string> unexpected;
    for (const auto& [label, e] : battery()) {
        for (const auto& pr : predicates()) {
            auto it = golden.find({label, pr.name});
            if (it == golden.end()) continue;  // predicate not in golden
            const char got = tri(ask(e, pr.key));
            const char want = it->second;
            if (got == want) continue;
            if (whitelist().count({label, pr.name})) continue;  // intentional
            unexpected.push_back(std::string(label) + '.' + pr.name + ": SymPP="
                                 + got + " SymPy=" + want);
        }
    }
    for (const auto& u : unexpected) UNSCOPED_INFO(u);
    CHECK(unexpected.empty());
}
