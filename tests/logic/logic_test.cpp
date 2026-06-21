// Propositional logic — satisfiable / to_cnf / to_dnf / simplify_logic,
// cross-checked against SymPy via the `logic` oracle op.

#include <map>
#include <optional>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/core/boolean.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/logic/logic.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
bool sympy_satisfiable(const Expr& e) {
    nlohmann::json req = {{"op", "logic"}, {"fn", "satisfiable"}, {"expr", e->str()}};
    auto r = sympp::testing::Oracle::instance().send(req);
    REQUIRE(r.ok);
    return r.result_bool();
}
bool sympy_equivalent(const Expr& a, const Expr& b) {
    nlohmann::json req = {{"op", "logic"}, {"fn", "equivalent"},
                          {"a", a->str()}, {"b", b->str()}};
    auto r = sympp::testing::Oracle::instance().send(req);
    REQUIRE(r.ok);
    return r.result_bool();
}
}  // namespace

TEST_CASE("connectives auto-evaluate", "[logic]") {
    auto a = symbol("a"), b = symbol("b");
    REQUIRE(bool_and(a, S::True()) == a);
    REQUIRE(bool_and(a, S::False()) == S::False());
    REQUIRE(bool_or(a, S::False()) == a);
    REQUIRE(bool_or(a, S::True()) == S::True());
    REQUIRE(bool_not(bool_not(a)) == a);
    REQUIRE(bool_and(a, bool_not(a)) == S::False());
    REQUIRE(bool_or(a, bool_not(a)) == S::True());
    REQUIRE(bool_and(a, a) == a);               // idempotent
    REQUIRE(bool_and(a, b) == bool_and(b, a));  // commutative (canonical order)
}

TEST_CASE("satisfiable agrees with SymPy (SAT/UNSAT)", "[logic][oracle]") {
    auto a = symbol("a"), b = symbol("b"), c = symbol("c");
    std::vector<Expr> formulas = {
        bool_and(a, bool_not(a)),                                  // UNSAT
        bool_or(a, bool_not(a)),                                   // valid
        implies(a, b),                                             // SAT
        bool_and(implies(a, b), bool_and(a, bool_not(b))),         // UNSAT
        bool_and(bool_xor(a, b), bool_xor(b, c)),                  // SAT
        equivalent(a, bool_not(a)),                                // UNSAT
    };
    for (const auto& f : formulas) {
        auto model = satisfiable(f);
        REQUIRE(model.has_value() == sympy_satisfiable(f));
    }
}

TEST_CASE("satisfiable returns a genuine model", "[logic]") {
    auto a = symbol("a"), b = symbol("b");
    auto f = bool_and(implies(a, b), a);  // forces a=T, b=T
    auto model = satisfiable(f);
    REQUIRE(model.has_value());
    REQUIRE((*model)["a"] == true);
    REQUIRE((*model)["b"] == true);
    REQUIRE_FALSE(satisfiable(bool_and(a, bool_not(a))).has_value());
}

TEST_CASE("to_cnf / to_dnf / simplify_logic preserve meaning (vs SymPy)", "[logic][oracle]") {
    auto a = symbol("a"), b = symbol("b"), c = symbol("c");
    std::vector<Expr> formulas = {
        bool_or(bool_and(a, b), bool_and(a, bool_not(b))),  // = a
        implies(a, bool_and(b, c)),
        bool_xor(a, bool_xor(b, c)),
        equivalent(a, b),
    };
    for (const auto& f : formulas) {
        REQUIRE(sympy_equivalent(to_cnf(f), f));
        REQUIRE(sympy_equivalent(to_dnf(f), f));
        REQUIRE(sympy_equivalent(simplify_logic(f), f));
    }
}

TEST_CASE("simplify_logic minimizes", "[logic][oracle]") {
    auto a = symbol("a"), b = symbol("b");
    // a&b | a&~b  ≡  a
    auto f = bool_or(bool_and(a, b), bool_and(a, bool_not(b)));
    REQUIRE(simplify_logic(f) == a);
    // (a|b) & (a|~b)  ≡  a
    auto g = bool_and(bool_or(a, b), bool_or(a, bool_not(b)));
    REQUIRE(simplify_logic(g) == a);
    REQUIRE(sympy_equivalent(simplify_logic(g), g));
}
