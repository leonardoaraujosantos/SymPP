// Combinatorics & group theory — permutations, permutation groups and integer
// partitions, cross-checked against sympy.combinatorics.

#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/combinatorics/combinatorics.hpp>
#include <sympp/core/integer.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
namespace cb = sympp::combinatorics;

namespace {
auto& oracle() { return sympp::testing::Oracle::instance(); }

std::string comb_oracle(nlohmann::json req) {
    req["op"] = "combinatorics";
    auto r = oracle().send(req);
    REQUIRE(r.ok);
    return r.result_str();
}

nlohmann::json gens_json(const cb::PermutationGroup& g) {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& p : g.generators()) arr.push_back(p.array_form());
    return arr;
}

void check_group(const cb::PermutationGroup& g, const std::string& name) {
    INFO(name);
    REQUIRE(std::to_string(g.order()) ==
            comb_oracle({{"fn", "group_order"}, {"gens", gens_json(g)}}));
    std::string ab = g.is_abelian() ? "true" : "false";
    REQUIRE(ab == comb_oracle({{"fn", "is_abelian"}, {"gens", gens_json(g)}}));
}
}  // namespace

TEST_CASE("permutation: composition, inverse, sign, order", "[combinatorics]") {
    cb::Permutation p({1, 2, 0});  // 3-cycle (0 1 2)
    REQUIRE(p.order() == 3);
    REQUIRE(p.sign() == 1);
    // p∘p⁻¹ = identity.
    REQUIRE(p.compose(p.inverse()) == cb::Permutation::identity(3));
    // p² = (0 2 1), p³ = id.
    REQUIRE(p.compose(p) == cb::Permutation({2, 0, 1}));
    REQUIRE(p.compose(p).compose(p).is_identity());

    cb::Permutation q({1, 0, 3, 2});  // (0 1)(2 3): two transpositions
    REQUIRE(q.sign() == 1);
    REQUIRE(q.order() == 2);
    REQUIRE(q.cyclic_form().size() == 2);

    cb::Permutation t({1, 0, 2});  // single transposition (0 1)
    REQUIRE(t.sign() == -1);
}

TEST_CASE("permutation: sign and order match SymPy", "[combinatorics][oracle]") {
    auto check = [&](std::vector<int> img) {
        cb::Permutation p(img);
        nlohmann::json j = img;
        REQUIRE(std::to_string(p.sign()) ==
                comb_oracle({{"fn", "perm_sign"}, {"p", j}}));
        REQUIRE(std::to_string(p.order()) ==
                comb_oracle({{"fn", "perm_order"}, {"p", j}}));
    };
    check({1, 2, 0});
    check({1, 0, 3, 2});
    check({2, 3, 4, 0, 1});
    check({0, 1, 2, 3, 4});  // identity
    check({4, 3, 2, 1, 0});
}

TEST_CASE("permutation groups: order and abelian-ness match SymPy",
          "[combinatorics][oracle]") {
    check_group(cb::symmetric_group(3), "S3");
    check_group(cb::symmetric_group(4), "S4");
    check_group(cb::cyclic_group(5), "C5");
    check_group(cb::cyclic_group(6), "C6");
    check_group(cb::dihedral_group(4), "D4");
    check_group(cb::dihedral_group(5), "D5");
    check_group(cb::alternating_group(4), "A4");

    // Spot-check the canonical orders directly.
    REQUIRE(cb::symmetric_group(4).order() == 24);   // 4!
    REQUIRE(cb::dihedral_group(5).order() == 10);     // 2n
    REQUIRE(cb::alternating_group(4).order() == 12);  // n!/2
    REQUIRE(cb::cyclic_group(7).is_abelian());
    REQUIRE_FALSE(cb::symmetric_group(3).is_abelian());
}

TEST_CASE("permutation group: membership", "[combinatorics]") {
    auto a4 = cb::alternating_group(4);
    // A 3-cycle is even → in A4; a transposition is odd → not in A4.
    REQUIRE(a4.contains(cb::Permutation({1, 2, 0, 3})));        // (0 1 2)
    REQUIRE_FALSE(a4.contains(cb::Permutation({1, 0, 2, 3})));  // (0 1)
}

TEST_CASE("integer partitions: counts match SymPy", "[combinatorics][oracle]") {
    for (int n : {0, 1, 2, 5, 10, 20, 50, 100}) {
        REQUIRE(cb::partition_count(n)->str() ==
                comb_oracle({{"fn", "partition"}, {"n", n}}));
    }
    // The enumerated partitions number p(n) and each sums to n.
    for (int n : {1, 4, 6, 7}) {
        auto parts = cb::integer_partitions(n);
        REQUIRE(static_cast<int>(parts.size()) ==
                std::stoi(cb::partition_count(n)->str()));
        for (const auto& part : parts) {
            int s = 0;
            for (int v : part) s += v;
            REQUIRE(s == n);
        }
    }
}

// ----- Group actions: orbits & Pólya/Burnside enumeration --------------------

namespace {
nlohmann::json gens_of(const cb::PermutationGroup& g) {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& p : g.generators()) arr.push_back(p.array_form());
    return arr;
}
std::string orbits_str(const cb::PermutationGroup& g) {
    std::string s;
    bool first_o = true;
    for (const auto& o : cb::orbits(g)) {
        if (!first_o) s += ";";
        first_o = false;
        bool first = true;
        for (int v : o) {
            if (!first) s += ",";
            first = false;
            s += std::to_string(v);
        }
    }
    return s;
}
}  // namespace

TEST_CASE("group orbits match SymPy", "[combinatorics][oracle]") {
    auto check = [&](const cb::PermutationGroup& g, const std::string& name) {
        INFO(name);
        REQUIRE(orbits_str(g) == comb_oracle({{"fn", "orbits"}, {"gens", gens_of(g)}}));
    };
    check(cb::cyclic_group(4), "C4");
    check(cb::dihedral_group(5), "D5");
    check(cb::symmetric_group(4), "S4");
    // A group acting with two orbits: ⟨(0 1)(2 3 4)⟩.
    cb::PermutationGroup mixed({cb::Permutation({1, 0, 3, 4, 2})});
    check(mixed, "(01)(234)");
    REQUIRE(cb::orbit(cb::cyclic_group(4), 0) == std::vector<int>{0, 1, 2, 3});
}

TEST_CASE("Pólya/Burnside coloring counts match SymPy", "[combinatorics][oracle]") {
    auto check = [&](const cb::PermutationGroup& g, int k, const std::string& name) {
        INFO(name);
        REQUIRE(cb::colorings_count(g, k)->str() ==
                comb_oracle({{"fn", "burnside"}, {"gens", gens_of(g)}, {"k", std::to_string(k)}}));
    };
    // Necklaces (cyclic) and bracelets (dihedral) — classic Burnside counts.
    check(cb::cyclic_group(4), 2, "C4 k=2");   // 6 binary necklaces
    check(cb::cyclic_group(4), 3, "C4 k=3");   // 24
    check(cb::cyclic_group(6), 2, "C6 k=2");   // 14
    check(cb::dihedral_group(4), 2, "D4 k=2"); // 6 binary bracelets
    check(cb::symmetric_group(3), 2, "S3 k=2");
    check(cb::dihedral_group(5), 3, "D5 k=3");

    // Spot value: 6 binary necklaces of length 4.
    REQUIRE(cb::colorings_count(cb::cyclic_group(4), 2) == integer(6));
}

// BSGS-1: Schreier–Sims gives correct order/membership for LARGE groups that the
// BFS closure cannot enumerate (S₁₂ has 479 001 600 elements).
TEST_CASE("permutation groups: Schreier–Sims order & membership (large groups)",
          "[combinatorics][oracle][bsgs]") {
    auto fact = [](long n) { long f = 1; for (long i = 2; i <= n; ++i) f *= i; return f; };
    for (long n : {6, 8, 10, 12}) {
        REQUIRE(cb::symmetric_group(static_cast<int>(n)).order() == fact(n));        // n!
        REQUIRE(cb::alternating_group(static_cast<int>(n)).order() == fact(n) / 2);  // n!/2
        REQUIRE(cb::cyclic_group(static_cast<int>(n)).order() == n);
        REQUIRE(cb::dihedral_group(static_cast<int>(n)).order() == 2 * n);
    }
    // Cross-check a couple against SymPy's own group order.
    check_group(cb::symmetric_group(7), "S7");
    check_group(cb::alternating_group(6), "A6");

    // Membership by sifting in A₈ (8! / 2 = 20160 elements): even ✓, odd ✗.
    auto a8 = cb::alternating_group(8);
    REQUIRE(a8.contains(cb::Permutation({1, 2, 0, 3, 4, 5, 6, 7})));        // (0 1 2) even
    REQUIRE_FALSE(a8.contains(cb::Permutation({1, 0, 2, 3, 4, 5, 6, 7})));  // (0 1) odd
}
