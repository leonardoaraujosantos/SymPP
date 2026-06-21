// Discrete transforms — cross-checked against sympy.discrete.

#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/discrete/discrete.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
namespace dsc = sympp::discrete;

namespace {
auto& oracle() { return sympp::testing::Oracle::instance(); }

std::vector<std::string> to_strs(const std::vector<Expr>& v) {
    std::vector<std::string> s;
    for (const auto& e : v) s.push_back(e->str());
    return s;
}
std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ';') { out.push_back(cur); cur.clear(); }
        else cur += c;
    }
    if (!cur.empty() || !s.empty()) out.push_back(cur);
    return out;
}
std::vector<std::string> ints(std::initializer_list<long> xs) {
    std::vector<std::string> s;
    for (long x : xs) s.push_back(std::to_string(x));
    return s;
}
// SymPP result vs SymPy result, element-wise SymPy-equivalence.
bool matches(const std::vector<Expr>& got, nlohmann::json req) {
    req["op"] = "discrete";
    auto r = oracle().send(req);
    REQUIRE(r.ok);
    auto want = split(r.result_str());
    if (got.size() != want.size()) return false;
    for (std::size_t i = 0; i < got.size(); ++i) {
        if (!oracle().equivalent(got[i]->str(), want[i])) return false;
    }
    return true;
}
std::vector<Expr> ivec(std::initializer_list<long> xs) {
    std::vector<Expr> v;
    for (long x : xs) v.push_back(integer(x));
    return v;
}
}  // namespace

TEST_CASE("FFT / IFFT vs SymPy", "[discrete][oracle]") {
    auto seq = ivec({1, 2, 3, 4});
    REQUIRE(matches(dsc::fft(seq), {{"fn", "fft"}, {"seq", ints({1, 2, 3, 4})}}));
    // round-trip
    REQUIRE(to_strs(dsc::ifft(dsc::fft(seq))) == to_strs(seq));
    // non-power-of-2 is zero-padded (matching SymPy).
    REQUIRE(matches(dsc::fft(ivec({1, 2, 3})), {{"fn", "fft"}, {"seq", ints({1, 2, 3})}}));
}

TEST_CASE("NTT / INTT vs SymPy", "[discrete][oracle]") {
    auto seq = ivec({1, 2, 3, 4});
    REQUIRE(to_strs(dsc::ntt(seq, integer(337))) == ints({10, 39, 335, 294}));
    REQUIRE(matches(dsc::ntt(seq, integer(337)),
                    {{"fn", "ntt"}, {"seq", ints({1, 2, 3, 4})}, {"p", "337"}}));
    REQUIRE(to_strs(dsc::intt(dsc::ntt(seq, integer(337)), integer(337))) == to_strs(seq));
}

TEST_CASE("convolution vs SymPy", "[discrete][oracle]") {
    auto c = dsc::convolution(ivec({1, 2, 3}), ivec({4, 5, 6}));
    REQUIRE(to_strs(c) == ints({4, 13, 28, 27, 18}));
    REQUIRE(matches(c, {{"fn", "convolution"},
                        {"seq", ints({1, 2, 3})}, {"seq2", ints({4, 5, 6})}}));
}

TEST_CASE("Mobius (subset) transform vs SymPy", "[discrete][oracle]") {
    auto seq = ivec({1, 2, 3, 4});
    REQUIRE(to_strs(dsc::mobius_transform(seq)) == ints({1, 3, 4, 10}));
    REQUIRE(matches(dsc::mobius_transform(seq), {{"fn", "mobius"}, {"seq", ints({1, 2, 3, 4})}}));
    REQUIRE(to_strs(dsc::inverse_mobius_transform(dsc::mobius_transform(seq))) == to_strs(seq));
}
