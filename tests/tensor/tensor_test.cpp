// Dense tensor algebra — contraction / tensor product cross-checked against
// sympy.tensor.array, plus metric raise/lower round-trips.

#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/tensor/tensor.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
namespace tn = sympp::tensor;

namespace {
auto& oracle() { return sympp::testing::Oracle::instance(); }

std::vector<std::string> data_strs(const tn::Tensor& t) {
    std::vector<std::string> s;
    for (const auto& e : t.data()) s.push_back(e->str());
    return s;
}
std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ';') { out.push_back(cur); cur.clear(); }
        else cur += c;
    }
    out.push_back(cur);
    return out;
}
// Element-wise SymPy-equivalence of a tensor's flat data against a ';'-joined ref.
bool flat_equiv(const tn::Tensor& t, const std::string& ref) {
    auto got = data_strs(t);
    auto want = split(ref);
    if (got.size() != want.size()) return false;
    for (std::size_t i = 0; i < got.size(); ++i) {
        if (!oracle().equivalent(got[i], want[i])) return false;
    }
    return true;
}
}  // namespace

TEST_CASE("tensor contraction matches SymPy (trace)", "[tensor][oracle]") {
    tn::Tensor m({2, 2}, {integer(1), integer(2), integer(3), integer(4)});
    auto tr = tn::contract(m, 0, 1);
    REQUIRE(tr.rank() == 0);
    REQUIRE(tr.data()[0] == integer(5));

    nlohmann::json req = {{"op", "tensor"}, {"fn", "contract"}};
    req["a"]["shape"] = nlohmann::json::array({2, 2});
    req["a"]["data"] = nlohmann::json::array({"1", "2", "3", "4"});
    req["axes"] = nlohmann::json::array({0, 1});
    auto r = oracle().send(req);
    REQUIRE(r.ok);
    REQUIRE(oracle().equivalent(tr.data()[0]->str(), r.result_str()));

    // Identity contracts to its dimension.
    tn::Tensor id({2, 2}, {integer(1), integer(0), integer(0), integer(1)});
    REQUIRE(tn::contract(id, 0, 1).data()[0] == integer(2));
}

TEST_CASE("tensor product matches SymPy", "[tensor][oracle]") {
    auto a = symbol("a"), b = symbol("b"), c = symbol("c"), d = symbol("d");
    tn::Tensor u({2}, {a, b}), v({2}, {c, d});
    auto p = tn::tensor_product(u, v);
    REQUIRE(p.shape() == std::vector<std::size_t>{2, 2});

    nlohmann::json req = {{"op", "tensor"}, {"fn", "tensorproduct"}};
    req["a"]["shape"] = nlohmann::json::array({2});
    req["a"]["data"] = nlohmann::json::array({"a", "b"});
    req["b"]["shape"] = nlohmann::json::array({2});
    req["b"]["data"] = nlohmann::json::array({"c", "d"});
    auto r = oracle().send(req);
    REQUIRE(r.ok);
    REQUIRE(flat_equiv(p, r.result_str()));
}

TEST_CASE("metric raising/lowering round-trips", "[tensor]") {
    // Vector v^a = (x, y); metric g = diag(2, 3); inverse = diag(1/2, 1/3).
    auto x = symbol("x"), y = symbol("y");
    tn::Tensor v({2}, {x, y});
    tn::Tensor g({2, 2}, {integer(2), integer(0), integer(0), integer(3)});
    tn::Tensor gi({2, 2}, {rational(1, 2), integer(0), integer(0), rational(1, 3)});

    auto lowered = tn::lower_index(v, 0, g);     // v_a = (2x, 3y)
    REQUIRE(lowered.data()[0] == mul(integer(2), x));
    auto restored = tn::raise_index(lowered, 0, gi);  // back to (x, y)
    REQUIRE(restored == v);
}
