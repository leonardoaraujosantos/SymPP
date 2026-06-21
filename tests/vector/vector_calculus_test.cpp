// Vector calculus & differential geometry — cross-checked against SymPy.

#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/vector/vector_calculus.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
namespace vc = sympp::vector;

namespace {
auto& oracle() { return sympp::testing::Oracle::instance(); }

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
bool field_matches(const vc::VectorField& got, const std::string& joined) {
    auto want = split(joined);
    if (got.size() != want.size()) return false;
    for (std::size_t i = 0; i < got.size(); ++i) {
        if (!oracle().equivalent(got[i]->str(), want[i])) return false;
    }
    return true;
}
std::string vcalc(const std::string& fn, nlohmann::json req) {
    req["op"] = "vectorcalc";
    req["fn"] = fn;
    auto r = oracle().send(req);
    REQUIRE(r.ok);
    return r.result_str();
}
}  // namespace

TEST_CASE("gradient / divergence / curl / laplacian vs SymPy", "[vector][oracle]") {
    auto x = symbol("x"), y = symbol("y"), z = symbol("z");
    std::vector<Expr> xyz = {x, y, z};

    Expr f = x * y * z + pow(x, integer(2));
    REQUIRE(field_matches(vc::gradient(f, xyz),
                          vcalc("gradient", {{"f", f->str()}, {"vars", {"x", "y", "z"}}})));

    vc::VectorField F = {x * y, y * z, z * x};
    std::vector<std::string> Fs = {(x * y)->str(), (y * z)->str(), (z * x)->str()};
    REQUIRE(oracle().equivalent(
        vc::divergence(F, xyz)->str(),
        vcalc("divergence", {{"field", Fs}, {"vars", {"x", "y", "z"}}})));
    REQUIRE(field_matches(vc::curl(F, xyz),
                          vcalc("curl", {{"field", Fs}, {"vars", {"x", "y", "z"}}})));
    REQUIRE(oracle().equivalent(
        vc::laplacian(f, xyz)->str(),
        vcalc("laplacian", {{"f", f->str()}, {"vars", {"x", "y", "z"}}})));
}

TEST_CASE("vector calculus: known closed forms", "[vector]") {
    auto x = symbol("x"), y = symbol("y"), z = symbol("z");
    std::vector<Expr> xyz = {x, y, z};
    REQUIRE(vc::divergence({x, y, z}, xyz) == integer(3));          // div(r) = 3
    REQUIRE(vc::curl({x, y, z}, xyz)[0] == integer(0));            // curl(grad) = 0
    auto rot = vc::curl({mul(integer(-1), y), x, integer(0)}, xyz);
    REQUIRE(rot[2] == integer(2));                                 // curl(-y, x, 0) = 2 ẑ
    REQUIRE(vc::laplacian(pow(x, integer(2)) + pow(y, integer(2)) + pow(z, integer(2)), xyz)
            == integer(6));
}

TEST_CASE("differential geometry: Ricci scalar", "[vector][diffgeom][oracle]") {
    auto th = symbol("theta"), ph = symbol("phi");

    // Unit 2-sphere: ds² = dθ² + sin²θ dφ²  →  Ricci scalar = 2.
    Matrix sphere{{integer(1), integer(0)},
                  {integer(0), pow(sin(th), integer(2))}};
    REQUIRE(vc::ricci_scalar(sphere, {th, ph}) == integer(2));
    nlohmann::json req;
    req["op"] = "diffgeom";
    req["fn"] = "ricci_scalar";
    req["coords"] = nlohmann::json::array({"theta", "phi"});
    req["metric"] = nlohmann::json::array({nlohmann::json::array({"1", "0"}),
                                           nlohmann::json::array({"0", "sin(theta)**2"})});
    auto resp = oracle().send(req);
    REQUIRE(resp.ok);
    REQUIRE(oracle().equivalent(vc::ricci_scalar(sphere, {th, ph})->str(),
                                resp.result_str()));

    // Flat plane in polar coordinates: ds² = dr² + r² dφ²  →  Ricci scalar = 0.
    auto r = symbol("r");
    Matrix flat{{integer(1), integer(0)}, {integer(0), pow(r, integer(2))}};
    REQUIRE(vc::ricci_scalar(flat, {r, ph}) == integer(0));
}
