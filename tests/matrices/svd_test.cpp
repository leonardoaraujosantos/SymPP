// Symbolic singular values — cross-checked against SymPy's
// Matrix.singular_values().

#include <sstream>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/matrices/matrix.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, sep)) out.push_back(item);
    return out;
}
std::vector<std::string> sympy_singular_values(const std::vector<std::vector<std::string>>& rows) {
    nlohmann::json req = {{"op", "matrix"}, {"fn", "singular_values"}, {"rows", rows}};
    auto r = sympp::testing::Oracle::instance().send(req);
    REQUIRE(r.ok);
    return split(r.result_str(), ';');
}
bool equiv(const Expr& a, const std::string& b) {
    return sympp::testing::Oracle::instance().equivalent(a->str(), b);
}
}  // namespace

TEST_CASE("singular values match SymPy (descending)", "[matrix][svd][oracle]") {
    auto check = [](const Matrix& m, const std::vector<std::vector<std::string>>& rows) {
        auto got = m.singular_values();
        auto want = sympy_singular_values(rows);
        REQUIRE(got.size() == want.size());
        for (std::size_t i = 0; i < got.size(); ++i) REQUIRE(equiv(got[i], want[i]));
    };

    check(Matrix{{integer(3), integer(0)}, {integer(0), integer(4)}},
          {{"3", "0"}, {"0", "4"}});
    check(Matrix{{integer(2), integer(0)}, {integer(0), integer(0)}},
          {{"2", "0"}, {"0", "0"}});
    check(Matrix{{integer(1), integer(1)}, {integer(0), integer(1)}},
          {{"1", "1"}, {"0", "1"}});
    // Non-square 3×2.
    check(Matrix{{integer(1), integer(2)}, {integer(3), integer(4)}, {integer(5), integer(6)}},
          {{"1", "2"}, {"3", "4"}, {"5", "6"}});
}

TEST_CASE("singular values: diagonal and known cases", "[matrix][svd]") {
    Matrix d{{integer(5), integer(0)}, {integer(0), integer(2)}};
    auto sv = d.singular_values();
    REQUIRE(sv.size() == 2);
    REQUIRE(sv[0] == integer(5));  // descending
    REQUIRE(sv[1] == integer(2));
}

TEST_CASE("full SVD reconstructs A = U S Vᵀ", "[matrix][svd][oracle]") {
    auto reconstructs = [](const Matrix& A) {
        auto [U, S, V] = A.svd();
        Matrix R = U * S * V.transpose();
        REQUIRE(R.rows() == A.rows());
        REQUIRE(R.cols() == A.cols());
        for (std::size_t i = 0; i < A.rows(); ++i) {
            for (std::size_t j = 0; j < A.cols(); ++j) {
                if (!sympp::testing::Oracle::instance().equivalent(R.at(i, j)->str(),
                                                                   A.at(i, j)->str())) {
                    return false;
                }
            }
        }
        // S is diagonal with descending positive entries.
        for (std::size_t k = 0; k + 1 < S.rows(); ++k) {
            REQUIRE(S.at(k, k + 1) == integer(0));
        }
        return true;
    };

    REQUIRE(reconstructs(Matrix{{integer(3), integer(0)}, {integer(0), integer(4)}}));
    REQUIRE(reconstructs(Matrix{{integer(2), integer(1)}, {integer(1), integer(2)}}));
    REQUIRE(reconstructs(Matrix{{integer(1), integer(1)}, {integer(0), integer(1)}}));
    REQUIRE(reconstructs(Matrix{{integer(1), integer(2)}, {integer(2), integer(4)}}));  // rank 1
}
