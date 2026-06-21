// Plane geometry — cross-checked against SymPy's geometry module.

#include <optional>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/geometry/geometry.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
Point2D P(long x, long y) { return {integer(x), integer(y)}; }

std::string geo(const std::string& fn, const std::vector<std::string>& coords) {
    nlohmann::json req = {{"op", "geometry"}, {"fn", fn}, {"coords", coords}};
    auto r = sympp::testing::Oracle::instance().send(req);
    REQUIRE(r.ok);
    return r.result_str();
}
bool equiv(const Expr& a, const std::string& sympy_str) {
    return sympp::testing::Oracle::instance().equivalent(a->str(), sympy_str);
}
}  // namespace

TEST_CASE("distance / midpoint", "[geometry][oracle]") {
    REQUIRE(distance(P(1, 2), P(4, 6)) == integer(5));
    REQUIRE(equiv(distance(P(0, 0), P(2, 3)), geo("distance", {"0", "0", "2", "3"})));
    auto m = midpoint(P(1, 2), P(3, 8));
    REQUIRE(m.x == integer(2));
    REQUIRE(m.y == integer(5));
}

TEST_CASE("collinearity and triangle area", "[geometry][oracle]") {
    REQUIRE(are_collinear(P(0, 0), P(1, 1), P(2, 2)));
    REQUIRE_FALSE(are_collinear(P(0, 0), P(1, 0), P(0, 1)));
    REQUIRE(geo("collinear", {"0", "0", "1", "1", "2", "2"}) == "true");

    // Right triangle with legs 3 and 4 → area 6.
    REQUIRE(triangle_area(P(0, 0), P(4, 0), P(0, 3)) == integer(6));
}

TEST_CASE("lines: slope, parallel, perpendicular, intersection", "[geometry]") {
    auto l1 = line_through(P(0, 0), P(2, 4));   // y = 2x  → slope 2
    REQUIRE(line_slope(l1) == integer(2));

    auto l2 = line_through(P(0, 1), P(1, 3));   // slope 2 → parallel to l1
    REQUIRE(lines_parallel(l1, l2));
    REQUIRE_FALSE(line_intersection(l1, l2).has_value());

    // x-axis ⟂ y-axis.
    REQUIRE(lines_perpendicular(line_through(P(0, 0), P(1, 0)),
                                line_through(P(0, 0), P(0, 1))));

    // The diagonals of the unit square meet at (1, 1).
    auto hit = line_intersection(line_through(P(0, 0), P(2, 2)),
                                 line_through(P(0, 2), P(2, 0)));
    REQUIRE(hit.has_value());
    REQUIRE(hit->x == integer(1));
    REQUIRE(hit->y == integer(1));
}

TEST_CASE("point-to-line distance", "[geometry]") {
    auto l = line_through(P(0, 0), P(1, 0));   // the x-axis
    REQUIRE(point_line_distance(P(3, 4), l) == integer(4));
}

TEST_CASE("polygon area and perimeter", "[geometry][oracle]") {
    std::vector<Point2D> unit_square = {P(0, 0), P(1, 0), P(1, 1), P(0, 1)};
    REQUIRE(polygon_area(unit_square) == integer(1));
    REQUIRE(polygon_perimeter(unit_square) == integer(4));
    REQUIRE(equiv(polygon_area(unit_square),
                  geo("polygon_area", {"0", "0", "1", "0", "1", "1", "0", "1"})));

    std::vector<Point2D> tri = {P(0, 0), P(4, 0), P(0, 3)};
    REQUIRE(polygon_area(tri) == integer(6));
}
