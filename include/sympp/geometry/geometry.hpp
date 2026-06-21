#pragma once

// Analytic plane geometry over exact symbolic coordinates — points, segments,
// lines and polygons, with the usual measurements and predicates.
//
//   Point2D a{integer(1), integer(2)}, b{integer(4), integer(6)};
//   distance(a, b);            // → 5
//   Line2D l = line_through(a, b);
//   line_slope(l);             // → 4/3
//
// Reference: sympy/geometry/{point,line,polygon}.py.

#include <optional>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

struct Point2D {
    Expr x, y;
};

// A line as the implicit form a·x + b·y + c = 0.
struct Line2D {
    Expr a, b, c;
};

// ----- Points & segments -----------------------------------------------------
[[nodiscard]] SYMPP_EXPORT Expr distance(const Point2D& p, const Point2D& q);
[[nodiscard]] SYMPP_EXPORT Point2D midpoint(const Point2D& p, const Point2D& q);
[[nodiscard]] SYMPP_EXPORT bool are_collinear(const Point2D& p, const Point2D& q,
                                              const Point2D& r);
[[nodiscard]] SYMPP_EXPORT Expr triangle_area(const Point2D& p, const Point2D& q,
                                              const Point2D& r);

// ----- Lines -----------------------------------------------------------------
[[nodiscard]] SYMPP_EXPORT Line2D line_through(const Point2D& p, const Point2D& q);
// Slope −a/b. Throws for a vertical line (b = 0).
[[nodiscard]] SYMPP_EXPORT Expr line_slope(const Line2D& l);
[[nodiscard]] SYMPP_EXPORT bool lines_parallel(const Line2D& l1, const Line2D& l2);
[[nodiscard]] SYMPP_EXPORT bool lines_perpendicular(const Line2D& l1, const Line2D& l2);
// Intersection point, or std::nullopt when the lines are parallel.
[[nodiscard]] SYMPP_EXPORT std::optional<Point2D> line_intersection(const Line2D& l1,
                                                                    const Line2D& l2);
[[nodiscard]] SYMPP_EXPORT Expr point_line_distance(const Point2D& p, const Line2D& l);

// ----- Polygons --------------------------------------------------------------
// Signed area (shoelace; positive for counter-clockwise vertex order, matching
// SymPy's Polygon.area) and perimeter.
[[nodiscard]] SYMPP_EXPORT Expr polygon_area(const std::vector<Point2D>& vertices);
[[nodiscard]] SYMPP_EXPORT Expr polygon_perimeter(const std::vector<Point2D>& vertices);

}  // namespace sympp
