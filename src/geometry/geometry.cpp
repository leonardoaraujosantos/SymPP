#include <sympp/geometry/geometry.hpp>

#include <stdexcept>
#include <vector>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

namespace {
[[nodiscard]] Expr sq(const Expr& e) { return pow(e, integer(2)); }
[[nodiscard]] bool is_zero_expr(const Expr& e) { return simplify(e) == S::Zero(); }
// (q - p) × (r - p) — the 2D cross product, zero iff p, q, r are collinear.
[[nodiscard]] Expr cross(const Point2D& p, const Point2D& q, const Point2D& r) {
    return (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
}
}  // namespace

Expr distance(const Point2D& p, const Point2D& q) {
    return simplify(sqrt(sq(q.x - p.x) + sq(q.y - p.y)));
}

Point2D midpoint(const Point2D& p, const Point2D& q) {
    return {simplify((p.x + q.x) * rational(1, 2)),
            simplify((p.y + q.y) * rational(1, 2))};
}

bool are_collinear(const Point2D& p, const Point2D& q, const Point2D& r) {
    return is_zero_expr(cross(p, q, r));
}

Expr triangle_area(const Point2D& p, const Point2D& q, const Point2D& r) {
    return simplify(abs(cross(p, q, r)) * rational(1, 2));
}

Line2D line_through(const Point2D& p, const Point2D& q) {
    Expr a = q.y - p.y;
    Expr b = p.x - q.x;
    Expr c = mul(integer(-1), a * p.x + b * p.y);
    return {simplify(a), simplify(b), simplify(c)};
}

Expr line_slope(const Line2D& l) {
    if (is_zero_expr(l.b)) throw std::runtime_error("line_slope: vertical line");
    return simplify(mul(integer(-1), l.a) * pow(l.b, integer(-1)));
}

bool lines_parallel(const Line2D& l1, const Line2D& l2) {
    return is_zero_expr(l1.a * l2.b - l2.a * l1.b);
}

bool lines_perpendicular(const Line2D& l1, const Line2D& l2) {
    return is_zero_expr(l1.a * l2.a + l1.b * l2.b);
}

std::optional<Point2D> line_intersection(const Line2D& l1, const Line2D& l2) {
    Expr det = l1.a * l2.b - l2.a * l1.b;
    if (is_zero_expr(det)) return std::nullopt;
    Expr inv = pow(det, integer(-1));
    Expr x = (l1.b * l2.c - l2.b * l1.c) * inv;
    Expr y = (l2.a * l1.c - l1.a * l2.c) * inv;
    return Point2D{simplify(x), simplify(y)};
}

Expr point_line_distance(const Point2D& p, const Line2D& l) {
    return simplify(abs(l.a * p.x + l.b * p.y + l.c) *
                    pow(sqrt(sq(l.a) + sq(l.b)), integer(-1)));
}

Expr polygon_area(const std::vector<Point2D>& v) {
    if (v.size() < 3) throw std::invalid_argument("polygon_area: need ≥ 3 vertices");
    Expr s = S::Zero();
    for (std::size_t i = 0; i < v.size(); ++i) {
        const Point2D& a = v[i];
        const Point2D& b = v[(i + 1) % v.size()];
        s = s + (a.x * b.y - b.x * a.y);
    }
    return simplify(s * rational(1, 2));
}

Expr polygon_perimeter(const std::vector<Point2D>& v) {
    if (v.size() < 2) throw std::invalid_argument("polygon_perimeter: need ≥ 2 vertices");
    Expr s = S::Zero();
    for (std::size_t i = 0; i < v.size(); ++i) {
        s = s + distance(v[i], v[(i + 1) % v.size()]);
    }
    return simplify(s);
}

}  // namespace sympp
