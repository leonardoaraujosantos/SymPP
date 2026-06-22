// SymPP micro-benchmarks — wall-clock timing of representative symbolic
// operations. Dependency-free (just <chrono>); run as `sympp_bench`.
//
// Output is a fixed-width table of median nanoseconds-per-op over a small
// repeat count. It is a smoke-grade performance harness for spotting
// regressions, not a statistically rigorous benchmark framework.

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include <sympp/sympp.hpp>

using namespace sympp;
using Clock = std::chrono::steady_clock;

namespace {

// Median wall-clock per call over `reps` repetitions (each repetition times a
// fresh evaluation so interning/caches do not hide real cost).
double median_ns(const std::function<void()>& fn, int reps) {
    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(reps));
    for (int i = 0; i < reps; ++i) {
        auto t0 = Clock::now();
        fn();
        auto t1 = Clock::now();
        samples.push_back(
            std::chrono::duration<double, std::nano>(t1 - t0).count());
    }
    std::sort(samples.begin(), samples.end());
    return samples[samples.size() / 2];
}

struct Case {
    std::string name;
    std::function<void()> fn;
    int reps;
};

}  // namespace

int main() {
    auto x = symbol("x");
    auto y = symbol("y");
    auto oo = S::Infinity();

    std::vector<Case> cases;
    auto add = [&](std::string n, std::function<void()> f, int reps = 200) {
        cases.push_back({std::move(n), std::move(f), reps});
    };

    add("build x^3 + 2xy", [&] {
        volatile auto e = pow(x, integer(3)) + integer(2) * x * y;
        (void)e;
    });
    add("diff x^10", [&] {
        volatile auto e = diff(pow(x, integer(10)), x);
        (void)e;
    });
    add("expand (x+1)^8", [&] {
        volatile auto e = expand(pow(x + integer(1), integer(8)));
        (void)e;
    });
    add("factor x^4 - 1", [&] {
        volatile auto e = factor(pow(x, integer(4)) - integer(1), x);
        (void)e;
    });
    add("simplify sin^2+cos^2", [&] {
        volatile auto e =
            simplify(pow(sin(x), integer(2)) + pow(cos(x), integer(2)));
        (void)e;
    });
    add("integrate log(x)", [&] {
        volatile auto e = integrate(log(x), x);
        (void)e;
    });
    add("limit sin(x)/x @0", [&] {
        volatile auto e = limit(sin(x) * pow(x, integer(-1)), x, S::Zero());
        (void)e;
    });
    add("limit (1+1/x)^x @oo", [&] {
        volatile auto e =
            limit(pow(integer(1) + pow(x, integer(-1)), x), x, oo);
        (void)e;
    });
    add("solve x^2-5x+6", [&] {
        auto r = solve(pow(x, integer(2)) - integer(5) * x + integer(6), x);
        volatile auto n = r.size();
        (void)n;
    });
    add("det 3x3", [&] {
        Matrix m = {{integer(2), integer(1), integer(0)},
                    {integer(1), integer(3), integer(1)},
                    {integer(0), integer(1), integer(4)}};
        volatile auto d = m.det();
        (void)d;
    });

    std::printf("%-26s %14s\n", "operation", "median ns/op");
    std::printf("%-26s %14s\n", "--------------------------", "--------------");
    for (const auto& c : cases) {
        double ns = median_ns(c.fn, c.reps);
        std::printf("%-26s %14.0f\n", c.name.c_str(), ns);
    }
    return 0;
}
