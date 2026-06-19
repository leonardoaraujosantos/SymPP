#include <sympp/core/assumption_context.hpp>

#include <utility>
#include <vector>

#include <sympp/core/basic.hpp>

namespace sympp {

namespace {

// Stack of active `assuming` scopes, innermost last. Each frame binds one key
// expression to its closed assumption mask. Thread-local so concurrent solver
// threads do not see each other's scopes.
thread_local std::vector<std::pair<Expr, AssumptionMask>> g_frames;

}  // namespace

bool assumption_context_active() noexcept { return !g_frames.empty(); }

std::optional<bool> assumption_context_get(const Expr& e,
                                           AssumptionKey k) noexcept {
    if (!e || g_frames.empty()) return std::nullopt;
    // Innermost scope wins: walk from the back.
    for (auto it = g_frames.rbegin(); it != g_frames.rend(); ++it) {
        if (it->first == e) {
            if (auto v = it->second.get(k); v.has_value()) return v;
        }
    }
    return std::nullopt;
}

assuming::assuming(const Expr& target, AssumptionMask facts) {
    g_frames.emplace_back(target, close_assumptions(std::move(facts)));
}

assuming::~assuming() {
    if (!g_frames.empty()) g_frames.pop_back();
}

}  // namespace sympp
