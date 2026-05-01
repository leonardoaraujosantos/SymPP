#include <sympp/core/expr_cache.hpp>

#include <cstddef>
#include <memory>
#include <mutex>
#include <utility>

#include <sympp/core/basic.hpp>

namespace sympp {

ExprCache& ExprCache::instance() {
    static ExprCache cache;
    return cache;
}

Expr ExprCache::intern(Expr candidate) {
    if (!candidate) return candidate;
    if (!enabled_) return candidate;

    auto h = candidate->hash();
    std::lock_guard lock(mtx_);

    auto range = buckets_.equal_range(h);
    for (auto it = range.first; it != range.second;) {
        auto live = it->second.lock();
        if (!live) {
            // Lazy GC: drop expired weak_ptrs as we walk the bucket.
            it = buckets_.erase(it);
            continue;
        }
        if (live->equals(*candidate)) {
            return live;  // cache hit — candidate will be deallocated
        }
        ++it;
    }

    // Cache miss — record candidate and return it.
    buckets_.emplace(h, std::weak_ptr<const Basic>{candidate});
    return candidate;
}

std::size_t ExprCache::cleanup() {
    std::lock_guard lock(mtx_);
    std::size_t removed = 0;
    for (auto it = buckets_.begin(); it != buckets_.end();) {
        if (it->second.expired()) {
            it = buckets_.erase(it);
            ++removed;
        } else {
            ++it;
        }
    }
    return removed;
}

std::size_t ExprCache::size() const noexcept {
    std::lock_guard lock(mtx_);
    return buckets_.size();
}

void ExprCache::set_enabled(bool enabled) noexcept {
    std::lock_guard lock(mtx_);
    enabled_ = enabled;
}

bool ExprCache::enabled() const noexcept {
    std::lock_guard lock(mtx_);
    return enabled_;
}

}  // namespace sympp
