#include <sympp/core/function.hpp>

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>

namespace sympp {

namespace {

constexpr std::size_t kFunctionHashSeed = 0xF1F2'F3F4'5566'7788ULL;

[[nodiscard]] std::size_t hash_combine(std::size_t a, std::size_t b) noexcept {
    return a ^ (b + 0x9E37'79B9'7F4A'7C15ULL + (a << 6) + (a >> 2));
}

}  // namespace

Function::Function(std::vector<Expr> args) : args_(std::move(args)), hash_(0) {
    // FunctionId-aware hashing happens in compute_hash, called by subclasses
    // once their static identity is established. We seed with a placeholder
    // so a subclass forgetting to call compute_hash gets a stable (if useless)
    // hash rather than UB.
    hash_ = kFunctionHashSeed;
}

void Function::compute_hash(FunctionId fid) noexcept {
    std::size_t h = kFunctionHashSeed;
    h = hash_combine(h, static_cast<std::size_t>(fid));
    for (const auto& a : args_) {
        h = hash_combine(h, a->hash());
    }
    hash_ = h;
}

bool Function::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Function) return false;
    const auto& o = static_cast<const Function&>(other);
    if (function_id() != o.function_id()) return false;
    if (args_.size() != o.args_.size()) return false;
    for (std::size_t i = 0; i < args_.size(); ++i) {
        if (!(args_[i] == o.args_[i])) return false;
    }
    return true;
}

std::string Function::str() const {
    std::string out{name()};
    out.push_back('(');
    for (std::size_t i = 0; i < args_.size(); ++i) {
        if (i > 0) out.append(", ");
        out.append(args_[i]->str());
    }
    out.push_back(')');
    return out;
}

}  // namespace sympp
