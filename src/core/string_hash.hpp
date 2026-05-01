#pragma once

// Internal helper — deterministic string hash.
//
// std::hash<std::string> is implementation-defined; libc++ (macOS Apple
// clang) and libstdc++ (Linux gcc / clang) produce different values for
// identical inputs. SymPP uses string-derived hashes to seed the Add /
// Mul canonical sort key, so an unstable hash makes simplified-form
// `str()` output platform-dependent — tests that pin a specific
// argument order break across stdlib implementations.
//
// FNV-1a 64-bit is small, fast, and bit-exact across platforms. Used
// inside Symbol::hash_ and UndefinedFunction::hash_.
//
// Reference: http://www.isthe.com/chongo/tech/comp/fnv/index.html

#include <cstddef>
#include <cstdint>
#include <string_view>

namespace sympp::detail {

[[nodiscard]] constexpr std::size_t fnv1a_64(std::string_view s) noexcept {
    std::uint64_t h = 0xcbf2'9ce4'8422'2325ULL;
    for (char c : s) {
        h ^= static_cast<std::uint64_t>(static_cast<unsigned char>(c));
        h *= 0x0000'0100'0000'01B3ULL;
    }
    return static_cast<std::size_t>(h);
}

}  // namespace sympp::detail
