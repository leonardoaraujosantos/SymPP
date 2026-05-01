#pragma once

// SymPy oracle — long-lived Python subprocess running oracle/oracle.py.
// See docs/05-validation-strategy.md for the protocol and rationale.
//
// Usage:
//     auto& oracle = sympp::testing::Oracle::instance();
//     auto resp = oracle.send({{"op", "ping"}});
//     REQUIRE(resp.ok);
//
// Or via convenience helpers:
//     REQUIRE(oracle.equivalent("sin(x)**2 + cos(x)**2", "1"));

#include <mutex>
#include <string>
#include <string_view>

#include <nlohmann/json.hpp>

namespace sympp::testing {

struct OracleResponse {
    bool ok = false;
    nlohmann::json raw;

    [[nodiscard]] std::string result_str() const;
    [[nodiscard]] bool result_bool() const;
    [[nodiscard]] std::string error() const;
    [[nodiscard]] std::string detail() const;
};

class Oracle {
public:
    Oracle(const Oracle&) = delete;
    Oracle& operator=(const Oracle&) = delete;
    Oracle(Oracle&&) = delete;
    Oracle& operator=(Oracle&&) = delete;

    static Oracle& instance();

    OracleResponse send(nlohmann::json request);

    // -- Convenience helpers (thin wrappers over send) --
    [[nodiscard]] std::string srepr(std::string_view expr);
    [[nodiscard]] std::string sympy_str(std::string_view expr);
    [[nodiscard]] bool equivalent(std::string_view a, std::string_view b);
    [[nodiscard]] std::string evalf(std::string_view expr, int prec = 15);
    [[nodiscard]] std::string diff(std::string_view expr,
                                   std::string_view var,
                                   int order = 1);
    [[nodiscard]] std::string integrate(std::string_view expr, std::string_view var);
    [[nodiscard]] std::string simplify(std::string_view expr);
    [[nodiscard]] std::string sympy_version();

private:
    Oracle();
    ~Oracle();

    void spawn();
    void shutdown_subprocess() noexcept;

    int pid_ = -1;
    int stdin_fd_ = -1;   // we write -> subprocess stdin
    int stdout_fd_ = -1;  // we read  <- subprocess stdout
    int next_id_ = 1;
    std::string read_buf_;  // accumulator for partial lines
    std::mutex mtx_;
};

}  // namespace sympp::testing
