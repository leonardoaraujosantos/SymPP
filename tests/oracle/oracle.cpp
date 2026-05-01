#include "oracle.hpp"

#include <array>
#include <cerrno>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#ifndef SYMPP_TEST_ORACLE_SCRIPT
#error "SYMPP_TEST_ORACLE_SCRIPT must be defined by the build system (path to oracle.py)"
#endif

#ifndef SYMPP_TEST_PYTHON_EXECUTABLE
#error "SYMPP_TEST_PYTHON_EXECUTABLE must be defined by the build system"
#endif

namespace sympp::testing {

// ---- OracleResponse -------------------------------------------------------

std::string OracleResponse::result_str() const {
    if (!raw.contains("result")) {
        return {};
    }
    const auto& r = raw.at("result");
    if (r.is_string()) return r.get<std::string>();
    return r.dump();
}

bool OracleResponse::result_bool() const {
    if (!raw.contains("result")) return false;
    const auto& r = raw.at("result");
    return r.is_boolean() && r.get<bool>();
}

std::string OracleResponse::error() const {
    return raw.value("error", std::string{});
}

std::string OracleResponse::detail() const {
    return raw.value("detail", std::string{});
}

// ---- Oracle ---------------------------------------------------------------

Oracle& Oracle::instance() {
    static Oracle instance;
    return instance;
}

Oracle::Oracle() {
    spawn();
}

Oracle::~Oracle() {
    shutdown_subprocess();
}

void Oracle::spawn() {
    int in_pipe[2];   // parent writes [1] -> child reads [0]
    int out_pipe[2];  // child writes  [1] -> parent reads [0]

    if (pipe(in_pipe) != 0 || pipe(out_pipe) != 0) {
        throw std::runtime_error(std::string{"pipe() failed: "} + std::strerror(errno));
    }

    pid_t pid = fork();
    if (pid < 0) {
        throw std::runtime_error(std::string{"fork() failed: "} + std::strerror(errno));
    }

    if (pid == 0) {
        // ---- Child ----
        // Wire stdin <- in_pipe[0], stdout -> out_pipe[1]
        dup2(in_pipe[0], STDIN_FILENO);
        dup2(out_pipe[1], STDOUT_FILENO);
        // Close all pipe fds in the child (after dup2)
        close(in_pipe[0]);
        close(in_pipe[1]);
        close(out_pipe[0]);
        close(out_pipe[1]);
        // Force unbuffered I/O on Python stdout
        setenv("PYTHONUNBUFFERED", "1", 1);

        const char* py = SYMPP_TEST_PYTHON_EXECUTABLE;
        const char* script = SYMPP_TEST_ORACLE_SCRIPT;
        // execvp argv must be writable char* -- duplicate
        std::array<char*, 3> argv{
            const_cast<char*>(py),
            const_cast<char*>(script),
            nullptr,
        };
        execvp(py, argv.data());
        // If we get here exec failed
        std::fprintf(stderr, "execvp failed for %s: %s\n", py, std::strerror(errno));
        std::_Exit(127);
    }

    // ---- Parent ----
    pid_ = pid;
    close(in_pipe[0]);
    close(out_pipe[1]);
    stdin_fd_ = in_pipe[1];
    stdout_fd_ = out_pipe[0];
}

void Oracle::shutdown_subprocess() noexcept {
    if (pid_ <= 0) return;
    try {
        // Best-effort polite shutdown
        nlohmann::json bye{{"op", "shutdown"}};
        auto payload = bye.dump() + "\n";
        ::write(stdin_fd_, payload.data(), payload.size());
    } catch (...) {
        // ignored — destructor must not throw
    }
    if (stdin_fd_ >= 0) close(stdin_fd_);
    if (stdout_fd_ >= 0) close(stdout_fd_);
    int status = 0;
    waitpid(pid_, &status, 0);
    pid_ = -1;
    stdin_fd_ = -1;
    stdout_fd_ = -1;
}

OracleResponse Oracle::send(nlohmann::json request) {
    std::lock_guard lock(mtx_);

    int id = next_id_++;
    request["id"] = id;
    auto payload = request.dump() + "\n";

    // ---- Write request ----
    const char* p = payload.data();
    std::size_t remaining = payload.size();
    while (remaining > 0) {
        ssize_t n = ::write(stdin_fd_, p, remaining);
        if (n < 0) {
            if (errno == EINTR) continue;
            throw std::runtime_error(std::string{"oracle write failed: "} + std::strerror(errno));
        }
        p += n;
        remaining -= static_cast<std::size_t>(n);
    }

    // ---- Read response (line-delimited) ----
    while (true) {
        auto nl = read_buf_.find('\n');
        if (nl != std::string::npos) {
            auto line = read_buf_.substr(0, nl);
            read_buf_.erase(0, nl + 1);

            OracleResponse resp;
            try {
                resp.raw = nlohmann::json::parse(line);
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string{"oracle returned invalid JSON: "} + e.what()
                                         + "\nline: " + line);
            }
            resp.ok = resp.raw.value("ok", false);
            return resp;
        }

        std::array<char, 4096> chunk{};
        ssize_t n = ::read(stdout_fd_, chunk.data(), chunk.size());
        if (n < 0) {
            if (errno == EINTR) continue;
            throw std::runtime_error(std::string{"oracle read failed: "} + std::strerror(errno));
        }
        if (n == 0) {
            throw std::runtime_error("oracle closed stdout unexpectedly");
        }
        read_buf_.append(chunk.data(), static_cast<std::size_t>(n));
    }
}

// ---- High-level helpers ---------------------------------------------------

std::string Oracle::srepr(std::string_view expr) {
    auto resp = send({{"op", "srepr"}, {"expr", expr}});
    if (!resp.ok) {
        throw std::runtime_error("oracle srepr failed: " + resp.error() + " — " + resp.detail());
    }
    return resp.result_str();
}

std::string Oracle::sympy_str(std::string_view expr) {
    auto resp = send({{"op", "parse"}, {"expr", expr}});
    if (!resp.ok) {
        throw std::runtime_error("oracle parse failed: " + resp.error() + " — " + resp.detail());
    }
    return resp.raw.value("str", std::string{});
}

bool Oracle::equivalent(std::string_view a, std::string_view b) {
    auto resp = send({{"op", "equivalent"}, {"a", a}, {"b", b}});
    if (!resp.ok) {
        throw std::runtime_error("oracle equivalent failed: " + resp.error() + " — "
                                 + resp.detail());
    }
    return resp.result_bool();
}

std::string Oracle::evalf(std::string_view expr, int prec) {
    auto resp = send({{"op", "evalf"}, {"expr", expr}, {"prec", prec}});
    if (!resp.ok) {
        throw std::runtime_error("oracle evalf failed: " + resp.error());
    }
    return resp.result_str();
}

std::string Oracle::diff(std::string_view expr, std::string_view var, int order) {
    auto resp = send({{"op", "diff"}, {"expr", expr}, {"var", var}, {"order", order}});
    if (!resp.ok) {
        throw std::runtime_error("oracle diff failed: " + resp.error());
    }
    return resp.result_str();
}

std::string Oracle::integrate(std::string_view expr, std::string_view var) {
    auto resp = send({{"op", "integrate"}, {"expr", expr}, {"var", var}});
    if (!resp.ok) {
        throw std::runtime_error("oracle integrate failed: " + resp.error());
    }
    return resp.result_str();
}

std::string Oracle::simplify(std::string_view expr) {
    auto resp = send({{"op", "simplify"}, {"expr", expr}});
    if (!resp.ok) {
        throw std::runtime_error("oracle simplify failed: " + resp.error());
    }
    return resp.result_str();
}

std::string Oracle::sympy_version() {
    auto resp = send({{"op", "ping"}});
    if (!resp.ok) {
        throw std::runtime_error("oracle ping failed");
    }
    return resp.raw.value("sympy_version", std::string{});
}

}  // namespace sympp::testing
