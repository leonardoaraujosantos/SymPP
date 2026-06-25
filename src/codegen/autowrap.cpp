#include <sympp/codegen/autowrap.hpp>

#include <dlfcn.h>
#include <unistd.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <sympp/printing/printing.hpp>

namespace sympp::codegen {

namespace {

// Resolve the C compiler: $CC if set and runnable, else first of cc/gcc/clang
// found via `command -v`. Returns empty string if none is usable.
std::string find_c_compiler() {
    std::vector<std::string> candidates;
    if (const char* cc = std::getenv("CC"); cc != nullptr && cc[0] != '\0') {
        candidates.emplace_back(cc);
    }
    candidates.insert(candidates.end(), {"cc", "gcc", "clang"});

    for (const auto& cand : candidates) {
        // `command -v` prints the resolved path and exits 0 when found.
        std::string probe = "command -v " + cand + " >/dev/null 2>&1";
        if (std::system(probe.c_str()) == 0) {
            return cand;
        }
    }
    return {};
}

// RAII deleter for a dlopen handle. The CompiledExpr captures a shared_ptr to
// keep the loaded shared object resident for the callable's lifetime.
struct DlHandle {
    void* handle = nullptr;
    explicit DlHandle(void* h) noexcept : handle(h) {}
    ~DlHandle() {
        if (handle != nullptr) {
            dlclose(handle);
        }
    }
    DlHandle(const DlHandle&) = delete;
    DlHandle& operator=(const DlHandle&) = delete;
};

// Unique temp path with the given suffix under TMPDIR (or /tmp).
std::string temp_path(const std::string& suffix) {
    const char* tmp = std::getenv("TMPDIR");
    std::string dir = (tmp != nullptr && tmp[0] != '\0') ? tmp : "/tmp";
    std::ostringstream os;
    os << dir << "/sympp_autowrap_" << ::getpid() << '_'
       << reinterpret_cast<std::uintptr_t>(&suffix) << '_'
       << std::rand() << suffix;
    return os.str();
}

}  // namespace

bool autowrap_available() noexcept {
    try {
        return !find_c_compiler().empty();
    } catch (...) {
        return false;
    }
}

CompiledExpr autowrap(const std::vector<Expr>& vars, const Expr& body) {
    const std::string cc = find_c_compiler();
    if (cc.empty()) {
        throw std::runtime_error(
            "autowrap: no C compiler found (set $CC or install cc/gcc/clang); "
            "gate on autowrap_available()");
    }

    const std::string fn_name = "sympp_autowrap_fn";
    const std::string call_name = "sympp_autowrap_call";

    // Emit the scalar-argument C function, plus a vector-argument trampoline so
    // we can call it uniformly through a const double* of arbitrary arity.
    std::ostringstream src;
    src << "#include <math.h>\n\n";
    src << printing::c_function(fn_name, body, vars) << "\n\n";
    src << "double " << call_name << "(const double* a) {\n";
    src << "    return " << fn_name << "(";
    for (std::size_t i = 0; i < vars.size(); ++i) {
        if (i > 0) src << ", ";
        src << "a[" << i << "]";
    }
    src << ");\n}\n";

    const std::string c_path = temp_path(".c");
    const std::string so_path = temp_path(".so");

    {
        std::ofstream out(c_path);
        if (!out) {
            throw std::runtime_error("autowrap: cannot write temp source " +
                                     c_path);
        }
        out << src.str();
    }

    // Compile into a shared object. Redirect compiler output to keep tests quiet.
    std::ostringstream cmd;
    cmd << cc << " -O2 -shared -fPIC -o " << so_path << ' ' << c_path
        << " -lm >/dev/null 2>&1";
    const int rc = std::system(cmd.str().c_str());
    if (rc != 0) {
        std::remove(c_path.c_str());
        throw std::runtime_error("autowrap: C compilation failed (compiler '" +
                                 cc + "')");
    }

    void* raw = dlopen(so_path.c_str(), RTLD_NOW | RTLD_LOCAL);
    // The handle pins the file in memory; we no longer need the on-disk copies.
    std::remove(c_path.c_str());
    std::remove(so_path.c_str());
    if (raw == nullptr) {
        throw std::runtime_error(std::string("autowrap: dlopen failed: ") +
                                 dlerror());
    }

    auto handle = std::make_shared<DlHandle>(raw);
    using CallFn = double (*)(const double*);
    // The double cast through void* silences the ISO C++ object/function-pointer
    // warning while remaining the portable POSIX idiom for dlsym.
    void* sym = dlsym(raw, call_name.c_str());
    if (sym == nullptr) {
        throw std::runtime_error("autowrap: dlsym failed for " + call_name);
    }
    auto fn = reinterpret_cast<CallFn>(reinterpret_cast<std::uintptr_t>(sym));

    const std::size_t arity = vars.size();
    return [handle, fn, arity](const std::vector<double>& xs) -> double {
        if (xs.size() < arity) {
            throw std::runtime_error("autowrap: too few arguments");
        }
        return fn(xs.data());
    };
}

CompiledExpr autowrap(const Expr& var, const Expr& body) {
    return autowrap(std::vector<Expr>{var}, body);
}

}  // namespace sympp::codegen
