# SymPP task runner — https://github.com/casey/just
#
#   just            # list recipes
#   just build      # configure + compile (library, tests, examples)
#   just test       # full suite (unit + integration/oracle)
#   just test-unit  # pure-C++ unit tests (no SymPy oracle)
#   just examples   # build & run every example
#
# Override defaults inline, e.g.:  just build_type=Debug test
#   build_dir   - CMake build directory (default: build)
#   build_type  - CMAKE_BUILD_TYPE     (default: Release)
#   gmp_lib     - extra -L/-rpath for GMP/MPFR if not on the default path

build_dir  := "build"
build_type := "Release"
gmp_lib    := ""

# Show the recipe list.
default:
    @just --list

# ---- build -----------------------------------------------------------------

# Configure CMake (idempotent; enables tests + examples).
configure:
    cmake -S . -B {{build_dir}} \
      -DCMAKE_BUILD_TYPE={{build_type}} \
      -DSYMPP_BUILD_TESTS=ON \
      -DSYMPP_BUILD_EXAMPLES=ON

# Compile the library, tests and examples.
build: configure
    cmake --build {{build_dir}} --parallel

# Alias for `build`.
compile: build

# Wipe the build directory and rebuild from scratch.
rebuild: clean build

# ---- tests -----------------------------------------------------------------

# Full test suite (unit + integration / SymPy-oracle).
test: build
    {{build_dir}}/tests/sympp_tests

# Unit tests only — pure C++, no SymPy oracle needed.
test-unit: build
    {{build_dir}}/tests/sympp_tests "~[oracle]"

# Integration tests — cross-checked against the SymPy oracle (needs python3 + sympy).
test-integration: build
    {{build_dir}}/tests/sympp_tests "[oracle]"

# Run a Catch2 tag/name filter, e.g.  just test-filter "[polys]"
test-filter pattern: build
    {{build_dir}}/tests/sympp_tests "{{pattern}}"

# Python-side oracle harness tests (protocol / clean shutdown).
test-oracle-py:
    python3 -m unittest discover -s oracle -p 'test_*.py' -v

# ---- examples --------------------------------------------------------------

# Build and run every example (quick-start + engineering-course series).
examples: build
    #!/usr/bin/env bash
    set -euo pipefail
    for exe in $(find {{build_dir}}/examples -maxdepth 2 -type f -executable | sort); do
        echo "========== $(basename "$exe") =========="
        "$exe"
        echo
    done

# Run one example by target name, e.g.  just example calculus_03_integration
example name: build
    #!/usr/bin/env bash
    set -euo pipefail
    exe=$(find {{build_dir}}/examples -type f -executable -name "{{name}}" | head -1)
    if [ -z "$exe" ]; then echo "no example named '{{name}}'"; exit 1; fi
    "$exe"

# ---- parity gate -----------------------------------------------------------

# Build & run the Tier 1 + Tier 2 acceptance probe (must print REMAINING=0).
parity: build
    #!/usr/bin/env bash
    set -euo pipefail
    cache={{build_dir}}/CMakeCache.txt
    json=$(find {{build_dir}}/_deps -name json.hpp -path '*nlohmann*' 2>/dev/null | head -1 | xargs -r dirname | xargs -r dirname)
    # GMP/MPFR include + lib dirs as CMake found them (portable across hosts).
    gmpinc=$(grep -E '^GMP_INCLUDE_DIR:' "$cache" | cut -d= -f2)
    gmpdir=$(dirname "$(grep -E '^GMP_LIBRARY:' "$cache" | cut -d= -f2)")
    incflag=""; [ -n "$gmpinc" ] && incflag="-I$gmpinc"
    g++ -std=c++20 -O2 -I include -I tests -I "$json" $incflag \
        -DPROBE_GROUP -DPROBE_PHYS -DPROBE_PERF \
        -DSYMPP_TEST_ORACLE_SCRIPT="\"$PWD/oracle/oracle.py\"" \
        -DSYMPP_TEST_PYTHON_EXECUTABLE="\"python3\"" \
        parity-workspace/acceptance_probe.cpp tests/oracle/oracle.cpp \
        -o {{build_dir}}/parity_probe {{build_dir}}/src/libsympp.so \
        -L"$gmpdir" -Wl,-rpath,"$gmpdir" {{gmp_lib}} -lgmpxx -lgmp -lmpfr
    LD_LIBRARY_PATH={{build_dir}}/src:"$gmpdir"${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} {{build_dir}}/parity_probe | tail -6

# ---- extras ----------------------------------------------------------------

# Build and run the benchmark harness (Release recommended).
bench:
    cmake -S . -B {{build_dir}} -DCMAKE_BUILD_TYPE=Release -DSYMPP_BUILD_BENCHMARKS=ON
    cmake --build {{build_dir}} --target sympp_bench --parallel
    {{build_dir}}/benchmarks/sympp_bench

# Generate Doxygen API docs (output under build/docs/html; needs Doxygen).
docs:
    #!/usr/bin/env bash
    set -euo pipefail
    if ! command -v doxygen >/dev/null 2>&1; then
        echo "error: Doxygen not found — 'just docs' needs it to build the API docs."
        echo "  install: apt-get install doxygen graphviz   # Debian/Ubuntu"
        echo "           brew install doxygen graphviz       # macOS"
        exit 1
    fi
    cmake -S . -B {{build_dir}} -DCMAKE_BUILD_TYPE={{build_type}} \
        -DSYMPP_BUILD_TESTS=ON -DSYMPP_BUILD_EXAMPLES=ON -DSYMPP_BUILD_DOCS=ON
    cmake --build {{build_dir}} --target sympp_docs
    echo "API docs generated at {{build_dir}}/docs/html/index.html"

# Install the library + CMake package config to a prefix.
install prefix="/usr/local": build
    cmake --install {{build_dir}} --prefix {{prefix}}

# Format all C++ sources in place (needs clang-format; uses .clang-format).
format:
    find include src tests examples benchmarks -name '*.hpp' -o -name '*.cpp' \
      | xargs clang-format -i

# Validate the OpenSpec change/spec set.
spec-validate:
    openspec validate --strict

# List active OpenSpec changes.
spec:
    openspec list

# Configure + build + full test suite — the CI pipeline.
ci: build test

# Remove the build directory.
clean:
    rm -rf {{build_dir}}
