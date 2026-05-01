# -----------------------------------------------------------------------------
# Test executable
# -----------------------------------------------------------------------------
add_executable(sympp_tests)

target_sources(sympp_tests
    PRIVATE
        oracle/oracle.hpp
        oracle/oracle.cpp
        calculus/diff_test.cpp
        calculus/series_limit_test.cpp
        core/arithmetic_test.cpp
        core/assumptions_test.cpp
        core/collection_test.cpp
        core/expand_test.cpp
        core/expr_cache_test.cpp
        core/float_test.cpp
        core/imaginary_unit_test.cpp
        core/piecewise_test.cpp
        core/integer_test.cpp
        core/rational_test.cpp
        core/singletons_test.cpp
        core/symbol_test.cpp
        core/traversal_test.cpp
        core/undefined_function_test.cpp
        functions/combinatorial_test.cpp
        functions/exponential_test.cpp
        functions/hyperbolic_test.cpp
        functions/integers_test.cpp
        functions/inverse_trig_test.cpp
        functions/miscellaneous_test.cpp
        functions/special_test.cpp
        functions/trigonometric_test.cpp
        integrals/integrate_test.cpp
        integrals/transforms_test.cpp
        matrices/matrix_test.cpp
        oracle/oracle_smoke_test.cpp
        polys/poly_test.cpp
        sets/sets_test.cpp
        simplify/simplify_test.cpp
        solvers/solve_test.cpp
)

target_include_directories(sympp_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(sympp_tests
    PRIVATE
        SymPP::sympp
        SymPP::warnings
        SymPP::sanitizers
        Catch2::Catch2WithMain
        nlohmann_json::nlohmann_json
)

# Catch2's REQUIRE(int_value) decomposer triggers -Wsign-conversion on
# integer expressions; suppress that one for tests only.
if(NOT MSVC)
    target_compile_options(sympp_tests PRIVATE -Wno-sign-conversion)
endif()

# Pass paths through to tests so the Oracle subprocess can find oracle.py
target_compile_definitions(sympp_tests PRIVATE
    SYMPP_TEST_ORACLE_SCRIPT="${CMAKE_SOURCE_DIR}/oracle/oracle.py"
    SYMPP_TEST_PYTHON_EXECUTABLE="${Python3_EXECUTABLE}"
)

include(Catch)
catch_discover_tests(sympp_tests)
