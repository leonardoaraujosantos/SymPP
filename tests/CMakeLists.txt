# -----------------------------------------------------------------------------
# Test executable
# -----------------------------------------------------------------------------
add_executable(sympp_tests)

target_sources(sympp_tests
    PRIVATE
        oracle/oracle.hpp
        oracle/oracle.cpp
        calculus/diff_test.cpp
        calculus/product_test.cpp
        calculus/series_limit_test.cpp
        core/arithmetic_test.cpp
        core/assumptions_test.cpp
        core/satask_test.cpp
        core/rewrite_test.cpp
        combinatorics/combinatorics_test.cpp
        crypto/crypto_test.cpp
        discrete/discrete_test.cpp
        core/lambdify_test.cpp
        core/lambdify_jit_test.cpp
        core/parity_probe.cpp
        core/collection_test.cpp
        core/expand_test.cpp
        core/expr_cache_test.cpp
        core/float_test.cpp
        core/imaginary_unit_test.cpp
        core/infinity_test.cpp
        core/piecewise_test.cpp
        core/integer_test.cpp
        core/rational_test.cpp
        core/singletons_test.cpp
        core/symbol_test.cpp
        core/traversal_test.cpp
        core/undefined_function_test.cpp
        functions/bessel_test.cpp
        functions/combinatorial_test.cpp
        functions/ntheory_test.cpp
        functions/orthopolys_test.cpp
        functions/exponential_test.cpp
        functions/hyperbolic_test.cpp
        functions/hypergeometric_test.cpp
        functions/integers_test.cpp
        functions/inverse_trig_test.cpp
        functions/miscellaneous_test.cpp
        functions/special_test.cpp
        functions/special_integrals_test.cpp
        functions/trigonometric_test.cpp
        integrals/integrate_test.cpp
        integrals/transforms_test.cpp
        integrals/meijerint_test.cpp
        integrals/hankel_test.cpp
        geometry/geometry_test.cpp
        logic/logic_test.cpp
        physics/physics_test.cpp
        physics/physics_quantum_test.cpp
        physics/physics_fock_test.cpp
        physics/physics_multimode_fock_test.cpp
        physics/physics_fermion_test.cpp
        physics/physics_normal_order_test.cpp
        stats/stats_test.cpp
        tensor/tensor_test.cpp
        vector/vector_calculus_test.cpp
        matlab/matlab_test.cpp
        matlab/matlab_assumptions_test.cpp
        matlab/matlab_ode_test.cpp
        matlab/matlab_parsing_test.cpp
        matlab/matlab_pde_test.cpp
        matlab/matlab_solvers_test.cpp
        matrices/matrix_test.cpp
        matrices/jordan_general_test.cpp
        matrices/svd_test.cpp
        ode/dsolve_test.cpp
        parsing/parser_test.cpp
        parsing/latex_parser_test.cpp
        physics/units_test.cpp
        printing/printing_test.cpp
        printing/pretty_test.cpp
        codegen/codegen_test.cpp
        oracle/oracle_smoke_test.cpp
        oracle/oracle_crosscheck_test.cpp
        polys/poly_test.cpp
        polys/domain_matrix_test.cpp
        polys/zassenhaus_test.cpp
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
    SYMPP_TEST_PARITY_DATA="${CMAKE_SOURCE_DIR}/tools/parity_expected.tsv"
)

include(Catch)
catch_discover_tests(sympp_tests)
