# -----------------------------------------------------------------------------
# sympp library
# -----------------------------------------------------------------------------
add_library(sympp)
add_library(SymPP::sympp ALIAS sympp)

target_sources(sympp
    PRIVATE
        core/add.cpp
        core/assumption_mask.cpp
        core/basic.cpp
        core/boolean.cpp
        core/expand.cpp
        core/expr_cache.cpp
        core/float.cpp
        core/function.cpp
        core/imaginary_unit.cpp
        core/integer.cpp
        core/mul.cpp
        core/number_arith.cpp
        core/number_symbol.cpp
        core/piecewise.cpp
        core/pow.cpp
        core/queries.cpp
        core/rational.cpp
        core/refine.cpp
        core/singletons.cpp
        core/symbol.cpp
        core/traversal.cpp
        core/undefined_function.cpp
        calculus/asymptotes.cpp
        calculus/diff.cpp
        calculus/euler_lagrange.cpp
        calculus/extrema.cpp
        calculus/limit.cpp
        calculus/order.cpp
        calculus/pade.cpp
        calculus/series.cpp
        calculus/summation.cpp
        functions/combinatorial.cpp
        functions/exponential.cpp
        functions/hyperbolic.cpp
        functions/integers.cpp
        functions/miscellaneous.cpp
        functions/special.cpp
        functions/trigonometric.cpp
        integrals/integrate.cpp
        integrals/quadrature.cpp
        integrals/transforms.cpp
        matrices/matrix.cpp
        matrices/matrix_symbol.cpp
        polys/poly.cpp
        polys/rootof.cpp
        simplify/simplify.cpp
        solvers/solve.cpp
)

target_sources(sympp
    PUBLIC FILE_SET HEADERS
        BASE_DIRS ${CMAKE_SOURCE_DIR}/include
        FILES
            ${CMAKE_SOURCE_DIR}/include/sympp/sympp.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/fwd.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/version.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/add.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/api.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/assumption_key.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/assumption_mask.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/basic.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/boolean.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/expand.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/expr_cache.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/expr_collections.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/float.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/function.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/function_id.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/imaginary_unit.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/integer.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/mul.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/number.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/number_arith.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/number_symbol.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/operators.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/piecewise.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/pow.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/queries.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/refine.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/rational.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/singletons.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/symbol.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/traversal.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/type_id.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/undefined_function.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/calculus/asymptotes.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/calculus/diff.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/calculus/euler_lagrange.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/calculus/extrema.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/calculus/limit.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/calculus/order.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/calculus/pade.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/calculus/series.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/calculus/summation.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/functions/combinatorial.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/functions/exponential.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/functions/hyperbolic.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/functions/integers.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/functions/miscellaneous.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/functions/special.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/functions/trigonometric.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/integrals/integrate.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/integrals/quadrature.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/integrals/transforms.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/matrices/matrix.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/matrices/matrix_symbol.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/polys/poly.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/polys/rootof.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/simplify/simplify.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/solvers/solve.hpp
)

# When building the library itself we need to export symbols (Windows DLL).
# Consumers see SYMPP_EXPORT as dllimport.
target_compile_definitions(sympp PRIVATE SYMPP_BUILDING_DLL)
if(NOT SYMPP_BUILD_SHARED)
    target_compile_definitions(sympp PUBLIC SYMPP_STATIC)
endif()

target_include_directories(sympp
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(sympp PUBLIC cxx_std_20)

target_link_libraries(sympp
    PUBLIC
        GMP::gmpxx
        GMP::gmp
        MPFR::mpfr
    PRIVATE
        SymPP::warnings
        SymPP::sanitizers
)

set_target_properties(sympp PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    POSITION_INDEPENDENT_CODE ON
)

# Generate version.hpp from version.hpp.in
configure_file(
    ${CMAKE_SOURCE_DIR}/include/sympp/version.hpp.in
    ${CMAKE_SOURCE_DIR}/include/sympp/version.hpp
    @ONLY
)
