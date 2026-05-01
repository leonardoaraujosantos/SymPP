# SymPP examples — small standalone consumer programs.
#
# When built from the SymPP root tree (SYMPP_BUILD_EXAMPLES=ON), each
# example links against the in-tree `sympp` target. The same CMakeLists
# also works as a standalone consumer project — just point it at an
# installed SymPP via find_package.

if(NOT TARGET SymPP::sympp)
    # Standalone build path.
    cmake_minimum_required(VERSION 3.25)
    project(sympp_examples CXX)
    set(CMAKE_CXX_STANDARD 20)
    find_package(SymPP REQUIRED)
endif()

set(_examples
    01_calculus
    02_polynomials
    03_linear_algebra
    04_solving
    05_codegen
)

foreach(ex IN LISTS _examples)
    add_executable(${ex} ${ex}.cpp)
    target_link_libraries(${ex} PRIVATE SymPP::sympp)
endforeach()
