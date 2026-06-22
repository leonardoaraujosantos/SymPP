# SymPP benchmarks — dependency-free wall-clock timing harness.
#
# Enabled with -DSYMPP_BUILD_BENCHMARKS=ON. Also builds standalone against an
# installed SymPP via find_package(SymPP).

if(NOT TARGET SymPP::sympp)
    cmake_minimum_required(VERSION 3.25)
    project(sympp_benchmarks CXX)
    set(CMAKE_CXX_STANDARD 20)
    find_package(SymPP REQUIRED)
endif()

add_executable(sympp_bench bench.cpp)
target_link_libraries(sympp_bench PRIVATE SymPP::sympp)
