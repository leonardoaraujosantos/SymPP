# SymPP examples — calculus (engineering-course series).
# Built from the SymPP tree (SYMPP_BUILD_EXAMPLES=ON); also works standalone
# against an installed SymPP via find_package.

if(NOT TARGET SymPP::sympp)
    cmake_minimum_required(VERSION 3.25)
    project(sympp_calculus_examples CXX)
    set(CMAKE_CXX_STANDARD 20)
    find_package(SymPP REQUIRED)
endif()

file(GLOB _srcs CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
foreach(_src IN LISTS _srcs)
    get_filename_component(_name ${_src} NAME_WE)
    add_executable(calculus_${_name} ${_src})
    target_link_libraries(calculus_${_name} PRIVATE SymPP::sympp)
endforeach()
