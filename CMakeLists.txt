cmake_minimum_required(VERSION 3.25)

project(SymPP
    VERSION 0.0.1
    DESCRIPTION "Modern C++ symbolic mathematics library"
    HOMEPAGE_URL "https://github.com/leonardoaraujo/SymPP"
    LANGUAGES CXX
)

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------
option(SYMPP_BUILD_SHARED       "Build shared library"           ON)
option(SYMPP_BUILD_TESTS        "Build test suite"               ON)
option(SYMPP_BUILD_EXAMPLES     "Build examples"                 OFF)
option(SYMPP_BUILD_BENCHMARKS   "Build benchmarks"               OFF)
option(SYMPP_ENABLE_ASAN        "Enable AddressSanitizer"        OFF)
option(SYMPP_ENABLE_UBSAN       "Enable UndefinedBehaviorSanitizer" OFF)
option(SYMPP_WARNINGS_AS_ERRORS "Treat warnings as errors"       ON)

# -----------------------------------------------------------------------------
# Project-wide settings
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile_commands.json" FORCE)

# Force shared/static based on option
if(SYMPP_BUILD_SHARED)
    set(BUILD_SHARED_LIBS ON)
else()
    set(BUILD_SHARED_LIBS OFF)
endif()

# Hide internal symbols by default
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

# -----------------------------------------------------------------------------
# Compiler warnings
# -----------------------------------------------------------------------------
add_library(sympp_warnings INTERFACE)
add_library(SymPP::warnings ALIAS sympp_warnings)

if(MSVC)
    target_compile_options(sympp_warnings INTERFACE /W4 /permissive-)
    if(SYMPP_WARNINGS_AS_ERRORS)
        target_compile_options(sympp_warnings INTERFACE /WX)
    endif()
else()
    target_compile_options(sympp_warnings INTERFACE
        -Wall -Wextra -Wpedantic
        -Wshadow -Woverloaded-virtual
        -Wcast-align -Wunused -Wmisleading-indentation
        -Wnull-dereference -Wdouble-promotion -Wformat=2
        -Wconversion -Wsign-conversion
    )
    # -Wnon-virtual-dtor disabled project-wide: Catch2's internal expression
    # types trip it on every REQUIRE(). Our own discipline (every abstract
    # base declares `virtual ~T() = default`) is enforced by code review.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(sympp_warnings INTERFACE
            -Wlogical-op -Wduplicated-cond -Wduplicated-branches
        )
    endif()
    if(SYMPP_WARNINGS_AS_ERRORS)
        target_compile_options(sympp_warnings INTERFACE -Werror)
    endif()
endif()

# -----------------------------------------------------------------------------
# Sanitizers
# -----------------------------------------------------------------------------
add_library(sympp_sanitizers INTERFACE)
add_library(SymPP::sanitizers ALIAS sympp_sanitizers)

if(SYMPP_ENABLE_ASAN AND NOT MSVC)
    target_compile_options(sympp_sanitizers INTERFACE -fsanitize=address -fno-omit-frame-pointer)
    target_link_options(sympp_sanitizers INTERFACE -fsanitize=address)
endif()
if(SYMPP_ENABLE_UBSAN AND NOT MSVC)
    target_compile_options(sympp_sanitizers INTERFACE -fsanitize=undefined)
    target_link_options(sympp_sanitizers INTERFACE -fsanitize=undefined)
endif()

# -----------------------------------------------------------------------------
# Required dependencies
# -----------------------------------------------------------------------------
find_package(GMP REQUIRED)
find_package(MPFR 4.0 REQUIRED)

message(STATUS "  GMP:                  ${GMP_VERSION} (${GMP_LIBRARY})")
message(STATUS "  MPFR:                 ${MPFR_VERSION} (${MPFR_LIBRARY})")

# -----------------------------------------------------------------------------
# Subdirectories
# -----------------------------------------------------------------------------
add_subdirectory(src)

if(SYMPP_BUILD_TESTS)
    enable_testing()
    include(CTest)
    include(FetchContent)

    # Catch2 — test framework
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.5.4
    )

    # nlohmann/json — JSON for oracle harness protocol
    FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
    )

    FetchContent_MakeAvailable(Catch2 nlohmann_json)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

    # Locate Python — oracle harness runs it as a subprocess at test time
    find_package(Python3 3.10 COMPONENTS Interpreter REQUIRED)

    add_subdirectory(tests)
endif()

if(SYMPP_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
message(STATUS "")
message(STATUS "SymPP ${PROJECT_VERSION} configuration:")
message(STATUS "  Build type:           ${CMAKE_BUILD_TYPE}")
message(STATUS "  Compiler:             ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Shared lib:           ${SYMPP_BUILD_SHARED}")
message(STATUS "  Tests:                ${SYMPP_BUILD_TESTS}")
message(STATUS "  Examples:             ${SYMPP_BUILD_EXAMPLES}")
message(STATUS "  Benchmarks:           ${SYMPP_BUILD_BENCHMARKS}")
message(STATUS "  AddressSanitizer:     ${SYMPP_ENABLE_ASAN}")
message(STATUS "  UndefinedBehaviorSan: ${SYMPP_ENABLE_UBSAN}")
message(STATUS "  Warnings as errors:   ${SYMPP_WARNINGS_AS_ERRORS}")
message(STATUS "")
