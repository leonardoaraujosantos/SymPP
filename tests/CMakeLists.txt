# -----------------------------------------------------------------------------
# Test executable
# -----------------------------------------------------------------------------
add_executable(sympp_tests)

target_sources(sympp_tests
    PRIVATE
        oracle/oracle.hpp
        oracle/oracle.cpp
        core/symbol_test.cpp
        oracle/oracle_smoke_test.cpp
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
