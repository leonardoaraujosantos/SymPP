# -----------------------------------------------------------------------------
# sympp library
# -----------------------------------------------------------------------------
add_library(sympp)
add_library(SymPP::sympp ALIAS sympp)

target_sources(sympp
    PRIVATE
        core/basic.cpp
        core/symbol.cpp
)

target_sources(sympp
    PUBLIC FILE_SET HEADERS
        BASE_DIRS ${CMAKE_SOURCE_DIR}/include
        FILES
            ${CMAKE_SOURCE_DIR}/include/sympp/sympp.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/fwd.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/version.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/api.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/basic.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/symbol.hpp
            ${CMAKE_SOURCE_DIR}/include/sympp/core/type_id.hpp
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
