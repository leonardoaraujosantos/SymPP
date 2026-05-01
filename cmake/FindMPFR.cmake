# FindMPFR.cmake — locate the GNU MPFR library.
#
# Defines:
#   MPFR_FOUND
#   MPFR_INCLUDE_DIRS
#   MPFR_LIBRARIES
#   MPFR_VERSION
#
# Imported target:
#   MPFR::mpfr  (links GMP::gmp transitively)

include(FindPackageHandleStandardArgs)

find_path(MPFR_INCLUDE_DIR
    NAMES mpfr.h
    HINTS
        ENV MPFR_DIR
        ${CMAKE_PREFIX_PATH}
        /usr/local
        /opt/homebrew
        /opt/local
    PATH_SUFFIXES include
)

find_library(MPFR_LIBRARY
    NAMES mpfr
    HINTS
        ENV MPFR_DIR
        ${CMAKE_PREFIX_PATH}
        /usr/local
        /opt/homebrew
        /opt/local
    PATH_SUFFIXES lib
)

if(MPFR_INCLUDE_DIR AND EXISTS "${MPFR_INCLUDE_DIR}/mpfr.h")
    file(STRINGS "${MPFR_INCLUDE_DIR}/mpfr.h" mpfr_version_line REGEX "^#define[ \t]+MPFR_VERSION_STRING[ \t]+\".*\"")
    string(REGEX REPLACE "^#define[ \t]+MPFR_VERSION_STRING[ \t]+\"([^\"]+)\".*" "\\1" MPFR_VERSION "${mpfr_version_line}")
endif()

find_package_handle_standard_args(MPFR
    REQUIRED_VARS MPFR_INCLUDE_DIR MPFR_LIBRARY
    VERSION_VAR MPFR_VERSION
)

if(MPFR_FOUND)
    set(MPFR_INCLUDE_DIRS ${MPFR_INCLUDE_DIR})
    set(MPFR_LIBRARIES ${MPFR_LIBRARY})

    if(NOT TARGET MPFR::mpfr)
        add_library(MPFR::mpfr UNKNOWN IMPORTED)
        set_target_properties(MPFR::mpfr PROPERTIES
            IMPORTED_LOCATION "${MPFR_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR}"
        )
        if(TARGET GMP::gmp)
            set_property(TARGET MPFR::mpfr APPEND PROPERTY
                INTERFACE_LINK_LIBRARIES GMP::gmp
            )
        endif()
    endif()
endif()

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)
