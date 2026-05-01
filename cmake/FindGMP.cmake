# FindGMP.cmake — locate GNU Multiple Precision library.
#
# Defines:
#   GMP_FOUND
#   GMP_INCLUDE_DIRS
#   GMP_LIBRARIES
#   GMPXX_LIBRARIES   (C++ wrapper, libgmpxx)
#
# Imported targets:
#   GMP::gmp
#   GMP::gmpxx

include(FindPackageHandleStandardArgs)

find_path(GMP_INCLUDE_DIR
    NAMES gmp.h
    HINTS
        ENV GMP_DIR
        ${CMAKE_PREFIX_PATH}
        /usr/local
        /opt/homebrew
        /opt/local
    PATH_SUFFIXES include
)

find_library(GMP_LIBRARY
    NAMES gmp
    HINTS
        ENV GMP_DIR
        ${CMAKE_PREFIX_PATH}
        /usr/local
        /opt/homebrew
        /opt/local
    PATH_SUFFIXES lib
)

find_library(GMPXX_LIBRARY
    NAMES gmpxx
    HINTS
        ENV GMP_DIR
        ${CMAKE_PREFIX_PATH}
        /usr/local
        /opt/homebrew
        /opt/local
    PATH_SUFFIXES lib
)

# Extract version from gmp.h if available
if(GMP_INCLUDE_DIR AND EXISTS "${GMP_INCLUDE_DIR}/gmp.h")
    file(STRINGS "${GMP_INCLUDE_DIR}/gmp.h" gmp_version_line REGEX "^#define[ \t]+__GNU_MP_VERSION[ \t]+[0-9]+")
    string(REGEX REPLACE "^#define[ \t]+__GNU_MP_VERSION[ \t]+([0-9]+).*" "\\1" GMP_VERSION_MAJOR "${gmp_version_line}")
    file(STRINGS "${GMP_INCLUDE_DIR}/gmp.h" gmp_version_line REGEX "^#define[ \t]+__GNU_MP_VERSION_MINOR[ \t]+[0-9]+")
    string(REGEX REPLACE "^#define[ \t]+__GNU_MP_VERSION_MINOR[ \t]+([0-9]+).*" "\\1" GMP_VERSION_MINOR "${gmp_version_line}")
    file(STRINGS "${GMP_INCLUDE_DIR}/gmp.h" gmp_version_line REGEX "^#define[ \t]+__GNU_MP_VERSION_PATCHLEVEL[ \t]+[0-9]+")
    string(REGEX REPLACE "^#define[ \t]+__GNU_MP_VERSION_PATCHLEVEL[ \t]+([0-9]+).*" "\\1" GMP_VERSION_PATCH "${gmp_version_line}")
    set(GMP_VERSION "${GMP_VERSION_MAJOR}.${GMP_VERSION_MINOR}.${GMP_VERSION_PATCH}")
endif()

find_package_handle_standard_args(GMP
    REQUIRED_VARS GMP_INCLUDE_DIR GMP_LIBRARY
    VERSION_VAR GMP_VERSION
)

if(GMP_FOUND)
    set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
    set(GMP_LIBRARIES ${GMP_LIBRARY})
    if(GMPXX_LIBRARY)
        set(GMPXX_LIBRARIES ${GMPXX_LIBRARY})
    endif()

    if(NOT TARGET GMP::gmp)
        add_library(GMP::gmp UNKNOWN IMPORTED)
        set_target_properties(GMP::gmp PROPERTIES
            IMPORTED_LOCATION "${GMP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
        )
    endif()

    if(GMPXX_LIBRARY AND NOT TARGET GMP::gmpxx)
        add_library(GMP::gmpxx UNKNOWN IMPORTED)
        set_target_properties(GMP::gmpxx PROPERTIES
            IMPORTED_LOCATION "${GMPXX_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES GMP::gmp
        )
    endif()
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY)
