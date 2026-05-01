#pragma once

// Symbol-visibility / DLL-import macros for SymPP's public surface.
//
// Usage: mark every public class declaration with SYMPP_EXPORT and every
// non-template public free function the same way. Templates and inline
// definitions don't need it.

#if defined(_WIN32) || defined(__CYGWIN__)
#    if defined(SYMPP_BUILDING_DLL)
#        define SYMPP_EXPORT __declspec(dllexport)
#    elif defined(SYMPP_STATIC)
#        define SYMPP_EXPORT
#    else
#        define SYMPP_EXPORT __declspec(dllimport)
#    endif
#    define SYMPP_LOCAL
#else
#    define SYMPP_EXPORT __attribute__((visibility("default")))
#    define SYMPP_LOCAL __attribute__((visibility("hidden")))
#endif
