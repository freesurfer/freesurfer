#pragma once
/*
  The following are taken from boost/predef
  Ideally we would use that library, but it only appeared in boost 1.55 (Nov 2013)
 */

#ifdef BOOST_VERSION
#error "Please use boost/predef in place of FS_COMP_* macros"
#endif

#undef FS_COMP_CLANG
#undef FS_COMP_GNUC

#if defined(__clang__)
#define FS_COMP_CLANG "Clang detected"
#endif

// Clang emulates GCC, so need to detect it first

#ifndef FS_COMP_CLANG
#if defined(__GNUC__)
#define FS_COMP_GNUC "GCC Detected"
#endif
#endif
