#ifndef string_fixed_h
#define string_fixed_h

#  define COMPILER_VERSION (__GNUC__ * 10000 \
                            + __GNUC_MINOR__ * 100 \
                            + __GNUC_PATCHLEVEL__)

#if (COMPILER_VERSION <= 2096)
// pick particular version of alloc.h!
#include "/usr/include/g++-3/alloc.h"
#endif

// then
#include <string>


#endif
