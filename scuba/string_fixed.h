#ifndef string_fixed_h
#define string_fixed_h

#include <ansidecl.h>


#if (GCC_VERSION <= 2096)
// pick particular version of alloc.h!
#include "/usr/include/g++-3/alloc.h"
#endif

// then
#include <string>


#endif
