/**
 * @brief common definitions we want to be widely available
 *
 */
/*
 * Original Author: Bevin Brett
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#pragma once

#if !defined(DARWIN) && !defined(__APPLE__)
#include <malloc.h>
// DON'T #include <mm_malloc.h> IT BREAKS icc AND IS NOT NEEDED
#endif

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <time.h>

#include "proto.h"  // malloc   is needed from here

extern const char *Progname;
    // various of the utility programs define this global variable to be their name

// you can do #if GCC_VERSION > 30200 for gcc 3.2.0
#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)
// from Graham Wideman:
// __func__ in gcc appears not to be a predefined macro that can be
// tested and also glued to string literals. Instead we have to use
// printf-style
#define __MYFUNCTION__ __func__

#ifdef  __FUNCTION__
#undef  __MYFUNCTION__
#define __MYFUNCTION__  __FUNCTION__
#endif

#ifdef  __FUNC__
#undef  __MYFUNCTION__
#define __MYFUNCTION__  __FUNC__
#endif

#endif

// Causes the compiler to not inline this function
//
#define NOINLINE __attribute__((noinline))

// defines the maximum number of threads used in OpenMP code
//
#define _MAX_FS_THREADS 128 

// assertions
//
void assertFailed(const char* file, int line, const char* tst);
#define cheapAssert(TST)        { if (!(TST)) assertFailed(__FILE__, __LINE__, #TST); }
#define costlyAssert(TST) //    { if (!(TST)) assertFailed(__FILE__, __LINE__, #TST); }

#define cheapAssertValidFno(_MRIS, _FNO) cheapAssert((0 <= _FNO) && (_FNO < _MRIS->nfaces))
#define cheapAssertValidVno(_MRIS, _VNO) cheapAssert((0 <= _VNO) && (_VNO < _MRIS->nvertices))

typedef enum LogicProblemResponse {
  LogicProblemResponse_old,
  LogicProblemResponse_fix
} LogicProblemResponse;

bool spendTimeCheckingForLogicProblem(const char* file, int line);
LogicProblemResponse copeWithLogicProblem2(
    bool* wasReported,          // set to whether reported
    const char* envvarFixer,    // user told to set this env var to fix, as part of the report
    const char* msg,            // the report
    const char* file, int line, // multiple reports at the same location are partially elided
    const char* function);
#define copeWithLogicProblem(ENV, MSG) copeWithLogicProblem2(NULL, (ENV), (MSG), __FILE__, __LINE__, __MYFUNCTION__)

// Regardless of whether the __real_malloc etc. or the __wrap_ ones, it is still desirable
// to know where in the program the allocations are happening.  This mechanism allows that to happen.
//
//#define DEBUG_MEMLEAK

#if defined(DEBUG_MEMLEAK)

void *mallocHere (              size_t size,                        const char* file, const char* function, int line);
void  freeHere   (void *ptr,                                        const char* file, const char* function, int line);
void* callocHere (size_t nmemb, size_t size,                        const char* file, const char* function, int line);
void *reallocHere(void *ptr,    size_t size,                        const char* file, const char* function, int line);
int posix_memalignHere(void **memptr, size_t alignment, size_t size,const char* file, const char* function, int line);
    // implemented in mgh_malloc.c but defined here to get them widely used

#define malloc(SIZE)                        mallocHere ((SIZE),                         __FILE__, __MYFUNCTION__, __LINE__)
#define free(PTR)                           freeHere   ((PTR),                          __FILE__, __MYFUNCTION__, __LINE__)
#define calloc(NMEMB,SIZE)                  callocHere ((NMEMB), (SIZE),                __FILE__, __MYFUNCTION__, __LINE__)
#define realloc(PTR,SIZE)                   reallocHere((PTR),   (SIZE),                __FILE__, __MYFUNCTION__, __LINE__)
#define posix_memalign(MEMPTR,ALIGN,SIZE)   posix_memalignHere((MEMPTR),(ALIGN),(SIZE), __FILE__, __MYFUNCTION__, __LINE__)

#endif

//#define freeAndNULL(PTR) { free((void*)(PTR)); (PTR) = NULL;}
#define freeAndNULL(PTR) { if((PTR)!=NULL){free((void*)(PTR)); (PTR) = NULL;} }


// Some trivial types
//
#ifndef uchar
typedef unsigned char unsigned_char;
#define uchar unsigned_char
#endif

typedef const float * ptr_to_const_float;


// Some trivial math functions needed lots
//
#pragma GCC diagnostic ignored "-Wunused-function"
static float squaref(float x) { return x*x; }

#pragma GCC diagnostic ignored "-Wunused-function"
static double squared(double x) { return x*x; }

typedef struct FloatXYZ {
    float x,y,z;
} FloatXYZ;



template <typename T, size_t SIZE>
struct FixedSizeArray {
    T&         operator[](size_t i)         { return v[i]; }
    T const&   operator[](size_t i)   const { return v[i]; }
    
    operator T       *()                    { return v;    }
    operator T const *()              const { return v;    }
    
    T       * data()                        { return v;    }
    T const * data()                  const { return v;    }
    
private:
    T v[SIZE];
};
