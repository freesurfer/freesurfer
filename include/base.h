/**
 * @file  base.h
 * @brief common definitions we want to be widely available
 *
 */
/*
 * Original Author: Bevin Brett
 * CVS Revision Info:
 *    $Author: bbrett $
 *    $Date: 2018/03/14 18:28:00 $
 *    $Revision: 1.0 $
 *
 * Copyright © 2018 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <mm_malloc.h>
#endif

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <time.h>

#include "proto.h"  // malloc   is needed from here


#if defined(__cplusplus)
extern "C" {
#endif



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


// defines the maximum number of threads used in OpenMP code
//
#define _MAX_FS_THREADS 128 


// Regardless of whether the __real_malloc etc. or the __wrap_ ones, it is still desirable
// to know where in the program the allocations are happening.  This mechanism allows that to happen.
//
#if 1

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


#if defined(__cplusplus)
};
#endif

