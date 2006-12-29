/**
 * @file  macros.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.34 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


/*------------------------------------------------------------------------
      File Name:  macros.h

         Author:  Bruce Fischl

        Created:  Jan. 1994

    Description:

------------------------------------------------------------------------*/
#ifndef MACROS_H
#define MACROS_H

#include "const.h"
#include "utils.h"
#include "proto.h"
#include "math.h"

#ifdef _MSDOS
#include <math.h>
#define exp2(f)     pow(2.0,(f))
#define log2(f)     (log(f) / log(2.0))
#ifndef M_E
#define M_E 2.718282 /* exp(1) */
#endif
#ifndef M_PI
#define M_PI  3.141593
#endif

#endif
#ifndef UCHAR
#define UCHAR        unsigned char
#endif

#ifndef UINT
#define UINT         unsigned int
#endif

#define ISOPTION(c)  ((c) == '-')

/* these are defined in hipl_format.h */
#ifndef MIN
#define MIN(a,b)     ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b)     ((a) >= (b) ? (a) : (b))
#endif

#define EVEN(n)      ((((n) / 2) * 2) == n)
#define ODD(n)       (!EVEN(n))
#define ISEVEN       EVEN
#define ISODD        ODD

#define RADIANS(deg) ((2.0 * M_PI * (double)(deg)) / (360.0))
#define DEGREES(rad) ((360.0 * (double)(rad)) / (2.0 * M_PI))
#define NDEGREES(rad) (DEGREES(normAngle(rad)))

#define ISSMALL(f)   (fabs(f) < 0.000001f)
#define ISTINY(f)    (fabs(f) < 0.00000001f)

#define FZERO(f)     (fabs(f) < 0.0000001F)
#define DZERO(d)     (fabs(d) < 1e-15)
#define iszero(f)   (FZERO(f))
#define FEQUAL(f1,f2) (FZERO(f1-f2))
#define DEQUAL(d1,d2) (DZERO(d1-d2))

#define ISINT(f)      ((float)((int)f) == f)

#ifndef SQR
#define SQR(a)   ((a)*(a))
#endif

#ifndef ISIGN
#define ISIGN(a)   (a >= 0 ? 1 : -1)
#endif
#include <stdlib.h>
#include <string.h>

#define STRALLOC(str)   ((char *)calloc(strlen(str)+1, sizeof(char)))
#define STRCPALLOC(str) strcpy(STRALLOC(str), str)

#ifdef Linux
#define exp2(f)     pow(2.0,(f))
#define log2(f)     (log(f) / log(2.0))
#endif

#ifdef IRIX
#define exp2(f)     pow(2.0,(f))
#define log2(f)     (log(f) / log(2.0))
#endif

#ifdef SunOS
#define exp2(f)     pow(2.0,(f))
#define log2(f)     (log(f) / log(2.0))
#define ceilf(f)    (int)ceil((double)f)
#define floorf(f)   (int)floor((double)f)
#endif

#define ISPOW2(n)   (exp2((float)nint(log2((float)n))) == (float)n)

#ifdef INT32B
#undef INT32B
#endif

#define INT32B  long

#ifndef SGN
#define SGN(x)  (((x) < 0) ? -1 : 1)
#endif

#ifdef SunOS
#define memmove  memcpy
#endif

#endif

// you can do #if GCC_VERSION > 30200 for gcc 3.2.0
#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)
#endif
