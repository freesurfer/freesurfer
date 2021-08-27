/**
 * @brief common macro definitions
 *
 */
/*
 * Original Author: Bruce Fischl
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

#ifndef MACROS_H
#define MACROS_H

#include "base.h"

#include "const.h"
#include "utils.h"
#include "proto.h"
#include "math.h"

#define DEFINE_LOG2 static double log2(double x) { return log(x) / log(2.0); }
	// defining log2 as a macro has problems with other header files that also define it as a function taking an integer arg 
	
#ifdef _MSDOS
#include <math.h>
#define exp2(f)     pow(2.0,(f))
DEFINE_LOG2
#ifndef M_E
#define M_E 2.718282 /* exp(1) */
#endif
#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

#endif
#ifndef UCHAR
#define UCHAR        unsigned char
#endif

#ifndef UINT
#define UINT         unsigned int
#endif

#ifndef FALSE
#define FALSE  	     0
#endif

#ifndef TRUE
#define TRUE         1
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

#ifndef FLT_EPSILON
#define FLT_EPSILON 1e-5
#endif

#define FZERO(f)     (fabs(f) < FLT_EPSILON)
#define FZEROTHR(f,thr)     (fabs(f) < thr)
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
#endif

#ifdef IRIX
#define exp2(f)     pow(2.0,(f))
DEFINE_LOG2
#endif

#ifdef SunOS
#define exp2(f)     pow(2.0,(f))
DEFINE_LOG2
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
