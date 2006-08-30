#if !defined(__MCONF_H)
#define __MCONF_H

#ifdef __cplusplus
extern "C" {
#endif



/*                                                      mconf.h
 *
 *     Common include file for math routines
 *
 *
 *
 * SYNOPSIS:
 *
 * #include "mconf.h"
 *
 *
 *
 * DESCRIPTION:
 *
 * This file contains definitions for error codes that are
 * passed to the common error handling routine mtherr()
 * (which see).
 *
 * The file also includes a conditional assembly definition
 * for the type of computer arithmetic (IEEE, DEC, Motorola
 * IEEE, or UNKnown).
 *
 * For Digital Equipment PDP-11 and VAX computers, certain
 * IBM systems, and others that use numbers with a 56-bit
 * significand, the symbol DEC should be defined.  In this
 * mode, most floating point constants are given as arrays
 * of octal integers to eliminate decimal to binary conversion
 * errors that might be introduced by the compiler.
 *
 * For computers, such as IBM PC, that follow the IEEE
 * Standard for Binary Floating Point Arithmetic (ANSI/IEEE
 * Std 754-1985), the symbol IBMPC should be defined.  These
 * numbers have 53-bit significands.  In this mode, constants
 * are provided as arrays of hexadecimal 16 bit integers.
 *
 * To accommodate other types of computer arithmetic, all
 * constants are also provided in a normal decimal radix
 * which one can hope are correctly converted to a suitable
 * format by the available C language compiler.  To invoke
 * this mode, the symbol UNK is defined.
 *
 * An important difference among these modes is a predefined
 * set of machine arithmetic constants for each.  The numbers
 * MACHEP (the machine roundoff error), MAXNUM (largest number
 * represented), and several other parameters are preset by
 * the configuration symbol.  Check the file const.c to
 * ensure that these values are correct for your computer.
 *
 */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1987, 1989, 1995 by Stephen L. Moshier
*/


/* Constant definitions for math error conditions
 */

#ifndef DOMAIN
#define DOMAIN          1       /* argument domain error */
#endif
#ifndef SING
#define SING            2       /* argument singularity */
#endif
#ifndef OVERFLOW
#define OVERFLOW        3       /* overflow range error */
#endif
#ifndef UNDERFLOW
#define UNDERFLOW       4       /* underflow range error */
#endif
#ifndef TLOSS
#define TLOSS           5       /* total loss of precision */
#endif
#ifndef PLOSS
#define PLOSS           6       /* partial loss of precision */
#endif

#define EDOM            33
#define ERANGE          34

/* Complex numeral.  */
typedef struct
       {
       double r;
       double i;
       }cmplx;

/* Long double complex numeral.  */
typedef struct
       {
       double r;
       double i;
       }cmplxl;

/* Type of computer arithmetic */

/* PDP-11, Pro350, VAX:
 */
/*define DEC 1*/

/* Intel IEEE, low order words come first:
#define IBMPC 1
 */

/* Motorola IEEE, high order words come first
 * (Sun 680x0 workstation):
 */
/*#define MIEEE 1*/

/* UNKnown arithmetic, invokes coefficients given in
 * normal decimal format.  Beware of range boundary
 * problems (MACHEP, MAXLOG, etc. in const.c) and
 * roundoff problems in pow.c:
 * (Sun SPARCstation)
 */
#define UNK 1

/* If you define UNK, then be sure to set BIGENDIAN properly. */
#define BIGENDIAN 0

/* Define this `volatile' if your compiler thinks
 * that floating point arithmetic obeys the associative
 * and distributive laws.  It will defeat some optimizations
 * (but probably not enough of them).
 *
 * #define VOLATILE volatile
 */
/*#define VOLATILE volatile*/

/* For 12-byte long doubles on an i386, pad a 16-bit short 0
 * to the end of real constants initialized by integer arrays.
 *
 * #define XPD 0,
 *
 * Otherwise, the type is 10 bytes long and XPD should be
 * defined blank (e.g., Microsoft C).
 *
 * #define XPD
 */
#define XPD 0,


/* Define to ask for infinity support, else undefine. */
/* #define INFINITIES 1 */

/* Define to ask for support of numbers that are Not-a-Number,
   else undefine.  This may automatically define INFINITY in some files. */
/* #define NANS 1 */

/* Define to support tiny denormal numbers, else undefine. */
#define DENORMAL 1

/* Define to distinguish between -0.0 and +0.0.  */
/* #define MINUSZERO 1 */

/* Define 1 for ANSI C atan2() function
   See atan.c and clog.c. */
#define ANSIC 1

/* Get ANSI function prototypes, if you want them. */
#ifdef __STDC__
#define ANSIPROT
#include "protos.h"
#else
int mtherr();
#endif

/* Variable for error reporting.  See mtherr.c.  */
extern int merror;


#ifdef __cplusplus
}
#endif

#endif
