#ifndef  MINC_BASIC_HEADER_FILE
#define  MINC_BASIC_HEADER_FILE

/* ----------------------------- MNI Header -----------------------------------
@NAME       : minc_basic.h
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Constants and macros for private use by MINC routines.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 28, 1992 (Peter Neelin)
@MODIFIED   : 
 * $Log: minc_basic.h,v $
 * Revision 6.2  2001/04/17 18:40:13  neelin
 * Modifications to work with NetCDF 3.x
 * In particular, changed NC_LONG to NC_INT (and corresponding longs to ints).
 * Changed NC_UNSPECIFIED to NC_NAT.
 * A few fixes to the configure script.
 *
 * Revision 6.1  1999/10/19 14:45:08  neelin
 * Fixed Log subsitutions for CVS
 *
 * Revision 6.0  1997/09/12 13:24:54  neelin
 * Release of minc version 0.6
 *
 * Revision 5.0  1997/08/21  13:25:53  neelin
 * Release of minc version 0.5
 *
 * Revision 4.0  1997/05/07  20:07:52  neelin
 * Release of minc version 0.4
 *
 * Revision 3.0  1995/05/15  19:33:12  neelin
 * Release of minc version 0.3
 *
 * Revision 2.0  1994/09/28  10:38:01  neelin
 * Release of minc version 0.2
 *
 * Revision 1.8  94/09/28  10:37:26  neelin
 * Pre-release
 * 
 * Revision 1.7  93/10/28  10:18:23  neelin
 * Added FILLVALUE_EPSILON for doing fillvalue checking in icv's.
 * 
 * Revision 1.6  93/08/11  12:06:37  neelin
 * Added RCS logging in source.
 * 
@COPYRIGHT  :
              Copyright 1993 Peter Neelin, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
@RCSID      : $Header: /private-cvsroot/minc/libsrc/minc_basic.h,v 6.2 2001/04/17 18:40:13 neelin Exp $ MINC (MNI)
---------------------------------------------------------------------------- */

#include <math.h>

static long longMIN(long lhs, long rhs) {
   return (lhs < rhs) ? lhs : rhs;
}

static long longMAX(long lhs, long rhs) {
   return (lhs > rhs) ? lhs : rhs;
}

static double doubleMIN(double lhs, double rhs) {
   return (lhs < rhs) ? lhs : rhs;
}

static double doubleMAX(double lhs, double rhs) {
   return (lhs > rhs) ? lhs : rhs;
}

#define MI_PRIV_DEFSIGN   0
#define MI_PRIV_SIGNED    1
#define MI_PRIV_UNSIGNED  2

#define MI_PRIV_GET 10
#define MI_PRIV_PUT 11

/* Epsilon for detecting fillvalues */
#define FILLVALUE_EPSILON (10.0 * FLT_EPSILON)

#define MI_TO_DOUBLE(dvalue, type, sign, ptr) \
   switch (type) { \
   case NC_BYTE : \
      switch (sign) { \
      case MI_PRIV_UNSIGNED : \
         dvalue = (double) *((unsigned char *) ptr); break; \
      case MI_PRIV_SIGNED : \
         dvalue = (double) *((signed char *) ptr); break; \
      } \
      break; \
   case NC_SHORT : \
      switch (sign) { \
      case MI_PRIV_UNSIGNED : \
         dvalue = (double) *((unsigned short *) ptr); break; \
      case MI_PRIV_SIGNED : \
         dvalue = (double) *((signed short *) ptr); break; \
      } \
      break; \
   case NC_INT : \
      switch (sign) { \
      case MI_PRIV_UNSIGNED : \
         dvalue = (double) *((unsigned int *) ptr); break; \
      case MI_PRIV_SIGNED : \
         dvalue = (double) *((signed int  *) ptr); break; \
      } \
      break; \
   case NC_FLOAT : \
      dvalue = (double) *((float *) ptr); \
      break; \
   case NC_DOUBLE : \
      dvalue = (double) *((double *) ptr); \
      break; \
   default: \
      fprintf(stderr,"%s:%d no default\n", __FILE__,__LINE__); exit(1); \
   } 

#define MI_FROM_DOUBLE(dvalue, type, sign, ptr) \
   switch (type) { \
   case NC_BYTE : \
      switch (sign) { \
      case MI_PRIV_UNSIGNED : \
         dvalue = doubleMAX(0, dvalue); \
         dvalue = doubleMIN(UCHAR_MAX, dvalue); \
         *((unsigned char *) ptr) = round(dvalue); \
         break; \
      case MI_PRIV_SIGNED : \
         dvalue = doubleMAX(SCHAR_MIN, dvalue); \
         dvalue = doubleMIN(SCHAR_MAX, dvalue); \
         *((signed char *) ptr) = round(dvalue); \
         break; \
      } \
      break; \
   case NC_SHORT : \
      switch (sign) { \
      case MI_PRIV_UNSIGNED : \
         dvalue = doubleMAX(0, dvalue); \
         dvalue = doubleMIN(USHRT_MAX, dvalue); \
         *((unsigned short *) ptr) = round(dvalue); \
         break; \
      case MI_PRIV_SIGNED : \
         dvalue = doubleMAX(SHRT_MIN, dvalue); \
         dvalue = doubleMIN(SHRT_MAX, dvalue); \
         *((signed short *) ptr) = round(dvalue); \
         break; \
      } \
      break; \
   case NC_INT : \
      switch (sign) { \
      case MI_PRIV_UNSIGNED : \
         dvalue = doubleMAX(0, dvalue); \
         dvalue = doubleMIN(UINT_MAX, dvalue); \
         *((unsigned int *) ptr) = round(dvalue); \
         break; \
      case MI_PRIV_SIGNED : \
         dvalue = doubleMAX(INT_MIN, dvalue); \
         dvalue = doubleMIN(INT_MAX, dvalue); \
         *((signed int *) ptr) = round(dvalue); \
         break; \
      } \
      break; \
   case NC_FLOAT : \
      dvalue = doubleMAX(-FLT_MAX,dvalue); \
      *((float *) ptr) = doubleMIN(FLT_MAX,dvalue); \
      break; \
   case NC_DOUBLE : \
      *((double *) ptr) = dvalue; \
      break; \
   default: \
      fprintf(stderr,"%s:%d no default\n", __FILE__,__LINE__); exit(1); \
   }

#endif
