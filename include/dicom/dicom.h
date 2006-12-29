/*
          Copyright (C) 1993, 1994, RSNA and Washington University

          The software and supporting documentation for the Radiological
          Society of North America (RSNA) 1993, 1994 Digital Imaging and
          Communications in Medicine (DICOM) Demonstration were developed
          at the
                  Electronic Radiology Laboratory
                  Mallinckrodt Institute of Radiology
                  Washington University School of Medicine
                  510 S. Kingshighway Blvd.
                  St. Louis, MO 63110
          as part of the 1993, 1994 DICOM Central Test Node project for, and
          under contract with, the Radiological Society of North America.

          THIS SOFTWARE IS MADE AVAILABLE, AS IS, AND NEITHER RSNA NOR
          WASHINGTON UNIVERSITY MAKE ANY WARRANTY ABOUT THE SOFTWARE, ITS
          PERFORMANCE, ITS MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR
          USE, FREEDOM FROM ANY COMPUTER DISEASES OR ITS CONFORMITY TO ANY
          SPECIFICATION. THE ENTIRE RISK AS TO QUALITY AND PERFORMANCE OF
          THE SOFTWARE IS WITH THE USER.

          Copyright of the software and supporting documentation is
          jointly owned by RSNA and Washington University, and free access
          is hereby granted as a license to use this software, copy this
          software and prepare derivative works based upon this software.
          However, any distribution of this software source code or
          supporting documentation or derivative works (source code and
          supporting documentation) must include the three paragraphs of
          the copyright notice.
*/
/* Copyright marker.  Copyright will be inserted above.  Do not remove */
/*
** @$=@$=@$=
*/
/*
**    DICOM 93
**       Electronic Radiology Laboratory
**     Mallinckrodt Institute of Radiology
**  Washington University School of Medicine
**
** Module Name(s):
** Author, Date: Stephen M. Moore, 15-Apr-93
** Intent:  This header defines public typedefs for the DICOM
**   software produced at the Mallinckrodt Institute of
**   Radiology.  It also defines unique identifiers
**   for standard classes and objects defined by the
**   standard.
** Last Update:  $Author: nicks $, $Date: 2006/12/29 02:09:01 $
** Source File:  $RCSfile: dicom.h,v $
** Revision:  $Revision: 1.5 $
** Status:  $State: Exp $
*/

#ifndef DICOM_IS_IN
#define DICOM_IS_IN 1

/* RKT - commented this out, don't know why */
/* #ifdef _MSC_VER */
#include "dicom_platform.h"
/* #endif */

#ifdef  __cplusplus
extern "C"
{
#endif

#ifndef _SITE_MACROS
  typedef unsigned long CONDITION;
  typedef unsigned short U_SHORT; /* normal unsigned short */
  typedef unsigned long U_LONG; /* normal unsigned long */
  typedef unsigned long MASK_32; /* For bit masks  */
  typedef unsigned long CTNBOOLEAN; /* for boolean ops  */

#if !defined(SHORTSIZE) || SHORTSIZE != 16
  /* The writers of this code assume that shorts are 16 bits long.
  ** If that is not the case, this system will not operate properly.
  ** This code will trip the compiler.  This code also tripes the
  ** compiler if you have not defined the macro SHORTSIZE.  You
  ** also want to define INTSIZE and LONGSIZE.
  */

  /* RKT - commented this out, our shorts are ints. but it seems to
     work anyway. */
  /*     short c; */
  /*     char c;  */  /* See note above */
#endif

  typedef unsigned short U16; /* unsigned, 16 bit */
  typedef short S16;  /* signed, 16 bit */

#if LONGSIZE==32
  /* note that under both 64 bit and 32 bit OS, unsigned int is 32 bit
     and int is 32 bit.  Thus we can avoid setting the above */
  typedef unsigned int U32;
  typedef int S32;
  // the following is unnecessary
  // #if LONGSIZE == 64 && INTSIZE == 32 /* Such as an Alpha */
  //  typedef unsigned int U32;
  //  typedef int S32;

  // #elif LONGSIZE == 32  /* Most 32 bit workstations */
  //  typedef unsigned long U32;
  //   typedef long S32;
#else    /* Something we do not support */

  /* The writers of this code assume that we can find a 32 bit integer
  ** defined for this system as an int or a long.  If that assumption
  ** is not true, this code will not operate properly.
  ** This code will trip the compiler.
  */

  /* RKT - commented this out, our shorts are ints. but it seems to
  work anyway. */
  /*     short c; */
  /*     char c;  */  /* See note above */

#endif

#endif

#define FORM_COND(facility, severity, value) \
 (CONDITION)((((unsigned long)value)<<16) | \
 (((unsigned long)facility) << 4) | ((unsigned long)severity))

#define SEV_SUCC 1
#define SEV_INFORM 3
#define SEV_WARN 5
#define SEV_ERROR 2
#define SEV_FATAL 4

#define CTN_SUCCESS(A) (((A)&0xf) == SEV_SUCC)
#define CTN_INFORM(A) (((A)&0xf) == SEV_INFORM)
#define CTN_WARNING(A) (((A)&0xf) == SEV_WARN)
#define CTN_ERROR(A) (((A)&0xf) == SEV_ERROR)
#define CTN_FATAL(A) (((A)&0xf) == SEV_FATAL)

#if 0
  /* We turn these on to force compiler errors to find dependencies
  ** on these older macros.  These are retired as of 2.8.6.
  */
#define SUCCESS(A) (zzzz)
#define INFORM(A) (zzzz)
#define WARNING(A) (zzzz)
#define ERROR(A) (zzzz)
#define FATAL(A) (zzzz)
#endif

#define FACILITY(A) ((unsigned long)(A)>>4) & 0xfff

#ifndef _FACILITY_CODES
#define FAC_DUL  1
#define FAC_IDBMB 2 /* IDB Multibyte */
#define FAC_IDX  3
#define FAC_LST  4
#define FAC_DIAG 5
#define FAC_COND 6
#define FAC_GQ  7
#define FAC_SRV  8
#define FAC_DCM  9
#define FAC_MSG  10
#define FAC_HUNK 11
#define FAC_DB  12
#define FAC_CFG  13
#define FAC_IAP  14
#define FAC_HIS  15
#define FAC_HAP  16
#define FAC_IE  17
#define FAC_UID  18
#define FAC_SQ  19
#define FAC_ICON 20
#define FAC_PRN  21
#define FAC_TBL  22 /* Table functions (relational database) */
#define FAC_DMAN 23 /* DICOM Management of application
  * connections */
#define FAC_UTL  24 /* Utility functions */
#define FAC_IDB  25 /* Image database */
#define FAC_MUT  26 /* Motif utilities */
#define FAC_IMAN 27 /* Image management */
#define FAC_ICPY 30 /* Image copy (structures for queueing) */
#define FAC_FIS  31 /* Fake information system */
#define FAC_SNP  32 /* TCP/IP snoop facility */
#define FAC_LUT  34 /* LUT facility */
#define FAC_IODV 35 /* IOD Verification */
#define FAC_THR  36 /* CTN Threading routines */
#define FAC_DDR  37 /* DICOM Directory Services */
#define FAC_ATH  38 /* Application thread usage */
#define FAC_IRS  39 /* Image recycle system */
#define FAC_TBLMB 40 /* Table functions (relational database) */
#define FAC_CHR  41 /* Character set encoding utilities */

#define FAC_MAXIMUM 50 /* Maximum number of facilities.  This can increase */

#define FAC_APP 0x0fff  /* for stand-alone programs  */
#endif

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif


#ifndef MAX
#define MAX(x, y) (((x) < (y)) ? (y) : (x))
#endif
#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif
#define IS_EVEN(i) (~(i) & 0x01)
#define DIM_OF(a) (sizeof(a) / sizeof(a[0]))
#define IN_RANGE(n, lo, hi) ((lo) <= n && (n) <= (hi))
#define STRUCT_OFFSET(s, f)  (off_t)(((s *)(0))->f)

#ifdef NO_STRERROR
  static char *
  strerror(int e)
  {
    static char string[256];

    sprintf(string, "Error number: %d", e);
    return string;
  }
#endif

#define DICOM_AS_LENGTH 4
#define DICOM_CS_LENGTH 16
#define DICOM_DS_LENGTH 16
#define DICOM_IS_LENGTH 12
#define DICOM_PN_LENGTH 64
#define DICOM_DA_LENGTH 8
#define DICOM_LO_LENGTH 64
#define DICOM_TM_LENGTH 16
#define DICOM_UI_LENGTH 64
#define DICOM_SH_LENGTH 16
#define DICOM_AE_LENGTH 16
#define DICOM_ST_LENGTH 1024
#define DICOM_LT_LENGTH 10240
#define DICOM_DT_LENGTH 26

#define VERSION_JUN1993 199306
#define VERSION_JUL1993 199307
#define VERSION_AUG1993 199308
#define VERSION_SEP1993 199309
#define VERSION_OCT1993 199310
#define VERSION_NOV1993 199311
#define VERSION_DEC1993 199312
#define VERSION_JAN1994 199401
#define VERSION_FEB1994 199402
#define VERSION_MAR1994 199403
#define VERSION_APR1994 199404
#define VERSION_MAY1994 199405
#define VERSION_JUN1994 199406
#define VERSION_JUL1994 199407
#define VERSION_AUG1994 199408
#define VERSION_SEP1994 199409
#define VERSION_OCT1994 199410
#define VERSION_NOV1994 199411
#define VERSION_DEC1994 199412
#define VERSION_JAN1995 199501
#define VERSION_FEB1995 199502
#define VERSION_MAR1995 199503
#define VERSION_APR1995 199504
#define VERSION_MAY1995 199505
#define VERSION_JUN1995 199506

#ifndef STANDARD_VERSION
#define STANDARD_VERSION VERSION_JUN1995
#endif

#ifdef MALLOC_TEST

#define malloc(a) (void *)COND_Malloc((a), __FILE__, __LINE__)
#define free(a) (void)COND_Free((a), __FILE__, __LINE__)

#endif

#define CTN_MALLOC(a) malloc((a))
#define CTN_FREE(a)   free((a))

#ifdef  __cplusplus
}
#endif

#endif
