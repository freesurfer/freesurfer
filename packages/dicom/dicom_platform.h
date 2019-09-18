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
**    DICOM 97
**       Electronic Radiology Laboratory
**     Mallinckrodt Institute of Radiology
**  Washington University School of Medicine
**
** Module Name(s):
** Author, Date: Stephen M. Moore, 6-Aug-1997
** Intent:  This header defines platform specific values for
**   MIR CTN software.  This version happens to be for
**   Windows NT.
** Last Update:  $Author: nicks $, $Date: 2006/12/29 02:09:01 $
** Source File:  $RCSfile: dicom_platform.h,v $
** Revision:  $Revision: 1.6 $
** Status:  $State: Exp $
*/

#ifndef DICOMPLATFORM_IS_IN
#define DICOMPLATFORM_IS_IN 1

#ifdef __cplusplus
extern "C"
{
#endif

#define LONGSIZE 32
#define INTSIZE 32
#define SHORTSIZE 16

#define TBL_SQLSERVER

  /* RKT - commeted this out, we don't want to use threads. */
  /* #define CTN_USE_THREADS */

#ifdef __cplusplus
}
#endif

#endif
