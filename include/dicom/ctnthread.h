/*
          Copyright (C) 1996 RSNA and Washington University

          The software and supporting documentation for the Radiological
          Society of North America (RSNA) 1993 - 1996 Digital Imaging and
          Communications in Medicine (DICOM) Demonstration were developed
          at the
                  Electronic Radiology Laboratory
                  Mallinckrodt Institute of Radiology
                  Washington University School of Medicine
                  510 S. Kingshighway Blvd.
                  St. Louis, MO 63110
          as part of the 1993 - 1996 DICOM Central Test Node project for, and
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
**    DICOM 96
**       Electronic Radiology Laboratory
**     Mallinckrodt Institute of Radiology
**  Washington University School of Medicine
**
** Module Name(s):
** Author, Date: Stephen M. Moore, 30-Jun-96
** Intent:  This is the include file for the CTN threads.
**   facility.  This facility provides functions for
**   simple thread operations needed to make some
**   parts of the code thread safe.
** Last Update:  $Author: nicks $, $Date: 2006/12/29 02:09:01 $
** Source File:  $RCSfile: ctnthread.h,v $
** Revision:  $Revision: 1.4 $
** Status:  $State: Exp $
*/

#ifndef CTN_THREADS_IS_IN
#define CTN_THREADS_IS_IN 1

#ifdef  __cplusplus
extern "C"
{
#endif

  /* Define the function prototypes for this set of routines.
  ** The first set defines initialization routines for using these
  ** services as a user or provider.
  */
  CONDITION THR_Init(void);
  CONDITION  THR_Shutdown(void);
  CONDITION  THR_ObtainMutex(int fac);
  CONDITION  THR_ReleaseMutex(int fac);

#define THR_ObtainMutexA(a) EEE;
#define THR_ReleaseMutexA(a) FFF;

#define THR_NORMAL   FORM_COND(FAC_THR, SEV_SUCC, 1)
#define THR_GENERICFAILURE  FORM_COND(FAC_THR, SEV_ERROR, 2)
#define THR_NOTINITIALIZED  FORM_COND(FAC_THR, SEV_ERROR, 3)

#ifdef  __cplusplus
}
#endif

#endif
