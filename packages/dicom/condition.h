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
** Author, Date: Stephen M. Moore, 14-Apr-1993
** Intent:  This module defines function prototypes for the
**   CONDITION facility which is used to record status
**   and error messages on a stack.
** Last Update:  $Author: nicks $, $Date: 2006/12/29 02:09:01 $
** Source File:  $RCSfile: condition.h,v $
** Revision:  $Revision: 1.4 $
** Status:  $State: Exp $
*/

#ifndef COND_IS_IN
#define COND_IS_IN 1

#include <stdio.h>

#ifdef  __cplusplus
extern "C"
{
#endif

  CONDITION COND_PushCondition(CONDITION cond, const char *controlString,...);
  CONDITION
  COND_ExtractConditions(CTNBOOLEAN(*callback) ());
  CONDITION
  COND_TopCondition(CONDITION * condition, char *text,
                    unsigned long maxlength);
  CONDITION COND_PopCondition(CTNBOOLEAN clearstack);
  CONDITION COND_EstablishCallback(void (*callback) ());
  void COND_DumpConditions(void);
  void COND_CopyText(char *txt, size_t length);
  void COND_WriteConditions(FILE * lfp);

  /*  Now define the fixed values for conditions returned by this
  **  package.  Note that FAC_COND is used to generate these
  **  conditions.  This should be defined in some global include
  **  file so that we can keep all of the facilities straight.
  */

#define COND_NORMAL /* Successful return */ \
 FORM_COND(FAC_COND, SEV_SUCC, 1)


#ifdef  __cplusplus
}
#endif

#endif
