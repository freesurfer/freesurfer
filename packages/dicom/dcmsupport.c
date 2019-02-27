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
**    DICOM 93
**       Electronic Radiology Laboratory
**     Mallinckrodt Institute of Radiology
**  Washington University School of Medicine
**
** Module Name(s):
**   DCM_ListToString
**   DCM_IsString
** Author, Date: Stephen M. Moore, 13-Jun-93
** Intent:  This file contains more DCM routines which are used
**   as support for the DCM facility and for applications.
**   These routines help parse strings and other data
**   values that are encoded in DICOM objects.
** Last Update:  $Author: nicks $, $Date: 2007/01/11 20:15:14 $
** Source File:  $RCSfile: dcmsupport.c,v $
** Revision:  $Revision: 1.7 $
** Status:  $State: Exp $
*/

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#ifndef MACOS
#include <stdlib.h>
#endif
#ifdef MALLOC_DEBUG
#include "malloc.h"
#endif

#include "dicom.h"
#include "condition.h"
#include "lst.h"
#include "dicom_objects.h"
#include "dcmprivate.h"

/* DCM_ListToString
**
** Purpose:
** Convert the list of strings into a single string separated by '\'
**
** Parameter Dictionary:
** list  Handle to the list of strings
** offset  The actual string starts at "offset" offset in
**   each individual structure chained in the list
** string  The single large string returned to the caller
**
** Return Values:
** DCM_NORMAL
** DCM_LISTFAILURE
** DCM_MALLOCFAILURE
**
** Notes:
**
** Algorithm:
** Description of the algorithm (optional) and any other notes.
*/
typedef struct {
  void *reserved[2];
  char *s;
}
GENERIC;

CONDITION
DCM_ListToString(LST_HEAD * list, long offset, char **string) {
  GENERIC
  * g;
  char
  *c,
  *p;
  long
  length;

  *string = NULL;
  if (list == NULL)
    return DCM_NORMAL;

  g = (GENERIC*)LST_Head(&list);
  if (g == NULL)
    return DCM_NORMAL;

  (void) LST_Position(&list, g);

  length = 0;
  while (g != NULL) {
    c = ((char *) g) + offset;
    length += strlen(c) + 1;
    g = (GENERIC*)LST_Next(&list);
  }

  p = (char*)CTN_MALLOC(length);
  if (p == NULL)
    return COND_PushCondition(DCM_MALLOCFAILURE,
                              DCM_Message(DCM_MALLOCFAILURE), 
                              length, "DCM_ListToString");

  *string = p;
  g = (GENERIC*)LST_Head(&list);
  if (g == NULL)
    return COND_PushCondition(DCM_LISTFAILURE, DCM_Message(DCM_LISTFAILURE),
                              "DCM_ListToString");
  (void) LST_Position(&list, g);

  length = 0;
  while (g != NULL) {
    c = ((char *) g) + offset;
    length = strlen(c);
    (void) memcpy(p, c, length);
    p += length;
    *p++ = '\\';
    g = (GENERIC*)LST_Next(&list);
  }
  *--p = '\0';
  return DCM_NORMAL;
}


/* DCM_IsString
**
** Purpose:
** Verify if the DICOM value representation is that of a string
**
** Parameter Dictionary:
** representation  One of the many DICOM value representations
**
** Return Values:
** TRUE
** FALSE
**
** Notes:
**
** Algorithm:
** Description of the algorithm (optional) and any other notes.
*/

CTNBOOLEAN
DCM_IsString(DCM_VALUEREPRESENTATION representation) {
  CTNBOOLEAN
  flag = FALSE;

  switch (representation) {
  case DCM_AE:  /* Application Entity */
  case DCM_AS:  /* Age string */
    flag = TRUE;
    break;
  case DCM_AT:  /* Attribute tag */
    break;
  case DCM_CS:  /* Control string */
  case DCM_DA:  /* Date */
    flag = TRUE;
    break;
  case DCM_DD:  /* Data set */
    break;
  case DCM_DS:  /* Decimal string */
  case DCM_DT:  /* Old date/time */
    flag = TRUE;
    break;
  case DCM_FD:  /* Floating double */
  case DCM_FL:  /* Float */
    break;
  case DCM_IS:  /* Integer string */
  case DCM_LO:  /* Long string */
  case DCM_LT:  /* Long text */
    flag = TRUE;
    break;
  case DCM_OB:  /* Other binary value (byte) */
  case DCM_OT:  /* Other binary value */
  case DCM_OW:  /* Other binary value (word) */
    break;
  case DCM_SH:  /* Short string */
    flag = TRUE;
    break;
  case DCM_SL:  /* Signed long */
  case DCM_SQ:  /* Sequence of items */
  case DCM_SS:  /* Signed short */
    break;
  case DCM_ST:  /* Short text */
  case DCM_TM:  /* Time */
    flag = TRUE;
    break;
  case DCM_UL:  /* Unsigned long */
  case DCM_US:  /* Unsigned short */
    /*case DCM_UNKNOWN:*/ /* Unknown/unspecified */
  case DCM_RET:  /* Retired */
  case DCM_CTX:  /* Context sensitive */
    break;
  case DCM_PN:  /* Person Name */
  case DCM_UI:  /* Unique identifier (UID) */
  case DCM_UT:  /* Unlimited Text */
    flag = TRUE;
    break;
  default: // DCM_UN and DCM_DLM are not handled  added by tosa
    break;
  };
  return flag;
}
