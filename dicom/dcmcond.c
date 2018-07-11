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
** Module Name(s): DCM_Message
** Author, Date: Stephen M. Moore, 27-Apr-93
** Intent:  Define the ASCIZ messages that go with each DCM
**   error number and provide a function for looking up
**   the error message.
** Last Update:  $Author: nicks $, $Date: 2007/01/11 20:15:14 $
** Source File:  $RCSfile: dcmcond.c,v $
** Revision:  $Revision: 1.7 $
** Status:  $State: Exp $
*/

// static char rcsid[] = "$Revision: 1.7 $ $RCSfile: dcmcond.c,v $";

#include <stdio.h>
#include <sys/types.h>
#include "dicom.h"
#include "lst.h"
#include "dicom_objects.h"
#include "dcmprivate.h"

typedef struct vector {
  CONDITION cond;
  const char *message;
}
VECTOR;

static VECTOR messageVector[] = {
  {
    DCM_NORMAL, "Normal return from DCM routine"
  },
  {DCM_FILEOPENFAILED, "DCM failed to open file: %s in %s"},
  {DCM_FILEACCESSERROR, "DCM failed to access file: %s in %s"},
  {DCM_OBJECTCREATEFAILED, "DCM failed to create object in %s"},
  {DCM_NULLOBJECT, "NULL object passed to routine %s"},
  {DCM_ILLEGALOBJECT, "Illegal object passed to routine %s"},
  {DCM_ELEMENTNOTFOUND, "Requested element (%x %x) not found in %s"},
  {DCM_ILLEGALSTREAMLENGTH,
   "DCM Illegal Stream Length (%ld) (Not enough data to define a full element) in %s"},
  {DCM_ELEMENTCREATEFAILED, "DCM failed to create element in %s (%04x %04x %d)"},
  {DCM_UNRECOGNIZEDGROUP, "DCM unrecognized group: %04x in %s"},
  {DCM_UNRECOGNIZEDELEMENT, "DCM unrecognized element: (%04x %04x) in %s"},
  {DCM_ELEMENTOUTOFORDER, "DCM group/element out of order (%04x %04x) in %s"},
  {DCM_LISTFAILURE, "DCM routine failed on list operation in %s"},
  {DCM_ILLEGALOPTION, "DCM illegal stream option: %s"},
  {DCM_ILLEGALADD, "DCM attempt to add illegal element: %x %x in %s"},
  {DCM_GETINCOMPLETE, "DCM Get Element incomplete in %s"},
  {DCM_ILLEGALCONTEXT, "DCM Illegal context value in %s"},
  {DCM_ILLEGALREPRESENTATION,
   "DCM Caller specified illegal representation for element (%04x %04x) in %s"},
  {DCM_UNEVENELEMENTLENGTH,
   "DCM attempt to add data element (%x %x) with uneven length (%ld)  in %s"},
  {DCM_ELEMENTLENGTHERROR,
   "DCM Data Element (%04x %04x) longer (%ld) than remaining length (%ld) of \
                                   data in stream or file in %s"},
  {DCM_GROUPNOTFOUND, "Requested group (%x) not found in %s"},
  {DCM_FILECREATEFAILED, "DCM failed to create file %s (%s)"},
  {DCM_FILEIOERROR, "DCM io error on file %s (%s)"},
  {DCM_INSERTFAILED,
   "DCM failed to insert new element (%04x %04x) in %s"},
  {DCM_CANNOTGETSEQUENCEVALUE,
   "DCM Cannot retrieve value of element with SQ representation (%08x) in (%s)"},
  {DCM_FILEDELETEFAILED, "DCM Failed to delete file %s in %s"},
  {DCM_MALLOCFAILURE, "DCM Failed to malloc %ld bytes in %s"},
  {DCM_NULLADDRESS, "DCM NULL address passed to routine %s"},
  {DCM_UNEXPECTEDREPRESENTATION,
   "DCM Routine %s expected %s representation for element %04x %04x"},
  {DCM_BADELEMENTINGROUP,
   "DCM Bad element (%04x %04x) found in group %04x in %s"},
  {DCM_CALLBACKABORTED, "DCM Callback aborted by user in %s"},
  {DCM_READSTREAMFAILED, "DCM Failed to read stream in %s"},
  {DCM_UNRECOGNIZEDVRCODE, "DCM Unrecognized VR code (%s) in %s"},
  {DCM_VRMISMATCH, "DCM Incorrect VR (%s) for attribute with tag %08x"},
  {DCM_EXPORTBUFFERTOOSMALL,
   "DCM Caller's export buffer length (%d) is too short in %s"},
  {DCM_BADOFFSET,
   "DCM Offset value (%d) larger than attribute length (%d) in %s"},
  {DCM_BADLENGTH,
   "DCM Combination of offset, length (%d %d) is longer than element length (%d) in %s"},
  {DCM_NOTASEQUENCE,
   "DCM Attempt to perform sequence operation on element (%04x %04x) not a sequence in %s"},
  {DCM_GENERALWARNING, "DCM General warning in %s: %s"},
  {DCM_UNEVENFRAGMENTLENGTH,
   "DCM attempt to add fragment with uneven length (%ld) in %s"},
  {0, NULL}

};


/* DCM_Message
**
** Purpose:
** Find the ASCIZ message that goes with an DCM error number and
** return a pointer to static memory containing that error message.
**
** Parameter Dictionary:
** condition The error condition for which the message is to be
**   returned
**
** Return Values:
** The error message if a valid error condition was reported else NULL.
**
** Algorithm:
** Description of the algorithm (optional) and any other notes.
*/

const char *
DCM_Message(CONDITION condition) {
  int
    index;

  for (index = 0; messageVector[index].message != NULL; index++)
    if (condition == messageVector[index].cond)
      return messageVector[index].message;

  return NULL;
}

void DCM_DumpVector() {
  int index;

  for (index = 0; index < (int) DIM_OF(messageVector); index++) {
    if (messageVector[index].message != NULL)
      printf("%8x %8ld %s\n", (unsigned int) messageVector[index].cond,
             messageVector[index].cond,
             messageVector[index].message);
  }
}
