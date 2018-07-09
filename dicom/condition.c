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
** Module Name(s): COND_PushCondition
**   COND_ExtractConditions
**   COND_TopCondition
**   COND_PopCondition
**   COND_DumpConditions
**   COND_CopyText
**   COND_WriteConditions
** Author, Date: Stephen M. Moore, 15-Apr-93
** Intent:  This file contains functions implementing a simple
**   error facility.  It was first written by Stephen Moore
**   (smm@wuerl.wustl.edu) to support PACS development at
**   the Mallinckrodt Institute of Radiology.  The function
**   names have been modified to have a slightly more
**   generic name, but the basic model and concepts are
**   the same.
**
**   The condition package maintains a stack of
**   <condition, message> pairs that callers can push or
**   pop.  When a routine returns an abnormal value, it
**   should push a condition onto the stack so that the
**   caller can examine the value at a later time.  Nested
**   routines may push a number of conditions onto the
**   stack providing more detailed information about why
**   a routine did not return normally.
**
**   The stack is maintained as a simple stack array.  If
**   it overflows, we dump the stack to stdout and reset it.
**
** Last Update:  $Author: nicks $, $Date: 2006/12/29 02:08:57 $
** Source File:  $RCSfile: condition.c,v $
** Revision:  $Revision: 1.8 $
** Status:  $State: Exp $
*/

/*
**
**  INCLUDE FILES
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "dicom.h"
#include "lst.h"
#include "condition.h"
#include "ctnthread.h"

typedef struct {
  CONDITION statusCode;
  char statusText[256];
}
EDB;

#define MAXEDB  100

static int stackPtr = -1;
static EDB EDBStack[MAXEDB];
static void (*ErrorCallback) (CONDITION, const char*) = NULL;
static void dumpstack(FILE * fp);


/*
**++
**  FUNCTIONAL DESCRIPTION:
**
**      COND_PushCondition
** This routine is used to log a condition on the stack.  The user
** passes an error code (currently uninterpreted), a format string
** and optional format arguments.  We use the vsprintf routine to
** interpret the user's format string and arguments, and we place the
** error condition and resultant string on the stack.
**
**  FORMAL PARAMETERS:
**
**      code:
**          The condition code of the error/warning.
**
**      controlString:
**          format string for vsprintf statement
**
**      [varargs]:
**          variable arguments to be used with controlString
**
**  RETURN VALUE:
**
**      code (as passed in by the user)
**
**  SIDE EFFECTS:
**
**      Places a new entry on the stack.  If the stack
** fills up, drop the last condition.
** Calls a user-established callback just before return.
**
*/
CONDITION
COND_PushCondition(CONDITION cond, const char *controlString,...) {
  va_list
  args;
  char
  buffer[1024];

  /*lint -e40 -e50 */
  va_start(args, controlString);
  if (controlString == NULL)
    controlString = "NULL Control string passedto PushCondition";
  (void) vsprintf(buffer, controlString, args);
  va_end(args);
  /*lint +e40 +e50 */

#ifdef CTN_USE_THREADS
  if (THR_ObtainMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_PushCondition unable to obtain mutex\n");
    return cond;
  }
#endif

  stackPtr++;
  EDBStack[stackPtr].statusCode = cond;
  buffer[256] = '\0';

  (void) strcpy(EDBStack[stackPtr].statusText, buffer);
  if (ErrorCallback != NULL)
    ErrorCallback(EDBStack[stackPtr].statusCode,
                  EDBStack[stackPtr].statusText);

  if (stackPtr >= MAXEDB - 2) {
    dumpstack(stderr);
    fprintf(stderr, "CONDITION Stack overflow\n");
    stackPtr = 0;
  }
#ifdef CTN_USE_THREADS
  if (THR_ReleaseMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_PushCondition unable to release mutex, exiting\n");
    exit(1);
  }
#endif

  return cond;

}


/*
**++
**  FUNCTIONAL DESCRIPTION:
**
**  COND_ExtractConditions
** This routine walks through the stack and passes the condition
** codes and text back to the user.  The caller supplies a
** callback routine.  We start at the top of the stack and
** call the user's callback for each message on the stack.  The
** user can terminate the process at any time by returning
** a zero from his callback.
**
**  FORMAL PARAMETERS:
**
**      callback:
**          User routine to call for each message on the stack.
**
**  RETURN VALUE:
**
**      1
**
**  SIDE EFFECTS:
**
**
**  DESIGN:
**
**      None
**--
*/

CONDITION
COND_ExtractConditions(CTNBOOLEAN(*callback) (CONDITION, const char*)) {
  int
  index,
  returnflag;

#ifdef CTN_USE_THREADS
  if (THR_ObtainMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_ExtractConditions unable to obtain mutex, exiting\n");
    exit(1);
  }
#endif

  for (index = stackPtr, returnflag = 1; index >= 0 && returnflag != 0;
       index--) {
    returnflag = callback(EDBStack[index].statusCode,
                          EDBStack[index].statusText);
  }

#ifdef CTN_USE_THREADS
  if (THR_ReleaseMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_ExtractConditions unable to release mutex, exiting\n");
    exit(1);
  }
#endif
  return COND_NORMAL;
}

/*
**++
**  FUNCTIONAL DESCRIPTION:
**
**      COND_TopCondition
** This routine is used to look at the top condition message on
** the stack.  The user passes pointers to areas to place
** the error message.  The function also returns the code
** for the top error message.  If the stack is empty, the
** success code (0) is returned.
**
**  FORMAL PARAMETERS:
**
**      code:
**          Pointer to the user's area to hold the error code
**
**      text
**          Pointer to the user's area to hold the error text
**
**      maxlength
**          Maximum buffer length in the user's text area
**
**  RETURN VALUE:
**
**      top error code on the stack
**
**  SIDE EFFECTS:
**
**
*/

CONDITION
COND_TopCondition(CONDITION * code, char *text, unsigned long maxlength) {
  CONDITION rtnValue;
#ifdef CTN_USE_THREADS
  if (THR_ObtainMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_TopCondition unable to obtain mutex, exiting\n");
    exit(1);
  }
#endif

  if (stackPtr >= 0) {
    *code = EDBStack[stackPtr].statusCode;
    (void) strncpy(text, EDBStack[stackPtr].statusText, maxlength - 1);
    text[maxlength - 1] = '\0';
    rtnValue = EDBStack[stackPtr].statusCode;
  } else {
    *code = COND_NORMAL;
    *text = '\0';
    rtnValue = COND_NORMAL;
  }

#ifdef CTN_USE_THREADS
  if (THR_ReleaseMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_TopCondition unable to release mutex, exiting\n");
    exit(1);
  }
#endif
  return rtnValue;
}

/*
**++
**  FUNCTIONAL DESCRIPTION:
**
**      COND_PopCondition
** This routine pops one or all conditions off the stack.
** The user passes a flag which indicates the operation.
** After the clear, the current top error code is returned.
** If the stack is empty at this point, the success code (0)
** is returned.
**
**  FORMAL PARAMETERS:
**
**      clearstack:
**          Flag which indicates if the entire stack is to be cleared.
**  0     Just pop the top error
**  non zero    Clear the entire stack
**
**  RETURN VALUE:
**
**      The new top error code.  0 if the stack is empty
**
**  SIDE EFFECTS:
**
**
*/

CONDITION
COND_PopCondition(CTNBOOLEAN clearstack) {
  CONDITION
  value;

#ifdef CTN_USE_THREADS
  if (THR_ObtainMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_PopCondition unable to obtain mutex, exiting\n");
    exit(1);
  }
#endif

  if (stackPtr >= 0)
    value = EDBStack[stackPtr].statusCode;
  else
    value = COND_NORMAL;

  if (clearstack) {
    stackPtr = -1;
  } else if (stackPtr <= 0) {
    stackPtr = -1;
  } else {
    stackPtr--;
  }

#ifdef CTN_USE_THREADS
  if (THR_ReleaseMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_PopCondition unable to release mutex, exiting\n");
    exit(1);
  }
#endif

  return value;
}

/*
**++
**  FUNCTIONAL DESCRIPTION:
**
**      COND_EstablishCallback
** Establishes a callback routine to be called whenever a
** new condition is placed on the stack.  There is no stack
** mechanism for these callbacks, so each new callback routine
** completely supersedes the previous one.
**
**  FORMAL PARAMETERS:
**
**      callback:
**          The new callback routine.  If NULL, this will
**     disable callbacks.
**
**  RETURN VALUE:
**
**      0
**
**  SIDE EFFECTS:
**
**
*/

CONDITION
COND_EstablishCallback(void (*callback) (CONDITION, const char*)) {
#ifdef CTN_USE_THREADS
  if (THR_ObtainMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_EstablishCallback unable to obtain mutex, exiting\n");
    exit(1);
  }
#endif

  ErrorCallback = callback;

#ifdef CTN_USE_THREADS
  if (THR_ReleaseMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_EstablishCallback unable to release mutex, exiting\n");
    exit(1);
  }
#endif
  return COND_NORMAL;
}


/* function name
**
** Purpose:
** Describe the purpose of the function
**
** Parameter Dictionary:
** Define the parameters to the function
**
** Return Values:
**
** Algorithm:
** Description of the algorithm (optional) and any other notes.
*/

void
COND_DumpConditions(void) {
#ifdef CTN_USE_THREADS
  if (THR_ObtainMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_DumpConditions unable to obtain mutex\n");
    return;
  }
#endif

  dumpstack(stderr);
  stackPtr = -1;

#ifdef CTN_USE_THREADS
  if (THR_ReleaseMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_DumpConditions unable to release mutex\n");
    return;
  }
#endif
}

static void
dumpstack(FILE * lfp) {
  int
  index;

  for (index = 0; index <= stackPtr; index++)
    fprintf(lfp, "%8x %s\n", (unsigned int) EDBStack[index].statusCode,
            EDBStack[index].statusText);
}

/*
**++
**  FUNCTIONAL DESCRIPTION:
**
**      COND_CopyText
** This function copies as much text as possible from the
** condition stack and places it in the caller's buffer.
**
**  FORMAL PARAMETERS:
**
**      txt
**          Pointer to the user's area to hold the error text
**
**      length
**          Maximum buffer length in the user's text area
**
**  RETURN VALUE:
**
**      none
**
**  SIDE EFFECTS:
**
*/

void
COND_CopyText(char *txt, size_t length) {
  size_t i;
  int j;

  txt[0] = '\0';

#ifdef CTN_USE_THREADS
  if (THR_ObtainMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_CopyText unable to obtain mutex\n");
    return;
  }
#endif

  j = stackPtr;
  while (length > 2 && j >= 0) {
    i = strlen(EDBStack[j].statusText);
    if (i > length)
      i = length - 2;
    strncpy(txt, EDBStack[j].statusText, i);
    txt[i++] = '\n';
    txt[i] = '\0';
    length -= i;
    txt += i;
    j--;
  }

#ifdef CTN_USE_THREADS
  if (THR_ReleaseMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_CopyText unable to release mutex, exiting\n");
    exit(1);
  }
#endif
}

/* COND_WriteConditions
**
** Purpose:
** Write the condition stack to a file, ie stdout or stderr.
**
** Parameter Dictionary:
** File * lfp, the file to which the stack is written.
**
** Return Values:
**
** Algorithm:
** A reiteration of the COND_DumpConditions except this takes an argument.
*/

void
COND_WriteConditions(FILE * lfp) {
#ifdef CTN_USE_THREADS
  if (THR_ObtainMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_WriteConditions unable to obtain mutex\n");
    return;
  }
#endif
  dumpstack(lfp);
  stackPtr = -1;

#ifdef CTN_USE_THREADS
  if (THR_ReleaseMutex(FAC_COND) != THR_NORMAL) {
    fprintf(stderr, "COND_WriteConditions unable to release mutex\n");
    return;
  }
#endif
}

