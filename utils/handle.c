/**
 * @file  handle.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
 *    $Revision: 1.8 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


/*----------------------------------------------------------------------

      File Name:    handle.c

      Description:  control access to pointers through handles.

----------------------------------------------------------------------*/


/*-----------------------------------------------------------------
              INCLUDE FILES
-----------------------------------------------------------------*/

#ifndef Linux
int not_used_000(void) ;

int not_used_000(void)
{
  int i;
  i=0;
  return(i);
}

#else

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "handle.h"
#include "machine.h"


/*-----------------------------------------------------------------
MACROS AND CONSTANTS
-----------------------------------------------------------------*/

#define DEFAULT_MAX_HANDLES   100
#define MAX_HANDLES         (0x7fff)
#define HANDLE_FREE         0
#define HANDLE_ALLOCATED      1

#ifndef ULONG
#define ULONG  unsigned long
#endif

/*----------------------------------------------------------------------
STRUCTURES
----------------------------------------------------------------------*/

typedef struct
{
  void  *ptr ;
  int status ;
}
HandleInfo ;

/*-----------------------------------------------------------------
PROTOTYPES
-----------------------------------------------------------------*/

static void growHandleTable(void) ;

/*-----------------------------------------------------------------
STATIC DATA
-----------------------------------------------------------------*/

static HandleInfo *handleTable ;
static int      maxHandles = 0L ;
static int      nhandles = 0L ;
static int      freeHandles = 0L ;

/*-----------------------------------------------------------------
FUNCTIONS
-----------------------------------------------------------------*/

/*----------------------------------------------------------------------
Parameters:

Description:
Allocate a new handle and assign the given pointer to it.
If the current handle space is used up, grow the handle table.

Returns:
the newly allocated handle.
----------------------------------------------------------------------*/
PTR_HANDLE
HandleAlloc(void *ptr)
{
  HandleInfo  *handleInfo ;
  PTR_HANDLE  handle ;

  if (nhandles >= maxHandles) growHandleTable() ;

  handle = ++nhandles ;
  handleInfo = handleTable + (handle - 1) ;
  handleInfo->ptr = ptr ;
  handleInfo->status = HANDLE_ALLOCATED ;

  return(handle) ;
}
/*----------------------------------------------------------------------
Parameters:

Description:
Free a previously allocated handle.

Returns:
nothing.
----------------------------------------------------------------------*/
void
HandleFree(PTR_HANDLE handle)
{
  HandleInfo  *handleInfo ;

  if (HandleOk(handle) <= 0)
    ESCAPE(ERROR_BADPARM, "HandleFree: bad handle %d", handle) ;

  handleInfo = handleTable + (handle - 1) ;

  freeHandles++ ;
  handleInfo->status = HANDLE_FREE ;
  handleInfo->ptr = NULL ;
}
/*----------------------------------------------------------------------
Parameters:

Description:
turn a handle into the pointer which it reprents.

Returns:
the pointer which the handle represents.
----------------------------------------------------------------------*/
void *
HandleToPtr(PTR_HANDLE handle)
{
  HandleInfo  *handleInfo ;

  if (HandleOk(handle) <= 0)
    ESCAPE(ERROR_BADPARM, "HandleToPtr: bad handle %d", handle) ;

  handleInfo = handleTable + (handle - 1) ;

  return(handleInfo->ptr) ;
}
/*----------------------------------------------------------------------
Parameters:

Description:
determine whether a handle is a valid one.

Returns:
1  if the handle is ok and allocated.
-1  if the handle is ok and not allocated
0  if the handle is out of range.
----------------------------------------------------------------------*/
int
HandleOk(PTR_HANDLE handle)
{
  HandleInfo  *handleInfo ;

  if ((handle <= (PTR_HANDLE)0) || (handle > nhandles)) return(0) ;

  handleInfo = handleTable + (handle - 1) ;
  if (handleInfo->status == HANDLE_FREE) return(-1) ;

  return(1) ;
}
/*----------------------------------------------------------------------
Parameters:

Description:
We have run out of room for more handles.  Free the current
handle table, and allocate one twice as large, moving all
the current information into it.

Returns:
nothing.
----------------------------------------------------------------------*/
static void
growHandleTable(void)
{
  HandleInfo  *newTable ;
  int     newMaxHandles = 0 ;

  if (!maxHandles) newMaxHandles = DEFAULT_MAX_HANDLES ;
  else if (maxHandles < MAX_HANDLES) newMaxHandles = maxHandles << 1 ;

  if (newMaxHandles <= 0)
    ESCAPE(ERROR_NO_MEMORY, "growHandleTable: too many handles") ;

  newTable = (HandleInfo *)calloc((long)newMaxHandles, sizeof(HandleInfo)) ;
  if (maxHandles > 0)   /* not the first time */
  {
    hmemmove((long huge *)newTable, (long huge *)handleTable,
             (ULONG)maxHandles * sizeof(HandleInfo)) ;
    free(handleTable) ;
  }
  maxHandles = newMaxHandles ;
  handleTable = newTable ;
}

#endif


