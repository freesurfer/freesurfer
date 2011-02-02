/**
 * @file  xUndoList.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/02 19:25:19 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifndef xUndoList_H
#define xUndoList_H

#include "xList.h"

typedef enum
{

  xUndL_tErr_NoErr = 0,
  xUndL_tErr_AllocationFailed,
  xUndL_tErr_InternalAllocationFailed,
  xUndL_tErr_ItemDeletionFailed,
  xUndL_tErr_ListDeletionFailed,
  xUndL_tErr_InvalidFunctionPtr,
  xUndL_tErr_InsertionFailed,
  xUndL_tErr_RetrievalFailed,
  xUndL_tErr_AllocationOfSwapListFailed,
  xUndL_tErr_SwapListInsertionFailed,
  xUndL_tErr_InvalidListPtr,
  xUndL_tErr_ClearFailed,
  xUndL_tErr_ListCountFailed,
  xUndL_tErr_ListResetFailed,
  xUndL_tErr_InvalidPrintFunction,
  xUndL_tErr_InvalidErrorCode,
  xUndL_knNumErrorCodes

} xUndL_tErr;

#define xUndL_kSignature 0xf0f9f8f7

typedef void* xUndL_tEntryPtr;

typedef void(*xUndL_tSwapFuncPtr)       ( xUndL_tEntryPtr, xUndL_tEntryPtr* );
typedef void(*xUndL_tDeleteEntryFuncPtr)( xUndL_tEntryPtr* );
typedef void(*xUndL_tPrintEntryFuncPtr) ( xUndL_tEntryPtr );

typedef struct
{

  long                       mSignature;

  xListRef                   mpList;
  xUndL_tSwapFuncPtr         mpSwapFunction;
  xUndL_tDeleteEntryFuncPtr  mpDeleteFunction;
  xUndL_tPrintEntryFuncPtr   mpPrintFunction;

}
xUndoList, *xUndoListRef;

xUndL_tErr xUndL_New ( xUndoListRef*             oppList,
                       xUndL_tSwapFuncPtr        ipSwapFunction,
                       xUndL_tDeleteEntryFuncPtr ipDeleteFunction );

xUndL_tErr xUndL_Delete ( xUndoListRef* ioppList );

xUndL_tErr xUndL_Clear ( xUndoListRef ipList );

xUndL_tErr xUndL_AddEntry ( xUndoListRef    ipList,
                            xUndL_tEntryPtr ipEntry );

xUndL_tErr xUndL_Restore ( xUndoListRef ipList );

xUndL_tErr xUndL_SetPrintFunction ( xUndoListRef             ipList,
                                    xUndL_tPrintEntryFuncPtr ipPrintFunction );

xUndL_tErr xUndL_Print ( xUndoListRef ipList );

xUndL_tErr xUndL_Verify ( xUndoListRef ipList );

char * xUndL_GetErrorString ( xUndL_tErr ieCode );

/* private */

xUndL_tErr xUndL_NewList_ ( xListRef* ippList );

xUndL_tErr xUndL_DeleteList_ ( xUndoListRef ipUndoList,
                               xListRef*    oppList );





#endif
