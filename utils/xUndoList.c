/**
 * @file  xUndoList.c
 * @brief general-purpose utils
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 1.12 $
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


#include <stdlib.h>
#include "xUndoList.h"
#include "xDebug.h"

static char *xUndL_ksaError [xUndL_knNumErrorCodes]  =
{

  "No error.",
  "Allocation of undo list structure failed.",
  "Allocation of internal list structure failed.",
  "Couldn't delete items from list.",
  "Couldn't delete list.",
  "Invalid swap or delete function ptr.",
  "Couldn't add item to list.",
  "Couldn't retrieve item from list.",
  "Couldn't allocate swap list (non-fatal).",
  "Couldn't add item to swap list (non-fatal).",
  "Invalid undo list ptr.",
  "Failed to clear the list.",
  "Failed to count items in the list.",
  "Failed to reset the list position.",
  "Invalid print function.",
  "Invalid error code."
};

xUndL_tErr xUndL_New ( xUndoListRef*               oppList,
                       xUndL_tSwapFuncPtr          ipSwapFunction,
                       xUndL_tDeleteEntryFuncPtr   ipDeleteFunction )
{

  xUndoListRef this     = NULL;
  xUndL_tErr   eResult  = xUndL_tErr_NoErr;

  // assume failure.
  *oppList = NULL;

  // check our functions.
  if ( NULL == ipSwapFunction
       || NULL == ipDeleteFunction )
  {
    eResult = xUndL_tErr_InvalidFunctionPtr;
    goto cleanup;
  }

  // allocate our private structure.
  this = (xUndoListRef) malloc ( sizeof(xUndoList) );
  if ( NULL == this )
  {
    eResult = xUndL_tErr_AllocationFailed;
    goto cleanup;
  }

  // set the sig.
  this->mSignature = xUndL_kSignature;

  // allocate us a list.
  eResult = xUndL_NewList_ ( &(this->mpList) );
  if ( xUndL_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // set the function ptrs.
  this->mpSwapFunction     = ipSwapFunction;
  this->mpDeleteFunction   = ipDeleteFunction;
  this->mpPrintFunction    = NULL;

  // return our list;
  *oppList = this;

cleanup:

  return eResult;
}

xUndL_tErr xUndL_NewList_ ( xListRef* oppList )
{

  xListRef    pList       = NULL;
  xList_tErr  eListResult = xList_tErr_NoErr;
  xUndL_tErr  eResult     = xUndL_tErr_NoErr;

  // assume failure.
  *oppList = NULL;

  // try to allocate our xList.
  eListResult = xList_New ( &pList );
  if ( xList_tErr_NoErr != eListResult )
  {
    eResult = xUndL_tErr_InternalAllocationFailed;
    goto cleanup;
  }

  // set the list.
  *oppList = pList;

cleanup:

  return eResult;
}

xUndL_tErr xUndL_Delete ( xUndoListRef* ioppList )
{

  xUndoListRef this    = NULL;
  xUndL_tErr   eResult = xUndL_tErr_NoErr;

  // get and verify ourself.
  this = *ioppList;
  eResult = xUndL_Verify ( this );
  if ( xUndL_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // delete the internal structure.
  eResult = xUndL_DeleteList_ ( this, &(this->mpList) );
  if ( xUndL_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // mangle our signature
  this->mSignature = 0x1;

  // delete the structure
  free ( this );

  // set outgoing ptr.
  *ioppList = NULL;

cleanup:

  return eResult;
}

xUndL_tErr xUndL_DeleteList_ ( xUndoListRef ipUndoList,
                               xListRef*    ioppList )
{

  xUndL_tErr       eResult     = xUndL_tErr_NoErr;
  xListRef         pList       = NULL;
  xList_tErr       eListResult = xList_tErr_NoErr;
  xUndL_tEntryPtr  pEntry      = NULL;

  // get the list.
  pList = *ioppList;

  // we need to pop every entry off the list...
  pEntry = NULL;
  eListResult = xList_tErr_NoErr;
  while ( xList_tErr_NoErr == eListResult )
  {

    eListResult = xList_PopItem ( pList, (void**)&pEntry );

    // if we got one, call our delete function on it.
    if ( pEntry )
    {
      ipUndoList->mpDeleteFunction ( &pEntry );
    }
  }

  // check for odd results.
  if ( xList_tErr_EndOfList != eListResult )
  {
    eResult = xUndL_tErr_ItemDeletionFailed;
    goto cleanup;
  }

  // now delete the list.
  eListResult = xList_Delete ( &pList );
  if ( xList_tErr_NoErr != eListResult )
  {
    eResult = xUndL_tErr_ListDeletionFailed;
    goto cleanup;
  }

cleanup:

  DebugCode
  if ( xList_tErr_NoErr != eListResult )
  {
    DebugPrint(
      ("xUndL_DeleteInternalList(): Error in xList function %d: %s\n",
       eListResult, xList_GetErrorString ( eListResult ) ) );
  }
  EndDebugCode;

  return eResult;
}

xUndL_tErr xUndL_Clear ( xUndoListRef ipList )
{

  xUndoListRef     this        = NULL;
  xUndL_tErr       eResult     = xUndL_tErr_NoErr;
  xList_tErr       eListResult = xList_tErr_NoErr;
  xUndL_tEntryPtr  pEntry      = NULL;

  // get and verify ourself.
  this = ipList;
  eResult = xUndL_Verify ( this );
  if ( xUndL_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // we need to pop every entry off the list...
  pEntry = NULL;
  eListResult = xList_tErr_NoErr;
  while ( xList_tErr_NoErr == eListResult )
  {

    eListResult = xList_PopItem ( this->mpList, (void**)&pEntry );

    // if we got one, call our delete function on it.
    if ( pEntry )
    {
      this->mpDeleteFunction ( &pEntry );
    }
  }

cleanup:

  return eResult;
}

xUndL_tErr xUndL_AddEntry ( xUndoListRef    ipList,
                            xUndL_tEntryPtr ipEntry )
{

  xUndoListRef this        = NULL;
  xUndL_tErr   eResult     = xUndL_tErr_NoErr;
  xList_tErr   eListResult = xList_tErr_NoErr;

  // get and verify ourself.
  this = ipList;
  eResult = xUndL_Verify ( this );
  if ( xUndL_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  /* push the entry on the list. we want to
  actually add the ptr to the ptr to the
  list, not thing it points to, so pass
  in a ptr to the in ptr. */
  eListResult = xList_PushItem ( this->mpList,
                                 ipEntry );
  if ( xList_tErr_NoErr != eListResult )
  {
    eResult = xUndL_tErr_InsertionFailed;
    goto cleanup;
  }

cleanup:

  DebugCode
  if ( xList_tErr_NoErr != eListResult )
  {
    DebugPrint( ("xUndL_AddEntry(): Error in xList function %d: %s\n",
                 eListResult, xList_GetErrorString ( eListResult ) ) );
  }
  EndDebugCode;

  return eResult;
}

xUndL_tErr xUndL_Restore ( xUndoListRef ipList )
{

  xUndoListRef    this        = NULL;
  xUndL_tErr      eResult     = xUndL_tErr_NoErr;
  xListRef        pSwapList   = NULL;
  xList_tErr      eListResult = xList_tErr_NoErr;
  xUndL_tEntryPtr pItem       = NULL;
  xUndL_tEntryPtr pSwapItem   = NULL;

  // get and verify ourself.
  this = ipList;
  eResult = xUndL_Verify ( this );
  if ( xUndL_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // make a swap list.
  eResult = xUndL_NewList_ ( &pSwapList );
  if ( xUndL_tErr_NoErr != eResult )
  {

    // if error, just mark the swap list empty, not the end of the world.
    eResult = xUndL_tErr_AllocationOfSwapListFailed;
    pSwapList = NULL;
  }

  // pop every item on the list. for each one...
  eListResult = xList_tErr_NoErr;
  while ( xList_tErr_EndOfList != eListResult )
  {

    // pop the item.
    eListResult = xList_PopItem ( this->mpList, (void**)&pItem );

    // if we got one..
    if ( NULL != pItem )
    {

      // call our swap function on it.
      this->mpSwapFunction ( pItem, (void**)&pSwapItem );

      // if we have a swap list..
      if ( NULL != pSwapList )
      {

        // add the swapped item to the swap list.
        eListResult = xList_PushItem ( pSwapList, pSwapItem );
        if ( xList_tErr_NoErr != eListResult )
        {
          eResult = xUndL_tErr_SwapListInsertionFailed;
          goto cleanup;
        }
      }

      // call the deletion function on the item.
      this->mpDeleteFunction ( &pItem );
    }
  }

  // clear the end of list flag
  if ( xList_tErr_EndOfList == eListResult )
  {
    eListResult = xList_tErr_NoErr;
  }
  else
  {
    goto cleanup;
  }

  // if we have a swap list...
  if ( NULL != pSwapList )
  {

    // delete the old list.
    eResult = xUndL_DeleteList_ ( this, &(this->mpList) );
    if ( xUndL_tErr_NoErr != eResult )
    {
      goto cleanup;
    }

    // make our swap list the new one.
    this->mpList = pSwapList;
  }

cleanup:

  DebugCode
  if ( xList_tErr_NoErr != eListResult )
  {
    DebugPrint( ("xUndL_RestoreList(): Error in xList function %d: %s\n",
                 eListResult, xList_GetErrorString ( eListResult ) ) );
  }
  EndDebugCode;

  return eResult;
}


xUndL_tErr xUndL_SetPrintFunction ( xUndoListRef             this,
                                    xUndL_tPrintEntryFuncPtr ipPrintFunction )
{

  xUndL_tErr eResult = xUndL_tErr_NoErr;

  // verify ourself.
  eResult = xUndL_Verify ( this );
  if ( xUndL_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // set the function.
  this->mpPrintFunction = ipPrintFunction;

cleanup:

  return eResult;
}

xUndL_tErr xUndL_Print ( xUndoListRef this )
{

  xUndL_tErr      eResult     = xUndL_tErr_NoErr;
  xList_tErr      eListResult = xList_tErr_NoErr;
  xUndL_tEntryPtr pItem       = NULL;
  int             nCount      = 0;

  // verify ourself.
  eResult = xUndL_Verify ( this );
  if ( xUndL_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // make sure we have a print function.
  if ( NULL == this->mpPrintFunction )
  {
    eResult = xUndL_tErr_InvalidPrintFunction;
    goto cleanup;
  }

  DebugPrint( ("Printing undo list....\n" ) );

  // get the count
  eListResult = xList_GetCount ( this->mpList, &nCount );
  if ( xList_tErr_NoErr != eListResult )
  {
    eResult = xUndL_tErr_ListCountFailed;
    goto cleanup;
  }

  DebugPrint( ("\nNum items = %d\n", nCount ) );

  // start at the beginning.
  eListResult = xList_ResetPosition ( this->mpList );
  if ( xList_tErr_NoErr != eListResult )
  {
    eResult = xUndL_tErr_ListResetFailed;
    goto cleanup;
  }

  // traverse the list.
  pItem = NULL;
  eListResult  = xList_tErr_NoErr;
  while ( eListResult != xList_tErr_EndOfList )
  {

    eListResult = xList_GetNextItemFromPosition ( this->mpList,
                  (void**)&pItem );
    if ( xList_tErr_NoErr != eListResult
         && xList_tErr_EndOfList != eListResult )
    {
      eResult = xUndL_tErr_RetrievalFailed;
      goto cleanup;
    }

    // if we got an item, print it.
    if ( pItem )
    {
      DebugPrint( ("\t" ) );
      this->mpPrintFunction ( pItem );
    }
  }

  DebugPrint( ("\tDone.\n\n" ) );

cleanup:

  DebugCode
  if ( xList_tErr_NoErr != eListResult )
  {
    DebugPrint( ("xUndL_PrintList(): Error in xList function %d: %s\n",
                 eListResult, xList_GetErrorString ( eListResult ) ) );
  }
  EndDebugCode;

  return eResult;
}


xUndL_tErr xUndL_Verify ( xUndoListRef ipList )
{

  xUndL_tErr eResult = xUndL_tErr_NoErr;

  // check for a null list ptr.
  if ( NULL == ipList )
  {
    eResult = xUndL_tErr_InvalidListPtr;
  }

  // check for our sig.
  if ( ipList->mSignature != xUndL_kSignature )
  {
    eResult = xUndL_tErr_InvalidListPtr;
  }

  // check for null function ptrs.
  if ( NULL == ipList->mpSwapFunction
       || NULL == ipList->mpDeleteFunction )
  {
    eResult = xUndL_tErr_InvalidFunctionPtr;
  }

  return eResult;
}

char * xUndL_GetErrorString ( xUndL_tErr ieCode )
{

  xUndL_tErr eCode = ieCode;

  if ( ieCode < xUndL_tErr_NoErr
       || ieCode >= xUndL_knNumErrorCodes )
  {

    eCode = xUndL_tErr_InvalidErrorCode;
  }

  return xUndL_ksaError [ eCode ];

}
