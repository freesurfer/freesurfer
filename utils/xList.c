/**
 * @file  xList.c
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 1.9 $
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
#include "xTypes.h"
#include "xList.h"
#include "xDebug.h"

static char *xList_ksaError [xList_knNumErrorCodes]  =
{

  "No error",
  "List allocation failed.",
  "Internal list allocation (RTS) failed.",
  "Internal list deletion (RTS) failed.",
  "List deleteion failed.",
  "Invalid list pointer (was NULL).",
  "List is full.",
  "Item not in list.",
  "Error finding item.",
  "Error getting list size.",
  "End of list.",
  "List is empty.",
  "Invalid error code."
};


xList_tErr xList_New ( xListRef* oppList )
{

  xListRef   this       = NULL;
  xList_tErr eResult    = xList_tErr_NoErr;

  // assume failure.
  *oppList = NULL;

  // allocate our list structure and check it.
  this = (xListRef) malloc ( sizeof(xList) );
  if ( NULL == this )
  {
    eResult = xList_tErr_AllocationFailed;
    goto cleanup;
  }

  // set the singnature.
  this->mSignature = xList_kSignature;

  // null all ptrs
  this->mpHead = NULL;
  this->mpTail = NULL;
  this->mpNext = NULL;

  /* no comparator by default */
  this->mComparator = NULL;

  // all good, set the outgoing ptr.
  *oppList = this;

cleanup:

  return eResult;
}

xList_tErr xList_Delete ( xListRef* ioppList )
{

  xListRef   this       = NULL;
  xList_tErr eResult    = xList_tErr_NoErr;

  // grab the list.
  this = *ioppList;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // clear the list first.
  xList_Clear ( this );

  // mangle our signature
  this->mSignature = 0x1;

  // delete our list structure.
  free ( this );

  // set list ptr to nil.
  *ioppList = NULL;

cleanup:

  return eResult;
}

xList_tErr xList_InsertItem ( xListRef this,
                              void* ipItemToInsert )
{

  xList_tErr   eResult    = xList_tErr_NoErr;
  xListNodeRef pNewNode   = NULL;
  tBoolean     bIsInList  = FALSE;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  /* make sure the item isn't already in here. */
  xList_IsInList ( this, ipItemToInsert, &bIsInList );
  if ( bIsInList )
  {
    goto cleanup;
  }

  // make a new node.
  pNewNode = (xListNodeRef) malloc ( sizeof(xListNode) );
  if ( NULL == pNewNode )
  {
    eResult = xList_tErr_InternalAllocationFailed;
    goto cleanup;
  }

  // set its data.
  pNewNode->mSignature = xList_kSignature;
  pNewNode->mpData     = ipItemToInsert;
  pNewNode->mpNext     = NULL;

  // insert it at the tail.
  if ( NULL == this->mpHead )
  {

    this->mpHead = pNewNode;
    this->mpNext = pNewNode;
    this->mpTail = pNewNode;

  }
  else
  {

    this->mpTail->mpNext = pNewNode;
    this->mpTail         = pNewNode;
  }

cleanup:

  return eResult;
}

xList_tErr xList_RemoveItem ( xListRef this,
                              void**   iopItemToRemove )
{

  xList_tErr   eResult   = xList_tErr_NoErr;
  xListNodeRef pCurNode  = NULL;
  xListNodeRef pBackNode = NULL;
  tBoolean         bFound    = FALSE;
  void*        pItemToRemove = NULL;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // get the item to remove.
  pItemToRemove = *iopItemToRemove;

  // scan through the list, keeping a back node.
  pCurNode = this->mpHead;
  while ( NULL != pCurNode
          && !bFound )
  {

    // compare the nodes. if we found it exit the loop.
    if ( xList_CompareItems_( this, pCurNode->mpData, pItemToRemove )
         == xList_tCompare_Match )
    {

      bFound = TRUE;

    }
    else
    {

      // next node.
      pBackNode = pCurNode;
      pCurNode  = pCurNode->mpNext;
    }
  }

  // error if not found.
  if ( !bFound )
  {
    eResult = xList_tErr_ItemNotInList;
    goto cleanup;
  }

  // remove this node.
  if ( NULL == pBackNode )
  {

    // node to be deleted is head.
    this->mpHead = pCurNode->mpNext;

  }
  else
  {

    // wrap list around it.
    pBackNode->mpNext = pCurNode->mpNext;

    // reattach tail if necessary.
    if ( this->mpTail == pCurNode )
    {
      this->mpTail = pBackNode;
    }
  }

  // return the item.
  *iopItemToRemove = pCurNode->mpData;

  // delete the node.
  free ( pCurNode );

cleanup:

  return eResult;
}


xList_tErr xList_IsInList ( xListRef this,
                            void* ipItemToFind, tBoolean* obpIsInList )
{

  xList_tErr    eResult    = xList_tErr_NoErr;
  tBoolean          bFound     = FALSE;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // try and find the item.
  bFound = ( NULL != xList_FindItem_ ( this, ipItemToFind ) );

  // set found var.
  *obpIsInList = bFound;

cleanup:

  return eResult;
}

xListNodeRef xList_FindItem_ ( xListRef this, void *ipItem )
{

  xListNodeRef  pCurNode   = NULL;

  // if no comparator, return null.
  if ( NULL == this->mComparator )
  {
    return NULL;
  }

  // scan through the list.
  pCurNode = this->mpHead;
  while ( NULL != pCurNode )
  {

    // compare the nodes. if we found it, return it.
    if ( xList_CompareItems_( this, pCurNode->mpData, ipItem )
         == xList_tCompare_Match )
    {

      return pCurNode;

    }
    else
    {

      // next node.
      pCurNode  = pCurNode->mpNext;
    }
  }

  // not found.
  return NULL;
}

xList_tErr xList_Clear ( xListRef this )
{

  xList_tErr    eResult    = xList_tErr_NoErr;
  xListNodeRef  pCurNode   = NULL;
  xListNodeRef  pDelNode   = NULL;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // go through the list......
  pCurNode = this->mpHead;
  while ( NULL != pCurNode )
  {

    // mark this node.
    pDelNode = pCurNode;

    // next node.
    pCurNode = pCurNode->mpNext;

    // check the node.
    if ( xList_kSignature != pDelNode->mSignature )
    {
      eResult = xList_tErr_InvalidListRef;
      goto cleanup;
    }

    // mangle sig.
    pDelNode->mSignature = 0x1;

    // delete it.
    free ( pDelNode );
  }

  this->mpHead = NULL;
  this->mpTail = NULL;
  this->mpNext = NULL;

cleanup:

  return eResult;
}

xList_tErr xList_GetCount ( xListRef this,
                            int*     opnCount )
{

  xList_tErr    eResult    = xList_tErr_NoErr;
  xListNodeRef  pCurNode   = NULL;
  int           nCount     = 0;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // go through the list..
  nCount = 0;
  pCurNode = this->mpHead;
  while ( NULL != pCurNode )
  {

    // inc count
    nCount ++;

    // next node.
    pCurNode = pCurNode->mpNext;
  }

  // return count.
  *opnCount = nCount;

cleanup:

  return eResult;
}

xList_tErr xList_GetFirstItem ( xListRef this,
                                void**   oppFirstItem )
{

  xList_tErr eResult    = xList_tErr_NoErr;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // assume failure.
  *oppFirstItem = NULL;

  // if there's a head...
  if ( this->mpHead )
  {

    // return the data ptr.
    *oppFirstItem = this->mpHead->mpData;

  }
  else
  {

    eResult = xList_tErr_ListEmpty;
  }

cleanup:

  return eResult;
}

xList_tErr xList_GetNextItem ( xListRef this,
                               void*    ipCurrentItem,
                               void**   oppNextItem )
{

  xList_tErr    eResult    = xList_tErr_NoErr;
  xListNodeRef  pCurNode   = NULL;
  void*         pNextItem  = NULL;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // assume failure.
  *oppNextItem = FALSE;

  // try to find the item.
  pCurNode = xList_FindItem_ ( this, ipCurrentItem );
  if ( NULL != pCurNode )
  {

    // if there's a next one...
    if ( NULL != pCurNode->mpNext )
    {

      // get its data.
      pNextItem = ((xListNodeRef)pCurNode->mpNext)->mpData;

    }
    else
    {

      // end of list.
      eResult = xList_tErr_EndOfList;
    }

  }
  else
  {

    // couldn't find item.
    eResult = xList_tErr_ItemNotInList;
  }

  // return the result.
  *oppNextItem = pNextItem;

cleanup:

  return eResult;
}

xList_tErr xList_ResetPosition ( xListRef this )
{

  xList_tErr        eResult    = xList_tErr_NoErr;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // set current ptr to head.
  this->mpNext = this->mpHead;

cleanup:

  return eResult;
}

xList_tErr xList_GetNextItemFromPosition ( xListRef this,
    void**   oppNextItem )
{

  xList_tErr    eResult    = xList_tErr_NoErr;
  void*         pNextItem  = NULL;

  // assume failure.
  *oppNextItem = FALSE;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // make sure we have an item.
  if ( NULL == this->mpNext )
  {
    eResult = xList_tErr_EndOfList;
    goto cleanup;
  }

  // verify this item.
  if ( xList_kSignature != this->mpNext->mSignature )
  {
    eResult = xList_tErr_InvalidListRef;
    goto cleanup;
  }

  // get its item.
  pNextItem = this->mpNext->mpData;

  // advance the ptr.
  this->mpNext = this->mpNext->mpNext;

  // return the ptr.
  *oppNextItem = pNextItem;

cleanup:

  return eResult;
}

xList_tErr xList_NextFromPos ( xListRef this,
                               void**   oppNextItem )
{

  return xList_GetNextItemFromPosition( this, oppNextItem );
}

xList_tErr xList_PushItem ( xListRef this,
                            void* ipItemToInsert )
{
  xList_tErr   eResult    = xList_tErr_NoErr;
  xListNodeRef pNewNode   = NULL;
  tBoolean     bIsInList  = FALSE;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  /* make sure the item isn't already in here. */
  xList_IsInList ( this, ipItemToInsert, &bIsInList );
  if ( bIsInList )
  {
    goto cleanup;
  }

  // make a new node.
  pNewNode = (xListNodeRef) malloc ( sizeof(xListNode) );
  if ( NULL == pNewNode )
  {
    eResult = xList_tErr_InternalAllocationFailed;
    goto cleanup;
  }

  // set its data.
  pNewNode->mSignature = xList_kSignature;
  pNewNode->mpData     = ipItemToInsert;
  pNewNode->mpNext     = NULL;

  // insert it at the head.
  if ( NULL == this->mpHead )
  {

    this->mpHead = pNewNode;
    this->mpNext = pNewNode;
    this->mpTail = pNewNode;

  }
  else
  {

    pNewNode->mpNext = this->mpHead;
    this->mpHead     = pNewNode;
  }

cleanup:

  return eResult;

}

xList_tErr xList_PopItem ( xListRef this, void** oppItem )
{

  xList_tErr    eResult  = xList_tErr_NoErr;
  xListNodeRef  pDelNode = NULL;
  void*         pItem    = NULL;

  // assume failure.
  *oppItem = FALSE;

  // verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  // make sure we have a head.
  if ( NULL == this->mpHead )
  {
    eResult = xList_tErr_EndOfList;
    goto cleanup;
  }

  // verify this item.
  if ( xList_kSignature != this->mpHead->mSignature )
  {
    eResult = xList_tErr_InvalidListRef;
    goto cleanup;
  }

  // get the head's data.
  pItem = this->mpHead->mpData;

  // remove the head.
  pDelNode = this->mpHead;
  this->mpHead = this->mpHead->mpNext;

  // delete the head.
  free ( pDelNode );

  // return the item.
  *oppItem = pItem;

cleanup:

  return eResult;
}

xList_tErr xList_SetComparator ( xListRef             this,
                                 xList_tCompare(*iComparator)(void*,void*) )
{
  xList_tErr    eResult  = xList_tErr_NoErr;

// verify the list.
  eResult = xList_Verify ( this );
  if ( xList_tErr_NoErr != eResult )
  {
    goto cleanup;
  }

  /* set the comparator */
  this->mComparator = iComparator;

  goto cleanup;

cleanup:

  return eResult;
}

xList_tCompare xList_CompareItems_ ( xListRef this,
                                     void*    pItemA,
                                     void*    pItemB )
{

  xList_tCompare eResult = xList_tCompare_Match;

  if ( this->mComparator )
  {
    eResult = this->mComparator( pItemA, pItemB );
  }
  else
  {
    if ( pItemA == pItemB)
    {
      eResult = xList_tCompare_Match;
    }
    else
    {
      eResult = xList_tCompare_GreaterThan;
    }
  }

  return eResult;
}

xList_tErr xList_Verify ( xListRef this )
{

  xList_tErr eResult    = xList_tErr_NoErr;

  // check pointer.
  if ( NULL == this )
  {
    eResult = xList_tErr_InvalidListRef;
    goto cleanup;
  }

  // check signature.
  if ( xList_kSignature != this->mSignature )
  {
    eResult = xList_tErr_InvalidListRef;
    goto cleanup;
  }

cleanup:

  return eResult;
}

char* xList_GetErrorString ( xList_tErr ieCode )
{

  xList_tErr eCode = ieCode;

  if ( ieCode < xList_tErr_NoErr
       || ieCode >= xList_knNumErrorCodes )
  {

    eCode = xList_tErr_InvalidErrorCode;
  }

  return xList_ksaError [ eCode ];
}

