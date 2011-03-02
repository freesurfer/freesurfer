/**
 * @file  xList.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.7 $
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


#ifndef xList_H
#define xList_H

#include "xTypes.h"

/* various error codes that will be
returned by all functions. these can be
translated into human readable strings
with the GetErrorString func below. note
that some of these code are not
necessarily errors. if you ask for the
next item and get a EndOfList error, it
just means you're at the end of the
list. */
typedef enum
{

  xList_tErr_NoErr = 0,                 // no error
  xList_tErr_AllocationFailed,          // local memory allocation failed
  xList_tErr_InternalAllocationFailed,  // rts list new failed
  xList_tErr_InternalDeletionFailed,    // rts list delete failed
  xList_tErr_DeletionFailed,            // local delete failed
  xList_tErr_InvalidListRef,            // null list arg ptr
  xList_tErr_ListFull,                  // couldn't insert item
  xList_tErr_ItemNotInList,             // item not found
  xList_tErr_ErrorFindingItem,          // couldn't look for item
  xList_tErr_ErrorGettingSize,          // couldn't get size
  xList_tErr_EndOfList,                 // end of list
  xList_tErr_ListEmpty,                 // list is empty
  xList_tErr_InvalidErrorCode,
  xList_knNumErrorCodes

} xList_tErr;

typedef enum
{

  xList_tCompare_LessThan    = -1,
  xList_tCompare_Match       = 0,
  xList_tCompare_GreaterThan = 1
} xList_tCompare;


/* we use a signature to detect corrupted
memory or multiple deletes on the same
ptr. */
#define xList_kSignature 0x8f7f6f5c

/* private list node struct. */
typedef struct
{

  long               mSignature;   // the signature
  void*              mpData;       // pointer to the data in the list
  void*              mpNext;       // ptr to the next node
}
xListNode, *xListNodeRef;

/* private list struct. */
typedef struct
{

  long               mSignature;   // the signature
  xListNodeRef       mpHead;       // first node
  xListNodeRef       mpTail;       // last node
  xListNodeRef       mpNext;       // the next node in a traversal

  xList_tCompare(*mComparator)(void*,void*);

}
xList, *xListRef;

/* public interface */

/* desc:
call to create a new list. pass in
an uninited xListRef and get back an
inited list. allocated a single node.
* returns:
AllocationFailed: couldn't allocate
necessary memory. oppList will be nil.
*/
xList_tErr xList_New ( xListRef* oppList );

/* desc:
deletes the list. note that this will
leave orphaned ptrs.
* returns:
InvalidListRef: null list ptr or bad
signature.
*/
xList_tErr xList_Delete ( xListRef* ioppList );

/* desc:
inserts an item onto the tail of the
list.
* returns:
InvalidListRef: null list ptr or bad
signature.
*/
xList_tErr xList_InsertItem ( xListRef this,
                              void*    ipItemToInsert );

/* desc:
removes an item from the list. note
that if a comparator is defined, the
ptr that is passed in may not be the
same as the ptr in the list; both can
point to different pieces of memory.
since the one in the list needs to be
removed, it will be returned in the
item paramter. it should then be
deleted.
* returns:
InvalidListRef: null list ptr or bad
signature.
InternalAllocationFailed: couldn't
allocate node.
*/
xList_tErr xList_RemoveItem ( xListRef this,
                              void**   iopItemToRemove );

/* desc:
searches list for item, returns true if
found.
* returns:
InvalidListRef: null list ptr or bad
signature.
ItemNotInList: couldn't find
ipItemToRemove.
*/
xList_tErr xList_IsInList ( xListRef this,
                            void*    ipItem,
                            tBoolean*    obpIsInList );

/* desc:
empties the list. only deletes nodes,
not the data they're pointing to.
* returns:
InvalidListRef: null list ptr or bad
signature.
*/
xList_tErr xList_Clear ( xListRef this );

/* desc:
returns the number of items in the
list.
* returns:
InvalidListRef: null list ptr or bad
signature for list or node.
*/
xList_tErr xList_GetCount ( xListRef this,
                            int*     opCount );

/* desc:
item by item acccess. gets first item
at head of list. for simple list
traversing, it's quicker to use the
ResetPosition and
GetNextItemFromPosition functions.
* returns:
InvalidListRef: null list ptr or bad
signature.
*/
xList_tErr xList_GetFirstItem ( xListRef this,
                                void**   oppFirstItem );

/* desc: given an item, searches the list
for the item, gets the next node, and
returns the data at that node.
* returns:
InvalidListRef: null list ptr or bad
signature.
*/
xList_tErr xList_GetNextItem ( xListRef this,
                               void*    ipCurrentItem,
                               void**   oppNextItem );

/* desc:
quicker for simple sequential access.
resets internal position to the head of
the list. call before traversing.
* returns:
InvalidListRef: null list ptr or bad
signature.
EndOfList: no next item. oppNextItem
will be null.
ItemNotInList: couldn't find
ipCurrentItem.
*/
xList_tErr xList_ResetPosition ( xListRef this );

/* desc:
gets the data at the current internal
position and advances the internal
ptr.
* returns:
InvalidListRef: null list ptr or bad
signature.
*/
xList_tErr xList_GetNextItemFromPosition ( xListRef this,
    void**   oppNextItem );
xList_tErr xList_NextFromPos ( xListRef this,
                               void**   oppNextItem );

/* desc:
for stack style access. pushes onto the
head.
* returns:
InvalidListRef: null list ptr or bad
signature for list or node.
EndOfList: no next item. oppNextItem
will be nil.
*/
xList_tErr xList_PushItem ( xListRef this,
                            void*    ipItemToInsert );

/* desc:
for stack style access. pops an item off
the head and returns its data
* returns:
InvalidListRef: null list ptr or bad
signature.
InternalAllocationFailed: couldn't
allocate node.
*/
xList_tErr xList_PopItem ( xListRef this, void** oppItem );

/* desc:
error checking. checks for null ptr and
for a valid signature. every function
calls this one, so you shouldn't have to,
but you can.
* returns:
InvalidListRef: null list ptr or bad
signature.
*/
xList_tErr xList_SetComparator ( xListRef             this,
                                 xList_tCompare(*iComparator)(void*,void*) );

xList_tCompare xList_CompareItems_ ( xListRef this,
                                     void*    pItemA,
                                     void*    pItemB );

xList_tErr xList_Verify ( xListRef this );

/* desc:
gives back an error string.
* returns:
ptr to a string.
*/
char* xList_GetErrorString ( xList_tErr ieCode );

/* private functions */

/* desc:
internal method for finding the node
containing a specific ptr.
* returns:
ptr to the node or nil if not found.
*/
xListNodeRef xList_FindItem_ ( xListRef this,
                               void*    ipItem );


#endif
