/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


/*----------------------------------------------------------------------------
//
//      File: heap.h (renamed from pgutil.h)
//      Data Structure: list and heap
//
//--------------------------------------------------------------------------*/

/*============================================================================
//                        Summary
//============================================================================
//
//  DATA STRUCTURE:
//  PGlist ---> {elementSize, size, data, capacity, capacityIncrement}
//
//  NOTE: the type of PGlist is a pointer.
//
//  Assume the following variable declarations:
//      PGlist list, list1, list2;
//      ListElement element;
//      ListElement *data;
//      int capacity, capacityIncrement, size,
//          elementSize, index, returnFlag;
//
//  LIST MANIPULATIONS:
//
//  pgList ........................ default constructor
//           list  = pgList(sizeof(ListElement));
//
//  pgList1 ....................... default constructor with capacity
//           list1 = pgList(sizeof(ListElement), capacity);
//
//  pgList2 ....................... default constructor with capacity
//                                  and capacity increment
//           list2 = pgList(sizeof(ListElement), capacity, capacityIncrement);
//
//  pgListOfSize .................. default constructor with a specified size
//           list  = pgListOfSize(size, sizeof(ListElement));
//
//  pgListDelete .................. delete the list
//           pgListDelete(list);
//
//  pgListAddElement .............. add an element to this list
//           pgListAddElement(list, &element);
//
//  pgListInsertElementAt ......... insert an element in the list at the given
//                                  index. Each list element's index greater
//                                  or equal to the specified index is
//                                  shifted upward than its previous value.
//           returnFlag = pgListInsertElementAt(list, index, &element);
//
//  pgListElementAt ............... retrieve an element at index
//           returnFlag = pgListElementAt(list, index, &element);
//
//  pgListSetElementAt ............ set the element at the specified index of
//                                  this list by copying the value of
//                                  given element.
//           returnFlag = pgListSetElementAt(list, index, &element);
//
//  pgListRemoveElementAt ......... Delete the element at the specified index.
//                                  The index of each element after the
//                                  specified index is decreased by 1.
//           returnFlag = pgListRemoveElementAt(list, index);
//
//  pgListRemoveAllElements ....... removes all elements from this list
//                                  and sets its size to zero
//           pgListRemoveAllElements(list);
//
//  pgListTrim .................... trim this list to its current size
//           pgListTrim(list);
//
//  pgListIsEmpty ................. PG_TRUE if this list has no elements
//           returnFlag = pgListIsEmpty(list);
//
//  pgListElementSize ............. the element size of each component
//           elementSize = pgListElementSize(list);
//
//  pgListSize .................... the current size of this list
//           size = pgListSize(list);
//
//  pgListData .................... the data of this list
//           data = pgListData(list);
//
//==========================================================================*/

#ifndef PGUTIL_TOOLS
#define PGUTIL_TOOLS
#include <stdio.h>

#ifndef PG_ERROR
#define PG_ERROR -1
#endif

typedef struct
{
  int   elementSize;
  int   size;
  void *data;
  int   capacity;
  int   capacityIncrement;
}
PGlistStruct;

typedef PGlistStruct *PGlist;

#define pgListIsEmpty(list)              (list->size == 0 ? PG_TRUE: PG_FALSE)
#define pgListSize(list)                 (list->size)
#define pgListData(list)                 (list->data)
#define pgListElementSize(list)          (list->elementSize)

PGlist  pgList(int elementSize);
PGlist  pgList1(int elementSize, int capacity);
PGlist  pgList2(int elementSize, int capacity, int capacityIncrement);
PGlist  pgListOfSize(int size, int elementSize);
void    pgListDelete(PGlist list);
void    pgListAddElement(PGlist list, void *element);
int pgListInsertElementAt(PGlist list, int index, void *element);
int pgListSetElementAt(PGlist list, int index, void *element);
int pgListElementAt(PGlist list, int index,
                    /* stores the result at */ void *element);
int pgListRemoveElementAt(PGlist list, int index);
void    pgListRemoveAllElements(PGlist list);
void    pgListTrim(PGlist list);

/*============================================================================
//                        PGutil private functions
//                        USER PLEASE DO NOT USE
// NOTE: these numbers should not be changed once the list is created
//============================================================================
//  pgListCapacity ................ the current capacity of this list
//  pgListCapacityIncrement ....... the current capacityIncrement of this list
//  pgListSetElementSize .......... set the element size of this list
//  pgListSetSize ................. set the size of this list
//  pgListSetData ................. set the data of this list
//  pgListSetCapacity ............. set the capacity of this list
//  pgListSetCapacityIncrement .... set the capacityIncrement of this list
//  pgListInfo .................... print the information about the list
//==========================================================================*/

/*--------------------------------------------------------------------------
// Private functions
//------------------------------------------------------------------------*/
#define pgListCapacity(list)                  (list->capacity)
#define pgListCapacityIncrement(list)         (list->capacityIncrement)
#define pgListSetElementSize(list, val)       list->elementSize = val
#define pgListSetSize(list, val)              list->size = val
#define pgListSetData(list, val)              list->data = (void *)val
#define pgListSetCapacity(list, val)          list->capacity = val
#define pgListSetCapacityIncrement(list, val) list->capacityIncrement = val

/* The following is the definition of HEAP data structure,
 which is built upon pgList above. It implements a Min_heap that stores
 floating point value with an associated attribute of id */
typedef struct
{
  double value;
  int id;
  int*  p; /* backpointer */
}
XheapElement;

typedef PGlist Xheap;

int   xhSize(Xheap H);   /* get size of heap */
Xheap xhInitEmpty();     /* an empty heap */
Xheap xhInit(XheapElement *array, int N); /* init from an array (0,N-1) */
void  xhDestroy(Xheap H); /* destroy the heap and free the memory */
int   xhUpHeap(int k, Xheap H);
int   xhDownHeap(int k, Xheap H);
int   xhInsert(double value, int id, int *p, Xheap H);
XheapElement xhRemove(Xheap H);
XheapElement xhReplace(double value, int id, int *p, Xheap H);
XheapElement xhDelete(int k, Xheap H);
XheapElement xhChange(int k, double value, int id, int *p, Xheap H);
XheapElement xhChangeValue(int k, double value, Xheap H);
XheapElement xhGet(int k, Xheap H); /* k must be 1, 2, ... N */
#define xhIsEmpty(H)  (xhSize(H) == 0)

#endif
