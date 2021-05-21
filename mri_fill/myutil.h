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
//      File: myutil.h
//      A MY utility library
//
//--------------------------------------------------------------------------*/

/*============================================================================
//                        MYutil Summary
//============================================================================
//
//  DATA STRUCTURE:
//  MYlist ---> {elementSize, size, data, capacity, capacityIncrement}
//
//  NOTE: the type of MYlist is a pointer.
//
//  Assume the following variable declarations:
//      MYlist list, list1, list2;
//      ListElement element;
//      ListElement *data;
//      int capacity, capacityIncrement, size,
//          elementSize, index, returnFlag;
//
//  LIST MANIPULATIONS:
//
//  myList ........................ default constructor
//           list  = myList(sizeof(ListElement));
//
//  myList1 ....................... default constructor with capacity
//           list1 = myList(sizeof(ListElement), capacity);
//
//  myList2 ....................... default constructor with capacity
//                                  and capacity increment
//           list2 = myList(sizeof(ListElement), capacity, capacityIncrement);
//
//  myListOfSize .................. default constructor with a specified size
//           list  = myListOfSize(size, sizeof(ListElement));
//
//  myListDelete .................. delete the list
//           myListDelete(list);
//
//  myListAddElement .............. add an element to this list
//           myListAddElement(list, &element);
//
//  myListInsertElementAt ......... insert an element in the list at the given
//                                  index. Each list element's index greater
//                                  or equal to the specified index is
//                                  shifted upward than its previous value.
//           returnFlag = myListInsertElementAt(list, index, &element);
//
//  myListElementAt ............... retrieve an element at index
//           returnFlag = myListElementAt(list, index, &element);
//
//  myListSetElementAt ............ set the element at the specified index of
//                                  this list by copying the value of
//                                  given element.
//           returnFlag = myListSetElementAt(list, index, &element);
//
//  myListRemoveElementAt ......... Delete the element at the specified index.
//                                  The index of each element after the
//                                  specified index is decreased by 1.
//           returnFlag = myListRemoveElementAt(list, index);
//
//  myListRemoveAllElements ....... removes all elements from this list
//                                  and sets its size to zero
//           myListRemoveAllElements(list);
//
//  myListTrim .................... trim this list to its current size
//           myListTrim(list);
//
//  myListIsEmpty ................. 1 if this list has no elements
//           returnFlag = myListIsEmpty(list);
//
//  myListElementSize ............. the element size of each component
//           elementSize = myListElementSize(list);
//
//  myListSize .................... the current size of this list
//           size = myListSize(list);
//
//  myListData .................... the data of this list
//           data = myListData(list);
//
//  DATA STRUCTURE:
//  MYstack ---> implemented on top of MYlist
//
//  NOTE: the type of MYstack is a pointer.
//
//  STACK MANIPULATIONS:
//  myStack ....................... default constructor
//  myStackPush ................... push an element into the stack
//  myStackPop .................... pop out the top element from the stack
//  myStackIsEmpty ................ 1 if this stack is empty
//  myStackRemoveAllElements ...... removes all elements from the stack
//                                  and sets its size to zero
//  myStackTrim ................... trim the stack to its current size
//  myStackDelete ................. delete the stack
//  myStackElementSize ............ the element size of each component
//
//  QUEUE DATA STRUCTURE:
//  MYqueue ---> implemented on top of MYlist
//
//  QUEUE  MANIPULATIONS:
//  myQueue ....................... default constructor
//  myQueue1 ...................... constructor 1
//  myQueue2 ...................... constructor 2
//  myQueueDelete ................. delete the queue and its memory
//  myQueueRemoveAllElements ...... removes all elements from the queue
//                                  without releasing memory
//  myQueuePush ................... push an element to the end of the queue
//  myQueuePop .................... pop out the first element from the queue
//  myQueueTrim ................... trim the queue to its current size
//  myQueueToArray ................ extract an array from the queue
//  myQueueSize ................... the current size of the queue
//  myQueueIsEmpty ................ 1 if this stack is empty
//  myQueueInfo ................... print the queue information
//  myQueueElementSize ............ returns the size of an element
//==========================================================================*/

#ifndef MYUTIL_TOOLS
#define MYUTIL_TOOLS
#include <stdio.h>

typedef struct {
  int   elementSize;
  int   size;
  void *data;
  int   capacity;
  int   capacityIncrement;
}
MYlistStruct;

typedef MYlistStruct *MYlist;
#define MYstack MYlist

typedef struct {
  int start;
  int end;
  MYlist list;
}
MYqueueStruct;

typedef MYqueueStruct *MYqueue;

#define MY_QUEUE_Q    2


#define myListIsEmpty(list)            ((list)->size == 0 ? 1: 0)
#define myListSize(list)               ((list)->size)
#define myListData(list)               ((list)->data)
#define myListElementSize(list)        ((list)->elementSize)

void* myMalloc(int size);
void* myRealloc(void* ptr, int size);
void myError(char error_text[]);

MYlist  myList(int elementSize);
MYlist  myList1(int elementSize, int capacity);
MYlist  myList2(int elementSize, int capacity, int capacityIncrement);
MYlist  myListOfSize(int size, int elementSize);
void    myListDelete(MYlist list);
void    myListAddElement(MYlist list, void *element);
void    myListAddInt(MYlist list, int element);
void    myListAddArray(MYlist list, void *array, int num);
int myListInsertElementAt(MYlist list, int index, void *element);
int myListSetElementAt(MYlist list, int index, void *element);
int myListElementAt(MYlist list, int index,
                    /* stores the result at */ void *element);
int myListRemoveElementAt(MYlist list, int index);
void    myListRemoveAllElements(MYlist list);
void    myListTrim(MYlist list);
void    myListInfo(MYlist list);

#define myStack(elementSize)       myList(elementSize)
#define myStack1(elementSize,capacity)  myList(elementSize,capacity)
#define myStack2(elementSize,capacity,capacityIncrement) myList(elementSize,capacity,capacityIncrement)
#define myStackPush(stack, element) myListAddElement(stack, element)
void    myStackPop(MYstack stack, /* stores the result at */ void *element);
#define myStackIsEmpty(stack)       myListIsEmpty(stack)
#define myStackRemoveAllElements(stack) myListRemoveAllElements(stack)
#define myStackTrim(stack)          myListTrim(stack)
#define myStackDelete(stack)        myListDelete(stack)
#define myStackElementSize(stack)  (stack->elementSize)

MYqueue myQueue2(int elementSize, int capacity, int capacityIncrement);
MYqueue myQueue1(int elementSize, int capacity);
MYqueue myQueue(int elementSize);
void    myQueueDelete(MYqueue queue);
void    myQueueRemoveAllElements(MYqueue queue);
void    myQueuePush(MYqueue queue, void* element);
void    myQueuePushArray(MYqueue queue, void* array, int num);
int myQueuePop(MYqueue queue, void* element);
void    myQueueTrim(MYqueue queue);
void*   myQueueToArray(MYqueue queue);
void    myQueueInfo(MYqueue queue);
#define myQueueSize(queue)         (queue->end - queue->start + 1)
#define myQueueIsEmpty(queue)      (myQueueSize(queue) == 0)
#define myQueueElementSize(queue)  (myListElementSize(queue->list))

/*============================================================================
//                        MYutil private functions
//                        USER PLEASE DO NOT USE
// NOTE: these numbers should not be changed once the list is created
//============================================================================
//  myListCapacity ................ the current capacity of this list
//  myListCapacityIncrement ....... the current capacityIncrement of this list
//  myListSetElementSize .......... set the element size of this list
//  myListSetSize ................. set the size of this list
//  myListSetData ................. set the data of this list
//  myListSetCapacity ............. set the capacity of this list
//  myListSetCapacityIncrement .... set the capacityIncrement of this list
//  myListInfo .................... print the information about the list
//==========================================================================*/

/*--------------------------------------------------------------------------
// Private functions
//------------------------------------------------------------------------*/
#define myListCapacity(list)                  ((list)->capacity)
#define myListCapacityIncrement(list)         ((list)->capacityIncrement)
#define myListSetElementSize(list, val)       (list)->elementSize = val
#define myListSetSize(list, val)              (list)->size = val
#define myListSetData(list, val)              (list)->data = (void *)val
#define myListSetCapacity(list, val)          (list)->capacity = val
#define myListSetCapacityIncrement(list, val) (list)->capacityIncrement = val

#endif
