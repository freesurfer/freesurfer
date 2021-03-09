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
//      File: myutil.c
//      A utility library
//
//--------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "myutil.h"

/*---------------------------------------------------------------------------
// Construct an empty MYlist with specified capacity and capacity increment
//-------------------------------------------------------------------------*/
MYlist myList2(int elementSize, int capacity, int capacityIncrement) {
  MYlist list;
  void * data;

  if (elementSize<1) {
    fprintf(stderr, "myList(): elementSize must be a postive integer!\n");
    exit(1);
  }
  if (capacity<0) {
    fprintf(stderr, "myList(): capacity must not be negative!\n");
    exit(1);
  }

  list = (MYlist)myMalloc(sizeof(MYlistStruct));

  myListSetElementSize(list, elementSize);
  myListSetSize(list, 0);
  myListSetCapacity(list, capacity);
  myListSetCapacityIncrement(list, capacityIncrement);
  if (capacity == 0)
    myListSetData(list, NULL);
  else {
    data = (void *)myMalloc(elementSize*capacity);
    myListSetData(list, data);
  }

  return(list);
}

/*---------------------------------------------------------------------------
// Construct an empty MYlist with specified capacity and default
// capacity increment as 100
//-------------------------------------------------------------------------*/
MYlist myList1(int elementSize, int capacity) {
  MYlist list;

  list = myList2(elementSize, capacity, 100);

  return(list);
}

/*---------------------------------------------------------------------------
// Construct an empty MYlist with default capacity as 0 and
// capacity increment as 100
//-------------------------------------------------------------------------*/
MYlist myList(int elementSize) {
  MYlist list;

  list = myList2(elementSize, 0, 100);

  return(list);
}


/*---------------------------------------------------------------------------
// Construct an empty MYlist with specified size, all the elements are set to
// zero
//-------------------------------------------------------------------------*/
MYlist myListOfSize(int size, int elementSize) {
  MYlist list;
  char *data;
  int i;
  int capacity, capacityIncrement;

  if (size<0) {
    fprintf(stderr, "myListOfSize(): size must not be negative!\n");
    exit(1);
  }

  capacity = size;
  capacityIncrement = 100;
  list = myList2(elementSize, capacity, capacityIncrement);
  myListSetSize(list, size);
  data = (char *)myListData(list);
  for (i=0; i<elementSize*size; i++)
    data[i] = 0;

  return(list);
}

/*---------------------------------------------------------------------------
// Delete this list
//-------------------------------------------------------------------------*/
void myListDelete(MYlist list) {
  void *data;

  data = myListData(list);
  free(data);
  free(list);
}

/*---------------------------------------------------------------------------
// Add an element to this list
//-------------------------------------------------------------------------*/
void myListAddElement(MYlist list, void *element) {
  int size, capacity, elementSize, capacityIncrement;
  void *data;

  size        = myListSize(list);
  capacity    = myListCapacity(list);
  elementSize = myListElementSize(list);
  data        = myListData(list);
  if (size >= capacity) {
    capacityIncrement = myListCapacityIncrement(list);
    capacity += capacityIncrement;
    myListSetCapacity(list, capacity);
    if (data == NULL) {
      /* initial list */
      data = (void *)myMalloc(elementSize * capacity);
    } else {
      /* allocate a larger list */
      data = (void *)myRealloc(data, elementSize*capacity);
    }
    myListSetData(list, data);
  }

  memcpy((char *)data+size*elementSize, (char *)element, elementSize);
  myListSetSize(list, size+1);
}

/*---------------------------------------------------------------------------
// Add an integer to this list (must be a list consists of only integers)
//-------------------------------------------------------------------------*/
void myListAddInt(MYlist list, int element) {
  int size, capacity, elementSize, capacityIncrement;
  int *data;

  size        = myListSize(list);
  capacity    = myListCapacity(list);
  elementSize = myListElementSize(list);
  data        = myListData(list);

  if (size >= capacity) {
    capacityIncrement = myListCapacityIncrement(list);
    capacity += capacityIncrement;
    myListSetCapacity(list, capacity);
    if (data == NULL) {
      /* initial list */
      data = (int *)myMalloc(elementSize * capacity);
    } else {
      /* allocate a larger list */
      data = (int *)myRealloc(data, elementSize*capacity);
    }
    myListSetData(list, data);
  }

  data[size] = element;
  myListSetSize(list, size+1);
}

/*---------------------------------------------------------------------------
// Add an array to this list
//-------------------------------------------------------------------------*/
void    myListAddArray(MYlist list, void *array, int num) {
  int size, capacity, elementSize, capacityIncrement, actualIncrement;
  void *data;

  size        = myListSize(list);
  capacity    = myListCapacity(list);
  elementSize = myListElementSize(list);
  data        = myListData(list);
  if (size + num > capacity) {
    capacityIncrement = myListCapacityIncrement(list);
    actualIncrement = (capacityIncrement > num)? capacityIncrement: num;
    capacity += actualIncrement;
    myListSetCapacity(list, capacity);

    if (data == NULL) {
      /* initial list */
      data = (void *)myMalloc(elementSize * capacity);
    } else {
      /* allocate a larger list */
      data = (void *)myRealloc(data, elementSize*capacity);
    }
    myListSetData(list, data);
  }

  memcpy((char *)data+size*elementSize, (char *)array, num*elementSize);
  myListSetSize(list, size+num);
}

/*---------------------------------------------------------------------------
// Insert an element into the list at the specified index
//-------------------------------------------------------------------------*/
int myListInsertElementAt(MYlist list, int index, void *element) {
  int size, elementSize;
  void *data;
  void *tempPtr;
  char *currentPtr, *nextPtr;
  int i;

  size        = myListSize(list);
  elementSize = myListElementSize(list);

  if (index<0 || index>size-1) {
    return(-1); /* out of bound error */
  }

  tempPtr = (void *)myMalloc(elementSize);
  myListAddElement(list, tempPtr);

  data        = myListData(list);

  for (i=size-1; i>=index; i--) {
    currentPtr = (char *)data+i*elementSize;
    nextPtr    = (char *)currentPtr + elementSize;
    memcpy(nextPtr, currentPtr, elementSize);
  }

  memcpy((char *)data+index*elementSize, (char *)element, elementSize);

  return(0);
}



/*---------------------------------------------------------------------------
// Retrieve an element from this list at a given index
//-------------------------------------------------------------------------*/
int myListElementAt(MYlist list, int index, void *element) {
  int size, elementSize;
  void *data;

  size        = myListSize(list);
  elementSize = myListElementSize(list);
  data        = myListData(list);

  if (index<0 || index>size-1) {
    return(-1); /* out of bound error */
  }
  memcpy((char *)element, (char *)data+index*elementSize, elementSize);

  return(0);
}

/*---------------------------------------------------------------------------
// Sets a list element at a given index
//-------------------------------------------------------------------------*/
int myListSetElementAt(MYlist list, int index, void *element) {
  int size, elementSize;
  void *data;

  size        = myListSize(list);
  elementSize = myListElementSize(list);
  data        = myListData(list);

  if (index<0 || index>size-1) {
    return(-1); /* out of bound error */
  }

  memcpy((char *)data+index*elementSize, (char *)element, elementSize);

  return(0);
}

/*---------------------------------------------------------------------------
// Removes all elements from this list and sets its size to zero
//-------------------------------------------------------------------------*/
int myListRemoveElementAt(MYlist list, int index) {
  int size, elementSize;
  void *data;
  char *currentPtr, *nextPtr;
  int i;

  size        = myListSize(list);
  elementSize = myListElementSize(list);
  data        = myListData(list);

  if (index<0 || index>size-1) {
    return(-1); /* out of bound error */
  }

  for (i=index; i<size-1; i++) {
    currentPtr = (char *)data+i*elementSize;
    nextPtr    = (char *)currentPtr + elementSize;
    memcpy(currentPtr, nextPtr, elementSize);
  }

  myListSetSize(list, size-1);

  return(0);
}

/*---------------------------------------------------------------------------
// Removes all elements from this list and sets its size to zero
//-------------------------------------------------------------------------*/
void myListRemoveAllElements(MYlist list) {
  myListSetSize(list, 0);
}

/*---------------------------------------------------------------------------
// Trim this list to current size
//-------------------------------------------------------------------------*/
void myListTrim(MYlist list) {
  void *data;
  int size, elementSize;

  size        = myListSize(list);
  elementSize = myListElementSize(list);
  data = myListData(list);

  data = (void *)myRealloc(data, elementSize*size);
  myListSetData(list, data);
  myListSetCapacity(list, size);
}

void myListInfo(MYlist list) {
  int elementSize, size, capacity, capacityIncrement;

  elementSize = myListElementSize(list);
  size        = myListSize(list);
  capacity    = myListCapacity(list);
  capacityIncrement = myListCapacityIncrement(list);

  printf("         elementSize = %d\n", elementSize);
  printf("                size = %d\n", size);
  printf("            capacity = %d\n", capacity);
  printf("   capacityIncrement = %d\n", capacityIncrement);
  printf("\n");
}


/*---------------------------------------------------------------------------
// pop out the top element from the stack
//-------------------------------------------------------------------------*/
void myStackPop(MYstack stack, void *element) {
  int size, elementSize;
  void *data;

  size        = myListSize(stack);
  elementSize = myListElementSize(stack);
  data        = myListData(stack);

  memcpy((char *)element, (char *)data+(size-1)*elementSize, elementSize);
  myListSetSize(stack, size-1);
}

/*---------------------------------------------------------------------------
// Construct an empty MYqueue with specified capacity and capacity increment
//-------------------------------------------------------------------------*/
MYqueue myQueue2(int elementSize, int capacity, int capacityIncrement) {
  MYqueue queue;

  queue = (MYqueue)myMalloc(sizeof(MYqueueStruct));
  queue->start = 0;
  queue->end   = -1;
  queue->list  = myList2(elementSize, capacity, capacityIncrement);

  return(queue);
}

/*---------------------------------------------------------------------------
// Construct an empty MYqueue with specified capacity and default
// capacity increment as 100
//-------------------------------------------------------------------------*/
MYqueue myQueue1(int elementSize, int capacity) {
  MYqueue queue;

  queue = myQueue2(elementSize, capacity, 100);

  return(queue);
}

/*---------------------------------------------------------------------------
// default constructor
//-------------------------------------------------------------------------*/
MYqueue myQueue(int elementSize) {
  MYqueue queue;

  queue = myQueue2(elementSize, 0, 100);

  return(queue);
}

/* private functinos */

void    myQueueEnsureSize(MYqueue queue) {
  myListSetSize(queue->list, myQueueSize(queue));
}

/*---------------------------------------------------------------------------
// delete the queue and its memory
//-------------------------------------------------------------------------*/
void    myQueueDelete(MYqueue queue) {
  myListDelete(queue->list);
  free(queue);
}

/*---------------------------------------------------------------------------
// removes all elements from the queue without releasing memory
//-------------------------------------------------------------------------*/
void    myQueueRemoveAllElements(MYqueue queue) {
  queue->start = 0;
  queue->end   = -1;
  myQueueEnsureSize(queue);
}

/* a private function for clean queue, ie. move things to the front */
void    myQueueMoveToFront(MYqueue queue) {
  void *data;
  void *s1, *s2;
  int elementSize, start, end;

  elementSize = myListElementSize(queue->list);
  data  = myListData(queue->list);
  start = queue->start;
  end   = queue->end;
  s2 = (char*)data + start*elementSize;
  s1 = data;
  memmove(s1, s2, (end - start + 1)*elementSize);
  queue->end   = end - start;
  queue->start = 0;
  myQueueEnsureSize(queue);
}


/*---------------------------------------------------------------------------
// push an element to the end of the queue
//-------------------------------------------------------------------------*/
void    myQueuePush(MYqueue queue, void* element) {
  int size, capacity;
  int q = MY_QUEUE_Q;
  /* this is the factor determines the frequency of cleaning */
  /* a factor of n denotes that (n-1)*memory(queue) will be
  wasted */

  size     = myQueueSize(queue);
  capacity = myListCapacity(queue->list);

  if (queue->end >= (capacity - 1) && q*size < capacity )
    /* move the block to front, and release more memory */
    myQueueMoveToFront(queue);

  /* just keep adding the element */
  queue->end = queue->end + 1;
  myListAddElement(queue->list, element);
}

/*---------------------------------------------------------------------------
// push an array to the end of the queue
//-------------------------------------------------------------------------*/
void    myQueuePushArray(MYqueue queue, void* array, int num) {
  int size, capacity;
  int q = MY_QUEUE_Q;
  /* this is the factor determines the frequency of cleaning */
  /* a factor of n denotes that (n-1)*memory(queue) will be
  wasted */

  size     = myQueueSize(queue);
  capacity = myListCapacity(queue->list);

  if (queue->end >= (capacity - 1) && q*size < capacity )
    /* move the block to front, and release more memory */
    myQueueMoveToFront(queue);

  /* just keep adding the array */
  queue->end = queue->end + num;
  myListAddArray(queue->list, array, num);
}

/*---------------------------------------------------------------------------
// pop out the first element from the queue
//-------------------------------------------------------------------------*/
int myQueuePop(MYqueue queue, void* element) {
  int start;

  start = queue->start;
  if (myQueueIsEmpty(queue) == 0) {
    myListElementAt(queue->list, start, element); /* get the first one */
    queue->start = start + 1;
    return(1); /* correctly got the result */
  } else {
    return(0); /* element is undefined */
  }
}

/*---------------------------------------------------------------------------
// trim the queue to its current size
//-------------------------------------------------------------------------*/
void    myQueueTrim(MYqueue queue) {
  myQueueMoveToFront(queue);

  myListTrim(queue->list);
}

/*---------------------------------------------------------------------------
// extract an array from the queue
//-------------------------------------------------------------------------*/
void*   myQueueToArray(MYqueue queue) {
  void *arr, *data;
  int size, elementSize;

  size = myQueueSize(queue);
  elementSize = myListElementSize(queue->list);

  if (size == 0) {
    arr = NULL;
  } else {
    data = myListData(queue->list);
    arr = (void *)myMalloc(elementSize*size);
    memcpy(arr, (char*)data + queue->start*elementSize, size*elementSize);
  }

  return(arr);
}


void    myQueueInfo(MYqueue queue) {
  printf("               start = %d\n", queue->start);
  printf("               end   = %d\n", queue->end);
  myListInfo(queue->list);
}


void* myMalloc(int size) {
  register void *value = (void *)malloc (size);
  if (value == 0)
    myError("Virtual memory exhausted!");

  return value;
}

void* myRealloc(void* ptr, int size) {
  register void *value = (void *)realloc (ptr, size);
  if (value == 0)
    myError("Virtual memory exhausted!");

  return value;
}

/* Prints error message */
void myError(char error_text[]) {
  fprintf(stderr,"Utility run-time error:\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"Now exiting to system.\n");
  exit(1);
}
