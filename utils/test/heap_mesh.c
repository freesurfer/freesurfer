/**
 * @brief list utils
 *
 */
/*
 * Original Author: Xaio Han
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

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "heap.h"

#define INCREMENT 1000
#define PG_OK 0

/*---------------------------------------------------------------------------
// Construct an empty PGlist with specified capacity and capacity increment
//-------------------------------------------------------------------------*/
PGlist pgList2(int elementSize, int capacity, int capacityIncrement)
{
  PGlist list;
  void * data;

  if (elementSize<1)
  {
    fprintf(stderr, "pgList(): elementSize must be a postive integer!\n");
    exit(1);
  }
  if (capacity<0)
  {
    fprintf(stderr, "pgList(): capacity must not be negative!\n");
    exit(1);
  }

  list = (PGlist)malloc(sizeof(PGlistStruct));

  pgListSetElementSize(list, elementSize);
  pgListSetSize(list, 0);
  pgListSetCapacity(list, capacity);
  pgListSetCapacityIncrement(list, capacityIncrement);
  if (capacity == 0)
    pgListSetData(list, NULL);
  else
  {
    data = (void *)malloc(elementSize*capacity);
    pgListSetData(list, data);
  }

  return(list);
}

/*---------------------------------------------------------------------------
// Construct an empty PGlist with specified capacity and default
// capacity increment as 100
//-------------------------------------------------------------------------*/
PGlist pgList1(int elementSize, int capacity)
{
  PGlist list;

  list = pgList2(elementSize, capacity, 100);

  return(list);
}

/*---------------------------------------------------------------------------
// Construct an empty PGlist with default capacity as 0 and
// capacity increment as 100
//-------------------------------------------------------------------------*/
PGlist pgList(int elementSize)
{
  PGlist list;

  list = pgList2(elementSize, 0, 100);

  return(list);
}


/*---------------------------------------------------------------------------
// Construct an empty PGlist with specified size, all the elements are set to
// zero
//-------------------------------------------------------------------------*/
PGlist pgListOfSize(int size, int elementSize)
{
  PGlist list;
  char *data;
  int i;
  int capacity, capacityIncrement;

  if (size<0)
  {
    fprintf(stderr, "pgListOfSize(): size must not be negative!\n");
    exit(1);
  }

  capacity = size;
  capacityIncrement = 100;
  list = pgList2(elementSize, capacity, capacityIncrement);
  pgListSetSize(list, size);
  data = (char *)pgListData(list);
  for (i=0; i<elementSize*size; i++)
    data[i] = 0;

  return(list);
}

/*---------------------------------------------------------------------------
// Delete this list
//-------------------------------------------------------------------------*/
void pgListDelete(PGlist list)
{
  void *data;

  data = pgListData(list);
  free(data);
  free(list);
}

/*---------------------------------------------------------------------------
// Add an element to this list
//-------------------------------------------------------------------------*/
void pgListAddElement(PGlist list, void *element)
{
  int size, capacity, elementSize, capacityIncrement;
  void *data;

  size        = pgListSize(list);
  capacity    = pgListCapacity(list);
  elementSize = pgListElementSize(list);
  data        = pgListData(list);
  if (size >= capacity)
  {
    capacityIncrement = pgListCapacityIncrement(list);
    if (data == NULL)
    {
      /* initial list */
      capacity += capacityIncrement;
      pgListSetCapacity(list, capacity);
      data = (void *)malloc(elementSize * capacity);
    }
    else
    {
      /* allocate a larger list */
      capacity += capacityIncrement;
      pgListSetCapacity(list, capacity);
      data = (void *)realloc(data, elementSize*capacity);
    }
    pgListSetData(list, data);
  }

  memmove((char *)data+size*elementSize, (char *)element, elementSize);
  pgListSetSize(list, size+1);
}

/*---------------------------------------------------------------------------
// Insert an element into the list at the specified index
//-------------------------------------------------------------------------*/
int pgListInsertElementAt(PGlist list, int index, void *element)
{
  int size, elementSize;
  void *data;
  void *tempPtr;
  char *currentPtr, *nextPtr;
  int i;

  size        = pgListSize(list);
  elementSize = pgListElementSize(list);

  if (index<0 || index>size-1)
  {
    return(PG_ERROR); /* out of bound error */
  }

  tempPtr = (void *)malloc(elementSize);
  pgListAddElement(list, tempPtr);

  data        = pgListData(list);

  for (i=size-1; i>=index; i--)
  {
    currentPtr = (char *)data+i*elementSize;
    nextPtr    = (char *)currentPtr + elementSize;
    memmove(nextPtr, currentPtr, elementSize);
  }

  memmove((char *)data+index*elementSize, (char *)element, elementSize);

  return(PG_OK);
}



/*---------------------------------------------------------------------------
// Retrieve an element from this list at a given index
//-------------------------------------------------------------------------*/
int pgListElementAt(PGlist list, int index, void *element)
{
  int size, elementSize;
  void *data;

  size        = pgListSize(list);
  elementSize = pgListElementSize(list);
  data        = pgListData(list);

  if (index<0 || index>size-1)
  {
    return(PG_ERROR); /* out of bound error */
  }
  memmove((char *)element, (char *)data+index*elementSize, elementSize);

  return(PG_OK);
}

/*---------------------------------------------------------------------------
// Sets a list element at a given index
//-------------------------------------------------------------------------*/
int pgListSetElementAt(PGlist list, int index, void *element)
{
  int size, elementSize;
  void *data;

  size        = pgListSize(list);
  elementSize = pgListElementSize(list);
  data        = pgListData(list);

  if (index<0 || index>size-1)
  {
    return(PG_ERROR); /* out of bound error */
  }

  memmove((char *)data+index*elementSize, (char *)element, elementSize);

  return(PG_OK);
}

/*---------------------------------------------------------------------------
// Removes all elements from this list and sets its size to zero
//-------------------------------------------------------------------------*/
int pgListRemoveElementAt(PGlist list, int index)
{
  int size, elementSize;
  void *data;
  char *currentPtr, *nextPtr;
  int i;

  size        = pgListSize(list);
  elementSize = pgListElementSize(list);
  data        = pgListData(list);

  if (index<0 || index>size-1)
  {
    return(PG_ERROR); /* out of bound error */
  }

  for (i=index; i<size-1; i++)
  {
    currentPtr = (char *)data+i*elementSize;
    nextPtr    = (char *)currentPtr + elementSize;
    memmove(currentPtr, nextPtr, elementSize);
  }

  pgListSetSize(list, size-1);

  return(PG_OK);
}

/*---------------------------------------------------------------------------
// Removes all elements from this list and sets its size to zero
//-------------------------------------------------------------------------*/
void pgListRemoveAllElements(PGlist list)
{
  pgListSetSize(list, 0);
}

/*---------------------------------------------------------------------------
// Trim this list to current size
//-------------------------------------------------------------------------*/
void pgListTrim(PGlist list)
{
  void *data;
  int size, elementSize;

  size        = pgListSize(list);
  elementSize = pgListElementSize(list);
  data = pgListData(list);

  data = (void *)realloc(data, elementSize*size);
  pgListSetData(list, data);
  pgListSetCapacity(list, size);
}

void pgListInfo(PGlist list)
{
  int elementSize, size, capacity, capacityIncrement;

  elementSize = pgListElementSize(list);
  size        = pgListSize(list);
  capacity    = pgListCapacity(list);
  capacityIncrement = pgListCapacityIncrement(list);

  printf("         elementSize = %d\n", elementSize);
  printf("                size = %d\n", size);
  printf("            capacity = %d\n", capacity);
  printf("   capacityIncrement = %d\n", capacityIncrement);
  printf("\n");
}

Xheap xhInitEmpty()
{
  Xheap H;
  XheapElement he;

  H = (Xheap)pgList2(sizeof(XheapElement),0,INCREMENT);

  /* assume that a[0] = smallest floating point such as -1e33
     serve as sentinel for stopping condition.
     heap starts from 1 */
  he.value = 1e-34;  /* not used in UpHeap anymore */
  he.id = -1;

  pgListAddElement(H, &he);

  return(H);
}

Xheap xhInit(XheapElement *array, int N)
{
  Xheap H;
  XheapElement *data, he;
  int i, *p;

  H = (Xheap)pgList2(sizeof(XheapElement),0,INCREMENT);

  he.value = 1e-34;
  he.id = -1;

  pgListAddElement(H, &he);
  pgListSetSize(H,N+1);

  data = (XheapElement *)pgListData(H);
  for (i = 1; i <= N; i++)
  {
    data[i] = array[i-1];
    p = data[i].p;
    *p = i;
  }
  /* down build */
  for (i = N/2; i >= 1; i--)
    xhDownHeap(i, H);

  return(H);
}

/* destroy the heap and free the memory */
void  xhDestroy(Xheap H)
{
  pgListDelete(H);
}


int xhUpHeap(int k, Xheap H)
{
  XheapElement *a, v;
  int k_father;
  int *p;

  a = (XheapElement *)pgListData(H);

  v = a[k];
  k_father = k/2;  /* integer divsion to retrieve its parent */
  while (k_father > 0 && a[k_father].value > v.value)
  {
    a[k] = a[k_father];
    p = a[k].p;
    *p = k;
    k = k_father;
    k_father = k/2;
  }
  a[k] = v;
  p = a[k].p;
  *p = k;

  return(k);
}

int xhDownHeap(int k, Xheap H)
{
  XheapElement *a, v;
  int N, k_minson;
  int *p;

  a = (XheapElement *)pgListData(H);
  N = xhSize(H);

  v = a[k];
  while ( k <= N/2 )
  {
    k_minson = k+k;
    if ( k_minson < N )
    {
      if (a[k_minson].value > a[k_minson+1].value)
        k_minson = k_minson + 1;  /* always locate the smallest son */
    }
    if ( v.value <= a[k_minson].value )
      break; /* break out the loop */
    a[k] = a[k_minson];
    p = a[k].p;
    *p = k;
    k = k_minson;      /* go down one level */
  }
  a[k] = v;
  p = a[k].p;
  *p = k;

  return(k);
}


int xhInsert(double value, int id, int *p, Xheap H)
{
  XheapElement *a, v;
  int N, k;

  a = (XheapElement *)pgListData(H);

  v.value = value;
  v.id = id;
  v.p = p;

  pgListAddElement(H, &v);
  N = xhSize(H);
  k = xhUpHeap(N, H);

  return(k);
}

/* remove the smallest element */
XheapElement xhRemove(Xheap H)
{
  XheapElement v, *a;
  int N;

  N = xhSize(H);
  a = (XheapElement *)pgListData(H);

  v = a[1];
  a[1] = a[N];
  pgListSetSize(H, N);
  /* the size of list is always 1 more than the size of heap
     since the heap starts at 1 */

  xhDownHeap(1,H);

  return(v);
}

/* replace the smallest value with a new value if the new value is smaller
   otherwise the new value is returned and the heap is unchanged */
XheapElement xhReplace(double value, int id, int* p, Xheap H)
{
  XheapElement *a, v;

  a = (XheapElement *)pgListData(H);

  if ( value < a[1].value )
  {
    v = a[1];
    a[1].value = value;
    a[1].id = id;
    a[1].p = p;
    xhDownHeap(1,H);
  }
  else
  {
    v.value = value;
    v.id = id;
    v.p = p;
  }
  return(v);
}

/* delete an item in the heap and its value is returned */
XheapElement xhDelete(int k, Xheap H)
{
  XheapElement *a, v;
  int N;

  N = xhSize(H);
  a = (XheapElement *)pgListData(H);

  v = a[k];

  a[k] = a[N];
  pgListSetSize(H, N);
  /* the size of list is always 1 more than the size of heap
     since the heap starts at 1 */

  xhDownHeap(k, H);

  return(v);
}

/* change the value of an item and its original value is returned */
XheapElement xhChange(int k, double value, int id, int *p, Xheap H)
{
  XheapElement *a, v;

  a = (XheapElement *)pgListData(H);
  v = a[k];

  if (value != a[k].value)
  {
    a[k].value = value;
    a[k].id = id;
    a[k].p = p;
    if ( value < v.value )
      xhUpHeap(k, H);
    else
      xhDownHeap(k, H);
  }

  return(v);
}

/* change the value of an item and its original value is returned */
XheapElement xhChangeValue(int k, double value, Xheap H)
{
  XheapElement *a, v;

  a = (XheapElement *)pgListData(H);
  v = a[k];

  if (value != a[k].value)
  {
    a[k].value = value;
    if ( value < v.value )
      xhUpHeap(k, H);
    else
      xhDownHeap(k, H);
  }

  return(v);
}

XheapElement xhGet(int k, Xheap H)
{
  XheapElement v;

  pgListElementAt(H, k, &v);

  return(v);
}

int xhSize(Xheap H)
{
  return( pgListSize(H) - 1 );
}

