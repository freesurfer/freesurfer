// MARS Project (Thomas Yeo and Mert Sabuncu), MIT, CSAIL (c) 2006-2008

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "min_heap.h"

// Note that indices goes from 0 to size-1.

static int Min_HeapExchangeElements(MIN_HEAP *MH, int i, int j);
static int Min_HeapifyDown(MIN_HEAP *MH, int index);
static int Min_HeapifyUp(MIN_HEAP *MH, int index);

MIN_HEAP
*Min_HeapAllocate(int max_size, int max_id_array_size)
{
  MIN_HEAP *MH;
  int i;

  MH = (MIN_HEAP *)calloc(1, sizeof(MIN_HEAP));
  if (MH == NULL) {
    fprintf(stderr, "Min_HeapAllocate: Unable to allocate memory for min heap!!\n");
    return (NULL);
  }

  MH->MaxHeapSize = max_size;
  MH->CurrHeapSize = 0;
  MH->MHE_array = (MHE *)calloc(max_size, sizeof(MHE));
  if (MH->MHE_array == NULL) {
    fprintf(stderr, "Min_HeapAllocate: Unable to allocate memory for min heapelement array!!\n");
    return (NULL);
  }

  if (max_id_array_size > 0) {
    MH->id_array = (int *)calloc(max_id_array_size, sizeof(int));
    if (MH->id_array == NULL) {
      fprintf(stderr, "*Min_HeapAllocate: Unable to allocate memory for id array!!\n");
      return (NULL);
    }

    for (i = 0; i < max_id_array_size; i++) MH->id_array[i] = -1;
  }
  else {
    fprintf(stderr, "Min_HeapAllocate: max_id_array_size is 0!!\n");
    return (NULL);
  }
  MH->max_id_array_size = max_id_array_size;

  return MH;
}

static int Min_HeapExchangeElements(MIN_HEAP *MH, int i, int j)
{
  double tempHeapKey;
  int tempID = 0;
  void *tempData;

  if (i < 0 || i >= MH->CurrHeapSize) {
    fprintf(stderr, "Min_HeapExchangeElements: index %d exceeds size of heap!!\n", i);
    return (ERROR);
  }
  if (j < 0 || j >= MH->CurrHeapSize) {
    fprintf(stderr, "Min_HeapExchangeElements: index %d exceeds size of heap!!\n", j);
    return (ERROR);
  }

  // note that i and j should have parent and child relationship
  if (i == j) {
    fprintf(stderr, "Min_HeapExchangeElements: The indices %d and %d are the same!!\n", i, j);
    return (ERROR);
  }
  else if (i < j) {
    if (j != 2 * i + 1 && j != 2 * i + 2) {
      fprintf(stderr, "Min_HeapExchangeElements: %d is not the child of %d!!\n", j, i);
      return (ERROR);
    }
  }
  else {
    if (i != 2 * j + 1 && i != 2 * j + 2) {
      fprintf(stderr, "Min_HeapExchangeElements: %d is not the child of %d!!\n", i, j);
      return (ERROR);
    }
  }

  // Copy i's stuff over to temporary storage
  tempHeapKey = MH->MHE_array[i].HeapKey;
  tempData = MH->MHE_array[i].Data;
  if (MH->max_id_array_size > 0) tempID = MH->MHE_array[i].id;

  // Copy j's stuff to i
  MH->MHE_array[i].HeapKey = MH->MHE_array[j].HeapKey;
  MH->MHE_array[i].Data = MH->MHE_array[j].Data;

  // Copy i's stuff to j via temporary storage
  MH->MHE_array[j].HeapKey = tempHeapKey;
  MH->MHE_array[j].Data = tempData;

  // Need to update id_array
  if (MH->max_id_array_size > 0) {
    MH->MHE_array[i].id = MH->MHE_array[j].id;
    MH->MHE_array[j].id = tempID;
    MH->id_array[MH->MHE_array[i].id] = i;
    MH->id_array[MH->MHE_array[j].id] = j;
  }
  return (NO_ERROR);
}

int Min_HeapEditKeyIndexID(MIN_HEAP *MH, int id, double newKey)
{
  double tempKey;

  if (id < 0 || id >= MH->max_id_array_size) {
    fprintf(stderr, "Min_HeapEditKeyIndexID: ID %d is out of range. id_array size is %d\n", id, MH->max_id_array_size);
    return (ERROR);
  }

  if (Min_HeapIdIsInHeap(MH, id) == NOT_IN_HEAP) {
    fprintf(stderr, "Min_HeapEditKeyIndexID: ID %d not in heap\n", id);
    return (ERROR);
  }

  tempKey = MH->MHE_array[MH->id_array[id]].HeapKey;
  MH->MHE_array[MH->id_array[id]].HeapKey = newKey;

  if (tempKey > newKey) {
    Min_HeapifyUp(MH, MH->id_array[id]);
  }
  else if (tempKey < newKey) {
    Min_HeapifyDown(MH, MH->id_array[id]);
  }
  else {
    // key same as before, nothing to do
  }
  return (NO_ERROR);
}

int Min_HeapQueryKeyIndexID(MIN_HEAP *MH, int id, double *key)
{
  if (Min_HeapIdIsInHeap(MH, id)) {
    *key = MH->MHE_array[MH->id_array[id]].HeapKey;
    return (NO_ERROR);
  }
  else {
    return (-1);
  }
}

int Min_HeapExtract(MIN_HEAP *MH, double *key, void **data, int *id)
{
  if (MH->CurrHeapSize == 0) {
    fprintf(stderr, "Min_HeapExtract: There's no element to be extracted!!\n");
    return (ERROR);
  }

  // Copy data out...
  *key = MH->MHE_array[0].HeapKey;
  *data = MH->MHE_array[0].Data;
  *id = MH->MHE_array[0].id;

  MH->id_array[*id] = -1;  // after extracting element, the pointer into the MHE_array has to be set to negative
  MH->CurrHeapSize--;

  if (MH->CurrHeapSize != 0) {
    MH->MHE_array[0].HeapKey = MH->MHE_array[MH->CurrHeapSize].HeapKey;
    MH->MHE_array[0].Data = MH->MHE_array[MH->CurrHeapSize].Data;
    MH->MHE_array[0].id = MH->MHE_array[MH->CurrHeapSize].id;
    MH->id_array[MH->MHE_array[MH->CurrHeapSize].id] = 0;
    Min_HeapifyDown(MH, 0);
  }

  return (NO_ERROR);
}

// Min_HeapifyDown assumes index might be bigger than its children but everyone above it is ok
static int Min_HeapifyDown(MIN_HEAP *MH, int index)
{
  int left, right, smallest, curr_index;

  if (index < 0 || index >= MH->CurrHeapSize) {
    fprintf(stderr, "Min_HeapifyDown: index %d exceeds size of heap!!\n", index);
    return (ERROR);
  }

  curr_index = index;

  while (1) {
    left = 2 * curr_index + 1;
    right = 2 * curr_index + 2;

    if (left < MH->CurrHeapSize) {
      if (MH->MHE_array[left].HeapKey < MH->MHE_array[curr_index].HeapKey)
        smallest = left;
      else
        smallest = curr_index;

      if (right < MH->CurrHeapSize && MH->MHE_array[right].HeapKey < MH->MHE_array[smallest].HeapKey) smallest = right;
    }
    else {
      smallest = curr_index;
    }

    if (smallest == curr_index)  // If heap property is not violated, break out of loop
      break;

    Min_HeapExchangeElements(MH, smallest, curr_index);

    curr_index = smallest;
  }

  return (NO_ERROR);
}

// Min_HeapifyUp assumes index's parent might be bigger than it, but everyone belong it is ok.
static int Min_HeapifyUp(MIN_HEAP *MH, int index)
{
  int parent, curr_index;

  if (index < 0 || index >= MH->CurrHeapSize) {
    fprintf(stderr, "Min_HeapifyUp: index %d exceeds size of heap!!\n", index);
    return (ERROR);
  }

  curr_index = index;

  while (1) {
    if (curr_index <= 0)  // If at root, there's no parents...
      break;

    parent = (curr_index - 1) / 2;  // Note that truncation gives you correct parent
    if (MH->MHE_array[parent].HeapKey > MH->MHE_array[curr_index].HeapKey) {
      Min_HeapExchangeElements(MH, parent, curr_index);
      curr_index = parent;
    }
    else {
      break;
    }
  }

  return (NO_ERROR);
}

int Min_HeapIdIsInHeap(MIN_HEAP *MH, int id)
{
  if (id < 0 || id >= MH->max_id_array_size) return (NOT_IN_HEAP);

  if (MH->id_array[id] < 0)
    return (NOT_IN_HEAP);
  else
    return (IS_IN_HEAP);
}

int Min_HeapInsert(MIN_HEAP *MH, double key, void *data, int id)
{
  if (id < 0 || id >= MH->max_id_array_size) {
    fprintf(stderr, "Min_HeapInsert: ID %d is out of range. id_array size is %d\n", id, MH->max_id_array_size);
    return (ERROR);
  }

  if (Min_HeapIdIsInHeap(MH, id) == IS_IN_HEAP) {
    fprintf(stderr, "Min_HeapInsert: ID %d is already found in heap!!!\n", id);
    return (ERROR);
  }

  if (MH->CurrHeapSize == MH->MaxHeapSize) {
    fprintf(
        stderr,
        "Min_HeapInsert: Heap runs out of space. Current version doesn't allocate more space when heap is full!!\n");
    return (ERROR);
  }
  else {
    MH->MHE_array[MH->CurrHeapSize].HeapKey = key;
    MH->MHE_array[MH->CurrHeapSize].Data = data;
    MH->MHE_array[MH->CurrHeapSize].id = id;
    MH->id_array[id] = MH->CurrHeapSize;

    MH->CurrHeapSize++;
    Min_HeapifyUp(MH, MH->CurrHeapSize - 1);
  }

  return (NO_ERROR);
}

int Min_HeapFree(MIN_HEAP *MH)
{
  free(MH->MHE_array);
  if (MH->id_array) free(MH->id_array);
  free(MH);
  return (NO_ERROR);
}

void Min_HeapInternalCheck(MIN_HEAP *MH, int PrintContent)
{
  int i, total = 0;
  int child1, child2;

  fprintf(stderr, "Performing Heap Consistency Check: Current heap size is: %d\n", MH->CurrHeapSize);

  // test id consistency
  for (i = 0; i < MH->max_id_array_size; i++) {
    if (MH->id_array[i] != -1) {
      total++;

      if (MH->id_array[i] < 0 || MH->id_array[i] >= MH->max_id_array_size) {
        fprintf(stderr,
                "Min_HeapInternalCheck: Index %d is out of range. Current heap size is %d\n",
                MH->id_array[i],
                MH->CurrHeapSize);
        fprintf(stderr, "Consistency Check Unsuccessful\n\n");
        return;
      }
      if (MH->MHE_array[MH->id_array[i]].id != i) {
        fprintf(stderr,
                "Min_HeapInternalCheck: i = %d, but MH->MHE_array[MH->id_array[i]].id = %d\n",
                i,
                MH->MHE_array[MH->id_array[i]].id);
        fprintf(stderr, "Consistency Check Unsuccessful\n\n");
        return;
      }
    }
  }

  // test total
  if (total != MH->CurrHeapSize) {
    fprintf(stderr, "Total (%d) in id_array not equal to heap size: %d\n", total, MH->CurrHeapSize);
    fprintf(stderr, "Consistency Check Unsuccessful\n\n");
    return;
  }
  // test heap property
  for (i = 0; i < MH->CurrHeapSize; i++) {
    child1 = 2 * i + 1;
    child2 = 2 * i + 2;

    if (PrintContent) fprintf(stderr, "Parent's key is %f, id %d\n", MH->MHE_array[i].HeapKey, MH->MHE_array[i].id);
    if (child1 < MH->CurrHeapSize) {
      if (MH->MHE_array[child1].HeapKey < MH->MHE_array[i].HeapKey) {
        fprintf(stderr,
                "Child index %d (key %f) < (key %f) Parent %d\n",
                child1,
                MH->MHE_array[child1].HeapKey,
                MH->MHE_array[i].HeapKey,
                i);
        fprintf(stderr, "Consistency Check Unsuccessful\n\n");
        return;
      }
      if (PrintContent)
        fprintf(stderr, "Child1 key is %f, id %d\n", MH->MHE_array[child1].HeapKey, MH->MHE_array[child1].id);
    }
    if (child2 < MH->CurrHeapSize) {
      if (MH->MHE_array[child2].HeapKey < MH->MHE_array[i].HeapKey) {
        fprintf(stderr,
                "Child index %d (key %f) < (key %f) Parent %d\n",
                child2,
                MH->MHE_array[child2].HeapKey,
                MH->MHE_array[i].HeapKey,
                i);
        fprintf(stderr, "Consistency Check Unsuccessful\n\n");
        return;
      }
      if (PrintContent)
        fprintf(stderr, "Child2 key is %f, id %d\n", MH->MHE_array[child2].HeapKey, MH->MHE_array[child2].id);
    }
  }
  fprintf(stderr, "Consistency Check Successful\n\n");
}

int Min_HeapGetCurrSize(MIN_HEAP *MH) { return (MH->CurrHeapSize); }
