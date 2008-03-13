// MARS Project (Thomas Yeo and Mert Sabuncu), MIT, CSAIL (c) 2006-2008

#ifndef min_heap
#define min_heap

#ifndef ERROR
#define NO_ERROR 1
#define ERROR 0
#endif

#define IS_IN_HEAP 1
#define NOT_IN_HEAP 0

//indices start from 0 to size-1
typedef struct {
      double HeapKey;  //The HeapKey defines Heap Order
      void *Data; //Each Heap entry has a corresponding data pointer
      int id; //this is a hack that allows user to find and modify a node of a heap really easily, for example id can be equal to 5 for a particular node. (note that id >= 0)Also, note that id should be unique to each node and is provided by the users. For example, it can correspond to vertex number. Currently, one need to set a maximum possible id size so that we can allocate the table accordingly.
} MIN_HEAP_ELEMENT, MHE;

typedef struct {
      MHE *MHE_array;
      int *id_array; //This allows the user to find and modify a node of a heap really easily. For example, id_array[5] = 10 means MHE_array[10].id = 5. Thus, if a user wants to modify the key of the node of id 5, This can be done through the function: Min_HeapEditKeyIndexID
      int CurrHeapSize;
      int MaxHeapSize;
      int max_id_array_size;
} MIN_HEAP;


MIN_HEAP *Min_HeapAllocate(int max_size, int max_id_array_size);
int Min_HeapEditKeyIndexID(MIN_HEAP *MH, int id, double newKey);
int Min_HeapQueryKeyIndexID(MIN_HEAP *MH, int id, double *key);
int Min_HeapExtract(MIN_HEAP *MH, double *key, void **data, int *id);
int Min_HeapInsert(MIN_HEAP *MH, double key, void *data, int id);
int Min_HeapFree(MIN_HEAP *MH);
int Min_HeapIdIsInHeap(MIN_HEAP *MH, int id);
void Min_HeapInternalCheck(MIN_HEAP *MH, int PrintContent);
int Min_HeapGetCurrSize(MIN_HEAP *MH);

#endif
