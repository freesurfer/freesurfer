#include "sort.h"

// Self-test by
//                  gcc -fopenmp -I../include -DSORT_TEST sort.c
//                  ./a.out
//
#ifdef SORT_TEST
int main() {
    size_t capacity = 100000000;
    int*   elts     = (int*)malloc(capacity*sizeof(int));
    size_t size;
    for (size = 0; size < capacity; size = size*3 + 7) {
        size_t i;
        for (i = 0; i < size; i++) elts[i] = i*(size^0x55555) ^ 0x13579;
        sort_int(elts,size,true);
        cheapAssert(sort_int_isSorted(elts, size));
    }
    return 0;
}
#endif
