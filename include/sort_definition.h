/**
 * @brief parallel qsort 
 *      with threads
 *      with comparison being a macro
 *
 */
/*
 * Original Author: Bevin Brett
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

// See sort.h for examples of usage
//
_Pragma("GCC diagnostic ignored \"-Wunused-function\"")
static void SORT_NAME_small(SORT_ELEMENT* elts, size_t size) {
    while (size > 1) {
        SORT_ELEMENT elt_0 = elts[0], least = elt_0;
        size_t least_i = 0, i;
        for (i = 1; i < size; i++) {
            if (SORT_LE(least, elts[i])) continue;
            least   = elts[i];
            least_i = i;
        }
        elts[0]       = least; 
        elts[least_i] = elt_0;
        size--;
        elts++;
    }
}

_Pragma("GCC diagnostic ignored \"-Wunused-function\"")
static size_t SORT_NAME_partition(SORT_ELEMENT* elts, size_t size) 
//
// Returns 0 only when all values are equal, otherwise partitions into two non-empty sets with all the elts
// in the [0..result) set being less than or equal to all the elts in the [result..size) set.
//
{
    if (size < 2) return 0;

    // Choose the splitting value
    //
    SORT_ELEMENT splitter;
    if (size < 100) {
        splitter = elts[size/2];
    } else {
        size_t step = size/8;
        SORT_ELEMENT sample[5];
        sample[0] = elts[step*0];
        sample[1] = elts[step*1];
        sample[2] = elts[step*2];
        sample[3] = elts[step*3];
        sample[4] = elts[step*4];
        SORT_NAME_small(sample, 5);
        splitter = sample[2];
    }

    // Partition
    // [0..lo)    are <= splitter
    // [hi..size) are not
    // in between is unknown
    //
Partition:;
    size_t lo = 0;
    size_t hi = size;
    while (lo < hi) {
        while (            SORT_LE(elts[lo  ], splitter)) { lo++; if (lo == hi) break; }
        while (lo < hi && !SORT_LE(elts[hi-1], splitter)) { hi--;                      }
        // now lo <= hi, [0..lo) <= splitter  [hi..size) > splitter
        if (lo == hi) break;
        // elts[lo] > splitter, elts[hi-1] <= splitter, so swap them and skip over them.
        SORT_ELEMENT temp = elts[lo]; elts[lo] = elts[hi-1]; elts[hi-1] = temp;
        lo++; hi--;
    } 
    
    // The bottom partition must at least contain the splitter.
    // The high partition might be empty, for two reasons
    //      (a) all the values are the same, or
    //      (b) by bad luck, the splitter is the highest value
    // There are two choices
    //      (a) find a lower splitter
    //      (b) filter at least some of the splitter values out of the low partition
    // The following does both
    //
    if (hi == size) {
        while (0 < hi) {
            SORT_ELEMENT temp = elts[hi-1];
            if   (SORT_LE(splitter, temp)) hi--;        // accept some of the values equal to splitter into [hi..size)
            else { splitter = temp; goto Partition; }   // found a lower splitter
        }
    }
    
    // Return the begin of the high partition
    //
    return hi;
}

_Pragma("GCC diagnostic ignored \"-Wunused-function\"")
static void SORT_NAME(SORT_ELEMENT* elts, size_t size, bool useThreads) {

#ifdef HAVE_OPENMP

    // It is only sensible to use threads when there are many elements to sort.
    // This is done by recursively partitioning the data until all the threads are busy. 
    // Currently allow up to 16 partitions
    //
    if (useThreads && size > 10000) {
#define PARTITIONS_CAPACITY 16
        size_t losA[PARTITIONS_CAPACITY], sizesA[PARTITIONS_CAPACITY];
        size_t losB[PARTITIONS_CAPACITY], sizesB[PARTITIONS_CAPACITY];
        size_t* losIn  = losA, *sizesIn  = sizesA;
        size_t* losOut = losB, *sizesOut = sizesB;
        size_t nPartitions = 1;
        losIn  [nPartitions-1] = 0;
        sizesIn[nPartitions-1] = size;
        while (2*nPartitions <= PARTITIONS_CAPACITY) {
            unsigned int p;
            #pragma omp parallel for schedule(guided)
            for (p = 0; p < nPartitions; p++) {
                size_t lo     = losIn[p];
                size_t middle = SORT_NAME_partition(&elts[lo], sizesIn[p]);
                losOut[2*p+0] = lo;        sizesOut[2*p+0] = middle;
                losOut[2*p+1] = lo+middle; sizesOut[2*p+1] = sizesIn[p]-middle;
            }
            // Have doubled the number of partitions
            nPartitions *= 2;
            // Swap the ins and the outs
            size_t* temp; 
            temp = losIn;   losIn   = losOut;   losOut   = temp;
            temp = sizesIn; sizesIn = sizesOut; sizesOut = temp;
        }
        unsigned int p;
        #pragma omp parallel for schedule(guided)
        for (p = 0; p < nPartitions; p++) {
            SORT_NAME(&elts[losIn[p]],sizesIn[p], false);
        }        
        return;
    }
#endif
    
    // It is only sensible to partition when there are several elements to sort
    //
    if (size > 16) {
        size_t middle = SORT_NAME_partition(elts, size);    // returns 0 only when all values are equal
        if (middle > 0) {
            SORT_NAME(elts +      0,        middle, false);
            SORT_NAME(elts + middle, size - middle, false);
        }
        return;
    }
    
    // When there are only a few elements, tight loops are best
    // This code avoids branch mispredicts, and avoids unnecessary mem accesses
    //
    SORT_NAME_small(elts, size);
}

_Pragma("GCC diagnostic ignored \"-Wunused-function\"")
static bool SORT_NAME_isSorted(SORT_ELEMENT* elts, size_t size) {
    for (; size > 1; size--) { 
        if (SORT_LE(elts[size-2], elts[size-1])) continue;
        return false;
    }
    return true;
}

#undef SORT_NAME
#undef SORT_NAME_partition
#undef SORT_NAME_isSorted
#undef SORT_NAME_small

#undef SORT_ELEMENT
#undef SORT_NAME
#undef SORT_LE
