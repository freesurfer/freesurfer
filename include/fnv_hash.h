/**
 * @brief A hash function that can be used to replace getting a sequence of random numbers
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

// See https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
//
static unsigned long fnv_init() { return 0x811c9dc5; }

static unsigned long fnv_add(unsigned long hash, const unsigned char* p, unsigned long size) {
    unsigned long c;
    for (c = 0; c < size; c++) {
        hash ^= p[c];
        hash *= 16777619L;
    }
    return hash;
}

#pragma GCC diagnostic ignored "-Wunused-function"
static float fnv_hash(int i, int k, int* random_counter, double LO, double HI) {

    unsigned int combined[3]; combined[0] = i; combined[1] = k; combined[2] = (*random_counter)++; 

    unsigned long hash = fnv_add(fnv_init(), (unsigned char*)&combined[0], sizeof(combined));

    double f = (double)(hash & (unsigned long)0xffffffff) 
             / (double)        (unsigned long)0xffffffff;

    return (float)((HI-LO)*f + LO);
}
