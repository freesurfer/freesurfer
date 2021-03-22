//
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

#include "base.h"

// Useful functions
//
_Pragma("GCC diagnostic ignored \"-Wunused-function\"")
static bool int_le(int lhs, int rhs) { return lhs <= rhs; }

// Usage
//
#define SORT_NAME           sort_int
#define SORT_NAME_partition sort_int_partition
#define SORT_NAME_isSorted  sort_int_isSorted
#define SORT_NAME_small     sort_int_small
   
#define SORT_ELEMENT    int
#define SORT_LE         int_le
#include "sort_definition.h"
