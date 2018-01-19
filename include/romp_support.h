/**
 * @file  romp_support.h
 * @brief prototypes and structures for getting reprodiucible results from and for timing omp loops.
 *
 */
/*
 * Original Author: Bevin Brett
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2017/12 $
 *    $Revision: 1.0 $
 *
 * Copyright © 2012 The General Hospital Corporation (Boston, MA) "MGH"
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

#pragma once

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "timer.h"
#include <stdio.h>


void ROMP_show_stats(FILE*);

// An annotated omp for loop looks like this...
//
#if 0
	ROMP_PF_begin
	#ifdef HAVE_OPENMP
	#pragma omp parallel for if_ROMP(experimental) ...
	#endif
  	for (n=0; n < nsegid; n++)
	{
	    ROMP_PFLB_begin
	    ...
	    ROMP_PFLB_continue
	    ...
	    ROMP_PFLB_end
	}
	ROMP_PF_end
    
#endif

// A possibly parallel for loop that has a reproducible set of reductions
// Macros and include files simplify coding it
//
#if 0
                
    See the end of romp_support.c
        
#endif


// Conditionalize a parallel for with
//
typedef enum ROMP_level { 
    ROMP_level_serial,                // always run this code serially
    ROMP_level_experimental,          // hasn't even been tested for correctness yet
    ROMP_level_fast,                  // is known to get differing results parallel and serial
    ROMP_level_assume_reproducible,   // is suspected of getting same results
    ROMP_level_shown_reproducible,    // is tested and shown to get same results 
    ROMP_level__size
    } ROMP_level;
extern ROMP_level romp_level;

// Surround a parallel for

#define ROMP_maxWatchedThreadNum 4

typedef struct ROMP_pf_static_struct { 
    void * volatile ptr; 
    const char*     file; 
    unsigned int    line; 
} ROMP_pf_static_struct;

int ROMP_if_parallel(ROMP_level, ROMP_pf_static_struct*);

 
typedef struct ROMP_pf_stack_struct  { 
    struct ROMP_pf_static_struct * staticInfo; 
    NanosecsTimer beginTime;
    Nanosecs      watchedThreadBeginCPUTimes[ROMP_maxWatchedThreadNum];
    int 	  tids_active;
    int 	  skip_pflb_timing;
    ROMP_level    saved_ROMP_level;
} ROMP_pf_stack_struct;

#define ROMP_main ROMP_main_started(__FILE__, __LINE__);
    
void ROMP_main_started(const char* file, int line);

void ROMP_pf_begin(
    ROMP_pf_static_struct * pf_static,
    ROMP_pf_stack_struct  * pf_stack);

void ROMP_pf_end(
    ROMP_pf_stack_struct  * pf_stack);

// Begin and end of the loop body
// Note: continue does need to be annotated
//       break and return and similar exits are illegal
//
typedef struct ROMP_pflb_stack_struct {
    ROMP_pf_stack_struct * pf_stack;
    NanosecsTimer beginTime;
    int tid;
} ROMP_pflb_stack_struct;

#if 0

#define if_ROMP(LEVEL)

#define if_ROMP2(CONDITION, LEVEL) \
    if ((CONDITION)) \
    // end of macro

#define ROMP_PF_begin \
    {

#define ROMP_PF_end \
    }

#define ROMP_PFLB_begin
#define ROMP_PFLB_end
#define ROMP_PFLB_continue \
    { continue; }
#define ROMP_PF_continue \
    ROMP_PFLB_continue
    
#else

#define if_ROMPLEVEL(LEVEL) \
    if (ROMP_pf_stack.staticInfo && \
        (ROMP_if_parallel(LEVEL,&ROMP_pf_static))) \
    // end of macro

#define if_ROMP2(CONDITION, LEVEL) \
    if ((CONDITION) && \
        ROMP_pf_stack.staticInfo && \
        (ROMP_if_parallel(ROMP_level_##LEVEL,&ROMP_pf_static))) \
    // end of macro

#define if_ROMP(LEVEL) if_ROMPLEVEL(ROMP_level_##LEVEL)

#define ROMP_PF_begin \
    { \
    static ROMP_pf_static_struct ROMP_pf_static = { 0L, __FILE__, __LINE__ }; \
    ROMP_pf_stack_struct  ROMP_pf_stack;  \
    ROMP_pf_begin(&ROMP_pf_static, &ROMP_pf_stack);

#define ROMP_PF_end \
    ROMP_pf_end(&ROMP_pf_stack); \
    }

#define ROMP_PFLB_begin \
    /* ROMP_pflb_stack_struct  ROMP_pflb_stack;  \
    if (!ROMP_pf_stack.skip_pflb_timing) ROMP_pflb_begin(&ROMP_pf_stack, &ROMP_pflb_stack); */ \
    // end of macro

#define ROMP_PFLB_end \
    /* if (!ROMP_pf_stack.skip_pflb_timing) ROMP_pflb_end(&ROMP_pflb_stack); */ \
    // end of macro

#define ROMP_PFLB_continue \
    { /* if (!ROMP_pf_stack.skip_pflb_timing) ROMP_PFLB_end; */ continue; } \
    // end of macro
    
#define ROMP_PF_continue \
    ROMP_PFLB_continue

#endif

void ROMP_pflb_begin(
    ROMP_pf_stack_struct   * pf_stack,
    ROMP_pflb_stack_struct * pflb_stack);

void ROMP_pflb_end(
    ROMP_pflb_stack_struct  * pflb_stack);


// Reproducible reductions
//
typedef struct ROMP_Distributor ROMP_Distributor;

struct ROMP_Distributor {

    #define ROMP_DISTRIBUTOR_PARTIAL_CAPACITY   128
    #define ROMP_DISTRIBUTOR_REDUCTION_CAPACITY   3

    double* originals[ROMP_DISTRIBUTOR_REDUCTION_CAPACITY];

    struct Partials {
        int lo;
        int hi;
        double  partialSum[ROMP_DISTRIBUTOR_REDUCTION_CAPACITY];
    } partials[ROMP_DISTRIBUTOR_PARTIAL_CAPACITY];

    int partialSize;
};

#define ROMP_PARTIALSUM(REDUCTION_INDEX) ROMP_distributor.partials[ROMP_index].partialSum[REDUCTION_INDEX]
 
void ROMP_Distributor_begin(ROMP_Distributor* distributor,
    int lo, int hi, 
    double* sumReducedDouble0, 
    double* sumReducedDouble1, 
    double* sumReducedDouble2); 

double ROMP_Distributor_end(ROMP_Distributor* distributor);
