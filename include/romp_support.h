/**
 * @brief prototypes and structures for getting reprodiucible results from and for timing omp loops.
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

#pragma once

// OPTIONS

// Uncomment this to compile in the code that collects the statistics 
//#define ROMP_SUPPORT_ENABLED


// For efficiency reasons, only the first few threads are studied
#define ROMP_maxWatchedThreadNum 4


// This code requires omp to work, but can be compiled without it
#ifdef HAVE_OPENMP
  #include <omp.h>
  int romp_omp_get_thread_num();
  #define omp_get_thread_num romp_omp_get_thread_num
#else
  // Make it easier to write code which is insensitive to OpenMP being present
  static inline int omp_get_max_threads() { return 1; }
  static inline int omp_get_thread_num() { return 0; }
  #define omp_set_num_threads(n)
#endif


// The source of the times
//
#include "timer.h"


// Optionally tell romp when the main program has started
// so it can time the period before the first parallel loop    
//
#define ROMP_main ROMP_main_started(__BASE_FILE__, __LINE__);
void ROMP_main_started(const char* file, int line);


// Write any collected stats to the file
//
void ROMP_show_stats(FILE* file);


// Return the number of times the code has gone from serial to parallel
// 
size_t ROMP_countGoParallel();
    // Useful during debugging to write conditional code looking for a problem 
    // that is being caused by parallelism

// omp for loops should be annotated with the following macros
// so that the romp support knows of their existence.  Other omp loops hjave no data collection.
//
// An annotated omp for loop that does not to reductions looks like this...
//
#if 0
	ROMP_PF_begin	    	    	    	    	    // times the whole loop
	#ifdef HAVE_OPENMP
	#pragma omp parallel for if_ROMP(experimental) ...  // decides whether it should be parallel
	#endif
  	for (n=0; n < nsegid; n++)
	{
	    ROMP_PFLB_begin 	    	    	    	    // expensive timing of each loop body
	    ...     	    	    	    	    	    // so these are usually defined as no code
	    ROMP_PFLB_continue	    	    	    	    // Needed since doesn't go thru the end
	    ...
	    ROMP_PFLB_end
	}
	ROMP_PF_end
#endif

// A parallel for loop that has a reproducible set of floating point reductions
// is shown at the end of romp_support.c
        


// The omp for loops should be conditionalized to only be parallel if the romp_level allows it
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


// To speed up collection, an annotated omp for loop
// has a static structure holding information about the loop
// and one on the stack holding current timing information
//
typedef struct ROMP_pf_static_struct { 
    void * volatile ptr; 
    const char*     file; 
    const char*     func; 
    unsigned int    line; 
} ROMP_pf_static_struct;

typedef struct ROMP_pf_stack_struct  { 
    struct ROMP_pf_static_struct * staticInfo; 
    Timer timer;
    long      watchedThreadBeginCPUTimes[ROMP_maxWatchedThreadNum];
    int 	  gone_parallel;
    ROMP_level    entry_level;
} ROMP_pf_stack_struct;


void ROMP_pf_begin(
    ROMP_pf_static_struct * pf_static,
    ROMP_pf_stack_struct  * pf_stack);

int ROMP_if_parallel1(ROMP_level);
int ROMP_if_parallel2(ROMP_level, ROMP_pf_stack_struct*);
    // it is the call of this function that tells the ROMP code what level the loop is
    // so there is more done here than just supplying the condition value to the omp if clause
 
void ROMP_pf_end(
    ROMP_pf_stack_struct  * pf_stack);


// Instrumentation of  the Begin and end of the loop body
// Because many loop bodies are too brief to hide the cost of doing these,
// the macros below that add these are usually defined as noops rather than adding them.
// 
// Note: continue is the same as ROMP_pf_end, hence use it
//       break and return and similar exits are illegal
//
typedef struct ROMP_pflb_stack_struct {
    ROMP_pf_stack_struct * pf_stack;
    Timer timer;
    int tid;
} ROMP_pflb_stack_struct;

void ROMP_pflb_begin(
    ROMP_pf_stack_struct   * pf_stack,
    ROMP_pflb_stack_struct * pflb_stack);

void ROMP_pflb_end(
    ROMP_pflb_stack_struct  * pflb_stack);


// The conditionalized macros that either do or don't add the variables and calls based on the above
//
#if !defined(ROMP_SUPPORT_ENABLED)

    #define if_ROMPLEVEL(LEVEL) \
	if (ROMP_if_parallel1(LEVEL)) \
	// end of macro

    #define if_ROMP2(CONDITION, LEVEL) \
	if ((CONDITION) && ROMP_if_parallel1(ROMP_level_##LEVEL)) \
	// end of macro

    #define ROMP_PF_begin \
	{

    #define ROMP_PF_end \
	}

    #define ROMP_PFLB_begin
    #define ROMP_PFLB_end
    #define ROMP_PFLB_continue \
	{ continue; }
	
#else

    #define if_ROMPLEVEL(LEVEL) \
	if (ROMP_pf_stack.staticInfo && \
            (ROMP_if_parallel2(LEVEL,&ROMP_pf_stack))) \
	// end of macro

    #define if_ROMP2(CONDITION, LEVEL) \
	if ((CONDITION) && \
            ROMP_pf_stack.staticInfo && \
            (ROMP_if_parallel2(ROMP_level_##LEVEL,&ROMP_pf_stack))) \
	// end of macro

    #define ROMP_PF_begin \
	{ \
	static ROMP_pf_static_struct ROMP_pf_static = { 0L, __BASE_FILE__, __func__, __LINE__ }; \
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
    
#endif


#define if_ROMP(LEVEL) if_ROMPLEVEL(ROMP_level_##LEVEL)

#define ROMP_PF_continue \
    ROMP_PFLB_continue


// When measuring where time is going, 
// it is convenient to be able to do some scopes other than parallel loops
// Such scopes look a single trip loop, so can exploit the above support
// 
#define ROMP_SCOPE_begin    ROMP_PF_begin
#define ROMP_SCOPE_end	    ROMP_PF_end


// Deterministic reductions
//
// omp reduction clauses do the reduction operations in arbitrary order
// but floating point add and mul are not precisely commutative, and are compiler and thead-count specific
// so they are causing small initial variations that can cascade into test system failures.
// The following implements are deterministic ordering of the operations independent of thread-count.
//
// For usage, see the example in romp_support.c in the test case at the end.
//
// Only sum is currently supported.
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

void ROMP_Distributor_end(ROMP_Distributor* distributor);
