#pragma once

#include <omp.h>

#include "timer.h"
#include <stdio.h>


void ROMP_show_stats(FILE*);

// An annotated omp for loop looks like this...
//
#if 0
	#ifdef HAVE_OPENMP
	ROMP_PF_begin
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

#if 1

#define if_ROMP(LEVEL)

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

#define if_ROMP(LEVEL) \
    if (ROMP_pf_stack.staticInfo && \
        (ROMP_if_parallel(ROMP_level_##LEVEL,&ROMP_pf_static))) \
    // end of macro

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


