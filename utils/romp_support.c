#include "romp_support.h"
#include "romp_support.h"
#include <malloc.h>


ROMP_level romp_level = ROMP_shown_reproducible;


typedef struct StaticData {
    ROMP_pf_static_struct* next;
    Nanosecs in_pf, in_pflb;
} StaticData;

ROMP_pf_static_struct* known_ROMP_pf;


static StaticData* initStaticData(ROMP_pf_static_struct * pf_static)
{
    StaticData* ptr = (StaticData*)pf_static->ptr;
    if (ptr) return ptr;
    #pragma omp critical
    {	// Might have been made by another thread
    	ptr = (StaticData*)pf_static->ptr;
    	if (!ptr) {
    	    ptr = (StaticData*)calloc(1, sizeof(StaticData));
    	    pf_static->ptr = ptr;
    	    ptr->next = known_ROMP_pf;
    	    known_ROMP_pf = pf_static;
    	}
    }
    return ptr;
}

void ROMP_pf_begin(
    ROMP_pf_static_struct * pf_static,
    ROMP_pf_stack_struct  * pf_stack) 
{
    pf_stack->staticInfo = pf_static;
    StaticData* staticData = initStaticData(pf_static);
    if (!staticData) return;
    TimerStartNanosecs(&pf_stack->beginTime);
    pf_stack->tids_active = 0;
}

void ROMP_pf_end(
    ROMP_pf_stack_struct  * pf_stack) 
{
    ROMP_pf_static_struct * pf_static = pf_stack->staticInfo;
    if (!pf_static) return;
    Nanosecs delta = TimerElapsedNanosecs(&pf_stack->beginTime);
    StaticData* staticData = (StaticData*)(pf_static->ptr);
    #pragma omp atomic
    staticData->in_pf.ns += delta.ns;
}

void ROMP_pflb_begin(
    ROMP_pf_stack_struct   * pf_stack,
    ROMP_pflb_stack_struct * pflb_stack)
{
    pflb_stack->pf_stack = pf_stack;
    ROMP_pf_static_struct * pf_static = pf_stack->staticInfo;
    if (!pf_static) { pflb_stack->pf_stack = NULL; return; }
    TimerStartNanosecs(&pflb_stack->beginTime);
    int tid = omp_get_thread_num();
    pflb_stack->tid = tid;
    if (tid >= 8*sizeof(int) return;
    int tidMask = 1<<tid;
    if (pf_stack->tids_active & tidMask) {
        fprintf(stderr, "Acitive tid in ROMP_pflb_begin\n");
        exit(1);
    }
    #pragma omp atomic
    pf_stack->tids_active ^= tidMask;
}

void ROMP_pflb_end(
    ROMP_pflb_stack_struct  * pflb_stack)
{
    if (!pflb_stack->pf_stack) return;
    ROMP_pf_stack_struct  * pf_stack = pflb_stack->pf_stack;
    ROMP_pf_static_struct * pf_static = pf_stack->staticInfo;
    Nanosecs delta = TimerElapsedNanosecs(&pf_stack->beginTime);
    StaticData* staticData = (StaticData*)(pf_static->ptr);
    #pragma omp atomic
    staticData->in_pflb.ns += delta.ns;
    int tid = omp_get_thread_num();
    if (pflb_stack->tid != tid) {
        fprintf(stderr, "Bad tid in ROMP_pflb_end %s:%d\n",
	    pf_static->file, pf_static->line);
        exit(1);
    }
    if (tid >= 8*sizeof(int) return;
    int tidMask = 1<<tid;
    if (!(pf_stack->tids_active & tidMask)) {
        fprintf(stderr, "Inactive tid in ROMP_pflb_end %s:%d\n",
	    pf_static->file, pf_static->line);
        exit(1);
    }
    #pragma omp atomic
    pf_stack->tids_active ^= tidMask;
    
}

void ROMP_show_stats(FILE* file)
{
    fprintf(file, "ROMP_show_stats\n");
    ROMP_pf_static_struct* pf;
    for (pf = known_ROMP_pf; pf; ) {
    	StaticData* sd = (StaticData*)(pf->ptr);
    	fprintf(file, "%s:%d %ld %ld\n", pf->file, pf->line, sd->in_pf.ns, sd->in_pflb.ns);
    	pf = sd->next;
    }
}



#if 0

// example
//
int main(int argc, char* argv[])
{
    fprintf(stdout, "%s:%d main()\n", __FILE__, __LINE__);
    static const int v_size = 1000;
    int* v = (int*)malloc(sizeof(int)*v_size);
    int i;
    double sum = 0;

    omp_set_num_threads(1);
    fprintf(stdout, "#threads:%d\n", omp_get_max_threads());

    ROMP_PF_begin
    #pragma omp parallel for if_ROMP(experimental) reduction(+:sum) if_ROMP(shown_reproducible)
    for (i = 0; i < v_size; i++) {
    	ROMP_PFLB_begin
    	sum += 1.0 / i;
    	ROMP_PFLB_end;
    }
    ROMP_PF_end

    ROMP_show_stats(stdout);
    return 0;
}

#endif
