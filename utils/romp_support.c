#include "romp_support.h"
#include "romp_support.h"
#include <malloc.h>
#include <stdlib.h>

ROMP_level romp_level = 
	ROMP_shown_reproducible;
	// should not be set ro ROMP_serial

typedef struct StaticData {
    ROMP_pf_static_struct* next;
    Nanosecs in_pf, in_pflb;
} StaticData;

ROMP_pf_static_struct* known_ROMP_pf;


static void rompExitHandler(void)
{
  static int count;
  if (count++ > 0) return;
  fprintf(stderr, "ROMP staticExitHandler called\n");
  ROMP_show_stats(stderr);
}

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
	
	static int atExitCount;
	if (atExitCount++ == 0) atexit(rompExitHandler);
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
    if (tid >= 8*sizeof(int)) return;
    int tidMask = 1<<tid;
    if (pf_stack->tids_active & tidMask) {
        fprintf(stderr, "Active tid in ROMP_pflb_begin\n");
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
    Nanosecs delta = TimerElapsedNanosecs(&pflb_stack->beginTime);
    StaticData* staticData = (StaticData*)(pf_static->ptr);
    #pragma omp atomic
    staticData->in_pflb.ns += delta.ns;
    int tid = omp_get_thread_num();
    if (pflb_stack->tid != tid) {
        fprintf(stderr, "Bad tid in ROMP_pflb_end %s:%d\n",
	    pf_static->file, pf_static->line);
        exit(1);
    }
    if (tid >= 8*sizeof(int)) return;
    int tidMask = 1<<tid;
    if (!(pf_stack->tids_active & tidMask)) {
        fprintf(stderr, "Inactive tid in ROMP_pflb_end %s:%d\n",
	    pf_static->file, pf_static->line);
        exit(1);
    }
    #pragma omp atomic
    pf_stack->tids_active ^= tidMask;
}

static const char* mainFile = NULL;
static int mainLine = 0;
static NanosecsTimer mainTimer;
void ROMP_show_stats(FILE* file)
{
    fprintf(file, "ROMP_show_stats\n");
    fprintf(file, "file, line, in pf, in pflb, body/elapsed\n");
    if (mainFile)  {
      Nanosecs mainDuration = TimerElapsedNanosecs(&mainTimer);
      fprintf(file, "%s, %d, %12ld, %12ld, 1.0\n", mainFile, mainLine, mainDuration.ns, mainDuration.ns);
    }
    ROMP_pf_static_struct* pf;
    for (pf = known_ROMP_pf; pf; ) {
    	StaticData* sd = (StaticData*)(pf->ptr);
    	fprintf(file, "%s, %d, %12ld, %12ld, %6.3g\n", pf->file, pf->line, sd->in_pf.ns, sd->in_pflb.ns,
	  (double)sd->in_pflb.ns/(double)sd->in_pf.ns);
    	pf = sd->next;
    }
}

void ROMP_main_started(const char* file, int line) {
    TimerStartNanosecs(&mainTimer);
    mainFile = file;
    mainLine = line;
}

#if 0

// example
//
int main(int argc, char* argv[])
{
    fprintf(stdout, "%s:%d main()\n", __FILE__, __LINE__);
    static const int v_size = 1000;
    int i;
    double sum = 0;

    omp_set_num_threads(1);
    fprintf(stdout, "#threads:%d\n", omp_get_max_threads());

    ROMP_PF_begin
    #pragma omp parallel for if_ROMP(experimental) reduction(+:sum)
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
