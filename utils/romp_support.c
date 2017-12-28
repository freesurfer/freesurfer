#include "romp_support.h"
#include "romp_support.h"
#include <malloc.h>
#include <stdlib.h>
#include <pthread.h>

ROMP_level romp_level = 
    ROMP_level_fast;              	// should not be set to ROMP_serial


typedef struct StaticData {
    ROMP_pf_static_struct* next;
    int skip_pflb_timing;
    Nanosecs in_pf, in_pfThread, in_pflb;
    ROMP_level level;
} StaticData;

ROMP_pf_static_struct* known_ROMP_pf;

int ROMP_if_parallel(ROMP_level level, ROMP_pf_static_struct* pf_static) 
{
    StaticData* ptr = (StaticData*)pf_static->ptr;
    if (ptr) ptr->level = level;
    
    ROMP_level current_level = romp_level;
    int result = level >= current_level;

    if (result) romp_level = ROMP_level__size;	// disable nested parallelism
    
    return result;
}

static const char* mainFile = NULL;
static int         mainLine = 0;

static const char* getMainFile() 
{
    if (!mainFile) {    
        static char commBuffer[1024];
        FILE* commFile = fopen("/proc/self/comm", "r");
        int commSize = 0;
        if (commFile) {
            commSize = fread(commBuffer, 1, 1023, commFile);
	    if (commSize > 0) commSize-=1; // drop the \n
	    int i = 0;
	    for (i = 0; i < commSize; i++) {
	        if (commBuffer[i] == '/') commBuffer[i] = '@';
	    }
            fclose(commFile);
        }
        commBuffer[commSize] = 0;
        if (commSize) mainFile = commBuffer;
    }
    return mainFile;
}

static void rompExitHandler(void)
{
    static int once;
    if (once++ > 0) return;
    fprintf(stderr, "ROMP staticExitHandler called\n");
     
    ROMP_show_stats(stderr);
  
    if (getMainFile()) {
        char ROMP_statsFileName[1024];
        ROMP_statsFileName[0] = 0;
        snprintf(ROMP_statsFileName, 1024, "/tmp/ROMP_statsFiles/%s.csv", mainFile);
        FILE* comFile = fopen(ROMP_statsFileName, "a");
        if (!comFile) {
            fprintf(stderr, "Could not create %s\n", ROMP_statsFileName);
        } else {
            ROMP_show_stats(comFile);
            fclose(comFile);
        }
    }
}

static NanosecsTimer mainTimer;

static void initMainTimer() {
    static int once;
    if (once++ == 0) {
        TimerStartNanosecs(&mainTimer);
	atexit(rompExitHandler);
    }
}

void ROMP_main_started(const char* file, int line) {
    initMainTimer();
    mainFile = file;
    mainLine = line;
}

static StaticData* initStaticData(ROMP_pf_static_struct * pf_static)
{
    StaticData* ptr = (StaticData*)pf_static->ptr;
    if (ptr) return ptr;
    #pragma omp critical
    {	// Might have been made by another thread
    	ptr = (StaticData*)pf_static->ptr;
    	if (!ptr) {
	    initMainTimer();
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

    if (romp_level == ROMP_level__size) {
    	// ignore nested parallel so the loop body times are not added to more than one scope
	pf_stack->staticInfo = NULL;
	return;
    }
    
    int i;
    for (i = 0; i < ROMP_maxWatchedThreadNum; i++) {
        pf_stack->watchedThreadBeginCPUTimes[i].ns = 0;
    }
    #pragma omp parallel for schedule(static,1)
    for (i = 0; i < ROMP_maxWatchedThreadNum; i++) {
        int tid = omp_get_thread_num();
	if (tid >= ROMP_maxWatchedThreadNum) continue;
	clockid_t clockid;
	int s = pthread_getcpuclockid(pthread_self(), &clockid);
	if (s != 0) {
	    fprintf(stderr, "%s:%d pthread_getcpuclockid failed", __FILE__, __LINE__);
	    exit(1);
	}
	struct timespec timespec;
        if (clock_gettime(clockid, &timespec) != 0) {
	    fprintf(stderr, "%s:%d clock_gettime failed", __FILE__, __LINE__);
	    exit(1);
	}
	pf_stack->watchedThreadBeginCPUTimes[tid].ns =
    	    timespec.tv_sec * 1000000000L + timespec.tv_nsec;
    }
    
    TimerStartNanosecs(&pf_stack->beginTime);
    pf_stack->skip_pflb_timing = staticData->skip_pflb_timing;
    pf_stack->tids_active      = 0;
    pf_stack->saved_ROMP_level = romp_level;
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

    int i;
    #pragma omp parallel for schedule(static,1)
    for (i = 0; i < ROMP_maxWatchedThreadNum; i++) {
        int tid = omp_get_thread_num();
	if (tid >= ROMP_maxWatchedThreadNum) continue;
	
	Nanosecs* startCPUTime = &pf_stack->watchedThreadBeginCPUTimes[tid];
	if (startCPUTime->ns == 0) {
	    fprintf(stderr, "%s:%d no start time for %d\n", __FILE__, __LINE__, tid);
	    continue;
	}
	if (startCPUTime->ns == -1) {
	    continue;
	}
	clockid_t clockid;
	int s = pthread_getcpuclockid(pthread_self(), &clockid);
	if (s != 0) {
	    fprintf(stderr, "%s:%d pthread_getcpuclockid failed", __FILE__, __LINE__);
	    exit(1);
	}
	struct timespec timespec;
        if (clock_gettime(clockid, &timespec) != 0) {
	    fprintf(stderr, "%s:%d clock_gettime failed", __FILE__, __LINE__);
	    exit(1);
	}
	Nanosecs threadCpuTime;
	threadCpuTime.ns = timespec.tv_sec * 1000000000L + timespec.tv_nsec - startCPUTime->ns;
	startCPUTime->ns = -1;
    
    	#pragma omp atomic
    	staticData->in_pfThread.ns += threadCpuTime.ns;
    }
    
    for (i = 0; i < ROMP_maxWatchedThreadNum; i++) {
        pf_stack->watchedThreadBeginCPUTimes[i].ns = 0;
    }

    if (pf_stack->saved_ROMP_level != ROMP_level__size) {
    	// exiting a non-nested parallel for
        romp_level = pf_stack->saved_ROMP_level;
    }
    
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
        fprintf(stderr, "Active tid in ROMP_pflb_begin %s:%d\n",
	    pf_static->file, pf_static->line);
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

    if (delta.ns < 1000) {	// loop body is too small to time this way...
      staticData->skip_pflb_timing = pf_stack->skip_pflb_timing = 1;
      return;
    }
}

void ROMP_show_stats(FILE* file)
{
    fprintf(file, "ROMP_show_stats\n");
    fprintf(file, "file, line, level, in pf, in pflb, in_pfThread, pflb/elapsed, pft/elapsed\n");

    if (getMainFile())  {
      	Nanosecs mainDuration = TimerElapsedNanosecs(&mainTimer);
      	fprintf(file, "%s, %d, 0, %12ld, %12ld, %12ld, %6.3g, %6.3g\n", 
      	    mainFile, mainLine, 
	    mainDuration.ns, 0L, 0L,
	    1.0,
	    1.0);
    }
    ROMP_pf_static_struct* pf;
    for (pf = known_ROMP_pf; pf; ) {
    	StaticData* sd = (StaticData*)(pf->ptr);
	if (sd->skip_pflb_timing) sd->in_pflb.ns = 0;
    	fprintf(file, "%s, %d, %d, %12ld, %12ld, %12ld, %6.3g, %6.3g\n", 
	    pf->file, pf->line,
	    sd->level,
	    sd->in_pf.ns, sd->in_pflb.ns, sd->in_pfThread.ns,
	    (double)sd->in_pflb.ns    /(double)sd->in_pf.ns,
	    (double)sd->in_pfThread.ns/(double)sd->in_pf.ns);
    	pf = sd->next;
    }
    fprintf(file, "ROMP_show_stats end\n");
}

#if 0

// example
//
int main(int argc, char* argv[])
{
    FILE* comFile = fopen("/proc/self/comm", "r");
    char comBuffer[1024];
    int comSize = 0;
    if (comFile) {
        comSize = fread(comBuffer, 1, 1023, comFile);
        fclose(comFile);
    }
    comBuffer[comSize] = 0;
    fprintf(stdout, "%s:%d main() of %s\n", __FILE__, __LINE__, comBuffer);
    
    static const int v_size = 30000;
    int i;
    double sum = 0;

    fprintf(stdout, "#threads:%d\n", omp_get_max_threads());

    long threadMask = 0;
    ROMP_PF_begin
    #pragma omp parallel for if_ROMP(fast) reduction(+:sum)
    for (i = 0; i < v_size; i++) {
    	ROMP_PFLB_begin
	
	threadMask |= 1 << omp_get_thread_num();
	
    	sum += 1.0 / i;
	
	int j;
	ROMP_PF_begin
    	#pragma omp parallel for if_ROMP(fast) reduction(+:sum)
    	for (j = 0; j < i; j++) {
    	    //ROMP_PFLB_begin
	    
	    sum += 1.0 / j;
	    
    	    ROMP_PFLB_end;
        }
        ROMP_PF_end
	
    	ROMP_PFLB_end;
    }
    ROMP_PF_end

    fprintf(stdout, "ThreadMask:%p\n", (void*)threadMask);
    
    ROMP_show_stats(stdout);
    return 0;
}

#endif
