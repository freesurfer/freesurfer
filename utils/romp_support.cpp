static const int debug = 0;

/**
 * @brief prototypes and structures for getting reproducible results from and for timing omp loops.
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
#include "romp_support.h"
#include "base.h"

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <stdlib.h>
#include <pthread.h>

static void __attribute__((constructor)) before_main() 
{
    int n = omp_get_max_threads();
    if (n <= _MAX_FS_THREADS) return;
    omp_set_num_threads(_MAX_FS_THREADS);
}

int romp_omp_get_thread_num()  { 
#undef omp_get_thread_num
    int n = omp_get_thread_num();
#define omp_get_thread_num romp_omp_get_thread_num
    if (n >= _MAX_FS_THREADS) {
        fprintf(stderr, "Freesurfer supports a maximum of %d threads, set OMP_NUM_THREADS accordingly\n", _MAX_FS_THREADS);
        exit(1);
    }
    return n;
}

// measurement is showing there is only a few percent difference in timing
// between fast and assume_reproducible
// so make the later the default 
//
ROMP_level romp_level =                 // should not be set to ROMP_serial 
    //ROMP_level_fast;
    ROMP_level_assume_reproducible;     // doesn't require the random seed to be     //ROMP_level_shown_reproducible;

typedef struct PerThreadScopeTreeData {
    ROMP_pf_static_struct*         key;
    struct PerThreadScopeTreeData* parent;
    struct PerThreadScopeTreeData* next_sibling;
    struct PerThreadScopeTreeData* first_child;
    long in_scope;
    long in_child_threads[ROMP_maxWatchedThreadNum];    
        // threads may be running in a different scope tree!
} PerThreadScopeTreeData;

static PerThreadScopeTreeData  scopeTreeRoots[ROMP_maxWatchedThreadNum];
static PerThreadScopeTreeData* scopeTreeToS  [ROMP_maxWatchedThreadNum];

static struct { Timer timer; bool inited = false; } tidStartTime[ROMP_maxWatchedThreadNum];
static void maybeInitTidStartTime(int tid) {
    if (tidStartTime[tid].inited) return;
    tidStartTime[tid].timer.reset();
    tidStartTime[tid].inited = true;
}

static PerThreadScopeTreeData* enterScope(PerThreadScopeTreeData* parent, ROMP_pf_static_struct* key) {
    PerThreadScopeTreeData** prev = &parent->first_child;
    while (*prev && (*prev)->key != key) { prev = &(*prev)->next_sibling; } // might need speeding up with a hash table
    if (!*prev) {
        PerThreadScopeTreeData* ptr = (PerThreadScopeTreeData *)calloc(1, sizeof(PerThreadScopeTreeData));
        ptr->key    = key;
        ptr->parent = parent;
        *prev = ptr;
    }
    return *prev;
}

typedef struct StaticData {
    ROMP_pf_static_struct* next;
    int skip_pflb_timing;
    long in_pf_sum, in_pfThread_sum, in_pflb_sum;
    ROMP_level level;
} StaticData;

ROMP_pf_static_struct* known_ROMP_pf;

static long cpuTimeUsed() {
#ifdef __APPLE__
    // not yet supported on mac
    return 0;
#else
    clockid_t clockid = clockid_t();
    if (pthread_getcpuclockid(pthread_self(), &clockid) != 0) {
        fprintf(stderr, "%s:%d pthread_getcpuclockid failed", __FILE__, __LINE__);
        exit(1);
    }
    struct timespec timespec;
    if (clock_gettime(clockid, &timespec) != 0) {
        fprintf(stderr, "%s:%d gettime failed", __FILE__, __LINE__);
        exit(1);
    }
    return timespec.tv_sec * 1000000000L + timespec.tv_nsec;
#endif
}

int ROMP_if_parallel1(ROMP_level level)
{
    return (level >= romp_level);               // sadly this allows nested parallelism
}                                               // sad only because it hasn't been analyzed

static size_t countGoParallel;
size_t ROMP_countGoParallel() { return countGoParallel; }

int ROMP_if_parallel2(ROMP_level level, ROMP_pf_stack_struct* pf_stack) 
{
    ROMP_pf_static_struct * pf_static = pf_stack->staticInfo;
    StaticData* ptr = (StaticData*)pf_static->ptr;
    if (ptr) ptr->level = level;

    ROMP_level current_level = romp_level;
    int result = level >= current_level;

    if (result) {
        countGoParallel++;
        if (debug) 
            fprintf(stderr, "ROMP_if_parallel2 tid:%d pf_stack:%p gone parallel\n",
                omp_get_thread_num(), pf_stack);
        pf_stack->gone_parallel = 1;
        romp_level = ROMP_level__size;  // disable nested parallelism because hard to analyze
    }
    
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
#if defined(ROMP_SUPPORT_ENABLED)
    static int once;
    if (once++ > 0) return;
    if (debug) fprintf(stderr, "ROMP staticExitHandler called\n");
     
    int tid;
    for (tid = 0; tid < ROMP_maxWatchedThreadNum; tid++) {
        PerThreadScopeTreeData* root = &scopeTreeRoots[tid];
        if (tidStartTime[tid].inited) root->in_scope = tidStartTime[tid].timer.nanoseconds();
    }
    
    ROMP_show_stats(stderr);
  
    if (getMainFile()) {
        char ROMP_statsFileName[1024];
        ROMP_statsFileName[0] = 0;
        snprintf(ROMP_statsFileName, 1024, "/tmp/ROMP_statsFiles/%s.csv", mainFile);
        FILE* comFile = fopen(ROMP_statsFileName, "a");
        if (!comFile) {
            fprintf(stderr, "Could not create %s\n", ROMP_statsFileName);
        } else {
            fprintf(stderr, "Created %s\n", ROMP_statsFileName);
            ROMP_show_stats(comFile);
            fclose(comFile);
        }
    }
#endif
}

static Timer mainTimer;

static void initMainTimer() {
    static int once;
    if (once++ == 0) {
        mainTimer.reset();
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
    
    // Don't use omp critical because caller might be within a critical already
    //
#ifdef HAVE_OPENMP
    static omp_lock_t lock; static volatile int lock_inited;
    if (!lock_inited) {
        #pragma omp critical
        if (!lock_inited) {
            omp_init_lock(&lock); 
            lock_inited = 1; 
        }
    }
    omp_set_lock(&lock);
#endif
    {   // Might have been made by another thread
        ptr = (StaticData*)pf_static->ptr;
        if (!ptr) {
            initMainTimer();
            ptr = (StaticData*)calloc(1, sizeof(StaticData));
            pf_static->ptr = ptr;
            ptr->next = known_ROMP_pf;
            known_ROMP_pf = pf_static;
        }   
    }
#ifdef HAVE_OPENMP
    omp_unset_lock(&lock);
#endif
    return ptr;
}

void ROMP_pf_begin(
    ROMP_pf_static_struct * pf_static,
    ROMP_pf_stack_struct  * pf_stack) 
{
    pf_stack->staticInfo = pf_static;
    StaticData* staticData = initStaticData(pf_static);
    if (!staticData) return;

    pf_stack->gone_parallel = 0;
    pf_stack->entry_level = romp_level;

    int tid = 
#ifdef HAVE_OPENMP
    omp_get_thread_num();
#else
    0;
#endif
    if (tid >= ROMP_maxWatchedThreadNum) return;
    
    PerThreadScopeTreeData* tos = scopeTreeToS[tid];
    if (!tos) { tos = &scopeTreeRoots[tid]; scopeTreeToS[tid] = tos; maybeInitTidStartTime(tid); }
    scopeTreeToS[tid] = enterScope(tos, pf_static);
    
    int i;
    for (i = 0; i < ROMP_maxWatchedThreadNum; i++) {
        pf_stack->watchedThreadBeginCPUTimes[i] = 0;
    }

    if (romp_level >= ROMP_level__size) {
        if (debug) fprintf(stderr, "%s:%d only getting child tid start times for tid:%d pf_stack:%p\n", __FILE__, __LINE__, tid, pf_stack);
        pf_stack->watchedThreadBeginCPUTimes[tid] = cpuTimeUsed();
    } else {
        if (debug) fprintf(stderr, "%s:%d get child tid start times for tid:%d pf_stack:%p\n", __FILE__, __LINE__, tid, pf_stack);

#ifdef HAVE_OPENMP
        #pragma omp parallel for schedule(static,1)
#endif
        for (i = 0; i < ROMP_maxWatchedThreadNum; i++) {
            int childTid = 
#ifdef HAVE_OPENMP
            omp_get_thread_num();
#else
            0;
#endif
            if (childTid >= ROMP_maxWatchedThreadNum) continue;
        pf_stack->watchedThreadBeginCPUTimes[childTid] = cpuTimeUsed();
            if (debug) {
#ifdef HAVE_OPENMP
                #pragma omp critical
#endif
                fprintf(stderr, "%s:%d start time gotten for tid:%d childTid:%d\n", __FILE__, __LINE__, tid, childTid);
            }
        }
    
        for (i = 0; i < ROMP_maxWatchedThreadNum; i++) {
            if (pf_stack->watchedThreadBeginCPUTimes[i] == 0) {
                if (debug) fprintf(stderr, "%s:%d no start time for %d\n", __FILE__, __LINE__, i);
            }
        }
    }
    
    pf_stack->timer.reset();
}


static void pf_end_one_thread(int tid, ROMP_pf_stack_struct * pf_stack, PerThreadScopeTreeData* tos)
{
    int childTid = 
#ifdef HAVE_OPENMP
    omp_get_thread_num();
#else
    0;
#endif
    if (childTid >= ROMP_maxWatchedThreadNum) return;

    long &startCPUTime = pf_stack->watchedThreadBeginCPUTimes[childTid];
    if (startCPUTime == -1) {
        return;
    }
    if (startCPUTime == 0) {
        if (debug) fprintf(stderr, "%s:%d no corresponding start time for tid:%d pf_stack:%p\n", __FILE__, __LINE__, tid, pf_stack);
        return;
    }
        
    long threadCpuTime;
    threadCpuTime = cpuTimeUsed() - startCPUTime;
    startCPUTime = -1;  // copes when another iteration of this loop is executed on the same thread
    
    if (debug) fprintf(stderr, "Adding to tid:%d childTid:%d ns:%ld\n", tid, childTid, threadCpuTime);
    tos->in_child_threads[childTid] += threadCpuTime;
}


void ROMP_pf_end(
    ROMP_pf_stack_struct  * pf_stack) 
{
    ROMP_pf_static_struct * pf_static = pf_stack->staticInfo;
    if (!pf_static) return;

    int tid = 
#ifdef HAVE_OPENMP
    omp_get_thread_num();
#else
    0;
#endif
    if (tid < ROMP_maxWatchedThreadNum) {

        PerThreadScopeTreeData* tos = scopeTreeToS[tid];
        if (tos) {

            long delta = pf_stack->timer.nanoseconds();
            tos->in_scope += delta;

            if (pf_stack->gone_parallel)
                if (debug) fprintf(stderr, "ROMP_pf_end tid:%d pf_stack:%p getting other thread times\n",
                    omp_get_thread_num(), pf_stack);

            if (!pf_stack->gone_parallel) {
                // Not the same as putting a condition on the else loop
                // The later is implemented by openmp choosing any available thread to execute the loop body, not the current thread
                pf_end_one_thread(tid, pf_stack, tos);
            } else {
                int i;
                #ifdef HAVE_OPENMP
                    #pragma omp parallel for schedule(static,1)
                #endif
                for (i = 0; i < ROMP_maxWatchedThreadNum; i++) pf_end_one_thread(tid, pf_stack, tos);
            }

            if (romp_level < ROMP_level__size) {
                int i;
                for (i = 0; i < ROMP_maxWatchedThreadNum; i++) {
                    long startCPUTime = pf_stack->watchedThreadBeginCPUTimes[i];
                    if (startCPUTime != 0 && startCPUTime != -1) {
                        if (debug) fprintf(stderr, "%s:%d no end time for %d\n", __FILE__, __LINE__, tid);
                    }
                } 
            }

            scopeTreeToS[tid] = tos->parent;
        }
    }

    romp_level = pf_stack->entry_level;
}


void ROMP_pflb_begin(
    ROMP_pf_stack_struct   * pf_stack,
    ROMP_pflb_stack_struct * pflb_stack)
{
}

void ROMP_pflb_end(
    ROMP_pflb_stack_struct  * pflb_stack)
{
}


static void node_show_stats(FILE* file, PerThreadScopeTreeData* node, unsigned int depth) {
    ROMP_pf_static_struct* pf = node->key;
    StaticData* sd = pf ? (StaticData*)(pf->ptr) : (StaticData*)(NULL);
    
    long inAllThreads; inAllThreads = 0;
    int tid;
    for (tid = 0; tid < ROMP_maxWatchedThreadNum; tid++) {
        inAllThreads += node->in_child_threads[tid];
    }
    
    {
        for (unsigned int d = 0; d < depth; d++) fprintf(file, "    ");
    }
    
    fprintf(file, "%s:%s, %d, %d, %12ld, %12ld, %12ld, %6.3g, %6.3g\n", 
        pf ? pf->file : "<file>", 
        pf ? pf->func : "<func>", 
        pf ? pf->line : 0,
        sd ? sd->level : 0,
    node->in_scope, 0L, inAllThreads, 
        0.0, (double)inAllThreads/(double)node->in_scope);

    PerThreadScopeTreeData* child;
    for (child = node->first_child; child; child = child->next_sibling) {
        node_show_stats(file, child, depth+1);
    }
}

void ROMP_show_stats(FILE* file)
{
    fprintf(file, "ROMP_show_stats %s\n", currentDateTime(false).c_str());
    fprintf(file, "file, line, level, in pf, in pflb, in_pfThread, pflb/elapsed, pft/elapsed\n");

    if (getMainFile())  {
        long mainDuration = mainTimer.nanoseconds();
        fprintf(file, "%s, %d, 0, %12ld, %12ld, %12ld, %6.3g, %6.3g\n", mainFile, mainLine, 
            mainDuration, 0L, 0L, 1.0, 1.0);
    }
        
    int tid;
    for (tid = 0; tid < ROMP_maxWatchedThreadNum; tid++) {
        PerThreadScopeTreeData* node = &scopeTreeRoots[tid];
        if (node) node_show_stats(file, node, 0);
    }

    fprintf(file, "ROMP_show_stats end\n");
}


void ROMP_Distributor_begin(ROMP_Distributor* distributor,
    int lo, int hi, 
    double* sumReducedDouble0, 
    double* sumReducedDouble1, 
    double* sumReducedDouble2) {
    
    distributor->originals[0] = sumReducedDouble0;
    distributor->originals[1] = sumReducedDouble1;
    distributor->originals[2] = sumReducedDouble2;

    distributor->partialSize = (hi - lo);
    if (distributor->partialSize <= 0) {
        distributor->partialSize = 0;
        return;
    }
    
    if (distributor->partialSize > ROMP_DISTRIBUTOR_PARTIAL_CAPACITY)
        distributor->partialSize = ROMP_DISTRIBUTOR_PARTIAL_CAPACITY;
    
    int i;
    for (i = 0; i < distributor->partialSize; i++) {
        int step = (hi - lo) / (distributor->partialSize - i);
        distributor->partials[i].lo = lo;
        distributor->partials[i].hi = (lo += step);
        int j;
        for (j = 0; j < ROMP_DISTRIBUTOR_REDUCTION_CAPACITY; j++) {
            distributor->partials[i].partialSum[j] = 0.0;
        }
    }
}

void ROMP_Distributor_end(ROMP_Distributor* distributor) {
    int j;
    for (j = 0; j < ROMP_DISTRIBUTOR_REDUCTION_CAPACITY; j++) {
        if (!distributor->originals[j]) continue;
        double sum = *distributor->originals[j];
        int i;
        for (i = 0; i < distributor->partialSize; i++) {
            sum += distributor->partials[i].partialSum[j];
        }   
        *distributor->originals[j] = sum;
    }
}



#if 0

// example
//
// pushd utils; rm -f a.out ; gcc -g -DHAVE_OPENMP -fopenmp -I../include romp_support.c timer.c ; ./a.out ; popd
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
    
    // Trivial timing 
    //
    if (0) {
        ROMP_PF_begin

        int i; double sum = 0;
        #pragma omp parallel for if_ROMP(assume_reproducible)
        for (i = 0; i < 4; i++) {
            int j;
            for (j = 0; j < 250000; j++) {   
                sum += i/(double)j;
            }
        }
        ROMP_PF_end
        
        rompExitHandler();
        return 0;
    }
    
    
    if (1) {
        ROMP_PF_begin

        int i; double sum = 0;
        #pragma omp parallel for if_ROMP(assume_reproducible)
        for (i = 0; i < 4; i++) {
            ROMP_PF_begin
            int j;
            for (j = 0; j < 250000; j++) {   
                sum += i/(double)j;
            }
            ROMP_PF_end
        }
        ROMP_PF_end
        
        return 0;
    }
    
    
    // Simple timing and control example
    //
    static const int v_size = 30000;
    int i;
    double sum = 0;

    fprintf(stdout, "#threads:%d\n", omp_get_max_threads());

    long threadMask = 0;
    ROMP_PF_begin
    #pragma omp parallel for if_ROMP(fast) reduction(+:sum)
    for (i = 0; i < v_size; i++) {
        ROMP_PFLB_begin
    
    threadMask |= 1 << 
#ifdef HAVE_OPENMP
        omp_get_thread_num();
#else
        0;
#endif  
        sum += 1.0 / i;
    
    int j;
    ROMP_PF_begin
        #pragma omp parallel for if_ROMP(assume_reproducible) reduction(+:sum)
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
    
    // Reproducible reduction example
    //
    ROMP_PF_begin
    
    double doubleToSum0_0,doubleToSum1_0;
    int numThreads;
    for (numThreads = 1; numThreads <= 10; numThreads++) {
    
        ROMP_PF_begin
        
        omp_set_num_threads(numThreads);
    
        double doubleToSum0 = 0.0, doubleToSum1 = 0.0;
        
        // the original loop control
        //
        #define ROMP_VARIABLE   originalVariable
        #define ROMP_LO         1
        #define ROMP_HI         10000
    
        // the original reductions
        //
        #define ROMP_SUMREDUCTION0  doubleToSum0
        #define ROMP_SUMREDUCTION1  doubleToSum1
        #define ROMP_LEVEL          assume_reproducible
        
#ifdef ROMP_SUPPORT_ENABLED
        const int romp_for_line = __LINE__;
#endif
        #include "romp_for_begin.h"
    
            #define doubleToSum0 ROMP_PARTIALSUM(0)
            #define doubleToSum1 ROMP_PARTIALSUM(1)
    
            // the original loop body
            //
            doubleToSum0 += 1.0 /             (double)(originalVariable)  ;
            doubleToSum1 += 1.0 / ( ROMP_HI - (double)(originalVariable) );

            #undef doubleToSum0
            #undef doubleToSum1
            
        #include "romp_for_end.h"

        if (numThreads == 1) {        
            doubleToSum0_0 = doubleToSum0;
            doubleToSum1_0 = doubleToSum1;
        } else {
            if (doubleToSum0_0 != doubleToSum0) fprintf(stderr, "diff 0\n");
            if (doubleToSum1_0 != doubleToSum1) fprintf(stderr, "diff 1\n");
        }
        
        printf("romp_for numThreads:%d doubleToSum0:%g doubleToSum1:%g\n",
            numThreads, doubleToSum0,doubleToSum1);
            
        ROMP_PF_end
    }
    ROMP_PF_end
    
    // Done
    //
    ROMP_show_stats(stdout);
    return 0;
}

#endif
