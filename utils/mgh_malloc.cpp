//  #define MEMORY_CORRUPTION_DEBUGGING     // uncomment this to turn the feature on



/**
 * @brief overrides or wraps malloc et. al. with code that allows us
 *        to measure their uses and debug their failures
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
#include "mgh_malloc.h"

#include "romp_support.h"


#ifdef HAVE_MCHECK
#include <mcheck.h>
#endif

#if !defined(__APPLE__)
#define USE_MPROTECT
#include <sys/mman.h>
#endif

#undef malloc
#undef free
#undef calloc
#undef realloc
#undef posix_memalign

typedef enum MallocAction {
    MA_insert = 1,
    MA_clear  = 2,
    MA_copy   = 4,
    MA_remove = 8
} MallocAction;

// Sometimes the code tracks where the allocation was done
//
typedef struct MallocStats {
    const char* file;
    const char* function;
    int         line;
    size_t      size;
    size_t      inserted;
    size_t      cleared;
    size_t      copied;
    size_t      removed;
    size_t      sumCurrentAllocs;
    size_t      maxSumCurrentAllocs;
} MallocStats;

static void clearMallocStatsAllocCounters();
static void showMallocStatsAllocCounters();



//  The link command needs to be something like this to enable this support
//
//      make CCLD="g++ -Wl,--wrap=free -Wl,--wrap=malloc -Wl,--wrap=calloc -Wl,--wrap=realloc -Wl,--wrap=posix_memalign"
//
//  In CMakeLists.txt for the entity being linked
//      SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,--wrap=free -Wl,--wrap=malloc -Wl,--wrap=calloc -Wl,--wrap=realloc -Wl,--wrap=posix_memalign")
// or   SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lmcheck") 

//
//  It may also need
//           -lmcheck


#if !defined(MEMORY_CORRUPTION_DEBUGGING) && !defined(DEBUG_MEMLEAK)

static void noteWhereAllocated(void* ptr, MallocStats* mallocStats) {
}

#else

//
// TBD get this working in conjunction with ROMP_SUPPORT


#include <mcheck.h>

// This is the GNU linker mechanism for intercepting calls to a function.
// When linked with the command line above, the original functions get renamed __real_name
// and the uses of the original functions get replaced with uses of __wrap_name
//
void *__real_malloc(size_t size);
void  __real_free(void *ptr);
void *__real_calloc(size_t nmemb, size_t size);
void *__real_realloc(void *ptr, size_t size);
int   __real_posix_memalign(void **memptr, size_t alignment, size_t size);

// The code here is only thread-safe because it is all done inside a critical section
// protected by this lock
//
#ifdef HAVE_OPENMP
static int lock_inited;
static omp_lock_t lock;
#endif

static void wrapper_lock_acq() {
#ifdef HAVE_OPENMP
   if (!lock_inited) { lock_inited = 1; omp_init_lock(&lock); }
   omp_set_lock(&lock);
#endif
}
static void wrapper_lock_rel() {
#ifdef HAVE_OPENMP
   omp_unset_lock(&lock);
#endif
}


#if defined(MEMORY_CORRUPTION_DEBUGGING)    
// In the extreme situation, rather than using __real_malloc etc., this code can
// use whole pages of memory to satisfy each request.  It can change the protection
// on one of these pages when it is 'free' to catch writes to free pages.  It also
// checksums the free page, and before it is reassigned, checks it has not been written while free.
// Sadly the o/s limits the number of pages that can be protected to just a few...
//
typedef struct Page {
    size_t size;                                // also available on not-first page
    size_t loGuard;
    size_t mostRecentlyAllocationClock;
    char   available[4096-sizeof(size_t)*3];
} Page;

size_t const pagesCapacity = 8000000;
size_t pagesSize = 0;
Page*  pages;

Page*  pageChains   [100];
size_t pageChainsGet[100];
size_t pageChainsPut[100];

size_t maxPageChainsIndexSeen = 10;
static void updatedMaxPageChainsIndexSeen() {
    fprintf(stderr, "pageChainsIndex:%ld\n", maxPageChainsIndexSeen);
}

static void getPageChainsIndexAndCount(size_t* index, size_t* count, size_t size) {
    size_t availablePerPage = sizeof(((Page*)0)->available);
    size_t neededPages = (size + availablePerPage - 1)/availablePerPage;
    size_t pageChainsIndex = 0;
    while ((1L << pageChainsIndex) < neededPages) pageChainsIndex++;
    if (maxPageChainsIndexSeen < pageChainsIndex) {
        maxPageChainsIndexSeen = pageChainsIndex;
        updatedMaxPageChainsIndexSeen();
    }
    *index = pageChainsIndex;
    *count = 1L << pageChainsIndex;
}

    
static size_t computePageHash(Page* loPage, size_t pageCount) {
    size_t count = pageCount*(4096/sizeof(size_t));
    size_t* p    = (size_t*)loPage;
    size_t hash  = 0, i;
    for (i = 0; i < count; i++) hash = (hash*123)^(*p++);
    return hash;
}

static void checkNotUsed(Page* nLoPage, size_t npcPageCount) {
    size_t loGuard = nLoPage->loGuard;
    nLoPage->loGuard = 0;
    if (loGuard != computePageHash(nLoPage, npcPageCount)) {
        *(int*)-1 = 0;                                                                  // has been changed when should not be available
    }
    nLoPage->loGuard = loGuard;
}

static size_t const memProtected_mostRecentlyAllocationClock = 12000221200L;
static Page*        memProtected_loPage = NULL;

static void canUse(Page* nLoPage, size_t npcPageCount) {
    
    if (nLoPage == memProtected_loPage) {
#ifdef USE_MPROTECT
        if (mprotect(nLoPage, npcPageCount*4096, PROT_READ|PROT_WRITE)) {   // the limit (see below) should not affect this
            perror("mprotect failed");                                      // 
            *(int*)-1 = 0;                                                  // if can protect, should be able to unprotect
        }
#endif
        memProtected_loPage = NULL;
    }
    
    checkNotUsed(nLoPage, npcPageCount);

}

static void shouldNotUse(Page* oLoPage, size_t opcPageCount) {

    oLoPage->loGuard = 0;
    oLoPage->loGuard = computePageHash(oLoPage, opcPageCount);

    if (oLoPage->mostRecentlyAllocationClock != memProtected_mostRecentlyAllocationClock) return;
    
#ifdef USE_MPROTECT
    if (mprotect(oLoPage, opcPageCount*4096, PROT_NONE)) {      // sadly there is a limit
        perror("mprotect failed");                              // on how variegated memory can be
        *(int*)-1 = 0;                // weird                  // but we can use it to find a repeatable problem
    }
#endif

    memProtected_loPage = oLoPage;
}

#endif

// The code also tracks all the mallocated storage in a huge hash table
// These can be either the pages above or storage managed by __real_malloc etc.
// This allows checking for double-frees, frees of non-malloc'ed storage, etc.
// It also puts a few guard bytes at the end of the storage that can be checked for overruns
//
typedef struct Malloced {
    size_t       size;
    size_t       padding;
    void*        ptr;
    MallocStats* whereAllocated;
} Malloced;

#define mallocedSizeLog2 23
#define mallocedSize     (1<<mallocedSizeLog2)
#define mallocedSizeMask (mallocedSize-1)

static Malloced* malloced = NULL;
static size_t countAlloc[_MAX_FS_THREADS];
static size_t countFrees[_MAX_FS_THREADS];
static size_t sumAlloc;
static size_t sumAllocLimit = 1000000;

static int doMallocAction(MallocAction action, void** ptr, size_t size) {

#if !defined(MEMORY_CORRUPTION_DEBUGGING)    
    if (action & MA_insert) countAlloc[omp_get_thread_num()]++;
    if (action & MA_remove) countFrees[omp_get_thread_num()]++;
    
    // Doesn't do the allocation
    //
    return 0;
#else
    if (!pages) {
        size_t unalignedPage = (size_t)    __real_malloc((pagesCapacity+1) * sizeof(Page));
        size_t alignedPage   = unalignedPage & ~(sizeof(Page)-1);
        if (unalignedPage > alignedPage) alignedPage += sizeof(Page);
        pages = (Page*)alignedPage;
    }

    char* op      = (action & MA_remove) ? (char*)*ptr : NULL;
    Page* oLoPage = (Page*)((size_t)op & ~(4096-1));
    Page* nLoPage = NULL;
    char* np      = NULL;

    size_t opci = 0, opcPageCount = 0; if (action & MA_remove) getPageChainsIndexAndCount(&opci, &opcPageCount, oLoPage->size);
    size_t npci = 0, npcPageCount = 0; if (action & MA_insert) getPageChainsIndexAndCount(&npci, &npcPageCount,          size);

    if (action & MA_insert) {                                                               // get the new pages
        size_t saveCountAlloc = countAlloc[omp_get_thread_num()]++;
        if (pageChains[npci]) {
            nLoPage = pageChains[npci];
            canUse(nLoPage, npcPageCount);
            pageChains[npci] = (Page*)nLoPage->size;
        } else {
            size_t newPagesSize = pagesSize + (1L << npci);
            if (newPagesSize > pagesCapacity) *(int*)-1 = 0;                                // out of memory
            nLoPage = &pages[pagesSize];
            pagesSize = newPagesSize;
        }
        pageChainsGet[npci]++;
        nLoPage->size    = size;
        nLoPage->loGuard = 0xFEEDBABEFEEDBABE;
        nLoPage->mostRecentlyAllocationClock = 
            ROMP_countGoParallel()*100000000L + saveCountAlloc*100 + omp_get_thread_num();
        np = (char*)&nLoPage->available;
    }

    if (action & MA_clear) {
        bzero(np, size);
    }

    if (action & MA_copy) {
        size_t minSize = (size < oLoPage->size) ? size : oLoPage->size; 
        memcpy(np, op, minSize);
    }

    if (action & MA_remove) {
        countFrees[omp_get_thread_num()]++;
        
        oLoPage->size = (size_t)pageChains[opci];                                           // put the old pages
        pageChains[opci] = oLoPage;
        pageChainsPut[opci]++;
        shouldNotUse(oLoPage, opcPageCount);
    }

    *ptr = np;
    
    return 1;           // DONE
#endif
}    


static Malloced* noteMallocAction(MallocAction action, void* ptr, size_t size, size_t padding) {

    if (!malloced) {
        malloced = (Malloced*)__real_calloc(mallocedSize,       sizeof(Malloced));
    }
    
    size_t stabs = 0;
    size_t hash = ((((size_t)ptr) >> 5)*75321) & mallocedSizeMask;
    Malloced* m = &malloced[hash];
    while (m->ptr != ptr) {
        if (m->ptr == 0) break;                                     // reached the end of the chain, some items on which may have been freed
        hash = (hash*327 + 1) & mallocedSizeMask;                   //      but their ptr's were left non-null so they did not act as end of chain!
        m = &malloced[hash];
        if (++stabs > 1000) *(int*)-1 = 0;  // table too full 
    }
    
    if (action & MA_remove) {
        if (m->ptr  != ptr) *(int*)-1 = 0;                          // not in chain 
        if (m->size == 0  ) *(int*)-1 = 0;                          // in chain but freed already
        char* pad = (char*)ptr + m->size - m->padding;              // 
        while (m->padding--) if (*pad++ != 0x5c) *(int*)-1 = 0;     // trailing padding has been corrupted
        sumAlloc -= m->size;
        m->whereAllocated = NULL;
        //m->ptr must not be zeroed' because that would end chains that others might be in
        m->size = 0;                        // free, but following items
    }                                       //  in chain still reachable
    
    if (action & MA_insert) {
        if (m->size != 0  ) *(int*)-1 = 0;  // already in chain
        if (size == 0) size = 1;            // special case, can malloc a 0 sized object!
        m->size = size;
        sumAlloc += m->size;
        m->ptr  = ptr;
        m->padding = padding;        
        char* pad = (char*)ptr + size - padding;
        while (padding--) *pad++ = 0x5c;
    }
    
    return m;
}

void round_up_size(size_t* size, size_t* padding) {
    *padding = 32 - (*size&31);
    *size += *padding;
}

static void showCurrentAllocs() {
    clearMallocStatsAllocCounters();
    
    size_t hash;
    for (hash = 0; hash < mallocedSize; hash++) {
        Malloced* m = &malloced[hash];
        if (!m->size || !m->whereAllocated) continue;
        m->whereAllocated->sumCurrentAllocs += m->size;
    }
    
    showMallocStatsAllocCounters();
}

static void noteWhereAllocated(void* ptr, MallocStats* mallocStats) {
    Malloced* malloced = noteMallocAction(0, ptr, 0, 0);
    malloced->whereAllocated = mallocStats;
    if (sumAllocLimit < sumAlloc) {
        sumAllocLimit = sumAlloc + 10000000L;
        fprintf(stderr,  "sumAlloc:%ld\n", sumAlloc);
        showCurrentAllocs();
        // Clear the counters
    }
}

// Here are the wrap functions that use the above mechanisms
// I have seen mris_fix_topology run almost to completion before the hash table filled up
//
void *__wrap_malloc(size_t size) {
    size_t padding;
    round_up_size(&size, &padding);
    wrapper_lock_acq();
    
    void* r = NULL; if (!doMallocAction(MA_insert, &r, size))                       r = __real_malloc(size);
    
    noteMallocAction(MA_insert, r, size, padding);
    wrapper_lock_rel();
    return r;
}


void *__wrap_calloc(size_t nmemb, size_t size) {
    size *= nmemb;
    size_t padding;
    round_up_size(&size, &padding);
    wrapper_lock_acq();
    void* r = NULL; if (!doMallocAction(MA_insert|MA_clear, &r, size))              r = __real_calloc(1, size);
    noteMallocAction(MA_insert, r, size, padding);
    wrapper_lock_rel();
    return r;
}

void *__wrap_realloc(void *ptr, size_t size) {
    if (!ptr) return __wrap_malloc(size);
    size_t padding;
    round_up_size(&size, &padding);
    wrapper_lock_acq();
    if (ptr) noteMallocAction(MA_remove, ptr, 0, 0);
    void* r = ptr;  if (!doMallocAction(MA_insert|MA_copy|MA_remove, &r, size))    r = __real_realloc(ptr, size);
    noteMallocAction(MA_insert, r, size, padding);
    wrapper_lock_rel();
    return r;
}

void  __wrap_free(void *ptr) {
    if (!ptr) return;
    wrapper_lock_acq();
    noteMallocAction(MA_remove, ptr, 0, 0);
    if (!doMallocAction(MA_remove, &ptr, 0))                                           __real_free(ptr);
    wrapper_lock_rel();
}

int __wrap_posix_memalign(void **memptr, size_t alignment, size_t size) {
    size_t padding;
    round_up_size(&size, &padding);
    wrapper_lock_acq();
    
    int   r = 0;    if (!doMallocAction(MA_insert, memptr, size))                   r = __real_posix_memalign(memptr, alignment, size);

    noteMallocAction(MA_insert, *memptr, size, padding);
    wrapper_lock_rel();
    return r;
}


#endif


// Regardless of whether the __real_malloc etc. or the __wrap_ ones are being used, it is still desirable
// to know where in the program the allocations are happening.  This mechanism allows that to happen.
//
#define mallocStatsSizeLog2 13
#define mallocStatsSize     (1<<mallocStatsSizeLog2)
#define mallocStatsSizeMask (mallocStatsSize-1)

static MallocStats* mallocStats = NULL;

static int stats_compare(const void* lhs_ptr, const void* rhs_ptr) {
   int lhs = *(int*)lhs_ptr;
   int rhs = *(int*)rhs_ptr;
   size_t lhsPriority = mallocStats[lhs].inserted + mallocStats[lhs].removed;
   size_t rhsPriority = mallocStats[rhs].inserted + mallocStats[rhs].removed;
   if (lhsPriority < rhsPriority) return +1;    // ascending order
   if (lhsPriority > rhsPriority) return -1;    // ascending order
   return 0;
}

static void mallocStatsExitHandler(void) {
    if (!mallocStats) return;

    size_t count = 0;
    size_t i;
    for (i = 0; i < mallocStatsSize; i++) {
        MallocStats* m = &mallocStats[i];
        if (!m->line) continue;
        count++;
    }

    int* indexs = (int*)malloc(count*sizeof(int));
    count = 0;
    for (i = 0; i < mallocStatsSize; i++) {
        MallocStats* m = &mallocStats[i];
        if (!m->line) continue;
        indexs[count++] = i;
    }

    qsort(indexs, count, sizeof(int), stats_compare);
       
    fprintf(stderr, "MallocStats\n   file, function, line, size, inserted, cleared, copied, removed\n");
    for (i = 0; i < count; i++) {
        MallocStats* m = &mallocStats[indexs[i]];
        fprintf(stderr, "%s, %s, %d, %g, %g, %g, %g, %g\n",
            m->file, m->function, m->line, (float)m->size, (float)m->inserted, (float)m->cleared, (float)m->copied, (float)m->removed);
    }
    fprintf(stderr, "MallocStats\n   file, function, line, size, inserted, cleared, copied, removed\n");
}

static MallocStats* noteMallocStatsAction(MallocAction action, size_t size, const char* file, const char* function, int line) {

    if (!mallocStats) {
        mallocStats = (MallocStats*)calloc(mallocStatsSize, sizeof(MallocStats));
        atexit(mallocStatsExitHandler);
    }
    
    size_t stabs = 0;
    size_t hash = (((size_t)line ^ (size_t)file ^ (size_t)function)*75321) & mallocStatsSizeMask;
    MallocStats* m = &mallocStats[hash];
    while (m->line != line || m->function != function || m->file != file) {
        if (m->line == 0) break; // not in chain
        hash = (hash*327 + 1) & mallocStatsSizeMask;
        m = &mallocStats[hash];
        if (++stabs > 1000) *(int*)-1 = 0;  // table too full 
    }
    
    if (m->line == 0) {     // There is a slight chance this might find the same empty cell as another thread - who cares?
        m->line = line; m->function = function; m->file = file;
    }
    
    if (action & MA_insert) { m->inserted++; m->size += size; }
    if (action & MA_clear)    m->cleared++;
    if (action & MA_copy)     m->copied++;
    if (action & MA_remove)   m->removed++;
    
    return m;
}

_Pragma("GCC diagnostic ignored \"-Wunused-function\"")
static void clearMallocStatsAllocCounters() {
    size_t hash;
    for (hash = 0; hash < mallocStatsSize; hash++) mallocStats[hash].sumCurrentAllocs = 0;
}

_Pragma("GCC diagnostic ignored \"-Wunused-function\"")
static void showMallocStatsAllocCounters() {
    fprintf(stdout, "\nshowMallocStatsAllocCounters\n");
    size_t hash;
    for (hash = 0; hash < mallocStatsSize; hash++) {
        MallocStats* m = &mallocStats[hash];
        if (!m->sumCurrentAllocs || m->maxSumCurrentAllocs < m->sumCurrentAllocs) continue;
        fprintf(stdout, "%10g %s:%d %s\n", (double)m->sumCurrentAllocs, m->file, m->line, m->function);
    }
    for (hash = 0; hash < mallocStatsSize; hash++) {
        MallocStats* m = &mallocStats[hash];
        if (!m->sumCurrentAllocs || m->maxSumCurrentAllocs > m->sumCurrentAllocs) continue;
        fprintf(stdout, "%10g %s:%d %s, up from %g\n", (double)m->sumCurrentAllocs, m->file, m->line, m->function, (double)m->maxSumCurrentAllocs);
        m->maxSumCurrentAllocs = m->sumCurrentAllocs;
    }
}

void *mallocHere (              size_t size, const char* file, const char* function, int line) {
    void* r = malloc(size);       MallocStats* m = noteMallocStatsAction(MA_insert, size, file, function, line);
    noteWhereAllocated(r,m);
    return r;
}

void  freeHere   (void *ptr,                 const char* file, const char* function, int line) {
    free(ptr);                                     noteMallocStatsAction(MA_remove, 0, file, function, line);
}

void* callocHere (size_t nmemb, size_t size, const char* file, const char* function, int line) {
    void* r = calloc(nmemb,size); MallocStats* m = noteMallocStatsAction(MA_insert|MA_clear, size, file, function, line);
    noteWhereAllocated(r,m);
    return r;
}

void *reallocHere(void *ptr,    size_t size, const char* file, const char* function, int line) {
    void* r = realloc(ptr,size);  MallocStats* m = noteMallocStatsAction(MA_insert|MA_remove, size, file, function, line);
    noteWhereAllocated(r,m);
    return r;
}

int posix_memalignHere(void **memptr, size_t alignment, size_t size, const char* file, const char* function, int line) {
    int r = posix_memalign(memptr, alignment, size); 
                                  MallocStats* m = noteMallocStatsAction(MA_insert, size, file, function, line);
    noteWhereAllocated(*memptr,m);
    return r;
}

