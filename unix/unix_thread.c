#ifndef IRIX
#ifndef Linux
/*
   @(#)unix_thread.c  1.3
   3/24/94
*/
/*------------------------------------------------------------------------
      File Name:

         Author:

        Created:

    Description:

------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include <lwp/lwp.h>
#include <lwp/stackdep.h>
#include <lwp/check.h>
#include <lwp/lwpmachdep.h>
#include <lwp/lwperror.h>
#include "lwpproto.h"
#include "proto.h"

#include "thread.h"
#include "mthread.h"
#include "error.h"

/*------------------------------------------------------------------------
                            CONSTANTS
------------------------------------------------------------------------*/

#define ERRSTR   (lwp_errstr()[lwp_geterr()])

/*------------------------------------------------------------------------
                            STRUCTURES
------------------------------------------------------------------------*/
typedef struct
{
   thread_t  thr ;
    char      *name ;
} MACH_THREAD ;

/*------------------------------------------------------------------------
                            STATIC DATA
------------------------------------------------------------------------*/

static int          max_pri = 1 ;
static int          iMaxThreads = 0 ;
static int          iNthreads = 0 ;
static MACH_THREAD  *pmachThreadTable ;

/*------------------------------------------------------------------------
                            STATIC PROTOTYPES
------------------------------------------------------------------------*/

#define ANSI 1
#if ANSI
void print_stat(char *name, thread_t *iTid, statvec_t *stat) ;
#else
void print_stat() ;
#endif

/*------------------------------------------------------------------------
                              FUNCTIONS
------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
#if ANSI
MachThreadInit(int iMaxThr, int stacksize, int npri)
#else
MachThreadInit(iMaxThr, stacksize, npri)
    int iMaxThr ;
    int stacksize ;
    int npri ;
#endif
{
   MACH_THREAD *mthr ;

   if (lwp_setstkcache(2 * stacksize * sizeof(stkalign_t), iNthreads+1) < 0)
      return(ErrorSet(ERROR_NO_MEMORY, "lwp_setstkcache(%d, %d) failed (%s)",
                   2 * stacksize * sizeof(stkalign_t), iNthreads+1, ERRSTR)) ;

   iMaxThreads = iMaxThr + 1 ;
   pmachThreadTable = (MACH_THREAD *)calloc(iMaxThreads, sizeof(MACH_THREAD)) ;
   if (!pmachThreadTable)
      return(ErrorSet(ERROR_NO_MEMORY, 
                      "MachThreadInit: could not allocate %d threads\n",
                      iMaxThreads)) ;

   max_pri = npri + 1 ;
   pod_setmaxpri(max_pri) ;

   /* thread 0 is reserved for main */
   mthr = &pmachThreadTable[iNthreads++] ;
   mthr->name = "main" ;
   lwp_self(&mthr->thr) ;
   
   return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
            The TID of the just started thread, < 0 on error.

------------------------------------------------------------------------*/
int
#ifdef ANSI
MachThreadStart(char *name, void (*func)(int iTid, void *parm), int iTid, 
                        void *parm, int priority)
#else
MachThreadStart(name, func, iTid, parm, priority)
    char *name ;
    void (*func)() ;
    int  iTid ;
    void *parm ;
    int  priority ;
#endif
{
   thread_t   *pthr ;
   stkalign_t *stack ;

   if (iNthreads >= iMaxThreads)
      return(ErrorSet(ERROR_NO_MEMORY, "MachThreadStart: no threads left\n")) ;

   stack = lwp_newstk() ;
   if (!stack)
      return(ErrorSet(ERROR_NO_MEMORY,
                              "MachThreadStart(%s) could not allocate stack\n", name));

   pthr = &pmachThreadTable[iNthreads].thr ;
   pmachThreadTable[iNthreads].name = name ;

   if (lwp_create(pthr, func, priority, LWPSUSPEND, stack, 2, iTid, parm) < 0)
      return(ErrorSet(-1, 
                              "MachThreadStart: could not create stack (%s)\n", ERRSTR));

   return(iNthreads++) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
#ifdef ANSI
MachThreadYield(int iTid)
#else
MachThreadYield(iTid)
    int iTid ;
#endif
{
   thread_t  *iMtid ;

   if (iTid >= iNthreads)
      return(ErrorSet(ERROR_NO_MEMORY, 
                      "MachThreadYield(%d) - invalid thread #\n", iTid)) ;

   if (iTid == TID_SELF)
      iMtid = &SELF ;
   else
      iMtid = &pmachThreadTable[iTid].thr ;

   if (lwp_yield(*iMtid) < 0)
      return(ErrorSet(-4, "lwp_yeild failed: %s\n", ERRSTR));

   return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
#ifdef ANSI
MachThreadKill(int iTid)
#else
MachThreadKill(iTid)
    int iTid ;
#endif
{
   thread_t  *iMtid ;

   if (iTid >= iNthreads)
      return(ErrorSet(ERROR_NO_MEMORY, 
                      "MachThreadKill(%d) - invalid thread #\n", iTid)) ;

   if (iTid == TID_SELF)
      iMtid = &SELF ;
   else
      iMtid = &pmachThreadTable[iTid].thr ;

   if (lwp_destroy(*iMtid) < 0)
      return(ErrorSet(-4,"lwp_destroy failed: %s\n",ERRSTR));
   return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
#ifdef ANSI
MachThreadSleep(int iTid, long usec)
#else
MachThreadSleep(iTid, usec)
    int  iTid ;
    long usec ;
#endif
{
   struct timeval timeout ;

   timeout.tv_sec = (int)(usec / 1000000) ;
   timeout.tv_usec = (int)(usec % 1000000) ;
   if (lwp_sleep(&timeout) < 0)
      return(ErrorSet(-1, "lwp_sleep failed (%s)\n",ERRSTR));

   return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
#ifdef ANSI
MachThreadSuspend(int iTid)
#else
MachThreadSuspend(iTid)
    int iTid ;
#endif
{
   thread_t  *iMtid ;

   if (iTid >= iNthreads)
     return(ErrorSet(ERROR_NO_MEMORY, 
                     "MachThreadSuspend(%d) - invalid thread #\n", iTid)) ;

   if (iTid == TID_SELF)
      iMtid = &SELF ;
   else
      iMtid = &pmachThreadTable[iTid].thr ;

   if (lwp_suspend(*iMtid) < 0)
    return(ErrorSet(-2,"MachThreadSuspend: lwp_suspend failed (%s)\n",ERRSTR));

   return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
            nothing

------------------------------------------------------------------------*/
void
#ifdef ANSI
print_stat(char *name, thread_t *pthr, statvec_t *stat)
#else
print_stat(name, pthr, stat)
    char *name ;
    thread_t *pthr ;
    statvec_t *stat ;
#endif
{
   printf("%s (%d @ %lx) :  priority %d, ", 
             name, pthr->thread_key, (long)(pthr->thread_id), stat->stat_prio);

   switch (stat->stat_blocked.any_kind)
   {
   case NO_TYPE:
      printf("not blocked\n") ;
      break ;
   case CV_TYPE:
      printf("condition variable blocked\n") ;
      break ;
   case MON_TYPE:
      printf("monitor blocked\n") ;
      break ;
   case LWP_TYPE:
      printf("thread or agent blocked\n") ;
      break ;
   default:
      printf("unknown (%d) blockage\n", stat->stat_blocked.any_kind) ;
      break ;
   }
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
#ifdef ANSI
MachThreadResume(int iTid)
#else
MachThreadResume(iTid)
    int iTid ;
#endif
{
   thread_t  *iMtid ;

   if (iTid >= iNthreads)
      return(ErrorSet(ERROR_NO_MEMORY, 
                      "MachThreadResume(%d) - invalid thread #\n", iTid)) ;

   if (iTid == TID_SELF)
      iMtid = &SELF ;
   else
      iMtid = &pmachThreadTable[iTid].thr ;

   if (lwp_resume(*iMtid) < 0)
     return(ErrorSet(-2, "MachThreadResume: lwp_resume failed (%s)\n",ERRSTR));

   return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
            the ID of the current thread on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
#ifdef ANSI
MachThreadGetTid(void)
#else
MachThreadGetTid()
#endif
{
   int      iTid ;
   thread_t self ;

   if (lwp_self(&self) < 0)
      ErrorSet(-2, "MachThreadGetTid: lwp_self failed (%s)\n", ERRSTR) ;
   for (iTid = 0 ; iTid < iNthreads ; iTid++)
   {
      if (SAMETHREAD(self, pmachThreadTable[iTid].thr))
         return(iTid) ;
   }
   return(TID_NONE) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
            nothing.

------------------------------------------------------------------------*/
void
#ifdef ANSI
MachThreadExit(int status)
#else
MachThreadExit(status)
    int status ;
#endif
{
   pod_exit(status) ;
}
#else
#include "mthread.h"

int MachThreadInit(int nthreads, int stacksize, int npriorities) {return(1);}
int MachThreadStart(char *name, void (*func)(int iTid, void *parm), int iTid, 
                    void *parm, int priority) {return(1) ;}
int MachThreadYield(int iTid) {return(1) ;}
int MachThreadKill(int iTid) {return(1) ;}
int MachThreadSleep(int iTid, long lUsec) {return(1) ;}
int MachThreadSuspend(int iTid) {return(1) ;}
int MachThreadResume(int iTid) {return(1) ;}
int MachThreadGetTid(void) {return(1) ;}
void MachThreadExit(int status) {}

#endif
#endif
