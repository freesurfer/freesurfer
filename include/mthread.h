/*
   @(#)mthread.h  1.2
   3/24/94
*/
/*------------------------------------------------------------------------
      File Name:  mthread.h

         Author:  Bruce Fischl

        Created:  Jan. 1994

    Description:  prototypes for machine-specific functions that are
                  required by the multi-threader.  These functions should
                  be implemented separately on each platform.

------------------------------------------------------------------------*/
#ifndef MACH_THREAD_H
#define MACH_THREAD_H

#ifdef ANSI
int MachThreadInit(int nthreads, int stacksize, int npriorities) ;
int MachThreadStart(char *name, void (*func)(int iTid, void *parm), int iTid, 
                    void *parm, int priority) ;
int MachThreadYield(int iTid) ;
int MachThreadKill(int iTid) ;
int MachThreadSleep(int iTid, long lUsec) ;
int MachThreadSuspend(int iTid) ;
int MachThreadResume(int iTid) ;
int MachThreadGetTid(void) ;
void MachThreadExit(int status) ;

#else

int MachThreadInit() ;
int MachThreadStart() ;
int MachThreadYield() ;
int MachThreadKill() ;
int MachThreadSleep() ;
int MachThreadSuspend() ;
int MachThreadResume() ;
int MachThreadGetTid() ;
void MachThreadExit() ;

#endif

#endif
