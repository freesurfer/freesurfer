/**
 * @file  thread.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


/*
   @(#)thread.h 1.4
   8/10/95
*/
/*------------------------------------------------------------------------
      File Name:  thread.h

         Author:  Bruce Fischl

        Created:  Jan. 1994

    Description:

------------------------------------------------------------------------*/
#ifndef THREAD_H
#define THREAD_H

#include "queue.h"

#ifdef ANSI
int ThreadResume(int iTid, int iSignal) ;
int ThreadInit(int machine_id, int nthreads, int stacksize, int npriorities) ;
int ThreadEnd(int tid) ;
int ThreadStart(char *name, void (*func)(int iTid, void *parm),
                void *parm,int priority);
int ThreadSuspend(int tid, int iSignal) ;
int ThreadSleep(int tid, long usec) ;
int ThreadSignal(int tid, int signal) ;
int ThreadYield(void) ;
int  ThreadEnqueue(int tid, void *msg) ;
void *ThreadDequeue(int mode) ;
int  ThreadGetTid(void) ;
int  ThreadGetMid(void) ;
int  ThreadCheckMailbox(void) ;
void ThreadExit(int status) ;

#else

int ThreadInit() ;
int ThreadEnd() ;
int ThreadStart() ;
int ThreadSuspend() ;
int ThreadSleep() ;
int ThreadSignal() ;
int ThreadYield() ;
int  ThreadEnqueue() ;
void *ThreadDequeue() ;
int  ThreadGetTid() ;
int  ThreadGetMid() ;
int  ThreadCheckMailbox() ;

#endif
/* 0 and positive thread id #s are reserved for real threads */
#define TID_NONE      -3
#define TID_ALL       -2
#define TID_SELF      -1

#define SIG_ALL       0xffffffff
#define SIG_ANY       0xffffffff
#define SIG_Q_PENDING 0x00000001


#define MIN_PRIORITY  1
#define SECONDS(sec) ((long)sec * 1000000)
#define MSEC(msec)   ((long)msec * 1000)

#define THREAD_DONE   iDone
extern int iDone ;

#endif
