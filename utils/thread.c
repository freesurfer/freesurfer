/**
 * @file  thread.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 1.4 $
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
   @(#)thread.c 1.3
   4/7/94
*/
/*------------------------------------------------------------------------
      File Name: thread.c

         Author: Bruce Fischl

        Created: Jan. 1994

    Description: Generic multi-threading utility.  Does not contain any
                   platform-specific code.

------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "thread.h"
#include "mthread.h"
#include "queue.h"
#include "error.h"
#include "mailbox.h"

/*------------------------------------------------------------------------
                            CONSTANTS
------------------------------------------------------------------------*/

#define QSIZE   20

/*------------------------------------------------------------------------
                            STRUCTURES
------------------------------------------------------------------------*/

typedef struct
{
  char  *name ;
  int   iMtid ;          /* machine specific id of this thread */
  int   iSusSignal ;
  QUEUE *inQ ;          /* input queue for thread */
}
THREAD ;

/*------------------------------------------------------------------------
                            GLOBAL DATA
------------------------------------------------------------------------*/

int iDone = 0 ;   /* will cause shutdown when set to 1 */

/*------------------------------------------------------------------------
                            STATIC DATA
------------------------------------------------------------------------*/

static int      iMaxThreads = 0 ;   /* max allowable # of threads */
static int      iNthreads = 0 ;     /* number of current threads */
static THREAD   *pthrTable ;
static int      iMachineId = 0 ;

/*------------------------------------------------------------------------
                              FUNCTIONS
------------------------------------------------------------------------*/


/*------------------------------------------------------------------------
       Parameters:

      Description:
         initialize the thread utility, setting up a hard maximum on
         the number of allowable threads.

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
ThreadInit(int mid, int iMaxThr, int stacksize, int npriorities)
{
  THREAD *thread ;
  int    iError ;

  iMaxThreads = iMaxThr + 1 ;
  pthrTable = (THREAD *)calloc(iMaxThreads, sizeof(THREAD)) ;
  if (!pthrTable)
    return(ErrorSet(ERROR_NO_MEMORY,
                    "ThreadInit(%d, %d) - could not allocate table\n",
                    mid, iMaxThr)) ;

  iMachineId = mid ;
  iError = MachThreadInit(iMaxThr, stacksize, npriorities) ;
  if (iError < 0)
    return(iError) ;

  /* thread #0 is reserved for main */
  thread = &pthrTable[iNthreads++] ;
  thread->name = "main" ;
  thread->iMtid = MachThreadGetTid() ;
  thread->inQ = Qalloc(QSIZE) ;
  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
         Kill the thread specified by the iTid field.

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
ThreadEnd(int iTid)
{
  int  iMtid ;

  if (iTid >= iNthreads)
    return(ErrorSet(ERROR_NO_MEMORY,
                    "ThreadEnd(%d) - invalid thread #\n", iTid)) ;

  if (iTid == TID_SELF)
    iTid = ThreadGetTid() ;

  iMtid = pthrTable[iTid].iMtid ;

  return(MachThreadKill(iMtid)) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
         Start a new thread.  The entry point for the thread is specified
         by the 'func' parameter.  It will be passed it's thread id (iTid),
         and a parameter of its own choosing.

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
ThreadStart(char *name, void (*func)(int iTid, void *parm), void *parm,
            int priority)
{
  int     iMtid, iTid ;
  THREAD  *pthr ;

  if (iNthreads >= iMaxThreads)
    return(ErrorSet(ERROR_NO_MEMORY,
                    "ThreadStart() - no more threads available\n")) ;

  iTid = iNthreads++ ;
  pthr = &pthrTable[iTid] ;
  pthr->name = name ;
  pthr->inQ = Qalloc(QSIZE) ;
  if (!pthr->inQ)
    return(ErrorSet(ERROR_NO_MEMORY,
                    "ThreadStart(%s): could not allocate Q\n", name)) ;

  iMtid = MachThreadStart(name, func, iTid, parm, priority) ;
  if (iMtid < 0)
  {
    Qfree(pthr->inQ) ;
    iNthreads-- ;
    return(iMtid) ;
  }
  pthr->iMtid = iMtid ;

  ThreadResume(iTid, SIG_ALL) ;

  return(iTid) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
         Suspend the thread specified by iTid.  The thread will not gain
         control of the machine again until either ThreadResume is called
         explicitly, or ThreadSignal (not yet implemented) is called with
         the signal that this thread is waiting for.  Note that the iSignal
         field is a bit field, and any overlap between the incoming signal
         and that which is being waited for will result in the thread being
         returned to the active pool.

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
ThreadSuspend(int iTid, int iSignal)
{
  int     iMtid ;
  THREAD  *pthr ;

  if (iTid >= iNthreads)
    return(ErrorSet(ERROR_NO_MEMORY,
                    "ThreadSuspend(%d) - invalid thread\n", iTid)) ;

  if (iTid == TID_SELF)
    iTid = ThreadGetTid() ;

  pthr = &pthrTable[iTid] ;
  pthr->iSusSignal |= iSignal ;
  iMtid = pthr->iMtid ;

  return(MachThreadSuspend(iMtid)) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
        Wake up a thread that was previously slept or suspended.
        iSignal is currently unused.

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
ThreadResume(int iTid, int iSignal)
{
  THREAD  *pthr ;
  int     iMtid ;

  if (iTid >= iNthreads)
    return(ErrorSet(ERROR_NO_MEMORY,
                    "ThreadResume(%d) - invalid thread\n", iTid)) ;

  if (iTid == TID_ALL)
  {
    for (iTid = 0 ; iTid < iNthreads ; iTid++)
      ThreadResume(iTid, iSignal) ;

    return(0) ;
  }
  else
  {
    if (iTid == TID_SELF)
      iTid = ThreadGetTid() ;

    pthr = &pthrTable[iTid] ;
    iMtid = pthr->iMtid ;
    pthr->iSusSignal = 0 ;

    return(MachThreadResume(iMtid)) ;
  }
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
         Put the thread specified by iTid to sleep for 'usec' micro
         seconds.  The thread will wake up after the specified amount of
         time, or if ThreadResume is called in the interim.

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
ThreadSleep(int iTid, long usec)
{
  int  iMtid ;

  if (iTid >= iNthreads)
    return(ErrorSet(ERROR_NO_MEMORY,
                    "ThreadSleep(%d) - invalid thread #\n", iTid)) ;

  if (iTid == TID_SELF)
    iMtid = TID_SELF ;
  else
    iMtid = pthrTable[iTid].iMtid ;

  return(MachThreadSleep(iMtid, usec)) ;
}
/*------------------------------------------------------------------------
       Parameters:

       Description:
         Signal a thread that an event that is has (possibly) been waiting
         for has occurred.  If the thread iSusSignal field has any bits in
         common with the incoming signal than the thread will be woken up.

    Return Values:
          1 if the thread was woken up, 0 if not.

------------------------------------------------------------------------*/
int
ThreadSignal(int iTid, int iSignal)
{
  THREAD  *pthr ;

  if (iTid >= iNthreads)
    return(ErrorSet(ERROR_NO_MEMORY,
                    "ThreadSignal(%d) - invalid thread #\n", iTid)) ;

  if (iTid == TID_SELF)
    iTid = ThreadGetTid() ;

  pthr = &pthrTable[iTid] ;
  if (pthr->iSusSignal & iSignal)
  {
    ThreadResume(iTid, iSignal) ;
    return(1) ;
  }
  else
    return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
         Yield control of the machine.  The caller will not return from
         this function until other processes (if any) on the machine have
         run, and this thread is rescheduled.

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
ThreadYield(void)
{
  int  iMtid, iTid ;

  /* ThreadCheckMailbox() ;*/
  iTid = ThreadGetTid() ;
  iMtid = pthrTable[iTid].iMtid ;

  return(MachThreadYield(iMtid)) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
        Enqueue a message for the thread specified by iTid.  If iTid is
        set to TID_ALL, the message will be iteratively enqueued on all
        threads incoming queues.  If the threads in question are in a
        state of SIG_Q_PENDING (i.e. they have called Qget and are waiting
        for data), they will be woken up, and will return from the call.

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
ThreadEnqueue(int iTid, void *msg)
{
  THREAD *thr ;

  if (iTid >= iNthreads)
    return(ErrorSet(ERROR_NO_MEMORY,
                    "ThreadEnqueue(%d) - invalid thread #\n", iTid)) ;

  if (iTid == TID_ALL)
  {
    for (iTid = 0 ; iTid < iNthreads ; iTid++)
      ThreadEnqueue(iTid, msg) ;
  }
  else
  {
    thr = &pthrTable[iTid] ;
    if (msg)
      Qput(thr->inQ, msg) ;
    if (thr->iSusSignal & SIG_Q_PENDING)
      ThreadResume(iTid, SIG_Q_PENDING) ;
  }

  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
         Read an element off of the the thread's in queue.  If no data
         is present, and the caller has specified Q_WAIT_FOR_DATA,
         then suspend the thread until data is available.

    Return Values:
            NULL if nothing is in the Q, pointer to the object at the
               head of the Q otherwise.

------------------------------------------------------------------------*/
void *
ThreadDequeue(int mode)
{
  int    iTid ;
  THREAD *thr ;

  iTid = ThreadGetTid() ;

  if ((iTid >= iNthreads) || (iTid < 0))
  {
    ErrorSet(ERROR_NO_MEMORY,"ThreadDequeue(%d) - invalid thread # (max %d)\n",
             iTid, iNthreads) ;
    return(NULL) ;
  }

  thr = &pthrTable[iTid] ;

  if ((mode == Q_WAIT_FOR_DATA) && (Qempty(thr->inQ)))
  {
    ThreadSuspend(iTid, SIG_Q_PENDING) ;
  }

  return(Qget(thr->inQ, mode)) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
        Get the id of the currently running thread.

    Return Values:
            the TID on success (> 0), < 0 on error.

------------------------------------------------------------------------*/
int
ThreadGetTid(void)
{
  int  iMtid, iTid ;

  iMtid = MachThreadGetTid() ;

  for (iTid = 0 ; iTid  < iNthreads ; iTid++)
    if (pthrTable[iTid].iMtid == iMtid)
      return(iTid) ;

  return(TID_NONE) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
        Kill all threads and exit.  This call will not return.

    Return Values:

------------------------------------------------------------------------*/
void
ThreadExit(int status)
{
  MachThreadExit(status) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
        get the machine id.

    Return Values:
           the ID of the local machine.

------------------------------------------------------------------------*/
int
ThreadGetMid(void)
{
  return(iMachineId) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
        Check the threads incoming queue, dequeing all elements and
        invoking the appropriate mailbox handler for each.

    Return Values:
         the number of elements processed.

------------------------------------------------------------------------*/
int
ThreadCheckMailbox(void)
{
  MSG  *msg ;
  int  iNcalls = 0 ;

  msg = ThreadDequeue(Q_DONT_WAIT) ;
  while (msg)
  {
    iNcalls++ ;
    /*    MBinvoke(msg) ;*/
    msg = ThreadDequeue(Q_DONT_WAIT) ;
  }

  return(iNcalls) ;
}
