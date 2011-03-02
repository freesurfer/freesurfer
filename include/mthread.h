/**
 * @file  mthread.h
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
