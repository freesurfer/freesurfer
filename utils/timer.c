/**
 * @file  timer.c
 * @brief TimerStart and TimerStop routines for timing programs
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 1.6 $
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

#include <stdio.h>
#include <stdlib.h>

#include "proto.h"
#include "sys/timeb.h"
#include "timer.h"


// clock_gettime is not available on osx < 10.12 sierra
// mach_gettime is a custom replacement for that
#ifndef HAVE_CLOCK_GETTIME
#include <mach/clock.h>
#include <mach/mach.h>

int mach_gettime(clockid_t clk_id, struct timespec *tp)
{
  int ret;
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), clk_id, &cclock);
  ret = clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  tp->tv_sec = mts.tv_sec;
  tp->tv_nsec = mts.tv_nsec;
  return ret;
}
#endif


struct timeb *TimerStart(struct timeb *then)
{
/* according to the header ftime() is obsolete */
#if 0
  ftime(then) ;
  return(then) ;
#endif
  /* struct timeval { long tv_sec; long tv_usec; }; */
  struct timeval tv;
  gettimeofday(&tv, 0);
  then->time = tv.tv_sec;
  then->millitm = (tv.tv_usec + 500) / 1000;
  then->dstflag = 0;
  return then;
}


int TimerStop(struct timeb *then)
{
#if 0
  struct timeb now ;
  int          msec, now_msec, then_msec ;

  ftime(&now) ;
  then_msec = (int)then->time * 1000 + (int)then->millitm ;
  now_msec = (int)now.time * 1000 + (int)now.millitm ;
  msec = now_msec - then_msec ;
  return(msec) ;
#endif
  int msec;
  struct timeval tvnow, tvthen;
  gettimeofday(&tvnow, 0);
  tvthen.tv_sec = then->time;
  tvthen.tv_usec = ((long)(then->millitm)) * 1000;
  msec = 1000 * (tvnow.tv_sec - tvthen.tv_sec) + (tvnow.tv_usec - tvthen.tv_usec + 500) / 1000;
  return msec;
}


void TimerStartNanosecs(struct NanosecsTimer * nst)
{
#ifndef HAVE_CLOCK_GETTIME
  mach_gettime(CALENDAR_CLOCK, &nst->now);
#else
  clock_gettime(CLOCK_REALTIME, &nst->now);
#endif
}


struct Nanosecs TimerElapsedNanosecs(struct NanosecsTimer * nst) // returns delta in nanosecs
{
  struct timespec now;
#ifndef HAVE_CLOCK_GETTIME
  mach_gettime(CALENDAR_CLOCK, &now);
#else
  clock_gettime(CLOCK_REALTIME, &now);
#endif
  struct timespec * then = &nst->now;
  struct Nanosecs result;
  result.ns = (long)(now.tv_sec - then->tv_sec)*1000000000 + (long)(now.tv_nsec - then->tv_nsec);
  return result;
}


// Use this string for writing into files etc. where it might be compared by
// the test system.  The string is usually the current date and time but
// can be overridden with a string that can the same across runs.
//
const char* current_date_time() {
  time_t tt = time(&tt);
  const char* time_str = ctime(&tt);
  const char* override_time_str =
    getenv("FREESURFER_REPLACEMENT_FOR_CREATION_TIME_STRING");
  if (override_time_str) time_str = override_time_str;
  return time_str;
}
