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
