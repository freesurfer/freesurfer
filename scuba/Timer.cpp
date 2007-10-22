/**
 * @file  Timer.cpp
 * @brief Simple stopwatch class
 *
 * Just call Start() when you want and TimeNow() when you want the
 * time between calls. Start() resets the timer.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/22 04:39:30 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include "Timer.h"

Timer::Timer () {
}


Timer::~Timer () {
}


void
Timer::Start () {

  struct timeval tv;
  gettimeofday(&tv, 0);
  mStartTime.time = tv.tv_sec;
  mStartTime.millitm = (tv.tv_usec+500)/1000;
  mStartTime.dstflag = 0;
}

int
Timer::TimeNow () {

  int msec;
  struct timeval tvnow, tvthen;
  gettimeofday(&tvnow, 0);
  tvthen.tv_sec = mStartTime.time;
  tvthen.tv_usec = ((long) (mStartTime.millitm))*1000;
  msec = 1000*(tvnow.tv_sec - tvthen.tv_sec) + (tvnow.tv_usec - tvthen.tv_usec+500)/1000;
  return msec;
}
