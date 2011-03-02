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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
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
