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
