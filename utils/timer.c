#include <stdio.h>
#include <stdlib.h>

#include <sys/timeb.h>
#include "timer.h"
#include "proto.h"

struct timeb *
TimerStart(struct timeb *then)
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
  then->millitm = (tv.tv_usec+500)/1000;
  then->dstflag = 0;
  return then;
}

int
TimerStop(struct timeb *then)
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
  tvthen.tv_usec = ((long) (then->millitm))*1000;
  msec = 1000*(tvnow.tv_sec - tvthen.tv_sec) + (tvnow.tv_usec - tvthen.tv_usec+500)/1000;
  return msec;
}
