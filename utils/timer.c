#include <stdio.h>
#include <stdlib.h>

#include <sys/timeb.h>
#include "timer.h"
#include "proto.h"

struct timeb *
TimerStart(struct timeb *then)
{
  ftime(then) ;
  return(then) ;
}

int
TimerStop(struct timeb *then)
{
  struct timeb now ;
  int          msec, now_msec, then_msec ;

  ftime(&now) ;
  then_msec = then->time * 1000 + then->millitm ;
  now_msec = now.time * 1000 + now.millitm ;
  msec = now_msec - then_msec ;
  return(msec) ;
}
