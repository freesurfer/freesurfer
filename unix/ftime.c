#ifdef IRIX

#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/timeb.h>
#include <sys/param.h>

int 
ftime(struct timeb *t)
{
  struct tms tm ;
  struct timeval tv;
  struct timezone tz;
  int    msec ;
  int    sec ;

  tz.tz_minuteswest = 0;
  tz.tz_dsttime = 0;
  gettimeofday(&tv,&tz);
  msec = (int)((double)tv.tv_usec * 0.001 + 0.5) ;
  sec = (int)((double)tv.tv_sec+0.5) ;
  t->millitm = msec ;
  t->time = sec ;
  return(time(NULL)) ;
}

#endif
