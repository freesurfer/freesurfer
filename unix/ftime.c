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

  t->millitm = 1000 * times(&tm) / HZ ;
#if 0
  t->millitm = 1000 * tm.tms_utime / HZ ;
#endif
  t->time = 0 ;
  return(time(NULL)) ;
}

#endif
