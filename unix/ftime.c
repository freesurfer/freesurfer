#ifdef IRIX

#include <sys/types.h>
#include <sys/time.h>
#include <sys/timeb.h>

int 
ftime(struct timeb *t)
{
  t->millitm = 0 ;
  return(time(NULL)) ;
}

#endif
