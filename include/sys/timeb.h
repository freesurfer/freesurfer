/* @(#)timeb.h 2.6 88/08/19 SMI; from UCB 4.2 81/02/19 */

/*
 * Structure returned by ftime system call
 */

#ifndef _sys_timeb_h
#define _sys_timeb_h

#include <sys/time.h>

struct timeb
{
  time_t time;
  unsigned short millitm;
  short timezone;
  short dstflag;
};

#if defined(IRIX)
int ftime(struct timeb *t) ;
#endif

#endif /*!_sys_timeb_h*/
