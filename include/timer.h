#ifndef TIMER_H
#define TIMER_H

#include <sys/timeb.h>

struct timeb *TimerStart(struct timeb *then) ;
int TimerStop(struct timeb *then) ;

#ifdef Linux
/* don't know why this doesn't work on linux, but.... */
extern int ftime (struct timeb *__timebuf);
#endif

#endif
