#ifndef TIMER_H
#define TIMER_H

#include <sys/timeb.h>

struct timeb *TimerStart(struct timeb *then) ;
int TimerStop(struct timeb *then) ;

#endif
