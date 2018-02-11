/*

Gheorghe Postelnicu, 2007

simple timer class - basic time measurements

*/
#ifndef _h_simple_timer_h
#define _h_simple_timer_h


// use FS timer
extern "C"
{
#include "timer.h"
};

class SimpleTimer
{
public:
  SimpleTimer()
  {
    TimerStart(&m_start);
  }
  double elapsed() const
  {
    return double( TimerStop(const_cast<timeb*>(&m_start) ) ) / 1000.;
  }
  double elapsed_min() const
  {
    return this->elapsed() / 60.;
  }
  double elapsed_hour() const
  {
    return this->elapsed() / 3600. ;
  }
private:
  timeb m_start;
};

#endif
