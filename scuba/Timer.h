#ifndef Timer_h
#define Timer_h

extern "C" {
#include <sys/time.h>
#include <sys/timeb.h>
}

class Timer {
 public:
  Timer ();
  ~Timer ();

  void Start ();
  int  TimeNow ();

 protected:
  struct timeb mStartTime;

};

#endif
