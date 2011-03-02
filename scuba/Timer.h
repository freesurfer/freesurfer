/**
 * @file  Timer.h
 * @brief Simple stopwatch class
 *
 * Just call Start() when you want and TimeNow() when you want the
 * time between calls. Start() resets the timer.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.5 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


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

  // Call to start the timer.
  void Start ();

  // Returns ms between now and the time of the last call to Start().
  int  TimeNow ();

protected:
  struct timeb mStartTime;

};

#endif
