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
 *    $Author: kteich $
 *    $Date: 2007/10/22 04:39:30 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
