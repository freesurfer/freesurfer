/*! \file
  C file to provide some simple timers, modelled on those in CUTIL
  @author R. G. Edgar
  Commit Info : $Id: chronometer.c,v 1.1 2009/07/13 21:25:45 krish Exp $
*/
// Modelled after those in CUTIL, but not as elegant

#include <stdlib.h>

#include "chronometer.h"

//! Stopped state for timer
#define kStopped 0
//! Running state for timer
#define kRunning 1

// Private Prototypes

//! Get time difference in milliseconds between two timeval structures
static float GetMilliSeconds(const struct timeval t1, const struct timeval t2);

// ===========================================================

float GetMilliSeconds(const struct timeval t1, const struct timeval t2)
{
  // Returns the difference in ms between two timevals
  float secs, musecs, msecs;

  secs = t2.tv_sec - t1.tv_sec;
  musecs = t2.tv_usec - t1.tv_usec;

  msecs = (1e3 * secs) + (musecs / 1e3);

  return (msecs);
}

// ===========================================================

void InitChronometer(Chronometer *theChronometer)
{
  // Prepares a timer for use

  theChronometer->running = kStopped;
  theChronometer->accTime = 0;
  theChronometer->nStarts = 0;
}

// ===========================================================

void StartChronometer(Chronometer *theChronometer)
{
  // Starts the timer running

  // Sanity check
  if (theChronometer->running == kRunning) {
    return;
  }

  gettimeofday(&(theChronometer->tStart), NULL);
  theChronometer->running = kRunning;
  theChronometer->nStarts++;
}

// ===========================================================

void StopChronometer(Chronometer *theChronometer)
{
  // Stops the current timer

  struct timeval tStop;

  // Sanity check
  if (theChronometer->running == kStopped) {
    return;
  }

  // Get the stop time
  gettimeofday(&tStop, NULL);
  theChronometer->running = kStopped;

  // Add on the accumulated time
  theChronometer->accTime += GetMilliSeconds(theChronometer->tStart, tStop);
}

// ============================================================

void ResetChronometer(Chronometer *theChronometer)
{
  // For the C library, this is a wrapper

  InitChronometer(theChronometer);
}

// ============================================================

float GetChronometerValue(const Chronometer *theChronometer)
{
  // Returns the time in milliseconds
  struct timeval now;
  float msecs;

  msecs = theChronometer->accTime;

  if (theChronometer->running == kRunning) {
    // Have to add on a little more
    gettimeofday(&now, NULL);
    msecs += GetMilliSeconds(theChronometer->tStart, now);
  }

  return (msecs);
}

// ============================================================

float GetAverageChronometerValue(const Chronometer *theChronometer)
{
  // Returns the time in milliseconds averaged over the number of starts
  float msecs;

  msecs = GetChronometerValue(theChronometer);

  return (msecs / theChronometer->nStarts);
}
