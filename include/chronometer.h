/*! \file
  Header file for simple C timer library
  Original Author : Richard Edgar
  $Id: chronometer.h,v 1.1 2009/07/13 21:25:33 krish Exp $
*/

#ifndef CHRONOMETER_H
#define CHRONOMETER_H

// These two system libraries are required
#include <sys/time.h>
#include <time.h>

//! Timer Datatype
typedef struct _Chronometer {
  int running;           //!< Whether the timer is running
  float accTime;         //!< Accumulated time for this timer (ms)
  long nStarts;          //!< Count of times this timer has been started
  struct timeval tStart; //!< When this timer started
} Chronometer;

// Prototypes

//! Initialises a timer
void InitChronometer( Chronometer *theChronometer );
//! Starts a timer
void StartChronometer( Chronometer *theChronometer );
//! Stops a timer
void StopChronometer( Chronometer *theChronometer );
//! Resets a timer
void ResetChronometer( Chronometer *theChronometer );
//! Gets the time for a timer in milliseconds
float GetChronometerValue( const Chronometer *theChronometer );
//! Gets the average time per start
float GetAverageChronometerValue( const Chronometer *theChronometer );

#endif
