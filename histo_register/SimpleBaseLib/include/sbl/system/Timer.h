#ifndef _SBL_TIMER_H_
#define _SBL_TIMER_H_
#include <sbl/core/Dict.h>
namespace sbl {


/*! \file Timer.h
    \brief The Timer module provides utilities for measuring time 
    (for example, to measure the running time of an algorithm).
    These should not be used to measure extremely short periods of time.
*/


/// returns number of seconds (since boot-up?)
double getPerfTime();


/// wait for the given number of seconds (busy loop checking for interface events)
void delaySeconds( float seconds );


//-------------------------------------------
// INTERFACE TO PRECISE TIMING CODE
//-------------------------------------------


/// The Timer class is designed for timing parts of algorithms.
class Timer {
public:

    /// create (and, if requested, start) the timer
    explicit Timer( bool start = false );

    /// start the timer (assumes not already started)
    void start();

    /// stop the timer (assumes already started)
    void stop();

    /// get the total elapsed time
    float timeSum( bool stop = true );

    /// reset the total elapsed time to sero
    inline void reset() { m_timeSum = 0; }

    /// set the total elapsed time to the given value
    inline void setTimeSum( double timeSum ) { m_timeSum = timeSum; }

private:

    /// the current start time, if any
    double m_startTime;

    /// get the total elapsed time
    double m_timeSum;

    // disable copy constructor and assignment operator
    Timer( const Timer &x );
    Timer &operator=( const Timer &x );
};


//-------------------------------------------
// TIMER SET CLASS
//-------------------------------------------


/// The TimerSet class holds a named set of timers (useful for timing different parts of a complex algorithm).
class TimerSet {
public:

    /// start the specified timer
    void start( const String &timerName );

    /// stop the specified timer
    void stop( const String &timerName );

    /// get the total elapsed time of the specified timer
    float timeSum( const String &timerName );

    /// display the elapsed time of all of the timers
    void display( int indent, const String &caption );

private:

    /// timers by name
    StringDict<Timer> m_timerDict;

    // disable copy constructor and assignment operator
    TimerSet( const TimerSet &x );
    TimerSet &operator=( const TimerSet &x );
};


} // end namespace sbl
#endif // _SBL_TIMER_H_

