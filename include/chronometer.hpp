/*! \file
  Header file for simple C++ timer library
*/

#ifndef CHRONOMETER_H
#define CHRONOMETER_H
// These two system libraries are required
#include <sys/time.h>
#include <time.h>

#include <iostream>


namespace SciGPU
{

//! Namespace to hold utilities for SciGPU
namespace Utilities
{

//! Timer Datatype
/*!
  This is a very simple timer, which should provide
  accuracy on the order of 10 microseconds
*/
class Chronometer
{

public:
  //! Default constructor
  Chronometer( void ) : running(false),
    accTime(0),
    nStarts(0),
    tStart()
  {
  }

  //! Start a chronometer
  void Start( void )
  {
    // Sanity check
    if( running )
    {
      return;
    }

    running = true;
    nStarts++;
    gettimeofday( &(tStart), NULL );
  }

  //! Stop a chronometer
  void Stop( void )
  {
    // Sanity check
    if( !running )
    {
      return;
    }

    struct timeval tStop;
    gettimeofday( &tStop, NULL );

    running = false;
    accTime += GetMilliseconds( tStart, tStop );
  }

  //! Reset a chronometer
  void Reset( void )
  {
    running = false;
    accTime = 0;
    nStarts = 0;
  }


  //! Returns the current accumulated time
  float GetTime( void ) const
  {
    float msecs;

    msecs = accTime;

    if( running )
    {
      struct timeval now;

      gettimeofday( &now, NULL );
      msecs += GetMilliseconds( tStart, now );
    }

    return( msecs );
  }

  //! Returns the average accumulated time
  float GetAverageTime( void ) const
  {
    float msecs;

    msecs = GetTime();

    return( msecs / nStarts );
  }

  friend std::ostream& operator<<( std::ostream& os,
                                   const Chronometer& me );
private:
  bool running;          //!< Whether the timer is running
  float accTime;         //!< Accumulated time for this timer (ms)
  long nStarts;          //!< Count of times this timer has been started
  struct timeval tStart; //!< When this timer started


  //! Function to return the difference in milliseconds between two timevals
  float GetMilliseconds( const struct timeval t1, const struct timeval t2 ) const
  {
    float secs, musecs, msecs;

    secs = t2.tv_sec - t1.tv_sec;
    musecs = t2.tv_usec - t1.tv_usec;

    msecs = (1e3*secs) + (musecs/1e3);

    return( msecs );
  }
};

}
}


#endif
