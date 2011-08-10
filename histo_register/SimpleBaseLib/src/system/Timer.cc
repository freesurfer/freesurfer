// Licensed under MIT license; see license.txt.

#include <sbl/system/Timer.h>
#include <sbl/core/Command.h> // for checkCommandEvents() inside delaySeconds()
#ifdef WIN32
    #include <windows.h>
#else
    #include <sys/time.h>
#endif
namespace sbl {


/// returns number of seconds (since boot-up?)
double getPerfTime() {
#ifdef WIN32
    static double freqVal = 0;
    static bool freqSet = false;

    // get freq
    if (freqSet == false) {
        static LARGE_INTEGER freq;
        QueryPerformanceFrequency(&freq);
        freqVal = (double) freq.QuadPart;
        freqSet = true;
    }

    // get time
    LARGE_INTEGER curTime;
    QueryPerformanceCounter(&curTime);
    return ((double) curTime.QuadPart) / freqVal;
#else
    struct timeval tv;
    gettimeofday( &tv, NULL );
    double secs = (double) tv.tv_sec + 0.000001 * (double) tv.tv_usec;
    return secs;
#endif
}


/// wait for the given number of seconds (busy loop checking for interface events)
void delaySeconds( float seconds ) {
    double startTime = getPerfTime();
    while (getPerfTime() < startTime + seconds)
        checkCommandEvents();
}


//-------------------------------------------
// INTERFACE TO PRECISE TIMING CODE
//-------------------------------------------


/// create (and, by default, start) the timer
Timer::Timer( bool start ) {
    m_startTime = -1;
    m_timeSum = 0;
    if (start)
        Timer::start();
}


/// start the timer (assumes not already started)
void Timer::start() {
    if (m_startTime > 0)
        fatalError( "called timeSum::start after already started" );
    m_startTime = getPerfTime();
}


/// stop the timer (assumes already started)
void Timer::stop() {
    if (m_startTime < 0)
        fatalError( "called timeSum::stop without starting" );
    m_timeSum += getPerfTime() - m_startTime;
    m_startTime = -1;
}


/// get the total elapsed time
float Timer::timeSum( bool stop ) {
    if (stop && m_startTime > 0)
        Timer::stop();
    return (float) m_timeSum;
}


//-------------------------------------------
// TIMER SET CLASS
//-------------------------------------------


/// start the specified timer
void TimerSet::start( const String &timerName ) {
    Timer *timer = m_timerDict.find( timerName );
    if (timer)
        timer->start();
    else
        m_timerDict.add( timerName, new Timer );
}


/// stop the specified timer
void TimerSet::stop( const String &timerName ) {
    Timer *timer = m_timerDict.find( timerName );
    if (timer) 
        timer->stop();
    else
        warning( "timer not found: %s", timerName.c_str() );
}


/// get the total elapsed time of the specified timer
float TimerSet::timeSum( const String &timerName ) {
    float timeSum = 0;
    Timer *timer = m_timerDict.find( timerName );
    if (timer) 
        timeSum = timer->timeSum();
    else
        warning( "timer not found: %s", timerName.c_str() );
    return timeSum;
}


/// display the elapsed time of all of the timers
void TimerSet::display( int indent, const String &caption ) {
    Array<String> dispList;
    for (int i = 0; i < m_timerDict.count(); i++) 
        dispList.appendCopy( sprintF( "%s: %1.2f", m_timerDict.refArrayKey( i ).c_str(), m_timerDict.refArray( i ).timeSum() ));
    String dispText;
    if (dispList.count())
        dispText = caption + join( dispList, ", " );
    else
        dispText = caption + "[none]";
    disp( 1, dispText.c_str() );
}


} // end namespace sbl


