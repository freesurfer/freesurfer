#ifndef _SBL_TIME_UTIL_H_
#define _SBL_TIME_UTIL_H_
#include <sbl/core/String.h>
namespace sbl {


/*! \file TimeUtil.h
    \brief The TimeUtil module provides utilities for manipulating timestamps.
    These timestamps are represented using the unix standard of seconds since the epoch.
    In future versions the int timestamp types will be replaced with 64-bit ints and/or doubles.
*/


//-------------------------------------------
// CONVERT TIME INTERVALS
//-------------------------------------------


/// convert days to seconds
inline int daysToSecs( int days ) { return days * 24 * 60 * 60; }
inline float daysToSecs( float days ) { return days * 24.0f * 60.0f * 60.0f; }


/// convert hours to seconds
inline int hoursToSecs( int hours ) { return hours * 60 * 60; }
inline float hoursToSecs( float hours ) { return hours * 60.0f * 60.0f; }


/// convert minutes to seconds
inline int minutesToSecs( int minutes ) { return minutes * 60; }
inline float minutesToSecs( float minutes ) { return minutes * 60.0f; }


/// get the minute of the day 
inline int minuteOfDay( int hour, int minute ) { return hour * 60 + minute; }


//-------------------------------------------
// CREATE TIMESTAMPS
//-------------------------------------------


/// returns the current timestamp
int now();


/// creates timestamp corresponding to date (midnight)
int dateUTC( int year, int month, int day );
int dateUTC( const String &dateStr );


/// creates timestamp from date and time
int dateTimeLocal( int year, int month, int day, int hour, int minute, int second );
int dateTimeUTC( int year, int month, int day, int hour, int minute, int second );


// parse a date-time string into a timestamp, assuming YMD-HMS order
int parseTimestamp( const String &data, const String &format, bool local, int &msec );


//-------------------------------------------
// EXTRACT COMPONENT OF TIMESTAMP
//-------------------------------------------


/// get day of week (0=Sunday, 1=Monday, ...)
int weekDayUTC( int timestamp );
int weekDayLocal( int timestamp );


/// get hour of day
int hourOfDayLocal( int timestamp );


/// get minute of day
int minuteOfDayLocal( int timestamp );


/// get second of day
int secondOfDayLocal( int timestamp );


/// get fraction of day (from midnight to midnight)
float fracOfDayUTC( int timestamp );
float fracOfDayLocal( int timestamp );


/// returns timestamp corresponding to start of day (UTC) 
int datePartUTC( int timestamp );


//-------------------------------------------
// FORMAT TIMESTAMPS AS STRINGS
//-------------------------------------------


/// convert timestamp to local-time string 
String tsToStrLocal( int timestamp, bool showSeconds = true );
inline String tsToStrLocal( double timestamp, bool showSeconds = true ) { return tsToStrLocal( (int) timestamp, showSeconds ); }
String tsToStrLocal24( int timestamp, bool showSeconds = true );
inline String tsToStrLocal24( double timestamp, bool showSeconds = true ) { return tsToStrLocal24( (int) timestamp, showSeconds ); }
String tsToStrDateLocal( int timestamp );


/// convert timestamp to local-time string using filename-safe delimiters
String tsToFileNameLocal( int timestamp );


/// convert timestamp to UTC string
String tsToStrUTC( int timestamp );
String tsToStrDateUTC( int timestamp );


//-------------------------------------------
// A SIMPLE DATE CLASS
//-------------------------------------------


/// The Date class represents a date value, allowing simple and efficient operations on dates.
class Date {
public:
    
    // basic constructors
    Date( int year, int month, int day );
    Date( int timestamp );

    /// date components
    inline int year() const { return m_year; }
    inline int month() const { return m_month; }
    inline int day() const { return m_day; }

    /// start of day / end of day timestamp
    inline int startTimestamp() const { return m_startTimestamp; }
    inline int endTimestamp() const { return m_startTimestamp + 24 * 60 * 60 - 1; }

    /// move date forward by one day
    void increment();

private:

    /// date components
    int m_year;
    int m_month;
    int m_day;

    /// the UTC starting timestamp (midnight)
    int m_startTimestamp;
};


} // end namespace sbl
#endif // _SBL_TIME_UTIL_H_

