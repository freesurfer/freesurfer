// Licensed under MIT license; see license.txt.

#include <sbl/system/TimeUtil.h>
#include <sbl/core/StringUtil.h>
#include <time.h>
#include <stdio.h>
#ifdef WIN32
    #include <windows.h>
    #include <external/win_timegm.h>
#else
    #include <sys/time.h>
#endif
namespace sbl {


// fix(clean): order in this file doesn't match header file


//-------------------------------------------
// DATE/TIME UTILS
//-------------------------------------------


/// convert timestamp to local-time string 
String tsToStrLocal( int timestamp, bool showSeconds ) {
    time_t tsLong = timestamp;
    char cstr[ 1000 ];
    if (timestamp > 0) {
        if (showSeconds)
            strftime( cstr, 1000, "%Y/%m/%d %I:%M:%S%p", localtime( &tsLong ));
        else
            strftime( cstr, 1000, "%Y/%m/%d %I:%M%p", localtime( &tsLong ));
    } else {
        sprintf( cstr, "[invalid timestamp: %d]", timestamp );
    }
    return String( cstr );
}


/// convert timestamp to local-time string 
String tsToStrLocal24( int timestamp, bool showSeconds ) {
    time_t tsLong = timestamp;
    char cstr[ 1000 ];
    if (timestamp > 0) {
        if (showSeconds)
            strftime( cstr, 1000, "%Y/%m/%d %H:%M:%S", localtime( &tsLong ));
        else 
            strftime( cstr, 1000, "%Y/%m/%d %H:%M", localtime( &tsLong ));
    } else {
        sprintf( cstr, "[invalid timestamp: %d]", timestamp );
    }
    return String( cstr );
}


/// convert timestamp to local-time string 
String tsToStrDateLocal( int timestamp ) {
    time_t tsLong = timestamp;
    char cstr[ 1000 ];
    strftime( cstr, 1000, "%Y-%m-%d", localtime( &tsLong ));
    return String( cstr );
}


/// convert timestamp to local-time string using filename-safe delimiters
String tsToFileNameLocal( int timestamp ) {
    time_t tsLong = timestamp;
    char cstr[ 1000 ];
    strftime( cstr, 1000, "%Y-%m-%d-%H-%M-%S", localtime( &tsLong ));
    return String( cstr );
}


/// convert timestamp to UTC string
String tsToStrUTC( int timestamp ) {
    time_t tsLong = timestamp;
    char cstr[ 1000 ];
    if (timestamp > 0) {
        strftime( cstr, 1000, "%Y/%m/%d %H:%M:%S", gmtime( &tsLong ));
    } else {
        sprintf( cstr, "[invalid timestamp: %d]", timestamp );
    }
    return String( cstr );
}


/// convert timestamp to UTC string
String tsToStrDateUTC( int timestamp ) {
    time_t tsLong = timestamp;
    char cstr[ 1000 ];
    strftime( cstr, 1000, "%Y-%m-%d", gmtime( &tsLong ));
    return String( cstr );
}


/// get day of week (0=Sunday, 1=Monday, ...)
int weekDayUTC( int timestamp ) {
    time_t tsLong = timestamp;
    return gmtime( &tsLong )->tm_wday;
}


/// get day of week (0=Sunday, 1=Monday, ...)
int weekDayLocal( int timestamp ) {
    time_t tsLong = timestamp;
    struct tm *timeTuple = localtime( &tsLong );
    return timeTuple->tm_wday;
}


/// get hour of day
int hourOfDayLocal( int timestamp ) {
    time_t tsLong = timestamp;
    return localtime( &tsLong )->tm_hour;
}


/// get minute of day
int minuteOfDayLocal( int timestamp ) {
    time_t tsLong = timestamp;
    struct tm *timeTuple = localtime( &tsLong );
    return timeTuple->tm_hour * 60 + timeTuple->tm_min;
}


/// get second of day
int secondOfDayLocal( int timestamp ) {
    time_t tsLong = timestamp;
    struct tm *timeTuple = localtime( &tsLong );
    return (timeTuple->tm_hour * 60 + timeTuple->tm_min) * 60 + timeTuple->tm_sec;
}


/// get fraction of day (from midnight to midnight)
float fracOfDayUTC( int timestamp ) {
    time_t tsLong = timestamp;
    struct tm *timeTuple = gmtime( &tsLong );
    int secondOfDay = (timeTuple->tm_hour * 60 + timeTuple->tm_min) * 60 + timeTuple->tm_sec;
    return (float) secondOfDay / (24.0f * 60.0f * 60.0f);
}


/// get fraction of day (from midnight to midnight)
float fracOfDayLocal( int timestamp ) {
    time_t tsLong = timestamp;
    struct tm *timeTuple = localtime( &tsLong );
    int secondOfDay = (timeTuple->tm_hour * 60 + timeTuple->tm_min) * 60 + timeTuple->tm_sec;
    return (float) secondOfDay / (24.0f * 60.0f * 60.0f);
}


/// creates timestamp corresponding to date (midnight)
int dateUTC( int year, int month, int day ) {
    if (year < 1980 || year > 2037 || month < 1 || month > 12 || day < 1 || day > 31) {
        warning( "invalid date" );
        return 0;
    }
    struct tm timeTuple;
    timeTuple.tm_year = year - 1900;
    timeTuple.tm_mon = month - 1;
    timeTuple.tm_mday = day;
    timeTuple.tm_hour = 0;
    timeTuple.tm_min = 0;
    timeTuple.tm_sec = 0;
    timeTuple.tm_isdst = 0;
    return (int) timegm( &timeTuple );
}


/// returns timestamp corresponding to start of day (UTC) 
int dateUTC( const String &dateStr ) {
    int date = 0;
    Array<String> split = dateStr.split( "-" );
    if (split.count() == 3) 
        date = dateUTC( split[ 0 ].toInt(), split[ 1 ].toInt(), split[ 2 ].toInt() );
    return date;
}


/// returns timestamp corresponding to start of day (UTC) 
int datePartUTC( int timestamp ) {
    time_t tsLong = timestamp;
    struct tm *timeTuple = gmtime( &tsLong );    
    return dateUTC( timeTuple->tm_year + 1900, timeTuple->tm_mon + 1, timeTuple->tm_mday );
}


/// creates timestamp from date and time
int dateTimeLocal( int year, int month, int day, int hour, int minute, int second ) {
    if (year < 1980 || year > 2037 || month < 1 || month > 12 || day < 1 || day > 31) {
        warning( "invalid date" );
        return 0;
    }
    if (hour < 0 || hour > 23 || minute < 0 || minute > 59 || second < 0 || second > 59) {
        warning( "invalid time (hour: %d, minute: %d, second: %d)", hour, minute, second );
        return 0;
    }
    struct tm timeTuple;
    timeTuple.tm_year = year - 1900;
    timeTuple.tm_mon = month - 1;
    timeTuple.tm_mday = day;
    timeTuple.tm_hour = hour;
    timeTuple.tm_min = minute;
    timeTuple.tm_sec = second;
    timeTuple.tm_isdst = -1;
    return (int) mktime( &timeTuple );
}


/// creates timestamp from date and time
int dateTimeUTC( int year, int month, int day, int hour, int minute, int second ) {
    if (year < 1980 || year > 2037 || month < 1 || month > 12 || day < 1 || day > 31) {
        warning( "invalid date" );
        return 0;
    }
    if (hour < 0 || hour > 23 || minute < 0 || minute > 59 || second < 0 || second > 59) {
        warning( "invalid time (hour: %d, minute: %d, second: %d)", hour, minute, second );
        return 0;
    }
    struct tm timeTuple;
    timeTuple.tm_year = year - 1900;
    timeTuple.tm_mon = month - 1;
    timeTuple.tm_mday = day;
    timeTuple.tm_hour = hour;
    timeTuple.tm_min = minute;
    timeTuple.tm_sec = second;
    timeTuple.tm_isdst = 0;
    return (int) timegm( &timeTuple );
}


/// parse a date-time string into a timestamp, assuming YMD-HMS order
int parseTimestamp( const String &data, const String &format, bool local, int &msec ) {
    msec = 0;
    int year = 0, month = 0, day = 0, hour = 0, minute = 0, second = 0;
    // fix(safety): check number of positions in format    
    sscanf( data.c_str(), format.c_str(), &year, &month, &day, &hour, &minute, &second, &msec );
    if (local)
        return dateTimeLocal( year, month, day, hour, minute, second );
    else
        return dateTimeUTC( year, month, day, hour, minute, second );
}


/// returns the current timestamp
int now() {
    return (int) time( NULL );
}


// basic constructors
Date::Date( int year, int month, int day ) { 
    m_year = year; 
    m_month = month; 
    m_day = day; 
    m_startTimestamp = dateUTC( year, month, day );
}


// basic constructor
Date::Date( int timestamp ) {
    time_t tsLong = timestamp;
    struct tm *timeTuple = gmtime( &tsLong );    
    m_year = timeTuple->tm_year + 1900;
    m_month = timeTuple->tm_mon + 1;
    m_day = timeTuple->tm_mday;
    m_startTimestamp = dateUTC( m_year, m_month, m_day );
}


/// move date forward by one day
void Date::increment() {
    m_startTimestamp += 24 * 60 * 60; 
    time_t tsLong = m_startTimestamp;
    struct tm *timeTuple = gmtime( &tsLong );    
    m_year = timeTuple->tm_year + 1900;
    m_month = timeTuple->tm_mon + 1;
    m_day = timeTuple->tm_mday;
}


} // end namespace sbl

