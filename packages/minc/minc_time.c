/*
 * Original Author: David MacDonald, modified to compile within freesurfer/utils by Bevin Brett
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include "minc_internals.h"
#include <math.h>
#include <string.h>
#include <unistd.h>

#define VIOAPI

#include  <sys/types.h>

#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif
#if HAVE_UNISTD_H
#include  <unistd.h>
#endif

#ifndef CLK_TCK
#define CLK_TCK CLOCKS_PER_SEC
#endif

// #include "timer.h"

#ifndef lint
//static char rcsid[] = "$Header: /private-cvsroot/minc/volume_io/Prog_utils/time.c,v 1.21.2.3 2005/07/13 19:59:19 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_clock_ticks_per_second
@INPUT      : 
@OUTPUT     : 
@RETURNS    : number clock ticks per second
@DESCRIPTION: Returns the number of clock ticks per second in a system
              independent fashion
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jul 3, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static double  get_clock_ticks_per_second( )
{
    static volatile bool    initialized = false;
    static volatile double  clock_ticks_per_second;

    if( !initialized )
    {
#if HAVE_SYSCONF
        clock_ticks_per_second = (double) sysconf( _SC_CLK_TCK );
#else
        clock_ticks_per_second = (double) CLK_TCK;
#endif
        initialized = true;
    }

    return( clock_ticks_per_second );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : current_cpu_seconds
@INPUT      : 
@OUTPUT     : 
@RETURNS    : # seconds
@DESCRIPTION: Returns the number of cpu seconds used by the program to date.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

double  current_cpu_seconds( )
{
    static volatile bool    first_call = true;
    static volatile clock_t first;
    
    double secs;

    if (first_call)
    {
        first = clock();
        secs = (double) first / get_clock_ticks_per_second();
        first_call = false;
    }
    else
    {
    	clock_t current = clock();
        secs = (double) (current - first) / get_clock_ticks_per_second();
    }
    return (secs);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : current_realtime_seconds
@INPUT      : 
@OUTPUT     : 
@RETURNS    : # seconds
@DESCRIPTION: Returns the number of seconds since the first invocation of this
            : function.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

double  current_realtime_seconds( )
{
    static volatile bool first_call = true;
    static volatile time_t first;
    double secs;

    if( first_call )
    {
        first = time(NULL);
        secs = 0.0;
        first_call = false;
    }
    else
    {
        time_t current = time(NULL);
        secs = (double) (current - first);
    }

    return( secs );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : format_time
@INPUT      : format
            : seconds
@OUTPUT     : str
@RETURNS    : 
@DESCRIPTION: Decides what time unit to use and displays the seconds value
            : in str, using format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

char* format_time(
    const char* format,
    double      seconds )
{
 #define units_size 7
    static   char   *units [units_size] = { "us",   "ms",   "sec", "min", "hrs", "days", "years" };
    static   double  scales[units_size] = { 1000.0, 1000.0, 60.0,  60.0,  24.0,  365.0,  0.0     };
    char     buffer[1024];

    bool negative = seconds < 0.0;
    if( negative )  seconds = -seconds;

    seconds *= 1.0e6;

    int      i;
    for( i = 0; i < units_size - 1; i++ )
    {
        if( seconds > 2.0 * scales[i] )
        {
            seconds /= scales[i];
        }
        else
        {
            break;
        }
    }

    seconds = (double) round( 10.0 * seconds ) / 10.0;

    if( negative )  seconds = -seconds;

    (void) sprintf( buffer, format, seconds, units[i] );

    return( strdup( buffer ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : print_time
@INPUT      : format
            : seconds
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Prints out the time in suitable units.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void  print_time(
    const char* format,
    double      seconds )
{
    char* str = format_time( format, seconds );

    printf( "%s", str );

    free( str );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_clock_time
@INPUT      : 
@OUTPUT     : time_str
@RETURNS    : 
@DESCRIPTION: Stores the current time of day in the "time_str".
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

char* get_clock_time( void )
{
    time_t           clock_time;
    struct  tm       *time_tm;
    char             *str;

    (void) time( &clock_time );

    time_tm = localtime( &clock_time );

    str = asctime( time_tm );

    return( strdup( str ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : sleep_program
@INPUT      : seconds
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Make the program sleep for the specified number of seconds.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void sleep_program( double seconds )
{
#if HAVE_SELECT
    struct  timeval  timeout;

    timeout.tv_sec = (long) seconds;
    timeout.tv_usec = (long) (1.0e6 * (seconds - (double) timeout.tv_sec) + 0.5);

    (void) select( 0, NULL, NULL, NULL, &timeout );
#else
    sleep((unsigned int) seconds);
#endif
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_date
@INPUT      : 
@OUTPUT     : date_str
@RETURNS    : 
@DESCRIPTION: Fills in the date into the string.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

const char* get_date() 
{
    time_t           clock_time;
    struct  tm       *time_tm;
    char             *str;

    (void) time( &clock_time );

    time_tm = localtime( &clock_time );

    str = asctime( time_tm );

    return( strdup( str ) );
}
