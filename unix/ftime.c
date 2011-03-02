/**
 * @file  ftime.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:41 $
 *    $Revision: 1.7 $
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


#ifdef IRIX

#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/timeb.h>
#include <sys/param.h>

int
ftime(struct timeb *t)
{
  struct tms tm ;
  struct timeval tv;
  struct timezone tz;
  int    msec ;
  int    sec ;

  tz.tz_minuteswest = 0;
  tz.tz_dsttime = 0;
  gettimeofday(&tv,&tz);
  msec = (int)((double)tv.tv_usec * 0.001 + 0.5) ;
  sec = (int)((double)tv.tv_sec+0.5) ;
  t->millitm = msec ;
  t->time = sec ;
  return(time(NULL)) ;
}

#endif

#ifdef Darwin

#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/timeb.h>
#include <sys/param.h>

int
ftime(struct timeb *t)
{
  struct timeval tv;
  struct timezone tz;
  int    msec ;
  int    sec ;

  tz.tz_minuteswest = 0;
  tz.tz_dsttime = 0;
  gettimeofday(&tv,&tz);
  msec = (int)((double)tv.tv_usec * 0.001 + 0.5) ;
  sec = (int)((double)tv.tv_sec+0.5) ;
  t->millitm = msec ;
  t->time = sec ;
  return(time(NULL)) ;
}

#endif
