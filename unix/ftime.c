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
 *    $Date: 2006/12/29 02:09:18 $
 *    $Revision: 1.6 $
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
