/**
 * @file  timeb.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:02 $
 *    $Revision: 1.2 $
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


/* @(#)timeb.h 2.6 88/08/19 SMI; from UCB 4.2 81/02/19 */

/*
 * Structure returned by ftime system call
 */

#ifndef _sys_timeb_h
#define _sys_timeb_h

#include <sys/time.h>

struct timeb
{
  time_t time;
  unsigned short millitm;
  short timezone;
  short dstflag;
};

#if defined(IRIX)
int ftime(struct timeb *t) ;
#endif

#endif /*!_sys_timeb_h*/
