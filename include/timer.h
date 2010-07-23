/**
 * @file  timer.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2010/07/23 21:07:45 $
 *    $Revision: 1.4 $
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


#ifndef TIMER_H
#define TIMER_H

#if defined(__cplusplus)
extern "C" {
#endif

#include <sys/timeb.h>

struct timeb *TimerStart(struct timeb *then) ;
int TimerStop(struct timeb *then) ;

#ifdef Linux
/* don't know why this doesn't work on linux, but.... */
extern int ftime (struct timeb *__timebuf);
#endif

#if defined(__cplusplus)
};
#endif

#endif
