/**
 * @file  runfuncs.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
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


/*
  @(#)runfuncs.h  1.1
  4/4/94
*/
/*------------------------------------------------------------------------
      File Name:  runfuncs.h

         Author:  Bruce Fischl

        Created:  Jan. 1993

    Description:

------------------------------------------------------------------------*/
#ifndef RUNFUNCS_H
#define RUNFUNCS_H

typedef int (*fwd_func)(unsigned char **outPtr) ;
typedef void (*inv_func)(unsigned char **outPtr, unsigned char val) ;

void runFuncInit(fwd_func *fwd_array, inv_func *inv_array) ;


#endif
