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
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.3 $
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
