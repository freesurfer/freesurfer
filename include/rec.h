/**
 * @file  rec.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.5 $
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


#ifndef REC_H
#define REC_H

#include "matrix.h"

typedef struct
{
  int    ntimepts ;
  int    ptime ;
  int    neeg_channels ;
  int    nmeg_channels ;
  float  *latencies ;
  MATRIX *m_data ;
}
REC ;


/* each column is a time course (except 1st)*/
REC *RecRead(char *fname, int iop_neeg, int iop_nmeg) ;

/* flag == 1 - load only EEG channels
   flag == 2 - load only MEG channels
*/
REC *RecReadPartially(char *fname, int iop_neeg, int iop_nmeg, int flag) ;
#endif
