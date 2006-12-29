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
 *    $Date: 2006/12/29 02:09:00 $
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
