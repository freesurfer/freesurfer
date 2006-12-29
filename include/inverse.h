/**
 * @file  inverse.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.3 $
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


#ifndef INVERSE_H
#define INVERSE_H

#include "matrix.h"
#include "rec.h"

typedef struct
{
  int    version ;
  int    neeg_channels ;
  int    nmeg_channels ;
  int    nchan ;                 /* sum of the previous two */
  int    ndipoles_per_location ; /* 1 for orientation, 3 otherwise */
  int    ndipoles ;              /* # of dipole locations */
  int    ndipole_files ;
  float  pthresh ;               /*threshold that was used to generate priors*/
  int    *dipole_vertices ;      /* vertex numbers where dipoles are situated*/
  VECTOR *spatial_priors ;       /* fMRI at the moment */
  int    *bad_sensors ;
  VECTOR *pvals ;                /* statistics from fMRI */
  float  *dipole_normalization ;
  MATRIX *m_iop ;                /* actual inverse operator */
  MATRIX *m_forward ;            /* forward solution */
}
INVERSE_OPERATOR, IOP ;

IOP *IOPRead(char *fname, int hemi) ;
int IOPWrite(IOP *iop, char *fname) ;
int IOPFree(IOP **piop) ;
int IOPNormalize(IOP *iop) ;
MATRIX *IOPapplyInverseOperator(IOP *iop, REC *rec, MATRIX *m_sol) ;


#define INVERSE_RH    1
#define INVERSE_LH    2

#endif
