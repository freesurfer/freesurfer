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
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.4 $
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
