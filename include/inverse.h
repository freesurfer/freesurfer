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
} INVERSE_OPERATOR, IOP ;

IOP *IOPRead(char *fname, int hemi) ;
int IOPWrite(IOP *iop, char *fname) ;
int IOPFree(IOP **piop) ;
int IOPNormalize(IOP *iop) ;
MATRIX *IOPapplyInverseOperator(IOP *iop, REC *rec, MATRIX *m_sol) ;


#define INVERSE_RH    1
#define INVERSE_LH    2

#endif
