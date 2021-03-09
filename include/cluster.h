/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef CLUSTER_H
#define CLUSTER_H

#include "matrix.h"

typedef struct
{
  MATRIX  *m_scatter ;      /* within-cluster scatter matrix */
  MATRIX  *m_inverse ;      /* inverse of scatter (covariance) matrix */
  VECTOR  *v_means ;        /* vector of means */
  VECTOR  *v_seed ;         /* used for distance-to-cluster calculations */
  MATRIX  *m_evectors ;     /* columns are eigenvectors of scatter matrix */
  float   *evalues ;        /* eigenvalues of scatter matrix */
  int     nobs ;            /* # of observation in this cluster */
  int     nsamples ;        /* # of samples in this cluster */
  int     cno ;             /* index of this cluster (for diagnostics) */
  float   det ;             /* determinant of scatter matrix */
  float   norm ;            /* sqrt(det(scatter)) for Gaussian */
  int     ill_conditioned;  /* was the scatter matrix singular */
}
CLUSTER ;

typedef struct
{
  int     nclusters ;       /* current # of clusters */
  int     ninputs ;         /* dimensionality of the space */
  int     max_clusters ;    /* # of desired clusters */
  CLUSTER *clusters ;       /* array of clusters */
  MATRIX  *m_sw ;           /* within-cluster scatter matrix */
  MATRIX  *m_sb ;           /* between-cluster scatter matrix */
  int     nobs ;            /* # of observations so far */
  int     nsamples ;        /* total # of samples */
  int     normalize ;       /* normalize samples if true */
  float   *means ;          /* means each input dimension */
  float   *stds ;           /* sigma of the input dimensions */
}
CLUSTER_SET ;

CLUSTER_SET  *CSinit(int max_clusters, int ninputs, int normalize) ;
int          CSfree(CLUSTER_SET **pcs) ;
int          CSnewObservation(CLUSTER_SET *cs, VECTOR *v_obs) ;
int          CSdivide(CLUSTER_SET *cs) ;
int          CSprint(CLUSTER_SET *cs, FILE *fp) ;
int          CScluster(CLUSTER_SET *cs,
                       int (*get_observation_func)
                       (VECTOR *v_obs, int no, void *parm), void *parm) ;
int          CSreset(CLUSTER_SET *cs) ;
int          CScomputeStatistics(CLUSTER_SET *cs) ;
int          CScomputeDimensionStatistics(CLUSTER_SET *cs,
    int (*get_observation_func)
    (VECTOR *v_obs, int no, void *parm),
    void *parm) ;
int          CSrenormalize(CLUSTER_SET *cs) ;
int          CSwriteInto(FILE *fp, CLUSTER_SET *cs) ;
CLUSTER_SET  *CSreadFrom(FILE *fp, CLUSTER_SET *cs) ;

/* return values for CSdivide. Negative ones indicate error */
#define CS_CONTINUE              0
#define CS_DONE                  1

/* values for the normalize parm of CSinit */
#define CLUSTER_DONT_NORMALIZE   0
#define CLUSTER_NORMALIZE        1

#endif
