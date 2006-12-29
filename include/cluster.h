/**
 * @file  cluster.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.6 $
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
