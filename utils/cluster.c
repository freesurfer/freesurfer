/*
 *       FILE NAME:   cluster.c
 *
 *       DESCRIPTION: clustering utilities
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        5/16/97
 *
*/

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <math.h>

#include "matrix.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "proto.h"
#include "cluster.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static int   clusterInit(CLUSTER *cluster, int ninputs) ;
static int   clusterFree(CLUSTER *cluster) ;
static int   clusterNewObservation(CLUSTER *cluster, VECTOR *v_obs) ;
static float clusterDistance(CLUSTER *cluster, VECTOR *v_obs) ;
static int   clusterComputeStatistics(CLUSTER *cluster) ;
static int   clusterPrint(CLUSTER *cluster, FILE *fp) ;
static float clusterVariance(CLUSTER *cluster) ;
static CLUSTER *clusterDivide(CLUSTER *csrc, CLUSTER *cdst) ;
static int   normalizeObservation(CLUSTER_SET *cs, VECTOR *v_obs) ;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Print the contents of a CLUSTER_SET to a file.
------------------------------------------------------*/
int
CSprint(CLUSTER_SET *cs, FILE *fp)
{
  int  c ;

  fprintf(fp, "cluster set %swith %d inputs, %d obs, and %d clusters, "
          "looking for %d\n",
          cs->normalize ? "normalized " : "", cs->ninputs, cs->nsamples, 
          cs->nclusters, cs->max_clusters)  ;
#if 0
  fprintf(fp, "within-cluster scatter matrix:\n") ;
  MatrixPrint(fp, cs->m_sw) ;
  fprintf(fp, "between-cluster scatter matrix:\n") ;
  MatrixPrint(fp, cs->m_sb) ;
#endif
  for (c = 0 ; c < cs->nclusters ; c++)
    clusterPrint(cs->clusters+c, fp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           allocate and initialize a CLUSTER_SET structure
           in preparation for clustering a data set.
------------------------------------------------------*/
CLUSTER_SET *
CSinit(int max_clusters, int ninputs, int normalize)
{
  CLUSTER_SET  *cs ;
  int          c ;

  cs = (CLUSTER_SET *)calloc(1, sizeof(CLUSTER_SET)) ;
  if (!cs)
    ErrorExit(ERROR_NO_MEMORY,"ClusterInit(%d, %d); could not allocate struct",
              max_clusters, ninputs) ;

  cs->normalize = normalize ;
  cs->ninputs = ninputs ;
  cs->nclusters = 1 ;
  cs->max_clusters = max_clusters ;
  cs->clusters = (CLUSTER *)calloc(max_clusters, sizeof(CLUSTER)) ;
  if (!cs->clusters)
    ErrorExit(ERROR_NO_MEMORY,
              "ClusterInit(%d, %d); could not allocate clusters",
              max_clusters, ninputs) ;

  for (c = 0 ; c < max_clusters ; c++)
  {
    cs->clusters[c].cno = c ;
    clusterInit(cs->clusters+c, ninputs) ;
  }

  cs->m_sw = MatrixAlloc(ninputs, ninputs, MATRIX_REAL) ;
  if (!cs->m_sw)
    ErrorExit(ERROR_NO_MEMORY,
              "ClusterInit(%d, %d); could not allocate Sw",
              max_clusters, ninputs) ;
  cs->m_sb = MatrixAlloc(ninputs, ninputs, MATRIX_REAL) ;
  if (!cs->m_sb)
    ErrorExit(ERROR_NO_MEMORY,
              "ClusterInit(%d, %d); could not allocate Sb",
              max_clusters, ninputs) ;
  cs->nobs = 0 ;
  cs->means = (float *)calloc(ninputs, sizeof(float)) ;
  if (!cs->means)
    ErrorExit(ERROR_NO_MEMORY,
              "ClusterInit(%d, %d); could not allocate means",
              max_clusters, ninputs) ;
  cs->stds =  (float *)calloc(ninputs, sizeof(float)) ;
  if (!cs->stds)
    ErrorExit(ERROR_NO_MEMORY,
              "ClusterInit(%d, %d); could not allocate std",
              max_clusters, ninputs) ;

  for (c = 0 ; c < ninputs ; c++)
    cs->stds[c] = 1.0f ;

  return(cs) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Free a CLUSTER_SET structure and all memory that
          it has allocated.
------------------------------------------------------*/
int
CSfree(CLUSTER_SET **pcs)
{
  CLUSTER_SET  *cs ;
  int          c ;

  cs = *pcs ;
  *pcs = NULL ;

  if (cs->means)
    free(cs->means) ;
  if (cs->stds)
    free(cs->stds) ;

  if (cs->m_sw)
    MatrixFree(&cs->m_sw) ;
  if (cs->m_sb)
    MatrixFree(&cs->m_sb) ;
  
  for (c = 0 ; c < cs->max_clusters ; c++)
    clusterFree(cs->clusters+c) ;
  free(cs) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           allocate the internal variables of the given cluster
           and initialize them to reasonable values.
------------------------------------------------------*/
static int 
clusterInit(CLUSTER *cluster, int ninputs)
{
  cluster->m_scatter = MatrixAlloc(ninputs, ninputs, MATRIX_REAL) ;
  if (!cluster->m_scatter)
    ErrorExit(ERROR_NO_MEMORY, 
              "clusterInit(%d): could not allocate scatter matrix", ninputs) ;
  cluster->v_means = VectorAlloc(ninputs, MATRIX_REAL) ;
  if (!cluster->v_means)
    ErrorExit(ERROR_NO_MEMORY, 
              "clusterInit(%d): could not allocate mean vector", ninputs) ;
  cluster->v_seed = VectorAlloc(ninputs, MATRIX_REAL) ;
  if (!cluster->v_seed)
    ErrorExit(ERROR_NO_MEMORY, 
              "clusterInit(%d): could not allocate centroid vector", ninputs) ;
  cluster->evalues = (float *)calloc(ninputs+1, sizeof(float)) ;
  if (!cluster->evalues)
    ErrorExit(ERROR_NO_MEMORY, 
              "clusterInit(%d): could not allocate eigen values", ninputs) ;
  cluster->nobs = cluster->nsamples = 0 ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Free a CLUSTER structure and all memory that
          it has allocated.
------------------------------------------------------*/
static int
clusterFree(CLUSTER *cluster)
{
  if (cluster->m_scatter)
    MatrixFree(&cluster->m_scatter) ;
  if (cluster->v_means)
    VectorFree(&cluster->v_means) ;
  if (cluster->evalues)
    free(cluster->evalues) ;
  if (cluster->v_seed)
    VectorFree(&cluster->v_seed) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Add this observation to the appropriate cluster.
------------------------------------------------------*/
int
CSnewObservation(CLUSTER_SET *cs, VECTOR *v_obs)
{
  int    c, min_cluster ;
  float  min_dist, dist ;

  if (cs->normalize)
    normalizeObservation(cs, v_obs) ;

  min_dist = clusterDistance(cs->clusters+0, v_obs) ;
  min_cluster = 0 ;
  cs->nobs++ ;

  for (c = 1 ; c < cs->nclusters ; c++)
  {
    dist = clusterDistance(cs->clusters+c, v_obs) ;
    if (dist < min_dist)
    {
      min_dist = dist ;
      min_cluster = c ;
    }
  }
  clusterNewObservation(cs->clusters+min_cluster, v_obs) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           A pass has just been made through the data set, divide
           the cluster with the maximum variance.
------------------------------------------------------*/
int
CSdivide(CLUSTER_SET *cs)
{
  int    c, max_cluster ;
  float  max_variance, variance ;

  if (cs->nclusters >= cs->max_clusters)
    return(CS_DONE) ;

  max_cluster = 0 ;
  max_variance = clusterVariance(cs->clusters+0) ;
  for (c = 1 ; c < cs->nclusters ; c++)
  {
    variance = clusterVariance(cs->clusters+c) ;
    if (variance > max_variance)
    {
      max_variance = variance ;
      max_cluster = c ;
    }
  }
  fprintf(stderr, "dividing cluster %d\n", max_cluster) ;
  clusterDivide(cs->clusters+max_cluster, cs->clusters+cs->nclusters++) ;
  return(cs->nclusters < cs->max_clusters ? CS_CONTINUE : CS_DONE) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Add this observation to the appropriate cluster.
------------------------------------------------------*/
static int
clusterNewObservation(CLUSTER *cluster, VECTOR *v_obs)
{
  float  covariance ;
  int    row, col ;

  cluster->nobs++ ;
  for (row = 1 ; row <= v_obs->rows ; row++)
    VECTOR_ELT(cluster->v_means, row) += VECTOR_ELT(v_obs, row) ;

  for (row = 1 ; row <= cluster->m_scatter->rows ; row++)
  {
    for (col = 1 ; col <= row ; col++)
    {
      covariance = cluster->m_scatter->rptr[row][col] + 
        VECTOR_ELT(v_obs,row) * VECTOR_ELT(v_obs, col) ;
      cluster->m_scatter->rptr[row][col] = covariance;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Determine the distance from this cluster to the observation
          vector.
------------------------------------------------------*/
static float
clusterDistance(CLUSTER *cluster, VECTOR *v_obs)
{
  float  dist ;

#if 0
  dist = VectorNormalizedDot(cluster->v_seed, v_obs) ;
#else
  dist = VectorDistance(cluster->v_seed, v_obs) ;
#endif
  return(dist) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          a round of observations is done - complete the calculation
          of means and scatter matrices for this cluster
------------------------------------------------------*/
static int
clusterComputeStatistics(CLUSTER *cluster)
{
  float   mean_a, mean_b, covariance ;
  int     nobs, row, col ;

  nobs = cluster->nobs ;
  if (nobs)
  {
    /* compute means */
    for (row = 1 ; row <= cluster->v_means->rows ; row++)
      VECTOR_ELT(cluster->v_means, row) /= (float)nobs ;

    for (row = 1 ; row <= cluster->v_means->rows ; row++)
    {
      mean_a = VECTOR_ELT(cluster->v_means, row) ;
      for (col = 1 ; col <= row ; col++)
      {
        mean_b = cluster->v_means->rptr[col][1] ;
        covariance = cluster->m_scatter->rptr[row][col] / nobs - mean_a*mean_b;
        cluster->m_scatter->rptr[row][col] = covariance ;
        cluster->m_scatter->rptr[col][row] = covariance ;
      }
    }
    cluster->m_inverse = MatrixInverse(cluster->m_scatter, cluster->m_inverse) ;
#if 0
    fprintf(stderr, "scatter matrix:\n") ;
    MatrixPrint(stderr, cluster->m_scatter) ;
    fprintf(stderr, "inverse scatter matrix:\n") ;
    MatrixPrint(stderr, cluster->m_inverse) ;
#endif
  }

  cluster->m_evectors = 
    MatrixEigenSystem(cluster->m_scatter, cluster->evalues, 
                      cluster->m_evectors) ;
  VectorCopy(cluster->v_means, cluster->v_seed) ;
  cluster->nsamples = cluster->nobs ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Print the contents of a cluster to a file.
------------------------------------------------------*/
static int
clusterPrint(CLUSTER *cluster, FILE *fp)
{
  /*  int    i ;*/

  fprintf(fp, "cluster %d has %d observations. Seed:", 
          cluster->cno, cluster->nsamples) ;
  MatrixPrintTranspose(fp, cluster->v_seed) ;
#if 0
  fprintf(fp, "scatter matrix:\n") ;
  MatrixPrint(fp, cluster->m_scatter) ;
  if (cluster->m_evectors)
  {
    fprintf(fp, "eigen vectors:\n") ;
    MatrixPrint(fp, cluster->m_evectors) ;
    fprintf(fp, "eigen values: ") ;
    for (i = 0 ; i < cluster->m_evectors->rows ; i++)
      fprintf(fp, "%2.3f  ", cluster->evalues[i]) ;
    fprintf(fp, "\n") ;
  }
#endif
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the variance of this cluster to decide whether (and how)
          it should be split.
------------------------------------------------------*/
static float
clusterVariance(CLUSTER *cluster)
{
  float  variance ;
  int    i ;

  for (variance = 0.0f, i = 0 ; i < cluster->m_scatter->rows ; i++)
    variance += cluster->evalues[i] ;

  return(i) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Divide this cluster in two by offsetting the mean vectors in
          the direction of the largest eigenvector.
------------------------------------------------------*/
static CLUSTER *
clusterDivide(CLUSTER *csrc, CLUSTER *cdst)
{
  VECTOR  *v_e, *v_mT ;
  float   len ;

  fprintf(stderr, "initial cluster centroid at: ") ;
  MatrixPrintTranspose(stderr, csrc->v_means) ;

  v_e = MatrixColumn(csrc->m_evectors, NULL, 1) ; /* extract eigenvector */

  fprintf(stderr, "evector: ") ;
  MatrixPrintTranspose(stderr, v_e) ;

  /*
    split the cluster by moving the two new seed values in the direction
    parallel to and antiparallel to the eigenvector associated with the
    largest eigenvalue of the scatter matrix. By moving the seed point only
    a little, we keep the other clusters relatively intact.
    */
  len = sqrt(csrc->evalues[0]) ;
  VectorScalarMul(v_e, .1f*len, v_e) ;
  VectorAdd(csrc->v_means, v_e, cdst->v_seed) ;
  VectorScalarMul(v_e, -1.0f, v_e) ;
  VectorAdd(csrc->v_means, v_e, csrc->v_seed) ;
  VectorFree(&v_e) ;

  fprintf(stderr, "cluster %d split into two:\n", csrc->cno) ;
  fprintf(stderr, "mean 1: ") ;
  v_mT = MatrixTranspose(csrc->v_seed, NULL) ;
  MatrixPrint(stderr, v_mT) ;
  MatrixFree(&v_mT) ;
  fprintf(stderr, "mean 2: ") ;
  v_mT = MatrixTranspose(cdst->v_seed, NULL) ;
  MatrixPrint(stderr, v_mT) ;
  MatrixFree(&v_mT) ;
  
#if 0
  clusterPrint(csrc, stderr) ;
  clusterPrint(cdst, stderr) ;
#endif
  return(cdst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Perform the clustering. This takes place in two stages:

          first, start with a single cluster and divide it based on
          the direction of maximum variance (eigenvector associated
          with the largest eigenvalue). After the initial clusters
          have been determined, then use an optimal procedure based
          on maxizing the eigenvectors of Sw-1 Sb, where Sw-1 is the
          inverse of the within-cluster scatter matrix, and Sb is the
          between cluster scatter matrix.
------------------------------------------------------*/
int
CScluster(CLUSTER_SET *cs, int (*get_observation_func)
          (VECTOR *v_obs, int no, void *parm), void *parm)
{
  int         obs_no ;
  VECTOR      *v_obs ;

  v_obs = VectorAlloc(cs->ninputs, MATRIX_REAL) ;

  if (cs->normalize)
    CScomputeDimensionStatistics(cs, get_observation_func, parm) ;

  do
  {
    CSreset(cs) ;
    obs_no = 0 ;
    while ((*get_observation_func)(v_obs, obs_no++, parm) == NO_ERROR)
      CSnewObservation(cs, v_obs) ;
    CScomputeStatistics(cs) ;
  } while (CSdivide(cs) == CS_CONTINUE) ;

  /* Go through the data one more time to rebuild scatter matrices and means 
     with final # of clusters */
  CSreset(cs) ;
  obs_no = 0 ;
  while ((*get_observation_func)(v_obs, obs_no++, parm) == NO_ERROR)
    CSnewObservation(cs, v_obs) ;
  CScomputeStatistics(cs) ;
  if (cs->normalize)
    CSrenormalize(cs) ;

  /* now we have the proper # of clusters, use a scale-invariant metric
     to find the optimum grouping
     */

  VectorFree(&v_obs) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Go through all the clusters and clear their mean and
          scatter matrices.
------------------------------------------------------*/
int
CSreset(CLUSTER_SET *cs)
{
  int      c ;
  CLUSTER  *cluster ;
  
  for (c = 0 ; c < cs->nclusters ; c++)
  {
    cluster = cs->clusters+c ;
    cluster->nobs = 0 ;
    MatrixClear(cluster->m_scatter) ;
    VectorClear(cluster->v_means) ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the covariances and means of all the clusters 
          (performed after each complete pass through the data).
------------------------------------------------------*/
int
CScomputeStatistics(CLUSTER_SET *cs)
{
  int    c ;

  for (c = 0 ; c < cs->nclusters ; c++)
    clusterComputeStatistics(cs->clusters+c) ;

  cs->nsamples = cs->nobs ;
  cs->nobs = 0 ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the means and standard deviations of the input
          dimensions so that they can be normalized.
------------------------------------------------------*/
int
CScomputeDimensionStatistics(CLUSTER_SET *cs, int (*get_observation_func)
          (VECTOR *v_obs, int no, void *parm), void *parm)
{
  int         obs_no = 0, i ;
  VECTOR      *v_obs ;
  float       v, mean ;

  v_obs = VectorAlloc(cs->ninputs, MATRIX_REAL) ;
  memset(cs->means, 0, cs->ninputs*sizeof(float)) ;
  memset(cs->stds, 0, cs->ninputs*sizeof(float)) ;
  while ((*get_observation_func)(v_obs, obs_no, parm) == NO_ERROR)
  {
    obs_no++ ;
    for (i = 0 ; i < cs->ninputs ; i++)
    {
      v = VECTOR_ELT(v_obs, i+1) ;
      cs->means[i] += v ;
      cs->stds[i] += (v*v) ;
    }
  }

  for (i = 0 ; i < cs->ninputs ; i++)
  {
    mean = cs->means[i] ;
    mean /= (obs_no+1) ;
    cs->means[i] = mean ;
    cs->stds[i] = sqrt(cs->stds[i] / obs_no - mean*mean) ;
  }
  
  VectorFree(&v_obs) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Normalize the observation using the mean and standard
          deviation of each dimension.
------------------------------------------------------*/
static int
normalizeObservation(CLUSTER_SET *cs, VECTOR *v_obs)
{
  int  i ;
  float  mean, std ;

  for (i = 1 ; i <= v_obs->rows ; i++)
  {
    mean = cs->means[i-1] ; std = cs->stds[i-1] ;
    VECTOR_ELT(v_obs,i) = (VECTOR_ELT(v_obs,i) - mean) / std ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Undo the effects of input normalization so that the
          clusters are represented in terms of the original
          measurement spaces.
------------------------------------------------------*/
int
CSrenormalize(CLUSTER_SET *cs)
{
  int      c, row ;
  CLUSTER  *cluster ;
  float    mean, sigma ;

  for (c = 0 ; c < cs->nclusters ; c++)
  {
    cluster = cs->clusters+c ;
    for (row = 1 ; row <= cs->ninputs ; row++)
    {
      mean = cs->means[row-1] ;   /* zero-based */
      sigma = cs->stds[row-1] ;   /* zero-based */
      VECTOR_ELT(cluster->v_seed, row) *= sigma ;
      VECTOR_ELT(cluster->v_seed, row) += mean ;
      VECTOR_ELT(cluster->v_means, row) *= sigma ;
      VECTOR_ELT(cluster->v_means, row) += mean ;
    }
  }
  return(NO_ERROR) ;
}

