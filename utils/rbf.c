/*
 *       FILE NAME:   rbf.c
 *
 *       DESCRIPTION: Radial Basis Function classification utilities.
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        5/19/97
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

#include "proto.h"
#include "error.h"
#include "diag.h"
#include "cluster.h"
#include "rbf.h"
#include "macros.h"
#include "utils.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define TRATE_DECREASE  0.7f
#define TRATE_INCREASE  1.05f
#define ERROR_RATIO     1.04f
#define MIN_TRATE       0.00001f
#define DEFAULT_TRATE   0.05f
#define MOMENTUM        0.8f
#define MAX_EPOCHS      4000   /* maximum # of training iterations */

/* 
   if the sse doesn't change by more than MIN_DELTA_SSE for MAX_SMALL
   epochs in a row, assume that training has asymptoted.
*/
#define MIN_DELTA_SSE   MIN_TRATE/100.0f
#define MAX_SMALL       10

#define TRAIN_OUTPUTS   0x001
#define TRAIN_SPREADS   0x002
#define TRAIN_CENTERS   0x004
#define TRAIN_HIDDEN    (TRAIN_SPREADS|TRAIN_CENTERS)
#define TRAIN_ALL       (TRAIN_OUTPUTS|TRAIN_SPREADS|TRAIN_CENTERS)

/*-----------------------------------------------------
                      STRUCTURES
-------------------------------------------------------*/

typedef struct
{
  int  current_class ;
  int  (*get_observation_func)(VECTOR *v, int obs_no, void *parm,int *classno);
  void *parm ;
} CLUSTERING_PARM ;

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static int   rbf_get_obs_func(VECTOR *v_obs, int obs_no, void *vrbf) ;
static int   rbfGradientDescent(RBF *rbf, int (*get_observation_func)
               (VECTOR *v_obs, int no, void *parm, int *pclass), void *parm) ;
static int   rbfComputeHiddenActivations(RBF *rbf, VECTOR *v_obs) ;
static float rbfGaussian(MATRIX *m_sigma_inverse,VECTOR *v_means,
                         VECTOR *v_obs, VECTOR *v_z);
static int   rbfShowClusterCenters(RBF *rbf, FILE *fp) ;
static float rbfComputeErrors(RBF *rbf, int class, VECTOR *v_error) ;
static int   rbfAdjustOutputWeights(RBF *rbf, VECTOR *v_error) ;
static int   rbfAdjustHiddenCenters(RBF *rbf, VECTOR *v_error) ;
static int   rbfAdjustHiddenSpreads(RBF *rbf, VECTOR *v_error) ;
static int   rbfPrintActivations(RBF *rbf, VECTOR *v_obs, VECTOR *v_error, 
                                 int class, FILE *fp);
static int   rbfComputeOutputs(RBF *rbf) ;
static int   rbfNormalizeObservation(RBF *rbf, VECTOR *v_in, VECTOR *v_out) ;
static float rbfTrain(RBF *rbf, int (*get_observation_func)
               (VECTOR *v_obs, int no, void *parm, int *pclass), void *parm,
                      int which) ;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          allocate and initialize an rbf classifier.
------------------------------------------------------*/
RBF *
RBFinit(int ninputs, int noutputs, int max_clusters[], char *names[])
{
  RBF  *rbf ;
  int  i ;

  rbf = (RBF *)calloc(1, sizeof(RBF)) ;
  if (!rbf)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate rbf struct") ;
  rbf->ninputs = ninputs ;
  rbf->trate = DEFAULT_TRATE ;
  rbf->momentum = MOMENTUM ;
  rbf->base_momentum = MOMENTUM ;
  rbf->noutputs = noutputs ;
  rbf->v_outputs = RVectorAlloc(noutputs, MATRIX_REAL) ;
  if (!rbf->v_outputs)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_outputs") ;

  if (names) for (i = 0 ; i < noutputs ; i++)
  {
    rbf->class_names[i] = (char *)calloc(strlen(names[i]+1), sizeof(char)) ;
    if (!rbf->class_names[i])
      ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate class name %d",
                i) ;
    strcpy(rbf->class_names[i], names[i]) ;
  }
  else for (i = 0 ; i < noutputs ; i++)
  {
    char name[200] ;

    sprintf(name, "class %d", i) ;
    rbf->class_names[i] = (char *)calloc(strlen(name+1), sizeof(char)) ;
    if (!rbf->class_names[i])
      ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate class name %d",
                i) ;
    strcpy(rbf->class_names[i], name) ;
  }

  for (rbf->nhidden = i = 0 ; i < noutputs ; i++)
  {
    rbf->cs[i] = CSinit(max_clusters[i], ninputs, CLUSTER_DONT_NORMALIZE) ;
    if (!rbf->cs[i])
      ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate cs %d", i) ;
    rbf->nhidden += max_clusters[i] ;
  }

  
  rbf->v_biases = RVectorAlloc(noutputs, MATRIX_REAL) ;
  if (!rbf->v_biases)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_biases") ;
  rbf->m_wij = MatrixAlloc(rbf->nhidden, noutputs, MATRIX_REAL) ;
  if (!rbf->m_wij)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate m_wij") ;
  rbf->v_hidden = RVectorAlloc(rbf->nhidden, MATRIX_REAL) ;
  if (!rbf->v_hidden)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_hidden") ;
  rbf->v_z = (VECTOR **)calloc(rbf->nhidden, sizeof(VECTOR *)) ;
  if (!rbf->v_z)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_z") ;
  for (i = 0 ; i < rbf->nhidden ; i++)
  {
    rbf->v_z[i] = VectorAlloc(ninputs, MATRIX_REAL) ;
    if (!rbf->v_z[i])
      ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_z[%d]", i) ;
  }

  rbf->clusters = (CLUSTER **)calloc(rbf->nhidden, sizeof(CLUSTER *)) ;
  if (!rbf->clusters)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate cluster pointers");
  return(rbf) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          train an RBF classifier.
------------------------------------------------------*/
int
RBFtrain(RBF *rbf, int (*get_observation_func)
               (VECTOR *v_obs, int no, void *parm, int *pclass), void *parm,
         float momentum)
{
  int              i, class, c, cno ;
  CLUSTERING_PARM  cp ;
  CLUSTER_SET      *cs ;

  if (!FZERO(momentum))
    rbf->base_momentum = rbf->momentum = momentum ;

  /* do allocation of trainig-specific stuff */
  rbf->m_delta_sigma_inv = (MATRIX **)calloc(rbf->nhidden, sizeof(MATRIX *));
  if (!rbf->m_delta_sigma_inv)
    ErrorExit(ERROR_NO_MEMORY,"RBFtrain: could not allocate delta sigma") ;
  rbf->v_delta_means = (VECTOR **)calloc(rbf->nhidden, sizeof(VECTOR *));
  if (!rbf->v_delta_means)
    ErrorExit(ERROR_NO_MEMORY,"RBFtrain: could not allocate delta means") ;
  for (i = 0 ; i < rbf->nhidden ; i++)
  {
    rbf->m_delta_sigma_inv[i] = 
      MatrixAlloc(rbf->ninputs, rbf->ninputs, MATRIX_REAL) ;
    if (!rbf->m_delta_sigma_inv[i])
      ErrorExit(ERROR_NO_MEMORY,"RBFtrain: could not allocate delta sigma %d", 
                i) ;
    rbf->v_delta_means[i] = VectorAlloc(rbf->ninputs, MATRIX_REAL) ;
    if (!rbf->v_delta_means[i])
      ErrorExit(ERROR_NO_MEMORY,"RBFtrain: could not allocate delta means %d", 
                i) ;
  }
  rbf->v_delta_biases = RVectorAlloc(rbf->noutputs, MATRIX_REAL) ;
  if (!rbf->v_delta_biases)
    ErrorExit(ERROR_NO_MEMORY, "RBFtrain: could not allocate v_delta_biases") ;
  rbf->m_delta_wij = MatrixAlloc(rbf->nhidden, rbf->noutputs, MATRIX_REAL) ;
  if (!rbf->m_delta_wij)
    ErrorExit(ERROR_NO_MEMORY, "RBFtrain: could not allocate m_delta_wij") ;

  cp.parm = parm ;
  cp.get_observation_func = get_observation_func ;

  for (i = 0 ; i < rbf->noutputs ; i++)
  {
    cp.current_class = i ;
    if (CScluster(rbf->cs[i], rbf_get_obs_func, (void *)&cp) != NO_ERROR)
      return(Gerror) ;
  }

  /* fill in cluster pointers in rbf struct for convenience sake */
  for (cno = class = 0 ; class < rbf->noutputs ; class++)
  {
    cs = rbf->cs[class] ;
    for (c = 0 ; c < cs->nclusters ; c++, cno++)
      rbf->clusters[cno] = cs->clusters+c ;
  }

  /* now that the initial cluster positions have been established,
     train the RBF using gradient descent.
     */
  fprintf(stdout, "initial inverse covariance matrix: \n") ;
  MatrixPrint(stdout, rbf->clusters[0]->m_inverse) ;
  if (rbfGradientDescent(rbf, get_observation_func, parm) != NO_ERROR)
    return(Gerror) ;

  fprintf(stdout, "final inverse covariance matrix: \n") ;
  MatrixPrint(stdout, rbf->clusters[0]->m_inverse) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           free the memory allocated by an RBF structure.
------------------------------------------------------*/
int
RBFfree(RBF **prbf)
{
  RBF  *rbf ;
  int  i ;

  rbf = *prbf ;
  *prbf = NULL ;

  MatrixFree(&rbf->m_wij) ;
  VectorFree(&rbf->v_outputs) ;
  VectorFree(&rbf->v_hidden) ;
  VectorFree(&rbf->v_biases) ;
  
  free(rbf->clusters) ;
  if (rbf->m_delta_sigma_inv)  /* free training-specific stuff */
  {
    MatrixFree(&rbf->m_delta_wij) ;
    VectorFree(&rbf->v_delta_biases) ;
    for (i = 0 ; i < rbf->nhidden ; i++)
    {
      if (rbf->m_delta_sigma_inv[i])
        MatrixFree(&rbf->m_delta_sigma_inv[i]) ;
      if (rbf->v_delta_means[i])
        VectorFree(&rbf->v_delta_means[i]) ;
    }
    free(rbf->v_delta_means) ;
    free(rbf->m_delta_sigma_inv) ;
  }

  for (i = 0 ; i < rbf->nhidden ; i++)
    VectorFree(&rbf->v_z[i]) ;
  free(rbf->v_z) ;
  for (i = 0 ; i < rbf->noutputs ; i++)
  {
    if (rbf->class_names[i])
      free(rbf->class_names[i]) ;
  }
  for (i = 0 ; i < rbf->ninputs ; i++)
    CSfree(&rbf->cs[i]) ;

  free(rbf) ;
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          print out the contents of an RBF to a file.
------------------------------------------------------*/
int
RBFprint(RBF *rbf, FILE *fp)
{
  int  i ;

  fprintf(fp,"rbf with %d inputs %d hidden, and %d outputs\n",
          rbf->ninputs,rbf->nhidden,rbf->noutputs);
  for (i = 0 ; i < rbf->noutputs ; i++)
  {
    if (rbf->class_names[i])
      fprintf(fp, "class %s:\n", rbf->class_names[i]) ;
    else
      fprintf(fp, "class %d:\n", i) ;
    CSprint(rbf->cs[i], stderr) ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Function called by clustering software. It must sift out
          only those inputs which have the appropriate class.
------------------------------------------------------*/
static int
rbf_get_obs_func(VECTOR *v_obs, int obs_no, void *vcp)
{
  CLUSTERING_PARM *cp ;
  int             ret, classno ;
  
  cp = (CLUSTERING_PARM *)vcp ;

  do
  {
    ret = (*cp->get_observation_func)(v_obs, obs_no++, cp->parm, &classno) ;
  } while ((ret == NO_ERROR) && (classno != cp->current_class)) ;

  return(ret) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Function called by clustering software. It must sift out
          only those inputs which have the appropriate class.
------------------------------------------------------*/
static int
rbfGradientDescent(RBF *rbf, int (*get_observation_func)
               (VECTOR *v_obs, int no, void *parm, int *pclass), void *parm)
{
  VECTOR           *v_obs, *v_error ;
  int              obs_no, class ;
  float            sse ;

  v_obs = VectorAlloc(rbf->ninputs, MATRIX_REAL) ;
  v_error = VectorAlloc(rbf->noutputs, MATRIX_REAL) ;
#if 0
/*rbf->momentum = rbf->base_momentum = 0.0f ;*/
parm = (void *)fopen("c2.dat", "r") ;
#if 0
rbf->nhidden = 1 ;
RVECTOR_ELT(rbf->v_biases,1) = 0.0f ;
RVECTOR_ELT(rbf->v_biases,2) = 1.0f ;
rbf->m_wij->rptr[1][1] = 1.0f ;
rbf->m_wij->rptr[1][2] = -1.0f ;
#endif
VECTOR_ELT(rbf->clusters[0]->v_means,1) = 0.2f ;
VECTOR_ELT(rbf->clusters[0]->v_means,2) = 0.2f ;
VECTOR_ELT(rbf->clusters[1]->v_means,1) = 0.8f ;
VECTOR_ELT(rbf->clusters[1]->v_means,2) = 0.8f ;
rbf->clusters[0]->m_inverse->rptr[1][1] = .10f ;
rbf->clusters[0]->m_inverse->rptr[2][2] = .10f ;
rbf->clusters[0]->m_inverse->rptr[1][2] = 0.0f ;
rbf->clusters[0]->m_inverse->rptr[2][1] = 0.0f ;
rbf->clusters[1]->m_inverse->rptr[1][1] = .10f ;
rbf->clusters[1]->m_inverse->rptr[2][2] = .10f ;
rbf->clusters[1]->m_inverse->rptr[1][2] = 0.0f ;
rbf->clusters[1]->m_inverse->rptr[2][1] = 0.0f ;
#endif
  
  rbfShowClusterCenters(rbf, stderr) ;
  rbfTrain(rbf, get_observation_func, parm, TRAIN_OUTPUTS) ;
  rbfTrain(rbf, get_observation_func, parm, TRAIN_SPREADS) ;
  /*  rbfTrain(rbf, get_observation_func, parm, TRAIN_ALL) ;*/
  obs_no = 0 ;
  sse = 0.0f ;
  while ((*get_observation_func)(v_obs, obs_no++, parm, &class) == NO_ERROR)
  {
    rbfComputeHiddenActivations(rbf, v_obs) ;
    rbfComputeOutputs(rbf) ;
    if (Gdiag & DIAG_SHOW)
    {
      sse += rbfComputeErrors(rbf, class, v_error) ;
      rbfPrintActivations(rbf, v_obs, v_error, class, stderr) ;
    }
  }
  sse = sqrt(sse/obs_no) ;
  VectorFree(&v_error) ;
  VectorFree(&v_obs) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Given the observation vector v_obs, compute the activation
          of the hidden units (RBFs).
------------------------------------------------------*/
static int
rbfComputeHiddenActivations(RBF *rbf, VECTOR *v_obs)
{
  int          c ;
  CLUSTER      *cluster ;
  float        total ;

  for (total = 0.0f, c = 0 ; c < rbf->nhidden ; c++)
  {
    cluster = rbf->clusters[c] ;
    total += RVECTOR_ELT(rbf->v_hidden,c+1) = 
      rbfGaussian(cluster->m_inverse, cluster->v_means, v_obs, rbf->v_z[c]);
  }
#if 1
  /* normalize hidden layer activations */
  if (FZERO(total))  /* shouldn't ever happen */
    total = 1.0f ;
  for (c = 0 ; c < rbf->nhidden ; c++)
    RVECTOR_ELT(rbf->v_hidden,c+1) /= total ;
#endif

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Given the observation vector v_obs, compute the activation
          of the hidden units (RBFs).
------------------------------------------------------*/
static int
rbfShowClusterCenters(RBF *rbf, FILE *fp)
{
  int          c ;
  CLUSTER      *cluster ;

  for (c = 0 ; c < rbf->nhidden ; c++)
  {
    cluster = rbf->clusters[c] ;
    fprintf(stderr, "cluster %d has %d observations. Center:", 
            c, cluster->nsamples) ;
    MatrixPrintTranspose(stderr, cluster->v_means) ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the value of a gaussian with mean v_means and
          inverse covariance matrix m_sigma_inverse given the
          input v_obs.
------------------------------------------------------*/
static float
rbfGaussian(MATRIX *m_sigma_inverse,VECTOR *v_means,VECTOR *v_obs, VECTOR *v_z)
{
  float   val = 0.0f ;
  VECTOR  *v_zT ;
  MATRIX  *m_tmp1, *m_tmp2 ;

  VectorSubtract(v_obs, v_means, v_z) ;
  v_zT = VectorTranspose(v_z, NULL) ;

  m_tmp1 = MatrixMultiply(m_sigma_inverse, v_z, NULL) ;
  m_tmp2 = MatrixMultiply(v_zT, m_tmp1, NULL) ;
  val = exp(-.5f * m_tmp2->rptr[1][1]) ;

  MatrixFree(&m_tmp1) ;
  MatrixFree(&m_tmp2) ;
  VectorFree(&v_zT) ;
  return(val) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the value of a gaussian with mean v_means and
          inverse covariance matrix m_sigma_inverse given the
          input v_obs.
------------------------------------------------------*/
static float
rbfComputeErrors(RBF *rbf, int class, VECTOR *v_error)
{
  int  i ;
  float target, error, total_error ;

  for (total_error = 0.0f, i = 1 ; i <= rbf->noutputs ; i++)
  {
    if (class == (i-1))
      target = 1.0f ;
    else
      target = 0.0f ;
    error = target - RVECTOR_ELT(rbf->v_outputs, i) ;
    VECTOR_ELT(v_error,i) = error ;
    total_error += error*error ;
  }
  return(total_error) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Adjust the weights in the (linear) output layer
          using the LMS algorithm.
------------------------------------------------------*/
static int
rbfAdjustOutputWeights(RBF *rbf, VECTOR *v_error)
{
  int    i, j ;
  float  dE_dwi, Gi, delta_wij, delta_bias, one_minus_momentum, trate,
         momentum ;

  momentum = rbf->momentum ;
  trate = rbf->trate ;
  one_minus_momentum = trate * (1.0f - rbf->momentum) ;

  /* i refers to the hidden unit, while j is the output unit */
  for (j = 1 ; j <= rbf->noutputs ; j++)
  {
    delta_bias = VECTOR_ELT(v_error,j) * one_minus_momentum + 
      momentum * RVECTOR_ELT(rbf->v_delta_biases,j) ;
    RVECTOR_ELT(rbf->v_delta_biases,j) = delta_bias ;
    RVECTOR_ELT(rbf->v_biases,j) += delta_bias ;
  }

  for (i = 1 ; i <= rbf->nhidden ; i++)
  {
    dE_dwi = 0.0f ;  /* gradient for this weight */
    Gi = RVECTOR_ELT(rbf->v_hidden, i) ;

    /* adjust the weight connected to each output node for this hidden unit */
    for (j = 1 ; j <= rbf->noutputs ; j++)
    {
      delta_wij = momentum * rbf->m_delta_wij->rptr[i][j] + 
        one_minus_momentum * VECTOR_ELT(v_error,j) * Gi ;
      rbf->m_wij->rptr[i][j] += delta_wij ;
      rbf->m_delta_wij->rptr[i][j] = delta_wij ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Print the RBF activations to a file for debugging.
------------------------------------------------------*/
static int
rbfPrintActivations(RBF *rbf,VECTOR *v_obs,VECTOR *v_error,int class,FILE *fp)
{
  int i ;

/*    fprintf(fp, "rbf: ") ;*/
  for (i = 1 ; i <= v_obs->rows ; i++)
    fprintf(fp, "%+2.2f ", VECTOR_ELT(v_obs, i)) ;
  fprintf(fp, " --> ") ;
  for (i = 1 ; i <= rbf->nhidden ; i++)
    fprintf(fp, "%+2.2f ", RVECTOR_ELT(rbf->v_hidden, i)) ;
  fprintf(fp, " --> ") ;
  for (i = 1 ; i <= rbf->noutputs ; i++)
  {
    fprintf(fp, "%+2.2f (%+2.2f) ", RVECTOR_ELT(rbf->v_outputs, i),
            VECTOR_ELT(v_error,i)) ;
  }
  if (class >= 0)
    fprintf(fp, " : %d \n", class) ;
  else
    fprintf(fp, "\n") ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Adjust the centers (means) of the Gaussians in
          the hidden layer.
------------------------------------------------------*/
static int
rbfAdjustHiddenCenters(RBF *rbf, VECTOR *v_error)
{
  int     i, j ;
  float   Gi, one_minus_momentum, trate, momentum, total, delta ;
  VECTOR  *v_dE_dui, *v_means, *v_delta_means, *v_z ;
  CLUSTER *cluster ;

  momentum = rbf->momentum ;
  trate = rbf->trate ;
  one_minus_momentum = trate * (1.0f - rbf->momentum) ;

  /* derivative of SSE function w.r.t. the ith mean vector */
  v_dE_dui = VectorAlloc(rbf->ninputs, MATRIX_REAL) ;

  /* i refers to the hidden unit, while j is the output unit */
  for (i = 1 ; i <= rbf->nhidden ; i++)
  {
    cluster = rbf->clusters[i-1] ;
    Gi = RVECTOR_ELT(rbf->v_hidden, i) ;
    v_delta_means = rbf->v_delta_means[i-1] ;
    v_means = cluster->v_means ;
    v_z = rbf->v_z[i-1] ;

    /* adjust the weight connected to each output node for this hidden unit */
    for (total = 0.0f, j = 1 ; j <= rbf->noutputs ; j++)
    {
      delta = rbf->m_wij->rptr[i][j] * VECTOR_ELT(v_error, j) ;
      total += delta ;
    }
    MatrixMultiply(cluster->m_inverse, v_z, v_dE_dui) ;
    VectorScalarMul(v_dE_dui, Gi * one_minus_momentum * total, v_dE_dui) ;

    /* at this point v_dE_dui is the  gradient (except for the momentum turn), 
       now turn it into a 'weight' change.
       */
    VectorScalarMul(v_delta_means, momentum, v_delta_means) ;
    VectorAdd(v_delta_means, v_dE_dui, v_delta_means) ;
    VectorAdd(v_delta_means, v_means, v_means) ;
  }

  VectorFree(&v_dE_dui) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
         Adjust the spreads (covariance matrices) of the
         Gaussians in the hidden layer.
------------------------------------------------------*/
static int
rbfAdjustHiddenSpreads(RBF *rbf, VECTOR *v_error)
{
  int     i, j ;
  float   Gi, one_minus_momentum, trate, momentum, total, delta ;
  VECTOR  *v_z ;
  MATRIX  *m_dE_dsigma_inv, *m_sigma_inv, *m_delta_sigma_inv ;
  CLUSTER *cluster ;

  momentum = rbf->momentum ;
  trate = rbf->trate ;
  one_minus_momentum = -0.5f * trate * (1.0f - rbf->momentum) ;

  /* derivative of SSE function w.r.t. the ith mean vector */
  m_dE_dsigma_inv = MatrixAlloc(rbf->ninputs, rbf->ninputs, MATRIX_REAL) ;

  /* i refers to the hidden unit, while j is the output unit */
  for (i = 1 ; i <= rbf->nhidden ; i++)
  {
    cluster = rbf->clusters[i-1] ;
    Gi = RVECTOR_ELT(rbf->v_hidden, i) ;
    m_delta_sigma_inv = rbf->m_delta_sigma_inv[i-1] ;
    m_sigma_inv = cluster->m_inverse ;
    v_z = rbf->v_z[i-1] ;

    /* adjust the weight connected to each output node for this hidden unit */
    for (total = 0.0f, j = 1 ; j <= rbf->noutputs ; j++)
    {
      delta = rbf->m_wij->rptr[i][j] * VECTOR_ELT(v_error, j) ;
      total += delta ;
    }
    VectorOuterProduct(v_z, v_z, m_dE_dsigma_inv) ;
    MatrixScalarMul(m_dE_dsigma_inv, Gi * one_minus_momentum * total, 
                    m_dE_dsigma_inv) ;

    /* at this point m_dE_dsigma_inv is the  gradient (except for the 
       momentum turn), now turn it into a 'weight' change.
       */
    MatrixScalarMul(m_delta_sigma_inv, momentum, m_delta_sigma_inv) ;
    MatrixAdd(m_delta_sigma_inv, m_dE_dsigma_inv, m_delta_sigma_inv) ;
    MatrixAdd(m_delta_sigma_inv, m_sigma_inv, m_sigma_inv) ;
  }

  MatrixFree(&m_dE_dsigma_inv) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the outputs of the RBF (assumes 
          rbfComputeHiddenActivations has already been called).
------------------------------------------------------*/
static int
rbfComputeOutputs(RBF *rbf)
{
#if 0
  int    i ;
  float  min_out, val, total ;
#endif

  MatrixMultiply(rbf->v_hidden, rbf->m_wij, rbf->v_outputs) ;
  VectorAdd(rbf->v_biases, rbf->v_outputs, rbf->v_outputs) ;

#if 0
  /* normalize the outputs to sum to 1 and be in [0,1] */
  min_out = RVECTOR_ELT(rbf->v_outputs, 1) ;
  for (i = 2 ; i <= rbf->noutputs ; i++)
  {
    val = RVECTOR_ELT(rbf->v_outputs, i) ;
    if (val < min_out)
      min_out = val ;
  }
  for (total = 0.0f, i = 1 ; i <= rbf->noutputs ; i++)
  {
    val = RVECTOR_ELT(rbf->v_outputs, i) - min_out ;
    total += val ;
    RVECTOR_ELT(rbf->v_outputs, i) = val ;
  }
  if (FZERO(total))   /* should never happen */
    total = 1.0f ;
  for (i = 1 ; i <= rbf->noutputs ; i++)
    RVECTOR_ELT(rbf->v_outputs, i) /= total ;
#endif

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:
          The class with the largest activation.

        Description
          Classify an observation vector.
------------------------------------------------------*/
int
RBFclassify(RBF *rbf, VECTOR *v_obs)
{
  int    class, c ;
  float  max_val, val ;

  rbfComputeHiddenActivations(rbf, v_obs) ;
  rbfComputeOutputs(rbf) ;
  max_val = RVECTOR_ELT(rbf->v_outputs, 1) ;
  class = 0 ;
  for (c = 2 ; c <= rbf->noutputs ; c++)
  {
    val = RVECTOR_ELT(rbf->v_outputs, c) ;
    if (val > max_val)
    {
      max_val = val ;
      class = c-1 ;
    }
  }
  return(class) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
rbfNormalizeObservation(RBF *rbf, VECTOR *v_in, VECTOR *v_out)
{
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static float
rbfTrain(RBF *rbf, int (*get_observation_func)
         (VECTOR *v_obs, int no, void *parm, int *pclass), void *parm,int which)
{
  VECTOR           *v_obs, *v_error ;
  int              obs_no, class, epoch, nsmall = 0, nnegative = 0 ;
  float            error, sse, old_sse = 0.0f, delta_sse, rms ;

  v_obs = VectorAlloc(rbf->ninputs, MATRIX_REAL) ;
  v_error = VectorAlloc(rbf->noutputs, MATRIX_REAL) ;
  rbf->trate = DEFAULT_TRATE ;
  fprintf(stderr, "\n") ;
  for (epoch = 0 ; epoch < MAX_EPOCHS ; epoch++)
  {
    obs_no = 0 ;
    sse = 0.0f ;
    while ((*get_observation_func)(v_obs, obs_no++, parm, &class) == NO_ERROR)
    {
      rbfNormalizeObservation(rbf, v_obs, v_obs) ;
      rbfComputeHiddenActivations(rbf, v_obs) ;
      rbfComputeOutputs(rbf) ;
      error = rbfComputeErrors(rbf, class, v_error) ;
      sse += error ;
      if (which & TRAIN_CENTERS)
        rbfAdjustHiddenCenters(rbf, v_error) ;
      if (which & TRAIN_SPREADS)
        rbfAdjustHiddenSpreads(rbf, v_error) ;
      if (which & TRAIN_OUTPUTS)
        rbfAdjustOutputWeights(rbf, v_error) ;
    }
    rms = sqrt(sse / (float)obs_no) ;
    fprintf(stderr, "\rerror: %2.5f (trate %2.4f)", rms, rbf->trate) ;
    delta_sse = sse - old_sse ;
    if (epoch > 1)
    {
      if (sse > old_sse * ERROR_RATIO)  /* increased by a fair amount */
      {
        rbf->momentum = 0.0f ;
        rbf->trate *= TRATE_DECREASE ;
      }
      else
      {
        rbf->momentum = rbf->base_momentum ;
        if (sse < old_sse)  /* error decreased, increase training rate */
          rbf->trate *= TRATE_INCREASE ;
        else                /* error increased a little */
          rbf->trate *= TRATE_DECREASE ;
      }
    }
    if (delta_sse < 0)
      nnegative++ ;
    else
      nnegative = 0 ;

    if (fabs(delta_sse) < MIN_DELTA_SSE)
    {
      if (nsmall++ > MAX_SMALL)
      {
        if (nnegative > MAX_SMALL/2)
          continue ;
        break ;
      }
    }
    else
      nsmall = 0 ;

    old_sse = sse ;
    if (rbf->trate < MIN_TRATE)
      rbf->trate = MIN_TRATE ;
  }
  fprintf(stderr, "- training done in %d epochs\n", epoch) ;
  return(sse) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
RBFwrite(RBF *rbf, char *fname)
{
  int   i ;
  FILE  *fp ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, 
                (ERROR_NO_FILE, "RBFwrite(%s): could not open file", fname)) ;

  fprintf(fp, "%d %d %d\n", rbf->noutputs, rbf->ninputs, rbf->nhidden) ;
  fprintf(fp, "\n# class names and max # of clusters\n") ;
  for (i = 0 ; i < rbf->noutputs ; i++)
    fprintf(fp, "%d %s\n", rbf->cs[i]->max_clusters, rbf->class_names[i]) ;

  fprintf(fp, "# weights:\n") ;
  MatrixAsciiWriteInto(fp, rbf->m_wij) ;
  fprintf(fp, "# biases:\n") ;
  MatrixAsciiWriteInto(fp, rbf->v_biases) ;
  for (i = 0 ; i < rbf->noutputs ; i++)
  {
    fprintf(fp, "CLASS: %s\n", rbf->class_names[i]) ;
    CSwriteInto(fp, rbf->cs[i]) ;
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
RBF *
RBFread(char *fname)
{
  int   i, c, cno, class,noutputs,ninputs,nhidden,max_clusters[MAX_OUTPUTS];
  char  *names[MAX_OUTPUTS], *cp, line[100] ;
  FILE  *fp ;
  RBF   *rbf ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, 
                (ERROR_NO_FILE, "RBFread(%s): could not open file", fname)) ;

  if (fscanf(fp, "%d %d %d\n", &noutputs, &ninputs, &nhidden) != 3)
    ErrorReturn(NULL, (ERROR_BADFILE, "RBFread(%s): could not scan parms",
                       fname)) ;

  for (i = 0 ; i < noutputs ; i++)
  {
    cp = fgetl(line, 199, fp) ;
    if (sscanf(cp, "%d", &max_clusters[i]) != 1)
      ErrorReturn(NULL, (ERROR_BADFILE, 
                         "RBFread(%s): could not read class # of clusters",
                         fname,i));
    cp = StrSkipNumber(cp) ;
    names[i] = (char *)calloc(strlen(cp)+1, sizeof(char)) ;
    strcpy(names[i], cp) ;
  }
  rbf = RBFinit(ninputs, noutputs, max_clusters, names) ;
  for (i = 0 ; i < noutputs ; i++)
    free(names[i]) ;

  MatrixAsciiReadFrom(fp, rbf->m_wij) ;
  MatrixAsciiReadFrom(fp, rbf->v_biases) ;
  for (i = 0 ; i < rbf->noutputs ; i++)
  {
    fgetl(line, 199, fp) ;   /* skip class name */
    CSreadFrom(fp, rbf->cs[i]) ;
  }

  /* fill in cluster pointers in rbf struct for convenience sake */
  for (cno = class = 0 ; class < rbf->noutputs ; class++)
  {
    for (c = 0 ; c < rbf->cs[class]->nclusters ; c++, cno++)
      rbf->clusters[cno] = rbf->cs[class]->clusters+c ;
  }

  fclose(fp) ;
  return(rbf) ;
}

