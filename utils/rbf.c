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

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define TRATE_DECREASE  0.7f
#define TRATE_INCREASE  1.05f
#define ERROR_RATIO     1.04f
#define MIN_TRATE       0.001f
#define DEFAULT_TRATE   0.05f
#define MOMENTUM        0.8f
/* 
   if the sse doesn't change by more than MIN_DELTA_SSE for MAX_SMALL
   epochs in a row, assume that training has asymptoted.
*/
#define MIN_DELTA_SSE   MIN_TRATE  
#define MAX_SMALL       5


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
                         VECTOR *v_obs);
static int   rbfShowClusterCenters(RBF *rbf, FILE *fp) ;
static float rbfComputeErrors(RBF *rbf, int class, VECTOR *v_error) ;
static int   rbfAdjustOutputWeights(RBF *rbf, VECTOR *v_error) ;

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
  for (rbf->nhidden = i = 0 ; i < ninputs ; i++)
  {
    rbf->cs[i] = CSinit(max_clusters[i], ninputs, CLUSTER_NORMALIZE) ;
    if (!rbf->cs[i])
      ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate cs %d", i) ;
    rbf->nhidden += max_clusters[i] ;
  }

  
  rbf->v_biases = RVectorAlloc(noutputs, MATRIX_REAL) ;
  if (!rbf->v_biases)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_biases") ;
  rbf->v_delta_biases = RVectorAlloc(noutputs, MATRIX_REAL) ;
  if (!rbf->v_delta_biases)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_delta_biases") ;
  rbf->m_wij = MatrixAlloc(rbf->nhidden, noutputs, MATRIX_REAL) ;
  if (!rbf->m_wij)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate m_wij") ;
  rbf->m_delta_wij = MatrixAlloc(rbf->nhidden, noutputs, MATRIX_REAL) ;
  if (!rbf->m_delta_wij)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate m_delta_wij") ;
  rbf->v_hidden = RVectorAlloc(noutputs, MATRIX_REAL) ;
  if (!rbf->v_hidden)
    ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_hidden") ;
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
               (VECTOR *v_obs, int no, void *parm, int *pclass), void *parm)
{
  int              i ;
  CLUSTERING_PARM  cp ;

  cp.parm = parm ;
  cp.get_observation_func = get_observation_func ;

  for (i = 0 ; i < rbf->noutputs ; i++)
  {
    cp.current_class = i ;
    CScluster(rbf->cs[i], rbf_get_obs_func, (void *)&cp) ;
  }

  /* now that the initial cluster positions have been established,
     train the RBF using gradient descent.
     */
  rbfGradientDescent(rbf, get_observation_func, parm) ;

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
  MatrixFree(&rbf->m_delta_wij) ;
  VectorFree(&rbf->v_outputs) ;
  VectorFree(&rbf->v_hidden) ;
  VectorFree(&rbf->v_biases) ;
  VectorFree(&rbf->v_delta_biases) ;

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

  fprintf(fp,"rbf with %d inputs and %d outputs\n",rbf->ninputs,rbf->noutputs);
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
#define MAX_EPOCHS      2000

static int
rbfGradientDescent(RBF *rbf, int (*get_observation_func)
               (VECTOR *v_obs, int no, void *parm, int *pclass), void *parm)
{
  VECTOR           *v_obs, *v_error ;
  int              obs_no, class, i, epoch, nsmall = 0 ;
  CLUSTERING_PARM  cp ;
  float            error, sse, old_sse = 0.0f, delta_sse ;

  cp.parm = parm ;
  cp.get_observation_func = get_observation_func ;

  v_obs = VectorAlloc(rbf->ninputs, MATRIX_REAL) ;
  v_error = VectorAlloc(rbf->noutputs, MATRIX_REAL) ;
  rbfShowClusterCenters(rbf, stderr) ;

  fprintf(stderr, "\n") ;
  for (epoch = 0 ; epoch < MAX_EPOCHS ; epoch++)
  {
    obs_no = 0 ;
    sse = 0.0f ;
    while ((*get_observation_func)(v_obs, obs_no++, parm, &class) == NO_ERROR)
    {
      rbfComputeHiddenActivations(rbf, v_obs) ;
      MatrixMultiply(rbf->v_hidden, rbf->m_wij, rbf->v_outputs) ;
/*      VectorAdd(rbf->v_biases, rbf->v_outputs, rbf->v_outputs) ;*/
      error = rbfComputeErrors(rbf, class, v_error) ;
      sse += (error*error) ;
      rbfAdjustOutputWeights(rbf, v_error) ;

#if 0     
      fprintf(stderr, "rbf: ") ;
      for (i = 1 ; i <= v_obs->rows ; i++)
        fprintf(stderr, "%+2.2f ", VECTOR_ELT(v_obs, i)) ;
      fprintf(stderr, " --> ") ;
      for (i = 1 ; i <= rbf->nhidden ; i++)
        fprintf(stderr, "%+2.2f ", RVECTOR_ELT(rbf->v_hidden, i)) ;
      fprintf(stderr, " --> ") ;
      for (i = 1 ; i <= rbf->noutputs ; i++)
      {
        fprintf(stderr, "%+2.2f (%+2.2f) ", RVECTOR_ELT(rbf->v_outputs, i),
                VECTOR_ELT(v_error,i)) ;
      }
      fprintf(stderr, " : %d \n", class) ;
#endif
    }
    fprintf(stderr, "\rerror: %2.5f (trate %2.4f)", sse, rbf->trate) ;
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
    if (fabs(delta_sse) < MIN_DELTA_SSE)
    {
      if (nsmall++ > MAX_SMALL)
        break ;
    }
    else
      nsmall = 0 ;

    old_sse = sse ;
    if (rbf->trate < MIN_TRATE)
      rbf->trate = MIN_TRATE ;
  }
  fprintf(stderr, " done after %d epochs\n", epoch) ;
  obs_no = 0 ;
  while ((*get_observation_func)(v_obs, obs_no++, parm, &class) == NO_ERROR)
  {
    rbfComputeHiddenActivations(rbf, v_obs) ;
    MatrixMultiply(rbf->v_hidden, rbf->m_wij, rbf->v_outputs) ;
    rbfComputeErrors(rbf, class, v_error) ;
    rbfAdjustOutputWeights(rbf, v_error) ;
    
    fprintf(stderr, "rbf: ") ;
    for (i = 1 ; i <= v_obs->rows ; i++)
      fprintf(stderr, "%+2.2f ", VECTOR_ELT(v_obs, i)) ;
    fprintf(stderr, " --> ") ;
    for (i = 1 ; i <= rbf->nhidden ; i++)
      fprintf(stderr, "%+2.2f ", RVECTOR_ELT(rbf->v_hidden, i)) ;
    fprintf(stderr, " --> ") ;
    for (i = 1 ; i <= rbf->noutputs ; i++)
    {
      fprintf(stderr, "%+2.2f (%+2.2f) ", RVECTOR_ELT(rbf->v_outputs, i),
              VECTOR_ELT(v_error,i)) ;
    }
    fprintf(stderr, " : %d \n", class) ;
  }
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
  int          class, c, rbf_no ;
  CLUSTER_SET  *cs ;
  CLUSTER      *cluster ;

  for (rbf_no = class = 0 ; class < rbf->noutputs ; class++)
  {
    cs = rbf->cs[class] ;
    for (c = 0 ; c < cs->nclusters ; c++, rbf_no++)
    {
      cluster = cs->clusters+c ;
      RVECTOR_ELT(rbf->v_hidden,rbf_no+1) = 
        rbfGaussian(cluster->m_inverse, cluster->v_means, v_obs) ;
    }
  }
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
  int          class, c, rbf_no ;
  CLUSTER_SET  *cs ;
  CLUSTER      *cluster ;

  for (rbf_no = class = 0 ; class < rbf->noutputs ; class++)
  {
    cs = rbf->cs[class] ;
    for (c = 0 ; c < cs->nclusters ; c++, rbf_no++)
    {
      cluster = cs->clusters+c ;
      fprintf(stderr, "cluster %d has %d observations. Center:", 
              rbf_no, cluster->nsamples) ;
      MatrixPrintTranspose(stderr, cluster->v_seed) ;
    }
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
rbfGaussian(MATRIX *m_sigma_inverse,VECTOR *v_means,VECTOR *v_obs)
{
  float   val = 0.0f ;
  VECTOR  *v_z, *v_zT ;
  MATRIX  *m_tmp1, *m_tmp2 ;

  v_z = VectorSubtract(v_obs, v_means, NULL) ;
  v_zT = VectorTranspose(v_z, NULL) ;

  m_tmp1 = MatrixMultiply(m_sigma_inverse, v_z, NULL) ;
  m_tmp2 = MatrixMultiply(v_zT, m_tmp1, NULL) ;
  val = exp(-.5f * m_tmp2->rptr[1][1]) ;

  MatrixFree(&m_tmp1) ;
  MatrixFree(&m_tmp2) ;
  VectorFree(&v_z) ;
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
    total_error += error ;
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
  float  dE_dwi, Gi, delta_wij, /*delta_bias, */one_minus_momentum, trate,
         momentum ;

  momentum = rbf->momentum ;
  one_minus_momentum = 1.0f - rbf->momentum ;
  trate = rbf->trate ;

  /* i refers to the hidden unit, while j is the output unit */
#if 0
  for (j = 1 ; j <= rbf->noutputs ; j++)
  {
    delta_bias = trate * VECTOR_ELT(v_error,j) * one_minus_momentum + 
      momentum * RVECTOR_ELT(rbf->v_delta_biases,j) ;
    RVECTOR_ELT(rbf->v_delta_biases,j) = delta_bias ;
    RVECTOR_ELT(rbf->v_biases,j) += delta_bias ;
  }
#endif

  for (i = 1 ; i <= rbf->nhidden ; i++)
  {
    dE_dwi = 0.0f ;  /* gradient for this weight */
    Gi = RVECTOR_ELT(rbf->v_hidden, i) ;

    /* adjust the weight connected to each output node for this hidden unit */
    for (j = 1 ; j <= rbf->noutputs ; j++)
    {
      delta_wij = momentum * rbf->m_delta_wij->rptr[i][j] + 
        one_minus_momentum * trate * VECTOR_ELT(v_error,j) * Gi ;
      rbf->m_wij->rptr[i][j] += delta_wij ;
      rbf->m_delta_wij->rptr[i][j] = delta_wij ;
    }
  }

  return(NO_ERROR) ;
}

