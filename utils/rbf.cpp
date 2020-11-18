/*
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

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cluster.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "proto.h"
#include "rbf.h"
#include "utils.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define TRATE_DECREASE 0.8f
#define TRATE_INCREASE 1.05f
#define ERROR_RATIO 1.05f
#define MIN_TRATE 0.000001f
#define DEFAULT_TRATE 0.001f
#define MOMENTUM 0.8f
#define MAX_EPOCHS 4000 /* maximum # of training iterations */

/*
   if the sse doesn't change by more than MIN_DELTA_SSE for MAX_SMALL
   epochs in a row, assume that training has asymptoted.
*/
#define MIN_DELTA_SSE 0.001f /* 0.01 % change in sse */
#define MAX_SMALL 5
#define MAX_POSITIVE 5

#define TRAIN_OUTPUTS 0x001
#define TRAIN_SPREADS 0x002
#define TRAIN_CENTERS 0x004
#define TRAIN_HIDDEN (TRAIN_SPREADS | TRAIN_CENTERS)
#define TRAIN_ALL (TRAIN_OUTPUTS | TRAIN_SPREADS | TRAIN_CENTERS)

/*-----------------------------------------------------
                      STRUCTURES
-------------------------------------------------------*/

typedef struct
{
  int current_class;
  int (*get_observation_func)(VECTOR *v, int obs_no, void *parm, int same_class, int *classno);
  void *parm;
} CLUSTERING_PARM;

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static int rbf_get_obs_func(VECTOR *v_obs, int obs_no, void *vrbf);
static int rbfGradientDescent(
    RBF *rbf, int (*get_observation_func)(VECTOR *v_obs, int no, void *parm, int same_class, int *pclass), void *parm);
static int rbfComputeHiddenActivations(RBF *rbf, VECTOR *v_obs);
static float rbfGaussian(MATRIX *m_sigma_inverse, VECTOR *v_means, VECTOR *v_obs, VECTOR *v_z, float norm);
static int rbfShowClusterCenters(RBF *rbf, FILE *fp);
static int rbfAdjustOutputWeights(RBF *rbf, VECTOR *v_error);
static int rbfAdjustHiddenCenters(RBF *rbf, VECTOR *v_error);
static int rbfAdjustHiddenSpreads(RBF *rbf, VECTOR *v_error);
static int rbfComputeOutputs(RBF *rbf);
static int rbfNormalizeObservation(RBF *rbf, VECTOR *v_in, VECTOR *v_out);
static float rbfTrain(RBF *rbf,
                      int (*get_observation_func)(VECTOR *v_obs, int no, void *parm, int same_class, int *pclass),
                      void *parm,
                      int which);
static float rbfComputeCurrentError(
    RBF *rbf, int (*get_observation_func)(VECTOR *v_obs, int no, void *parm, int same_class, int *pclass), void *parm);
static int rbfCalculateOutputWeights(
    RBF *rbf, int (*get_observation_func)(VECTOR *v_obs, int no, void *parm, int same_class, int *pclass), void *parm);
static int rbfAllocateTrainingParameters(RBF *rbf);

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          allocate and initialize an rbf classifier.
------------------------------------------------------*/
RBF *RBFinit(int ninputs, int noutputs, int max_clusters[], const char *names[])
{
  RBF *rbf;
  int i;

  rbf = (RBF *)calloc(1, sizeof(RBF));
  if (!rbf) ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate rbf struct");
  rbf->ninputs = ninputs;
  rbf->trate = DEFAULT_TRATE;
  rbf->momentum = MOMENTUM;
  rbf->base_momentum = MOMENTUM;
  rbf->noutputs = noutputs;
  rbf->min_inputs = (float *)calloc(ninputs, sizeof(float));
  if (!rbf->min_inputs) ErrorExit(ERROR_NOMEMORY, "RBFinit: could not allocate min_inputs");
  rbf->max_inputs = (float *)calloc(ninputs, sizeof(float));
  if (!rbf->max_inputs) ErrorExit(ERROR_NOMEMORY, "RBFinit: could not allocate max_inputs");
  rbf->v_outputs = RVectorAlloc(noutputs, MATRIX_REAL);
  if (!rbf->v_outputs) ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_outputs");

  if (names)
    for (i = 0; i < noutputs; i++) {
      rbf->class_names[i] = (char *)calloc(strlen(names[i]) + 1, sizeof(char));
      if (!rbf->class_names[i]) ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate class name %d", i);
      strcpy(rbf->class_names[i], names[i]);
    }
  else
    for (i = 0; i < noutputs; i++) {
      char name[200];

      sprintf(name, "class %d", i);
      rbf->class_names[i] = (char *)calloc(strlen(name + 1), sizeof(char));
      if (!rbf->class_names[i]) ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate class name %d", i);
      strcpy(rbf->class_names[i], name);
    }

  for (rbf->nhidden = i = 0; i < noutputs; i++) {
    if (max_clusters[i] > 0) {
      rbf->cs[i] = CSinit(max_clusters[i], ninputs, CLUSTER_DONT_NORMALIZE);
      if (!rbf->cs[i]) ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate cs %d", i);
      rbf->nhidden += max_clusters[i];
    }
  }

  rbf->m_wij = MatrixAlloc(rbf->nhidden + 1, noutputs, MATRIX_REAL);
  if (!rbf->m_wij) ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate m_wij");
  rbf->v_hidden = RVectorAlloc(rbf->nhidden + 1, MATRIX_REAL);
  if (!rbf->v_hidden) ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_hidden");
  RVECTOR_ELT(rbf->v_hidden, rbf->nhidden + 1) = 1.0f; /* bias input */

  rbf->v_z = (VECTOR **)calloc(rbf->nhidden, sizeof(VECTOR *));
  if (!rbf->v_z) ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_z");
  for (i = 0; i < rbf->nhidden; i++) {
    rbf->v_z[i] = VectorAlloc(ninputs, MATRIX_REAL);
    if (!rbf->v_z[i]) ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate v_z[%d]", i);
  }

  rbf->clusters = (CLUSTER **)calloc(rbf->nhidden, sizeof(CLUSTER *));
  if (!rbf->clusters) ErrorExit(ERROR_NO_MEMORY, "RBFinit: could not allocate cluster pointers");
  return (rbf);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          train an RBF classifier.
------------------------------------------------------*/
int RBFtrain(RBF *rbf,
             int (*get_observation_func)(VECTOR *v_obs, int no, void *parm, int same_class, int *pclass),
             void *parm,
             float momentum)
{
  int i, classnum, c, cno;
  CLUSTERING_PARM cp;
  CLUSTER_SET *cs;

  if ((momentum < 1.0f) && (momentum >= 0.0f)) rbf->base_momentum = rbf->momentum = momentum;

  /* must examine training set before allocating training set data, because
     we need to know how many observations are in the training set.
     */
  RBFexamineTrainingSet(rbf, get_observation_func, parm);

  /* do allocation of training-specific stuff */
  if (rbfAllocateTrainingParameters(rbf) != NO_ERROR) return (Gerror);

  cp.parm = parm;
  cp.get_observation_func = get_observation_func;

  for (i = 0; i < rbf->noutputs; i++) {
    if (!rbf->cs[i]) continue;

    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "finding %d clusters for class %s...", rbf->cs[i]->max_clusters, rbf->class_names[i]);
    cp.current_class = i;
    if (CScluster(rbf->cs[i], rbf_get_obs_func, (void *)&cp) != NO_ERROR)
      ErrorReturn(Gerror, (Gerror, "RBFtrain: clustering failed for class %s", rbf->class_names[i]));
    {
      char fname[100];
      FILE *fp;
      sprintf(fname, "cluster.%s.dat", rbf->class_names[i]);
      fp = fopen(fname, "w");
      CSwriteInto(fp, rbf->cs[i]);
      fclose(fp);
    }
    if (Gdiag & DIAG_SHOW) fprintf(stderr, "done.\n");
  }

  /* fill in cluster pointers in rbf struct for convenience sake */
  for (cno = classnum = 0; classnum < rbf->noutputs; classnum ++) {
    cs = rbf->cs[classnum];
    if (cs)
      for (c = 0; c < cs->nclusters; c++, cno++) rbf->clusters[cno] = cs->clusters + c;
  }

  /* now that the initial cluster positions have been established,
     train the RBF using gradient descent.
     */
  if (rbfGradientDescent(rbf, get_observation_func, parm) != NO_ERROR) return (Gerror);

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          perform more training on an already trained RBF classifier.
------------------------------------------------------*/
int RBFretrain(RBF *rbf,
               int (*get_observation_func)(VECTOR *v_obs, int no, void *parm, int same_class, int *pclass),
               void *parm,
               float momentum)
{
  if ((momentum < 1.0f) && (momentum >= 0.0f)) rbf->base_momentum = rbf->momentum = momentum;

  /* must examine training set before allocating training set data, because
     we need to know how many observations are in the training set.
     */
  RBFexamineTrainingSet(rbf, get_observation_func, parm);

  /* do allocation of training-specific stuff */
  if (rbfAllocateTrainingParameters(rbf) != NO_ERROR) return (Gerror);

  /* now that the initial cluster positions have been established,
     train the RBF using gradient descent.
     */
  if (rbfGradientDescent(rbf, get_observation_func, parm) != NO_ERROR) return (Gerror);

  if (Gdiag & DIAG_SHOW) fprintf(stderr, "\ntraining network...");
  rbfTrain(rbf, get_observation_func, parm, TRAIN_ALL);
  if (Gdiag & DIAG_SHOW) fprintf(stderr, "training complete, calculating output weights...");

  if (rbfCalculateOutputWeights(rbf, get_observation_func, parm) != NO_ERROR)
    rbfTrain(rbf, get_observation_func, parm, TRAIN_OUTPUTS);
  if (Gdiag & DIAG_SHOW) {
    fprintf(stderr, "done.\n");
    rbfShowClusterCenters(rbf, stderr);
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           free the memory allocated by an RBF structure.
------------------------------------------------------*/
int RBFfree(RBF **prbf)
{
  RBF *rbf;
  int i;

  rbf = *prbf;
  *prbf = NULL;

  if (!rbf) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "RBFfree: null pointer"));

  if (rbf->min_inputs) free(rbf->min_inputs);
  if (rbf->max_inputs) free(rbf->max_inputs);
  MatrixFree(&rbf->m_wij);
  VectorFree(&rbf->v_outputs);
  VectorFree(&rbf->v_hidden);

  free(rbf->clusters);
  if (rbf->m_delta_sigma_inv) /* free training-specific stuff */
  {
    if (rbf->observed) free(rbf->observed);
    if (rbf->m_delta_wij) MatrixFree(&rbf->m_delta_wij);
    for (i = 0; i < rbf->nhidden; i++) {
      if (rbf->m_delta_sigma_inv[i]) MatrixFree(&rbf->m_delta_sigma_inv[i]);
      if (rbf->v_delta_means[i]) VectorFree(&rbf->v_delta_means[i]);
    }
    if (rbf->v_delta_means) free(rbf->v_delta_means);
    if (rbf->m_delta_sigma_inv) free(rbf->m_delta_sigma_inv);
  }

  if (rbf->v_z) {
    for (i = 0; i < rbf->nhidden; i++) VectorFree(&rbf->v_z[i]);
    free(rbf->v_z);
  }

  for (i = 0; i < rbf->noutputs; i++) {
    if (rbf->class_names[i]) free(rbf->class_names[i]);
  }
  for (i = 0; i < rbf->ninputs; i++)
    if (rbf->cs[i]) CSfree(&rbf->cs[i]);

  free(rbf);
  return (NO_ERROR);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          print out the contents of an RBF to a file.
------------------------------------------------------*/
int RBFprint(RBF *rbf, FILE *fp)
{
  int i, c;
  CLUSTER_SET *cs;
  CLUSTER *cluster;

  fprintf(fp, "rbf with %d inputs %d hidden, and %d outputs\n", rbf->ninputs, rbf->nhidden, rbf->noutputs);
  for (i = 0; i < rbf->noutputs; i++) {
    if (rbf->class_names[i])
      fprintf(fp, "class %s:\n", rbf->class_names[i]);
    else
      fprintf(fp, "class %d:\n", i);
    cs = rbf->cs[i];
    if (!cs) continue;
    for (c = 0; c < cs->nclusters; c++) {
      cluster = cs->clusters + c;
      fprintf(fp, "cluster %d has %d observations. Center: ", c, cluster->nsamples);
      MatrixPrintTranspose(fp, cluster->v_means);
    }
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Function called by clustering software. It must sift out
          only those inputs which have the appropriate class.

          Note that obs_no refers to the total # of observations,
          while class_obs refers to the # of observations for the
          given class.
------------------------------------------------------*/
static int rbf_get_obs_func(VECTOR *v_obs, int desired_class_obs_no, void *vcp)
{
  CLUSTERING_PARM *cp;
  int ret, classno, obs_no, class_obs_no;
  static int last_class = -1;
  static int last_class_obs = -1;
  static int last_obs = -1;

  cp = (CLUSTERING_PARM *)vcp;

  desired_class_obs_no++; /* it is a count - not an index */
  classno = cp->current_class;
  if (classno == last_class && desired_class_obs_no > last_class_obs) {
    /* start at one past previous observation, not at start */
    class_obs_no = last_class_obs;
    obs_no = last_obs + 1;
  }
  else
    class_obs_no = obs_no = 0;

  do {
    ret = (*cp->get_observation_func)(v_obs, obs_no++, cp->parm, 1, &classno);
    if ((ret == NO_ERROR) && (classno == cp->current_class)) class_obs_no++;
  } while ((ret == NO_ERROR) && (class_obs_no < desired_class_obs_no));

  if (ret == NO_ERROR) {
    last_class = classno;
    last_obs = obs_no - 1;
    last_class_obs = desired_class_obs_no;
  }
  else
    last_class = last_obs = last_class_obs = -1;

  return (ret);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Function called by clustering software. It must sift out
          only those inputs which have the appropriate class.
------------------------------------------------------*/
static int rbfGradientDescent(
    RBF *rbf, int (*get_observation_func)(VECTOR *v_obs, int no, void *parm, int same_class, int *pclass), void *parm)
{
  if (Gdiag & DIAG_SHOW) {
    rbfShowClusterCenters(rbf, stderr);
    fprintf(stderr, "calculating output weights...");
  }

  if (rbfCalculateOutputWeights(rbf, get_observation_func, parm) != NO_ERROR)
    rbfTrain(rbf, get_observation_func, parm, TRAIN_OUTPUTS);
  if (Gdiag & DIAG_SHOW) fprintf(stderr, "\ntraining network...");
  /*  rbfTrain(rbf, get_observation_func, parm, TRAIN_ALL) ;*/
  if (Gdiag & DIAG_SHOW) {
    fprintf(stderr, "training complete.\n");
    rbfShowClusterCenters(rbf, stderr);
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Given the observation vector v_obs, compute the activation
          of the hidden units (RBFs).
------------------------------------------------------*/
static int rbfComputeHiddenActivations(RBF *rbf, VECTOR *v_obs)
{
  int c;
  CLUSTER *cluster;
  float total;

  for (total = 0.0f, c = 0; c < rbf->nhidden; c++) {
    cluster = rbf->clusters[c];
    total += RVECTOR_ELT(rbf->v_hidden, c + 1) =
        rbfGaussian(cluster->m_inverse, cluster->v_means, v_obs, rbf->v_z[c], cluster->norm);
  }
  /* normalize hidden layer activations */
  if (FZERO(total)) /* shouldn't ever happen */
    total = 1.0f;
  for (c = 0; c < rbf->nhidden; c++) RVECTOR_ELT(rbf->v_hidden, c + 1) /= total;

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Given the observation vector v_obs, compute the activation
          of the hidden units (RBFs).
------------------------------------------------------*/
static int rbfShowClusterCenters(RBF *rbf, FILE *fp)
{
  int i, c;
  CLUSTER *cluster;
  CLUSTER_SET *cs;

  for (i = 0; i < rbf->noutputs; i++) {
    cs = rbf->cs[i];
    if (!cs) continue;
    fprintf(fp, "classnum %s:\n", rbf->class_names[i]);
    for (c = 0; c < cs->nclusters; c++) {
      cluster = cs->clusters + c;
      fprintf(fp, "cluster %d has %d observations. Center:", c, cluster->nsamples);

      MatrixPrintTranspose(fp, cluster->v_means);
    }
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the value of a gaussian with mean v_means and
          inverse covariance matrix m_sigma_inverse given the
          input v_obs.
------------------------------------------------------*/
static float rbfGaussian(MATRIX *m_sigma_inverse, VECTOR *v_means, VECTOR *v_obs, VECTOR *v_z, float norm)
{
  float val = 0.0f;
  static VECTOR *v_zT = NULL;
  static MATRIX *m_tmp1 = NULL, *m_tmp2 = NULL;

  VectorSubtract(v_obs, v_means, v_z);
  v_zT = VectorTranspose(v_z, v_zT);

  m_tmp1 = MatrixMultiply(m_sigma_inverse, v_z, m_tmp1);
  m_tmp2 = MatrixMultiply(v_zT, m_tmp1, m_tmp2);
  val = -.5f * m_tmp2->rptr[1][1];
  if (val < -15.0f) /* prevent underflow */
    val = 0.0f;
  else if (val > 0.0f) /* error - shouldn't happen, but can because of lms */
    val = 0.0f;
  else
    val = norm * exp(val);

  return (val);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the value of a gaussian with mean v_means and
          inverse covariance matrix m_sigma_inverse given the
          input v_obs.
------------------------------------------------------*/
float RBFcomputeErrors(RBF *rbf, int classnum, VECTOR *v_error)
{
  int i;
  float target, error, total_error;

  for (total_error = 0.0f, i = 1; i <= rbf->noutputs; i++) {
    if (classnum == (i - 1))
      target = 1.0f;
    else
      target = 0.0f;
    error = target - RVECTOR_ELT(rbf->v_outputs, i);

    VECTOR_ELT(v_error, i) = error;
    total_error += error * error;
  }
  return (total_error);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Adjust the weights in the (linear) output layer
          using the LMS algorithm.
------------------------------------------------------*/
static int rbfAdjustOutputWeights(RBF *rbf, VECTOR *v_error)
{
  int i, j;
  float Gi, delta_wij, one_minus_momentum, trate, momentum;
  // float dE_dwi;

  momentum = rbf->momentum;
  trate = rbf->trate;
  one_minus_momentum = trate * (1.0f - rbf->momentum);

  for (i = 1; i <= rbf->nhidden + 1; i++) {
    // dE_dwi = 0.0f; /* gradient for this weight */
    Gi = RVECTOR_ELT(rbf->v_hidden, i);

    /* adjust the weight connected to each output node for this hidden unit */
    for (j = 1; j <= rbf->noutputs; j++) {
      delta_wij = momentum * rbf->m_delta_wij->rptr[i][j] + one_minus_momentum * VECTOR_ELT(v_error, j) * Gi;
      rbf->m_wij->rptr[i][j] += delta_wij;
      rbf->m_delta_wij->rptr[i][j] = delta_wij;
    }
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Print the RBF activations to a file for debugging.
------------------------------------------------------*/
int RBFprintActivations(RBF *rbf, VECTOR *v_obs, VECTOR *v_error, int classnum, FILE *fp)
{
  int i;

  /*    fprintf(fp, "rbf: ") ;*/
  for (i = 1; i <= v_obs->rows; i++) fprintf(fp, "%+2.2f ", VECTOR_ELT(v_obs, i));
  fprintf(fp, " --> ");
  for (i = 1; i <= rbf->nhidden; i++) fprintf(fp, "%+2.2f ", RVECTOR_ELT(rbf->v_hidden, i));
  fprintf(fp, " --> ");
  for (i = 1; i <= rbf->noutputs; i++) {
    fprintf(fp, "%+2.2f ", RVECTOR_ELT(rbf->v_outputs, i));
    if (v_error) fprintf(fp, "(%+2.2f) ", VECTOR_ELT(v_error, i));
  }
  if (classnum >= 0)
    fprintf(fp, " : %d \n", classnum);
  else
    fprintf(fp, "\n");
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Adjust the centers (means) of the Gaussians in
          the hidden layer.
------------------------------------------------------*/
static int rbfAdjustHiddenCenters(RBF *rbf, VECTOR *v_error)
{
  int i, j;
  float Gi, one_minus_momentum, trate, momentum, total, delta;
  VECTOR *v_means, *v_delta_means, *v_z;
  CLUSTER *cluster;
  static VECTOR *v_dE_dui = NULL;

  momentum = rbf->momentum;
  trate = rbf->trate;
  one_minus_momentum = trate * (1.0f - rbf->momentum);

  /* i refers to the hidden unit, while j is the output unit */
  for (i = 1; i <= rbf->nhidden; i++) {
    cluster = rbf->clusters[i - 1];
    Gi = RVECTOR_ELT(rbf->v_hidden, i);
    v_delta_means = rbf->v_delta_means[i - 1];
    v_means = cluster->v_means;
    v_z = rbf->v_z[i - 1];

    /* adjust the weight connected to each output node for this hidden unit */
    for (total = 0.0f, j = 1; j <= rbf->noutputs; j++) {
      delta = rbf->m_wij->rptr[i][j] * VECTOR_ELT(v_error, j);
      total += delta;
    }
    v_dE_dui = MatrixMultiply(cluster->m_inverse, v_z, v_dE_dui);
    VectorScalarMul(v_dE_dui, Gi * one_minus_momentum * total, v_dE_dui);

    /* at this point v_dE_dui is the  gradient (except for the momentum turn),
       now turn it into a 'weight' change.
       */
    VectorScalarMul(v_delta_means, momentum, v_delta_means);
    VectorAdd(v_delta_means, v_dE_dui, v_delta_means);
    VectorAdd(v_delta_means, v_means, v_means);
  }

  VectorFree(&v_dE_dui);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
         Adjust the spreads (covariance matrices) of the
         Gaussians in the hidden layer.
------------------------------------------------------*/
static int rbfAdjustHiddenSpreads(RBF *rbf, VECTOR *v_error)
{
  int i, j;
  float Gi, one_minus_momentum, trate, momentum, total, delta, det;
  VECTOR *v_z;
  MATRIX *m_sigma_inv, *m_delta_sigma_inv;
  CLUSTER *cluster;
  static MATRIX *m_sigma_inv_save = NULL, *m_dE_dsigma_inv = NULL;

  momentum = rbf->momentum;
  trate = rbf->trate;
  one_minus_momentum = -0.5f * trate * (1.0f - rbf->momentum);

  /* i refers to the hidden unit, while j is the output unit */
  for (i = 1; i <= rbf->nhidden; i++) {
    cluster = rbf->clusters[i - 1];
    Gi = RVECTOR_ELT(rbf->v_hidden, i);
    m_delta_sigma_inv = rbf->m_delta_sigma_inv[i - 1];
    m_sigma_inv = cluster->m_inverse;
    v_z = rbf->v_z[i - 1];
    det = MatrixDeterminant(m_sigma_inv);
    if (det < 0.0f) continue; /* shouldn't happen, something wrong - don't bother with it */

    /* adjust the weight connected to each output node for this hidden unit */
    for (total = 0.0f, j = 1; j <= rbf->noutputs; j++) {
      delta = rbf->m_wij->rptr[i][j] * VECTOR_ELT(v_error, j);
      total += delta;
    }
    m_dE_dsigma_inv = VectorOuterProduct(v_z, v_z, m_dE_dsigma_inv);
    MatrixScalarMul(m_dE_dsigma_inv, Gi * one_minus_momentum * total, m_dE_dsigma_inv);

    /* at this point m_dE_dsigma_inv is the  gradient (except for the
       momentum turn), now turn it into a 'weight' change.
       */
    MatrixScalarMul(m_delta_sigma_inv, momentum, m_delta_sigma_inv);
    MatrixAdd(m_delta_sigma_inv, m_dE_dsigma_inv, m_delta_sigma_inv);
    m_sigma_inv_save = MatrixCopy(m_sigma_inv, m_sigma_inv_save);
    MatrixAdd(m_delta_sigma_inv, m_sigma_inv, m_sigma_inv);
    det = MatrixDeterminant(m_sigma_inv);
    while (det < 0.0f) /* step was too large, matrix should be positive def. */
    {
      MatrixCopy(m_sigma_inv_save, m_sigma_inv);
      MatrixScalarMul(m_delta_sigma_inv, 0.1f, m_delta_sigma_inv);
      MatrixAdd(m_delta_sigma_inv, m_sigma_inv, m_sigma_inv);
      det = MatrixDeterminant(m_sigma_inv);
    }
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the outputs of the RBF (assumes
          rbfComputeHiddenActivations has already been called).
------------------------------------------------------*/
static int rbfComputeOutputs(RBF *rbf)
{

  MatrixMultiply(rbf->v_hidden, rbf->m_wij, rbf->v_outputs);


  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:
          The class with the largest activation.

        Description
          Classify an observation vector.
------------------------------------------------------*/
int RBFclassify(RBF *rbf, VECTOR *v_obs)
{
  int classnum, c;
  float max_val, val;

  rbfComputeHiddenActivations(rbf, v_obs);
  rbfComputeOutputs(rbf);
  max_val = RVECTOR_ELT(rbf->v_outputs, 1);
  classnum = 0;
  for (c = 2; c <= rbf->noutputs; c++) {
    val = RVECTOR_ELT(rbf->v_outputs, c);
    if (val > max_val) {
      max_val = val;
      classnum = c - 1;
    }
  }
  return (classnum);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int rbfNormalizeObservation(RBF *rbf, VECTOR *v_in, VECTOR *v_out) { return (NO_ERROR); }
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static float rbfTrain(RBF *rbf,
                      int (*get_observation_func)(VECTOR *v_obs, int no, void *parm, int same_class, int *pclass),
                      void *parm,
                      int which)
{
  VECTOR *v_obs, *v_error;
  int obs_no, classnum, epoch = 0, nsmall = 0, nnegative = 0, nobs, positive = 0;
  float error, sse = 0.0f, delta_sse, rms, old_sse;
  RBF *rbf_save = NULL;

  v_obs = VectorAlloc(rbf->ninputs, MATRIX_REAL);
  v_error = VectorAlloc(rbf->noutputs, MATRIX_REAL);
  rbf->trate = DEFAULT_TRATE;
  rbf->momentum = rbf->base_momentum;

  /* first compute initial error on training set */
  old_sse = rbfComputeCurrentError(rbf, get_observation_func, parm);

  fprintf(stderr, "\n");
  if (!FZERO(old_sse))
    for (epoch = 0; epoch < MAX_EPOCHS; epoch++) {
      rbf_save = RBFcopyWeights(rbf, rbf_save);
      memset(rbf->observed, 0, rbf->nobs * sizeof(unsigned char));
      for (sse = 0.0f, nobs = 0; nobs < rbf->nobs; nobs++) {
        do {
          obs_no = nint(randomNumber(0.0, (double)(rbf->nobs - 1)));
        } while (rbf->observed[obs_no] != 0);
        if ((*get_observation_func)(v_obs, obs_no, parm, 0, &classnum) != NO_ERROR)
          ErrorExit(ERROR_BADPARM,
                    "rbfTrain: observation function failed to"
                    " obtain training sample # %d",
                    obs_no);
        rbf->observed[obs_no] = 1; /* mark this sample as used this epoch */

        rbfNormalizeObservation(rbf, v_obs, v_obs);
        rbfComputeHiddenActivations(rbf, v_obs);
        rbfComputeOutputs(rbf);
        error = RBFcomputeErrors(rbf, classnum, v_error);
        sse += error;
        if (which & TRAIN_CENTERS) rbfAdjustHiddenCenters(rbf, v_error);
        if (which & TRAIN_SPREADS) rbfAdjustHiddenSpreads(rbf, v_error);
        if (which & TRAIN_OUTPUTS) rbfAdjustOutputWeights(rbf, v_error);
      }
      delta_sse = sse - old_sse;
      if (sse > old_sse * ERROR_RATIO) /* increased by a fair amount */
      {
        rbf->momentum = 0.0f;
        RBFcopyWeights(rbf_save, rbf); /* restore old weights */
        rbf->trate *= TRATE_DECREASE;
        sse = old_sse;
      }
      else /* accept new network, error either went down or up a little */
      {
        rbf->momentum = rbf->base_momentum;
        if (sse < old_sse) /* error decreased, increase training rate */
          rbf->trate *= TRATE_INCREASE;
        else /* error increased a little */
          rbf->trate *= TRATE_DECREASE;
        old_sse = sse;
      }

      if (FZERO(sse)) break;

      rms = sqrt(sse / (float)rbf->nobs);
      fprintf(stderr, "\repoch %d, error: %2.5f (trate %2.6f) (nsmall %d)", epoch, rms, rbf->trate, nsmall);
      if (delta_sse < 0) {
        nnegative++;
        positive = 0;
      }
      else {
        nnegative = 0;
        if (delta_sse > 0) positive++;
      }

      if (positive > MAX_POSITIVE) break;

      if (fabs(delta_sse / old_sse) < MIN_DELTA_SSE) /* a really small step */
      {
        if (nsmall++ > MAX_SMALL) /* MAX_SMALL really small steps in a row */
        {
          if (nnegative > MAX_SMALL) /* all of them were downhill - continue */
            continue;
          break; /* too many small steps in a row, assume we have asymptoted */
        }
      }
      else
        nsmall = 0;

      if (rbf->trate < MIN_TRATE) rbf->trate = MIN_TRATE;
    }
  if (rbf_save) RBFfree(&rbf_save);
  fprintf(stderr, "- training done in %d epochs\n", epoch);
  return (sse);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int RBFwrite(RBF *rbf, char *fname)
{
  FILE *fp;
  int error;

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "RBFwrite(%s): could not open file", fname));

  error = RBFwriteInto(rbf, fp);
  fclose(fp);
  return (error);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
RBF *RBFread(char *fname)
{
  FILE *fp;
  RBF *rbf;

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_NO_FILE, "RBFread(%s): could not open file", fname));

  rbf = RBFreadFrom(fp);
  fclose(fp);
  return (rbf);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Run through the training set once to find the min and max
          values of the various input dimensions, as well as the
          total # of training samples.
------------------------------------------------------*/
int RBFexamineTrainingSet(RBF *rbf,
                          int (*get_observation_func)(VECTOR *v_obs, int no, void *parm, int same_class, int *pclass),
                          void *parm)
{
  int obs_no = 0, classnum, row, c, i;
  VECTOR *v_obs;
  float *means, *stds, v, mean;

  v_obs = VectorAlloc(rbf->ninputs, MATRIX_REAL);

  if ((*get_observation_func)(v_obs, obs_no++, parm, 0, &classnum) != NO_ERROR)
    ErrorExit(ERROR_BADPARM, "rbfExamineTrainingSet: no samples in set");

  /* not really a good idea to be messing with the CS internal parameters,
     but it's much easier than the alternative, because all the other CS stuff
     is class-specific, but the means and stds for normalizing the inputs
     should be across all classes.
     */
  means = rbf->cs[0]->means;
  stds = rbf->cs[0]->stds;

  /* initialize min and max to 1st elt in training set. */
  for (row = 1; row <= rbf->ninputs; row++) {
    v = VECTOR_ELT(v_obs, row);
    rbf->min_inputs[row - 1] = v;
    rbf->max_inputs[row - 1] = v;
    means[row - 1] = v;
    stds[row - 1] = v * v;
  }

  while ((*get_observation_func)(v_obs, obs_no++, parm, 0, &classnum) == NO_ERROR) {
    for (i = 0; i < rbf->ninputs; i++) {
      v = VECTOR_ELT(v_obs, i + 1);
      if (v > rbf->max_inputs[i]) rbf->max_inputs[i] = v;
      if (v < rbf->min_inputs[i]) rbf->min_inputs[i] = v;
      means[i] += v;
      stds[i] += (v * v);
    }
  }

  rbf->nobs = obs_no - 1;
  for (i = 0; i < rbf->ninputs; i++) {
    mean = means[i];
    if (obs_no) mean /= (float)rbf->nobs;
    means[i] = mean;
    stds[i] = sqrt(stds[i] / rbf->nobs - mean * mean);
  }

  /* copy means from 0th CLUSTER_SET to all other CLUSTER_SETs */
  for (c = 0; c < rbf->noutputs; c++) {
    for (i = 0; i < rbf->ninputs; i++) {
      if (rbf->cs[c]) {
        rbf->cs[c]->means[i] = means[i];
        rbf->cs[c]->stds[i] = stds[i];
      }
    }
  }

  VectorFree(&v_obs);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Save the weights of one RBF in another, for possible
          restoration.
------------------------------------------------------*/
RBF *RBFcopyWeights(RBF *rbf_src, RBF *rbf_dst)
{
  int i;

  if (!rbf_dst) {
    int max_clusters[MAX_OUTPUTS], cno, classnum, c;
    CLUSTER_SET *cs_src, *cs_dst;

    for (i = 0; i < rbf_src->noutputs; i++) max_clusters[i] = rbf_src->cs[i] ? rbf_src->cs[i]->max_clusters : 0;

    rbf_dst = RBFinit(rbf_src->ninputs, rbf_src->noutputs, max_clusters, const_cast<const char**>(rbf_src->class_names));

    /* fill in cluster pointers in rbf struct for convenience sake */
    for (cno = classnum = 0; classnum < rbf_src->noutputs; classnum ++) {
      cs_src = rbf_src->cs[classnum];
      cs_dst = rbf_dst->cs[classnum];
      if (cs_src)
        for (c = 0; c < cs_src->nclusters; c++, cno++) rbf_dst->clusters[cno] = cs_dst->clusters + c;
    }
    rbf_dst->m_delta_sigma_inv = (MATRIX **)calloc(rbf_dst->nhidden, sizeof(MATRIX *));
    if (!rbf_dst->m_delta_sigma_inv) ErrorExit(ERROR_NO_MEMORY, "RBFcopyWeights: could not alloc delta sigma");
    rbf_dst->v_delta_means = (VECTOR **)calloc(rbf_dst->nhidden, sizeof(VECTOR *));
    if (!rbf_dst->v_delta_means) ErrorExit(ERROR_NO_MEMORY, "RBFcopyWeights: could not alloc delta means");
  }
  rbf_dst->noutputs = rbf_src->noutputs;
  rbf_dst->ninputs = rbf_src->ninputs;
  rbf_dst->nhidden = rbf_src->nhidden;
  rbf_dst->nobs = rbf_src->nobs;
  rbf_dst->current_class = rbf_src->current_class;

  /* save output weights */
  MatrixCopy(rbf_src->m_wij, rbf_dst->m_wij);
  rbf_dst->m_delta_wij = MatrixCopy(rbf_src->m_delta_wij, rbf_dst->m_delta_wij);

  /* save hidden layer weights */
  for (i = 0; i < rbf_src->nhidden; i++) {
    rbf_dst->clusters[i]->m_inverse = MatrixCopy(rbf_src->clusters[i]->m_inverse, rbf_dst->clusters[i]->m_inverse);
    VectorCopy(rbf_src->clusters[i]->v_means, rbf_dst->clusters[i]->v_means);
    rbf_dst->m_delta_sigma_inv[i] = MatrixCopy(rbf_src->m_delta_sigma_inv[i], rbf_dst->m_delta_sigma_inv[i]);
    rbf_dst->v_delta_means[i] = VectorCopy(rbf_src->v_delta_means[i], rbf_dst->v_delta_means[i]);
  }

  return (rbf_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the sse on the training set.
------------------------------------------------------*/
static float rbfComputeCurrentError(
    RBF *rbf, int (*get_observation_func)(VECTOR *v_obs, int no, void *parm, int same_class, int *pclass), void *parm)
{
  float error, sse, obs_no;
  // float rms;
  VECTOR *v_obs, *v_error;
  int classnum;

  v_obs = VectorAlloc(rbf->ninputs, MATRIX_REAL);
  v_error = VectorAlloc(rbf->noutputs, MATRIX_REAL);
  for (sse = 0.0f, obs_no = 0; obs_no < rbf->nobs; obs_no++) {
    if ((*get_observation_func)(v_obs, obs_no, parm, 0, &classnum) != NO_ERROR)
      ErrorExit(ERROR_BADPARM,
                "rbfTrain: observation function failed to"
                " obtain training sample # %d",
                obs_no);
    rbfNormalizeObservation(rbf, v_obs, v_obs);
    rbfComputeHiddenActivations(rbf, v_obs);
    rbfComputeOutputs(rbf);
    error = RBFcomputeErrors(rbf, classnum, v_error);
    sse += error;
  }
  // rms = sqrt(sse / (float)rbf->nobs);

  VectorFree(&v_obs);
  VectorFree(&v_error);
  return (sse);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Calculate the optimal MSE output weights using the
          pseudo-inverse, where the observations correspond to
          the vector of hidden activations (G), and the targets
          are the vector d. Thus, calculate  Gt G in one pass
          through the training data, and invert it. Then calculate
          Gt D where D is a matrix, the jth column of which is the
          desired output of the jth node for each training pattern.
          The pth row is then the desired output vector for the
          pth pattern.
------------------------------------------------------*/
static int rbfCalculateOutputWeights(
    RBF *rbf, int (*get_observation_func)(VECTOR *v_obs, int no, void *parm, int same_class, int *pclass), void *parm)
{
  VECTOR *v_obs, *v_hidden;
  int classnum, obs_no, i, k;
  MATRIX *m_Gt_G, *m_Gt_G_inv, *m_Gt_D;

  v_obs = VectorAlloc(rbf->ninputs, MATRIX_REAL);

  m_Gt_G = MatrixAlloc(rbf->nhidden + 1, rbf->nhidden + 1, MATRIX_REAL);
  m_Gt_D = MatrixAlloc(rbf->nhidden + 1, rbf->noutputs, MATRIX_REAL);
  v_hidden = rbf->v_hidden;

  /* go through the training set again to calculate (Gt G) and Gt D
   where D is the vector of desired values */
  for (obs_no = 0; obs_no < rbf->nobs; obs_no++) {
    if ((*get_observation_func)(v_obs, obs_no, parm, 0, &classnum) != NO_ERROR)
      ErrorExit(ERROR_BADPARM,
                "rbfCalculateOutputWeights: observation "
                " function failed to obtain training sample # %d",
                obs_no);
    rbfNormalizeObservation(rbf, v_obs, v_obs);
    rbfComputeHiddenActivations(rbf, v_obs);

    for (i = 1; i <= rbf->nhidden + 1; i++) {
      /* only contributs to this class (all other outputs are 0) */
      m_Gt_D->rptr[i][classnum + 1] += RVECTOR_ELT(v_hidden, i);
      for (k = 1; k <= rbf->nhidden + 1; k++) m_Gt_G->rptr[i][k] += RVECTOR_ELT(v_hidden, i) * RVECTOR_ELT(v_hidden, k);
    }
  }


  if (MatrixSingular(m_Gt_G)) {
    fprintf(stderr, "regularizing matrix...\n");
    /*    MatrixRegularize(m_Gt_G, m_Gt_G) ;*/
    m_Gt_G_inv = MatrixSVDInverse(m_Gt_G, NULL);
  }
  else
    m_Gt_G_inv = MatrixInverse(m_Gt_G, NULL);
  if (!m_Gt_G_inv) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "rbfCalculateOutputWeights: matrix is singular"));

  MatrixMultiply(m_Gt_G_inv, m_Gt_D, rbf->m_wij);


  MatrixFree(&m_Gt_G);
  MatrixFree(&m_Gt_D);
  MatrixFree(&m_Gt_G_inv);
  VectorFree(&v_obs);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Allocate training-specific parameters.
------------------------------------------------------*/
static int rbfAllocateTrainingParameters(RBF *rbf)
{
  int i;

  rbf->observed = (unsigned char *)calloc(rbf->nobs, sizeof(unsigned char));
  if (!rbf->observed) ErrorExit(ERROR_NO_MEMORY, "rbfAllocateTrainingParms: could not allocate observed array");
  rbf->m_delta_sigma_inv = (MATRIX **)calloc(rbf->nhidden, sizeof(MATRIX *));
  if (!rbf->m_delta_sigma_inv) ErrorExit(ERROR_NO_MEMORY, "rbfAllocateTrainingParms: could not allocate delta sigma");
  rbf->v_delta_means = (VECTOR **)calloc(rbf->nhidden, sizeof(VECTOR *));
  if (!rbf->v_delta_means) ErrorExit(ERROR_NO_MEMORY, "rbfAllocateTrainingParms: could not allocate delta means");
  for (i = 0; i < rbf->nhidden; i++) {
    rbf->m_delta_sigma_inv[i] = MatrixAlloc(rbf->ninputs, rbf->ninputs, MATRIX_REAL);
    if (!rbf->m_delta_sigma_inv[i])
      ErrorExit(ERROR_NO_MEMORY, "rbfAllocateTrainingParms: could not allocate delta sigma %d", i);
    rbf->v_delta_means[i] = VectorAlloc(rbf->ninputs, MATRIX_REAL);
    if (!rbf->v_delta_means[i])
      ErrorExit(ERROR_NO_MEMORY, "rbfAllocateTrainingParms: could not allocate delta means %d", i);
  }
  rbf->m_delta_wij = MatrixAlloc(rbf->nhidden + 1, rbf->noutputs, MATRIX_REAL);
  if (!rbf->m_delta_wij) ErrorExit(ERROR_NO_MEMORY, "rbfAllocateTrainingParms: could not allocate m_delta_wij");

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Write a classifier into a previously opened file.
------------------------------------------------------*/
int RBFwriteInto(RBF *rbf, FILE *fp)
{
  int i;

  fprintf(fp, "%d %d %d\n", rbf->noutputs, rbf->ninputs, rbf->nhidden);
  fprintf(fp, "\n# classnum names and max # of clusters\n");
  for (i = 0; i < rbf->noutputs; i++)
    fprintf(fp, "%d %s\n", rbf->cs[i] ? rbf->cs[i]->max_clusters : 0, rbf->class_names[i]);

  fprintf(fp, "# weights:\n");
  MatrixAsciiWriteInto(fp, rbf->m_wij);
  for (i = 0; i < rbf->noutputs; i++) {
    fprintf(fp, "CLASS: %s\n", rbf->class_names[i]);
    if (rbf->cs[i]) CSwriteInto(fp, rbf->cs[i]);
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
         Read a classifier from a previously opened file
------------------------------------------------------*/
RBF *RBFreadFrom(FILE *fp)
{
  int i, c, cno, classnum, noutputs, ninputs, nhidden, max_clusters[MAX_OUTPUTS];
  char *cp, line[100];
  char *names[MAX_OUTPUTS];
  RBF *rbf;

  if (fscanf(fp, "%d %d %d\n", &noutputs, &ninputs, &nhidden) != 3)
    ErrorReturn(NULL, (ERROR_BADFILE, "RBFread: could not scan parms"));

  for (i = 0; i < noutputs; i++) {
    cp = fgetl(line, 199, fp);
    if (sscanf(cp, "%d", &max_clusters[i]) != 1)
      ErrorReturn(NULL, (ERROR_BADFILE, "RBFread(%s): could not read class # of clusters", i));
    cp = StrSkipNumber(cp);
    names[i] = (char *)calloc(strlen(cp) + 1, sizeof(char));
    strcpy(names[i], cp);
  }
  rbf = RBFinit(ninputs, noutputs, max_clusters, const_cast<const char**>(names));
  for (i = 0; i < noutputs; i++) free(names[i]);

  MatrixAsciiReadFrom(fp, rbf->m_wij);
  for (i = 0; i < rbf->noutputs; i++) {
    fgetl(line, 199, fp); /* skip class name */
    if (rbf->cs[i]) CSreadFrom(fp, rbf->cs[i]);
  }

  /* fill in cluster pointers in rbf struct for convenience sake */
  for (cno = classnum = 0; classnum < rbf->noutputs; classnum ++) {
    if (rbf->cs[classnum])
      for (c = 0; c < rbf->cs[classnum]->nclusters; c++, cno++) rbf->clusters[cno] = rbf->cs[classnum]->clusters + c;
  }

  return (rbf);
}

