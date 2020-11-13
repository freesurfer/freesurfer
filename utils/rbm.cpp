/**
 * @brief utilities for Restricted Boltzmann Machines
 *
 * Reference:
 * Hinton, G. E., Osindero, S., and Teh, Y. W. (2006a).
 * A fast learning algorithm for deep belief nets.
 * Neural Computation, 18(7):1527{1554.
 */
/*
 * Original Author: Bruce Fischl
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

#include "rbm.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "utils.h"
#include "voxlist.h"

#include "romp_support.h"

int dump_hidden(RBM *rbm, char *fname);
int dump_visible(RBM *rbm, char *fname);

static int compare_pvals(const void *v1, const void *v2)
{
  LABEL_PVAL *p1, *p2;

  p1 = (LABEL_PVAL *)v1;
  p2 = (LABEL_PVAL *)v2;

  if (p1->pval < p2->pval) return (-1);
  if (p1->pval > p2->pval) return (+1);
  return (0);  // equal
}

int RBMsortLabelProbabilities(RBM *rbm)
{
  int l;

  for (l = 0; l < rbm->nlabels; l++) {
    rbm->label_pvals[l].label = l;
    rbm->label_pvals[l].pval = rbm->labels[l];
  }
  qsort(rbm->label_pvals, rbm->nlabels, sizeof(LABEL_PVAL), compare_pvals);

  return (NO_ERROR);
}

static MRI *weights_to_mri(RBM *rbm)
{
  MRI *mri;
  int v, h, k1, k2;

  if (rbm->input_type == RBM_INPUT_IMAGE) {
    mri = MRIallocSequence(rbm->ksize, rbm->ksize, 1, MRI_FLOAT, rbm->nhidden);
    for (h = 0; h < rbm->nhidden; h++) {
      for (v = k1 = 0; k1 < rbm->ksize; k1++)
        for (k2 = 0; k2 < rbm->ksize; k2++, v++) MRIsetVoxVal(mri, k1, k2, 0, h, rbm->weights[v][h]);
    }
    MRIcopyHeader(rbm->mri_inputs, mri);
  }
  else {
    mri = MRIalloc(rbm->nvisible, rbm->nhidden, 1, MRI_FLOAT);
    for (v = 0; v < rbm->nvisible; v++)
      for (h = 0; h < rbm->nhidden; h++) MRIsetVoxVal(mri, v, h, 0, 0, rbm->weights[v][h]);
  }
  return (mri);
}

static MRI *layer_weights_to_mri(DBN *dbn, int layer)
{
  MRI *mri, *mri_previous_layer, *mri_tmp = NULL;
  int v, h;
  RBM *rbm_first = dbn->rbms[0], *rbm;

  if (layer == 0) return (weights_to_mri(rbm_first));

  if (rbm_first->input_type == RBM_INPUT_IMAGE) {
    /*
      each frame in the previous layer if a hidden node for it, but a visible node for us. Create
      an image that has as many frames as we have hidden nodes, and each frame is a linear combination of all
      the frames in the previous layer with the weights given by the connection strength.
    */
    mri_previous_layer = layer_weights_to_mri(dbn, layer - 1);
    rbm = dbn->rbms[layer];
    mri = MRIallocSequence(rbm_first->ksize, rbm_first->ksize, 1, MRI_FLOAT, rbm->nhidden);
    for (h = 0; h < rbm->nhidden; h++) {
      for (v = 0; v < rbm->nvisible; v++) {
        mri_tmp = MRIcopyFrame(mri_previous_layer, mri_tmp, v, 0);
        MRIscalarMul(mri_tmp, mri_tmp, rbm->weights[v][h]);
        MRIaddToFrame(mri, mri_tmp, mri, h, h);
      }
    }
    MRIfree(&mri_tmp);
    MRIfree(&mri_previous_layer);
    MRIcopyHeader(rbm->mri_inputs, mri);
  }
  else {
    rbm = dbn->rbms[layer];
    mri = MRIalloc(rbm->nvisible, rbm->nhidden, 1, MRI_FLOAT);
    for (v = 0; v < rbm->nvisible; v++)
      for (h = 0; h < rbm->nhidden; h++) MRIsetVoxVal(mri, v, h, 0, 0, rbm->weights[v][h]);
  }
  return (mri);
}
static int dump_gradients(
    RBM *rbm, double *dvisible_bias, double *dvariance, double *dhidden_bias, double **dw, RBM_PARMS *parms, int step)
{
  int v, h;

  printf("visible state:  ");
  for (v = 0; v < rbm->nvisible; v++) printf("%2.0f ", rbm->visible[v]);
  printf("\n");
  printf("hidden state:   ");
  for (h = 0; h < rbm->nhidden; h++) printf("%d ", (int)rbm->hidden_state[h]);
  printf("\n");

  printf("visible bias:   ");
  for (v = 0; v < rbm->nvisible; v++) printf("%+2.5f ", rbm->visible_bias[v]);
  printf("\n");
  printf("variance:       ");
  for (v = 0; v < rbm->nvisible; v++) printf("%+2.5f ", exp(rbm->variance[v]));
  printf("\n");
  printf("hidden bias:    ");
  for (h = 0; h < rbm->nhidden; h++) printf("%+2.5f ", rbm->hidden_bias[h]);
  printf("\n");

  printf("weights:\n");
  for (v = 0; v < rbm->nvisible; v++) {
    for (h = 0; h < rbm->nhidden; h++) printf("%+2.5f ", rbm->weights[v][h]);
    if (v < rbm->nvisible - 1) printf("\n");
  }
  printf("\nGRADIENTS\n");
  printf("visible bias:    ");
  for (v = 0; v < rbm->nvisible; v++) printf("%+2.5f ", dvisible_bias[v]);
  printf("\n");
  printf("variance:        ");
  for (v = 0; v < rbm->nvisible; v++) printf("%+2.5f ", dvariance[v]);
  printf("\n");
  printf("hidden bias:     ");
  for (h = 0; h < rbm->nhidden; h++) printf("%+2.5f ", dhidden_bias[h]);
  printf("\n");

  printf("weights:\n");
  for (v = 0; v < rbm->nvisible; v++) {
    for (h = 0; h < rbm->nhidden; h++) printf("%+2.5f ", dw[v][h]);
    printf("\n");
  }
  fflush(stdout);

  return (NO_ERROR);
}

double RBMfreeEnergy(RBM *rbm, double *visible)
{
  double free_energy[MAX_RBM_LABELS], total;
  int v, h, l;

  if (rbm->nlabels > 0) {
    for (total = 0.0, l = 0; l < rbm->nlabels; l++) {
      RBMsetLabel(rbm, l);
      RBMactivateForward(rbm, visible);
      RBMactivateBackward(rbm);
      free_energy[l] = 0;
      if (rbm->type == RBM_TYPE_BINARY_INPUTS) {
        for (v = 0; v < rbm->nvisible; v++) free_energy[l] -= visible[v] * rbm->visible_bias[v];
        for (h = 0; h < rbm->nhidden; h++) {
          double pj = rbm->hidden[h];
          if (DZERO(pj)) pj = .000001;

          free_energy[l] -= pj * rbm->act[h];
          free_energy[l] += (pj * log(pj) + (1 - pj) * log(1 - pj));
        }
      }
      else  // Gaussian inputs
      {
        double var;

        for (v = 0; v < rbm->nvisible; v++) {
          var = exp(rbm->variance[v]);
          free_energy[l] += SQR(rbm->visible_bias[v] - visible[v]) / var;
        }
        for (h = 0; h < rbm->nhidden; h++) free_energy[l] -= log(1 + exp(rbm->act[h]));
      }
      total += exp(-free_energy[l]);
    }
    if (DZERO(total)) total = 0.00000001;
    for (l = 0; l < rbm->nlabels; l++) rbm->labels[l] = exp(-free_energy[l]) / total;
  }
  else {
    free_energy[0] = 0;
    for (v = 0; v < rbm->nvisible; v++) free_energy[0] -= visible[v] * rbm->visible_bias[v];
    for (h = 0; h < rbm->nhidden; h++) {
      double pj = rbm->hidden[h];
      if (DZERO(pj)) pj = .000001;

      free_energy[0] -= pj * rbm->act[h];
      free_energy[0] += (pj * log(pj) + (1 - pj) * log(1 - pj));
    }
  }
  return (free_energy[0]);
}

int RBMwriteNetwork(RBM *rbm, int n, RBM_PARMS *parms, int layer)
{
  MRI *mri;
  char fname[STRLEN];

  mri = weights_to_mri(rbm);
  if (layer < 0) {
    if (n < 0) {
      int req = snprintf(fname, STRLEN, "%s.wts.mgz", parms->base_name);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname, STRLEN, "%s.%3.3d.wts.mgz", parms->base_name, n);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
  }
  else {
    if (n < 0) {
      int req = snprintf(fname, STRLEN, "%s.layer%d.wts.mgz", parms->base_name, layer);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname, STRLEN, "%s.%3.3d.layer%d.wts.mgz", parms->base_name, n, layer);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
  }

  printf("saving weights to %s\n", fname);
  MRIwrite(mri, fname);
  MRIfree(&mri);

  return (NO_ERROR);
}

int RBMprintNetworkActivations(RBM *rbm, FILE *fp, int n, RBM_PARMS *parms)
{
  int v, h;

  printf("visible state:  ");
  for (v = 0; v < rbm->nvisible; v++) printf("%2.0f ", rbm->visible[v]);
  printf("\n");
  printf("hidden state:   ");
  for (h = 0; h < rbm->nhidden; h++) printf("%d ", (int)rbm->hidden_state[h]);
  printf("\n");

  return (NO_ERROR);
}

RBM *RBMcopy(RBM *rbm_src, RBM *rbm_dst)
{
  int l, v, h;

  if (rbm_dst == NULL)
    rbm_dst = RBMalloc(rbm_src->type, rbm_src->nvisible, rbm_src->nhidden, rbm_src->nlabels, rbm_src->input_type);

  rbm_dst->ksize = rbm_src->ksize;
  rbm_dst->mri_inputs = rbm_src->mri_inputs;

  for (l = 0; l < rbm_dst->nlabels; l++) {
    rbm_dst->labels[l] = rbm_src->labels[l];
    rbm_dst->lact[l] = rbm_src->lact[l];
    rbm_dst->label_bias[l] = rbm_src->label_bias[l];
    rbm_dst->label_states[l] = rbm_src->label_states[l];
    if (rbm_src->sorted_labels && rbm_dst->sorted_labels) rbm_dst->sorted_labels[l] = rbm_src->sorted_labels[l];
    rbm_dst->label_pvals[l].label = rbm_src->label_pvals[l].label;
    rbm_dst->label_pvals[l].pval = rbm_src->label_pvals[l].pval;
    for (h = 0; h < rbm_dst->nhidden; h++) rbm_dst->label_weights[l][h] = rbm_src->label_weights[l][h];
  }
  for (h = 0; h < rbm_dst->nhidden; h++) {
    rbm_dst->hidden_bias[h] = rbm_src->hidden_bias[h];
    rbm_dst->act[h] = rbm_src->act[h];
    rbm_dst->hidden[h] = rbm_src->hidden[h];
    rbm_dst->hidden_state[h] = rbm_src->hidden_state[h];
    rbm_dst->active_pvals[h] = rbm_src->active_pvals[h];
    rbm_dst->active[h] = rbm_src->active[h];
  }

  for (v = 0; v < rbm_dst->nvisible; v++) {
    rbm_dst->visible[v] = rbm_src->visible[v];
    rbm_dst->visible_bias[v] = rbm_src->visible_bias[v];
    rbm_dst->variance[v] = rbm_src->variance[v];

    for (h = 0; h < rbm_dst->nhidden; h++) rbm_dst->weights[v][h] = rbm_src->weights[v][h];
  }

  return (rbm_dst);
}

RBM *RBMalloc(int type, int nvisible, int nhidden, int nlabels, int input_type)
{
  RBM *rbm;
  int v, h;
  double wt_lim;

  rbm = (RBM *)calloc(1, sizeof(RBM));
  if (rbm == NULL) ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate RBM");

  rbm->input_type = input_type;
  rbm->type = type;
  rbm->nvisible = nvisible;
  rbm->nhidden = nhidden;
  rbm->nlabels = nlabels;

  wt_lim = 1.0 / (rbm->nhidden);
  if (nlabels > 0) {
    rbm->label_bias = (double *)calloc(rbm->nlabels, sizeof(double));
    if (!rbm->label_bias) ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate label array");

    rbm->lact = (double *)calloc(nlabels, sizeof(rbm->labels[0]));
    if (!rbm->lact) ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate label act array");

    rbm->labels = (double *)calloc(nlabels, sizeof(rbm->labels[0]));
    if (!rbm->labels) ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate label  array");
    rbm->label_weights = (double **)calloc(nlabels, sizeof(double *));
    if (!rbm->label_weights) ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate label weight array");

    rbm->label_states = (int *)calloc(nlabels, sizeof(rbm->label_states[0]));
    if (!rbm->label_states) ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate label state array");

    rbm->label_pvals = (LABEL_PVAL *)calloc(nlabels, sizeof(rbm->label_pvals[0]));
    if (!rbm->label_pvals) ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate label pval array");

    for (v = 0; v < nlabels; v++) {
      rbm->label_weights[v] = (double *)calloc(nhidden, sizeof(double));
      if (rbm->label_weights[v] == NULL)
        ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate label weights[%d]\n", v);
      for (h = 0; h < nhidden; h++) {
        rbm->label_weights[v][h] = randomNumber(-2 * wt_lim, wt_lim);
      }
      rbm->label_bias[v] = randomNumber(-.01, .01);
    }
  }

  if (rbm->input_type == RBM_INPUT_VALUE)
    rbm->ksize = 1;  // input is just a single value per input node
  else
    rbm->ksize = (int)nint(sqrt((float)nvisible));  // input is a convoution kernel centered on node
  rbm->visible = (double *)calloc(nvisible, sizeof(double));
  rbm->variance = (double *)calloc(nvisible, sizeof(double));
  rbm->act = (double *)calloc(nhidden, sizeof(double));
  rbm->hidden = (double *)calloc(nhidden, sizeof(double));
  rbm->active = (double *)calloc(nhidden, sizeof(double));
  rbm->active_pvals = (double *)calloc(nhidden, sizeof(double));
  rbm->hidden_state = (double *)calloc(nhidden, sizeof(double));
  if (!rbm->act || !rbm->visible || !rbm->hidden || !rbm->hidden_state || !rbm->active_pvals)
    ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate hidden or input layers");

  rbm->weights = (double **)calloc(nvisible, sizeof(double *));
  if (!rbm->weights) ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate weight array");

  rbm->hidden_bias = (double *)calloc(nhidden, sizeof(double));
  if (!rbm->hidden_bias) ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate hidden bias array");

  rbm->visible_bias = (double *)calloc(nvisible, sizeof(double));
  if (!rbm->visible_bias) ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate visible bias array");

  for (v = 0; v < nvisible; v++) {
    rbm->weights[v] = (double *)calloc(nhidden, sizeof(double));
    if (rbm->weights[v] == NULL) ErrorExit(ERROR_NOMEMORY, "RBMalloc: could not allocate weights[%d]\n", v);
    for (h = 0; h < nhidden; h++) {
      rbm->weights[v][h] = randomNumber(-50 * wt_lim, wt_lim);
      //      rbm->weights[v][h] = randomNumber(-1.2*5*wt_lim, 5*wt_lim) ;
    }
    rbm->visible_bias[v] = randomNumber(-.01, .01);
    if (rbm->type == RBM_TYPE_BINARY_INPUTS)
      rbm->variance[v] = log(1);
    else
      rbm->variance[v] = log(4 * 4);  // zi = log(sigma^2)
  }
  for (h = 0; h < nhidden; h++) rbm->hidden_bias[h] = 0;  // was -4

  return (rbm);
}

int RBMactivateForward(RBM *rbm, double *visible)
{
  int h;

  if (visible) memcpy(rbm->visible, visible, rbm->nvisible * sizeof(visible[0]));

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) shared(rbm, visible) schedule(static, 1)
#endif
  for (h = 0; h < rbm->nhidden; h++) {
    ROMP_PFLB_begin
    
    double act, r, var, delta;
    int v;

    act = rbm->hidden_bias[h];
    for (v = 0; v < rbm->nvisible; v++) {
      if (!devFinite(act)) DiagBreak();

      var = exp(rbm->variance[v]);
      delta = (rbm->visible[v] / var) * rbm->weights[v][h];
      ;
      if (!devFinite(delta)) DiagBreak();
      act += delta;
      if (!devFinite(act)) DiagBreak();
    }
    if (rbm->nlabels > 0) {
      for (v = 0; v < rbm->nlabels; v++) {
        act += rbm->labels[v] * rbm->label_weights[v][h];
        if (!devFinite(act)) DiagBreak();
      }
    }

    rbm->act[h] = act;
    act = rbm->hidden[h] = SIGMOID(act);

    r = randomNumber(0.0, 1.0);
    if (act > r)
      act = 1;
    else
      act = 0;
    rbm->hidden_state[h] = (double)nint(act);
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  return (NO_ERROR);
}

int RBMactivateBackward(RBM *rbm)
{
  int v;

  if (rbm->type == RBM_TYPE_CONTINUOUS_INPUTS) {

    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(experimental) shared(rbm) schedule(static, 1)
#endif
    for (v = 0; v < rbm->nvisible; v++) {
      ROMP_PFLB_begin
      
      double act;
      int h;
      // double var;

      // var = exp(rbm->variance[v]);
      act = rbm->visible_bias[v];
      for (h = 0; h < rbm->nhidden; h++) act += rbm->weights[v][h] * rbm->hidden_state[h];

      rbm->visible[v] = act;
      
      ROMP_PFLB_end
    }
    ROMP_PF_end
  }
  else  // inputs are binary
  {
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(experimental) shared(rbm) schedule(static, 1)
#endif
    for (v = 0; v < rbm->nvisible; v++) {
      ROMP_PFLB_begin
      
      double act, r;
      int h;

      act = rbm->visible_bias[v];
      for (h = 0; h < rbm->nhidden; h++) act += rbm->weights[v][h] * rbm->hidden_state[h];
      act = SIGMOID(act);

      r = randomNumber(0.0, 1.0);
      if (act > r)
        act = 1;
      else
        act = 0;
      rbm->visible[v] = act;
      
      ROMP_PFLB_end
    }
    ROMP_PF_end
  }
  if (rbm->nlabels > 0) {
    double r, total;
    int label;
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(experimental) shared(rbm) schedule(static, 1)
#endif
    for (v = 0; v < rbm->nlabels; v++) {
      ROMP_PFLB_begin
      
      double act;
      int h;

      act = rbm->label_bias[v];
      for (h = 0; h < rbm->nhidden; h++) {
        act += rbm->label_weights[v][h] * rbm->hidden_state[h];
      }
      rbm->lact[v] = act;
      rbm->labels[v] = exp(act);
      if (!devFinite(rbm->labels[v])) rbm->labels[v] = 100;
      
      ROMP_PFLB_end
    }
    ROMP_PF_end
    
    for (total = 0.0, v = 0; v < rbm->nlabels; v++) total += rbm->labels[v];

    if (total > 0)
      for (v = 0; v < rbm->nlabels; v++) {
        rbm->labels[v] /= total;
        if (rbm->labels[v] > 1) DiagBreak();
      }
    RBMsortLabelProbabilities(rbm);
    r = randomNumber(0, 1.0);
    total = 0.0;
    label = 0;
    for (v = 0; v < rbm->nlabels - 1; v++) {
      total += rbm->label_pvals[v].pval;
      if (r < total) {
        label = rbm->label_pvals[v].label;
        break;
      }
    }
    memset(rbm->label_states, 0, sizeof(rbm->label_states[0]) * rbm->nlabels);
    rbm->label_states[label] = 1;
  }
  return (NO_ERROR);
}

int RBMfree(RBM **prbm)
{
  RBM *rbm;
  int v;

  rbm = *prbm;
  *prbm = NULL;

  free(rbm->act);
  free(rbm->hidden_bias);
  free(rbm->visible_bias);
  for (v = 0; v < rbm->nvisible; v++) free(rbm->weights[v]);

  if (rbm->label_states) free(rbm->label_states);
  if (rbm->label_pvals) free(rbm->label_pvals);

  if (rbm->active) free(rbm->active);
  if (rbm->label_weights)
    for (v = 0; v < rbm->nlabels; v++) free(rbm->label_weights[v]);
  if (rbm->labels) free(rbm->labels);
  if (rbm->lact) free(rbm->lact);
  if (rbm->label_bias) free(rbm->label_bias);
  free(rbm->weights);
  free(rbm->visible);
  free(rbm->hidden);
  free(rbm->variance);
  free(rbm->hidden_state);
  free(rbm);
  return (NO_ERROR);
}

int RBMwrite(RBM *rbm, char *fname) { return (NO_ERROR); }

RBM *RBMread(char *fname)
{
  RBM *rbm = NULL;

  return (rbm);
}
static int Ncd = 5;

int RBMfillVisible(RBM *rbm, MRI *mri_inputs, double *visible, int x0, int y0, int z0, int f0, int ksize)
{
  int xk, yk, xi, yi, whalf, v;
  float val;

  if (rbm->input_type == RBM_INPUT_IMAGE) {
    whalf = (ksize - 1) / 2;
    for (v = 0, xk = -whalf; xk <= whalf; xk++) {
      xi = mri_inputs->xi[x0 + xk];
      for (yk = -whalf; yk <= whalf; yk++, v++) {
        yi = mri_inputs->yi[y0 + yk];
        val = MRIgetVoxVal(mri_inputs, xi, yi, 0, f0);
        visible[v] = val;
      }
    }
  }
  else {
    int f;

    for (f = 0; f < mri_inputs->nframes; f++) visible[f] = MRIgetVoxVal(mri_inputs, x0, y0, z0, f);
  }

  return (NO_ERROR);
}

int RBMmostLikelyLabel(RBM *rbm)
{
  int best_label, label;
  double pbest;

  pbest = rbm->labels[0];
  best_label = 0;
  for (label = 1; label < rbm->nlabels; label++)
    if (rbm->labels[label] > pbest) {
      pbest = rbm->labels[label];
      best_label = label;
    }

  return (best_label);
}

int RBMsetLabel(RBM *rbm, int label)
{
  if (label < 0 || label >= rbm->nlabels) return (ERROR_BADPARM);
  memset(rbm->labels, 0, rbm->nlabels * sizeof(rbm->labels[0]));
  rbm->labels[label] = 1;
  return (NO_ERROR);
}

double RBMvoxlistRMS(RBM *rbm, VOXLIST *vl, RBM_PARMS *parms, int *indices, int index, int num)
{
  int i, x, y, z, f, n, ind, v, nvox, h;
  double rms, *visible;
  MRI *mri_inputs;

  mri_inputs = vl->mri;
  visible = (double *)calloc(rbm->nvisible, sizeof(double));
  if (indices) {
    for (h = 0; h < rbm->nhidden; h++) rbm->active[h] = 0;
    for (rms = 0.0, nvox = 0, ind = index; ind < index + num; ind++) {
      i = indices[ind];
      x = vl->xi[i];
      y = vl->yi[i];
      z = vl->zi[i];
      f = vl->fi[i];
      RBMsetLabel(rbm, nint(vl->vsrc[i]));
      RBMfillVisible(rbm, mri_inputs, visible, x, y, z, f, parms->ksize);
      RBMactivateForward(rbm, visible);
      RBMactivateBackward(rbm);
      for (n = 0; n < Ncd; n++) {
        RBMactivateForward(rbm, NULL);
        RBMactivateBackward(rbm);
      }
      for (h = 0; h < rbm->nhidden; h++) rbm->active[h] += rbm->hidden_state[h];
      for (v = 0; v < rbm->nvisible; v++, nvox++) {
        rms += SQR(visible[v] - rbm->visible[v]);
      }
    }
    for (h = 0; h < rbm->nhidden; h++) rbm->active[h] /= (double)num;
  }
  else {
    for (rms = 0.0, nvox = i = 0; i < vl->nvox; i++) {
      x = vl->xi[i];
      y = vl->yi[i];
      z = vl->zi[i];
      f = vl->fi[i];
      RBMsetLabel(rbm, nint(vl->vsrc[i]));
      RBMfillVisible(rbm, mri_inputs, visible, x, y, z, f, parms->ksize);
      RBMactivateForward(rbm, visible);
      RBMactivateBackward(rbm);
      for (n = 0; n < Ncd; n++) {
        RBMactivateForward(rbm, NULL);
        RBMactivateBackward(rbm);
      }
      for (v = 0; v < rbm->nvisible; v++, nvox++) rms += SQR(visible[v] - rbm->visible[v]);
    }
  }
  rms = sqrt(rms / nvox);
  free(visible);
  return (rms);
}

#define HISTO_BINS 100
#include "histo.h"
static double sigma = 10;
double RBMvoxlistHistoRMS(RBM *rbm, VOXLIST *vl, RBM_PARMS *parms, int *indices, int index, int num)
{
  int i, x, y, z, f, n, ind, v, min_ind, max_ind;
  double rms, *visible;
  HISTOGRAM **histo_data, **histo_recon;
  MRI *mri_inputs;
  static int callno = 0;
  char fname[STRLEN];

  mri_inputs = vl->mri;
  histo_recon = (HISTOGRAM **)calloc(rbm->nvisible, sizeof(HISTOGRAM *));
  histo_data = (HISTOGRAM **)calloc(rbm->nvisible, sizeof(HISTOGRAM *));
  for (v = 0; v < rbm->nvisible; v++) {
    histo_recon[v] = HISTOinit(NULL, HISTO_BINS, 0, 1);
    histo_data[v] = HISTOinit(NULL, HISTO_BINS, 0, 1);
  }
  visible = (double *)calloc(rbm->nvisible, sizeof(double));
  if (indices) {
    min_ind = index;
    max_ind = index + num - 1;
  }
  else {
    min_ind = 0;
    max_ind = vl->nvox - 1;
  }

  for (ind = min_ind; ind <= max_ind; ind++) {
    if (indices)
      i = indices[ind];
    else
      i = ind;
    x = vl->xi[i];
    y = vl->yi[i];
    z = vl->zi[i];
    f = vl->fi[i];
    RBMsetLabel(rbm, nint(vl->vsrc[i]));
    RBMfillVisible(rbm, mri_inputs, visible, x, y, z, f, parms->ksize);
    RBMactivateForward(rbm, visible);
    RBMactivateBackward(rbm);
    for (n = 0; n < Ncd; n++) {
      RBMactivateForward(rbm, NULL);
      RBMactivateBackward(rbm);
    }
    for (v = 0; v < rbm->nvisible; v++) {
      HISTOaddSample(histo_data[v], visible[v], 0, 0);
      HISTOaddSample(histo_recon[v], rbm->visible[v], 0, 0);
    }
  }

  for (rms = 0.0, v = 0; v < rbm->nvisible; v++) {
    HISTOsmooth(histo_data[v], histo_data[v], sigma);
    HISTOsmooth(histo_recon[v], histo_recon[v], sigma);
    if (!((callno + 1) % parms->write_iterations)) {
      int req = snprintf(fname, STRLEN, "hist.data.%3.3d.%2.2d.log", callno, v);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      HISTOplot(histo_data[v], fname);
      req = snprintf(fname, STRLEN, "hist.recon.%3.3d.%2.2d.log", callno, v);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      HISTOplot(histo_recon[v], fname);
    }
    rms += HISTOksDistance(histo_data[v], histo_recon[v]);
    //    rms += HISTOrmsDifference(histo_data[v], histo_recon[v]) ;
    HISTOfree(&histo_recon[v]);
    HISTOfree(&histo_data[v]);
  }
  rms /= rbm->nvisible;
  free(visible);
  free(histo_recon);
  free(histo_data);
  callno++;
  return (rms);
}

int RBMcomputeGradients(RBM *rbm,
                        VOXLIST *vl,
                        double **dw,
                        double *dvisible_bias,
                        double *dhidden_bias,
                        double *dvariance,
                        double *dlabel_bias,
                        double **dlabel_weights,
                        RBM_PARMS *parms,
                        int *indices,
                        int index)
{
  int i, x, y, z, f, n, v, h, ind, current_label;
  double *visible, Q0, Qn, V0, Vn, *hidden0, scale, *db_sparsity, *active;
  MRI *mri_inputs = vl->mri;

  for (v = 0; v < rbm->nvisible; v++) memset(dw[v], 0, rbm->nhidden * sizeof(dw[v][0]));
  active = (double *)calloc(rbm->nhidden, sizeof(double));
  hidden0 = (double *)calloc(rbm->nhidden, sizeof(double));
  db_sparsity = (double *)calloc(rbm->nhidden, sizeof(double));
  visible = (double *)calloc(rbm->nvisible, sizeof(double));

  if (rbm->nlabels > 0) {
    memset(dlabel_bias, 0, rbm->nlabels * sizeof(dlabel_bias[0]));
    for (v = 0; v < rbm->nlabels; v++) memset(dlabel_weights[v], 0, rbm->nhidden * sizeof(dlabel_weights[v][0]));
  }
  memset(dvisible_bias, 0, rbm->nvisible * sizeof(dvisible_bias[0]));
  memset(dvariance, 0, rbm->nvisible * sizeof(dvariance[0]));
  memset(dhidden_bias, 0, rbm->nhidden * sizeof(dhidden_bias[0]));
  memset(db_sparsity, 0, rbm->nhidden * sizeof(db_sparsity[0]));

  for (ind = index; ind < index + parms->mini_batch_size; ind++)  // do one mini batch
  {
    i = indices[ind];
    x = vl->xi[i];
    y = vl->yi[i];
    z = vl->zi[i];
    f = vl->fi[i];
    if (y > 0) DiagBreak();
    RBMsetLabel(rbm, current_label = nint(vl->vsrc[i]));
    RBMfillVisible(rbm, mri_inputs, visible, x, y, z, f, parms->ksize);

    RBMactivateForward(rbm, visible);
    memcpy(hidden0, rbm->hidden, rbm->nhidden * sizeof(rbm->hidden[0]));
    RBMactivateBackward(rbm);
    if (parms->debug > 1) {
      printf("index %d: %2.3f\n", i, visible[0]);
      RBMprintNetworkActivations(rbm, stdout, 0, parms);
    }
    for (n = 0; n < Ncd; n++) {
      RBMactivateForward(rbm, NULL);
      RBMactivateBackward(rbm);
    }

    if (parms->debug > 1) {
      printf("final activations\n");
      RBMprintNetworkActivations(rbm, stdout, 0, parms);
    }
    if (rbm->nlabels > 0)  // compute change in label biases and weights
    {
      double l0, ln;

      if (current_label == 3) DiagBreak();
      for (v = 0; v < rbm->nlabels; v++) {
        l0 = (v == current_label);

        ln = rbm->labels[v];
        dlabel_bias[v] += l0 - ln;
        for (h = 0; h < rbm->nhidden; h++) {
          Q0 = hidden0[h];
          Qn = rbm->hidden[h];
          dlabel_weights[v][h] += (Q0 * l0 - Qn * ln);
          dlabel_weights[v][h] -= parms->weight_decays[0] * rbm->label_weights[v][h];
        }
      }
    }

    // compute visible bias and variance change and
    for (v = 0; v < rbm->nvisible; v++) {
      double dvar_data, dvar_model, var;

      var = exp(rbm->variance[v]);
      V0 = visible[v];
      Vn = rbm->visible[v];
      dvisible_bias[v] += (V0 - Vn) / var;

      // compute variance update
      dvar_data = 0.5 * SQR(V0 - rbm->visible_bias[v]);
      dvar_model = 0.5 * SQR(Vn - rbm->visible_bias[v]);

      // compute weight update
      for (h = 0; h < rbm->nhidden; h++) {
        Q0 = hidden0[h];
        Qn = rbm->hidden[h];
        dw[v][h] += (Q0 * V0 - Qn * Vn) / var;
        dw[v][h] -= parms->weight_decays[0] * rbm->weights[v][h];

        dvar_data -= Q0 * rbm->weights[v][h];
        dvar_model -= Qn * rbm->weights[v][h];
      }
      dvar_data *= (V0);
      dvar_model *= (Vn);
      dvariance[v] += (dvar_data - dvar_model);
    }
    // compute hidden bias update
    for (h = 0; h < rbm->nhidden; h++) {
      Q0 = hidden0[h];
      Qn = rbm->hidden[h];
      db_sparsity[h] += Qn;
      dhidden_bias[h] += Q0 - Qn;
      if (rbm->hidden_state[h]) active[h]++;
    }
  }

  scale = 1.0 / parms->mini_batch_size;
  for (v = 0; v < rbm->nvisible; v++) {
    dvisible_bias[v] *= (scale * .1);
    dvariance[v] = exp(-rbm->variance[v]) * dvariance[v] * scale;
    if (dvariance[v] > .1)
      dvariance[v] = .1;
    else if (dvariance[v] < -.1)
      dvariance[v] = -.1;

    for (h = 0; h < rbm->nhidden; h++) dw[v][h] *= scale;
  }
  if (rbm->nlabels > 0) {
    for (v = 0; v < rbm->nlabels; v++) {
      dlabel_bias[v] *= (scale);
      for (h = 0; h < rbm->nhidden; h++) dlabel_weights[v][h] *= (scale);
    }
  }

  for (h = 0; h < rbm->nhidden; h++) {
    double delta;

    dhidden_bias[h] *= (scale * .01);
    db_sparsity[h] *= scale;
    //    dhidden_bias[h] += parms->l_sparsity*(parms->sparsity - db_sparsity[h]) ;

    active[h] *= scale;  // frequency with which this node was active
    rbm->active_pvals[h] = parms->sparsity_decay * rbm->active_pvals[h] + (1 - parms->sparsity_decay) * active[h];
    delta = parms->l_sparsity[0] * (parms->sparsity[0] - rbm->active_pvals[h]);
    dhidden_bias[h] += delta;
    delta /= rbm->nvisible;
    for (v = 0; v < rbm->nvisible; v++) dw[v][h] += delta;
  }

  free(active);
  free(visible);
  free(hidden0);
  free(db_sparsity);
  return (NO_ERROR);
}

int RBMtrainFromImage(RBM *rbm, MRI *mri_inputs, MRI *mri_labels, RBM_PARMS *parms)
{
  VOXLIST *vl;

  rbm->mri_inputs = mri_inputs;
  if (parms->write_iterations > 0) {
    char fname[STRLEN];
    int req = snprintf(fname, STRLEN, "%s.V%3.3d.mgz", parms->base_name, 0);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("writing snapshot to %s\n", fname);
    MRIwrite(mri_inputs, fname);
  }

  vl = VLSTcreate(mri_labels, 1, 255, NULL, 0, 0);
  vl->mri = mri_inputs;
  RBMtrainFromVoxlistImage(rbm, vl, parms);
  return (NO_ERROR);
}

int RBMtrainFromVoxlistImage(RBM *rbm, VOXLIST *vl, RBM_PARMS *parms)
{
  double training_rate, last_rms, pct_diff, rms, momentum, delta, min_rms, **dw, *dvisible_bias, *dhidden_bias,
      *dvariance, *last_dvariance, *dlabel_bias = NULL, **dlabel_weights = NULL, **last_dw, *last_dvisible_bias,
                                   *last_dhidden_bias, *last_dlabel_bias = NULL, **last_dlabel_weights = NULL, var;

  // double sparsity, mean;
  int v, h, step, nbad, *indices, index, b, held_out_index;

  if (!FZERO(parms->variance))
    var = parms->variance;
  else  // estimate it from data
  {
    // mean =
    VLSTmean(vl, NULL, &var);
    var /= parms->nclasses;
    printf("setting initial variances to %2.3f\n", var);
  }
  for (v = 0; v < rbm->nvisible; v++) rbm->variance[v] = log(var);

  printf("training on %d voxels, mini batch size %d (%d), held out %d\n",
         vl->nvox,
         parms->mini_batch_size,
         parms->batches_per_step,
         parms->held_out);
  indices = compute_permutation(vl->nvox, NULL);
  parms->held_out = MIN(vl->nvox - 1, parms->held_out);
  held_out_index = MAX(0, vl->nvox - (parms->held_out));
  //  held_out_index = 0 ;

  dw = (double **)calloc(rbm->nvisible, sizeof(double *));
  last_dw = (double **)calloc(rbm->nvisible, sizeof(double *));
  dvisible_bias = (double *)calloc(rbm->nvisible, sizeof(double));
  dvariance = (double *)calloc(rbm->nvisible, sizeof(double));
  last_dvariance = (double *)calloc(rbm->nvisible, sizeof(double));
  last_dvisible_bias = (double *)calloc(rbm->nvisible, sizeof(double));
  dhidden_bias = (double *)calloc(rbm->nhidden, sizeof(double));
  last_dhidden_bias = (double *)calloc(rbm->nhidden, sizeof(double));
  for (v = 0; v < rbm->nvisible; v++) {
    dw[v] = (double *)calloc(rbm->nhidden, sizeof(double));
    last_dw[v] = (double *)calloc(rbm->nhidden, sizeof(double));
    if (!dw[v] || !last_dw[v]) ErrorExit(ERROR_NOMEMORY, "RBMtrain: could not allocate weight gradients");
  }

  if (rbm->nlabels > 0) {
    last_dlabel_bias = (double *)calloc(rbm->nlabels, sizeof(double));
    dlabel_bias = (double *)calloc(rbm->nlabels, sizeof(double));
    dlabel_weights = (double **)calloc(rbm->nlabels, sizeof(double *));
    last_dlabel_weights = (double **)calloc(rbm->nlabels, sizeof(double *));

    for (v = 0; v < rbm->nlabels; v++) {
      dlabel_weights[v] = (double *)calloc(rbm->nhidden, sizeof(double));
      last_dlabel_weights[v] = (double *)calloc(rbm->nhidden, sizeof(double));
      if (!dlabel_weights[v] || !last_dlabel_weights[v])
        ErrorExit(ERROR_NOMEMORY, "RBMtrain: could not allocate weight gradients");
    }
  }
  // sparsity = parms->sparsity[0];
  training_rate = parms->training_rates[0];
  momentum = parms->momentum[0];

  min_rms = last_rms = RBMvoxlistRMS(rbm, vl, parms, indices, held_out_index, parms->held_out);
  printf("iter %3.3d: rms = %2.3f\n", 0, last_rms);
  for (nbad = step = 0; step < parms->nsteps; step++) {
    for (b = 0; b < parms->batches_per_step; b++) {
      index = (int)nint(randomNumber(0, held_out_index - (parms->mini_batch_size + 1)));
      RBMcomputeGradients(
          rbm, vl, dw, dvisible_bias, dhidden_bias, dvariance, dlabel_bias, dlabel_weights, parms, indices, index);

      if (parms->debug && ((!((step + 1) % parms->write_iterations)) || (step == 0)))
        dump_gradients(rbm, dvisible_bias, dvariance, dhidden_bias, dw, parms, step);
      if (rbm->nlabels > 0) {
        for (v = 0; v < rbm->nlabels; v++) {
          delta = training_rate * dlabel_bias[v] + momentum * last_dlabel_bias[v];
          rbm->label_bias[v] += delta;
          last_dlabel_bias[v] = delta;
          for (h = 0; h < rbm->nhidden; h++) {
            delta = training_rate * dlabel_weights[v][h] + momentum * last_dlabel_weights[v][h];
            rbm->label_weights[v][h] += delta;
            last_dlabel_weights[v][h] = delta;
          }
        }
      }
      for (v = 0; v < rbm->nvisible; v++) {
        delta = training_rate * dvisible_bias[v] + momentum * last_dvisible_bias[v];
        rbm->visible_bias[v] += delta;
        last_dvisible_bias[v] = delta;

        if (parms->learn_variance) {
          delta = training_rate * dvariance[v] + momentum * last_dvariance[v];
          rbm->variance[v] += delta;
          last_dvariance[v] = delta;
        }
        for (h = 0; h < rbm->nhidden; h++) {
          delta = training_rate * dw[v][h] + momentum * last_dw[v][h];
          rbm->weights[v][h] += delta;
          last_dw[v][h] = delta;
        }
      }

      for (h = 0; h < rbm->nhidden; h++) {
        delta = training_rate * dhidden_bias[h] + momentum * last_dhidden_bias[h];
        rbm->hidden_bias[h] += delta;
        last_dhidden_bias[h] = delta;
      }
    }
    rms = RBMvoxlistRMS(rbm, vl, parms, indices, held_out_index, parms->held_out);
    if (!((step + 1) % parms->write_iterations)) RBMwriteNetwork(rbm, step + 1, parms, -1);
    pct_diff = 100 * (last_rms - rms) / (.5 * (last_rms + rms));
    printf("iter %3.3d: rms = %2.5f (%2.5f%%)\n", step + 1, rms, pct_diff);
    if (last_rms < rms) {
      //      training_rate *= .99 ;
      //      printf("error increased - decreasing training rate to %f\n", training_rate) ;
      //      memset(last_dvisible_bias, 0, rbm->nvisible*sizeof(last_dvisible_bias[0])) ;
      //      memset(last_dhidden_bias, 0, rbm->nhidden*sizeof(last_dhidden_bias[0])) ;
      //      for (v = 0 ; v < rbm->nvisible ; v++)
      //	memset(last_dw[v], 0, rbm->nhidden*sizeof(last_dw[v][0])) ;
    }
    else
      training_rate *= 1.000;

    if (rms < min_rms) {
      min_rms = rms;
      nbad = 0;
    }
    else if (nbad++ > parms->max_no_progress) {
      printf("stopping learning due to lack of progress\n");
      break;
    }
    last_rms = rms;
  }

  if (rbm->nlabels > 0) {
    for (v = 0; v < rbm->nlabels; v++) {
      free(dlabel_weights[v]);
      free(last_dlabel_weights[v]);
    }
    free(last_dlabel_bias);
    free(last_dlabel_weights);
  }

  for (v = 0; v < rbm->nvisible; v++) {
    free(dw[v]);
    free(last_dw[v]);
  }
  free(dw);
  free(last_dw);
  free(dvisible_bias);
  free(last_dvisible_bias);
  free(dvariance);
  free(dhidden_bias);
  free(last_dhidden_bias);
  free(last_dvariance);
  return (NO_ERROR);
}
int RBMcountHiddenActive(RBM *rbm)
{
  int hidden_active, h;

  for (h = hidden_active = 0; h < rbm->nhidden; h++) {
    if (rbm->hidden_state[h] == 1) hidden_active++;
  }
  return (hidden_active);
}

MRI *RBMaverageActiveHiddenReceptiveFields(RBM *rbm, MRI *mri_receptive_fields, MRI *mri_inputs, int x0, int y0, int z0)
{
  int k1, k2, xi, yi, zi, h, whalf;
  double val;

  if (mri_receptive_fields == NULL)
    mri_receptive_fields = MRIallocSequence(rbm->ksize, rbm->ksize, 1, MRI_FLOAT, rbm->nhidden);

  whalf = (rbm->ksize - 1) / 2;
  for (h = 0; h < rbm->nhidden; h++) {
    if (rbm->hidden_state[h] == 0) continue;

    zi = z0;
    for (k1 = 0; k1 < rbm->ksize; k1++) {
      xi = mri_inputs->xi[x0 + k1 - whalf];
      for (k2 = 0; k2 < rbm->ksize; k2++) {
        yi = mri_inputs->yi[y0 + k2 - whalf];
        val = MRIgetVoxVal(mri_inputs, xi, yi, zi, 0);
        val += MRIgetVoxVal(mri_receptive_fields, k1, k2, 0, h);
        MRIsetVoxVal(mri_receptive_fields, k1, k2, 0, h, val);
      }
    }
  }
  return (mri_receptive_fields);
}

static double threshold = 0.9;
MRI *RBMreconstruct(RBM *rbm, MRI *mri_inputs, MRI *mri_reconstructed, MRI **pmri_labeled, RBM_PARMS *parms)
{
  int h, x, y, z, center, n, nvox, hidden_active, *hidden_counts, whalf, f;
  float rms, V0, Vn;
  HISTOGRAM *histo;
  HISTOGRAM2D *histo_labels;
  MRI *mri_receptive_fields = NULL, *mri_labeled = NULL;
  char fname[STRLEN];

  histo_labels = HISTO2Dinit(NULL, rbm->nhidden, rbm->nlabels, 0, rbm->nhidden - 1, 0, rbm->nlabels - 1);
  if (mri_reconstructed == NULL) mri_reconstructed = MRIclone(mri_inputs, NULL);

  if (rbm->nlabels > 0) {
    mri_labeled = MRIallocSequence(mri_inputs->width, mri_inputs->height, mri_inputs->depth, MRI_FLOAT, 3);
    MRIcopyHeader(mri_inputs, mri_labeled);
  }

  mri_receptive_fields = MRIallocSequence(rbm->ksize, rbm->ksize, 1, MRI_FLOAT, rbm->nhidden);
  hidden_counts = (int *)calloc(rbm->nhidden, sizeof(int));
  histo = HISTOinit(NULL, rbm->nhidden + 1, 0, rbm->nhidden);
  if (rbm->input_type == RBM_INPUT_IMAGE) {
    center = (parms->ksize * parms->ksize - 1) / 2;
    whalf = (parms->ksize - 1) / 2;
    for (rms = 0.0, f = 0; f < mri_inputs->nframes; f++)
      for (x = whalf; x < mri_inputs->width - whalf; x++) {
        if (!(((x + 1) % 100))) {
          printf("x = %d of %d\n", x, mri_inputs->width);
          MRIwrite(mri_reconstructed, "r.mgz");
        }
        for (y = whalf; y < mri_inputs->height - whalf; y++)
          for (z = 0; z < mri_inputs->depth; z++) {
            if (x == Gx && y == Gy && z == Gz)  // x = 135,y = 681
              DiagBreak();
            V0 = MRIgetVoxVal(mri_inputs, x, y, z, f);

            /*
                        RBMfillVisible(rbm, mri_inputs, rbm->visible, x, y, z, f, parms->ksize);
                        RBMactivateForward(rbm, rbm->visible) ;
                        RBMactivateBackward(rbm) ;
            */
            for (n = 0; n < Ncd; n++) {
              RBMfillVisible(rbm, mri_inputs, rbm->visible, x, y, z, f, parms->ksize);
              RBMactivateForward(rbm, NULL);
              RBMactivateBackward(rbm);
            }
            Vn = rbm->visible[center];
            rms += SQR(Vn - V0);
            MRIsetVoxVal(mri_reconstructed, x, y, z, f, Vn);
            if (rbm->nlabels > 0) {
              int label = RBMmostLikelyLabel(rbm);
              MRIsetVoxVal(mri_labeled, x, y, z, 0, label);
              MRIsetVoxVal(mri_labeled, x, y, z, 1, rbm->labels[label]);
              if (rbm->labels[label] > threshold) MRIsetVoxVal(mri_labeled, x, y, z, 2, label);

              for (h = 0; h < rbm->nhidden; h++)
                if (rbm->hidden_state[h]) HISTO2DaddSample(histo_labels, h, label, -1, -1, -1, -1);
            }
            hidden_active = RBMcountHiddenActive(rbm);
            HISTOaddSample(histo, hidden_active, -1, -1);
            RBMaverageActiveHiddenReceptiveFields(rbm, mri_receptive_fields, mri_inputs, x, y, z);
            for (h = 0; h < rbm->nhidden; h++) hidden_counts[h] += rbm->hidden_state[h];
          }
      }
    int req = snprintf(fname, STRLEN, "%s.hidden.plt", parms->base_name);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("saving %s\n", fname);
    HISTOplot(histo, fname);
    req = snprintf(fname, STRLEN, "%s.hidden_labels.plt", parms->base_name);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("saving %s\n", fname);
    HISTO2Dplot(histo_labels, fname);
  }
  else {
    int v;

    for (rms = 0.0, f = 0; f < mri_inputs->nframes; f++)
      for (x = 0; x < mri_inputs->width; x++) {
        if (!((x + 1) % 100)) {
          printf("x = %d of %d\n", x, mri_inputs->width);
          MRIwrite(mri_reconstructed, "r.mgz");
        }
        for (y = 0; y < mri_inputs->height; y++)
          for (z = 0; z < mri_inputs->depth; z++) {
            RBMfillVisible(rbm, mri_inputs, rbm->visible, x, y, z, f, parms->ksize);
            RBMactivateForward(rbm, rbm->visible);
            RBMactivateBackward(rbm);
            for (n = 0; n < Ncd; n++) {
              RBMactivateForward(rbm, NULL);
              RBMactivateBackward(rbm);
            }
            for (v = 0; v < rbm->nvisible; v++) {
              V0 = MRIgetVoxVal(mri_inputs, x, y, z, v);
              Vn = rbm->visible[v];
              rms += SQR(Vn - V0);
              MRIsetVoxVal(mri_reconstructed, x, y, z, v, Vn);
            }
          }
      }
  }

  for (h = 0; h < rbm->nhidden; h++) {
    if (hidden_counts[h] > 0) {
      float scale = 1.0 / hidden_counts[h];
      MRIscalarMulFrame(mri_receptive_fields, mri_receptive_fields, scale, h);
    }
    else
      printf("hidden node %d never active\n", h);
  }
  MRIwrite(mri_receptive_fields, "rfs.mgz");
  nvox = mri_inputs->width * mri_inputs->height * mri_inputs->depth * mri_inputs->nframes;
  printf("final RMS = %2.3f\n", sqrt(rms / nvox));
  HISTOfree(&histo);
  MRIfree(&mri_receptive_fields);
  HISTO2Dfree(&histo_labels);
  free(hidden_counts);
  if (pmri_labeled) *pmri_labeled = mri_labeled;
  return (mri_reconstructed);
}
DBN *DBNalloc(int type, int nlayers, int nvisible, int *nhidden, int nlabels, int input_type)
{
  DBN *dbn;
  int itype, nl, layer;

  dbn = (DBN *)calloc(1, sizeof(DBN));
  dbn->nlayers = nlayers;
  if (dbn == NULL) ErrorExit(ERROR_NOMEMORY, "DBNalloc: could not allocate DBN");
  dbn->rbms = (RBM **)calloc(dbn->nlayers, sizeof(dbn->rbms[0]));
  if (dbn->rbms == NULL) ErrorExit(ERROR_NOMEMORY, "DBNalloc: could not allocate DBN RBM array");

  for (layer = 0; layer < dbn->nlayers; layer++) {
    if (layer == 0)
      itype = input_type;  // only first layer is continuous (image) inputs
    else {
      nvisible = nhidden[layer - 1];
      type = RBM_TYPE_BINARY_INPUTS;
      itype = RBM_INPUT_VALUE;
    }
    if (layer == dbn->nlayers - 1)
      nl = nlabels;  // only final layer has labels
    else
      nl = 0;
    dbn->rbms[layer] = RBMalloc(type, nvisible, nhidden[layer], nl, itype);
  }

  return (dbn);
}
int DBNfree(DBN **pdbn)
{
  DBN *dbn;
  int layer;

  dbn = *pdbn;
  *pdbn = NULL;
  for (layer = 0; layer < dbn->nlayers; layer++) RBMfree(&dbn->rbms[layer]);
  free(dbn->rbms);
  free(dbn);
  return (NO_ERROR);
}
DBN *DBNread(char *fname) { return (NULL); }
int DBNwrite(DBN *dbn, char *fname) { return (NO_ERROR); }

int DBNactivateForward(DBN *dbn, double *visible, int nlayers)
{
  int l;

  if (nlayers <= 0) nlayers = dbn->nlayers;

  RBMactivateForward(dbn->rbms[0], visible);
  for (l = 1; l < nlayers; l++) RBMactivateForward(dbn->rbms[l], dbn->rbms[l - 1]->hidden_state);

  return (NO_ERROR);
}

int DBNactivateBackward(DBN *dbn, int last_layer, int first_layer)
{
  int l;

  if (last_layer < 0) last_layer = dbn->nlayers - 1;
  if (first_layer < 0) first_layer = 0;

  RBMactivateBackward(dbn->rbms[last_layer]);
  for (l = last_layer; l > first_layer; l--) {
    memcpy(dbn->rbms[l - 1]->hidden_state,
           dbn->rbms[l]->visible,
           dbn->rbms[l]->nvisible * sizeof(dbn->rbms[0]->visible[0]));
    RBMactivateBackward(dbn->rbms[l - 1]);
  }
  return (NO_ERROR);
}

int DBNtrainFromImage(DBN *dbn, MRI *mri_inputs, MRI *mri_labels, RBM_PARMS *parms)
{
  VOXLIST *vl;

  if (parms->write_iterations > 0) {
    char fname[STRLEN];
    int req = snprintf(fname, STRLEN, "%s.V%3.3d.mgz", parms->base_name, 0);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("writing snapshot to %s\n", fname);
    MRIwrite(mri_inputs, fname);
  }

  vl = VLSTcreate(mri_labels, 1, 255, NULL, 0, 0);
  vl->mri = mri_inputs;
  DBNtrainFromVoxlistImage(dbn, vl, parms);
  return (NO_ERROR);
}

int DBNtrainFromVoxlistImage(DBN *dbn, VOXLIST *vl, RBM_PARMS *parms)
{
  double training_rate, last_rms, pct_diff, rms, momentum, delta, min_rms, **dw, *dvisible_bias, *dhidden_bias,
      *dvariance, *last_dvariance, *dlabel_bias = NULL, **dlabel_weights = NULL, **last_dw, *last_dvisible_bias,
                                   *last_dhidden_bias, *last_dlabel_bias = NULL, **last_dlabel_weights = NULL, var;

  // double sparsity, mean;
  int l, v, h, step, nbad, *indices, index, b, held_out_index;
  RBM *rbm;

  if (!FZERO(parms->variance))
    var = parms->variance;
  else  // estimate it from data
  {
    // mean =
    VLSTmean(vl, NULL, &var);
    var /= parms->nclasses;
    printf("setting initial variances to %2.3f\n", var);
  }

  for (v = 0; v < dbn->rbms[0]->nvisible; v++) dbn->rbms[0]->variance[v] = log(var);

  printf("training DBN on %d voxels, mini batch size %d (%d), held out %d\n",
         vl->nvox,
         parms->mini_batch_size,
         parms->batches_per_step,
         parms->held_out);
  indices = compute_permutation(vl->nvox, NULL);
  parms->held_out = MIN(vl->nvox - 1, parms->held_out);
  held_out_index = MAX(0, vl->nvox - (parms->held_out));
  //  held_out_index = 0 ;

  for (l = 0; l < dbn->nlayers; l++) {
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!! TRAINING LAYER %d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", l);

    // sparsity = parms->sparsity[l];
    momentum = parms->momentum[l];
    training_rate = parms->training_rates[l];

    rbm = dbn->rbms[l];
    rbm->mri_inputs = vl->mri;
    dw = (double **)calloc(rbm->nvisible, sizeof(double *));
    last_dw = (double **)calloc(rbm->nvisible, sizeof(double *));
    dvisible_bias = (double *)calloc(rbm->nvisible, sizeof(double));
    dvariance = (double *)calloc(rbm->nvisible, sizeof(double));
    last_dvariance = (double *)calloc(rbm->nvisible, sizeof(double));
    last_dvisible_bias = (double *)calloc(rbm->nvisible, sizeof(double));
    dhidden_bias = (double *)calloc(rbm->nhidden, sizeof(double));
    last_dhidden_bias = (double *)calloc(rbm->nhidden, sizeof(double));
    for (v = 0; v < rbm->nvisible; v++) {
      dw[v] = (double *)calloc(rbm->nhidden, sizeof(double));
      last_dw[v] = (double *)calloc(rbm->nhidden, sizeof(double));
      if (!dw[v] || !last_dw[v]) ErrorExit(ERROR_NOMEMORY, "RBMtrain: could not allocate weight gradients");
    }

    if (rbm->nlabels > 0) {
      last_dlabel_bias = (double *)calloc(rbm->nlabels, sizeof(double));
      dlabel_bias = (double *)calloc(rbm->nlabels, sizeof(double));
      dlabel_weights = (double **)calloc(rbm->nlabels, sizeof(double *));
      last_dlabel_weights = (double **)calloc(rbm->nlabels, sizeof(double *));

      for (v = 0; v < rbm->nlabels; v++) {
        dlabel_weights[v] = (double *)calloc(rbm->nhidden, sizeof(double));
        last_dlabel_weights[v] = (double *)calloc(rbm->nhidden, sizeof(double));
        if (!dlabel_weights[v] || !last_dlabel_weights[v])
          ErrorExit(ERROR_NOMEMORY, "RBMtrain: could not allocate weight gradients");
      }
    }

    min_rms = last_rms = DBNvoxlistRMS(dbn, l, vl, parms, indices, held_out_index, parms->held_out);
    printf("iter %3.3d: rms = %2.3f\n", 0, last_rms);
    for (nbad = step = 0; step < parms->nsteps; step++) {
      for (b = 0; b < parms->batches_per_step; b++) {
        index = (int)nint(randomNumber(0, held_out_index - (parms->mini_batch_size + 1)));
        DBNcomputeGradients(
            dbn, l, vl, dw, dvisible_bias, dhidden_bias, dvariance, dlabel_bias, dlabel_weights, parms, indices, index);

        if (parms->debug && ((!((step + 1) % parms->write_iterations)) || (step == 0)))
          dump_gradients(rbm, dvisible_bias, dvariance, dhidden_bias, dw, parms, step);
        if (rbm->nlabels > 0) {
          for (v = 0; v < rbm->nlabels; v++) {
            delta = training_rate * dlabel_bias[v] + momentum * last_dlabel_bias[v];
            rbm->label_bias[v] += delta;
            last_dlabel_bias[v] = delta;
            for (h = 0; h < rbm->nhidden; h++) {
              delta = training_rate * dlabel_weights[v][h] + momentum * last_dlabel_weights[v][h];
              rbm->label_weights[v][h] += delta;
              last_dlabel_weights[v][h] = delta;
            }
          }
        }
        for (v = 0; v < rbm->nvisible; v++) {
          delta = training_rate * dvisible_bias[v] + momentum * last_dvisible_bias[v];
          rbm->visible_bias[v] += delta;
          last_dvisible_bias[v] = delta;

          if (parms->learn_variance) {
            delta = training_rate * dvariance[v] + momentum * last_dvariance[v];
            rbm->variance[v] += delta;
            last_dvariance[v] = delta;
          }
          for (h = 0; h < rbm->nhidden; h++) {
            delta = training_rate * dw[v][h] + momentum * last_dw[v][h];
            rbm->weights[v][h] += delta;
            last_dw[v][h] = delta;
          }
        }

        for (h = 0; h < rbm->nhidden; h++) {
          delta = training_rate * dhidden_bias[h] + momentum * last_dhidden_bias[h];
          rbm->hidden_bias[h] += delta;
          last_dhidden_bias[h] = delta;
        }
      }
      rms = DBNvoxlistRMS(dbn, l, vl, parms, indices, held_out_index, parms->held_out);
      if (!((step + 1) % parms->write_iterations)) RBMwriteNetwork(rbm, step + 1, parms, l);
      pct_diff = 100 * (last_rms - rms) / (.5 * (last_rms + rms));
      printf("iter %3.3d: rms = %2.5f (%2.5f%%)\n", step + 1, rms, pct_diff);
      if (step == Gdiag_no) DiagBreak();
      if (last_rms < rms) {
        //      training_rate *= .99 ;
        //	printf("error increased - decreasing training rate to %f\n", training_rate) ;
        //      memset(last_dvisible_bias, 0, rbm->nvisible*sizeof(last_dvisible_bias[0])) ;
        //      memset(last_dhidden_bias, 0, rbm->nhidden*sizeof(last_dhidden_bias[0])) ;
        //      for (v = 0 ; v < rbm->nvisible ; v++)
        //	memset(last_dw[v], 0, rbm->nhidden*sizeof(last_dw[v][0])) ;
      }
      else
        training_rate *= 1.000;

      if (rms < min_rms) {
        min_rms = rms;
        nbad = 0;
      }
      else if (nbad++ > parms->max_no_progress) {
        printf("stopping learning due to lack of progress\n");
        break;
      }
      last_rms = rms;
    }

    if (rbm->nlabels > 0) {
      for (v = 0; v < rbm->nlabels; v++) {
        free(dlabel_weights[v]);
        free(last_dlabel_weights[v]);
      }
      free(last_dlabel_bias);
      free(last_dlabel_weights);
    }

    for (v = 0; v < rbm->nvisible; v++) {
      free(dw[v]);
      free(last_dw[v]);
    }
    free(dw);
    free(last_dw);
    free(dvisible_bias);
    free(last_dvisible_bias);
    free(dvariance);
    free(dhidden_bias);
    free(last_dhidden_bias);
    free(last_dvariance);
  }
  return (NO_ERROR);
}
int DBNwriteNetwork(DBN *dbn, int n, RBM_PARMS *parms)
{
  int layer;
  MRI *mri;
  char fname[STRLEN];

  RBMwriteNetwork(dbn->rbms[0], n, parms, 0);
  for (layer = 1; layer < dbn->nlayers; layer++) {
    mri = layer_weights_to_mri(dbn, layer);
    if (n < 0) {
      int req = snprintf(fname, STRLEN, "%s.layer%d.wts.mgz", parms->base_name, layer);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname, STRLEN, "%s.%3.3d.layer%d.wts.mgz", parms->base_name, n, layer);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    printf("saving weights to %s\n", fname);
    MRIwrite(mri, fname);
    MRIfree(&mri);
  }
  return (NO_ERROR);
}

int DBNcomputeGradients(DBN *dbn,
                        int layer,
                        VOXLIST *vl,
                        double **dw,
                        double *dvisible_bias,
                        double *dhidden_bias,
                        double *dvariance,
                        double *dlabel_bias,
                        double **dlabel_weights,
                        RBM_PARMS *parms,
                        int *indices,
                        int index)
{
  int i, x, y, z, f, n, v, h, ind, current_label;
  double *visible, Q0, Qn, V0, Vn, *hidden0, scale, *db_sparsity, *active;
  MRI *mri_inputs = vl->mri;
  RBM *rbm, *rbm_first_layer;

  rbm_first_layer = dbn->rbms[0];
  rbm = dbn->rbms[layer];
  for (v = 0; v < rbm->nvisible; v++) memset(dw[v], 0, rbm->nhidden * sizeof(dw[v][0]));
  hidden0 = (double *)calloc(rbm->nhidden, sizeof(double));
  active = (double *)calloc(rbm->nhidden, sizeof(double));
  db_sparsity = (double *)calloc(rbm->nhidden, sizeof(double));
  visible = (double *)calloc(rbm->nvisible, sizeof(double));

  if (rbm->nlabels > 0) {
    memset(dlabel_bias, 0, rbm->nlabels * sizeof(dlabel_bias[0]));
    for (v = 0; v < rbm->nlabels; v++) memset(dlabel_weights[v], 0, rbm->nhidden * sizeof(dlabel_weights[v][0]));
  }
  memset(dvisible_bias, 0, rbm->nvisible * sizeof(dvisible_bias[0]));
  memset(dvariance, 0, rbm->nvisible * sizeof(dvariance[0]));
  memset(dhidden_bias, 0, rbm->nhidden * sizeof(dhidden_bias[0]));
  memset(db_sparsity, 0, rbm->nhidden * sizeof(db_sparsity[0]));

  for (ind = index; ind < index + parms->mini_batch_size; ind++) {
    i = indices[ind];
    x = vl->xi[i];
    y = vl->yi[i];
    z = vl->zi[i];
    f = vl->fi[i];
    if (y > 0) DiagBreak();
    RBMsetLabel(rbm, current_label = nint(vl->vsrc[i]));
    RBMfillVisible(rbm_first_layer, mri_inputs, rbm_first_layer->visible, x, y, z, f, parms->ksize);
    DBNactivateForward(dbn, rbm_first_layer->visible, layer + 1);

    memcpy(visible, rbm->visible, rbm->nvisible * sizeof(rbm->visible[0]));
    memcpy(hidden0, rbm->hidden, rbm->nhidden * sizeof(rbm->hidden[0]));
    RBMactivateBackward(rbm);
    if (parms->debug > 1) {
      printf("index %d: %2.3f\n", i, visible[0]);
      RBMprintNetworkActivations(rbm, stdout, 0, parms);
    }
    for (n = 0; n < Ncd; n++) {
      RBMactivateForward(rbm, NULL);
      RBMactivateBackward(rbm);
    }

    if (parms->debug > 1) {
      printf("final activations\n");
      RBMprintNetworkActivations(rbm, stdout, 0, parms);
    }
    if (rbm->nlabels > 0)  // compute change in label biases and weights
    {
      double l0, ln;

      if (current_label == 3) DiagBreak();
      for (v = 0; v < rbm->nlabels; v++) {
        l0 = (v == current_label);

        ln = rbm->labels[v];
        dlabel_bias[v] += l0 - ln;
        for (h = 0; h < rbm->nhidden; h++) {
          Q0 = hidden0[h];
          Qn = rbm->hidden[h];
          dlabel_weights[v][h] += (Q0 * l0 - Qn * ln);
          dlabel_weights[v][h] -= parms->weight_decays[layer] * rbm->label_weights[v][h];
        }
      }
    }

    // compute visible bias and variance change and
    for (v = 0; v < rbm->nvisible; v++) {
      double dvar_data, dvar_model, var;

      var = exp(rbm->variance[v]);
      V0 = visible[v];
      Vn = rbm->visible[v];
      dvisible_bias[v] += (V0 - Vn) / var;

      // compute variance update
      dvar_data = 0.5 * SQR(V0 - rbm->visible_bias[v]);
      dvar_model = 0.5 * SQR(Vn - rbm->visible_bias[v]);

      // compute weight update
      for (h = 0; h < rbm->nhidden; h++) {
        Q0 = hidden0[h];
        Qn = rbm->hidden[h];
        dw[v][h] += (Q0 * V0 - Qn * Vn) / var;
        dw[v][h] -= parms->weight_decays[layer] * rbm->weights[v][h];

        dvar_data -= Q0 * rbm->weights[v][h];
        dvar_model -= Qn * rbm->weights[v][h];
      }
      dvar_data *= (V0);
      dvar_model *= (Vn);
      dvariance[v] += (dvar_data - dvar_model);
    }
    // compute hidden bias update
    for (h = 0; h < rbm->nhidden; h++) {
      Q0 = hidden0[h];
      Qn = rbm->hidden[h];
      db_sparsity[h] += Qn;
      dhidden_bias[h] += Q0 - Qn;
      if (rbm->hidden_state[h]) active[h]++;
    }
  }

  scale = 1.0 / parms->mini_batch_size;
  for (v = 0; v < rbm->nvisible; v++) {
    dvisible_bias[v] *= scale;
    dvariance[v] = exp(-rbm->variance[v]) * dvariance[v] * scale;
    if (dvariance[v] > .1)
      dvariance[v] = .1;
    else if (dvariance[v] < -.1)
      dvariance[v] = -.1;

    for (h = 0; h < rbm->nhidden; h++) dw[v][h] *= scale;
  }
  if (rbm->nlabels > 0) {
    for (v = 0; v < rbm->nlabels; v++) {
      dlabel_bias[v] *= (scale * .01);
      for (h = 0; h < rbm->nhidden; h++) dlabel_weights[v][h] *= (scale * .01);
    }
  }

  for (h = 0; h < rbm->nhidden; h++) {
    double delta;

    active[h] *= scale;
    dhidden_bias[h] *= scale;
    db_sparsity[h] *= scale;
    //    dhidden_bias[h] += parms->l_sparsity*(parms->sparsity - db_sparsity[h]) ;
    rbm->active_pvals[h] = parms->sparsity_decay * rbm->active_pvals[h] + (1 - parms->sparsity_decay) * active[h];
    delta = parms->l_sparsity[layer] * (parms->sparsity[layer] - rbm->active_pvals[h]);

    dhidden_bias[h] += delta;
    delta /= rbm->nvisible;
    for (v = 0; v < rbm->nvisible; v++) dw[v][h] += delta;
  }

  free(active);
  free(visible);
  free(hidden0);
  free(db_sparsity);
  return (NO_ERROR);
}

MRI *DBNreconstruct(DBN *dbn, MRI *mri_inputs, MRI *mri_reconstructed, MRI **pmri_labeled, RBM_PARMS *parms)
{
  int x, y, z, center, n, nvox, whalf, f, h, hidden_active;
  float rms = 0.0, V0, Vn;
  MRI *mri_labeled = NULL;
  RBM *rbm_first, *rbm_last;
  HISTOGRAM *histo[MAX_DBN_LAYERS];
  HISTOGRAM2D *histo_labels;
  char fname[STRLEN];

  rbm_first = dbn->rbms[0];
  rbm_last = dbn->rbms[dbn->nlayers - 1];

  histo_labels =
      HISTO2Dinit(NULL, rbm_last->nhidden, rbm_last->nlabels, 0, rbm_last->nhidden - 1, 0, rbm_last->nlabels - 1);
  for (n = 0; n < dbn->nlayers; n++) histo[n] = HISTOinit(NULL, rbm_first->nhidden + 1, 0, rbm_first->nhidden);

  if (mri_reconstructed == NULL) mri_reconstructed = MRIclone(mri_inputs, NULL);

  if (dbn->rbms[dbn->nlayers - 1]->nlabels > 0) {
    mri_labeled = MRIallocSequence(mri_inputs->width, mri_inputs->height, mri_inputs->depth, MRI_FLOAT, 3);
    MRIcopyHeader(mri_inputs, mri_labeled);
  }
  if (rbm_first->input_type == RBM_INPUT_IMAGE) {
    center = (parms->ksize * parms->ksize - 1) / 2;
    whalf = (parms->ksize - 1) / 2;
    for (rms = 0.0, f = 0; f < mri_inputs->nframes; f++)
      for (x = whalf; x < mri_inputs->width - whalf; x++) {
        if (!(((x + 1) % 100))) {
          printf("x = %d of %d\n", x, mri_inputs->width);
          MRIwrite(mri_reconstructed, "r.mgz");
        }
        for (y = whalf; y < mri_inputs->height - whalf; y++)
          for (z = 0; z < mri_inputs->depth; z++) {
            if (x == Gx && y == Gy && z == Gz)  // x = 135,y = 681
              DiagBreak();
            V0 = MRIgetVoxVal(mri_inputs, x, y, z, f);

            RBMfillVisible(rbm_first, mri_inputs, rbm_first->visible, x, y, z, f, parms->ksize);
            DBNactivateForward(dbn, rbm_first->visible, dbn->nlayers);
            DBNactivateBackward(dbn, -1, -1);
            for (n = 0; n < Ncd; n++) {
              RBMfillVisible(rbm_first, mri_inputs, rbm_first->visible, x, y, z, f, parms->ksize);
              DBNactivateForward(dbn, NULL, dbn->nlayers);
              DBNactivateBackward(dbn, -1, -1);
            }
            Vn = rbm_first->visible[center];
            rms += SQR(Vn - V0);
            MRIsetVoxVal(mri_reconstructed, x, y, z, f, Vn);
            if (rbm_last->nlabels > 0) {
              int label = RBMmostLikelyLabel(rbm_last), out_label;
              ;
              out_label = label;
              if (out_label == 2)
                out_label = 4;
              else if (out_label == 3)
                out_label = 0;
              if (out_label == 3) DiagBreak();
              MRIsetVoxVal(mri_labeled, x, y, z, 0, out_label);
              MRIsetVoxVal(mri_labeled, x, y, z, 1, rbm_last->labels[label]);
              if (rbm_last->labels[label] > threshold) MRIsetVoxVal(mri_labeled, x, y, z, 2, out_label);
              for (h = 0; h < rbm_last->nhidden; h++)
                if (rbm_last->hidden_state[h]) HISTO2DaddSample(histo_labels, h, label, -1, -1, -1, -1);
            }
            for (n = 0; n < dbn->nlayers; n++) {
              hidden_active = RBMcountHiddenActive(dbn->rbms[n]);
              HISTOaddSample(histo[n], hidden_active, -1, -1);
            }
          }
      }
    for (n = 0; n < dbn->nlayers; n++) {
      int req = snprintf(fname, STRLEN, "%s.hidden.layer%d.plt", parms->base_name, n);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("saving %s\n", fname);
      HISTOplot(histo[n], fname);
    }
    int req = snprintf(fname, STRLEN, "%s.hidden_labels.plt", parms->base_name);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("saving %s\n", fname);
    HISTO2Dplot(histo_labels, fname);
  }

  nvox = mri_inputs->width * mri_inputs->height * mri_inputs->depth * mri_inputs->nframes;
  printf("final RMS = %2.3f\n", sqrt(rms / nvox));
  if (pmri_labeled) *pmri_labeled = mri_labeled;
  return (mri_reconstructed);
}


double DBNvoxlistRMS(DBN *dbn, int layer, VOXLIST *vl, RBM_PARMS *parms, int *indices, int index, int num)
{
  int i, x, y, z, f, n, ind, v, nvox;
  double rms, *visible;
  MRI *mri_inputs;
  RBM *rbm, *rbm_first_layer;

  if (layer < 0) layer = dbn->nlayers - 1;
  rbm_first_layer = dbn->rbms[0];
  rbm = dbn->rbms[layer];

  mri_inputs = vl->mri;
  visible = (double *)calloc(rbm->nvisible, sizeof(double));
  if (indices) {
    for (rms = 0.0, nvox = 0, ind = index; ind < index + num; ind++) {
      i = indices[ind];
      x = vl->xi[i];
      y = vl->yi[i];
      z = vl->zi[i];
      f = vl->fi[i];
      RBMsetLabel(rbm, nint(vl->vsrc[i]));
      RBMfillVisible(rbm_first_layer, mri_inputs, rbm_first_layer->visible, x, y, z, f, parms->ksize);
      DBNactivateForward(dbn, rbm_first_layer->visible, layer + 1);
      memcpy(visible, rbm->visible, rbm->nvisible * sizeof(rbm->visible[0]));
      RBMactivateBackward(rbm);
      for (n = 0; n < Ncd; n++) {
        RBMactivateForward(rbm, NULL);
        RBMactivateBackward(rbm);
      }
      for (v = 0; v < rbm->nvisible; v++, nvox++) {
        rms += SQR(visible[v] - rbm->visible[v]);
      }
    }
  }
  else {
    for (rms = 0.0, nvox = i = 0; i < vl->nvox; i++) {
      x = vl->xi[i];
      y = vl->yi[i];
      z = vl->zi[i];
      f = vl->fi[i];
      RBMsetLabel(rbm, nint(vl->vsrc[i]));
      RBMfillVisible(rbm, mri_inputs, visible, x, y, z, f, parms->ksize);
      RBMactivateForward(rbm, visible);
      RBMactivateBackward(rbm);
      for (n = 0; n < Ncd; n++) {
        RBMactivateForward(rbm, NULL);
        RBMactivateBackward(rbm);
      }
      for (v = 0; v < rbm->nvisible; v++, nvox++) rms += SQR(visible[v] - rbm->visible[v]);
    }
  }
  rms = sqrt(rms / nvox);
  free(visible);
  return (rms);
}
int CDBNtrainFromImage(CDBN *cdbn, MRI *mri_inputs, MRI *mri_labels, RBM_PARMS *parms)
{
  VOXLIST *vl;

  if (parms->write_iterations > 0) {
    char fname[STRLEN];
    int req = snprintf(fname, STRLEN, "%s.V%3.3d.mgz", parms->base_name, 0);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("writing snapshot to %s\n", fname);
    MRIwrite(mri_inputs, fname);
  }

  vl = VLSTcreate(mri_labels, 1, 255, NULL, 0, 0);
  vl->mri = mri_inputs;
  CDBNtrainFromVoxlistImage(cdbn, vl, parms, mri_inputs, mri_labels);
  return (NO_ERROR);
}

CDBN *CDBNalloc(int type, int nlayers, int *ksizes, int *ngroups, int nlabels, MRI *mri_inputs)
{
  CDBN *cdbn;
  int nl, layer, nvisible;

  cdbn = (CDBN *)calloc(1, sizeof(CDBN));
  cdbn->nlayers = nlayers;
  if (cdbn == NULL) ErrorExit(ERROR_NOMEMORY, "CDBNalloc: could not allocate CDBN");
  cdbn->rbms = (RBM **)calloc(cdbn->nlayers, sizeof(cdbn->rbms[0]));
  if (cdbn->rbms == NULL) ErrorExit(ERROR_NOMEMORY, "CDBNalloc: could not allocate CDBN RBM array");

  cdbn->mri_outputs = (MRI **)calloc(cdbn->nlayers, sizeof(MRI *));
  if (cdbn->mri_outputs == NULL) ErrorExit(ERROR_NOMEMORY, "CDBNalloc: could not allocate CDBN mri_outputs array");

  for (layer = 0; layer < cdbn->nlayers; layer++) {
    cdbn->mri_outputs[layer] =
        MRIallocSequence(mri_inputs->width, mri_inputs->height, mri_inputs->depth, MRI_FLOAT, ngroups[layer]);
    MRIcopyHeader(mri_inputs, cdbn->mri_outputs[layer]);

    nvisible = ksizes[layer] * ksizes[layer];
    if (layer > 0) {
      nvisible *= ngroups[layer - 1];
      type = RBM_TYPE_BINARY_INPUTS;
    }

    if (layer == cdbn->nlayers - 1)
      nl = nlabels;  // only final layer has labels
    else
      nl = 0;
    cdbn->rbms[layer] = RBMalloc(type, nvisible, ngroups[layer], nl, RBM_INPUT_IMAGE);
    cdbn->rbms[layer]->ksize = ksizes[layer];
    cdbn->rbms[layer]->mri_inputs = layer == 0 ? mri_inputs : cdbn->mri_outputs[layer - 1];
  }

  return (cdbn);
}

static int reset_hidden_nodes(CDBN *cdbn, int layer, double min_active, double max_active)
{
  int h, num_off, num_on;
  RBM *rbm;

  rbm = cdbn->rbms[layer];

  for (num_off = num_on = h = 0; h < rbm->nhidden; h++) {
    if (rbm->active[h] < min_active || rbm->active[h] > max_active)  // always on or off
    {
      if (rbm->active[h] < min_active) {
        num_off++;
      } else {
        num_on++;
      }
    }
  }
  if (num_off + num_on > 0)
    printf("resetting %d hidden nodes that are always (%d) or never (%d) active\n", num_off + num_on, num_on, num_off);
  if (0 && num_off > .9 * rbm->nhidden) {
    printf("hidden biases: \n");
    for (h = 0; h < rbm->nhidden; h++) printf(" %2.2f : ", rbm->hidden_bias[h]);
    printf("\n");
  }

  return (num_off + num_on);
}

int CDBNtrainFromVoxlistImage(CDBN *cdbn, VOXLIST *vl, RBM_PARMS *parms, MRI *mri_inputs, MRI *mri_labels)
{
  double sparsity, training_rate, last_rms, pct_diff, rms, momentum, delta, min_rms, **dw, *dvisible_bias,
      *dhidden_bias, *dvariance, *last_dvariance,
      *dlabel_bias = NULL, **dlabel_weights = NULL, **last_dw, *last_dvisible_bias, *last_dhidden_bias,
      *last_dlabel_bias = NULL, **last_dlabel_weights = NULL, mean, var, label_rms, min_label_rms, last_label_rms,
      label_pct_diff;

  // double saved_rms;
  int layer, v, h, step, nbad, *indices, index, b, held_out_index, new_min = 0;
  RBM *rbm, *rbm_min, *rbm_min_label, *rbm_save;
  MRI *mri_layer_inputs;

  if (!FZERO(parms->variance))
    var = parms->variance;
  else  // estimate it from data
  {
    mean = VLSTmean(vl, NULL, &var);
    var /= parms->nclasses;
    //    MRIaddScalar(mri_inputs, mri_inputs, -mean) ;
    printf("setting initial variances to %2.3f and subtracting %2.1f from inputs\n", var, mean);
  }

  for (v = 0; v < cdbn->rbms[0]->nvisible; v++) cdbn->rbms[0]->variance[v] = log(var);

  printf("training CDBN on %d voxels, mini batch size %d (%d), held out %d\n",
         vl->nvox,
         parms->mini_batch_size,
         parms->batches_per_step,
         parms->held_out);
  indices = compute_permutation(vl->nvox, NULL);
  parms->held_out = MIN(vl->nvox - 1, parms->held_out);
  held_out_index = MAX(0, vl->nvox - (parms->held_out));
  //  held_out_index = 0 ;

  for (layer = 0; layer < cdbn->nlayers; layer++) {
    sparsity = parms->sparsity[layer];
    momentum = parms->momentum[layer];
    training_rate = parms->training_rates[layer];
    if (layer > 0)  // create inputs to lth layer from outputs of l-1st
    {
      CDBNcreateOutputs(cdbn, parms, mri_inputs, layer - 1, layer - 1, NULL);
      mri_layer_inputs = cdbn->mri_outputs[layer - 1];
    }
    else
      mri_layer_inputs = mri_inputs;

    printf(
        "************* TRAINING LAYER %d, sparsity = %2.2f wt %2.2f, momentum = %2.1f, dt = %2.4f, weight decay "
        "%2.5f ************\n",
        layer,
        sparsity,
        parms->l_sparsity[layer],
        momentum,
        training_rate,
        parms->weight_decays[layer]);
    rbm = cdbn->rbms[layer];
    rbm_min = RBMcopy(rbm, NULL);
    rbm_min_label = RBMcopy(rbm, NULL);
    rbm->mri_inputs = vl->mri = mri_layer_inputs;
    dw = (double **)calloc(rbm->nvisible, sizeof(double *));
    last_dw = (double **)calloc(rbm->nvisible, sizeof(double *));
    dvisible_bias = (double *)calloc(rbm->nvisible, sizeof(double));
    dvariance = (double *)calloc(rbm->nvisible, sizeof(double));
    last_dvariance = (double *)calloc(rbm->nvisible, sizeof(double));
    last_dvisible_bias = (double *)calloc(rbm->nvisible, sizeof(double));
    dhidden_bias = (double *)calloc(rbm->nhidden, sizeof(double));
    last_dhidden_bias = (double *)calloc(rbm->nhidden, sizeof(double));
    for (v = 0; v < rbm->nvisible; v++) {
      dw[v] = (double *)calloc(rbm->nhidden, sizeof(double));
      last_dw[v] = (double *)calloc(rbm->nhidden, sizeof(double));
      if (!dw[v] || !last_dw[v]) ErrorExit(ERROR_NOMEMORY, "RBMtrain: could not allocate weight gradients");
    }

    if (rbm->nlabels > 0) {
      last_dlabel_bias = (double *)calloc(rbm->nlabels, sizeof(double));
      dlabel_bias = (double *)calloc(rbm->nlabels, sizeof(double));
      dlabel_weights = (double **)calloc(rbm->nlabels, sizeof(double *));
      last_dlabel_weights = (double **)calloc(rbm->nlabels, sizeof(double *));

      for (v = 0; v < rbm->nlabels; v++) {
        dlabel_weights[v] = (double *)calloc(rbm->nhidden, sizeof(double));
        last_dlabel_weights[v] = (double *)calloc(rbm->nhidden, sizeof(double));
        if (!dlabel_weights[v] || !last_dlabel_weights[v])
          ErrorExit(ERROR_NOMEMORY, "RBMtrain: could not allocate weight gradients");
      }
    }

    min_rms = last_rms =
        CDBNvoxlistRMS(cdbn, layer, vl, parms, indices, held_out_index, parms->held_out, &min_label_rms);
    last_label_rms = min_label_rms;
    if (rbm->nlabels > 0)
      printf("iter %3.3d: rms = %2.3f, label_rms = %2.3f\n", 0, last_rms, last_label_rms);
    else
      printf("iter %3.3d: rms = %2.3f\n", 0, last_rms);
    for (nbad = step = 0; step < parms->nsteps; step++) {
      for (b = 0; b < parms->batches_per_step; b++) {
        index = (int)nint(randomNumber(0, held_out_index - (parms->mini_batch_size + 1)));
        CDBNcomputeGradients(cdbn,
                             layer,
                             vl,
                             dw,
                             dvisible_bias,
                             dhidden_bias,
                             dvariance,
                             dlabel_bias,
                             dlabel_weights,
                             parms,
                             indices,
                             index);

        if (parms->debug && ((!((step + 1) % parms->write_iterations)) || (step == 0)))
          dump_gradients(rbm, dvisible_bias, dvariance, dhidden_bias, dw, parms, step);
        if (rbm->nlabels > 0) {
          for (v = 0; v < rbm->nlabels; v++) {
            delta = training_rate * dlabel_bias[v] + momentum * last_dlabel_bias[v];
            rbm->label_bias[v] += delta;
            last_dlabel_bias[v] = delta;
            for (h = 0; h < rbm->nhidden; h++) {
              delta = training_rate * dlabel_weights[v][h] + momentum * last_dlabel_weights[v][h];
              rbm->label_weights[v][h] += delta;
              last_dlabel_weights[v][h] = delta;
            }
          }
        }
        for (v = 0; v < rbm->nvisible; v++) {
          delta = training_rate * dvisible_bias[v] + momentum * last_dvisible_bias[v];
          rbm->visible_bias[v] += delta;
          last_dvisible_bias[v] = delta;
          rbm->visible_bias[v] = rbm->visible_bias[0];  // only 1 visible vias for CDBN
          if (parms->learn_variance) {
            delta = training_rate * dvariance[v] + momentum * last_dvariance[v];
            rbm->variance[v] += delta;
            last_dvariance[v] = delta;
          }
          for (h = 0; h < rbm->nhidden; h++) {
            delta = training_rate * dw[v][h] + momentum * last_dw[v][h];
            rbm->weights[v][h] += delta;
            last_dw[v][h] = delta;
          }
        }

        for (h = 0; h < rbm->nhidden; h++) {
          delta = training_rate * dhidden_bias[h] + momentum * last_dhidden_bias[h];
          rbm->hidden_bias[h] += delta;
          last_dhidden_bias[h] = delta;
        }
      }
      rms = CDBNvoxlistRMS(cdbn, layer, vl, parms, indices, held_out_index, parms->held_out, &label_rms);
      if (step > 0) reset_hidden_nodes(cdbn, layer, .001, .999);
      if (!((step + 1) % parms->write_iterations)) CDBNwriteNetwork(cdbn, step + 1, parms, layer);
      pct_diff = 100 * (last_rms - rms) / (.5 * (last_rms));
      new_min = (rms < min_rms);

      if (rbm->nlabels > 0) {
        label_pct_diff = 100 * (last_label_rms - label_rms) / (.5 * (last_label_rms));
        printf("iter %3.3d: rms = %2.5f (%2.3f%%), label rms = %2.5f (%2.3f%%) %s\n",
               step + 1,
               rms,
               pct_diff,
               label_rms,
               label_pct_diff,
               new_min ? "****" : "");
      }
      else
        printf("iter %3.3d: rms = %2.5f (%2.3f%%) %s\n", step + 1, rms, pct_diff, new_min ? "****" : "");
      if (step == Gdiag_no) DiagBreak();
      if (last_rms < rms) {
        //      training_rate *= .99 ;
        //	printf("error increased - decreasing training rate to %f\n", training_rate) ;
        //      memset(last_dvisible_bias, 0, rbm->nvisible*sizeof(last_dvisible_bias[0])) ;
        //      memset(last_dhidden_bias, 0, rbm->nhidden*sizeof(last_dhidden_bias[0])) ;
        //      for (v = 0 ; v < rbm->nvisible ; v++)
        //	memset(last_dw[v], 0, rbm->nhidden*sizeof(last_dw[v][0])) ;
      }
      else
        training_rate *= 1.000;

      if (label_rms < min_label_rms)  // save best label rms and RBM
      {
        min_label_rms = label_rms;
        nbad = 0;
        rbm_min_label = RBMcopy(rbm, rbm_min_label);
        rbm_save = cdbn->rbms[layer];
        cdbn->rbms[layer] = rbm_min_label;
        rms = CDBNvoxlistRMS(cdbn, layer, vl, parms, indices, held_out_index, parms->held_out, &label_rms);
        cdbn->rbms[layer] = rbm_save;
        if (label_rms < 0.12) DiagBreak();
      }
      if (rms < min_rms)  // save best RBM and it's RMS
      {
        min_rms = rms;
        nbad = 0;
        rbm_min = RBMcopy(rbm, rbm_min);
        rbm_save = cdbn->rbms[layer];
        cdbn->rbms[layer] = rbm_min;
        rms = CDBNvoxlistRMS(cdbn, layer, vl, parms, indices, held_out_index, parms->held_out, &label_rms);
        cdbn->rbms[layer] = rbm_save;
      }
      else if (nbad++ > parms->max_no_progress) {
        printf("stopping learning due to lack of progress\n");
        break;
      }
      last_rms = rms;
      last_label_rms = label_rms;
    }
    // saved_rms = rms;
    rbm_save = RBMcopy(cdbn->rbms[layer], NULL);

    rbm = RBMcopy(rbm_min, cdbn->rbms[layer]);  // restore best one
    if (rbm_min_label) {
      rbm = RBMcopy(rbm_min_label, cdbn->rbms[layer]);  // restore best one
      last_rms = CDBNvoxlistRMS(cdbn, layer, vl, parms, indices, held_out_index, parms->held_out, &last_label_rms);
    }

    if (rbm->nlabels > 0 && parms->discriminative_training) {
      printf("************* STARTING DISCRIMINATIVE TRAINING *****************\n");
      training_rate = parms->label_trate;
      memset(dlabel_bias, 0, rbm->nlabels * sizeof(dlabel_bias[0]));
      memset(last_dlabel_bias, 0, rbm->nlabels * sizeof(dlabel_bias[0]));
      for (v = 0; v < rbm->nlabels; v++) {
        memset(dlabel_weights[v], 0, rbm->nhidden * sizeof(dlabel_weights[v][0]));
        memset(last_dlabel_weights[v], 0, rbm->nhidden * sizeof(dlabel_weights[v][0]));
      }
      memset(dhidden_bias, 0, rbm->nhidden * sizeof(dhidden_bias[0]));
      memset(last_dhidden_bias, 0, rbm->nhidden * sizeof(dhidden_bias[0]));
      for (v = 0; v < rbm->nvisible; v++) {
        memset(dw[v], 0, rbm->nhidden * sizeof(dw[v][0]));
        memset(last_dw[v], 0, rbm->nhidden * sizeof(last_dw[v][0]));
      }
      for (nbad = 0; step < parms->nsteps; step++) {
        for (b = 0; b < parms->batches_per_step; b++) {
          index = (int)nint(randomNumber(0, held_out_index - (parms->mini_batch_size + 1)));
          CDBNcomputeDiscriminativeGradients(
              cdbn, layer, vl, dw, dhidden_bias, dvariance, dlabel_bias, dlabel_weights, parms, indices, index);
          //	  CDBNcomputeLabelGradients(cdbn, layer, vl, dw, dhidden_bias, dvariance, dlabel_bias, dlabel_weights,
          //				    parms, indices, index) ;

          if (parms->debug && ((!((step + 1) % parms->write_iterations)) || (step == 0)))
            dump_gradients(rbm, dvisible_bias, dvariance, dhidden_bias, dw, parms, step);

          for (v = 0; v < rbm->nlabels; v++) {
            delta = training_rate * dlabel_bias[v] + momentum * last_dlabel_bias[v];
            if (fabs(delta) > 1e-5) DiagBreak();
            rbm->label_bias[v] += delta;
            last_dlabel_bias[v] = delta;
            for (h = 0; h < rbm->nhidden; h++) {
              delta = training_rate * dlabel_weights[v][h] + momentum * last_dlabel_weights[v][h];
              if (fabs(delta) > 1e-5) DiagBreak();
              rbm->label_weights[v][h] += delta;
              last_dlabel_weights[v][h] = delta;
            }
          }

          for (v = 0; v < rbm->nvisible; v++) {
            for (h = 0; h < rbm->nhidden; h++) {
              delta = training_rate * dw[v][h] + momentum * last_dw[v][h];
              if (fabs(delta) > 1e-5) DiagBreak();
              rbm->weights[v][h] += delta;
              last_dw[v][h] = delta;
            }
          }

          for (h = 0; h < rbm->nhidden; h++) {
            delta = training_rate * dhidden_bias[h] + momentum * last_dhidden_bias[h];
            if (fabs(delta) > 1e-5) DiagBreak();
            rbm->hidden_bias[h] += delta;
            last_dhidden_bias[h] = delta;
          }
        }
        rms = CDBNvoxlistRMS(cdbn, layer, vl, parms, indices, held_out_index, parms->held_out, &label_rms);
        //	rms = CDBNvoxlistRMS(cdbn, layer, vl, parms, indices, index, parms->mini_batch_size, &label_rms) ;

        if (step > 0) reset_hidden_nodes(cdbn, layer, .001, .999);
        if (!((step + 1) % parms->write_iterations)) CDBNwriteNetwork(cdbn, step + 1, parms, layer);
        pct_diff = 100 * (last_rms - rms) / (.5 * (last_rms));
        label_pct_diff = 100 * (last_label_rms - label_rms) / (.5 * (last_label_rms));
        printf("iter %3.3d: rms = %2.5f (%2.3f%%), label rms = %2.5f (%2.3f%%)\n",
               step + 1,
               rms,
               pct_diff,
               label_rms,
               label_pct_diff);

        if (step == Gdiag_no) DiagBreak();
        if (last_label_rms < label_rms) {
          //      training_rate *= .99 ;
          //	  printf("error increased - decreasing training rate to %f\n", training_rate) ;
          //      memset(last_dvisible_bias, 0, rbm->nvisible*sizeof(last_dvisible_bias[0])) ;
          //      memset(last_dhidden_bias, 0, rbm->nhidden*sizeof(last_dhidden_bias[0])) ;
          //      for (v = 0 ; v < rbm->nvisible ; v++)
          //	memset(last_dw[v], 0, rbm->nhidden*sizeof(last_dw[v][0])) ;
        }
        else
          training_rate *= 1.000;

        if (label_rms < min_label_rms) {
          min_label_rms = label_rms;
          nbad = 0;
        }
        if (rms < min_rms) {
          min_rms = rms;
          nbad = 0;
          rbm_min = RBMcopy(rbm, rbm_min);
        }
        else if (nbad++ > parms->max_no_progress) {
          printf("stopping learning due to lack of progress\n");
          break;
        }
        last_rms = rms;
        last_label_rms = label_rms;
      }
    }

    if (rbm->nlabels > 0) {
      for (v = 0; v < rbm->nlabels; v++) {
        free(dlabel_weights[v]);
        free(last_dlabel_weights[v]);
      }
      free(last_dlabel_bias);
      free(last_dlabel_weights);
    }
    for (v = 0; v < rbm->nvisible; v++) {
      free(dw[v]);
      free(last_dw[v]);
    }
    free(dw);
    free(last_dw);
    free(dvisible_bias);
    free(last_dvisible_bias);
    free(dvariance);
    free(dhidden_bias);
    free(last_dhidden_bias);
    free(last_dvariance);
  }

  return (NO_ERROR);
}

double CDBNvoxlistRMS(
    CDBN *cdbn, int layer, VOXLIST *vl, RBM_PARMS *parms, int *indices, int index, int num, double *plabel_rms)
{
  int i, x, y, z, n, ind, h, v, nvox, current_label, l;
  // int f;
  double rms, *visible, label_rms;
  MRI *mri_inputs;
  RBM *rbm;

  if (layer < 0) layer = cdbn->nlayers - 1;
  rbm = cdbn->rbms[layer];

  mri_inputs = vl->mri;
  visible = (double *)calloc(rbm->nvisible, sizeof(double));
  if (indices) {
    for (h = 0; h < rbm->nhidden; h++) rbm->active[h] = 0;
    for (label_rms = rms = 0.0, nvox = 0, ind = index; ind < index + num; ind++) {
      i = indices[ind];
      x = vl->xi[i];
      y = vl->yi[i];
      z = vl->zi[i];
      // f = vl->fi[i];
      RBMsetLabel(rbm, current_label = nint(vl->vsrc[i]));
      CDBNfillVisible(cdbn, mri_inputs, rbm->visible, x, y, z, rbm->ksize);
      RBMactivateForward(rbm, rbm->visible);
      memcpy(visible, rbm->visible, rbm->nvisible * sizeof(rbm->visible[0]));
      RBMactivateBackward(rbm);
      for (n = 0; n < Ncd; n++) {
        RBMactivateForward(rbm, NULL);
        RBMactivateBackward(rbm);
      }
      for (h = 0; h < rbm->nhidden; h++) rbm->active[h] += rbm->hidden_state[h];
      for (v = 0; v < rbm->nvisible; v++, nvox++) {
        rms += SQR(visible[v] - rbm->visible[v]);
      }
      for (l = 0; l < rbm->nlabels; l++, nvox++) rms += SQR(l == current_label ? 1 : 0 - rbm->labels[l]);
      if (rbm->nlabels > 0) label_rms += SQR(1 - rbm->labels[current_label]);
    }
    for (h = 0; h < rbm->nhidden; h++) rbm->active[h] /= (double)num;
    if (rbm->nlabels > 0) {
      label_rms = sqrt(label_rms / num);
      if (plabel_rms)
        *plabel_rms = label_rms;
      else
        printf("label rms %2.3f\n", label_rms);
    }
  }
  else {
    for (label_rms = rms = 0.0, nvox = i = 0; i < vl->nvox; i++) {
      x = vl->xi[i];
      y = vl->yi[i];
      z = vl->zi[i];
      // f = vl->fi[i];
      RBMsetLabel(rbm, nint(vl->vsrc[i]));
      CDBNfillVisible(cdbn, mri_inputs, visible, x, y, z, parms->ksize);
      RBMactivateForward(rbm, visible);
      RBMactivateBackward(rbm);
      for (n = 0; n < Ncd; n++) {
        RBMactivateForward(rbm, NULL);
        RBMactivateBackward(rbm);
      }
      for (v = 0; v < rbm->nvisible; v++, nvox++) rms += SQR(visible[v] - rbm->visible[v]);
      if (rbm->nlabels > 0) label_rms += SQR(1 - rbm->labels[nint(vl->vsrc[i])]);
    }
    if (rbm->nlabels > 0) {
      label_rms = sqrt(label_rms / num);
      if (plabel_rms)
        *plabel_rms = label_rms;
      else
        printf("label rms %2.3f\n", label_rms);
    }
  }
  rms = sqrt(rms / nvox);
  free(visible);
  return (rms);
}
int CDBNcomputeGradients(CDBN *cdbn,
                         int layer,
                         VOXLIST *vl,
                         double **dw,
                         double *dvisible_bias,
                         double *dhidden_bias,
                         double *dvariance,
                         double *dlabel_bias,
                         double **dlabel_weights,
                         RBM_PARMS *parms,
                         int *indices,
                         int index)
{
  int i, x, y, z, n, v, h, ind, current_label;
  // int f;
  double *visible, Q0, Qn, V0, Vn, *hidden0, scale, *db_sparsity, *active;
  MRI *mri_inputs = vl->mri;
  RBM *rbm;

  rbm = cdbn->rbms[layer];
  for (v = 0; v < rbm->nvisible; v++) memset(dw[v], 0, rbm->nhidden * sizeof(dw[v][0]));
  hidden0 = (double *)calloc(rbm->nhidden, sizeof(double));
  active = (double *)calloc(rbm->nhidden, sizeof(double));
  db_sparsity = (double *)calloc(rbm->nhidden, sizeof(double));
  visible = (double *)calloc(rbm->nvisible, sizeof(double));

  if (rbm->nlabels > 0) {
    memset(dlabel_bias, 0, rbm->nlabels * sizeof(dlabel_bias[0]));
    for (v = 0; v < rbm->nlabels; v++) memset(dlabel_weights[v], 0, rbm->nhidden * sizeof(dlabel_weights[v][0]));
  }
  memset(dvisible_bias, 0, rbm->nvisible * sizeof(dvisible_bias[0]));
  if (parms->learn_variance) memset(dvariance, 0, rbm->nvisible * sizeof(dvariance[0]));
  memset(dhidden_bias, 0, rbm->nhidden * sizeof(dhidden_bias[0]));
  memset(db_sparsity, 0, rbm->nhidden * sizeof(db_sparsity[0]));

  for (ind = index; ind < index + parms->mini_batch_size; ind++) {
    i = indices[ind];
    x = vl->xi[i];
    y = vl->yi[i];
    z = vl->zi[i];
    // f = vl->fi[i];
    if (y > 0) DiagBreak();
    RBMsetLabel(rbm, current_label = nint(vl->vsrc[i]));
    CDBNfillVisible(cdbn, mri_inputs, rbm->visible, x, y, z, rbm->ksize);
    RBMactivateForward(rbm, rbm->visible);

    memcpy(visible, rbm->visible, rbm->nvisible * sizeof(rbm->visible[0]));
    memcpy(hidden0,
           rbm->hidden,
           rbm->nhidden * sizeof(rbm->hidden[0]));  // Bengio says use hidden_state here, but Lee and Ng posterior
    RBMactivateBackward(rbm);
    if (parms->debug > 1) {
      printf("index %d: %2.3f\n", i, visible[0]);
      RBMprintNetworkActivations(rbm, stdout, 0, parms);
    }
    for (n = 0; n < Ncd; n++) {
      RBMactivateForward(rbm, NULL);
      RBMactivateBackward(rbm);
    }

    if (parms->debug > 1) {
      printf("final activations\n");
      RBMprintNetworkActivations(rbm, stdout, 0, parms);
    }
    if (rbm->nlabels > 0)  // compute change in label biases and weights
    {
      double l0, ln;

      if (current_label == 3) DiagBreak();
      for (v = 0; v < rbm->nlabels; v++) {
        l0 = (v == current_label);

        ln = rbm->labels[v];
        dlabel_bias[v] += l0 - ln;

        for (h = 0; h < rbm->nhidden; h++) {
          Q0 = hidden0[h];
          Qn = rbm->hidden[h];
          dlabel_weights[v][h] += (Q0 * l0 - Qn * ln);
          dlabel_weights[v][h] -= parms->weight_decays[layer] * rbm->label_weights[v][h];
        }
      }
    }

    // compute visible bias and variance change and
    for (v = 0; v < rbm->nvisible; v++) {
      double dvar_data = 0, dvar_model = 0, var;

      var = exp(rbm->variance[v]);
      V0 = visible[v];
      Vn = rbm->visible[v];
      //      dvisible_bias[v] += (V0 - Vn)/var ;
      dvisible_bias[0] += (V0 - Vn) / var;

      // compute variance update
      // ATH: is this supposed to be missing?

      // compute weight update
      for (h = 0; h < rbm->nhidden; h++) {
        Q0 = hidden0[h];
        Qn = rbm->hidden[h];
        dw[v][h] += (Q0 * V0 - Qn * Vn) / var;
        dw[v][h] -= parms->weight_decays[layer] * rbm->weights[v][h];

        dvar_data -= Q0 * rbm->weights[v][h];
        dvar_model -= Qn * rbm->weights[v][h];
      }
      if (parms->learn_variance) {
        dvar_data = 0.5 * SQR(V0 - rbm->visible_bias[v]);
        dvar_model = 0.5 * SQR(Vn - rbm->visible_bias[v]);

        for (h = 0; h < rbm->nhidden; h++) {
          Q0 = hidden0[h];
          Qn = rbm->hidden[h];
          dvar_data -= Q0 * rbm->weights[v][h];
          dvar_model -= Qn * rbm->weights[v][h];
        }

        dvar_data *= (V0);
        dvar_model *= (Vn);
        dvariance[v] += (dvar_data - dvar_model);
      }
    }
    // compute hidden bias update
    for (h = 0; h < rbm->nhidden; h++) {
      Q0 = hidden0[h];
      Qn = rbm->hidden[h];
      db_sparsity[h] += Qn;
      dhidden_bias[h] += Q0 - Qn;
      if (rbm->hidden_state[h]) active[h]++;
    }
  }

  scale = 1.0 / parms->mini_batch_size;
  for (v = 0; v < rbm->nvisible; v++) {
    dvisible_bias[v] *= scale;
    if (parms->learn_variance) {
      dvariance[v] = exp(-rbm->variance[v]) * dvariance[v] * scale;
      if (dvariance[v] > .1)
        dvariance[v] = .1;
      else if (dvariance[v] < -.1)
        dvariance[v] = -.1;
    }

    for (h = 0; h < rbm->nhidden; h++) dw[v][h] *= scale;
  }
  if (rbm->nlabels > 0) {
    for (v = 0; v < rbm->nlabels; v++) {
      dlabel_bias[v] *= (scale * 1);                                           // was .01
      for (h = 0; h < rbm->nhidden; h++) dlabel_weights[v][h] *= (scale * 1);  // was .01
    }
  }

  for (h = 0; h < rbm->nhidden; h++) {
    double delta;

    active[h] *= scale;
    dhidden_bias[h] *= scale;
    db_sparsity[h] *= scale;
    //    dhidden_bias[h] += parms->l_sparsity*(parms->sparsity - db_sparsity[h]) ;
    rbm->active_pvals[h] = parms->sparsity_decay * rbm->active_pvals[h] + (1 - parms->sparsity_decay) * active[h];
    delta = parms->l_sparsity[layer] * (parms->sparsity[layer] - rbm->active_pvals[h]);
    dhidden_bias[h] += delta;
    delta /= rbm->nvisible;
    for (v = 0; v < rbm->nvisible; v++) dw[v][h] += delta;
  }

  free(active);
  free(visible);
  free(hidden0);
  free(db_sparsity);
  return (NO_ERROR);
}

static int compute_label_bias = 1;
static int compute_label_weights = 0;
static int compute_hidden_bias = 0;
static int compute_weights = 0;
#if 1
int CDBNcomputeDiscriminativeGradients(CDBN *cdbn,
                                       int layer,
                                       VOXLIST *vl,
                                       double **dw,
                                       double *dhidden_bias,
                                       double *dvariance,
                                       double *dlabel_bias,
                                       double **dlabel_weights,
                                       RBM_PARMS *parms,
                                       int *indices,
                                       int index)
{
  int l, i, x, y, z, n, v, h, ind, current_label;
  // int  f;
  double *visible, scale, var;
  MRI *mri_inputs = vl->mri;
  RBM *rbm;

  rbm = cdbn->rbms[layer];
  memset(dlabel_bias, 0, rbm->nlabels * sizeof(dlabel_bias[0]));
  for (v = 0; v < rbm->nlabels; v++) memset(dlabel_weights[v], 0, rbm->nhidden * sizeof(dlabel_weights[v][0]));
  if (parms->learn_variance) memset(dvariance, 0, rbm->nvisible * sizeof(dvariance[0]));
  memset(dhidden_bias, 0, rbm->nhidden * sizeof(dhidden_bias[0]));

  if (rbm->nlabels == 0) return (NO_ERROR);

  for (v = 0; v < rbm->nvisible; v++) memset(dw[v], 0, rbm->nhidden * sizeof(dw[v][0]));
  visible = (double *)calloc(rbm->nvisible, sizeof(double));

  for (ind = index; ind < index + parms->mini_batch_size; ind++) {
    i = indices[ind];
    x = vl->xi[i];
    y = vl->yi[i];
    z = vl->zi[i];
    // f = vl->fi[i];
    if (y > 0) DiagBreak();

    current_label = nint(vl->vsrc[i]);
    CDBNfillVisible(cdbn, mri_inputs, visible, x, y, z, rbm->ksize);

    RBMactivateForward(rbm, visible);
    RBMactivateBackward(rbm);
    for (n = 0; n < Ncd * 0; n++) {
      RBMactivateForward(rbm, NULL);
      RBMactivateBackward(rbm);
    }

    for (l = 0; l < rbm->nlabels; l++) {
      if (current_label == Gdiag_no) DiagBreak();

      dlabel_bias[l] = (l == current_label) - rbm->labels[l];

      // compute parameter updates
      for (h = 0; h < rbm->nhidden; h++) {
        if (compute_hidden_bias) {
          dhidden_bias[h] -= rbm->hidden[h] * rbm->labels[l];
        }
        if (compute_weights) {
          for (v = 0; v < rbm->nvisible; v++) {
            var = exp(rbm->variance[v]);
            dw[v][h] -= rbm->hidden[h] * rbm->labels[l] * visible[v] / var;
          }
        }
        if (compute_label_weights) dlabel_weights[l][h] -= rbm->hidden[h] * rbm->labels[l];
        if (l == current_label)  // add in positive term
        {
          if (compute_hidden_bias) {
            dhidden_bias[h] += rbm->hidden[h];
          }
          if (compute_weights) {
            for (v = 0; v < rbm->nvisible; v++) {
              var = rbm->variance[v];
              dw[v][h] += rbm->hidden[h] * visible[v] / var;
            }
          }
          if (compute_label_weights) dlabel_weights[l][h] += rbm->hidden[h];
        }
      }
    }
  }

  scale = 1.0 / parms->mini_batch_size;
  for (v = 0; v < rbm->nvisible; v++) {
    for (h = 0; h < rbm->nhidden; h++) dw[v][h] *= scale;
  }
  for (l = 0; l < rbm->nlabels; l++) {
    dlabel_bias[l] *= (scale * 1);                                           // was .01
    for (h = 0; h < rbm->nhidden; h++) dlabel_weights[l][h] *= (scale * 1);  // was .01
  }

  // sparsity stuff
  for (h = 0; h < rbm->nhidden; h++) {
    //    double delta ;

    dhidden_bias[h] *= scale;
//    rbm->active_pvals[h] = parms->sparsity_decay*rbm->active_pvals[h] + (1-parms->sparsity_decay)*active[h];
  }

  free(visible);
  return (NO_ERROR);
}

// compute the gradient of E = sum(SQR(true label - prob of true label)^2)
int CDBNcomputeLabelGradients(CDBN *cdbn,
                              int layer,
                              VOXLIST *vl,
                              double **dw,
                              double *dhidden_bias,
                              double *dvariance,
                              double *dlabel_bias,
                              double **dlabel_weights,
                              RBM_PARMS *parms,
                              int *indices,
                              int index)
{
  int l, i, x, y, z, n, v, h, ind, current_label;
  // int f;
  double *visible, scale, var, softmax_deriv[MAX_RBM_LABELS];
  MRI *mri_inputs = vl->mri;
  RBM *rbm;

  rbm = cdbn->rbms[layer];
  memset(dlabel_bias, 0, rbm->nlabels * sizeof(dlabel_bias[0]));
  for (v = 0; v < rbm->nlabels; v++) memset(dlabel_weights[v], 0, rbm->nhidden * sizeof(dlabel_weights[v][0]));
  if (parms->learn_variance) memset(dvariance, 0, rbm->nvisible * sizeof(dvariance[0]));
  memset(dhidden_bias, 0, rbm->nhidden * sizeof(dhidden_bias[0]));

  if (rbm->nlabels == 0) return (NO_ERROR);

  for (v = 0; v < rbm->nvisible; v++) memset(dw[v], 0, rbm->nhidden * sizeof(dw[v][0]));
  visible = (double *)calloc(rbm->nvisible, sizeof(double));

  for (ind = index; ind < index + parms->mini_batch_size; ind++) {
    i = indices[ind];
    x = vl->xi[i];
    y = vl->yi[i];
    z = vl->zi[i];
    // f = vl->fi[i];
    if (y > 0) DiagBreak();
    RBMsetLabel(rbm, current_label = nint(vl->vsrc[i]));
    CDBNfillVisible(cdbn, mri_inputs, visible, x, y, z, rbm->ksize);

    RBMactivateForward(rbm, visible);
    RBMactivateBackward(rbm);
    if (parms->debug > 1) {
      printf("index %d: %2.3f\n", i, visible[0]);
      RBMprintNetworkActivations(rbm, stdout, 0, parms);
    }
    for (n = 0; n < Ncd * 0; n++) {
      RBMactivateForward(rbm, NULL);
      RBMactivateBackward(rbm);
    }

    if (parms->debug > 1) {
      printf("final activations\n");
      RBMprintNetworkActivations(rbm, stdout, 0, parms);
    }

    if (current_label == 3) DiagBreak();

    // compute derivative of softmax term
    for (l = 0; l < rbm->nlabels; l++) {
      int l1;
      double num, denom;

      softmax_deriv[l] = 0;
      for (num = denom = 0.0, l1 = 0; l1 < rbm->nlabels; l1++) {
        if (l1 != l) num += exp(rbm->lact[l1]);
        denom += exp(rbm->lact[l1]);
      }
      softmax_deriv[l] = exp(rbm->act[l] * num) / (denom * denom);
    }
    // compute hidden weight and bias update
    for (h = 0; h < rbm->nhidden; h++) {
      int l0;
      for (l = 0; l < rbm->nlabels; l++) {
        l0 = (l == current_label);

        if (compute_hidden_bias) {
          dhidden_bias[h] += (l0 - rbm->labels[l]) * rbm->label_weights[l][h] * rbm->hidden[h] * (1 - rbm->hidden[h]);
          if (!devFinite(dhidden_bias[h])) DiagBreak();
        }
        if (compute_weights) {
          for (v = 0; v < rbm->nvisible; v++) {
            var = exp(rbm->variance[v]);
            dw[v][h] += (l0 - rbm->labels[l]) * rbm->label_weights[l][h] * rbm->hidden[h] * (1 - rbm->hidden[h]) *
                        visible[v] / var;
          }
        }
      }
    }
    // compute label bias update
    for (l = 0; l < rbm->nlabels; l++) {
      int l0;
      l0 = (l == current_label);

      if (compute_label_bias) dlabel_bias[l] += (l0 - rbm->labels[l]) * softmax_deriv[l];
      if (!devFinite(dlabel_bias[h])) DiagBreak();
      if (compute_label_weights) {
        for (h = 0; h < rbm->nhidden; h++)
          dlabel_weights[l][h] += (l0 - rbm->labels[l]) * (rbm->hidden[h] + rbm->label_weights[l][h] * l0);
      }
    }
  }

  scale = 1.0 / parms->mini_batch_size;
  for (v = 0; v < rbm->nvisible; v++) {
    for (h = 0; h < rbm->nhidden; h++) dw[v][h] *= scale;
  }
  for (l = 0; l < rbm->nlabels; l++) {
    dlabel_bias[l] *= (scale * 1);                                           // was .01
    for (h = 0; h < rbm->nhidden; h++) dlabel_weights[l][h] *= (scale * 1);  // was .01
  }

  // sparsity stuff
  for (h = 0; h < rbm->nhidden; h++) {
    //    double delta ;

    dhidden_bias[h] *= scale;
  }

  free(visible);
  return (NO_ERROR);
}
#endif
int CDBNactivateForward(CDBN *cdbn, double *visible, int nlayers)
{
  int l;

  if (nlayers <= 0) nlayers = cdbn->nlayers;

  RBMactivateForward(cdbn->rbms[0], visible);
  for (l = 1; l < nlayers; l++) RBMactivateForward(cdbn->rbms[l], cdbn->rbms[l - 1]->hidden_state);

  return (NO_ERROR);
}

int CDBNactivateBackward(CDBN *cdbn, int last_layer, int first_layer)
{
  int l;

  if (last_layer < 0) last_layer = cdbn->nlayers - 1;
  if (first_layer < 0) first_layer = 0;

  RBMactivateBackward(cdbn->rbms[last_layer]);
  for (l = last_layer; l > first_layer; l--) {
    memcpy(cdbn->rbms[l - 1]->hidden_state,
           cdbn->rbms[l]->visible,
           cdbn->rbms[l]->nvisible * sizeof(cdbn->rbms[0]->visible[0]));
    RBMactivateBackward(cdbn->rbms[l - 1]);
  }
  return (NO_ERROR);
}

int CDBNfillVisible(CDBN *cdbn, MRI *mri_inputs, double *visible, int x0, int y0, int z0, int ksize)
{
  int xk, yk, xi, yi, whalf, v, f;
  float val;

  whalf = (ksize - 1) / 2;
  for (v = f = 0; f < mri_inputs->nframes; f++)
    for (xk = -whalf; xk <= whalf; xk++) {
      xi = mri_inputs->xi[x0 + xk];
      for (yk = -whalf; yk <= whalf; yk++, v++) {
        yi = mri_inputs->yi[y0 + yk];
        val = MRIgetVoxVal(mri_inputs, xi, yi, 0, f);
        visible[v] = val;
      }
    }

  return (NO_ERROR);
}
MRI *CDBNcreateOutputs(
    CDBN *cdbn, RBM_PARMS *parms, MRI *mri_inputs, int first_layer, int last_layer, MRI **pmri_labeled)
{
  int layer, x, y, z, h, *hidden_counts;
  MRI *mri_layer_inputs, *mri_outputs = NULL, *mri_labeled = NULL;
  RBM *rbm;
  static int callno = 0;
  callno++;
  if (cdbn->rbms[cdbn->nlayers - 1]->nlabels > 0 && pmri_labeled) {
    mri_labeled = MRIallocSequence(mri_inputs->width, mri_inputs->height, mri_inputs->depth, MRI_FLOAT, 5);
    MRIcopyHeader(mri_inputs, mri_labeled);
  }
  for (layer = first_layer, mri_layer_inputs = mri_inputs; layer <= last_layer; layer++) {
    rbm = cdbn->rbms[layer];
    hidden_counts = (int *)calloc(rbm->nhidden, sizeof(int));
    mri_outputs = cdbn->mri_outputs[layer];
    for (x = 0; x < mri_inputs->width; x++) {
      if (x && !(x % 100)) printf("layer %d of %d: x = %d of %d\n", layer, last_layer, x, mri_inputs->width);
      for (y = 0; y < mri_inputs->height; y++)
        for (z = 0; z < mri_inputs->depth; z++) {
          if (x == Gx && y == Gy) DiagBreak();
          if (rbm->nlabels > 0) memset(rbm->labels, 0, rbm->nlabels * sizeof(rbm->labels[0]));
          CDBNfillVisible(cdbn, mri_layer_inputs, rbm->visible, x, y, z, rbm->ksize);
          RBMactivateForward(rbm, NULL);
          RBMactivateBackward(rbm);
          if (pmri_labeled) CDBNfillVisible(cdbn, mri_layer_inputs, rbm->visible, x, y, z, rbm->ksize);
          RBMactivateForward(rbm, NULL);
          for (h = 0; h < rbm->nhidden; h++) {
            MRIsetVoxVal(mri_outputs, x, y, z, h, rbm->hidden_state[h]);
            hidden_counts[h] += rbm->hidden_state[h];
          }
          if ((rbm->nlabels > 0) && pmri_labeled) {
            int label = RBMmostLikelyLabel(rbm), out_label;
            out_label = label;
            if (out_label == 2)
              out_label = 4;
            else if (out_label == 3)
              out_label = 0;
            if (out_label == 3) DiagBreak();
            MRIsetVoxVal(mri_labeled, x, y, z, 0, out_label);
            MRIsetVoxVal(mri_labeled, x, y, z, 1, rbm->labels[label]);
            if (rbm->labels[label] > threshold) MRIsetVoxVal(mri_labeled, x, y, z, 2, out_label);

            CDBNfillVisible(cdbn, mri_layer_inputs, rbm->visible, x, y, z, rbm->ksize);
            RBMfreeEnergy(rbm, rbm->visible);
            label = RBMmostLikelyLabel(rbm);
            out_label = label;
            if (out_label == 2)
              out_label = 4;
            else if (out_label == 3)
              out_label = 0;
            if (out_label == 3) DiagBreak();
            MRIsetVoxVal(mri_labeled, x, y, z, 3, rbm->labels[label]);
            MRIsetVoxVal(mri_labeled, x, y, z, 4, out_label);
          }
        }
    }
    mri_layer_inputs = mri_outputs;
    {
      char fname[STRLEN];
      FILE *fp;
      int nvox, always, never;

      int req = snprintf(fname, STRLEN, "%s.%3.3d.layer%d.hidden_counts.txt",
			 parms->base_name, callno, layer);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fp = fopen(fname, "w");
      always = never = 0;
      nvox = mri_inputs->height * mri_inputs->width * mri_inputs->depth;
      for (h = 0; h < rbm->nhidden; h++) {
        if (!hidden_counts[h])
          never++;
        else if (hidden_counts[h] == nvox)
          always++;
        fprintf(fp, "%d %d\n", h, hidden_counts[h]);
      }
      fclose(fp);
      printf("writing hidden counts to %s (%d/%2.1f%% always active, %d/%2.1f%% never)\n",
             fname,
             always,
             100.0 * always / rbm->nhidden,
             never,
             100.0f * never / rbm->nhidden);
    }

    free(hidden_counts);
  }
  if (pmri_labeled) *pmri_labeled = mri_labeled;

  return (mri_outputs);
}
MRI *cdbn_layer_weights(CDBN *cdbn, int layer)
{
  MRI *mri = NULL, *mri_prev, *mri_counts;
  int width, whalf_prev, x, y, xk, yk, v, h, count, xp, yp, hp;
  RBM *rbm, *rbm_prev;
  float val, val_prev;

  if (layer <= 0)
    return (weights_to_mri(cdbn->rbms[0]));
  else if (layer >= 1) {
    rbm = cdbn->rbms[layer];
    rbm_prev = cdbn->rbms[layer - 1];
    mri_prev = cdbn_layer_weights(cdbn, layer - 1);
    whalf_prev = (mri_prev->width - 1) / 2;
    width = rbm->ksize + 2 * whalf_prev;

    mri = MRIallocSequence(width, width, 1, MRI_FLOAT, rbm->nhidden);
    MRIcopyHeader(cdbn->rbms[0]->mri_inputs, mri);
    mri_counts = MRIallocSequence(width, width, 1, MRI_INT, 1);
    MRIcopyHeader(mri, mri_counts);

    // v is the visible unit in this layer, which is the hidden unit (or frame) in the previous one
    for (h = 0; h < rbm->nhidden; h++) {
      for (v = hp = 0; hp < rbm_prev->nhidden; hp++) {
        for (x = 0; x < rbm->ksize; x++) {
          for (y = 0; y < rbm->ksize; y++, v++) {
            for (xp = 0; xp < mri_prev->width; xp++) {
              xk = mri->xi[x + xp];  // this is the location in the larger kernel
              for (yp = 0; yp < mri_prev->height; yp++) {
                yk = mri->yi[y + yp];  // y position in larger kernel
                if (xk == Gx && yk == Gy) DiagBreak();
                if (h == 0 && x == rbm->ksize - 1 && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON && 0)
                  printf("x = %d, xp = %d, xi = %d    y = %d, yp = %d, yk = %d, v = %d\n", x, xp, xk, y, yp, yk, v);
                count = MRIgetVoxVal(mri_counts, xk, yk, 0, 0);
                count++;
                MRIsetVoxVal(mri_counts, xk, yk, 0, 0, count);
                val = MRIgetVoxVal(mri, xk, yk, 0, h);
                val_prev = MRIgetVoxVal(mri_prev, xp, yp, 0, hp);
                val_prev *= rbm->weights[v][h];
                MRIsetVoxVal(mri, xk, yk, 0, h, val + val_prev);
              }
            }  // end of interior kernel
            DiagBreak();
          }
        }  // end of current kernel
        DiagBreak();
      }  // end of all inputs to this node
      for (x = 0; x < mri_counts->width; x++)
        for (y = 0; y < mri_counts->height; y++) {
          count = MRIgetVoxVal(mri_counts, x, y, 0, 0);
          if (count > 0) {
            val = MRIgetVoxVal(mri, x, y, 0, h);
            val /= (float)count;
            MRIsetVoxVal(mri, x, y, 0, h, val);
          }
        }
      MRIclear(mri_counts);
    }
    MRIfree(&mri_prev);
  }
  return (mri);
}

int CDBNwriteNetwork(CDBN *cdbn, int n, RBM_PARMS *parms, int layer)
{
  MRI *mri;
  char fname[STRLEN];

  mri = cdbn_layer_weights(cdbn, layer);
  if (n < 0) {
    int req = snprintf(fname, STRLEN, "%s.layer%d.wts.mgz", parms->base_name, layer);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  } else {
    int req = snprintf(fname, STRLEN, "%s.%3.3d.layer%d.wts.mgz", parms->base_name, n, layer);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  }
  printf("saving weights to %s\n", fname);
  MRIwrite(mri, fname);
  MRIfree(&mri);
  return (NO_ERROR);
}

int dump_visible(RBM *rbm, char *fname)
{
  FILE *fp = fopen(fname, "w");
  int v;

  for (v = 0; v < rbm->nvisible; v++) fprintf(fp, "%f\n", rbm->visible[v]);

  fclose(fp);
  return (NO_ERROR);
}
int dump_hidden(RBM *rbm, char *fname)
{
  FILE *fp = fopen(fname, "w");
  int v;

  for (v = 0; v < rbm->nhidden; v++) fprintf(fp, "%f\n", rbm->hidden[v]);

  fclose(fp);
  return (NO_ERROR);
}
