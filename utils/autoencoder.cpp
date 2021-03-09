/**
 * @brief header file for creating and training a stacked autoencoder for
feature extraction.
 *
H.-C. Shin, M. R. Orton, D. J. Collins, S. J. Doran, and M. O. Leach,
"Stacked Autoencoders for
Unsupervised Feature Learning and Multiple Organ Detectionin a Pilot Study
Using 4D Patient Data,"
IEEE Transaction on Pattern Analysis and Machine Intelligence, 2012.

 */
/*
 * Original Author: Bruce Fischl
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

#include <stdio.h>
#include "romp_support.h"

#include "diag.h"
#include "error.h"
#include "fio.h"
#include "macros.h"
#include "utils.h"
#include "voxlist.h"

#include "autoencoder.h"

static double AEcomputeHiddenRMS(AE *ae, SAE_INTEGRATION_PARMS *parms);
static int aeApplyGradient(AE *ae, SAE_INTEGRATION_PARMS *parms, double dt);
static double AEcomputeRMS(AE *ae);
static int AEsaveState(AE *ae);
static int AErestoreState(AE *ae);
static AE *AEalloc(AE *prev, int ninputs, int nhidden, int noutputs);
static int AEactivateLayer(AE *ae, VECTOR *v_input);
static double AEtrain(AE *ae, SAE_INTEGRATION_PARMS *parms);
static double AEaccumulateGradient(AE *ae, SAE_INTEGRATION_PARMS *parms);
static int AEclearGradients(AE *ae);
static int AEwrite(AE *ae, FILE *fp);
static AE *AEread(FILE *fp, AE *prev);
static int CSAEcomputeGradient(CSAE *csae,
                               AE *ae,
                               VOXEL_LIST *vl,
                               MRI *mri,
                               int start_index,
                               int end_index,
                               int *indices,
                               SAE_INTEGRATION_PARMS *parms);
static double aeApplyAccumulatedGradient(AE *ae, SAE_INTEGRATION_PARMS *parms);
static int reset_constant_nodes(AE *ae, double thresh);

AE *SAEfindLastLayer(SAE *sae, AE *ae)
{
  if (ae == NULL)
    return (SAEfindLastLayer(sae, sae->first));
  else if (ae->next)
    return (SAEfindLastLayer(sae, ae->next));
  else
    return (ae);
}

SAE *SAEalloc(int whalf, int nlevels, int type, double scale)
{
  SAE *sae;
  int ninputs, wsize, nhidden, noutputs;

  wsize = 2 * whalf + 1;
  ninputs = wsize * wsize * nlevels;
  if (!(type & AUTOENCODER_2D)) ninputs *= wsize;  // 3D
  nhidden = nint(scale * ninputs);
  printf("allocating SAE with inputs/hidden %d/%d\n", ninputs, nhidden);

  sae = (SAE *)calloc(1, sizeof(SAE));
  if (sae == NULL) ErrorReturn(NULL, (ERROR_NOFILE, "SAEalloc(%s, %d): could not alloc sae", whalf));
  sae->whalf = whalf;
  sae->nencoders = 1;
  sae->nlevels = nlevels;
  sae->type = type;
  noutputs = type & FOCUSED_AUTOENCODER ? 1 : ninputs;
  sae->first = AEalloc(NULL, ninputs, nhidden, noutputs);
  sae->scale = scale;
  if (whalf == 0) {
    sae->first->v_input->rptr[1][1] = 0;
    SAEactivateNetwork(sae);
    AEdump(sae->first);
    sae->first->v_input->rptr[1][1] = 1;
    SAEactivateNetwork(sae);
    AEdump(sae->first);
  }
  return (sae);
}

static AE *AEalloc(AE *prev, int ninputs, int nhidden, int noutputs)
{
  AE *ae;
  int i, j, k;
  double norm, w, wt_lim = .1;

  ae = (AE *)calloc(1, sizeof(AE));
  if (ae == NULL) ErrorReturn(NULL, (ERROR_NOFILE, "AEaalloc(%d, %d): could not alloc ae", ninputs, nhidden));

  // the input and ouptut layers for an AE is the hidden layer of the previous
  // AE
  if (prev)
    ae->v_output = prev->v_hidden;
  else {
    ae->v_output = VectorAlloc(noutputs, MATRIX_REAL);
    if (ae->v_output == NULL)
      ErrorReturn(NULL, (ERROR_NOFILE, "AEaalloc(%d, %d): could not alloc v_output", ninputs, nhidden));
  }

  ae->v_input = VectorAlloc(ninputs, MATRIX_REAL);  // will be a copy of previous layer to allow calculation of errors
  if (ae->v_input == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "AEaalloc(%d, %d): could not alloc v_input", ninputs, nhidden));

  ae->m_input_to_hidden = MatrixAlloc(nhidden, ninputs, MATRIX_REAL);
  if (ae->m_input_to_hidden == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "AEaalloc(%d, %d): could not alloc m_input_to_hidden", ninputs, nhidden));
  ae->m_hidden_to_output = MatrixAlloc(noutputs, nhidden, MATRIX_REAL);
  if (ae->m_hidden_to_output == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "AEaalloc(%d, %d): could not alloc m_hidden_output_to", ninputs, nhidden));
  ae->v_output = VectorAlloc(noutputs, MATRIX_REAL);
  if (ae->v_output == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "AEaalloc(%d, %d): could not alloc v_output", noutputs, nhidden));
  ae->v_hidden_bias = VectorAlloc(nhidden, MATRIX_REAL);
  if (ae->v_hidden_bias == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "AEaalloc(%d, %d): could not alloc v_hidden_bias", ninputs, nhidden));
  ae->hidden_class_labels = (int *)calloc(nhidden, sizeof(int));
  if (ae->hidden_class_labels == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "AEaalloc(%d, %d): could not alloc hidden class labels", ninputs, nhidden));
  ae->v_hidden_net = VectorAlloc(nhidden, MATRIX_REAL);
  if (ae->v_hidden_net == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "AEaalloc(%d, %d): could not alloc v_hidden_net", ninputs, nhidden));
  ae->v_hidden = VectorAlloc(nhidden, MATRIX_REAL);
  if (ae->v_hidden == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "AEaalloc(%d, %d): could not alloc v_hidden", ninputs, nhidden));
  ae->v_output_bias = VectorAlloc(noutputs, MATRIX_REAL);
  if (ae->v_output_bias == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "AEaalloc(%d, %d): could not alloc v_output_bias", ninputs, nhidden));
  ae->prev = prev;
  if (ae->prev) ae->prev->next = ae;
  ae->next = NULL;

  for (norm = 0.0, i = 1; i <= ninputs; i++)
    for (j = 1; j <= nhidden; j++) {
      w = randomNumber(-1, 1);
      w = randomNumber(-wt_lim, wt_lim);
      norm += (w * w);
      *MATRIX_RELT(ae->m_input_to_hidden, j, i) = w;
    }

  for (norm = 0.0, j = 1; j <= nhidden; j++)
    for (k = 1; k <= noutputs; k++) {
      w = randomNumber(-1, 1);
      w = randomNumber(-wt_lim, wt_lim);
      norm += (w * w);
      *MATRIX_RELT(ae->m_hidden_to_output, k, j) = w;
    }

  for (norm = 0.0, j = 1; j <= nhidden; j++) {
    w = randomNumber(-1, 1);
    norm += (w * w);
    VECTOR_ELT(ae->v_hidden_bias, j) = w;
  }
  norm = sqrt(norm);
  for (j = 1; j <= nhidden; j++) VECTOR_ELT(ae->v_hidden_bias, j) /= norm;

  for (norm = 0.0, k = 1; k <= noutputs; k++) {
    w = randomNumber(-1, 1);
    norm += (w * w);
    VECTOR_ELT(ae->v_output_bias, k) = w;
  }
  norm = sqrt(norm);
  for (k = 1; k <= noutputs; k++) VECTOR_ELT(ae->v_output_bias, k) /= norm;

  if (ninputs == 1)  // diags - set weights to correct values
  {
    //    double f0, f1;

    //  f0 = SIGMOID(0);
    // f1 = SIGMOID(1);
    //    ae->v_hidden_bias->rptr[1][1] = 0 ;
    //    ae->m_input_to_hidden->rptr[1][1] = 1 ;
    //    ae->m_hidden_to_output->rptr[1][1] = 1 / (f1-f0) ;
    //    ae->v_output_bias->rptr[1][1] = 0 - (1 / (f1-f0)) * f0 ;
  }

  ae->average_act = (double *)calloc(ae->v_hidden_bias->rows, sizeof(double));
  if (ae->average_act == NULL) ErrorExit(ERROR_NOMEMORY, "AEalloc: could not allocate average_act");

  return (ae);
}

static void AEfree(AE **pae)
{
  AE *ae = *pae;

  *pae = NULL;
  MatrixFree(&ae->m_input_to_hidden);
  MatrixFree(&ae->m_hidden_to_output);
  VectorFree(&ae->v_output_bias);
  VectorFree(&ae->v_hidden_bias);
  VectorFree(&ae->v_hidden);
  VectorFree(&ae->v_hidden_net);
  if (ae->v_error) VectorFree(&ae->v_error);
  free(ae->hidden_class_labels);
  if (ae->prev == NULL) {
    VectorFree(&ae->v_input);
    VectorFree(&ae->v_output);
  }
  if (ae->zero_indices) free(ae->zero_indices);
  if (ae->saved_inputs) free(ae->saved_inputs);
  if (ae->m_grad_input_to_hidden)  // free training stuff
  {
    VectorFree(&ae->v_grad_hidden_bias);
    VectorFree(&ae->v_grad_output_bias);
    VectorFree(&ae->v_previous_step_hidden_bias);
    VectorFree(&ae->v_previous_step_output_bias);

    MatrixFree(&ae->m_grad_input_to_hidden);
    MatrixFree(&ae->m_grad_hidden_to_output);
    MatrixFree(&ae->m_previous_step_input_to_hidden);
    MatrixFree(&ae->m_previous_step_hidden_to_output);
  }
  free(ae->average_act);
}

void SAEfree(SAE **psae)
{
  SAE *sae = *psae;
  AE *ae, *next;

  *psae = NULL;

  ae = sae->first;
  do {
    next = ae->next;
    AEfree(&ae);
    ae = next;
  } while (ae);

  free(sae);
}
AE *SAEaddLayer(SAE *sae, float scale)
{
  AE *ae, *last;
  int nhidden, noutputs;

  last = SAEfindLastLayer(sae, sae->first);
  noutputs = last->v_hidden_bias->rows;
  nhidden = nint(noutputs * scale);
  //  noutputs = sae->type == FOCUSED_AUTOENCODER ? 1 : last->v_hidden_bias->rows ;
  ae = AEalloc(last, noutputs, nhidden, noutputs);
  sae->nencoders++;
  printf("stacked layer #%d added with %d hidden units\n", sae->nencoders, nhidden);
  return (ae);
}

SAE *SAEtrainLayer(SAE *sae, AE *layer, MRI **mri, double tol)
{
  int x, y, z;

  for (x = 0; x < mri[0]->width; x++)
    for (y = 0; y < mri[0]->height; y++)
      for (z = 0; z < mri[0]->depth; z++) {
        SAEfillInputVector(mri, sae->nlevels, x, y, z, sae->whalf, sae->first->v_input);
        SAEactivateNetwork(sae);
      }
  return (sae);
}
double SAEcomputeTotalRMS(SAE *sae, MRI **mri)
{
  int x, y, z, nvox;
  double rms, total_rms;

  total_rms = 0.0;
  nvox = 0;
  for (x = 0; x < mri[0]->width; x++)
    for (y = 0; y < mri[0]->height; y++)
      for (z = 0; z < mri[0]->depth; z++) {
        if (FZERO(MRIgetVoxVal(mri[0], x, y, z, 0))) continue;
        SAEfillInputVector(mri, sae->nlevels, x, y, z, sae->whalf, sae->first->v_input);
        SAEactivateNetwork(sae);
        rms = SAEcomputeRMS(sae);
        if (!devFinite(rms)) DiagBreak();
        if (rms > 1000) DiagBreak();
        total_rms += rms;
        nvox++;
      }
  return (total_rms / nvox);
}

VECTOR *SAEactivateNetwork(SAE *sae)
{
  AEactivateLayer(sae->first, sae->first->v_input);
  return (sae->first->v_output);
}

static int AEactivateLayer(AE *ae, VECTOR *v_input)
{
  int row, ninputs = ae->v_input->rows;
  double o, net;

  MatrixCopy(v_input, ae->v_input);
  if (ae->noise_fraction > 0) {
    int i, max_i;
    max_i = nint(ninputs * ae->noise_fraction);
    ae->zero_indices = compute_permutation(ninputs, ae->zero_indices);
    if (ae->saved_inputs == NULL) ae->saved_inputs = (double *)calloc(max_i, sizeof(double));

    for (i = 0; i < max_i; i++) {
      ae->saved_inputs[i] = VECTOR_ELT(ae->v_input, ae->zero_indices[i] + 1);
      VECTOR_ELT(ae->v_input, ae->zero_indices[i] + 1) = 0;
    }
  }
  MatrixMultiply(ae->m_input_to_hidden, ae->v_input, ae->v_hidden_net);
  if (ae->noise_fraction > 0) {
    int i, max_i;
    max_i = nint(ninputs * ae->noise_fraction);
    for (i = 0; i < max_i; i++)
      VECTOR_ELT(ae->v_input, ae->zero_indices[i] + 1) = ae->saved_inputs[i];  // restore inputs for training
  }
  MatrixAdd(ae->v_hidden_net, ae->v_hidden_bias, ae->v_hidden_net);
  for (row = 1; row <= ae->v_hidden->rows; row++)  // apply hidden layer nonlinearity
  {
    net = VECTOR_ELT(ae->v_hidden_net, row);
    if (!devFinite(net)) DiagBreak();
    o = SIGMOID(net);
    if (!devFinite(o)) DiagBreak();
    VECTOR_ELT(ae->v_hidden, row) = o;
  }
  if (ae->next) {
    AEactivateLayer(ae->next, ae->v_hidden);  // recurse - this hidden unit is next layers input
    MatrixCopy(ae->next->v_output,
               ae->v_hidden);  // output of next layer is hidden activation of this one since it autoencodes
  }

  MatrixMultiply(ae->m_hidden_to_output, ae->v_hidden, ae->v_output);
  MatrixAdd(ae->v_output, ae->v_output_bias, ae->v_output);
  return (NO_ERROR);
}

VECTOR *SAEactivateLastHiddenLayer(SAE *sae, MRI *mri) { return (sae->first->v_output); }

double SAEcomputeRMS(SAE *sae)
{
  //  VECTOR *v_output;
  int row;
  double error, rms;

  // v_output =
  SAEactivateNetwork(sae);

  if (sae->type == NORMAL_AUTOENCODER) {
    for (rms = 0.0, row = 1; row <= sae->first->v_input->rows; row++) {
      error = VECTOR_ELT(sae->first->v_input, row) - VECTOR_ELT(sae->first->v_output, row);
      rms += error * error;
    }
  }
  else {
    row = (sae->first->v_input->rows + 1) / 2;
    error = VECTOR_ELT(sae->first->v_input, row) - VECTOR_ELT(sae->first->v_output, 1);
    rms = error * error;
    row = 1;  // for rms calculation
  }
  return (sqrt(rms / row));
}
static double AEcomputeRMS(AE *ae)
{
  int row;
  double error, rms;

  if (ae->v_input->rows == ae->v_output->rows)  // NORMAL_AUTOENCODER
  {
    for (rms = 0.0, row = 1; row <= ae->v_input->rows; row++) {
      error = VECTOR_ELT(ae->v_input, row) - VECTOR_ELT(ae->v_output, row);
      rms += error * error;
    }
  }
  else {
    row = (ae->v_input->rows + 1) / 2;
    error = VECTOR_ELT(ae->v_input, row) - VECTOR_ELT(ae->v_output, 1);
    rms = error * error;
    row = 1;  // for rms calculation
  }

  return (sqrt(rms / row));
}

// NOTE: input MRI must be type float and scaled to be 0->1 (that is, scale down MRI_UCHAR vols by 255)
double SAEtrainFromMRI(SAE *sae, MRI **mri_pyramid, SAE_INTEGRATION_PARMS *parms)
{
  double error = 0.0, rms, last_rms, total_rms, last_total_rms, pct_decrease, running_last_rms, G_rms,
         G_last_rms = 1e10;
  int x, y, z, iter = 0, visited, ind, nvox, calls = 0;
  short *x_ind, *y_ind, *z_ind;
  double tol;
  // double dt,
  double acceptance_sigma, proposal_sigma;
  char *out_fname;
  MRI *mri = mri_pyramid[0];
  AE *ae_train;  // the deepest layer - which is what we are training now

  getVolGeom(mri_pyramid[0], &sae->vg);
  acceptance_sigma = parms->acceptance_sigma;
  proposal_sigma = parms->proposal_sigma;
  tol = parms->tol;
  // dt = parms->dt;
  out_fname = parms->out_fname;

  if (mri_pyramid[0]->type != MRI_FLOAT)
    ErrorExit(ERROR_BADPARM, "SAEtrainFromMRI: input type must be MRI_FLOAT scaled to [0->1]");
  nvox = mri->width * mri->height * mri->depth;

  x_ind = (short *)calloc(nvox, sizeof(short));
  y_ind = (short *)calloc(nvox, sizeof(short));
  z_ind = (short *)calloc(nvox, sizeof(short));
  if (!x_ind || !y_ind || !z_ind) ErrorExit(ERROR_NOMEMORY, "SAEtrainFromMRI: could not allocate permutation indices");

  ae_train = SAEfindLastLayer(sae, NULL);
  last_total_rms = 1;
  //  last_total_rms = SAEcomputeTotalRMS(sae, mri_pyramid) ;
  printf("%3.3d: rms = %2.4f\n", iter, last_total_rms);
  if (Gx >= 0) {
    int wsize, ind;
    wsize = sae->whalf * 2 + 1;
    ind = (wsize * wsize * wsize) / 2 + 1;
    SAEfillInputVector(mri_pyramid, sae->nlevels, Gx, Gy, Gz, sae->whalf, sae->first->v_input);
    SAEactivateNetwork(sae);
    G_last_rms = SAEcomputeRMS(sae);
    printf("voxel (%d, %d, %d), I = %2.2f, o = %2.2f, rms = %2.4f\n",
           Gx,
           Gy,
           Gz,
           sae->first->v_input->rptr[ind][1],
           sae->first->v_output->rptr[ind][1],
           G_last_rms);
  }
  do {
    MRIcomputeVoxelPermutation(mri, x_ind, y_ind, z_ind);
    for (running_last_rms = total_rms = 0.0, visited = ind = 0; ind < nvox; ind++) {
      if (ind && !(ind % MAX((nvox / 500), 2))) {
        double running_pct_dec;
        running_pct_dec = 100 * (running_last_rms / visited - total_rms / visited) /
                          (running_last_rms / visited + total_rms / visited);
        pct_decrease = 100 * (last_total_rms - total_rms / visited) / (last_total_rms + total_rms / visited);
        printf("%d of %d (%2.2f%%): rms = %2.4f (%2.3f%%) (%2.3f%%)\n",
               ind,
               nvox,
               100 * ind / (float)nvox,
               total_rms / visited,
               pct_decrease,
               running_pct_dec);
      }
      x = x_ind[ind];
      y = y_ind[ind];
      z = z_ind[ind];
      if (x < sae->whalf || y < sae->whalf || z < sae->whalf || x >= mri->width - sae->whalf ||
          y >= mri->height - sae->whalf || z >= mri->depth - sae->whalf)
        continue;
      if (!FZERO(MRIgetVoxVal(mri, x, y, z, 0)))  // there is real data at this spot
      {
        visited++;
        if (++calls == Gdiag_no) DiagBreak();
        SAEfillInputVector(mri_pyramid, sae->nlevels, x, y, z, sae->whalf, sae->first->v_input);
        SAEactivateNetwork(sae);
        last_rms = SAEcomputeRMS(sae);
        running_last_rms += last_rms;
        AEtrain(ae_train, parms);
        SAEactivateNetwork(sae);
        rms = SAEcomputeRMS(sae);
        total_rms += rms;
        if (rms > last_rms) DiagBreak();
        if (Gx >= 0) {
          // int ind;
          // float in, out;
          // ind = (sae->first->v_output->rows + 1) / 2;

          SAEfillInputVector(mri_pyramid, sae->nlevels, Gx, Gy, Gz, sae->whalf, sae->first->v_input);
          SAEactivateNetwork(sae);
          G_rms = SAEcomputeRMS(sae);
          if (G_rms > G_last_rms) DiagBreak();
          G_last_rms = G_rms;
          // in = sae->first->v_input->rptr[ind][1];
          // out = sae->first->v_output->rptr[ind][1];
          DiagBreak();
        }
      }
      if (out_fname && !(ind % MAX((nvox / 50), 2))) {
        char fname[STRLEN], path[STRLEN];

        if (sae->first->v_input->rows != 1 || !(ind % MAX((nvox / 500), 2))) {
          FileNameRemoveExtension(out_fname, path);
          int cx = snprintf(fname, STRLEN, "%s.%3.3d.%2.2d.ae", path, ind, iter);
	  if( (cx<0) || (cx>STRLEN) ) {
	    std::cerr << __FUNCTION__
		      << ": snprintf returned error on line "
		      << __LINE__ << std::endl;
	  }
          printf("writing SAE after %dth iteration to %s\n", iter, fname);
          SAEwrite(sae, fname);
        }
      }
    }
    total_rms /= visited;
    last_total_rms = running_last_rms / visited;
    pct_decrease = 100 * (last_total_rms - total_rms) / (last_total_rms + total_rms);
    last_total_rms = total_rms;
    printf("%3.3d: rms = %2.4f (%2.3f%%)\n", ++iter, total_rms, pct_decrease);
    if (out_fname) {
      char fname[STRLEN], path[STRLEN];
      FileNameRemoveExtension(out_fname, path);
      int cx = snprintf(fname, STRLEN, "%s.%2.2d.ae", path, iter);
      if( (cx<0) || (cx>STRLEN) ) {
	std::cerr << __FUNCTION__
		  << ": snprintf returned error on line "
		  << __LINE__ << std::endl;
      }
      printf("writing SAE after %dth iteration to %s\n", iter, fname);
      SAEwrite(sae, fname);
    }
    if (parms->integration_type == INTEGRATE_BOLTZMANN_MACHINE) {
      if (iter < 30) pct_decrease = 2 * tol;
      parms->acceptance_sigma *= .8;
      parms->proposal_sigma *= .8;
      printf("setting acceptance/proposal sigma to %2.3f/%2.3f\n", parms->acceptance_sigma, parms->proposal_sigma);
    }
  } while (pct_decrease > tol);

  if (Gx >= 0) {
    // int wsize, ind,
    int i, j;
    // float in, out, total_rms, init_total_rms;
    // wsize = sae->whalf * 2 + 1;
    // ind = (wsize * wsize * wsize) / 2 + 1;

    // init_total_rms =
    SAEcomputeTotalRMS(sae, mri_pyramid);
    for (j = 0; j < 10; j++) {
      // total_rms =
      SAEcomputeTotalRMS(sae, mri_pyramid);
      SAEfillInputVector(mri_pyramid, sae->nlevels, Gx, Gy, Gz, sae->whalf, sae->first->v_input);
      SAEactivateNetwork(sae);
      last_rms = SAEcomputeRMS(sae);

      for (i = 0; i < 100; i++) {
        AEtrain(sae->first, parms);
        SAEactivateNetwork(sae);
        rms = SAEcomputeRMS(sae);
      }
    }

    G_rms = SAEcomputeRMS(sae);
    if (G_rms > G_last_rms) DiagBreak();
    G_last_rms = G_rms;
    // in = sae->first->v_input->rptr[ind][1];
    // out = sae->first->v_output->rptr[ind][1];
    DiagBreak();
  }

  parms->acceptance_sigma = acceptance_sigma;
  parms->proposal_sigma = proposal_sigma;
  free(x_ind);
  free(y_ind);
  free(z_ind);
  return (error);
}

// NOTE: input MRI must be type float and scaled to be 0->1 (that is, scale down MRI_UCHAR vols by 255)
double SAEtrainFromVoxlist(SAE *sae, VOXEL_LIST *vl, MRI **mri_pyramid, SAE_INTEGRATION_PARMS *parms)
{
  double error = 0.0, rms, last_rms, total_rms, last_total_rms, pct_decrease, running_last_rms, G_rms,
         G_last_rms = 1e10;
  int x, y, z, iter = 0, visited, ind, calls = 0, *indices, i;
  double tol;
  // double dt,
  double acceptance_sigma, proposal_sigma;
  char *out_fname;
  MRI *mri = mri_pyramid[0];
  AE *ae_train;  // the deepest layer - which is what we are training now

  getVolGeom(mri_pyramid[0], &sae->vg);
  acceptance_sigma = parms->acceptance_sigma;
  proposal_sigma = parms->proposal_sigma;
  tol = parms->tol;
  // dt = parms->dt;
  out_fname = parms->out_fname;

  if (mri_pyramid[0]->type != MRI_FLOAT)
    ErrorExit(ERROR_BADPARM, "SAEtrainFromMRI: input type must be MRI_FLOAT scaled to [0->1]");

  ae_train = SAEfindLastLayer(sae, NULL);
  last_total_rms = 1;
  last_total_rms = SAEcomputeTotalRMS(sae, mri_pyramid);
  printf("%3.3d: rms = %2.4f\n", iter, last_total_rms);
  if (Gx >= 0) {
    int wsize, ind;
    wsize = sae->whalf * 2 + 1;
    ind = (wsize * wsize * wsize) / 2 + 1;
    SAEfillInputVector(mri_pyramid, sae->nlevels, Gx, Gy, Gz, sae->whalf, sae->first->v_input);
    SAEactivateNetwork(sae);
    G_last_rms = SAEcomputeRMS(sae);
    printf("voxel (%d, %d, %d), I = %2.2f, o = %2.2f, rms = %2.4f\n",
           Gx,
           Gy,
           Gz,
           sae->first->v_input->rptr[ind][1],
           sae->first->v_output->rptr[ind][1],
           G_last_rms);
  }
  do {
    indices = compute_permutation(vl->nvox, NULL);
    for (running_last_rms = total_rms = 0.0, visited = ind = 0; ind < vl->nvox; ind++) {
      if (ind && !(ind % MAX((vl->nvox / 100), 2))) {
        double running_pct_dec;
        running_pct_dec = 100 * (running_last_rms / visited - total_rms / visited) /
                          (running_last_rms / visited + total_rms / visited);
        pct_decrease = 100 * (last_total_rms - total_rms / visited) / (last_total_rms + total_rms / visited);
        printf("%d of %d (%2.2f%%): rms = %2.4f (%2.3f%%) (%2.3f%%)\n",
               ind,
               vl->nvox,
               100 * ind / (float)vl->nvox,
               total_rms / visited,
               pct_decrease,
               running_pct_dec);
      }
      i = indices[ind];
      x = vl->xi[i];
      y = vl->yi[i];
      z = vl->zi[i];
      if (x < sae->whalf || y < sae->whalf || (z < sae->whalf && mri->depth > 1) || x >= mri->width - sae->whalf ||
          y >= mri->height - sae->whalf || (z >= mri->depth - sae->whalf && mri->depth > 1))
        continue;

      visited++;
      if (++calls == Gdiag_no) DiagBreak();
      SAEfillInputVector(mri_pyramid, sae->nlevels, x, y, z, sae->whalf, sae->first->v_input);
      SAEactivateNetwork(sae);
      last_rms = SAEcomputeRMS(sae);
      running_last_rms += last_rms;
      AEtrain(ae_train, parms);
      SAEactivateNetwork(sae);
      rms = SAEcomputeRMS(sae);
      total_rms += rms;
      if (rms > last_rms) DiagBreak();
      if (Gx >= 0) {
        // int ind;
        // float in, out;
        // ind = (sae->first->v_output->rows + 1) / 2;

        SAEfillInputVector(mri_pyramid, sae->nlevels, Gx, Gy, Gz, sae->whalf, sae->first->v_input);
        SAEactivateNetwork(sae);
        G_rms = SAEcomputeRMS(sae);
        if (G_rms > G_last_rms) DiagBreak();
        G_last_rms = G_rms;
        // in = sae->first->v_input->rptr[ind][1];
        // out = sae->first->v_output->rptr[ind][1];
        DiagBreak();
      }

      if (out_fname && !(ind % MAX((vl->nvox / 50), 2))) {
        char fname[STRLEN], path[STRLEN];

        if (sae->first->v_input->rows != 1 || !(ind % MAX((vl->nvox / 500), 2))) {
          FileNameRemoveExtension(out_fname, path);
          int cx = snprintf(fname, STRLEN, "%s.%2.2d.%4.4d.ae", path, iter, ind);
	  if( (cx<0) || (cx>STRLEN) ) {
	    std::cerr << __FUNCTION__
		      << ": snprintf returned error on line "
		      << __LINE__ << std::endl;
	  }
          printf("writing SAE after %dth iteration to %s\n", iter, fname);
          SAEwrite(sae, fname);
        }
      }
    }
    total_rms /= visited;
    last_total_rms = running_last_rms / visited;
    pct_decrease = 100 * (last_total_rms - total_rms) / (last_total_rms + total_rms);
    last_total_rms = total_rms;
    printf("%3.3d: rms = %2.4f (%2.3f%%)\n", ++iter, total_rms, pct_decrease);
    if (out_fname) {
      char fname[STRLEN], path[STRLEN];
      FileNameRemoveExtension(out_fname, path);
      int cx = snprintf(fname, STRLEN, "%s.%2.2d.ae", path, iter);
      if( (cx<0) || (cx>STRLEN) ) {
	std::cerr << __FUNCTION__
		  << ": snprintf returned error on line "
		  << __LINE__ << std::endl;
      }
      printf("writing SAE after %dth iteration to %s\n", iter, fname);
      SAEwrite(sae, fname);
    }
    if (parms->integration_type == INTEGRATE_BOLTZMANN_MACHINE) {
      if (iter < 30) pct_decrease = 2 * tol;
      parms->acceptance_sigma *= .8;
      parms->proposal_sigma *= .8;
      printf("setting acceptance/proposal sigma to %2.3f/%2.3f\n", parms->acceptance_sigma, parms->proposal_sigma);
    }
  } while (pct_decrease > tol);

  if (Gx >= 0) {
    // int wsize, ind;
    int i, j;

    // float in, out, total_rms, init_total_rms;
    // wsize = sae->whalf * 2 + 1;
    // ind = (wsize * wsize * wsize) / 2 + 1;

    // init_total_rms =
    SAEcomputeTotalRMS(sae, mri_pyramid);
    for (j = 0; j < 10; j++) {
      // total_rms =
      SAEcomputeTotalRMS(sae, mri_pyramid);
      SAEfillInputVector(mri_pyramid, sae->nlevels, Gx, Gy, Gz, sae->whalf, sae->first->v_input);
      SAEactivateNetwork(sae);
      last_rms = SAEcomputeRMS(sae);

      for (i = 0; i < 100; i++) {
        AEtrain(sae->first, parms);
        SAEactivateNetwork(sae);
        rms = SAEcomputeRMS(sae);
      }
    }

    G_rms = SAEcomputeRMS(sae);
    if (G_rms > G_last_rms) DiagBreak();
    G_last_rms = G_rms;
    // in = sae->first->v_input->rptr[ind][1];
    // out = sae->first->v_output->rptr[ind][1];
    DiagBreak();
  }

  parms->acceptance_sigma = acceptance_sigma;
  parms->proposal_sigma = proposal_sigma;
  return (error);
}

static int CSAEcomputeGradient(CSAE *csae,
                               AE *ae,
                               VOXEL_LIST *vl,
                               MRI *mri,
                               int start_index,
                               int end_index,
                               int *indices,
                               SAE_INTEGRATION_PARMS *parms)
{
  int ind, i, j, k, whalf, x, y, z, ninputs, nhidden, noutputs, num_indices;

  AEclearGradients(ae);
  whalf = (ae->ksize - 1) / 2;
  for (i = start_index; i <= end_index; i++) {
    ind = indices[i];
    x = vl->xi[ind];
    y = vl->yi[ind];
    z = vl->zi[ind];
    parms->class_label = vl->vsrc[ind];
    if (x < whalf || y < whalf || (z < whalf && mri->depth > 1) || x >= mri->width - whalf ||
        y >= mri->height - whalf || (z >= mri->depth - whalf && mri->depth > 1))
      continue;

    if (x == Gx && y == Gy && z == Gz) DiagBreak();

    CSAEfillInputs(csae, mri, ae->v_input, x, y, z, ae->ksize);
    AEactivateLayer(ae, ae->v_input);
    AEaccumulateGradient(ae, parms);
  }
  num_indices = end_index - start_index + 1;
  MatrixScalarMul(ae->m_grad_input_to_hidden, 1.0 / num_indices, ae->m_grad_input_to_hidden);
  MatrixScalarMul(ae->m_grad_hidden_to_output, 1.0 / num_indices, ae->m_grad_hidden_to_output);
  MatrixScalarMul(ae->v_grad_hidden_bias, 1.0 / num_indices, ae->v_grad_hidden_bias);
  MatrixScalarMul(ae->v_grad_output_bias, 1.0 / num_indices, ae->v_grad_output_bias);

  ninputs = ae->v_input->rows;
  noutputs = ae->v_output->rows;
  nhidden = ae->v_hidden->rows;

  // compute sparsity and apply it to hidden bias and weights
  for (j = 1; j <= nhidden; j++) {
    ae->average_act[j - 1] /= (double)num_indices;  // make it the average activation level
    VECTOR_ELT(ae->v_grad_hidden_bias, j) += parms->sparsity_trate * (ae->average_act[j - 1] - parms->sparsity_target);
    for (i = 1; i <= ninputs; i++)
      *MATRIX_RELT(ae->m_grad_input_to_hidden, j, i) +=
          (parms->sparsity_trate / ninputs) * (ae->average_act[j - 1] - parms->sparsity_target);
  }

  // compute weight decay
  for (i = 1; i <= ninputs; i++)
    for (j = 1; j <= nhidden; j++)
      *MATRIX_RELT(ae->m_grad_input_to_hidden, j, i) +=
          parms->weight_decays[parms->layer] * *MATRIX_RELT(ae->m_input_to_hidden, j, i);
  for (j = 1; j <= nhidden; j++)
    for (k = 1; k <= noutputs; k++)
      *MATRIX_RELT(ae->m_grad_hidden_to_output, k, j) +=
          parms->weight_decays[parms->layer] * *MATRIX_RELT(ae->m_hidden_to_output, k, j);

  return (NO_ERROR);
}

static int AEclearGradients(AE *ae)
{
  int j;

  if (ae->m_grad_input_to_hidden == NULL) return (NO_ERROR);
  for (j = 0; j < ae->v_hidden->rows; j++) ae->average_act[j] = 0;
  MatrixClear(ae->m_grad_input_to_hidden);
  MatrixClear(ae->m_grad_hidden_to_output);
  MatrixClear(ae->v_grad_hidden_bias);
  MatrixClear(ae->v_grad_output_bias);
  return (NO_ERROR);
}

/*

 Simplify[D[1 / (1+Exp[-(w*x+b)]),w]]

           b + w x
          E        x
Out[7]= ---------------
              b + w x 2
        (1 + E       )
 Simplify[D[1 / (1+Exp[-(w*x+b)]),b]]

            b + w x
           E
Out[8]= ---------------
              b + w x 2
        (1 + E       )


*/

static double AEaccumulateGradient(AE *ae, SAE_INTEGRATION_PARMS *parms)
{
  double rms = 0, error, wt;
  int i, j, k, ninputs, nhidden, noutputs;

  wt = (1 - parms->class_weight);
  ninputs = ae->v_input->rows;
  noutputs = ae->v_output->rows;
  nhidden = ae->v_hidden->rows;
  if (ae->v_error == NULL) {
    ae->v_error = VectorClone(ae->v_output);
    ae->m_grad_input_to_hidden = MatrixClone(ae->m_input_to_hidden);
    ae->m_grad_hidden_to_output = MatrixClone(ae->m_hidden_to_output);
    ae->v_grad_hidden_bias = VectorClone(ae->v_hidden);
    ae->v_grad_output_bias = VectorClone(ae->v_output);

    ae->m_previous_step_input_to_hidden = MatrixClone(ae->m_input_to_hidden);
    ae->m_previous_step_hidden_to_output = MatrixClone(ae->m_hidden_to_output);
    ae->v_previous_step_hidden_bias = VectorClone(ae->v_hidden);
    ae->v_previous_step_output_bias = VectorClone(ae->v_output);
  }
  if (noutputs == ninputs)  // normal autoencoder
  {
    for (k = 1; k <= ae->v_input->rows; k++) {
      error = VECTOR_ELT(ae->v_error, k) = VECTOR_ELT(ae->v_output, k) - VECTOR_ELT(ae->v_input, k);
      rms += error * error;
    }
  }
  else {
    int ind;
    ind = (ninputs + 1) / 2;
    rms = VECTOR_ELT(ae->v_error, 1) = VECTOR_ELT(ae->v_output, 1) - VECTOR_ELT(ae->v_input, ind);
    rms *= rms;
  }

// NOTE: compute all GRADIENTS first. Will multiply by -dt at the end to get right direction of update

//  compute output weights grad
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) firstprivate(noutputs, nhidden) shared(ae) schedule(static, 1)
#endif
  for (k = 1; k <= noutputs; k++) {
    ROMP_PFLB_begin
    double error, hidden;
    int j;

    error = VECTOR_ELT(ae->v_error, k);
    VECTOR_ELT(ae->v_grad_output_bias, k) += error;

    for (j = 1; j <= nhidden; j++) {
      hidden = VECTOR_ELT(ae->v_hidden, j);
      *MATRIX_RELT(ae->m_grad_hidden_to_output, k, j) += error * hidden;
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

// compute hidden bias grad
  ROMP_PF_begin
#ifdef HAVE_OPENMP
#pragma omp parallel for if_ROMP(experimental) firstprivate(noutputs, nhidden) shared(ae) schedule(static, 1)
#endif
  for (j = 1; j <= nhidden; j++) {
    ROMP_PFLB_begin
    
    double dE_dbj, fprime, error, wjk, hidden;
    int k;

    hidden = VECTOR_ELT(ae->v_hidden, j);
    fprime = D_SIGMOID(hidden);
    for (dE_dbj = 0.0, k = 1; k <= noutputs; k++) {
      error = VECTOR_ELT(ae->v_error, k);
      wjk = *MATRIX_RELT(ae->m_hidden_to_output, k, j);
      dE_dbj += error * wjk * fprime;
    }
    VECTOR_ELT(ae->v_grad_hidden_bias, j) += wt * dE_dbj;
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
// compute hidden weight grad
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) firstprivate(noutputs, nhidden) shared(ae) schedule(static, 1)
#endif
  for (i = 1; i <= noutputs; i++) {
    ROMP_PFLB_begin
    
    int j;
    double Ii;

    Ii = VECTOR_ELT(ae->v_input, i);
    for (j = 1; j <= nhidden; j++) {
      double hidden, fprimej, dE_dwij, wjk;
      int k;

      hidden = VECTOR_ELT(ae->v_hidden, j);
      fprimej = D_SIGMOID(hidden);

      for (dE_dwij = 0.0, k = 1; k <= noutputs; k++) {
        error = VECTOR_ELT(ae->v_error, k);
        wjk = *MATRIX_RELT(ae->m_hidden_to_output, k, j);
        dE_dwij += error * wjk;
      }
      dE_dwij *= Ii * fprimej;
      *MATRIX_RELT(ae->m_grad_input_to_hidden, j, i) += wt * dE_dwij;
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  if (parms->class_label >= 0) {
    if (parms->class_label == Gdiag_no) DiagBreak();

    for (j = 1; j <= nhidden; j++) {
      double fprimej, dE_dwij, hidden, target, dE_dbj, Ii;
      int i;

      if (j == Gdiag_no) DiagBreak();

      if (ae->hidden_class_labels[j - 1] <= 0) continue;  // no class label specified for this node

      hidden = VECTOR_ELT(ae->v_hidden, j);
      fprimej = D_SIGMOID(hidden);

      target = (ae->hidden_class_labels[j - 1] == parms->class_label) ? 1 : 0;
      error = hidden - target;
      for (dE_dwij = 0.0, i = 1; i <= noutputs; i++) {
        Ii = VECTOR_ELT(ae->v_input, i);
        dE_dwij = error * fprimej * Ii;
        *MATRIX_RELT(ae->m_grad_input_to_hidden, j, i) += parms->class_weight * dE_dwij;
      }

      dE_dbj = error * fprimej;

      VECTOR_ELT(ae->v_grad_hidden_bias, j) += parms->class_weight * dE_dbj;  // XXX added division by noutputs
    }
  }

  for (j = 0; j < nhidden;
       j++)  // keep track of average hidden node activation for use in sparsity gradient calculated later
    ae->average_act[j] += VECTOR_ELT(ae->v_hidden, j + 1);

  return (rms);
}

static double aeApplyAccumulatedGradient(AE *ae, SAE_INTEGRATION_PARMS *parms)
{
  double rms = 0, dt, momentum, Egrad, Erandom;
  int noutputs;
  // int ninputs, nhidden;

  dt = parms->orig_dt;
  momentum = parms->momentum;
  // ninputs = ae->v_input->rows;
  noutputs = ae->v_output->rows;
  // nhidden = ae->v_hidden->rows;

  if (parms->integration_type != INTEGRATE_CONJUGATE_GRADIENT)  // don't scale grads by -dt for conjugate gradient
  {
    MatrixScalarMul(ae->v_grad_output_bias, -dt, ae->v_grad_output_bias);
    MatrixScalarMul(ae->v_grad_hidden_bias, -dt, ae->v_grad_hidden_bias);
    MatrixScalarMul(ae->m_grad_hidden_to_output, -dt, ae->m_grad_hidden_to_output);
    MatrixScalarMul(ae->m_grad_input_to_hidden, -dt, ae->m_grad_input_to_hidden);
  }
  if (!FZERO(momentum))  // use momentum
  {
    MatrixScalarMul(ae->v_previous_step_output_bias, momentum, ae->v_previous_step_output_bias);
    MatrixAdd(ae->v_grad_output_bias, ae->v_previous_step_output_bias, ae->v_grad_output_bias);

    MatrixScalarMul(ae->v_previous_step_hidden_bias, momentum, ae->v_previous_step_hidden_bias);
    MatrixAdd(ae->v_grad_hidden_bias, ae->v_previous_step_hidden_bias, ae->v_grad_hidden_bias);

    MatrixScalarMul(ae->m_previous_step_input_to_hidden, momentum, ae->m_previous_step_input_to_hidden);
    MatrixAdd(ae->m_grad_input_to_hidden, ae->m_previous_step_input_to_hidden, ae->m_grad_input_to_hidden);

    MatrixScalarMul(ae->m_previous_step_hidden_to_output, momentum, ae->m_previous_step_hidden_to_output);
    MatrixAdd(ae->m_grad_hidden_to_output, ae->m_previous_step_hidden_to_output, ae->m_grad_hidden_to_output);
  }

  AEsaveState(ae);
  if (parms->integration_type == INTEGRATE_GRADIENT_DESCENT) {
    MatrixAdd(ae->v_grad_output_bias, ae->v_output_bias, ae->v_output_bias);
    MatrixAdd(ae->v_grad_hidden_bias, ae->v_hidden_bias, ae->v_hidden_bias);
    MatrixAdd(ae->m_grad_hidden_to_output, ae->m_hidden_to_output, ae->m_hidden_to_output);
    MatrixAdd(ae->m_grad_input_to_hidden, ae->m_input_to_hidden, ae->m_input_to_hidden);

    MatrixCopy(ae->v_grad_output_bias, ae->v_previous_step_output_bias);
    MatrixCopy(ae->v_grad_hidden_bias, ae->v_previous_step_hidden_bias);
    MatrixCopy(ae->m_grad_hidden_to_output, ae->m_previous_step_hidden_to_output);
    MatrixCopy(ae->m_grad_input_to_hidden, ae->m_previous_step_input_to_hidden);
    rms = sqrt(rms / noutputs);
  }
  else if (parms->integration_type == INTEGRATE_CONJUGATE_GRADIENT) {
    double beta, best_dt, rms, best_rms;
    // double orig_dt;
    static int callno = 0;

    if (++callno == Gdiag_no) DiagBreak();

    if (parms->v_dir_hidden_bias == NULL)  // allocate everything the first time
    {
      parms->v_prev_grad_change_hidden_bias = MatrixClone(ae->v_hidden_bias);
      parms->v_prev_grad_change_output_bias = MatrixClone(ae->v_output_bias);
      parms->m_prev_grad_change_input_to_hidden = MatrixClone(ae->m_input_to_hidden);
      parms->m_prev_grad_change_hidden_to_output = MatrixClone(ae->m_hidden_to_output);
      parms->v_prev_grad_hidden_bias = MatrixClone(ae->v_hidden_bias);
      parms->v_prev_grad_output_bias = MatrixClone(ae->v_output_bias);
      parms->m_prev_grad_input_to_hidden = MatrixClone(ae->m_input_to_hidden);
      parms->m_prev_grad_hidden_to_output = MatrixClone(ae->m_hidden_to_output);
      parms->v_dir_hidden_bias = MatrixClone(ae->v_hidden_bias);
      parms->v_dir_output_bias = MatrixClone(ae->v_output_bias);
      parms->m_dir_input_to_hidden = MatrixClone(ae->m_input_to_hidden);
      parms->m_dir_hidden_to_output = MatrixClone(ae->m_hidden_to_output);
      parms->norm_hidden_bias = parms->norm_output_bias = parms->norm_hidden_to_output = parms->norm_input_to_hidden =
          1.0;
    }

    if (DZERO(parms->norm_hidden_bias) || !devFinite(parms->norm_hidden_bias)) DiagBreak();
    beta = VectorDot(parms->v_prev_grad_change_hidden_bias, parms->v_prev_grad_hidden_bias) / parms->norm_hidden_bias;
    if (!devFinite(beta)) DiagBreak();
    MatrixScalarMul(parms->v_dir_hidden_bias, beta, parms->v_dir_hidden_bias);
    MatrixSubtract(parms->v_dir_hidden_bias, ae->v_grad_hidden_bias, parms->v_dir_hidden_bias);
    //    MatrixAdd(parms->v_dir_hidden_bias, ae->v_hidden_bias, ae->v_hidden_bias) ;

    AEsaveState(ae);
    best_dt = 0;
    best_rms = AEcomputeRMS(ae);
    // for (orig_dt = dt, dt = dt * .1; dt <= 1000; dt *= 2) {
    for (dt = dt * .1; dt <= 1000; dt *= 2) {
      aeApplyGradient(ae, parms, dt);
      AEactivateLayer(ae, ae->v_input);
      rms = AEcomputeRMS(ae);
      if (rms < best_rms) {
        best_rms = rms;
        best_dt = dt;
      }
      AErestoreState(ae);
    }
    aeApplyGradient(ae, parms, best_dt);
    parms->dt = best_dt;
  }
  else if (parms->integration_type == INTEGRATE_BOLTZMANN_MACHINE) {
    static MATRIX *m_hidden_to_output_delta = NULL, *m_input_to_hidden_delta = NULL, *v_hidden_bias_delta = NULL,
                  *v_output_bias_delta = NULL;
    double acceptance_val;

    MatrixAdd(ae->v_grad_output_bias, ae->v_output_bias, ae->v_output_bias);
    MatrixAdd(ae->v_grad_hidden_bias, ae->v_hidden_bias, ae->v_hidden_bias);
    MatrixAdd(ae->m_grad_hidden_to_output, ae->m_hidden_to_output, ae->m_hidden_to_output);
    MatrixAdd(ae->m_grad_input_to_hidden, ae->m_input_to_hidden, ae->m_input_to_hidden);

    AEactivateLayer(ae, ae->v_input);
    Egrad = AEcomputeRMS(ae);
    AErestoreState(ae);  // undo gradient steop
    m_hidden_to_output_delta =
        MatrixDRand48ZeroMean(ae->m_hidden_to_output->rows, ae->m_hidden_to_output->cols, m_hidden_to_output_delta);
    m_input_to_hidden_delta =
        MatrixDRand48ZeroMean(ae->m_input_to_hidden->rows, ae->m_input_to_hidden->cols, m_input_to_hidden_delta);
    v_hidden_bias_delta = MatrixDRand48ZeroMean(ae->v_hidden_bias->rows, ae->v_hidden_bias->cols, v_hidden_bias_delta);
    v_output_bias_delta = MatrixDRand48ZeroMean(ae->v_output_bias->rows, ae->v_output_bias->cols, v_output_bias_delta);

    MatrixScalarMul(m_hidden_to_output_delta, dt * parms->proposal_sigma, m_hidden_to_output_delta);
    MatrixScalarMul(m_input_to_hidden_delta, dt * parms->proposal_sigma, m_input_to_hidden_delta);
    MatrixScalarMul(v_hidden_bias_delta, dt * parms->proposal_sigma, v_hidden_bias_delta);
    MatrixScalarMul(v_output_bias_delta, dt * parms->proposal_sigma, v_output_bias_delta);

    MatrixAdd(v_output_bias_delta, ae->v_output_bias, ae->v_output_bias);
    MatrixAdd(v_hidden_bias_delta, ae->v_hidden_bias, ae->v_hidden_bias);
    MatrixAdd(m_hidden_to_output_delta, ae->m_hidden_to_output, ae->m_hidden_to_output);
    MatrixAdd(m_input_to_hidden_delta, ae->m_input_to_hidden, ae->m_input_to_hidden);
    AEactivateLayer(ae, ae->v_input);
    Erandom = AEcomputeRMS(ae);
    acceptance_val = exp((Egrad - Erandom) / parms->acceptance_sigma);
    if (randomNumber(0.0, 1.0) > acceptance_val)  // take gradient step
    {
      AErestoreState(ae);  // undo gradient steop
      MatrixAdd(ae->v_grad_output_bias, ae->v_output_bias, ae->v_output_bias);
      MatrixAdd(ae->v_grad_hidden_bias, ae->v_hidden_bias, ae->v_hidden_bias);
      MatrixAdd(ae->m_grad_hidden_to_output, ae->m_hidden_to_output, ae->m_hidden_to_output);
      MatrixAdd(ae->m_grad_input_to_hidden, ae->m_input_to_hidden, ae->m_input_to_hidden);
    }
    else
      DiagBreak();
  }

  return (rms);
}

static int aeApplyGradient(AE *ae, SAE_INTEGRATION_PARMS *parms, double dt)
{
  static MATRIX *m_grad_hidden_to_output, *m_grad_input_to_hidden;
  static VECTOR *v_grad_hidden_bias, *v_grad_output_bias;

  if (v_grad_output_bias == NULL) {
    v_grad_output_bias = VectorCopy(ae->v_grad_output_bias, NULL);
    v_grad_hidden_bias = VectorCopy(ae->v_grad_hidden_bias, NULL);
    m_grad_hidden_to_output = MatrixCopy(ae->m_grad_hidden_to_output, NULL);
    m_grad_input_to_hidden = MatrixCopy(ae->m_grad_input_to_hidden, NULL);
  }

  v_grad_output_bias = MatrixScalarMul(parms->v_dir_output_bias, dt, v_grad_output_bias);
  v_grad_hidden_bias = MatrixScalarMul(parms->v_dir_hidden_bias, dt, v_grad_hidden_bias);
  m_grad_hidden_to_output = MatrixScalarMul(parms->m_dir_hidden_to_output, dt, m_grad_hidden_to_output);
  m_grad_input_to_hidden = MatrixScalarMul(parms->m_dir_input_to_hidden, dt, m_grad_input_to_hidden);

  MatrixAdd(v_grad_output_bias, ae->v_output_bias, ae->v_output_bias);
  MatrixAdd(v_grad_hidden_bias, ae->v_hidden_bias, ae->v_hidden_bias);
  MatrixAdd(m_grad_hidden_to_output, ae->m_hidden_to_output, ae->m_hidden_to_output);
  MatrixAdd(m_grad_input_to_hidden, ae->m_input_to_hidden, ae->m_input_to_hidden);

  return (NO_ERROR);
}

static double AEtrain(AE *ae, SAE_INTEGRATION_PARMS *parms)
{
  double rms = 0, error, dt, momentum, Egrad, Erandom, wt;
  int i, j, k, ninputs, nhidden, noutputs;

  dt = parms->orig_dt;
  momentum = parms->momentum;
  wt = (1 - parms->class_weight);
  ninputs = ae->v_input->rows;
  noutputs = ae->v_output->rows;
  nhidden = ae->v_hidden->rows;
  if (ae->v_error == NULL) {
    ae->v_error = VectorClone(ae->v_output);
    ae->m_grad_input_to_hidden = MatrixClone(ae->m_input_to_hidden);
    ae->m_grad_hidden_to_output = MatrixClone(ae->m_hidden_to_output);
    ae->v_grad_hidden_bias = VectorClone(ae->v_hidden);
    ae->v_grad_output_bias = VectorClone(ae->v_output);

    ae->m_previous_step_input_to_hidden = MatrixClone(ae->m_input_to_hidden);
    ae->m_previous_step_hidden_to_output = MatrixClone(ae->m_hidden_to_output);
    ae->v_previous_step_hidden_bias = VectorClone(ae->v_hidden);
    ae->v_previous_step_output_bias = VectorClone(ae->v_output);
  }
  if (noutputs == ninputs)  // normal autoencoder
  {
    for (k = 1; k <= ae->v_input->rows; k++) {
      error = VECTOR_ELT(ae->v_error, k) = VECTOR_ELT(ae->v_output, k) - VECTOR_ELT(ae->v_input, k);
      rms += error * error;
    }
  }
  else {
    int ind;
    ind = (ninputs + 1) / 2;
    rms = VECTOR_ELT(ae->v_error, 1) = VECTOR_ELT(ae->v_output, 1) - VECTOR_ELT(ae->v_input, ind);
    rms *= rms;
  }

// NOTE: compute all GRADIENTS first. Will multiply by -dt at the end to get right direction of update

// update compute output weights grad
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) firstprivate(noutputs, nhidden, dt) shared(ae) schedule(static, 1)
#endif
  for (k = 1; k <= noutputs; k++) {
    ROMP_PFLB_begin
    
    double error, hidden;
    int j;

    error = VECTOR_ELT(ae->v_error, k);
    VECTOR_ELT(ae->v_grad_output_bias, k) += error;
    if (!devFinite(VECTOR_ELT(ae->v_grad_output_bias, k))) DiagBreak();
    for (j = 1; j <= nhidden; j++) {
      hidden = VECTOR_ELT(ae->v_hidden, j);
      *MATRIX_RELT(ae->m_grad_hidden_to_output, k, j) += error * hidden;  // *fprime fprime?
      if (!devFinite(*MATRIX_RELT(ae->m_grad_hidden_to_output, k, j))) DiagBreak();
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

// compute hidden bias grad
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) firstprivate(noutputs, nhidden, dt) shared(ae) schedule(static, 1)
#endif
  for (j = 1; j <= nhidden; j++) {
    ROMP_PFLB_begin
    
    double dE_dbj, fprime, error, wjk, hidden;
    int k;

    hidden = VECTOR_ELT(ae->v_hidden, j);
    fprime = D_SIGMOID(hidden);
    if (!devFinite(fprime) || !devFinite(hidden)) DiagBreak();
    for (dE_dbj = 0.0, k = 1; k <= noutputs; k++) {
      error = VECTOR_ELT(ae->v_error, k);
      wjk = *MATRIX_RELT(ae->m_hidden_to_output, k, j);
      dE_dbj += error * wjk * fprime;
      if (!devFinite(dE_dbj)) DiagBreak();
    }
    VECTOR_ELT(ae->v_grad_hidden_bias, j) += wt * dE_dbj;  // XXX added division by noutputs
    if (!devFinite(VECTOR_ELT(ae->v_grad_hidden_bias, j))) DiagBreak();
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

// compute hidden weight grad
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) firstprivate(noutputs, nhidden, dt) shared(ae) schedule(static, 1)
#endif
  for (i = 1; i <= noutputs; i++) {
    ROMP_PFLB_begin
    
    int j;
    double Ii;

    Ii = VECTOR_ELT(ae->v_input, i);
    for (j = 1; j <= nhidden; j++) {
      double netj, fprimej, dE_dwij, wjk, o;
      int k;

      netj = VECTOR_ELT(ae->v_hidden_net, j);
      o = SIGMOID(netj);
      fprimej = D_SIGMOID(o);
      if (!devFinite(fprimej) || !devFinite(o)) DiagBreak();

      for (dE_dwij = 0.0, k = 1; k <= noutputs; k++) {
        error = VECTOR_ELT(ae->v_error, k);
        wjk = *MATRIX_RELT(ae->m_hidden_to_output, k, j);
        dE_dwij += error * wjk;
      }
      dE_dwij *= Ii * fprimej;
      if (dE_dwij < 0) DiagBreak();
      *MATRIX_RELT(ae->m_grad_input_to_hidden, j, i) += wt * dE_dwij;  // XXX added division by noutputs
      if (!devFinite(*MATRIX_RELT(ae->m_grad_input_to_hidden, j, i))) DiagBreak();
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

  if (parms->class_label >= 0) {
    if (parms->class_label == 3)
      DiagBreak();
    else if (parms->class_label == 0)
      DiagBreak();

    for (j = 1; j <= nhidden; j++) {
      double fprimej, dE_dwij, hidden, target, dE_dbj, Ii;
      int i;

      if (j == Gdiag_no) DiagBreak();

      if (ae->hidden_class_labels[j - 1] <= 0) continue;  // no class label specified for this node

      hidden = VECTOR_ELT(ae->v_hidden, j);
      fprimej = D_SIGMOID(hidden);
      if (!devFinite(fprimej) || !devFinite(hidden)) DiagBreak();

      target = (ae->hidden_class_labels[j - 1] == parms->class_label) ? 1 : 0;
      error = hidden - target;
      for (dE_dwij = 0.0, i = 1; i <= noutputs; i++) {
        Ii = VECTOR_ELT(ae->v_input, i);
        dE_dwij = error * fprimej * Ii;
        *MATRIX_RELT(ae->m_grad_input_to_hidden, j, i) += parms->class_weight * dE_dwij;
        if (!devFinite(*MATRIX_RELT(ae->m_grad_input_to_hidden, j, i))) DiagBreak();
      }

      dE_dbj = error * fprimej;

      VECTOR_ELT(ae->v_grad_hidden_bias, j) += parms->class_weight * dE_dbj;  // XXX added division by noutputs
      if (!devFinite(VECTOR_ELT(ae->v_grad_hidden_bias, j))) DiagBreak();
    }
  }

  if (parms->integration_type != INTEGRATE_CONJUGATE_GRADIENT)  // don't scale grads by -dt for conjugate gradient
  {
    // compute sparsity
    for (j = 0; j < nhidden; j++) {
      ae->average_act[j] = .999 * ae->average_act[j] + (1 - .999) * VECTOR_ELT(ae->v_hidden, j + 1);
      VECTOR_ELT(ae->v_grad_hidden_bias, j + 1) +=
          parms->sparsity_trate * (ae->average_act[j] - parms->sparsity_target);
    }

    // compute weight decay
    for (i = 1; i <= ninputs; i++)
      for (j = 1; j <= nhidden; j++)
        *MATRIX_RELT(ae->m_grad_input_to_hidden, j, i) +=
            parms->weight_decays[parms->layer] * *MATRIX_RELT(ae->m_input_to_hidden, j, i);
    for (j = 1; j <= nhidden; j++)
      for (k = 1; k <= noutputs; k++)
        *MATRIX_RELT(ae->m_grad_hidden_to_output, k, j) +=
            parms->weight_decays[parms->layer] * *MATRIX_RELT(ae->m_hidden_to_output, k, j);

    MatrixScalarMul(ae->v_grad_output_bias, -dt, ae->v_grad_output_bias);
    MatrixScalarMul(ae->v_grad_hidden_bias, -dt, ae->v_grad_hidden_bias);
    MatrixScalarMul(ae->m_grad_hidden_to_output, -dt, ae->m_grad_hidden_to_output);
    MatrixScalarMul(ae->m_grad_input_to_hidden, -dt, ae->m_grad_input_to_hidden);
  }
  if (!FZERO(momentum))  // use momentum
  {
    MatrixScalarMul(ae->v_previous_step_output_bias, momentum, ae->v_previous_step_output_bias);
    MatrixAdd(ae->v_grad_output_bias, ae->v_previous_step_output_bias, ae->v_grad_output_bias);

    MatrixScalarMul(ae->v_previous_step_hidden_bias, momentum, ae->v_previous_step_hidden_bias);
    MatrixAdd(ae->v_grad_hidden_bias, ae->v_previous_step_hidden_bias, ae->v_grad_hidden_bias);

    MatrixScalarMul(ae->m_previous_step_input_to_hidden, momentum, ae->m_previous_step_input_to_hidden);
    MatrixAdd(ae->m_grad_input_to_hidden, ae->m_previous_step_input_to_hidden, ae->m_grad_input_to_hidden);

    MatrixScalarMul(ae->m_previous_step_hidden_to_output, momentum, ae->m_previous_step_hidden_to_output);
    MatrixAdd(ae->m_grad_hidden_to_output, ae->m_previous_step_hidden_to_output, ae->m_grad_hidden_to_output);
  }

  AEsaveState(ae);
  if (parms->integration_type == INTEGRATE_GRADIENT_DESCENT) {
    MatrixAdd(ae->v_grad_output_bias, ae->v_output_bias, ae->v_output_bias);
    MatrixAdd(ae->v_grad_hidden_bias, ae->v_hidden_bias, ae->v_hidden_bias);
    MatrixAdd(ae->m_grad_hidden_to_output, ae->m_hidden_to_output, ae->m_hidden_to_output);
    MatrixAdd(ae->m_grad_input_to_hidden, ae->m_input_to_hidden, ae->m_input_to_hidden);

    MatrixCopy(ae->v_grad_output_bias, ae->v_previous_step_output_bias);
    MatrixCopy(ae->v_grad_hidden_bias, ae->v_previous_step_hidden_bias);
    MatrixCopy(ae->m_grad_hidden_to_output, ae->m_previous_step_hidden_to_output);
    MatrixCopy(ae->m_grad_input_to_hidden, ae->m_previous_step_input_to_hidden);
    rms = sqrt(rms / noutputs);
  }
  else if (parms->integration_type == INTEGRATE_CONJUGATE_GRADIENT) {
    double beta, best_dt, rms, best_rms;
    // double orig_dt;
    static int callno = 0;

    if (++callno == Gdiag_no) DiagBreak();

    if (parms->v_dir_hidden_bias == NULL)  // allocate everything the first time
    {
      parms->v_prev_grad_change_hidden_bias = MatrixClone(ae->v_hidden_bias);
      parms->v_prev_grad_change_output_bias = MatrixClone(ae->v_output_bias);
      parms->m_prev_grad_change_input_to_hidden = MatrixClone(ae->m_input_to_hidden);
      parms->m_prev_grad_change_hidden_to_output = MatrixClone(ae->m_hidden_to_output);
      parms->v_prev_grad_hidden_bias = MatrixClone(ae->v_hidden_bias);
      parms->v_prev_grad_output_bias = MatrixClone(ae->v_output_bias);
      parms->m_prev_grad_input_to_hidden = MatrixClone(ae->m_input_to_hidden);
      parms->m_prev_grad_hidden_to_output = MatrixClone(ae->m_hidden_to_output);
      parms->v_dir_hidden_bias = MatrixClone(ae->v_hidden_bias);
      parms->v_dir_output_bias = MatrixClone(ae->v_output_bias);
      parms->m_dir_input_to_hidden = MatrixClone(ae->m_input_to_hidden);
      parms->m_dir_hidden_to_output = MatrixClone(ae->m_hidden_to_output);
      parms->norm_hidden_bias = parms->norm_output_bias = parms->norm_hidden_to_output = parms->norm_input_to_hidden =
          1.0;
    }

    if (DZERO(parms->norm_hidden_bias) || !devFinite(parms->norm_hidden_bias)) DiagBreak();
    beta = VectorDot(parms->v_prev_grad_change_hidden_bias, parms->v_prev_grad_hidden_bias) / parms->norm_hidden_bias;
    if (!devFinite(beta)) DiagBreak();
    MatrixScalarMul(parms->v_dir_hidden_bias, beta, parms->v_dir_hidden_bias);
    MatrixSubtract(parms->v_dir_hidden_bias, ae->v_grad_hidden_bias, parms->v_dir_hidden_bias);
    //    MatrixAdd(parms->v_dir_hidden_bias, ae->v_hidden_bias, ae->v_hidden_bias) ;

    AEsaveState(ae);
    best_dt = 0;
    best_rms = AEcomputeRMS(ae);
    // for (orig_dt = dt, dt = dt * .1; dt <= 1000; dt *= 2) {
    for (dt = dt * .1; dt <= 1000; dt *= 2) {
      aeApplyGradient(ae, parms, dt);
      AEactivateLayer(ae, ae->v_input);
      rms = AEcomputeRMS(ae);
      if (rms < best_rms) {
        best_rms = rms;
        best_dt = dt;
      }
      AErestoreState(ae);
    }
    aeApplyGradient(ae, parms, best_dt);
    parms->dt = best_dt;
  }
  else if (parms->integration_type == INTEGRATE_BOLTZMANN_MACHINE) {
    static MATRIX *m_hidden_to_output_delta = NULL, *m_input_to_hidden_delta = NULL, *v_hidden_bias_delta = NULL,
                  *v_output_bias_delta = NULL;
    double acceptance_val;

    MatrixAdd(ae->v_grad_output_bias, ae->v_output_bias, ae->v_output_bias);
    MatrixAdd(ae->v_grad_hidden_bias, ae->v_hidden_bias, ae->v_hidden_bias);
    MatrixAdd(ae->m_grad_hidden_to_output, ae->m_hidden_to_output, ae->m_hidden_to_output);
    MatrixAdd(ae->m_grad_input_to_hidden, ae->m_input_to_hidden, ae->m_input_to_hidden);

    AEactivateLayer(ae, ae->v_input);
    Egrad = AEcomputeRMS(ae);
    AErestoreState(ae);  // undo gradient steop
    m_hidden_to_output_delta =
        MatrixDRand48ZeroMean(ae->m_hidden_to_output->rows, ae->m_hidden_to_output->cols, m_hidden_to_output_delta);
    m_input_to_hidden_delta =
        MatrixDRand48ZeroMean(ae->m_input_to_hidden->rows, ae->m_input_to_hidden->cols, m_input_to_hidden_delta);
    v_hidden_bias_delta = MatrixDRand48ZeroMean(ae->v_hidden_bias->rows, ae->v_hidden_bias->cols, v_hidden_bias_delta);
    v_output_bias_delta = MatrixDRand48ZeroMean(ae->v_output_bias->rows, ae->v_output_bias->cols, v_output_bias_delta);

    MatrixScalarMul(m_hidden_to_output_delta, dt * parms->proposal_sigma, m_hidden_to_output_delta);
    MatrixScalarMul(m_input_to_hidden_delta, dt * parms->proposal_sigma, m_input_to_hidden_delta);
    MatrixScalarMul(v_hidden_bias_delta, dt * parms->proposal_sigma, v_hidden_bias_delta);
    MatrixScalarMul(v_output_bias_delta, dt * parms->proposal_sigma, v_output_bias_delta);

    MatrixAdd(v_output_bias_delta, ae->v_output_bias, ae->v_output_bias);
    MatrixAdd(v_hidden_bias_delta, ae->v_hidden_bias, ae->v_hidden_bias);
    MatrixAdd(m_hidden_to_output_delta, ae->m_hidden_to_output, ae->m_hidden_to_output);
    MatrixAdd(m_input_to_hidden_delta, ae->m_input_to_hidden, ae->m_input_to_hidden);
    AEactivateLayer(ae, ae->v_input);
    Erandom = AEcomputeRMS(ae);
    acceptance_val = exp((Egrad - Erandom) / parms->acceptance_sigma);
    if (randomNumber(0.0, 1.0) > acceptance_val)  // take gradient step
    {
      AErestoreState(ae);  // undo gradient steop
      MatrixAdd(ae->v_grad_output_bias, ae->v_output_bias, ae->v_output_bias);
      MatrixAdd(ae->v_grad_hidden_bias, ae->v_hidden_bias, ae->v_hidden_bias);
      MatrixAdd(ae->m_grad_hidden_to_output, ae->m_hidden_to_output, ae->m_hidden_to_output);
      MatrixAdd(ae->m_grad_input_to_hidden, ae->m_input_to_hidden, ae->m_input_to_hidden);
    }
    else
      DiagBreak();
  }

  return (rms);
}

VECTOR *SAEfillInputVector(MRI **mri, int nlevels, int x0, int y0, int z0, int whalf, VECTOR *v_input)
{
  int i, n, xk, yk, zk;
  double x, y, z, xi, yi, zi, scale, val, zmin, zmax;

  if (v_input == NULL) {
    int wsize = 2 * whalf + 1;

    v_input = VectorAlloc(wsize * wsize * wsize * nlevels, MATRIX_REAL);
  }

  if (mri[0]->depth > 1) {
    zmin = -whalf;
    zmax = whalf;
  }
  else
    zmin = zmax = 0;
  for (n = 0; n < nlevels; n++) {
    scale = pow(2.0, n);
    x = (x0 / scale);
    y = (y0 / scale);
    z = (z0 / scale);
    for (xk = -whalf, i = 1; xk <= whalf; xk++) {
      xi = (x + xk);
      if (xi < 0)
        xi = 0;
      else if (xi > mri[n]->width - 1)
        xi = mri[n]->width - 1;
      for (yk = -whalf; yk <= whalf; yk++) {
        yi = (y + yk);
        if (yi < 0)
          yi = 0;
        else if (yi > mri[n]->height - 1)
          yi = mri[n]->height - 1;
        for (zk = -zmin; zk <= zmax; zk++, i++) {
          zi = (z + zk);
          if (zi < 0)
            zi = 0;
          else if (zi > mri[n]->depth - 1)
            zi = mri[n]->depth - 1;
          //	  VECTOR_ELT(v_input, i) = MRIgetVoxVal(mri[n], xi, yi, zi, 0) ;
          MRIsampleVolume(mri[n], xi, yi, zi, &val);
          VECTOR_ELT(v_input, i) = val;
        }
      }
    }
  }

  return (v_input);
}

MRI *SAEvectorToMRI(VECTOR *v_input, int nlevels, int whalf, MRI *mri)
{
  int xk, yk, zk, i, wsize, n, zsize;

  wsize = (2 * whalf) + 1;

  if (mri == NULL) mri = MRIallocSequence(wsize, wsize, wsize, MRI_FLOAT, nlevels);

  if (mri->depth == 1)
    zsize = 1;
  else
    zsize = wsize;
  for (n = 0; n < nlevels; n++)
    for (xk = 0, i = 1; xk < wsize; xk++)
      for (yk = 0; yk < wsize; yk++)
        for (zk = 0; zk < zsize; zk++, i++) MRIsetVoxVal(mri, xk, yk, zk, n, VECTOR_ELT(v_input, i));

  return (mri);
}

MRI *SAEinputWeightsToMRI(SAE *sae, MRI *mri)
{
  int xk, yk, zk, i, wsize, frame, hidden, n;

  wsize = (2 * sae->whalf) + 1;

  if (mri == NULL) {
    if (sae->type & AUTOENCODER_2D)
      mri = MRIallocSequence(wsize, wsize, 1, MRI_FLOAT, sae->first->v_hidden->rows * sae->nlevels);
    else
      mri = MRIallocSequence(wsize, wsize, wsize, MRI_FLOAT, sae->first->v_hidden->rows * sae->nlevels);
  }

  //  useVolGeomToMRI(&sae->vg, mri);

  for (frame = hidden = 0; hidden < sae->first->v_hidden->rows; hidden++)
    for (i = 1, n = 0; n < sae->nlevels; n++, frame++)
      for (xk = 0; xk < wsize; xk++)
        for (yk = 0; yk < wsize; yk++)
          for (zk = 0; zk < mri->depth; zk++, i++)
            MRIsetVoxVal(mri, xk, yk, zk, frame, *MATRIX_RELT(sae->first->m_input_to_hidden, hidden + 1, i));

  return (mri);
}

int SAEwrite(SAE *sae, char *fname)
{
  FILE *fp;
  AE *ae;
  int i, n;

  fp = fopen(fname, "wb");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "SAEwrite(%s): could not open file", fname));

  n = fprintf(fp, "%d %d %lf %d %d\n", sae->whalf, sae->nencoders, sae->scale, sae->type, sae->nlevels);
  if (n != 5) {
    fprintf(stderr, "File(%s) Warning: expected number of values to read does not match those read\n", fname);
  }
  writeVolGeom(fp, &sae->vg);
  AEwrite(sae->first, fp);
  ae = sae->first;
  for (i = 1; i < sae->nencoders && ae->next; i++) {
    AEwrite(ae->next, fp);
    ae = ae->next;
  }
  fclose(fp);
  return (NO_ERROR);
}
static int AEwrite(AE *ae, FILE *fp)
{
  MatrixWriteInto(fp, ae->m_input_to_hidden);
  MatrixWriteInto(fp, ae->v_output_bias);
  MatrixWriteInto(fp, ae->v_hidden_bias);
  MatrixWriteInto(fp, ae->m_hidden_to_output);
  MatrixWriteInto(fp, ae->v_input);
  MatrixWriteInto(fp, ae->v_hidden);
  MatrixWriteInto(fp, ae->v_output);
  fwriteInt(ae->ksize, fp);
  return (NO_ERROR);
}
static AE *AEread(FILE *fp, AE *prev)
{
  AE *ae;
  MATRIX *m_input_to_hidden, *v_output_bias;

  if (feof(fp)) return (NULL);
  m_input_to_hidden = MatrixReadFrom(fp, NULL);
  v_output_bias = MatrixReadFrom(fp, NULL);
  ae = AEalloc(prev, m_input_to_hidden->cols, m_input_to_hidden->rows, v_output_bias->rows);
  MatrixCopy(v_output_bias, ae->v_output_bias);
  MatrixCopy(m_input_to_hidden, ae->m_input_to_hidden);
  MatrixFree(&m_input_to_hidden);
  VectorFree(&v_output_bias);

  MatrixReadFrom(fp, ae->v_hidden_bias);
  MatrixReadFrom(fp, ae->m_hidden_to_output);
  MatrixReadFrom(fp, ae->v_input);
  MatrixReadFrom(fp, ae->v_hidden);
  MatrixReadFrom(fp, ae->v_output);
  ae->ksize = freadInt(fp);

  return (ae);
}

SAE *SAEread(char *fname)
{
  FILE *fp;
  int whalf, nencoders, i, type, nlevels, n;
  SAE *sae;
  AE *ae;
  double scale;

  fp = fopen(fname, "rb");
  if (fp == NULL) ErrorReturn(NULL, (ERROR_NOFILE, "SAEread(%s): could not open file", fname));

  n = fscanf(fp, "%d %d %lf %d %d\n", &whalf, &nencoders, &scale, &type, &nlevels);
  if (n != 5) {
    fprintf(stderr, "File(%s) Warning: expected number of values to read does not match those read\n", fname);
  }

  sae = SAEalloc(whalf, nlevels, type, scale);
  readVolGeom(fp, &sae->vg);
  sae->nencoders = nencoders;
  ae = sae->first = AEread(fp, NULL);
  for (i = 1; i < nencoders; i++) {
    ae->next = AEread(fp, ae);
    ae->next->prev = ae;
    ae = ae->next;
    printf("layer %d read with %d hidden units\n", i + 1, ae->v_hidden->rows);
  }
  fclose(fp);
  return (sae);
}

void AEdump(AE *ae)
{
  if (ae->v_input->rows == 1) {
    printf("f(I=%2.2f * wij=%2.2f + bj=%2.2f = %2.2f) = %2.2f * wjk=%2.2f + bk=%2.2f = %2.2f, E=%2.2f\n",
           ae->v_input->rptr[1][1],
           ae->m_input_to_hidden->rptr[1][1],
           ae->v_hidden_bias->rptr[1][1],
           ae->v_hidden_net->rptr[1][1],
           ae->v_hidden->rptr[1][1],
           ae->m_hidden_to_output->rptr[1][1],
           ae->v_output_bias->rptr[1][1],
           ae->v_output->rptr[1][1],
           ae->v_output->rptr[1][1] - ae->v_input->rptr[1][1]);
    if (ae->m_grad_input_to_hidden)
      printf("grad: wij=%2.4f, bj=%2.4f, wjk=%2.4f, bk=%2.4f\n",
             ae->m_grad_input_to_hidden->rptr[1][1],
             ae->v_grad_hidden_bias->rptr[1][1],
             ae->m_grad_hidden_to_output->rptr[1][1],
             ae->v_grad_output_bias->rptr[1][1]);

    return;
  }
  printf("v_input\n");
  MatrixPrint(stdout, ae->v_input);
  printf("v_output\n");
  MatrixPrint(stdout, ae->v_output);
  printf("v_error\n");
  MatrixPrint(stdout, ae->v_error);

  printf("v_hidden_bias\n");
  MatrixPrint(stdout, ae->v_hidden_bias);
  printf("v_hidden_net\n");
  MatrixPrint(stdout, ae->v_hidden_net);
  printf("v_hidden\n");
  MatrixPrint(stdout, ae->v_hidden);
  printf("m_input_to_hidden\n");
  MatrixPrint(stdout, ae->m_input_to_hidden);

  printf("v_output_bias\n");
  MatrixPrint(stdout, ae->v_output_bias);
  printf("m_hidden_to_output\n");
  MatrixPrint(stdout, ae->m_hidden_to_output);

  if (ae->m_grad_input_to_hidden && 0) {
    printf("m_grad_input_to_hidden\n");
    MatrixPrint(Gstdout, ae->m_grad_input_to_hidden);
    printf("v_grad_hidden_bias\n");
    MatrixPrint(Gstdout, ae->v_grad_hidden_bias);
    printf("m_grad_hidden_to_output\n");
    MatrixPrint(Gstdout, ae->m_grad_hidden_to_output);
    printf("v_grad_output_bias\n");
    MatrixPrint(Gstdout, ae->v_grad_output_bias);
  }
}
static int AEsaveState(AE *ae)
{
  ae->m_saved_input_to_hidden = MatrixCopy(ae->m_input_to_hidden, ae->m_saved_input_to_hidden);
  ae->m_saved_hidden_to_output = MatrixCopy(ae->m_hidden_to_output, ae->m_saved_hidden_to_output);
  ae->v_saved_hidden_bias = MatrixCopy(ae->v_hidden_bias, ae->v_saved_hidden_bias);
  ae->v_saved_output_bias = MatrixCopy(ae->v_output_bias, ae->v_saved_output_bias);
  return (NO_ERROR);
}

static int AErestoreState(AE *ae)
{
  MatrixCopy(ae->m_saved_input_to_hidden, ae->m_input_to_hidden);
  MatrixCopy(ae->m_saved_hidden_to_output, ae->m_hidden_to_output);
  MatrixCopy(ae->v_saved_hidden_bias, ae->v_hidden_bias);
  MatrixCopy(ae->v_saved_output_bias, ae->v_output_bias);
  return (NO_ERROR);
}

// Convolutional SAE stuff

int CSAEwrite(CSAE *csae, char *fname) { return (SAEwrite(csae->sae, fname)); }

CSAE *CSAEalloc(int type, int nlayers, int *ksizes, int *ngroups, MRI *mri_inputs)
{
  int whalf = (ksizes[0] - 1) / 2, layer;
  double scale;
  AE *ae;
  CSAE *csae;

  csae = (CSAE *)calloc(1, sizeof(CSAE));
  scale = (float)ngroups[0] / (ksizes[0] * ksizes[0]);
  if (!(type & AUTOENCODER_2D)) scale /= (float)ksizes[0];
  csae->sae = SAEalloc(whalf, 1, type, scale);
  ae = csae->aes[layer = 0] = csae->sae->first;
  ae->ksize = ksizes[0];
  while (ae->next) {
    csae->aes[++layer] = ae->next;
    ae = ae->next;
  }
  for (layer = 0; layer < nlayers; layer++) {
    csae->mri_outputs[layer] =
        MRIallocSequence(mri_inputs->width, mri_inputs->height, mri_inputs->depth, MRI_FLOAT, ngroups[layer]);
    MRIcopyHeader(mri_inputs, csae->mri_outputs[layer]);
  }
  return (csae);
}
int CSAEfillInputs(CSAE *csae, MRI *mri_inputs, VECTOR *v_visible, int x0, int y0, int z0, int ksize)
{
  int xk, yk, xi, yi, whalf, v, f;
  float val;

  whalf = (ksize - 1) / 2;
  for (v = 1, f = 0; f < mri_inputs->nframes; f++)
    for (xk = -whalf; xk <= whalf; xk++) {
      xi = mri_inputs->xi[x0 + xk];
      for (yk = -whalf; yk <= whalf; yk++, v++) {
        yi = mri_inputs->yi[y0 + yk];
        val = MRIgetVoxVal(mri_inputs, xi, yi, 0, f);
        VECTOR_ELT(v_visible, v) = val;
      }
    }

  return (NO_ERROR);
}
MRI *CSAEcreateOutputs(CSAE *csae, MRI *mri_inputs, int first_layer, int last_layer)
{
  int layer, x, y, z, h;
  // int whalf;
  MRI *mri_layer_inputs, *mri_outputs = NULL;
  AE *ae, *next;
  static int callno = 0;

  callno++;
  for (layer = first_layer, mri_layer_inputs = mri_inputs; layer <= last_layer; layer++) {
    ae = csae->aes[layer];
    mri_outputs = csae->mri_outputs[layer];
    // whalf = (ae->ksize - 1) / 2;
    for (x = 0; x < mri_inputs->width; x++) {
      if (x && !(x % 100)) printf("layer %d of %d: x = %d of %d\n", layer, last_layer, x, mri_inputs->width);
      for (y = 0; y < mri_inputs->height; y++)
        for (z = 0; z < mri_inputs->depth; z++) {
          if (x == Gx && y == Gy) DiagBreak();
          CSAEfillInputs(csae, mri_layer_inputs, ae->v_input, x, y, z, ae->ksize);
          next = ae->next;
          ae->next = NULL;
          AEactivateLayer(ae, ae->v_input);
          ae->next = next;
          for (h = 0; h < ae->v_hidden->rows; h++)
            MRIsetVoxVal(mri_outputs, x, y, z, h, VECTOR_ELT(ae->v_hidden, h + 1));
        }
    }
    mri_layer_inputs = mri_outputs;
  }

  return (mri_outputs);
}
AE *CSAEaddLayer(CSAE *csae, int ksize, int nhidden)
{
  AE *ae, *last;
  int ninputs;
  // float scale;

  last = SAEfindLastLayer(csae->sae, csae->sae->first);
  ninputs = ksize * ksize * last->v_hidden_bias->rows;
  // scale = ninputs / nhidden;
  ae = AEalloc(last, ninputs, nhidden, ninputs);
  ae->ksize = ksize;
  printf("stacked layer #%d added with %d hidden units\n", csae->sae->nencoders, nhidden);
  csae->aes[csae->sae->nencoders++] = ae;

  return (ae);
}
// NOTE: input MRI must be type float and scaled to be 0->1 (that is, scale down MRI_UCHAR vols by 255)
double CSAEtrainLayerFromVoxlist(CSAE *csae, int layer, VOXEL_LIST *vl, MRI **mri_pyramid, SAE_INTEGRATION_PARMS *parms)
{
  double total_rms, last_total_rms, pct_decrease, running_last_rms = 0, running_rms;
  int iter = 0, ind, *indices, end_index, end_mini_batch, always, never, nbad, i;
  // double dt;
  double acceptance_sigma, proposal_sigma, tol, min_rms;
  char *out_fname;
  AE *ae_train;  // the deepest layer - which is what we are training now

  parms->layer = layer;
  getVolGeom(mri_pyramid[0], &csae->sae->vg);
  acceptance_sigma = parms->acceptance_sigma;
  proposal_sigma = parms->proposal_sigma;
  tol = parms->tol;
  // dt = parms->dt;
  out_fname = parms->out_fname;

  if (mri_pyramid[0]->type != MRI_FLOAT)
    ErrorExit(ERROR_BADPARM, "CSAEtrainLayerFromVoxlist: input type must be MRI_FLOAT scaled to [0->1]");

  ae_train = SAEfindLastLayer(csae->sae, NULL);
  //  last_total_rms = CSAEcomputeTotalRMS(csae, layer, mri_pyramid) ;

  last_total_rms = 0;
  end_index = vl->nvox - (parms->held_out + 1);
  min_rms = -1;
  nbad = 0;
  do {
    indices = compute_permutation(vl->nvox, NULL);
    if (iter == 0) {
      //      min_rms = last_total_rms = CSAEcomputeVoxlistRMS(csae, parms, layer, mri_pyramid, vl, indices,
      //      end_index+1, vl->nvox, &always, &never) ;
      min_rms = last_total_rms =
          CSAEcomputeVoxlistRMS(csae, parms, layer, mri_pyramid, vl, indices, 0, end_index, &always, &never);
      printf("%3.3d: rms = %2.4f\n", iter, last_total_rms);
    }
    if (nbad < parms->max_no_progress - 2) reset_constant_nodes(ae_train, 1000.0);
    for (i = ind = 0; ind <= end_index; ind += parms->mini_batch_size) {
      end_mini_batch = MIN(ind + parms->mini_batch_size - 1, vl->nvox - parms->held_out);
      if (i && (i % (nint(200.0 / parms->mini_batch_size)) == 0))
        running_last_rms =
            CSAEcomputeVoxlistRMS(csae, parms, layer, mri_pyramid, vl, indices, ind, end_mini_batch, &always, &never);

      CSAEcomputeGradient(csae, ae_train, vl, mri_pyramid[0], ind, end_mini_batch, indices, parms);
      aeApplyAccumulatedGradient(ae_train, parms);

      if (i && (i % (nint(200.0 / parms->mini_batch_size)) == 0)) {
        double running_pct_dec;
        running_rms =
            CSAEcomputeVoxlistRMS(csae, parms, layer, mri_pyramid, vl, indices, ind, end_mini_batch, &always, &never);
        total_rms = CSAEcomputeVoxlistRMS(csae, parms, layer, mri_pyramid, vl, indices, 0, end_index, &always, &never);
        running_pct_dec = 100 * (running_last_rms - running_rms) / (running_last_rms);
        pct_decrease = 100 * (last_total_rms - total_rms) / (last_total_rms);
        printf("%d of %d (%2.2f%%): rms = %2.4f (%2.3f%%) (%2.3f%%)\n",
               end_mini_batch + 1,
               end_index,
               100 * ind / (float)end_index,
               total_rms,
               pct_decrease,
               running_pct_dec);
      }

      if (out_fname && !(++i % nint(1000.0 / parms->mini_batch_size))) {
        char fname[STRLEN], path[STRLEN];

        if (csae->sae->first->v_input->rows != 1 || !(ind % MAX((vl->nvox / 10), 2))) {
          FileNameRemoveExtension(out_fname, path);
          int cx = snprintf(fname, STRLEN, "%s.layer%d.%2.2d.%4.4d.ae", path, layer, iter, ind);
	  if( (cx<0) || (cx>STRLEN) ) {
	    std::cerr << __FUNCTION__
		      << ": snprintf returned error on line "
		      << __LINE__ << std::endl;
	  }
          printf("writing CSAE after %dth iteration to %s\n", iter, fname);
          CSAEwrite(csae, fname);
        }
      }
    }
//    total_rms /= visited ; last_total_rms = running_last_rms / visited ;
    //    total_rms = CSAEcomputeVoxlistRMS(csae, parms, layer, mri_pyramid, vl, indices, end_index+1, vl->nvox,
    //    &always, &never) ;
    total_rms = CSAEcomputeVoxlistRMS(csae, parms, layer, mri_pyramid, vl, indices, 0, end_index, &always, &never);
    if (always || never) printf("\thidden nodes always (%d) or never (%d) active\n", always, never);
    pct_decrease = 100 * (last_total_rms - total_rms) / (last_total_rms + total_rms);
    last_total_rms = total_rms;
    printf("%3.3d: rms = %2.4f (%2.3f%%)\n", ++iter, total_rms, pct_decrease);
    if (out_fname) {
      char fname[STRLEN], path[STRLEN];
      FileNameRemoveExtension(out_fname, path);
      int cx = snprintf(fname, STRLEN, "%s.layer%d.%2.2d.ae", path, layer, iter);
      if( (cx<0) || (cx>STRLEN) ) {
	std::cerr << __FUNCTION__
		  << ": snprintf returned error on line "
		  << __LINE__ << std::endl;
      }
      printf("writing CSAE after %dth iteration to %s\n", iter, fname);
      CSAEwrite(csae, fname);
    }
    if (parms->integration_type == INTEGRATE_BOLTZMANN_MACHINE) {
      if (iter < 30) pct_decrease = 2 * tol;
      parms->acceptance_sigma *= .8;
      parms->proposal_sigma *= .8;
      printf("setting acceptance/proposal sigma to %2.3f/%2.3f\n", parms->acceptance_sigma, parms->proposal_sigma);
    }
    if (total_rms >= min_rms) {
      nbad++;
      printf("RMS %2.4f >=  min RMS %2.4f --> no progress = %d (max %d)\n",
             total_rms,
             min_rms,
             nbad,
             parms->max_no_progress);
    }
    else {
      printf("new min RMS %2.4f found (previous %2.4f)\n", total_rms, min_rms);
      min_rms = total_rms;
      nbad = 0;
    }
  } while ((nbad < parms->max_no_progress) && iter < parms->max_iter);


  parms->acceptance_sigma = acceptance_sigma;
  parms->proposal_sigma = proposal_sigma;
  return (total_rms);
}

MRI *CSAElayerWeightsToMRI(CSAE *csae, int layer)
{
  MRI *mri = NULL, *mri_prev, *mri_counts;
  int width, whalf_prev, x, y, xk, yk, v, h, count, xp, yp, hp;
  AE *ae, *ae_prev;
  float val, val_prev;

  if (layer >= csae->sae->nencoders)
    ErrorExit(ERROR_BADPARM, "CSAlayerWeightsToMRI: invalid layer %d (%d encoders)\n", layer, csae->sae->nencoders);

  if (layer <= 0)
    return (SAEinputWeightsToMRI(csae->sae, NULL));
  else if (layer >= 1) {
    ae = csae->aes[layer];
    ae_prev = csae->aes[layer - 1];
    whalf_prev = (ae_prev->ksize - 1) / 2;
    width = ae->ksize + 2 * whalf_prev;
    mri_prev = CSAElayerWeightsToMRI(csae, layer - 1);
    mri = MRIallocSequence(width, width, 1, MRI_FLOAT, ae->v_hidden->rows);
    MRIcopyHeader(csae->sae->mri_inputs, mri);
    mri_counts = MRIallocSequence(width, width, 1, MRI_INT, 1);
    MRIcopyHeader(mri, mri_counts);

    // v is the visible unit in this layer, which is the hidden unit (or frame) in the previous one
    for (h = 0; h < ae->v_hidden->rows; h++) {
      for (v = hp = 0; hp < ae_prev->v_hidden->rows; hp++) {
        for (x = 0; x < ae->ksize; x++) {
          for (y = 0; y < ae->ksize; y++, v++) {
            for (xp = 0; xp < ae_prev->ksize; xp++) {
              xk = mri->xi[x + xp];  // this is the location in the larger kernel
              for (yp = 0; yp < ae_prev->ksize; yp++) {
                yk = mri->yi[y + yp];  // y position in larger kernel
                if (xk == Gx && yk == Gy) DiagBreak();
                if (h == 0 && x == ae->ksize - 1 && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON && 0)
                  printf("x = %d, xp = %d, xi = %d    y = %d, yp = %d, yk = %d, v = %d\n", x, xp, xk, y, yp, yk, v);
                count = MRIgetVoxVal(mri_counts, xk, yk, 0, 0);
                count++;
                MRIsetVoxVal(mri_counts, xk, yk, 0, 0, count);
                val = MRIgetVoxVal(mri, xk, yk, 0, h);
                val_prev = MRIgetVoxVal(mri_prev, xp, yp, 0, hp);
                val_prev *= *MATRIX_RELT(ae->m_input_to_hidden, h + 1, v + 1);
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
MRI *SAElayerWeightsToMRI(SAE *sae, int layer)
{
  MRI *mri = NULL, *mri_prev, *mri_counts;
  int width, whalf_prev, x, y, xk, yk, v, h, count, xp, yp, hp;
  AE *ae, *ae_prev;
  float val, val_prev;

  if (layer >= sae->nencoders)
    ErrorExit(ERROR_BADPARM, "CSAlayerWeightsToMRI: invalid layer %d (%d encoders)\n", layer, sae->nencoders);

  if (layer <= 0)
    return (SAEinputWeightsToMRI(sae, NULL));
  else if (layer >= 1) {
    int l;
    for (ae_prev = ae = sae->first, l = 0; l < layer; l++) {
      ae_prev = ae;
      ae = ae->next;
    }
    mri_prev = SAElayerWeightsToMRI(sae, layer - 1);
    whalf_prev = (mri_prev->width - 1) / 2;
    width = ae->ksize + 2 * whalf_prev;

    mri = MRIallocSequence(width, width, 1, MRI_FLOAT, ae->v_hidden->rows);
    useVolGeomToMRI(&sae->vg, mri);
    mri_counts = MRIallocSequence(width, width, 1, MRI_INT, 1);
    MRIcopyHeader(mri, mri_counts);

    // v is the visible unit in this layer, which is the hidden unit (or frame) in the previous one
    for (h = 0; h < ae->v_hidden->rows; h++) {
      for (v = hp = 0; hp < ae_prev->v_hidden->rows; hp++) {
        for (x = 0; x < ae->ksize; x++) {
          for (y = 0; y < ae->ksize; y++, v++) {
            for (xp = 0; xp < mri_prev->width; xp++) {
              xk = mri->xi[x + xp];  // this is the location in the larger kernel
              for (yp = 0; yp < mri_prev->height; yp++) {
                yk = mri->yi[y + yp];  // y position in larger kernel
                if (xk == Gx && yk == Gy) DiagBreak();
                if (h == 0 && x == ae->ksize - 1 && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON && 0)
                  printf("x = %d, xp = %d, xi = %d    y = %d, yp = %d, yk = %d, v = %d\n", x, xp, xk, y, yp, yk, v);
                count = MRIgetVoxVal(mri_counts, xk, yk, 0, 0);
                count++;
                MRIsetVoxVal(mri_counts, xk, yk, 0, 0, count);
                val = MRIgetVoxVal(mri, xk, yk, 0, h);
                val_prev = MRIgetVoxVal(mri_prev, xp, yp, 0, hp);
                val_prev *= *MATRIX_RELT(ae->m_input_to_hidden, h + 1, v + 1);
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
double CSAEcomputeTotalRMS(CSAE *csae, int layer, MRI **mri)
{
  int x, y, z, nvox;
  double rms, total_rms;
  AE *ae;

  ae = csae->aes[layer];

  total_rms = 0.0;
  nvox = 0;
  for (x = 0; x < mri[0]->width; x++)
    for (y = 0; y < mri[0]->height; y++)
      for (z = 0; z < mri[0]->depth; z++) {
        if (FZERO(MRIgetVoxVal(mri[0], x, y, z, 0))) continue;
        CSAEfillInputs(csae, mri[0], ae->v_input, x, y, z, ae->ksize);
        AEactivateLayer(ae, ae->v_input);
        rms = AEcomputeRMS(ae);
        if (!devFinite(rms)) DiagBreak();
        if (rms > 1000) DiagBreak();
        total_rms += rms;
        nvox++;
      }
  return (total_rms / nvox);
}
static double AEcomputeHiddenRMS(AE *ae, SAE_INTEGRATION_PARMS *parms)
{
  int j;
  double target, rms, hidden, nhidden;

  nhidden = ae->v_hidden->rows;
  for (rms = 0.0, j = 0; j < nhidden; j++) {
    target = (ae->hidden_class_labels[j] == parms->class_label) ? 1 : 0;
    hidden = VECTOR_ELT(ae->v_hidden, j + 1);
    rms += SQR(target - hidden);
  }
  return (sqrt(rms / nhidden));
}

double CSAEcomputeVoxlistRMS(CSAE *csae,
                             SAE_INTEGRATION_PARMS *parms,
                             int layer,
                             MRI **mri,
                             VOXEL_LIST *vl,
                             int *indices,
                             int start_index,
                             int end_index,
                             int *always,
                             int *never)
{
  int x, y, z, nvox, i, ind, *histo, h, nhidden, iz, num_indices;
  double class_rms, rms, total_rms, total_class_rms;
  AE *ae;
  static double *last_rms = NULL;
  static double *last_class_rms = NULL;
  static int last_num_indices = 0;
  static int last_start_index = -1;

  num_indices = end_index - start_index + 1;
  if (last_num_indices != num_indices && last_rms) {
    free(last_rms);
    free(last_class_rms);
    last_rms = last_class_rms = NULL;
  }
  if (last_rms == NULL) {
    last_rms = (double *)calloc(num_indices, sizeof(double));
    last_class_rms = (double *)calloc(num_indices, sizeof(double));
  }
  if (last_start_index != start_index) {
    last_start_index = start_index;
    memset(last_rms, 0, num_indices * sizeof(double));
    memset(last_class_rms, 0, num_indices * sizeof(double));
  }
  last_num_indices = num_indices;

  ae = csae->aes[layer];
  nhidden = ae->v_hidden->rows;
  histo = (int *)calloc(nhidden, sizeof(int));

  total_class_rms = total_rms = 0.0;
  nvox = 0;

  for (i = start_index; i <= end_index; i++) {
    iz = i - start_index;
    ind = indices[i];
    x = vl->xi[ind];
    y = vl->yi[ind];
    z = vl->zi[ind];
    parms->class_label = vl->vsrc[ind];

    CSAEfillInputs(csae, mri[0], ae->v_input, x, y, z, ae->ksize);
    AEactivateLayer(ae, ae->v_input);
    for (h = 0; h < nhidden; h++)
      if ((VECTOR_ELT(ae->v_hidden, h + 1)) > .5) histo[h]++;
    rms = AEcomputeRMS(ae);
    total_rms += rms;
    class_rms = AEcomputeHiddenRMS(ae, parms);
    total_class_rms += class_rms;

    if (!devFinite(rms)) DiagBreak();
    if (rms > 1000) DiagBreak();
    nvox++;
    if (last_class_rms[iz] > 0 && class_rms > last_class_rms[iz]) DiagBreak();
    if (last_rms[iz] > 0 && rms > last_rms[iz]) DiagBreak();
    last_rms[iz] = rms;
    last_class_rms[iz] = class_rms;
  }

  total_rms = (1.0 - parms->class_weight) * total_rms + parms->class_weight * total_class_rms;
  if (always) {
    *always = *never = 0;
    for (h = 0; h < nhidden; h++) {
      if (histo[h] == 0)
        (*never)++;
      else if (histo[h] == nvox)
        (*always)++;
    }
  }

  free(histo);
  return (total_rms / nvox);
}
CSAE *CSAEread(char *fname)
{
  CSAE *csae;
  SAE *sae;
  int layer;
  AE *ae;

  sae = SAEread(fname);
  if (sae == NULL) return (NULL);
  csae = (CSAE *)calloc(1, sizeof(CSAE));
  csae->sae = sae;
  ae = csae->aes[layer = 0] = csae->sae->first;
  while (ae->next) {
    csae->aes[++layer] = ae->next;
    ae = ae->next;
  }
  for (layer = 0; layer < sae->nencoders; layer++) {
    ae = csae->aes[layer];
    csae->mri_outputs[layer] =
        MRIallocSequence(csae->sae->vg.width, csae->sae->vg.height, csae->sae->vg.depth, MRI_FLOAT, ae->v_hidden->rows);
    useVolGeomToMRI(&csae->sae->vg, csae->mri_outputs[layer]);
  }

  return (csae);
}

static int reset_constant_nodes(AE *ae, double thresh)
{
  static int calls = 0;
  int i, j, nreset = 0, ninputs;
  float mean, std, val;

  if (++calls == Gdiag_no) DiagBreak();
  ninputs = ae->v_input->rows;
  for (j = 1; j <= ae->v_hidden->rows; j++) {
    for (mean = std = 0.0, i = 1; i <= ninputs; i++) {
      val = *MATRIX_RELT(ae->m_input_to_hidden, j, i);
      mean += val;
      std += val * val;
    }
    mean /= ninputs;
    std = sqrt(std / ninputs - (mean * mean));
    if (FZERO(std) || (fabs(mean / std) > thresh)) {
      printf("*******************************  resetting hidden node  %d  *******************************\n", j);
      for (i = 1; i <= ninputs; i++) {
        nreset++;
        val = randomNumber(-.1, .1);
        *MATRIX_RELT(ae->m_input_to_hidden, j, i) = val;
      }
    }
  }
  return (nreset);
}
