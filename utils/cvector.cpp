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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "const.h"
#include "cvector.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "mrisurf.h"
#include "sig.h"

int cvector_compute_variance(float *var, float *mean, int norm, int num)
{
  int i;

  for (i = 0; i < num; i++) {
    var[i] = var[i] / (float)norm - mean[i] * mean[i];
    if (var[i] < 0) /* numerical instability */
      var[i] = 0;
  }
  return (NO_ERROR);
}
int cvector_normalize(float *v, float norm, int num)
{
  int i;

  for (i = 0; i < num; i++) v[i] /= norm;
  return (NO_ERROR);
}

int cvector_accumulate_square(float *v, float *vtotal, int num)
{
  int i;

  for (i = 0; i < num; i++) vtotal[i] += v[i] * v[i];
  return (NO_ERROR);
}
int cvector_accumulate(float *v, float *vtotal, int num)
{
  int i;

  for (i = 0; i < num; i++) vtotal[i] += v[i];
  return (NO_ERROR);
}

double cvector_compute_t_test(float *c1_mean,
                              float *c1_var,
                              float *c2_mean,
                              float *c2_var,
                              int num_class1,
                              int num_class2,
                              float *pvals,
                              int num,
                              int *pvno)
{
  int i;
  double t, numer, denom, p, max_p;

  max_p = 0;
  for (i = 0; i < num; i++) {
    if (i == Gdiag_no) DiagBreak();
    numer = (c1_mean[i] - c2_mean[i]);
    denom = sqrt((c1_var[i] / num_class1) + (c2_var[i] / num_class2));
    if (FZERO(denom)) {
      t = 0;
      p = 0.0f;
    }
    else {
      if (num_class1 == 1 || num_class2 == 1) {
        int dof_out = num_class1 + num_class2;
        /*
          if one of the means is not a population but a single subject
          then use standard deviation rather than standard error.
        */
        denom = sqrt((dof_out - 1) * ((c1_var[i] / num_class1) + (c2_var[i] / num_class2)));
      }
      else
        denom = sqrt((c1_var[i] / num_class1) + (c2_var[i] / num_class2));
      t = numer / denom;
      p = sigt(t, num_class1 + num_class2 - 2);
      p = log10(p);
    }
    if (t > 0) p *= -1;
    pvals[i] = p;

    if (fabs(p) > fabs(max_p)) {
      if (pvno) *pvno = i;
      max_p = p;
    }
  }
  return (max_p);
}
double cvector_compute_mean_diff(float *c1_mean, float *c2_mean, float *vmean_diff, int num, int *pvno)
{
  int i;
  double max_diff;

  max_diff = 0;
  for (i = 0; i < num; i++) {
    if (i == Gdiag_no) DiagBreak();
    vmean_diff[i] = (c1_mean[i] - c2_mean[i]);
    if (fabs(vmean_diff[i]) > fabs(max_diff)) {
      if (pvno) *pvno = i;
      max_diff = vmean_diff[i];
    }
  }
  return (max_diff);
}

int cvector_subtract(float *v1, float *v2, float *vdst, int num)
{
  int i;

  for (i = 0; i < num; i++) vdst[i] = v1[i] - v2[i];
  return (NO_ERROR);
}
int cvector_mark_low_prob_vertices(float *pvals, float pthresh, MRI_SURFACE *mris)
{
  int i, num;

  for (num = i = 0; i < mris->nvertices; i++) {
    if (pvals[i] < pthresh) {
      num++;
      mris->vertices[i].marked = 1;
    }
  }
  printf("%d vertices p < %2.2e\n", num, pthresh);
  return (NO_ERROR);
}

double cvector_compute_dist_free_snr(float **c1_thickness,
                                     int num_class1,
                                     float **c2_thickness,
                                     int num_class2,
                                     float *c1_mean,
                                     float *c2_mean,
                                     float *vsnr,
                                     int num,
                                     int *pi)
{
  int i, max_i, n, correct, total;
  double max_snr, mean, snr;

  max_i = -1;
  for (max_snr = 0.0, i = 0; i < num; i++) {
    mean = (c1_mean[i] + c2_mean[i]) / 2;
    snr = 0;
    correct = 0;
    if (c1_mean[i] > c2_mean[i]) {
      for (n = 0; n < num_class1; n++)
        if (c1_thickness[n][i] > mean) correct++;

      for (n = 0; n < num_class2; n++)
        if (c2_thickness[n][i] < mean) correct++;
    }
    else {
      for (n = 0; n < num_class1; n++)
        if (c1_thickness[n][i] < mean) correct++;

      for (n = 0; n < num_class2; n++)
        if (c2_thickness[n][i] > mean) snr++;
    }

    total = num_class1 + num_class2;
    snr = (double)correct / (double)total;
    vsnr[i] = snr;
    if (snr > max_snr) {
      max_i = i;
      max_snr = snr;
    }
  }
  *pi = max_i;
  return (max_snr);
}
double cvector_compute_snr(
    float *c1_mean, float *c2_mean, float *vvar, float *vsnr, int num, int *pi, float bonferroni, int stat_type)
{
  double snr;

  switch (stat_type) {
    default:
    case STAT_T:
    case STAT_F:
      snr = cvector_compute_snr_F(c1_mean, c2_mean, vvar, vsnr, num, pi, bonferroni);
      break;
  }
  return (snr);
}
double cvector_compute_snr_F(
    float *c1_mean, float *c2_mean, float *vvar, float *snr, int num, int *pi, float bonferroni)
{
  int i, max_i;
  float f, max_snr;

  max_i = -1;
  for (max_snr = i = 0; i < num; i++) {
    f = (c1_mean[i] - c2_mean[i]);
    f *= f;
    if (!iszero(vvar[i])) f /= (vvar[i]);

    f += bonferroni;
    if (c2_mean[i] > c1_mean[i]) f *= -1; /* make it a signed quantity */

    if (fabs(f) > max_snr) {
      max_snr = fabs(f);
      max_i = i;
    }
    snr[i] = f;
  }
  *pi = max_i;
  return (max_snr);
}
double cvector_len(float *v, int num)
{
  int i;
  double len;

  for (len = 0.0, i = 0; i < num; i++) len += v[i] * v[i];
  len /= (double)num;
  len = sqrt(len);
  return (len);
}

float *cvector_alloc(int num)
{
  float *v;

  v = (float *)calloc(num, sizeof(float));
  if (!v) ErrorExit(ERROR_NOMEMORY, "cvector_alloc(%d): calloc failed", num);
  return (v);
}

int cvector_clear(float *v, int num)
{
  int i;

  for (i = 0; i < num; i++) v[i] = 0;
  return (NO_ERROR);
}
int cvector_add_variances(
    float *c1_var, float *c2_var, int num_class1, int num_class2, float *vtotal_var, int nvertices)
{
  int i, total_dof;

  total_dof = num_class1 + num_class2;
  for (i = 0; i < nvertices; i++) vtotal_var[i] = (c1_var[i] * num_class1 + c2_var[i] * num_class2) / total_dof;

  return (NO_ERROR);
}

int cvector_multiply_variances(
    float *c1_var, float *c2_var, int num_class1, int num_class2, float *vtotal_var, int nvertices)
{
  int i, total_dof;

  total_dof = num_class1 + num_class2;
  for (i = 0; i < nvertices; i++) vtotal_var[i] = (c1_var[i] * num_class1 * c2_var[i] * num_class2) / total_dof;

  return (NO_ERROR);
}

int cvector_track_best_snr(float *vsnr,
                           float *vbest_snr,
                           float *vbest_avgs,
                           float *c1_mean,
                           float *c2_mean,
                           float *c1_best_mean,
                           float *c2_best_mean,
                           float *c1_var,
                           float *c2_var,
                           float *c1_best_var,
                           float *c2_best_var,
                           float **c1_avg_thickness,
                           float **c2_avg_thickness,
                           float **c1_best_thicknesses,
                           int nc1,
                           float **c2_best_thicknesses,
                           int nc2,
                           int avgs,
                           int num,
                           float fthresh,
                           int *pnum_found)
{
  int i, n;

  *pnum_found = 0;
  for (i = 0; i < num; i++) {
    if (fabs(vsnr[i]) > fabs(vbest_snr[i])) {
      vbest_snr[i] = vsnr[i];
      vbest_avgs[i] = avgs;
      c1_best_mean[i] = c1_mean[i];
      c1_best_var[i] = c1_var[i];
      c2_best_mean[i] = c2_mean[i];
      c2_best_var[i] = c2_var[i];
      if (vsnr[i] >= fthresh) *pnum_found += 1;
      for (n = 0; n < nc1; n++) c1_best_thicknesses[n][i] = c1_avg_thickness[n][i];
      for (n = 0; n < nc2; n++) c2_best_thicknesses[n][i] = c2_avg_thickness[n][i];
    }
  }

  return (NO_ERROR);
}
int cvector_track_best_stats(float *vpvals,
                             float *vbest_pvals,
                             float *vbest_avgs,
                             float *c1_mean,
                             float *c2_mean,
                             float *c1_best_mean,
                             float *c2_best_mean,
                             float *c1_var,
                             float *c2_var,
                             float *c1_best_var,
                             float *c2_best_var,
                             float **c1_avg_thickness,
                             float **c2_avg_thickness,
                             float **c1_best_thicknesses,
                             int nc1,
                             float **c2_best_thicknesses,
                             int nc2,
                             int avgs,
                             int num,
                             float fthresh,
                             int *pnum_found)
{
  int i, n;

  *pnum_found = 0;
  for (i = 0; i < num; i++) {
    if (fabs(vpvals[i]) > fabs(vbest_pvals[i])) {
      vbest_pvals[i] = vpvals[i];
      vbest_avgs[i] = avgs;
      c1_best_mean[i] = c1_mean[i];
      c1_best_var[i] = c1_var[i];
      c2_best_mean[i] = c2_mean[i];
      c2_best_var[i] = c2_var[i];
      if (fabs(vpvals[i]) >= fthresh) *pnum_found += 1;
      for (n = 0; n < nc1; n++) c1_best_thicknesses[n][i] = c1_avg_thickness[n][i];
      for (n = 0; n < nc2; n++) c2_best_thicknesses[n][i] = c2_avg_thickness[n][i];
    }
  }

  return (NO_ERROR);
}

int cvector_copy(float *v1, float *v2, int num)
{
  int i;

  for (i = 0; i < num; i++) v2[i] = v1[i];

  return (NO_ERROR);
}

int cvector_extract_best_avg(float *vbest_avgs, float *vsrc, float *vdst, int navgs, int num)
{
  int i;

  for (i = 0; i < num; i++) {
    if (nint(vbest_avgs[i]) == navgs) vdst[i] = vsrc[i];
  }
  return (NO_ERROR);
}
double cvector_average_in_label(float *v, LABEL *area, int num)
{
  int i;
  double avg;

  for (avg = 0.0, i = 0; i < area->n_points; i++) {
    if (!std::isfinite(v[area->lv[i].vno])) DiagBreak();
    avg += v[area->lv[i].vno];
  }
  avg /= (double)area->n_points;
  return (avg);
}

typedef struct
{
  float snr;
  int index;
} SORT_ELT;
int compare_sort_elts(const void *vse1, const void *vse2);
int compare_sort_elts(const void *vse1, const void *vse2)
{
  const SORT_ELT *se1, *se2;

  se1 = (const SORT_ELT *)vse1;
  se2 = (const SORT_ELT *)vse2;
  return (se2->snr > se1->snr);
}

int *cvector_sort(float *vbest_snr, int nvertices)
{
  int *sorted_vertices, i;
  SORT_ELT *se_table;

  se_table = (SORT_ELT *)calloc(nvertices, sizeof(SORT_ELT));
  if (!se_table) ErrorExit(ERROR_NOMEMORY, "cvector_sort: could not allocated %d vector", nvertices);
  sorted_vertices = (int *)calloc(nvertices, sizeof(int));
  if (!sorted_vertices) ErrorExit(ERROR_NOMEMORY, "cvector_sort: could not allocated %d vector", nvertices);

  for (i = 0; i < nvertices; i++) {
    se_table[i].index = i;
    se_table[i].snr = vbest_snr[i];
  }

  qsort(se_table, nvertices, sizeof(SORT_ELT), compare_sort_elts);
  for (i = 0; i < nvertices; i++) sorted_vertices[i] = se_table[i].index;
  free(se_table);
  return (sorted_vertices);
}
double cvector_compute_pvalues(float *c1_mean,
                               float *c1_var,
                               float *c2_mean,
                               float *c2_var,
                               int num_class1,
                               int num_class2,
                               float *pvals,
                               int num,
                               int stat_type,
                               int *pvno)
{
  switch (stat_type) {
    default:
    case STAT_T:
      return (cvector_compute_t_test(c1_mean, c1_var, c2_mean, c2_var, num_class1, num_class2, pvals, num, pvno));
      break;
    case STAT_MEAN:
      return (cvector_compute_mean_diff(c1_mean, c2_mean, pvals, num, pvno));
  }
  return (NO_ERROR);
}

int cvector_scalar_mul(float *v, float m, int num)
{
  int i;

  for (i = 0; i < num; i++) v[i] *= m;

  return (NO_ERROR);
}

int cvector_set(float *v, float val, int num)
{
  int i;

  for (i = 0; i < num; i++) v[i] = val;

  return (NO_ERROR);
}
