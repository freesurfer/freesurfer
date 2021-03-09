/**
 * @brief utilities for computing flash intensities from T1/PD pairs
 *
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

#include "flash.h"
#include <math.h>
#include "diag.h"
#include "error.h"
#include "macros.h"

extern const char *Progname;

#define MAX_FLASH_VOLUMES 50
#define MIN_T1 10
#define MAX_T1 10000
#define T1_STEP_SIZE 5
#define T1_TO_INDEX(T1) (nint((T1 - MIN_T1) / T1_STEP_SIZE))

typedef struct
{
  double *flash; /* forward model f(T1,TR,alpha) */
  double TR;
  double alpha;
  double step_size; /* T1 increment between table elements */
  double min_T1;
  double max_T1;
  int size; /* number of elements in flash and norm */
} FLASH_LOOKUP_TABLE, FLT;

static int build_lookup_table(double tr, double flip_angle, double te, double min_T1, double max_T1, double step);
static double lookup_flash_value(double TR, double flip_angle, double PD, double T1);
static FLT *find_lookup_table(double TR, double flip_angle);
static FLT lookup_tables[MAX_FLASH_VOLUMES];
static double *norms = NULL;
static int ntables = 0;

static double FLASHforwardModelLookup(double T1, double PD, double TR, double flip_angle);

double dFlash_dT1(double T1, double PD, double TR, double flip_angle, double TE)
{
  double e1, numer, denom;

  e1 = exp(TR / T1);
  numer = e1 * PD * TR * (cos(flip_angle) - 1) * sin(flip_angle);
  denom = T1 * (e1 - cos(flip_angle));
  denom *= denom;
  if (DZERO(denom)) denom = 0.0001;
  return (numer / denom);
}

double dFlash_dPD(double T1, double PD, double TR, double flip_angle, double TE)
{
  double e1, numer, denom;

  e1 = exp(TR / T1);
  numer = (e1 - 1) * sin(flip_angle);
  denom = e1 - cos(flip_angle);
  if (DZERO(denom)) denom = 0.0001;
  return (numer / denom);
}

double FLASHforwardModel(double T1, double PD, double TR, double flip_angle, double TE)
{
  double FLASH, E1;
  double CFA, SFA;

  CFA = cos(flip_angle);
  SFA = sin(flip_angle);
  E1 = exp(-TR / T1);

  FLASH = PD * SFA;
  if (!DZERO(T1)) FLASH *= (1 - E1) / (1 - CFA * E1);
  return (FLASH);
}
MRI *MRIparameterMapsToFlash(MRI *mri_src, MRI *mri_dst, double *TRs, double *TEs, double *FAs, int nflash)
{
  int x, y, z, n;
  double T1, PD;
  double val;

  if (!mri_dst) mri_dst = MRIallocSequence(mri_src->width, mri_src->height, mri_src->depth, mri_src->type, nflash);

  for (x = 0; x < mri_src->width; x++) {
    for (y = 0; y < mri_src->height; y++) {
      for (z = 0; z < mri_src->depth; z++) {
        MRIsampleVolumeFrame(mri_src, x, y, z, 0, &val);
        T1 = val;
        MRIsampleVolumeFrame(mri_src, x, y, z, 1, &val);
        PD = val;

        for (n = 0; n < nflash; n++) {
          val = FLASHforwardModel(T1, PD, TRs[n], FAs[n], TEs[n]);
          MRISseq_vox(mri_dst, x, y, z, n) = (short)nint(val);
        }
      }
    }
  }

  return (mri_dst);
}

int compute_T1_PD(int nvolumes, float *image_vals, double *TRs, double *FAs, double *TEs, double *pT1, double *pPD)
{
  double best_T1, best_PD, norm_im, norm_pred, sse, T1, pred_vals[MAX_FLASH_VOLUMES], error, upper_T1,
      lower_T1, mid_T1, upper_sse, lower_sse, mid_sse, upper_norm, mid_norm, lower_norm, range;
  // double best_sse;
  int i, j, upper_j, lower_j, mid_j, niter;

  if (!norms) {
    norms = (double *)calloc(nint((MAX_T1 - MIN_T1) / T1_STEP_SIZE) + 1, sizeof(double));
    if (!norms) ErrorExit(ERROR_NOMEMORY, "%s: could not allocate norm table", Progname);
    for (j = 0, T1 = MIN_T1; T1 < MAX_T1; T1 += T1_STEP_SIZE, j++) {
      for (norm_pred = 0.0, i = 0; i < nvolumes; i++) {
        pred_vals[i] = FLASHforwardModelLookup(T1, 1.0, TRs[i], FAs[i]);
        norm_pred += (pred_vals[i] * pred_vals[i]);
      }
      norms[j] = sqrt(norm_pred);
      if (FZERO(norms[j])) {
        printf("norms[%d] is zero!\n", j);
        DiagBreak();
        exit(0);
      }
    }
  }

  for (norm_im = i = 0; i < nvolumes; i++) norm_im += (image_vals[i] * image_vals[i]);
  norm_im = sqrt(norm_im); /* length of image vector */
  if (FZERO(norm_im)) {
    *pT1 = MIN_T1;
    *pPD = MIN_T1;
    return (ERROR_BADPARM);
  }
  for (i = 0; i < nvolumes; i++) image_vals[i] /= norm_im; /* normalize them */

  mid_T1 = (MAX_T1 - MIN_T1) / 2;
  mid_j = T1_TO_INDEX(mid_T1);
  range = (MAX_T1 - MIN_T1) / 2;

  /* compute sse for mid T1 */
  mid_norm = norms[mid_j];
  for (mid_sse = 0.0, i = 0; i < nvolumes; i++) {
    pred_vals[i] = FLASHforwardModelLookup(mid_T1, 1.0, TRs[i], FAs[i]);
    pred_vals[i] /= mid_norm; /* normalize them */
    error = (pred_vals[i] - image_vals[i]);
    mid_sse += (error * error);
  }

  best_T1 = mid_T1;
  best_PD = norm_im / mid_norm;
  // best_sse = mid_sse;
  niter = 0;
  if (FZERO(mid_norm)) {
    printf("mid norm=0 at %d (%2.1f)\n", mid_j, mid_T1);
    DiagBreak();
    exit(0);
  }
  do {
    upper_T1 = mid_T1 + 0.5 * range;
    lower_T1 = mid_T1 - 0.5 * range;
    if (upper_T1 > MAX_T1) upper_T1 = MAX_T1;
    if (lower_T1 < MIN_T1) lower_T1 = MIN_T1;
    upper_j = T1_TO_INDEX(upper_T1);
    lower_j = T1_TO_INDEX(lower_T1);
    upper_norm = norms[upper_j];
    lower_norm = norms[lower_j];
    if (FZERO(upper_norm)) {
      printf("upper norm=0 at %d (%2.1f)\n", upper_j, upper_T1);
      DiagBreak();
      exit(0);
    }
    if (FZERO(lower_norm)) {
      printf("lower norm=0 at %d (%2.1f)\n", lower_j, lower_T1);
      DiagBreak();
      exit(0);
    }
    for (lower_sse = upper_sse = 0.0, i = 0; i < nvolumes; i++) {
      pred_vals[i] = FLASHforwardModelLookup(upper_T1, 1.0, TRs[i], FAs[i]);
      pred_vals[i] /= upper_norm; /* normalize them */
      error = (pred_vals[i] - image_vals[i]);
      upper_sse += (error * error);

      pred_vals[i] = FLASHforwardModelLookup(lower_T1, 1.0, TRs[i], FAs[i]);
      pred_vals[i] /= lower_norm; /* normalize them */
      error = (pred_vals[i] - image_vals[i]);
      lower_sse += (error * error);
    }

    if (lower_sse <= mid_sse && lower_sse <= upper_sse) /* make lower new mid */
    {
      mid_sse = lower_sse;
      mid_j = lower_j;
      mid_norm = lower_norm;
      mid_T1 = lower_T1;
      best_T1 = lower_T1;
      best_PD = norm_im / lower_norm;
      // best_sse = lower_sse;
    }
    else if (upper_sse < mid_sse) /* make upper new mid */
    {
      mid_sse = upper_sse;
      mid_j = upper_j;
      mid_norm = upper_norm;
      mid_T1 = upper_T1;
      best_T1 = upper_T1;
      best_PD = norm_im / upper_norm;
      // best_sse = upper_sse;
    }
    if (!std::isfinite(best_PD)) {
      printf("best_PD is not finite at %d (%2.1f)\n", mid_j, mid_T1);
      DiagBreak();
      exit(0);
    }
    range /= 2;
    niter++;
  } while (upper_j - lower_j > 3);

  for (i = 0; i < nvolumes; i++) image_vals[i] *= norm_im; /* restore them */

  for (sse = 0.0, i = 0; i < nvolumes; i++) {
    pred_vals[i] = FLASHforwardModelLookup(best_T1, best_PD, TRs[i], FAs[i]);
    error = (pred_vals[i] - image_vals[i]);
    sse += (error * error);
  }

  *pT1 = best_T1;
  *pPD = best_PD;
  return (NO_ERROR);
}
static double FLASHforwardModelLookup(double T1, double PD, double TR, double flip_angle)
{
  double FLASH;

  FLASH = lookup_flash_value(TR, flip_angle, PD, T1);
  return (FLASH);
}

static int build_lookup_table(double tr, double flip_angle, double te, double min_T1, double max_T1, double step)
{
  FLT *flt;
  int i;
  double T1;

  flt = find_lookup_table(tr, flip_angle);
  if (flt != NULL) return (NO_ERROR); /* already created one */

  if (ntables >= MAX_FLASH_VOLUMES)
    ErrorExit(ERROR_NOMEMORY, "%s: MAX_FLASH_VOLUMES %d exceeded", Progname, MAX_FLASH_VOLUMES);

  flt = &lookup_tables[ntables++];
  flt->TR = tr;
  flt->alpha = flip_angle;
  flt->step_size = step;

  flt->size = (int)((max_T1 - min_T1) / step) + 1;
  flt->flash = (double *)calloc(flt->size, sizeof(*flt->flash));
  if (!flt->flash) ErrorExit(ERROR_NOMEMORY, "%s: could not allocated %dth lookup table", Progname, ntables);

  flt->min_T1 = min_T1;
  flt->max_T1 = max_T1;

  for (i = 0, T1 = MIN_T1; T1 <= MAX_T1; i++, T1 += step)
    flt->flash[i] = FLASHforwardModel(T1, 1.0, tr, flip_angle, te);
  return (NO_ERROR);
}

static FLT *find_lookup_table(double TR, double flip_angle)
{
  int i;

  for (i = 0; i < ntables; i++)
    if (FEQUAL(lookup_tables[i].TR, TR) && FEQUAL(lookup_tables[i].alpha, flip_angle)) break;

  if (i >= ntables) return (NULL);
  return (&lookup_tables[i]);
}

int FlashBuildLookupTables(int nvolumes, double *TRs, double *FAs, double *TEs)
{
  int n;

  for (n = 0; n < nvolumes; n++) build_lookup_table(TRs[n], FAs[n], TEs[n], MIN_T1, MAX_T1, T1_STEP_SIZE);
  return (NO_ERROR);
}
static double lookup_flash_value(double TR, double flip_angle, double PD, double T1)
{
  int index;
  FLT *flt;
  double FLASH;

  flt = find_lookup_table(TR, flip_angle);
  if (!flt) return (FLASHforwardModel(T1, PD, TR, flip_angle, 3.0));
  /* should propagate TE here! */

  index = T1_TO_INDEX(T1);
  if (index < 0) index = 0;
  if (index >= flt->size) index = flt->size - 1;
  FLASH = PD * flt->flash[index];
  return (FLASH);
}
double FLASHforwardModelT2star(double T1, double PD, double T2star, double TR, double flip_angle, double TE)
{
  double FLASH, E1;
  double CFA, SFA;

  CFA = cos(flip_angle);
  SFA = sin(flip_angle);
  E1 = exp(-TR / T1);

  FLASH = PD * SFA;
  if (!DZERO(T1)) FLASH *= (1 - E1) / (1 - CFA * E1);
  if (!DZERO(T2star)) FLASH *= exp(-TE / T2star);
  return (FLASH);
}
