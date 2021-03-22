/**
 * @brief utilities for reading/writing a Connectome MATrix structure
 *
 * Reading and writing and utilities for the Connectome Matrix (CMAT)
 * structure.
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

/*-----------------------------------------------------
  INCLUDE FILES
  -------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>

#include "cmat.h"
#include "error.h"

CMAT *CMATread(const char *fname)
{
  CMAT *cmat;
  int nlabels, i, j, ind1, ind2;
  FILE *fp;

  fp = fopen(fname, "r");

  if (!fp) {
    ErrorReturn(NULL, (ERROR_NOFILE, "CMATread(%s): could not open file", fname));
  }

  if (fscanf(fp, "CMAT - %d\n", &nlabels) != 1) {
    ErrorPrintf(ERROR_BAD_FILE, "CMATread(%s): could not read parameter(s)", fname);
  }
  cmat = CMATalloc(nlabels, NULL);

  for (i = 0; i < nlabels; i++) {
    if (fscanf(fp, "%d\n", &cmat->labels[i]) != 1) {
      ErrorPrintf(ERROR_BAD_FILE, "CMATread(%s): could not read parameter(s)", fname);
    }
  }

  for (i = 0; i < cmat->nlabels - 1; i++) {
    for (j = i + 1; j < cmat->nlabels; j++) {
      if (fscanf(fp, "%lf", &cmat->weights[i][j]) != 1) {
        ErrorPrintf(ERROR_BAD_FILE, "CMATread(%s): could not read parameter(s)", fname);
      }
    }
    if (fscanf(fp, "\n") != 1) {
      ErrorPrintf(ERROR_BAD_FILE, "CMATread(%s): could not read parameter(s)", fname);
    }
  }
  for (i = 0; i < cmat->nlabels - 1; i++) {
    for (j = i + 1; j < cmat->nlabels; j++) {
      if (fscanf(fp, "%d %d\n", &ind1, &ind2) != 2) {
        ErrorPrintf(ERROR_BAD_FILE, "CMATread(%s): could not read parameter(s)", fname);
      }
      if (feof(fp)) {
        break;
      }
      cmat->splines[ind1][ind2] = LabelReadFrom(NULL, fp);
      if (cmat->coords == LABEL_COORDS_NONE) {
        cmat->coords = cmat->splines[ind1][ind2]->coords;
        printf("reading cmat in coords %d\n", cmat->coords);
      }
      if (feof(fp)) break;
    }
    if (feof(fp)) break;
  }
  fclose(fp);
  return (cmat);
}

int CMATwrite(CMAT *cmat, const char *fname)
{
  FILE *fp;
  int i, j;

  fp = fopen(fname, "w");

  fprintf(fp, "CMAT - %d\n", cmat->nlabels);
  for (i = 0; i < cmat->nlabels; i++) fprintf(fp, "%d\n", cmat->labels[i]);

  for (i = 0; i < cmat->nlabels - 1; i++) {
    for (j = i + 1; j < cmat->nlabels; j++) fprintf(fp, "%f", cmat->weights[i][j]);
    fprintf(fp, "\n");
  }

  for (i = 0; i < cmat->nlabels - 1; i++)
    for (j = i + 1; j < cmat->nlabels; j++) {
      if (cmat->splines[i][j] == NULL) continue;
      fprintf(fp, "%d %d\n", i, j);
      LabelWriteInto(cmat->splines[i][j], fp);
    }
  fclose(fp);
  return (NO_ERROR);
}

CMAT *CMATalloc(int nlabels, int *labels)
{
  CMAT *cmat;
  int i;

  cmat = (CMAT *)calloc(1, sizeof(CMAT));
  if (cmat == NULL) ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat", nlabels);

  cmat->coords = LABEL_COORDS_NONE;
  cmat->nlabels = nlabels;
  cmat->labels = (int *)calloc(nlabels, sizeof(int));
  if (cmat->labels == NULL) ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat->labels", nlabels);
  cmat->splines = (LABEL ***)calloc(nlabels, sizeof(LABEL **));
  if (cmat->splines == NULL) ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat->splines", nlabels);
  cmat->weights = (double **)calloc(nlabels, sizeof(double *));
  if (cmat->weights == NULL) ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat->weights", nlabels);
  for (i = 0; i < nlabels; i++) {
    if (labels) cmat->labels[i] = labels[i];
    cmat->splines[i] = (LABEL **)calloc(nlabels, sizeof(LABEL *));
    cmat->weights[i] = (double *)calloc(nlabels, sizeof(double));
    if (cmat->weights[i] == NULL)
      ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat->weights[%d]", nlabels, i);
    if (cmat->splines[i] == NULL)
      ErrorExit(ERROR_NOMEMORY, "CMATalloc(%d): could not allocate cmat->splines[%d]", nlabels, i);
  }

  return (cmat);
}

int CMATfree(CMAT **pcmat)
{
  CMAT *cmat;
  int i, j;

  cmat = *pcmat;
  *pcmat = NULL;

  free(cmat->labels);
  for (i = 0; i < cmat->nlabels; i++) {
    for (j = i + 1; j < cmat->nlabels; j++) {
      if (cmat->splines[i]) LabelFree(&cmat->splines[i][j]);
    }
    free(cmat->splines[i]);
    free(cmat->weights[i]);
  }

  free(cmat->splines);
  free(cmat->weights);
  free(cmat);
  return (NO_ERROR);
}
CMAT *CMATtransform(CMAT *csrc, TRANSFORM *xform, MRI *mri_src, MRI *mri_dst, CMAT *cdst)
{
  int i, j;

  cdst = CMATalloc(csrc->nlabels, csrc->labels);
  for (i = 0; i < csrc->nlabels; i++) {
    cdst->weights[i] = csrc->weights[i];
    for (j = i + 1; j < csrc->nlabels; j++) {
      if (csrc->splines[i][j]) {
        cdst->splines[i][j] = LabelTransform(csrc->splines[i][j], xform, mri_src, NULL);
      }
    }
  }

  return (cdst);
}

int CMATtoVoxel(CMAT *cmat, MRI *mri)
{
  int i, j;

  cmat->coords = LABEL_COORDS_VOXEL;
  for (i = 0; i < cmat->nlabels; i++) {
    cmat->weights[i] = cmat->weights[i];
    for (j = i + 1; j < cmat->nlabels; j++) {
      if (cmat->splines[i][j]) {
        LabelToVoxel(cmat->splines[i][j], mri, cmat->splines[i][j]);
      }
    }
  }

  return (NO_ERROR);
}
int CMATtoTKreg(CMAT *cmat, MRI *mri)
{
  int i, j;

  cmat->coords = LABEL_COORDS_TKREG_RAS;
  for (i = 0; i < cmat->nlabels; i++) {
    cmat->weights[i] = cmat->weights[i];
    for (j = i + 1; j < cmat->nlabels; j++) {
      if (cmat->splines[i][j]) {
        LabelFromScannerRAS(cmat->splines[i][j], mri, cmat->splines[i][j]);
      }
    }
  }

  return (NO_ERROR);
}
int CMATtoScannerRAS(CMAT *cmat, MRI *mri)
{
  int i, j;

  cmat->coords = LABEL_COORDS_SCANNER_RAS;
  for (i = 0; i < cmat->nlabels; i++) {
    cmat->weights[i] = cmat->weights[i];
    for (j = i + 1; j < cmat->nlabels; j++) {
      if (cmat->splines[i][j]) {
        LabelToScannerRAS(cmat->splines[i][j], mri, cmat->splines[i][j]);
      }
    }
  }

  return (NO_ERROR);
}
