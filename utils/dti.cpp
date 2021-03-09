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


#include "dti.h"
#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <time.h>
#include <unistd.h>
#include "DICOMRead.h"
#include "diag.h"
#include "fio.h"
#include "fmriutils.h"
#include "fsenv.h"
#include "mri.h"
#include "mri2.h"
#include "utils.h"
#include "version.h"

/* --------------------------------------------- */
int DTIfree(DTI **pdti)
{
  DTI *dti;
  dti = *pdti;
  if (dti->GradFile) free(dti->GradFile);
  if (dti->bValue) MatrixFree(&dti->bValue);
  if (dti->GradDir) MatrixFree(&dti->GradDir);
  if (dti->GradDirNorm) MatrixFree(&dti->GradDirNorm);
  if (dti->B) MatrixFree(&dti->B);
  free(*pdti);
  *pdti = NULL;
  return (0);
}

/*-----------------------------------------------------------------
  DTIparamsFromSiemensAscii() - reads in diffusion parameters from the
  Siemens ASCII header stored in fname. fname may be a siemens dicom
  file or an infodump file as produced by mri_probedicom run on a
  siemens dicom file.
  -----------------------------------------------------------------*/
int DTIparamsFromSiemensAscii(const char *fname, float *bValue, int *nDir, int *nB0)

{
  std::string tag;
  char *pc;

  if (!fio_FileExistsReadable(fname)) {
    printf("ERROR: cannot read %s\n", fname);
    return (1);
  }

  tag = "sDiffusion.alBValue[1]";
  pc = SiemensAsciiTag(fname, tag.c_str(), 0);
  if (pc == NULL) {
    printf("ERROR: cannot extract %s from %s\n", tag.c_str(), fname);
    return (1);
  }
  sscanf(pc, "%f", bValue);
  printf("bValue = %g\n", *bValue);
  free(pc);

  tag = "sWiPMemBlock.alFree[8]";
  pc = SiemensAsciiTag(fname, tag.c_str(), 0);
  if (pc == NULL) {
    printf("ERROR: cannot extract %s from %s\n", tag.c_str(), fname);
    return (1);
  }
  sscanf(pc, "%d", nB0);
  printf("nB0 = %d\n", *nB0);
  free(pc);

  tag = "sDiffusion.lDiffDirections";
  pc = SiemensAsciiTag(fname, tag.c_str(), 0);
  if (pc == NULL) {
    printf("ERROR: cannot extract %s from %s\n", tag.c_str(), fname);
    return (1);
  }
  sscanf(pc, "%d", nDir);
  printf("nDir = %d\n", *nDir);
  free(pc);

  return (0);
}
/*------------------------------------------------------------*/
int DTIloadGradients(DTI *dti, const char *GradFile)
{
  static char tmpstr[2000];
  FILE *fp;
  int c, r, n;
  FSENV *fsenv;

  fsenv = FSENVgetenv();

  if (GradFile) dti->GradFile = strcpyalloc(GradFile);
  if (dti->GradFile == NULL) {
    sprintf(tmpstr, "%s/diffusion/mgh-dti-seqpack/gradient_mgh_dti%02d.gdt", getenv("FREESURFER_HOME"), dti->nDir);
    // If it does not exist, try using %d instead of %02d
    fp = fopen(tmpstr, "r");
    if (fp == NULL)
      sprintf(tmpstr, "%s/diffusion/mgh-dti-seqpack/gradient_mgh_dti%0d.gdt", getenv("FREESURFER_HOME"), dti->nDir);
    else
      fclose(fp);
    dti->GradFile = strcpyalloc(tmpstr);
    printf("GradFile %s\n", dti->GradFile);
  }

  fp = fopen(dti->GradFile, "r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n", dti->GradFile);
    return (1);
  }

  dti->GradDir = MatrixAlloc(dti->nDir + dti->nB0, 3, MATRIX_REAL);

  // Set the first nB0 rows to be all 0s (no gradients)
  for (r = 1; r < dti->nB0 + 1; r++) {
    dti->GradDir->rptr[r][1] = 0;
    dti->GradDir->rptr[r][2] = 0;
    dti->GradDir->rptr[r][3] = 0;
  }

  for (r = dti->nB0 + 1; r <= dti->GradDir->rows; r++) {
    for (c = 1; c <= 3; c++) {
      n = fscanf(fp, "%f", &(dti->GradDir->rptr[r][c]));
      if (n != 1) {
        printf("ERROR: reading gradients from %s\n", dti->GradFile);
        printf("  row = %d, col = %d\n", r - 1, c);
        fclose(fp);
        return (1);
      }
    }
  }
  fclose(fp);

  FSENVfree(&fsenv);
  return (0);
}

/*--------------------------------------------------------*/
DTI *DTIstructFromSiemensAscii(const char *fname)
{
  int err;
  float bval;
  DTI *dti;
  int r, n;

  dti = (DTI *)calloc(sizeof(DTI), 1);

  err = DTIparamsFromSiemensAscii(fname, &bval, &dti->nDir, &dti->nB0);
  if (err) {
    free(dti);
    dti = NULL;
    return (NULL);
  }
  // Set the bValues. First nB0 = 0, next nDir = bval
  dti->bValue = MatrixAlloc(dti->nB0 + dti->nDir, 1, MATRIX_REAL);
  r = 1;
  for (n = 1; n <= dti->nB0; n++) {
    dti->bValue->rptr[r][1] = 0;
    r++;
  }
  for (n = 1; n <= dti->nDir; n++) {
    dti->bValue->rptr[r][1] = bval;
    r++;
  }

  err = DTIloadGradients(dti, NULL);
  if (err) {
    free(dti);
    dti = NULL;
    return (NULL);
  }
  DTInormGradDir(dti);

  DTIdesignMatrix(dti);

  return (dti);
}

/*--------------------------------------------------------*/
int DTInormGradDir(DTI *dti)
{
  int r, c;
  double len, maxlen;

  dti->GradDirNorm = MatrixAlloc(dti->GradDir->rows, dti->GradDir->cols, MATRIX_REAL);

  maxlen = 0;
  for (r = 1; r <= dti->GradDir->rows; r++) {
    len = 0;
    for (c = 1; c <= dti->GradDir->cols; c++) len += pow(dti->GradDir->rptr[r][c], 2.0);
    len = sqrt(len);
    if (maxlen < len) maxlen = len;
  }

  for (r = 1; r <= dti->GradDir->rows; r++) {
    len = 0;
    for (c = 1; c <= dti->GradDir->cols; c++) dti->GradDirNorm->rptr[r][c] = dti->GradDir->rptr[r][c] / maxlen;
  }

  return (0);
}
/*--------------------------------------------------------*/
int DTIdesignMatrix(DTI *dti)
{
  int r, xr;
  double bval;
  MATRIX *g;

  g = dti->GradDirNorm;
  dti->B = MatrixAlloc(g->rows, 7, MATRIX_REAL);

  xr = 1;
  for (r = 1; r <= g->rows; r++) {
    bval = dti->bValue->rptr[r][1];
    dti->B->rptr[xr][1] = bval * pow(g->rptr[r][1], 2.0);
    dti->B->rptr[xr][2] = 2 * bval * g->rptr[r][1] * g->rptr[r][2];
    dti->B->rptr[xr][3] = 2 * bval * g->rptr[r][1] * g->rptr[r][3];

    dti->B->rptr[xr][4] = bval * pow(g->rptr[r][2], 2.0);
    dti->B->rptr[xr][5] = 2 * bval * g->rptr[r][2] * g->rptr[r][3];

    dti->B->rptr[xr][6] = bval * pow(g->rptr[r][3], 2.0);

    dti->B->rptr[xr][7] = 1;
    xr++;
  }
  // MatrixWriteTxt("G.dat",dti->GradDirNorm);
  // MatrixWriteTxt("B.dat",dti->B);

  return (0);
}

/*---------------------------------------------------------*/
MRI *DTIbeta2Tensor(MRI *beta, MRI *mask, MRI *tensor)
{
  int c, r, s;
  double m, v;

  if (beta->nframes < 6) {
    printf("ERROR: beta must have at least 6 frames\n");
    return (NULL);
  }
  if (tensor == NULL) {
    tensor = MRIcloneBySpace(beta, MRI_FLOAT, 9);  // 9 = 3x3
    if (!tensor) return (NULL);
  }
  // should check consistency with spatial

  for (c = 0; c < beta->width; c++) {
    for (r = 0; r < beta->height; r++) {
      for (s = 0; s < beta->depth; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }

        // 0 1 2 --> 0 1 2
        // 1 3 4 --> 3 4 5
        // 2 4 5 --> 6 7 8

        // 0 -> 0
        v = MRIgetVoxVal(beta, c, r, s, 0);
        MRIsetVoxVal(tensor, c, r, s, 0, v);

        // 1 -> 1, 3
        v = MRIgetVoxVal(beta, c, r, s, 1);
        MRIsetVoxVal(tensor, c, r, s, 1, v);
        MRIsetVoxVal(tensor, c, r, s, 3, v);

        // 2 -> 2, 6
        v = MRIgetVoxVal(beta, c, r, s, 2);
        MRIsetVoxVal(tensor, c, r, s, 2, v);
        MRIsetVoxVal(tensor, c, r, s, 6, v);

        // 3 -> 4
        v = MRIgetVoxVal(beta, c, r, s, 3);
        MRIsetVoxVal(tensor, c, r, s, 4, v);

        // 4 -> 5, 7
        v = MRIgetVoxVal(beta, c, r, s, 4);
        MRIsetVoxVal(tensor, c, r, s, 5, v);
        MRIsetVoxVal(tensor, c, r, s, 7, v);

        // 5 -> 8
        v = MRIgetVoxVal(beta, c, r, s, 5);
        MRIsetVoxVal(tensor, c, r, s, 8, v);
      }
    }
  }

  return (tensor);
}
/*---------------------------------------------------------*/
int DTItensor2Eig(MRI *tensor, MRI *mask, MRI **evals, MRI **evec1, MRI **evec2, MRI **evec3)
{
  int c, r, s, a, b, n;
  double m;
  MATRIX *T, *Evec;
  float eval[3];

  if (tensor->nframes != 9) {
    printf("ERROR: tensor must have 9 frames\n");
    return (1);
  }
  if (*evals == NULL) {
    *evals = MRIcloneBySpace(tensor, MRI_FLOAT, 3);
    if (!*evals) return (1);
  }
  if (*evec1 == NULL) {
    *evec1 = MRIcloneBySpace(tensor, MRI_FLOAT, 3);
    if (!*evec1) return (1);
  }
  if (*evec2 == NULL) {
    *evec2 = MRIcloneBySpace(tensor, MRI_FLOAT, 3);
    if (!*evec2) return (1);
  }
  if (*evec3 == NULL) {
    *evec3 = MRIcloneBySpace(tensor, MRI_FLOAT, 3);
    if (!*evec3) return (1);
  }
  // should check consistency with spatial

  T = MatrixAlloc(3, 3, MATRIX_REAL);
  Evec = MatrixAlloc(3, 3, MATRIX_REAL);

  for (c = 0; c < tensor->width; c++) {
    for (r = 0; r < tensor->height; r++) {
      for (s = 0; s < tensor->depth; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }

        // Load up the tensor into a matrix struct
        // 0 1 2
        // 3 4 5
        // 6 7 8
        n = 0;
        for (a = 1; a <= 3; a++) {
          for (b = 1; b <= 3; b++) {
            T->rptr[a][b] = MRIgetVoxVal(tensor, c, r, s, n);
            n++;
          }
        }

        /* Do eigen-decomposition */
        MatrixEigenSystem(T, eval, Evec);
        DTIsortEV(eval, Evec);

        for (a = 0; a < 3; a++) {
          MRIsetVoxVal(*evals, c, r, s, a, eval[a]);
          MRIsetVoxVal(*evec1, c, r, s, a, Evec->rptr[a + 1][1]);
          MRIsetVoxVal(*evec2, c, r, s, a, Evec->rptr[a + 1][2]);
          MRIsetVoxVal(*evec3, c, r, s, a, Evec->rptr[a + 1][3]);
        }

      }  // slice
    }    // row
  }      // col

  MatrixFree(&T);
  MatrixFree(&Evec);

  return (0);
}

/*---------------------------------------------------------
  DTIsortEV() - sorts the eigenvalues and eigenvectors from
  max to min.
  ---------------------------------------------------------*/
int DTIsortEV(float *EigVals, MATRIX *EigVecs)
{
  int r;
  static MATRIX *EigVecsTmp = NULL;
  static float EigValsTmp[3];

  for (r = 0; r < 3; r++) EigValsTmp[r] = EigVals[r];
  EigVecsTmp = MatrixCopy(EigVecs, EigVecsTmp);

  if (EigVals[0] > EigVals[1] && EigVals[0] > EigVals[2]) {
    // 1st is max
    if (EigVals[1] > EigVals[2]) {
      // 1st > 2nd > 3rd -- nothing to do
      return (0);
    }
    else {
      // 1st > 3rd > 2nd -- swap 2nd and 3rd cols
      for (r = 1; r <= 3; r++) {
        EigVecs->rptr[r][2] = EigVecsTmp->rptr[r][3];
        EigVecs->rptr[r][3] = EigVecsTmp->rptr[r][2];
      }
      EigVals[2 - 1] = EigValsTmp[3 - 1];
      EigVals[3 - 1] = EigValsTmp[2 - 1];
      return (0);
    }
  }

  if (EigVals[1] > EigVals[0] && EigVals[1] > EigVals[2]) {
    // 2nd is max
    if (EigVals[0] > EigVals[2]) {
      // 2nd > 1st > 3rd -- swap 1st and 2nd
      for (r = 1; r <= 3; r++) {
        EigVecs->rptr[r][1] = EigVecsTmp->rptr[r][2];
        EigVecs->rptr[r][2] = EigVecsTmp->rptr[r][1];
      }
      EigVals[1 - 1] = EigValsTmp[2 - 1];
      EigVals[2 - 1] = EigValsTmp[1 - 1];
      return (0);
    }
    else {
      // 2nd > 3rd > 1st
      for (r = 1; r <= 3; r++) {
        EigVecs->rptr[r][1] = EigVecsTmp->rptr[r][2];
        EigVecs->rptr[r][2] = EigVecsTmp->rptr[r][3];
        EigVecs->rptr[r][3] = EigVecsTmp->rptr[r][1];
      }
      EigVals[1 - 1] = EigValsTmp[2 - 1];
      EigVals[2 - 1] = EigValsTmp[3 - 1];
      EigVals[3 - 1] = EigValsTmp[1 - 1];
      return (0);
    }
  }

  // 3rd is max if it gets here
  if (EigVals[0] > EigVals[1]) {
    // 3rd > 1st > 2nd
    for (r = 1; r <= 3; r++) {
      EigVecs->rptr[r][1] = EigVecsTmp->rptr[r][3];
      EigVecs->rptr[r][2] = EigVecsTmp->rptr[r][1];
      EigVecs->rptr[r][3] = EigVecsTmp->rptr[r][2];
    }
    EigVals[1 - 1] = EigValsTmp[3 - 1];
    EigVals[2 - 1] = EigValsTmp[1 - 1];
    EigVals[3 - 1] = EigValsTmp[2 - 1];
    return (0);
  }
  else {
    // 3rd > 2nd > 1st
    for (r = 1; r <= 3; r++) {
      EigVecs->rptr[r][1] = EigVecsTmp->rptr[r][3];
      EigVecs->rptr[r][2] = EigVecsTmp->rptr[r][2];
      EigVecs->rptr[r][3] = EigVecsTmp->rptr[r][1];
    }
    EigVals[1 - 1] = EigValsTmp[3 - 1];
    EigVals[2 - 1] = EigValsTmp[2 - 1];
    EigVals[3 - 1] = EigValsTmp[1 - 1];
    return (0);
  }

  printf("DTIsortEV(): ERROR: should never get here\n");
  for (r = 1; r <= 3; r++) printf("%g ", EigValsTmp[r]);
  printf("\n");

  return (1);
}

/*---------------------------------------------------------
  DTIbeta2LowB() - computes exp(-v) where v is the offset
  regression parameter. This should be the average volume
  when bvalue=0.
  ---------------------------------------------------------*/
MRI *DTIbeta2LowB(MRI *beta, MRI *mask, MRI *lowb)
{
  int c, r, s;
  double m, v;

  if (beta->nframes < 7) {
    printf("ERROR: beta must have at least 7 frames\n");
    return (NULL);
  }
  if (lowb == NULL) {
    lowb = MRIcloneBySpace(beta, MRI_FLOAT, 1);
    if (!lowb) return (NULL);
  }
  // should check consistency with spatial

  for (c = 0; c < beta->width; c++) {
    for (r = 0; r < beta->height; r++) {
      for (s = 0; s < beta->depth; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        v = MRIgetVoxVal(beta, c, r, s, 6);
        MRIsetVoxVal(lowb, c, r, s, 0, exp(-v));
      }
    }
  }

  return (lowb);
}
/*---------------------------------------------------------
  DTItensor2ADC() - computes apparent diffusion coefficient
  as the trace/3.
  ---------------------------------------------------------*/
MRI *DTItensor2ADC(MRI *tensor, MRI *mask, MRI *adc)
{
  int c, r, s;
  double m, v1, v2, v3, vadc;

  if (tensor->nframes != 9) {
    printf("ERROR: tensor must have at least 9 frames\n");
    return (NULL);
  }
  if (adc == NULL) {
    adc = MRIcloneBySpace(tensor, MRI_FLOAT, 1);
    if (!adc) return (NULL);
  }
  // should check consistency with spatial

  for (c = 0; c < tensor->width; c++) {
    for (r = 0; r < tensor->height; r++) {
      for (s = 0; s < tensor->depth; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        v1 = MRIgetVoxVal(tensor, c, r, s, 0);
        v2 = MRIgetVoxVal(tensor, c, r, s, 4);
        v3 = MRIgetVoxVal(tensor, c, r, s, 8);
        vadc = (v1 + v2 + v3) / 3;
        MRIsetVoxVal(adc, c, r, s, 0, vadc);
      }
    }
  }

  return (adc);
}
/*------------------------------------------------------------*/
MRI *DTIeigvals2FA(MRI *evals, MRI *mask, MRI *FA)
{
  int c, r, s;
  double m, v1, v2, v3, vmean, vsse, vnorm, v;

  if (evals->nframes != 3) {
    printf("ERROR: evals must have 3 frames\n");
    return (NULL);
  }
  if (FA == NULL) {
    FA = MRIcloneBySpace(evals, MRI_FLOAT, 1);
    if (!FA) return (NULL);
  }
  // should check consistency with spatial

  for (c = 0; c < evals->width; c++) {
    for (r = 0; r < evals->height; r++) {
      for (s = 0; s < evals->depth; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        v1 = MRIgetVoxVal(evals, c, r, s, 0);
        v2 = MRIgetVoxVal(evals, c, r, s, 1);
        v3 = MRIgetVoxVal(evals, c, r, s, 2);
        vmean = (v1 + v2 + v3) / 3.0;
        vsse = pow(v1 - vmean, 2.0) + pow(v2 - vmean, 2.0) + pow(v3 - vmean, 2.0);
        vnorm = pow(v1, 2.0) + pow(v2, 2.0) + pow(v3, 2.0);
        v = sqrt(1.5 * vsse / vnorm);  // correct formula?
        MRIsetVoxVal(FA, c, r, s, 0, v);
      }
    }
  }

  return (FA);
}
/*------------------------------------------------------------
  DTIeigvals2RA() - relative anisotropy
  ------------------------------------------------------------*/
MRI *DTIeigvals2RA(MRI *evals, MRI *mask, MRI *RA)
{
  int c, r, s;
  double m, v1, v2, v3, vmean, vsse, v;

  if (evals->nframes != 3) {
    printf("ERROR: evals must have 3 frames\n");
    return (NULL);
  }
  if (RA == NULL) {
    RA = MRIcloneBySpace(evals, MRI_FLOAT, 1);
    if (!RA) return (NULL);
  }
  // should check consistency with spatial

  for (c = 0; c < evals->width; c++) {
    for (r = 0; r < evals->height; r++) {
      for (s = 0; s < evals->depth; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        v1 = MRIgetVoxVal(evals, c, r, s, 0);
        v2 = MRIgetVoxVal(evals, c, r, s, 1);
        v3 = MRIgetVoxVal(evals, c, r, s, 2);
        vmean = (v1 + v2 + v3) / 3.0;
        if (vmean != 0) {
          vsse = pow(v1 - vmean, 2.0) + pow(v2 - vmean, 2.0) + pow(v3 - vmean, 2.0);
          v = sqrt(vsse / (3.0 * vmean));
        }
        else
          v = 0;
        MRIsetVoxVal(RA, c, r, s, 0, v);
      }
    }
  }

  return (RA);
}
/*------------------------------------------------------------
  DTIeigvals2VR() - volume ratio measure of anisotropy. Actually,
  1-VR is used so that it increases with anisotropy.
  ------------------------------------------------------------*/
MRI *DTIeigvals2VR(MRI *evals, MRI *mask, MRI *VR)
{
  int c, r, s;
  double m, v1, v2, v3, vmean, v;

  if (evals->nframes != 3) {
    printf("ERROR: evals must have 3 frames\n");
    return (NULL);
  }
  if (VR == NULL) {
    VR = MRIcloneBySpace(evals, MRI_FLOAT, 1);
    if (!VR) return (NULL);
  }
  // should check consistency with spatial

  for (c = 0; c < evals->width; c++) {
    for (r = 0; r < evals->height; r++) {
      for (s = 0; s < evals->depth; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        v1 = MRIgetVoxVal(evals, c, r, s, 0);
        v2 = MRIgetVoxVal(evals, c, r, s, 1);
        v3 = MRIgetVoxVal(evals, c, r, s, 2);
        vmean = (v1 + v2 + v3) / 3.0;
        if (vmean != 0)
          v = 1 - (v1 * v2 * v3) / pow(vmean, 3.0);
        else
          v = 0.0;
        MRIsetVoxVal(VR, c, r, s, 0, v);
      }
    }
  }

  return (VR);
}
/*----------------------------------------------------------------
  DTIfslBValFile() -- saves bvalues in a format that can be
  read in by FSL's dtifit with -b option. They put all the bvalues
  on one line.
  ----------------------------------------------------------------*/
int DTIfslBValFile(DTI *dti, const char *bvalfname)
{
  FILE *fp;
  int n;

  fp = fopen(bvalfname, "w");
  if (!fp) {
    printf("ERROR: opening %s for writing\n", bvalfname);
    return (1);
  }

  for (n = 1; n <= dti->bValue->rows; n++) fprintf(fp, "%f ", dti->bValue->rptr[n][1]);
  fprintf(fp, "\n");
  fclose(fp);

  return (0);
}
/*----------------------------------------------------------------
  DTIfslBVecFile() -- saves directions in a format that can be read in
  by FSL's dtifit with -r option. They put all the gradients on three
  lines (ie, there are 3 rows and nsamples columns).
  ----------------------------------------------------------------*/
int DTIfslBVecFile(DTI *dti, const char *bvecfname)
{
  FILE *fp;
  int n, c;

  fp = fopen(bvecfname, "w");
  if (!fp) {
    printf("ERROR: opening %s for writing\n", bvecfname);
    return (1);
  }

  for (c = 0; c < 3; c++) {
    // for(n=0; n < dti->nB0; n++)  fprintf(fp,"0 ");
    for (n = 0; n < dti->GradDir->rows; n++) fprintf(fp, "%f ", dti->GradDir->rptr[n + 1][c + 1]);
    fprintf(fp, "\n");
  }
  fclose(fp);

  return (0);
}
/*!/
  \fn MRI *DTIsynthDWI(MATRIX *X, MRI *beta, MRI *mask, MRI *synth);
  \brief Computes exp(-X*beta). Currently mask has no effect.
*/
MRI *DTIsynthDWI(MATRIX *X, MRI *beta, MRI *mask, MRI *synth)
{
  synth = fMRImatrixMultiply(beta, X, synth);
  if (synth == NULL) return (NULL);
  synth = MRIexp(synth, 1, -1, NULL, synth);
  return (synth);
}
/*----------------------------------------------------------*/
/*!/
  \fn MRI *DTIivc(MRI *evec, MRI *mask, MRI *ivc)
  \brief Computes intervoxel coherence
*/
MRI *DTIivc(MRI *evec, MRI *mask, MRI *ivc)
{
  int c, r, s, f, dc, dr, ds, err;
  double v1, v2, vsum, angle, anglesum, m;
  int nhits;

  if (ivc == NULL) {
    ivc = MRIcloneBySpace(evec, MRI_FLOAT, 1);
    if (ivc == NULL) {
      printf("ERROR: DTIivc: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(evec, ivc);
  }
  else {
    err = MRIdimMismatch(evec, ivc, 0);
    if (err) {
      printf("ERROR: DTIivc(): output dimension mismatch (%d)\n", err);
      return (NULL);
    }
    if (ivc->type != MRI_FLOAT) {
      printf("ERROR: DTIivc(): structure passed is not MRI_FLOAT\n");
      return (NULL);
    }
  }

  for (c = 1; c < ivc->width - 1; c++) {
    for (r = 1; r < ivc->height - 1; r++) {
      for (s = 1; s < ivc->depth - 1; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        anglesum = 0;
        nhits = 0;
        for (dc = -1; dc < 2; dc++) {
          for (dr = -1; dr < 2; dr++) {
            for (ds = -1; ds < 2; ds++) {
              if (dc == 0 && dr == 0 && ds == 0) continue;
              if (fabs(dc) + fabs(dr) + fabs(ds) == 3) continue;
              // if (fabs(dc)+ fabs(dr)+ fabs(ds)== 2) continue;
              // skip corners?
              if (mask) {
                m = MRIgetVoxVal(mask, c + dc, r + dr, s + ds, 0);
                if (m < 0.5) continue;
              }
              vsum = 0;
              for (f = 0; f < 3; f++) {
                v1 = MRIgetVoxVal(evec, c, r, s, f);
                v2 = MRIgetVoxVal(evec, c + dc, r + dr, s + ds, f);
                vsum += v1 * v2;
              }

              if (fabs(vsum) > 1) vsum = 1;
              angle = acos(fabs(vsum));
              anglesum += angle;
              // printf("%7.4lf %7.4lf  %7.4lf\n",vsum,angle,anglesum);
              nhits++;
            }  // ds
          }    // dr
        }      // ds
        if (nhits > 0) MRIsetVoxVal(ivc, c, r, s, 0, (M_PI / 2 - anglesum / nhits) / (M_PI / 2));
        // exit(1);
      }
    }
  }

  return (ivc);
}
/*!/
  \fn MRI *DTIloadBValues(const char *bvalfile)
  \brief Loads in bvalues from text file. It does not matter
    whether they are all on the same line or not.
*/
MATRIX *DTIloadBValues(const char *bvalfile)
{
  FILE *fp;
  double b;
  MATRIX *bvals;

  printf("Loading BValues from %s\n", bvalfile);
  fp = fopen(bvalfile, "r");
  if (fp == NULL) {
    printf("ERROR: DTIloadBValues: could not load %s\n", bvalfile);
    return (NULL);
  }

  // First just go through and count them
  int nbvalues = 0;
  while (!feof(fp)) {
    if (fscanf(fp, "%lf", &b) == 1) nbvalues++;
  }
  fclose(fp);
  std::cout << "Found " << nbvalues << " bvalues" << std::endl;
  if (nbvalues == 0) {
    fs::error() << "DTIloadBValues: no bvalues found in " << bvalfile;
    return nullptr;
  }

  // Alloc the matrix
  bvals = MatrixAlloc(nbvalues, 1, MATRIX_REAL);

  // Now read them in
  fp = fopen(bvalfile, "r");
  nbvalues = 0;
  while (!feof(fp)) {
    if (fscanf(fp, "%lf", &b) == 1) {
      bvals->rptr[nbvalues + 1][1] = b;
      nbvalues++;
    }
  }
  fclose(fp);

  // MatrixPrint(stdout,bvals);

  return (bvals);
}
/*---------------------------------------------------------------------*/
int DTIwriteBValues(MATRIX *bvals, const char *bvalfile)
{
  FILE *fp;
  int n;

  fp = fopen(bvalfile, "w");
  if (fp == NULL) {
    printf("ERROR: DTIwriteBValues(): could not open %s for writing\n", bvalfile);
    return (1);
  }
  for (n = 1; n < bvals->rows + 1; n++) fprintf(fp, "%f\n", bvals->rptr[n][1]);
  fclose(fp);
  return (0);
}

/*!/
  \fn MRI *DTIloadBVectors(const char *bvecfile)
  \brief Loads in gradient directions from text file. Each line
    has a different 3-component vector (not the same as FSL).
*/
MATRIX *DTIloadBVectors(const char *bvecfile)
{
  FILE *fp;
  double gx, gy, gz;
  int isFSL;
  MATRIX *bvecs;

  // Could add something to autodetect FSL format
  // which puts all gx on one line, then gy on the next, etc

  printf("Loading BVectors from %s\n", bvecfile);
  fp = fopen(bvecfile, "r");
  if (fp == NULL) {
    printf("ERROR: DTIloadBValues: could not load %s\n", bvecfile);
    return (NULL);
  }

  // First just go through and count them
  int nbvecs = 0;
  while (!feof(fp)) {
    if (fscanf(fp, "%lf %lf %lf", &gx, &gy, &gz) == 3) nbvecs++;
  }
  fclose(fp);
  std::cout << "Found " << nbvecs << " bvectors" << std::endl;
  if (nbvecs == 0) {
    fs::error() << "DTIloadBVectors: no bvectors found in " << bvecfile;
    return nullptr;
  }

  // Alloc the matrix
  bvecs = MatrixAlloc(nbvecs,3,MATRIX_REAL);

  isFSL = DTIisFSLBVec(bvecfile);

  // Now read them in
  fp = fopen(bvecfile,"r");
  if(!isFSL){
    printf("Detected BVec file as MGH formatted\n");
    nbvecs = 0;
    fscanf(fp,"%lf %lf %lf",&gx,&gy,&gz);
    while(!feof(fp)){
      bvecs->rptr[nbvecs+1][1] = gx;
      bvecs->rptr[nbvecs+1][2] = gy;
      bvecs->rptr[nbvecs+1][3] = gz;
      fscanf(fp,"%lf %lf %lf",&gx,&gy,&gz);
      nbvecs ++;
    }
  } else {
    printf("Detected BVec file as FSL formatted\n");
    for(nbvecs = 0; nbvecs < bvecs->rows; nbvecs++)
      fscanf(fp,"%f",&(bvecs->rptr[nbvecs+1][1]));
    for(nbvecs = 0; nbvecs < bvecs->rows; nbvecs++)
      fscanf(fp,"%f",&(bvecs->rptr[nbvecs+1][2]));
    for(nbvecs = 0; nbvecs < bvecs->rows; nbvecs++)
      fscanf(fp,"%f",&(bvecs->rptr[nbvecs+1][3]));
  }
  fclose(fp);

  // MatrixPrint(stdout,bvecs);

  return (bvecs);
}
/*---------------------------------------------------------------------*/
int DTIwriteBVectors(MATRIX *bvecs, const char *bvecfile)
{
  FILE *fp;
  int n;

  fp = fopen(bvecfile, "w");
  if (fp == NULL) {
    printf("ERROR: DTIwriteBVectors(): could not open %s for writing\n", bvecfile);
    return (1);
  }
  for (n = 1; n < bvecs->rows + 1; n++)
    fprintf(fp, "%16.14f %16.14f %16.14f\n", bvecs->rptr[n][1], bvecs->rptr[n][2], bvecs->rptr[n][3]);
  fclose(fp);
  return (0);
}

/*--------------------------------------------------------*/
DTI *DTIstructFromBFiles(const char *bvalfile, const char *bvecfile)
{
  MATRIX *bvals, *bvecs;
  DTI *dti;

  bvals = DTIloadBValues(bvalfile);
  if (bvals == NULL) return (NULL);

  bvecs = DTIloadBVectors(bvecfile);
  if (bvecs == NULL) return (NULL);

  if (bvals->rows != bvecs->rows) {
    printf("ERROR: DTIstructFromBFiles(): dimension mismatch\n");
    printf(" %s has %d rows, %s has %d rows\n", bvalfile, bvals->rows, bvecfile, bvecs->rows);
    return (NULL);
  }

  dti = (DTI *)calloc(sizeof(DTI), 1);
  dti->bValue = bvals;
  dti->GradDir = bvecs;

  DTInormGradDir(dti);
  DTIdesignMatrix(dti);

  // printf("B --------------\n");
  // MatrixPrint(stdout,dti->B);

  return (dti);
}

/*-------------------------------------------------------------------------*/
//  ep_bX#N    X = bvalue     N = nth acq for that bvalue
int DTIparsePulseSeqName(const char *pulseseq, double *bValue, int *nthDirection)
{
  int n;
  const char *pc;
  char tmpstr[100];

  if (strlen(pulseseq) < 7) return (1);

  pc = &(pulseseq[4]);
  n = 0;
  while (*pc != '#') {
    if (*pc == '\0') return (1);
    tmpstr[n] = *pc;
    n++;
    pc++;
  }
  tmpstr[n] = '\0';
  // printf("bvstring %s\n",tmpstr);
  sscanf(tmpstr, "%lf", bValue);

  pc++;
  n = 0;
  while (*pc != '\0') {
    tmpstr[n] = *pc;
    n++;
    pc++;
  }
  tmpstr[n] = '\0';
  // printf("ndstring %s\n",tmpstr);
  sscanf(tmpstr, "%d", nthDirection);

  // printf("bValue = %g, nthDir = %d\n",*bValue,*nthDirection);
  return (0);
}

/*-----------------------------------------------------------*/
/*!/
  \fn int DTIisFSLBVec(const char *fname)
  \brief Detects whether a bvec file is FSL format. The FSL
  format is to have three rows, and each col is a different
  direction. Works by determining whether the first line
  has more than 3 elements. Returns -1 if error, 1 if FSL,
  0 if not FSL.
*/
int DTIisFSLBVec(const char *fname)
{
  FILE *fp;
  char tmpstr[10000], *s;
  float f;
  int nread;

  fp = fopen(fname, "r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n", fname);
    return (-1);
  }

  s = fgets(tmpstr, 9999, fp);
  fclose(fp);

  nread = sscanf(s, "%f %f %f %f", &f, &f, &f, &f);
  if (nread == 3) return (0);
  return (1);
}
/*-----------------------------------------------------------*/
/*!/
  \fn MRI *DTIradialDiffusivity(MRI *evals, MRI *mask, MRI *RD)
  \brief Computes radial diffusivity, which is the average
  of the 2nd and 3rd eigenvalues.
*/
MRI *DTIradialDiffusivity(MRI *evals, MRI *mask, MRI *RD)
{
  int c, r, s;
  double m, v2, v3, vmean;

  if (evals->nframes != 3) {
    printf("ERROR: evals must have 3 frames\n");
    return (NULL);
  }
  if (RD == NULL) {
    RD = MRIcloneBySpace(evals, MRI_FLOAT, 1);
    if (!RD) return (NULL);
  }
  // should check consistency with spatial

  for (c = 0; c < evals->width; c++) {
    for (r = 0; r < evals->height; r++) {
      for (s = 0; s < evals->depth; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        v2 = MRIgetVoxVal(evals, c, r, s, 1);
        v3 = MRIgetVoxVal(evals, c, r, s, 2);
        vmean = (v2 + v3) / 2.0;
        MRIsetVoxVal(RD, c, r, s, 0, vmean);
      }
    }
  }

  return (RD);
}

int DTIbvecChangeSpace(MRI *vol, int desired_bvec_space)
{
  int b, i, f;
  MATRIX *Mdc, *M, *bvec;

  b = desired_bvec_space;
  if (vol->bvecs == NULL) {
    printf("ERROR: DTIbvecChangeSpace(): bvec is NULL\n");
    return (1);
  }

  if (vol->bvec_space == desired_bvec_space) return (0);

  if (b != BVEC_SPACE_SCANNER && b != BVEC_SPACE_VOXEL) {
    printf("ERROR: DTIbvecChangeSpace(): desired_bvec_space = %d, must be %d or %d\n",
           b,
           BVEC_SPACE_SCANNER,
           BVEC_SPACE_VOXEL);
    return (1);
  }

  /* Mdc maps from vox coords to scanner coords.*/
  Mdc = MRImatrixOfDirectionCosines(vol, NULL);
  if (b == BVEC_SPACE_SCANNER) M = Mdc;                     // voxel to scanner
  if (b == BVEC_SPACE_VOXEL) M = MatrixInverse(Mdc, NULL);  // scanner to voxel

  if (Gdiag_no > 0) {
    printf(
        "DTIbvecChangeSpace(): transforming gradient directions from %d to %d\n", vol->bvec_space, desired_bvec_space);
    printf("M------------------------------\n");
    MatrixPrint(stdout, M);
    printf("------------------------------\n");
  }

  // Apply M to each bvec
  bvec = MatrixAlloc(4, 1, MATRIX_REAL);
  bvec->rptr[4][1] = 1;
  for (f = 1; f <= vol->bvecs->rows; f++) {
    for (i = 1; i <= 3; i++) bvec->rptr[i][1] = vol->bvecs->rptr[f][i];
    MatrixMultiplyD(M, bvec, bvec);
    for (i = 1; i <= 3; i++) vol->bvecs->rptr[f][i] = bvec->rptr[i][1];
  }

  MatrixFree(&bvec);
  if (M != Mdc) MatrixFree(&M);
  MatrixFree(&Mdc);

  vol->bvec_space = desired_bvec_space;

  return (0);
}
