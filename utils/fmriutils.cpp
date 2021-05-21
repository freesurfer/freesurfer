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

/*!
  \file fmriutils.c
  \brief Multi-frame utilities

  Things to do:
  1. Add flag to turn use of weight on and off

*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double round(double x);
#include "MRIio_old.h"
#include "diag.h"
#include "float.h"
#include "fmriutils.h"
#include "fsglm.h"
#include "matrix.h"
#include "mri.h"
#include "mri2.h"
#include "numerics.h"
#include "pdf.h"
#include "randomfields.h"
#include "sig.h"
#include "utils.h"
#include "volcluster.h"
#include "romp_support.h"

#ifdef X
#undef X
#endif

/*--------------------------------------------------------*/
MRI *fMRImatrixMultiply(MRI *inmri, MATRIX *M, MRI *outmri)
{
  int c, r, s, fin, fout;
  int nframesout;
  float val;
  int nin, nout;
  int nin0, nout0;
  float *pin = NULL, *pout = NULL;

  if (inmri->nframes != M->cols) {
    printf("ERROR: fMRImatrixMultiply: input dimension mismatch\n");
    return (NULL);
  }
  if (inmri->type != MRI_FLOAT) {
    printf("ERROR: fMRImatrixMultiply: input is not MRI_FLOAT\n");
    return (NULL);
  }

  nframesout = M->rows;
  if (outmri == NULL) {
    outmri = MRIallocSequence(inmri->width, inmri->height, inmri->depth, MRI_FLOAT, nframesout);
    if (outmri == NULL) {
      printf("ERROR: fMRImatrixMultiply: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(inmri, outmri);
  }
  else {
    if (outmri->width != inmri->width || outmri->height != inmri->height || outmri->depth != inmri->depth ||
        outmri->nframes != nframesout) {
      printf("ERROR: fMRImatrixMultiply: output dimension mismatch\n");
      return (NULL);
    }
    if (outmri->type != MRI_FLOAT) {
      printf("ERROR: fMRImatrixMultiply: output is not MRI_FLOAT\n");
      return (NULL);
    }
  }

  MRIclear(outmri);
  for (fout = 0; fout < outmri->nframes; fout++) {
    nout0 = fout * outmri->depth;
    for (fin = 0; fin < inmri->nframes; fin++) {
      val = M->rptr[fout + 1][fin + 1];
      nin0 = fin * inmri->depth;
      for (s = 0; s < outmri->depth; s++) {
        nin = s + nin0;
        nout = s + nout0;
        for (r = 0; r < outmri->height; r++) {
          pin = (float *)inmri->slices[nin][r];
          pout = (float *)outmri->slices[nout][r];
          for (c = 0; c < outmri->width; c++) (*pout++) += val * (*pin++);
        }
      }
    }
  }

  return (outmri);
}

/*!
\fn MRI *fMRIcovariance(MRI *fmri, int Lag, float DOFAdjust, MRI *mask, MRI *covar)
\brief Computes the temporal covariance at Lag at each voxel where mask > 0.5. Does
not remove the mean unless DOFAdjust < 0.
\param fmri - input
\param Lag - covariance lag (0 for simple variance)
\param DOFAdjust - DOF = nframes - DOFAdjust. If DOFAdjust < 0, then removes the mean
and resets DOFAdjust=1.
\param mask - only compute where mask > 0.5 (or everywhere if mask is NULL)
\param covar - output (can be NULL).
*/
MRI *fMRIcovariance(MRI *fmri, int Lag, float DOFAdjust, MRI *mask, MRI *covar)
{
  int RemoveMean = 0;
  int DOF, DOFLag, c, r, s, f;
  double sumv1v2, val1, val2, valmean;
  MRI *mean = NULL;

  if (DOFAdjust < 0) {
    RemoveMean = 1;
    DOFAdjust = 1;
  }

  if (Lag < 0 || Lag >= fmri->nframes) {
    printf("ERROR: fMRIcovariance: Lag out of bound %d \n", Lag);
    return (NULL);
  }

  DOF = fmri->nframes - DOFAdjust;

  DOFLag = DOF - Lag;
  if (DOFLag <= 0) {
    printf("ERROR: fMRIcovariance: DOFLag = %d <= 0 (DOF=%d, Lag=%d)\n", DOFLag, DOF, Lag);
    return (NULL);
  }

  if (covar == NULL) {
    covar = MRIallocSequence(fmri->width, fmri->height, fmri->depth, MRI_FLOAT, 1);
    if (covar == NULL) {
      printf("ERROR: fMRIvariance: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(fmri, covar);
  }

  if (RemoveMean) mean = MRIframeMean(fmri, NULL);

  valmean = 0;
  for (c = 0; c < fmri->width; c++) {
    for (r = 0; r < fmri->height; r++) {
      for (s = 0; s < fmri->depth; s++) {
        if (mask) {
          if (MRIgetVoxVal(mask, c, r, s, 0) < 0.5) {
            MRIFseq_vox(covar, c, r, s, 0) = 0;
            continue;
          }
        }
        if (RemoveMean) valmean = MRIgetVoxVal(mean, c, r, s, 0);
        sumv1v2 = 0;
        for (f = 0; f < fmri->nframes - Lag; f++) {
          val1 = MRIgetVoxVal(fmri, c, r, s, f);
          val2 = MRIgetVoxVal(fmri, c, r, s, f + Lag);
          sumv1v2 += ((val1 - valmean) * (val2 - valmean));
        }
        MRIFseq_vox(covar, c, r, s, 0) = sumv1v2 / DOFLag;
      }
    }
  }

  if (mean) MRIfree(&mean);

  return (covar);
}

/*!
\fn MRI *fMRItemporalAR1(MRI *fmri, float DOFAdjust, MRI *mask, MRI *ar1)
\brief Computes the temporal AR(1) at each voxel where mask > 0.5. Does
not remove the mean unless DOFAdjust < 0.
\param fmri - input
\param DOFAdjust - DOF = nframes - DOFAdjust. If DOFAdjust < 0, then removes the mean
and resets DOFAdjust=1.
\param mask - only compute where mask > 0.5 (or everywhere if mask is NULL)
\param ar1 - output (can be NULL).
*/
MRI *fMRItemporalAR1(MRI *fmri, float DOFAdjust, MRI *mask, MRI *ar1)
{
  int c, r, s;
  double voxvar, voxcovar;
  MRI *var, *covar;

  var = fMRIcovariance(fmri, 0, DOFAdjust, mask, NULL);
  covar = fMRIcovariance(fmri, 1, DOFAdjust, mask, NULL);

  if (ar1 == NULL) {
    ar1 = MRIallocSequence(fmri->width, fmri->height, fmri->depth, MRI_FLOAT, 1);
    if (ar1 == NULL) {
      printf("ERROR: fMRItemporalAR1: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(fmri, ar1);
  }

  for (c = 0; c < fmri->width; c++) {
    for (r = 0; r < fmri->height; r++) {
      for (s = 0; s < fmri->depth; s++) {
        if (mask) {
          if (MRIgetVoxVal(mask, c, r, s, 0) < 0.5) {
            MRIFseq_vox(ar1, c, r, s, 0) = 0;
            continue;
          }
        }
        voxvar = MRIFseq_vox(var, c, r, s, 0);
        if (voxvar == 0) MRIFseq_vox(ar1, c, r, s, 0) = 0;
        voxcovar = MRIFseq_vox(covar, c, r, s, 0);
        MRIFseq_vox(ar1, c, r, s, 0) = voxcovar / voxvar;
      }
    }
  }

  return (ar1);
}

/*--------------------------------------------------------
  fMRIsumSquare() - computes the sum of the squares over the
  frames. If the Update flag is set, then the sum of the
  squares is added to that already in sumsqr. If sumsqr
  is NULL, it will be allocated.
  --------------------------------------------------------*/
MRI *fMRIsumSquare(MRI *fmri, int Update, MRI *sumsqr)
{
  int c, r, s, f, n;
  float v;
  float *pfmri = NULL, *psumsqr = NULL;

  if (sumsqr == NULL) {
    sumsqr = MRIallocSequence(fmri->width, fmri->height, fmri->depth, MRI_FLOAT, 1);
    if (sumsqr == NULL) {
      printf("ERROR: fMRIsumSquare: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(fmri, sumsqr);
  }
  else {
    if (sumsqr->width != fmri->width || sumsqr->height != fmri->height || sumsqr->depth != fmri->depth) {
      printf("ERROR: fMRIsumsqriance: output dimension mismatch\n");
      return (NULL);
    }
    if (sumsqr->type != MRI_FLOAT) {
      printf("ERROR: fMRIsumsqriance: structure passed is not MRI_FLOAT\n");
      return (NULL);
    }
  }

  if (Update) MRIclear(sumsqr);
  n = 0;
  for (f = 0; f < fmri->nframes; f++) {
    for (s = 0; s < fmri->depth; s++) {
      for (r = 0; r < fmri->height; r++) {
        pfmri = (float *)fmri->slices[n][r];
        psumsqr = (float *)sumsqr->slices[s][r];
        for (c = 0; c < fmri->width; c++) {
          v = (*pfmri++);
          (*psumsqr++) += (v * v);
        }
      }
      n++;
    }
  }

  return (sumsqr);
}

/*--------------------------------------------------------------------*/
MRI *fMRIcomputeT(MRI *ces, MATRIX *X, MATRIX *C, MRI *var, MRI *t)
{
  int c, r, s;
  MATRIX *Xt, *XtX, *iXtX, *CiXtX, *Ct, *CiXtXCt;
  float srf, cesval, std;

  if (C->rows != 1) {
    printf("ERROR: fMRIcomputeT: contrast matrix has more than 1 row.\n");
    return (NULL);
  }

  if (ces->nframes != 1) {
    printf(
        "ERROR: fMRIcomputeT: contrast effect size has "
        "more than 1 frame.\n");
    return (NULL);
  }

  if (t == NULL) {
    t = MRIallocSequence(ces->width, ces->height, ces->depth, MRI_FLOAT, 1);
    if (t == NULL) {
      printf("ERROR: fMRIcomputeT: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(ces, t);
  }
  else {
    if (t->width != ces->width || t->height != ces->height || t->depth != ces->depth) {
      printf("ERROR: fMRIcomputeT: output dimension mismatch\n");
      return (NULL);
    }
    if (t->type != MRI_FLOAT) {
      printf("ERROR: fMRIcomputeT: structure passed is not MRI_FLOAT\n");
      return (NULL);
    }
  }

  Xt = MatrixTranspose(X, NULL);
  XtX = MatrixMultiply(Xt, X, NULL);
  iXtX = MatrixInverse(XtX, NULL);
  CiXtX = MatrixMultiplyD(C, iXtX, NULL);
  Ct = MatrixTranspose(C, NULL);
  CiXtXCt = MatrixMultiplyD(CiXtX, Ct, NULL);
  srf = sqrt(CiXtXCt->rptr[1][1]);
  // printf("fMRIcomputeT: srf = %g\n",srf);

  for (c = 0; c < ces->width; c++) {
    for (r = 0; r < ces->height; r++) {
      for (s = 0; s < ces->depth; s++) {
        std = sqrt(MRIgetVoxVal(var, c, r, s, 0));
        if (std == 0)
          MRIFseq_vox(t, c, r, s, 0) = 0;
        else {
          cesval = MRIgetVoxVal(ces, c, r, s, 0);
          MRIFseq_vox(t, c, r, s, 0) = cesval / (srf * std);
        }
      }
    }
  }

  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&CiXtX);
  MatrixFree(&Ct);
  MatrixFree(&CiXtXCt);

  return (t);
}

/*--------------------------------------------------------*/
MRI *fMRIsigT(MRI *t, float DOF, MRI *sig)
{
  int c, r, s, f;
  float tval, sigtval;

  if (sig == NULL) {
    sig = MRIallocSequence(t->width, t->height, t->depth, MRI_FLOAT, t->nframes);
    if (sig == NULL) {
      printf("ERROR: fMRIsigT: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(t, sig);
  }
  else {
    if (t->width != sig->width || t->height != sig->height || t->depth != sig->depth || t->nframes != sig->nframes) {
      printf("ERROR: fMRIsigT: output dimension mismatch\n");
      return (NULL);
    }
    if (sig->type != MRI_FLOAT) {
      printf("ERROR: fMRIsigT: structure passed is not MRI_FLOAT\n");
      return (NULL);
    }
  }

  for (c = 0; c < t->width; c++) {
    for (r = 0; r < t->height; r++) {
      for (s = 0; s < t->depth; s++) {
        for (f = 0; f < t->nframes; f++) {
          tval = MRIFseq_vox(t, c, r, s, f);
          sigtval = sigt(tval, rint(DOF));
          if (tval < 0) sigtval *= -1;
          MRIFseq_vox(sig, c, r, s, f) = sigtval;
        }
      }
    }
  }

  return (sig);
}

/*--------------------------------------------------------------------*/
MRI *fMRIcomputeF(MRI *ces, MATRIX *X, MATRIX *C, MRI *var, MRI *F)
{
  int c, r, s, f, J;
  MATRIX *Xt, *XtX, *iXtX, *CiXtX, *Ct, *CiXtXCt, *iCiXtXCt;
  MATRIX *M, *cesvect, *cesvectt, *voxF;
  float cesval, voxvar;

  if (ces->nframes != C->rows) {
    printf(
        "ERROR: fMRIcomputeT: contrast effect size and contrast matrix "
        "have inconsistent dimensions.\n");
    return (NULL);
  }

  if (F == NULL) {
    F = MRIallocSequence(ces->width, ces->height, ces->depth, MRI_FLOAT, 1);
    if (F == NULL) {
      printf("ERROR: fMRIcomputeF: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(ces, F);
  }
  else {
    if (F->width != ces->width || F->height != ces->height || F->depth != ces->depth) {
      printf("ERROR: fMRIcomputeT: output dimension mismatch\n");
      return (NULL);
    }
    if (F->type != MRI_FLOAT) {
      printf("ERROR: fMRIcomputeT: structure passed is not MRI_FLOAT\n");
      return (NULL);
    }
  }

  Xt = MatrixTranspose(X, NULL);
  XtX = MatrixMultiplyD(Xt, X, NULL);
  iXtX = MatrixInverse(XtX, NULL);
  CiXtX = MatrixMultiplyD(C, iXtX, NULL);
  Ct = MatrixTranspose(C, NULL);
  CiXtXCt = MatrixMultiplyD(CiXtX, Ct, NULL);
  iCiXtXCt = MatrixInverse(CiXtXCt, NULL);
  J = C->rows;
  cesvect = MatrixAlloc(ces->nframes, 1, MATRIX_REAL);
  cesvectt = MatrixAlloc(1, ces->nframes, MATRIX_REAL);
  M = NULL;
  voxF = NULL;

  for (c = 0; c < ces->width; c++) {
    for (r = 0; r < ces->height; r++) {
      for (s = 0; s < ces->depth; s++) {
        voxvar = MRIgetVoxVal(var, c, r, s, 0);
        if (voxvar == 0)
          MRIFseq_vox(F, c, r, s, 0) = 0;
        else {
          for (f = 0; f < ces->nframes; f++) {
            cesval = MRIgetVoxVal(ces, c, r, s, f);
            cesvect->rptr[f + 1][1] = cesval;
            cesvectt->rptr[1][f + 1] = cesval;
          }
          M = MatrixMultiplyD(iCiXtXCt, cesvect, M);
          voxF = MatrixMultiplyD(cesvectt, M, voxF);

          MRIFseq_vox(F, c, r, s, 0) = (voxF->rptr[1][1]) / (J * voxvar);
        }
      }
    }
  }

  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&CiXtX);
  MatrixFree(&Ct);
  MatrixFree(&CiXtXCt);
  MatrixFree(&iCiXtXCt);
  MatrixFree(&cesvect);
  MatrixFree(&cesvectt);
  MatrixFree(&M);

  return (F);
}

/*--------------------------------------------------------*/
// DOF1 = dof of den (same as t DOF)
// DOF2 = dof of num (number of rows in C)
// Note: order is rev relative to fsfast's FTest.m
MRI *fMRIsigF(MRI *F, float DOFDen, float DOFNum, MRI *sig)
{
  int c, r, s, f;
  float Fval, sigFval;

  if (sig == NULL) {
    sig = MRIallocSequence(F->width, F->height, F->depth, MRI_FLOAT, F->nframes);
    if (sig == NULL) {
      printf("ERROR: fMRIsigF: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(F, sig);
  }
  else {
    if (F->width != sig->width || F->height != sig->height || F->depth != sig->depth || F->nframes != sig->nframes) {
      printf("ERROR: fMRIsigF: output dimension mismatch\n");
      return (NULL);
    }
    if (sig->type != MRI_FLOAT) {
      printf("ERROR: fMRIsigF: structure passed is not MRI_FLOAT\n");
      return (NULL);
    }
  }

  for (c = 0; c < F->width; c++) {
    for (r = 0; r < F->height; r++) {
      for (s = 0; s < F->depth; s++) {
        for (f = 0; f < F->nframes; f++) {
          Fval = MRIFseq_vox(F, c, r, s, f);
          sigFval = sc_cdf_fdist_Q(Fval, DOFNum, DOFDen);
          MRIFseq_vox(sig, c, r, s, f) = sigFval;
        }
      }
    }
  }

  return (sig);
}

/*--------------------------------------------------------
  fMRInskip() - skip the first nskip frames
  --------------------------------------------------------*/
MRI *fMRInskip(MRI *inmri, int nskip, MRI *outmri)
{
  int c, r, s, fin, fout;
  int nframesout;
  float val;

  if (inmri->nframes <= nskip) {
    printf("ERROR: fMRInskip: nskip >= nframes\n");
    return (NULL);
  }

  nframesout = inmri->nframes - nskip;
  if (outmri == NULL) {
    outmri = MRIallocSequence(inmri->width, inmri->height, inmri->depth, inmri->type, nframesout);
    if (outmri == NULL) {
      printf("ERROR: fMRInskip: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(inmri, outmri);
  }
  else {
    if (outmri->width != inmri->width || outmri->height != inmri->height || outmri->depth != inmri->depth ||
        outmri->nframes != nframesout) {
      printf("ERROR: fMRInskip: output dimension mismatch\n");
      return (NULL);
    }
    if (outmri->type != inmri->type) {
      printf("ERROR: fMRInskip: structure type mismatch\n");
      return (NULL);
    }
  }

  MRIclear(outmri);
  for (fout = 0; fout < outmri->nframes; fout++) {
    fin = fout + nskip;
    for (s = 0; s < outmri->depth; s++) {
      for (r = 0; r < outmri->height; r++) {
        for (c = 0; c < outmri->width; c++) {
          val = MRIgetVoxVal(inmri, c, r, s, fin);
          MRIsetVoxVal(outmri, c, r, s, fout, val);
        }
      }
    }
  }

  return (outmri);
}

/*--------------------------------------------------------
  fMRIndrop() - drop the last ndrop frames
  --------------------------------------------------------*/
MRI *fMRIndrop(MRI *inmri, int ndrop, MRI *outmri)
{
  int c, r, s, fin, fout;
  int nframesout;
  float val;

  if (inmri->nframes <= ndrop) {
    printf("ERROR: fMRIndrop: ndrop >= nframes\n");
    return (NULL);
  }

  nframesout = inmri->nframes - ndrop;
  if (outmri == NULL) {
    outmri = MRIallocSequence(inmri->width, inmri->height, inmri->depth, inmri->type, nframesout);
    if (outmri == NULL) {
      printf("ERROR: fMRIndrop: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(inmri, outmri);
  }
  else {
    if (outmri->width != inmri->width || outmri->height != inmri->height || outmri->depth != inmri->depth ||
        outmri->nframes != nframesout) {
      printf("ERROR: fMRIndrop: output dimension mismatch\n");
      return (NULL);
    }
    if (outmri->type != inmri->type) {
      printf("ERROR: fMRIndrop: structure type mismatch\n");
      return (NULL);
    }
  }

  MRIclear(outmri);
  for (fout = 0; fout < outmri->nframes; fout++) {
    fin = fout;
    for (s = 0; s < outmri->depth; s++) {
      for (r = 0; r < outmri->height; r++) {
        for (c = 0; c < outmri->width; c++) {
          val = MRIgetVoxVal(inmri, c, r, s, fin);
          MRIsetVoxVal(outmri, c, r, s, fout, val);
        }
      }
    }
  }

  return (outmri);
}

/*--------------------------------------------------------
  fMRIframe() - extract the nth frame. frame is 0-based.
  --------------------------------------------------------*/
MRI *fMRIframe(MRI *inmri, int frame, MRI *outmri)
{
  int c, r, s, nframesout;
  float val;

  if (inmri->nframes <= frame) {
    printf("ERROR: fMRIframe: frame >= nframes\n");
    return (NULL);
  }

  nframesout = 1;
  if (outmri == NULL) {
    outmri = MRIallocSequence(inmri->width, inmri->height, inmri->depth, inmri->type, nframesout);
    if (outmri == NULL) {
      printf("ERROR: fMRIframe: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(inmri, outmri);
    MRIcopyPulseParameters(inmri, outmri);
  }
  else {
    if (outmri->width != inmri->width || outmri->height != inmri->height || outmri->depth != inmri->depth ||
        outmri->nframes != nframesout) {
      printf("ERROR: fMRIframe: output dimension mismatch\n");
      return (NULL);
    }
    if (outmri->type != inmri->type) {
      printf("ERROR: fMRIframe: structure type mismatch\n");
      return (NULL);
    }
  }

  MRIclear(outmri);
  for (s = 0; s < outmri->depth; s++) {
    for (r = 0; r < outmri->height; r++) {
      for (c = 0; c < outmri->width; c++) {
        val = MRIgetVoxVal(inmri, c, r, s, frame);
        MRIsetVoxVal(outmri, c, r, s, 0, val);
      }
    }
  }

  return (outmri);
}
/*
  \fn MRI *fMRIinsertFrame(MRI *srcmri, int srcframe, MRI *fmri, int frame)
  \brief Inserts the specified frame from the source mri into the target frame
  of the fmri. If fmri is NULL, then it is allocated with frame+1 frames.
 */
MRI *fMRIinsertFrame(MRI *srcmri, int srcframe, MRI *fmri, int frame)
{
  int c, r, s;
  double v;

  if (fmri == NULL) {
    fmri = MRIallocSequence(srcmri->width, srcmri->height, srcmri->depth, srcmri->type, frame + 1);
    MRIcopyHeader(srcmri, fmri);
    MRIcopyPulseParameters(srcmri, fmri);
  }
  if (fmri->nframes <= frame) {
    printf("ERROR: fMRIinsertFrame() frame %d => nframes %d\n", frame, fmri->nframes);
    return (NULL);
  }
  if (srcmri->nframes <= srcframe) {
    printf("ERROR: fMRIinsertFrame() srcframe %d => nframes %d\n", srcframe, srcmri->nframes);
    return (NULL);
  }

  for (c = 0; c < srcmri->width; c++) {
    for (r = 0; r < srcmri->height; r++) {
      for (s = 0; s < srcmri->depth; s++) {
        v = MRIgetVoxVal(srcmri, c, r, s, srcframe);
        MRIsetVoxVal(fmri, c, r, s, frame, v);
      }
    }
  }

  return (fmri);
}

/*---------------------------------------------------------------*/
MATRIX *MRItoMatrix(MRI *mri, int c, int r, int s, int Mrows, int Mcols, MATRIX *M)
{
  int mr, mc, f;

  if (M == NULL)
    M = MatrixAlloc(Mrows, Mcols, MATRIX_REAL);
  else {
    if (M->rows != Mrows || M->cols != Mcols) {
      printf("ERROR: Matrix dim mismatch\n");
    }
  }

  if (mri->nframes != Mrows * Mcols) {
    printf("ERROR: MRItoMatrix: MRI frames = %d, does not equal\n", mri->nframes);
    printf("       matrix dim = %dx%d = %d", Mrows, Mcols, Mrows * Mcols);
    return (NULL);
  }

  f = 0;
  for (mr = 1; mr <= Mrows; mr++) {
    for (mc = 1; mc <= Mcols; mc++) {
      M->rptr[mr][mc] = MRIgetVoxVal(mri, c, r, s, f);
      f++;
    }
  }
  return (M);
}

/*---------------------------------------------------------------*/
MATRIX *MRItoSymMatrix(MRI *mri, int c, int r, int s, MATRIX *M)
{
  int mr, mc, f, Msize, nframesexp;

  if (M == NULL) {
    Msize = (int)(round((sqrt(8.0 * mri->nframes + 1.0) - 1.0) / 2.0));
    printf("Msize = %d\n", Msize);
    M = MatrixAlloc(Msize, Msize, MATRIX_REAL);
  }

  nframesexp = M->rows * (M->rows + 1) / 2;
  if (mri->nframes != nframesexp) {
    printf("ERROR: MRItoSymMatrix: MRI frames = %d, does not support sym\n", mri->nframes);
    return (NULL);
  }

  f = 0;
  for (mr = 1; mr <= M->rows; mr++) {
    for (mc = mr; mc <= M->cols; mc++) {
      M->rptr[mr][mc] = MRIgetVoxVal(mri, c, r, s, f);
      M->rptr[mc][mr] = MRIgetVoxVal(mri, c, r, s, f);
      f++;
    }
  }
  return (M);
}

/*---------------------------------------------------------------*/
int MRIfromMatrix(MRI *mri, int c, int r, int s, MATRIX *M, MRI *FrameMask)
{
  int mr, mc, f, nf;

  nf = mri->nframes;
  // Count the number of frames in frame mask
  if (FrameMask != NULL) {
    nf = 0;
    for (f = 1; f <= mri->nframes; f++)
      if (MRIgetVoxVal(FrameMask, c, r, s, f - 1) > 0.5) nf++;
  }

  if (nf != M->rows * M->cols) {
    printf("ERROR: MRIfromMatrix: MRI frames = %d, does not equal\n", nf);
    printf("       matrix dim = %dx%d = %d", M->rows, M->cols, M->rows * M->cols);
    return (1);
  }

  mr = 1;
  mc = 1;
  for (f = 1; f <= mri->nframes; f++) {
    if (FrameMask && MRIgetVoxVal(FrameMask, c, r, s, f - 1) < 0.5) {
      MRIsetVoxVal(mri, c, r, s, f - 1, 0.0);
      continue;
    }
    if (mc > M->cols) {
      mc = 1;
      mr++;
    }
    MRIsetVoxVal(mri, c, r, s, f - 1, M->rptr[mr][mc]);
    mc++;
  }

  // This is the old code without the FrameMask
  // f = 0;
  // for (mr=1; mr <= M->rows; mr++){
  //  for (mc=1; mc <= M->cols; mc++) {
  //    MRIsetVoxVal(mri,c,r,s,f,M->rptr[mr][mc]);
  //    f++;
  //  }
  //}
  return (0);
}

/*---------------------------------------------------------------*/
int MRIfromSymMatrix(MRI *mri, int c, int r, int s, MATRIX *M)
{
  int mr, mc, f, nframesexp;

  nframesexp = M->rows * (M->rows + 1) / 2;
  if (mri->nframes != nframesexp) {
    printf("ERROR: MRIfromSumMatrix: MRI frames = %d, does not equal\n", mri->nframes);
    printf("       matrix dim = %dx%d = %d", M->rows, M->cols, M->rows * M->cols);
    return (1);
  }

  f = 0;
  for (mr = 1; mr <= M->rows; mr++) {
    for (mc = mr; mc <= M->cols; mc++) {
      MRIsetVoxVal(mri, c, r, s, f, M->rptr[mr][mc]);
      f++;
    }
  }
  return (0);
}

/*------------------------------------------------------------------------
  MRInormWeights() - rescales each voxel so that the sum across all
  frames equals nframes. If sqrtFlag=1, then computes the sqrt(w)
  before normalzing.  If invFlag=1, then computes the 1/w before
  normalzing.  The sqrt and inv can be used if the weights are
  variances for WLMS. Can be done in-place. If mask, then ignores
  voxels where mask<0.5. Weights must be >  0.
  *------------------------------------------------------*/
MRI *MRInormWeights(MRI *w, int sqrtFlag, int invFlag, MRI *mask, MRI *wn)
{
  int c, r, s, f;
  double v, vsum, m;

  //-------------------------------------------
  if (wn == NULL) {
    wn = MRIallocSequence(w->width, w->height, w->depth, MRI_FLOAT, w->nframes);
    if (wn == NULL) {
      printf("ERROR: MRInormWeights(): could not alloc weights\n");
      return (NULL);
    }
    MRIcopyHeader(wn, w);
  }

  //-------------------------------------------
  for (c = 0; c < w->width; c++) {
    for (r = 0; r < w->height; r++) {
      for (s = 0; s < w->depth; s++) {
        if (mask != NULL) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }

        // First go through and compute the sum
        vsum = 0;
        for (f = 0; f < w->nframes; f++) {
          v = MRIgetVoxVal(w, c, r, s, f);
          if (v <= 0) {
            printf("ERROR: MRInormWeights: value less than or eq to 0.\n");
            printf("  c=%d, r=%d, s=%d, v=%g\n", c, r, s, v);
            // Should do a free here, I guess
            return (NULL);
          }
          if (sqrtFlag) v = sqrt(v);
          if (invFlag) v = 1 / v;
          vsum += v;
        }

        // So that the sum = nframes
        vsum /= w->nframes;

        // Now rescale
        for (f = 0; f < w->nframes; f++) {
          v = MRIgetVoxVal(w, c, r, s, f);
          if (sqrtFlag) v = sqrt(v);
          if (invFlag) v = 1 / v;
          v = v / vsum;
          MRIsetVoxVal(wn, c, r, s, f, v);
        }
      }
    }
  }

  return (wn);
}

/*---------------------------------------------------------------------
  MRIglmFitAndTest() - fits and tests glm on a voxel-by-voxel basis.
  There are also two other related functions, MRIglmFit() and
  MRIglmTest(), that accomplish the same thing except MRIglmFit() fits
  all the voxels and then MRIglmTest() tests all the voxels.
  MRIglmFitAndTest() fits and tests a voxel before moving on to the
  next voxel. MRIglmFitAndTest() will be computationally more
  efficient.  So why have MRIglmFit() and MRIglmTest()? So that the
  variance can be smoothed between the two if desired.
  --------------------------------------------------------------------*/
int MRIglmFitAndTest(MRIGLM *mriglm)
{
  int c, nc, nr, ns, nf, n;
  long nvoxtot;
  //int c, r, s, n, nc, nr, ns, nf, pctdone;
  //float m, Xcond;
  //long nvoxtot, nthvox;
  GLMMAT *glm = mriglm->glm;

  nc = mriglm->y->width;
  nr = mriglm->y->height;
  ns = mriglm->y->depth;
  nvoxtot = nc * nr * ns;
  nf = mriglm->y->nframes;
  mriglm->nregtot = MRIglmNRegTot(mriglm);

  mriglm->pervoxflag = 0;
  if (mriglm->w != NULL || mriglm->npvr != 0 || mriglm->FrameMask != NULL) mriglm->pervoxflag = 1;

  GLMcMatrices(glm);

  if (mriglm->FrameMask == NULL) {
    GLMallocX(glm, nf, mriglm->nregtot);
    GLMallocY(glm);
    if (mriglm->yffxvar) GLMallocYFFxVar(glm);
  }

  if (!mriglm->pervoxflag) {
    MatrixCopy(mriglm->Xg, glm->X);
    mriglm->XgLoaded = 1;
    GLMxMatrices(glm);
  }

  // If beta has not been allocated, assume that no one has been alloced
  if (mriglm->beta == NULL) {
    mriglm->beta = MRIallocSequence(nc, nr, ns, MRI_FLOAT, mriglm->nregtot);
    MRIcopyHeader(mriglm->y, mriglm->beta);
    mriglm->eres = MRIallocSequence(nc, nr, ns, MRI_FLOAT, nf);
    MRIcopyHeader(mriglm->y, mriglm->eres);
    mriglm->rvar = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
    MRIcopyHeader(mriglm->y, mriglm->rvar);
    if (mriglm->yhatsave) {
      mriglm->yhat = MRIallocSequence(nc, nr, ns, MRI_FLOAT, nf);
      MRIcopyHeader(mriglm->y, mriglm->yhat);
    }
    if (mriglm->condsave) {
      mriglm->cond = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
      MRIcopyHeader(mriglm->y, mriglm->cond);
    }

    for (n = 0; n < glm->ncontrasts; n++) {
      mriglm->gamma[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, glm->C[n]->rows);
      MRIcopyHeader(mriglm->y, mriglm->gamma[n]);
      if (glm->C[n]->rows == 1) {
        mriglm->gammaVar[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
        MRIcopyHeader(mriglm->y, mriglm->gammaVar[n]);
        if (glm->DoPCC) {
          mriglm->pcc[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
          MRIcopyHeader(mriglm->y, mriglm->pcc[n]);
        }
      }
      mriglm->F[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
      MRIcopyHeader(mriglm->y, mriglm->F[n]);
      mriglm->p[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
      MRIcopyHeader(mriglm->y, mriglm->p[n]);
      mriglm->z[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
      MRIcopyHeader(mriglm->y, mriglm->z[n]);
      if (glm->ypmfflag[n]) {
        mriglm->ypmf[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, nf);
        MRIcopyHeader(mriglm->y, mriglm->ypmf[n]);
      }
    }
  }

  //--------------------------------------------
  //pctdone = 0;
  //nthvox = 0;
  mriglm->n_ill_cond = 0;
  long n_ill_cond = 0;

  // Parallel does not work yet because need separate glm for each thread
  //#ifdef HAVE_OPENMP
  //#pragma omp parallel for if_ROMP(assume_reproducible) reduction(+ : n_ill_cond)
  //#endif
  for (c = 0; c < nc; c++) {
    int r,s,nthvox=0,m,n,pctdone=0;
    double Xcond;
    for (r = 0; r < nr; r++) {
      for (s = 0; s < ns; s++) {
        nthvox++;
        if (nthvox == (long)floor(.1 * nvoxtot)) {
          pctdone += 10;
          if (Gdiag_no > 0) {
            printf("%2d%% ", pctdone);
            fflush(stdout);
          }
          nthvox = 0;
        }

        // Check the mask -----------
        if (mriglm->mask != NULL) {
          m = MRIgetVoxVal(mriglm->mask, c, r, s, 0);
          if (m < 0.5) continue;
        }

        // Get data from mri and put in GLM
        MRIglmLoadVox(mriglm, c, r, s, 0, NULL);

        // Compute intermediate matrices
        GLMxMatrices(glm);

        // Compute condition
        if (mriglm->condsave) {
          Xcond = MatrixConditionNumber(glm->XtX);
          MRIsetVoxVal(mriglm->cond, c, r, s, 0, Xcond);
        }

        // Test condition
        if (glm->ill_cond_flag) {
          n_ill_cond++;
          continue;
        }

        GLMfit(glm);
        if (mriglm->yffxvar == NULL)
          GLMtest(glm);
        else
          GLMtestFFx(glm);

        // Pack data back into MRI
        MRIsetVoxVal(mriglm->rvar, c, r, s, 0, glm->rvar);
        MRIfromMatrix(mriglm->beta, c, r, s, glm->beta, NULL);
        MRIfromMatrix(mriglm->eres, c, r, s, glm->eres, mriglm->FrameMask);
        if (mriglm->yhatsave) MRIfromMatrix(mriglm->yhat, c, r, s, glm->yhat, mriglm->FrameMask);
        for (n = 0; n < glm->ncontrasts; n++) {
          MRIfromMatrix(mriglm->gamma[n], c, r, s, glm->gamma[n], NULL);
          if (glm->C[n]->rows == 1)
            MRIsetVoxVal(mriglm->gammaVar[n], c, r, s, 0, glm->gCVM[n]->rptr[1][1]);
          MRIsetVoxVal(mriglm->F[n], c, r, s, 0, glm->F[n]);
          MRIsetVoxVal(mriglm->p[n], c, r, s, 0, glm->p[n]);
          MRIsetVoxVal(mriglm->z[n], c, r, s, 0, glm->z[n]);
          if (glm->C[n]->rows == 1 && glm->DoPCC)
            MRIsetVoxVal(mriglm->pcc[n], c, r, s, 0, glm->pcc[n]);

          if (glm->ypmfflag[n])
            MRIfromMatrix(mriglm->ypmf[n], c, r, s, glm->ypmf[n], mriglm->FrameMask);
        }
      }
    }
  }
  if (Gdiag_no > 0) printf("\n");
  mriglm->n_ill_cond = n_ill_cond;
  // printf("n_ill_cond = %d\n",mriglm->n_ill_cond);
  return (0);
}

/*---------------------------------------------------------------------
  MRIglmFit() - fits glm (beta and rvar) on a voxel-by-voxel basis.
  Made to be followed by MRIglmTest(). See notes on MRIglmFitandTest()
  --------------------------------------------------------------------*/
int MRIglmFit(MRIGLM *mriglm)
{
  int c, r, s, nc, nr, ns, nf, pctdone;
  float m, Xcond;
  long nvoxtot, nthvox;

  nc = mriglm->y->width;
  nr = mriglm->y->height;
  ns = mriglm->y->depth;
  nf = mriglm->y->nframes;
  nvoxtot = nc * nr * ns;

  mriglm->nregtot = MRIglmNRegTot(mriglm);
  GLMallocX(mriglm->glm, nf, mriglm->nregtot);
  GLMallocY(mriglm->glm);
  if (mriglm->yffxvar) GLMallocYFFxVar(mriglm->glm);

  if (mriglm->w != NULL || mriglm->npvr != 0)
    mriglm->pervoxflag = 1;
  else
    mriglm->pervoxflag = 0;

  GLMcMatrices(mriglm->glm);

  if (!mriglm->pervoxflag) {
    MatrixCopy(mriglm->Xg, mriglm->glm->X);
    mriglm->XgLoaded = 1;
    GLMxMatrices(mriglm->glm);
  }

  // If beta has not been allocated, assume that no one has been alloced
  if (mriglm->beta == NULL) {
    mriglm->beta = MRIallocSequence(nc, nr, ns, MRI_FLOAT, mriglm->nregtot);
    MRIcopyHeader(mriglm->y, mriglm->beta);
    mriglm->eres = MRIallocSequence(nc, nr, ns, MRI_FLOAT, nf);
    MRIcopyHeader(mriglm->y, mriglm->eres);
    mriglm->rvar = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
    MRIcopyHeader(mriglm->y, mriglm->rvar);
    if (mriglm->yhatsave) {
      mriglm->yhat = MRIallocSequence(nc, nr, ns, MRI_FLOAT, nf);
      MRIcopyHeader(mriglm->y, mriglm->yhat);
    }
    if (mriglm->condsave) {
      mriglm->cond = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
      MRIcopyHeader(mriglm->y, mriglm->cond);
    }
  }

  //--------------------------------------------
  pctdone = 0;
  nthvox = 0;
  mriglm->n_ill_cond = 0;
  for (c = 0; c < nc; c++) {
    for (r = 0; r < nr; r++) {
      for (s = 0; s < ns; s++) {
        nthvox++;
        if (nthvox == (long)floor(.1 * nvoxtot)) {
          pctdone += 10;
          printf("%2d%% ", pctdone);
          fflush(stdout);
          nthvox = 0;
        }

        // Check the mask -----------
        if (mriglm->mask != NULL) {
          m = MRIgetVoxVal(mriglm->mask, c, r, s, 0);
          if (m < 0.5) continue;
        }

        // Get data from mri and put in GLM
        MRIglmLoadVox(mriglm, c, r, s, 0, NULL);

        // Compute intermediate matrices
        GLMxMatrices(mriglm->glm);

        // Compute condition
        if (mriglm->condsave) {
          Xcond = MatrixConditionNumber(mriglm->glm->XtX);
          MRIsetVoxVal(mriglm->cond, c, r, s, 0, Xcond);
        }

        // Test condition
        if (mriglm->glm->ill_cond_flag) {
          mriglm->n_ill_cond++;
          continue;
        }

        GLMfit(mriglm->glm);

        // Pack data back into MRI
        MRIsetVoxVal(mriglm->rvar, c, r, s, 0, mriglm->glm->rvar);
        MRIfromMatrix(mriglm->beta, c, r, s, mriglm->glm->beta, NULL);
        MRIfromMatrix(mriglm->eres, c, r, s, mriglm->glm->eres, mriglm->FrameMask);
        if (mriglm->yhatsave) MRIfromMatrix(mriglm->yhat, c, r, s, mriglm->glm->yhat, mriglm->FrameMask);
      }
    }
  }
  printf("\n");

  // printf("n_ill_cond = %d\n",mriglm->n_ill_cond);
  return (0);
}

/*---------------------------------------------------------------------
  MRIglmTest() - tests glm contrasts on a voxel-by-voxel basis.
  Made to be preceded by MRIglmFit(). See notes on MRIglmFitandTest()
  --------------------------------------------------------------------*/
int MRIglmTest(MRIGLM *mriglm)
{
  int c, r, s, n, nc, nr, ns, nf, pctdone;
  float m;
  long nvoxtot, nthvox;

  if (mriglm->glm->ncontrasts == 0) return (0);

  nc = mriglm->y->width;
  nr = mriglm->y->height;
  ns = mriglm->y->depth;
  nf = mriglm->y->nframes;
  nvoxtot = nc * nr * ns;

  // If gamma[0] not been allocated, assume that no one has been alloced
  if (mriglm->gamma[0] == NULL) {
    for (n = 0; n < mriglm->glm->ncontrasts; n++) {
      mriglm->gamma[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, mriglm->glm->C[n]->rows);
      MRIcopyHeader(mriglm->y, mriglm->gamma[n]);
      if (mriglm->glm->C[n]->rows == 1) {
        mriglm->gammaVar[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
        MRIcopyHeader(mriglm->y, mriglm->gammaVar[n]);
        if (mriglm->glm->DoPCC) {
          mriglm->pcc[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
          MRIcopyHeader(mriglm->y, mriglm->pcc[n]);
        }
      }
      mriglm->F[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
      MRIcopyHeader(mriglm->y, mriglm->F[n]);
      mriglm->p[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
      MRIcopyHeader(mriglm->y, mriglm->p[n]);
      mriglm->z[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, 1);
      MRIcopyHeader(mriglm->y, mriglm->z[n]);
      if (mriglm->glm->ypmfflag[n]) {
        mriglm->ypmf[n] = MRIallocSequence(nc, nr, ns, MRI_FLOAT, nf);
        MRIcopyHeader(mriglm->y, mriglm->ypmf[n]);
      }
    }
  }

  //--------------------------------------------
  pctdone = 0;
  nthvox = 0;
  for (c = 0; c < nc; c++) {
    for (r = 0; r < nr; r++) {
      for (s = 0; s < ns; s++) {
        nthvox++;
        if (nthvox == (long)floor(.1 * nvoxtot)) {
          pctdone += 10;
          printf("%2d%% ", pctdone);
          fflush(stdout);
          nthvox = 0;
        }

        // Check the mask -----------
        if (mriglm->mask != NULL) {
          m = MRIgetVoxVal(mriglm->mask, c, r, s, 0);
          if (m < 0.5) continue;
        }

        // Get data from mri and put in GLM
        MRIglmLoadVox(mriglm, c, r, s, 1, NULL);

        // Compute intermediate matrices
        GLMxMatrices(mriglm->glm);

        // Test
        if (mriglm->yffxvar == NULL)
          GLMtest(mriglm->glm);
        else
          GLMtestFFx(mriglm->glm);

        // Pack data back into MRI
        for (n = 0; n < mriglm->glm->ncontrasts; n++) {
          MRIfromMatrix(mriglm->gamma[n], c, r, s, mriglm->glm->gamma[n], NULL);
          if (mriglm->glm->C[n]->rows == 1) {
            MRIsetVoxVal(mriglm->gammaVar[n], c, r, s, 0, mriglm->glm->gCVM[n]->rptr[1][1]);
            if (mriglm->glm->DoPCC) MRIsetVoxVal(mriglm->pcc[n], c, r, s, 0, mriglm->glm->pcc[n]);
          }
          MRIsetVoxVal(mriglm->F[n], c, r, s, 0, mriglm->glm->F[n]);
          MRIsetVoxVal(mriglm->p[n], c, r, s, 0, mriglm->glm->p[n]);
          MRIsetVoxVal(mriglm->z[n], c, r, s, 0, mriglm->glm->z[n]);
          if (mriglm->glm->ypmfflag[n])
            MRIfromMatrix(mriglm->ypmf[n], c, r, s, mriglm->glm->ypmf[n], mriglm->FrameMask);
        }
      }
    }
  }
  printf("\n");

  // printf("n_ill_cond = %d\n",mriglm->n_ill_cond);
  return (0);
}

/* ---------------------------------------------------------------------------
   MRIglmLoadVox() - loads the data (X and y) for the voxel at crs
   into the GLM matrix struct and applies the weights (if applicable).
   X and y are allocated if they are NULL. mriglm->nregtot is also computed
   here. If Xg has already been loaded into X, then it is not loaded again
   unless mriglm->w is non-null.
   -------------------------------------------------------------------------*/
int MRIglmLoadVox(MRIGLM *mriglm, int c, int r, int s, int LoadBeta, GLMMAT *glm)
{
  int f, n, nthreg, nthf, nf;
  double v;
  static int nfprev = -1;
  if(glm == NULL) glm = mriglm->glm;

  nf = mriglm->y->nframes;
  // Count the number of frames in frame mask
  if (mriglm->FrameMask != NULL) {
    nf = 0;
    for (f = 1; f <= mriglm->y->nframes; f++)
      if (MRIgetVoxVal(mriglm->FrameMask, c, r, s, f - 1) > 0.5) nf++;
    if (nf == 0) printf("MRIglmLoadVox(): %d,%d,%d nf=0\n", c, r, s);
    // Free matrices if needed
    if (glm->X != NULL && nfprev != nf) MatrixFree(&(glm->X));
    if (glm->y != NULL && nfprev != nf) MatrixFree(&(glm->y));
    nfprev = nf;
  }

  // Alloc matrices if needed
  if (glm->X == NULL) {
    mriglm->nregtot = mriglm->Xg->cols + mriglm->npvr;
    glm->X = MatrixAlloc(nf, mriglm->nregtot, MATRIX_REAL);
  }
  if (glm->y == NULL) glm->y = MatrixAlloc(nf, 1, MATRIX_REAL);

  // Load y, Xg, and the per-vox reg --------------------------
  nthf = 0;
  for (f = 1; f <= mriglm->y->nframes; f++) {
    if (mriglm->FrameMask != NULL && MRIgetVoxVal(mriglm->FrameMask, c, r, s, f - 1) < 0.5) continue;
    nthf++;

    // Load y
    glm->y->rptr[nthf][1] = MRIgetVoxVal(mriglm->y, c, r, s, f - 1);

    // Load Xg->X the global design matrix if needed
    // For wg, this is a little bit of a hack. wg needs to be applied to Xg only once,
    // but it will get applied again and again. Including wg here forces Xg to be
    // freshly copied into X each time, then wg is applied.
    if (mriglm->w != NULL || mriglm->wg != NULL || !mriglm->XgLoaded || mriglm->FrameMask) {
      nthreg = 1;
      for (n = 1; n <= mriglm->Xg->cols; n++) {
        glm->X->rptr[nthf][nthreg] = mriglm->Xg->rptr[f][n];  // X=Xg
        nthreg++;
      }
    }
    else
      nthreg = mriglm->Xg->cols + 1;

    // Load the global per-voxel regressors matrix, X = [X pvr]
    for (n = 1; n <= mriglm->npvr; n++) {
      glm->X->rptr[nthf][nthreg] = MRIgetVoxVal(mriglm->pvr[n - 1], c, r, s, f - 1);
      nthreg++;
    }
  }
  mriglm->XgLoaded = 1;  // Set flag that Xg has been loaded

  // Weight X and y, X = w.*X, y = w.*y
  if ((mriglm->w != NULL || mriglm->wg != NULL) && !mriglm->skipweight) {
    nthf = 0;
    for (f = 1; f <= mriglm->y->nframes; f++) {
      if (mriglm->FrameMask != NULL && MRIgetVoxVal(mriglm->FrameMask, c, r, s, f - 1) < 0.5) continue;
      nthf++;
      if (mriglm->w != NULL)
        v = MRIgetVoxVal(mriglm->w, c, r, s, f - 1);
      else
        v = mriglm->wg->rptr[f][1];
      glm->y->rptr[nthf][1] *= v;
      for (n = 1; n <= glm->X->cols; n++) glm->X->rptr[nthf][n] *= v;
    }
  }

  // Load ffx variance, if there
  if (mriglm->yffxvar != NULL) {
    nthf = 0;
    for (f = 1; f <= mriglm->y->nframes; f++) {
      if (mriglm->FrameMask != NULL && MRIgetVoxVal(mriglm->FrameMask, c, r, s, f - 1) < 0.5) continue;
      nthf++;
      v = MRIgetVoxVal(mriglm->yffxvar, c, r, s, f - 1);
      glm->yffxvar->rptr[nthf][1] = v;
    }
    glm->ffxdof = mriglm->ffxdof;
  }

  // Beta
  if (LoadBeta) {
    for (f = 1; f <= glm->X->cols; f++) {
      v = MRIgetVoxVal(mriglm->beta, c, r, s, f - 1);
      glm->beta->rptr[f][1] = v;
    }
    v = MRIgetVoxVal(mriglm->rvar, c, r, s, 0);
    glm->rvar = v;
  }

  return (0);
}
/*----------------------------------------------------------------
  MRIglmNRegTot() - computes the total number of regressors based
  on the number of columns in Xg + number of per-voxel regressors
  ----------------------------------------------------------------*/
int MRIglmNRegTot(MRIGLM *mriglm)
{
  mriglm->nregtot = mriglm->Xg->cols + mriglm->npvr;
  return (mriglm->nregtot);
}

/*----------------------------------------------------------------
  MRItoVector() - copies all the frames from the given voxel
  in to a vector.
  ----------------------------------------------------------------*/
VECTOR *MRItoVector(MRI *mri, int c, int r, int s, VECTOR *v)
{
  int f;
  if (v == NULL) v = MatrixAlloc(mri->nframes, 1, MATRIX_REAL);

  for (f = 1; f <= v->rows; f++) v->rptr[f][1] = MRIgetVoxVal(mri, c, r, s, f - 1);
  return (v);
}

/*---------------------------------------------------------------
  MRIsetSign() - sets the sign of the invol based on the sign of the
  nth frame of the signvol. The values of the input are changed.  All
  frames of the input volume are affected.
  --------------------------------------------------------------*/
int MRIsetSign(MRI *invol, MRI *signvol, int frame)
{
  int c, r, s, f;
  double v, sgn;

  if (frame > signvol->nframes) {
    printf("ERROR: MRIsetSign(): input frame %d is too large", frame);
    return (1);
  }

  for (c = 0; c < invol->width; c++) {
    for (r = 0; r < invol->height; r++) {
      for (s = 0; s < invol->depth; s++) {
        sgn = MRIgetVoxVal(signvol, c, r, s, frame);
        for (f = 0; f < invol->nframes; f++) {
          v = MRIgetVoxVal(invol, c, r, s, f);
          if (sgn < 0.0) v = -1.0 * fabs(v);
          if (sgn > 0.0) v = +1.0 * fabs(v);
          MRIsetVoxVal(invol, c, r, s, f, v);
        }
      }
    }
  }
  return (0);
}

/*----------------------------------------------------------------
  MRI *MRIvolMax(MRI *vol, MRI *out) - the value at each voxel
  is the maximum over the frames at that voxel.
  --------------------------------------------------------------*/
MRI *MRIvolMax(MRI *invol, MRI *out)
{
  int c, r, s, f;
  double v, max;

  if (out == NULL) {
    out = MRIalloc(invol->width, invol->height, invol->depth, invol->type);
    if (out == NULL) return (NULL);
    MRIcopyHeader(invol, out);
  }
  if (out->width != invol->width || out->height != invol->height || out->depth != invol->depth) {
    printf("ERROR: MRIvolMax: dimension mismatch\n");
    return (NULL);
  }

  for (c = 0; c < invol->width; c++) {
    for (r = 0; r < invol->height; r++) {
      for (s = 0; s < invol->depth; s++) {
        max = MRIgetVoxVal(invol, c, r, s, 0);
        for (f = 1; f < invol->nframes; f++) {
          v = MRIgetVoxVal(invol, c, r, s, f);
          if (max < v) max = v;
        }
        MRIsetVoxVal(out, c, r, s, 0, max);
      }
    }
  }
  return (out);
}

/*----------------------------------------------------------------
  MRI *MRIvolMin(MRI *vol, MRI *out) - the value at each voxel
  is the minimum over the frames at that voxel. Note that nothing
  special is done with the sign (ie, a negative value can be the
  minimum).
  --------------------------------------------------------------*/
MRI *MRIvolMin(MRI *invol, MRI *out)
{
  int c, r, s, f;
  double v, min;

  if (out == NULL) {
    out = MRIalloc(invol->width, invol->height, invol->depth, invol->type);
    if (out == NULL) return (NULL);
    MRIcopyHeader(invol, out);
  }
  if (out->width != invol->width || out->height != invol->height || out->depth != invol->depth) {
    printf("ERROR: MRIvolMin: dimension mismatch\n");
    return (NULL);
  }

  for (c = 0; c < invol->width; c++) {
    for (r = 0; r < invol->height; r++) {
      for (s = 0; s < invol->depth; s++) {
        min = MRIgetVoxVal(invol, c, r, s, 0);
        for (f = 1; f < invol->nframes; f++) {
          v = MRIgetVoxVal(invol, c, r, s, f);
          if (min > v) min = v;
        }
        MRIsetVoxVal(out, c, r, s, 0, min);
      }
    }
  }
  return (out);
}

/*!
\fn MRI *MRIconjunct(MRI *invol, MRI *out)
\brief Performs "conjunction" by taking the min of the abs across frames
at each voxel. The value at the voxel is the min, including the true
sign of the min. Eg, if the two frames are:
   +2.1 and +3.4 --> +2.1
   -2.1 and -3.4 --> -2.1
   +2.1 and -3.4 --> +2.1
   -2.1 and +3.4 --> -2.1
See: Thomas Nichols, Matthew Brett, Jesper Andersson, Tor Wager &
Jean-Baptiste Poline.  NeuroImage, Volume 25, Issue 3, 15 April 2005,
Pages 653-660.
\param multi-frame input vol
\param out - can be NULL
*/
MRI *MRIconjunct(MRI *invol, MRI *out)
{
  int c, r, s, f;
  double v, min, minsign;

  if (out == NULL) {
    out = MRIalloc(invol->width, invol->height, invol->depth, MRI_FLOAT);
    if (out == NULL) return (NULL);
    MRIcopyHeader(invol, out);
  }
  if (out->width != invol->width || out->height != invol->height || out->depth != invol->depth) {
    printf("ERROR: MRIconjunct: dimension mismatch\n");
    return (NULL);
  }

  for (c = 0; c < invol->width; c++) {
    for (r = 0; r < invol->height; r++) {
      for (s = 0; s < invol->depth; s++) {
        v = MRIgetVoxVal(invol, c, r, s, 0);
        min = fabs(v);
        minsign = SIGN(v);
        for (f = 0; f < invol->nframes; f++) {
          v = MRIgetVoxVal(invol, c, r, s, f);
          if (min > fabs(v)) {
            min = fabs(v);
            minsign = SIGN(v);
          }
        }
        MRIsetVoxVal(out, c, r, s, 0, minsign * min);
      }
    }
  }
  return (out);
}

/*----------------------------------------------------------------
  MRI *MRIvolMaxIndex(MRI *vol, MRI *out) - the value at each voxel
  is the frame index of the maximum over the frames at that voxel.
  base is added to the index.
  --------------------------------------------------------------*/
MRI *MRIvolMaxIndex(MRI *invol, int base, MRI *mask, MRI *out)
{
  int c, r, s, f, index;
  double v, max, m;

  if (out == NULL) {
    out = MRIalloc(invol->width, invol->height, invol->depth, MRI_INT);
    if (out == NULL) return (NULL);
    MRIcopyHeader(invol, out);
  }
  if (out->width != invol->width || out->height != invol->height || out->depth != invol->depth) {
    printf("ERROR: MRIvolMax: dimension mismatch\n");
    return (NULL);
  }

  for (c = 0; c < invol->width; c++) {
    for (r = 0; r < invol->height; r++) {
      for (s = 0; s < invol->depth; s++) {
        MRIsetVoxVal(out, c, r, s, 0, 0);
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        max = MRIgetVoxVal(invol, c, r, s, 0);
        index = 0;
        for (f = 1; f < invol->nframes; f++) {
          v = MRIgetVoxVal(invol, c, r, s, f);
          if (max < v) {
            max = v;
            index = f;
          }
        }
        MRIsetVoxVal(out, c, r, s, 0, index + base);
      }
    }
  }
  return (out);
}

/*---------------------------------------------------------------
  MRIframeMax() - finds the maximum in the given frame. The max is
  returned. The CRS of the max are passed back as args. If mask is
  non-null, then only voxels in the mask are considered. If signflag
  is 0, then the maximum of the absolute is searched for (signed
  value still returned). If signflag = +1, then only the pos max
  is found. If signflag = -1, then only the neg max  is found.
  --------------------------------------------------------------*/
double MRIframeMax(MRI *vol, int frame, MRI *mask, int signflag, int *cmax, int *rmax, int *smax)
{
  int c, r, s, nhits;
  double v, vmax, m;

  if (frame > vol->nframes) {
    printf("ERROR: MRIframeMax(): input frame %d is too large", frame);
    return (1);
  }

  nhits = -1;
  vmax = 0.0;
  for (c = 0; c < vol->width; c++) {
    for (r = 0; r < vol->height; r++) {
      for (s = 0; s < vol->depth; s++) {
        if (mask != NULL) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        nhits++;
        v = MRIgetVoxVal(vol, c, r, s, frame);

        if (nhits == 0) {  // first hit
          vmax = v;
          *cmax = c;
          *rmax = r;
          *smax = s;
          continue;
        }

        switch (signflag) {
          case 0:  // absolute
            if (fabs(vmax) < fabs(v)) {
              vmax = v;
              *cmax = c;
              *rmax = r;
              *smax = s;
            }
            break;
          case 1:  // positive
            if (vmax < v) {
              vmax = v;
              *cmax = c;
              *rmax = r;
              *smax = s;
            }
            break;
          case -1:  // negative
            if (vmax > v) {
              vmax = v;
              *cmax = c;
              *rmax = r;
              *smax = s;
            }
            break;
        }  // end swtich
      }    // s
    }      // r
  }        // s
  return (vmax);
}

/*---------------------------------------------------------------
  MRIframeMean() - computes mean over frames of each voxel.
  --------------------------------------------------------------*/
MRI *MRIframeMean(MRI *vol, MRI *volmn)
{
  int c, r, s, f;
  double v;

  if (volmn == NULL) {
    volmn = MRIallocSequence(vol->width, vol->height, vol->depth, MRI_FLOAT, 1);
    MRIcopyHeader(vol, volmn);
  }

  for (c = 0; c < vol->width; c++) {
    for (r = 0; r < vol->height; r++) {
      for (s = 0; s < vol->depth; s++) {
        v = 0;
        for (f = 0; f < vol->nframes; f++) v += MRIgetVoxVal(vol, c, r, s, f);
        MRIsetVoxVal(volmn, c, r, s, 0, v / vol->nframes);
      }  // s
    }    // r
  }      // s
  return (volmn);
}

/*---------------------------------------------------------------
  MRIframeMedian() - computes median over frames of each voxel.
  --------------------------------------------------------------*/
MRI *MRIframeMedian(MRI *vol, MRI *volmn)
{
  int c, r, s, f;
  float *t = (float *)calloc(vol->nframes, sizeof(float));

  if (volmn == NULL) {
    volmn = MRIallocSequence(vol->width, vol->height, vol->depth, MRI_FLOAT, 1);
    MRIcopyHeader(vol, volmn);
  }

  for (c = 0; c < vol->width; c++) {
    for (r = 0; r < vol->height; r++) {
      for (s = 0; s < vol->depth; s++) {
        for (f = 0; f < vol->nframes; f++) t[f] = MRIgetVoxVal(vol, c, r, s, f);
        MRIsetVoxVal(volmn, c, r, s, 0, median(t, vol->nframes));
      }  // s
    }    // r
  }      // c
  free(t);
  return (volmn);
}

/*---------------------------------------------------------------
  MRIframeSum() - computes sum over frames of each voxel.
  --------------------------------------------------------------*/
MRI *MRIframeSum(MRI *vol, MRI *volsum)
{
  int c, r, s, f;
  double v;

  if (volsum == NULL) {
    volsum = MRIallocSequence(vol->width, vol->height, vol->depth, MRI_FLOAT, 1);
    MRIcopyHeader(vol, volsum);
  }

  for (c = 0; c < vol->width; c++) {
    for (r = 0; r < vol->height; r++) {
      for (s = 0; s < vol->depth; s++) {
        v = 0;
        for (f = 0; f < vol->nframes; f++) v += MRIgetVoxVal(vol, c, r, s, f);
        MRIsetVoxVal(volsum, c, r, s, 0, v);
      }  // s
    }    // r
  }      // s
  return (volsum);
}

/*---------------------------------------------------------------
  fMRIdetrend() - returns (I-inv(X'*X)*X')*y
  ---------------------------------------------------------------*/
MRI *fMRIdetrend(MRI *y, MATRIX *X)
{
  MATRIX *Xt, *XtX, *iXtX, *B;
  MRI *beta, *yhat, *res;

  if (X->rows != y->nframes) {
    printf("ERROR: dimension mismatch between X and input\n");
    return (NULL);
  }

  Xt = MatrixTranspose(X, NULL);
  XtX = MatrixMultiplyD(Xt, X, NULL);
  iXtX = MatrixInverse(XtX, NULL);
  if (iXtX == NULL) {
    printf("ERROR: could not compute psuedo inverse of X\n");
    exit(1);
  }
  B = MatrixMultiplyD(iXtX, Xt, NULL);

  beta = fMRImatrixMultiply(y, B, NULL);
  yhat = fMRImatrixMultiply(beta, X, NULL);
  res = MRIsubtract(y, yhat, NULL);

  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&B);
  MRIfree(&beta);
  MRIfree(&yhat);

  return (res);
}

/*---------------------------------------------------------------------
  fMRIspatialAR1() - computes spatial AR1, ie, the correlation between
  the time course at one voxel and that at an adjacent voxel. There
  will be six frame, 2 for each dim. If a mask is used, then a voxel
  and all it's neighbors must be in the mask in order to compute the
  AR1 at that point. This function assumess that the mean and any
  other trends have been removed. It works regardless of the DOF of
  the time series.
  --------------------------------------------------------------------*/
MRI *fMRIspatialAR1(MRI *src, MRI *mask, MRI *ar1)
{
  int c, r, s, f, dc, dr, ds, skip, nhits;
  MRI *srcsumsq, *srctmp;
  double m, c1sum, c2sum, r1sum, r2sum, s1sum, s2sum;
  double v0, vc1, vc2, vr1, vr2, vs1, vs2, sumsq0;
  double car1, rar1, sar1;
  int freetmp;

  freetmp = 0;
  if (src->type != MRI_FLOAT) {
    srctmp = MRISeqchangeType(src, MRI_FLOAT, 0, 0, 0);
    freetmp = 1;
  }
  else {
    srctmp = src;
    freetmp = 0;
  }

  // alloc vol with 6 frames
  if (ar1 == NULL) {
    ar1 = MRIcloneBySpace(src, MRI_FLOAT, 6);
    if (ar1 == NULL) {
      printf("ERROR: could not alloc\n");
      return (NULL);
    }
  }

  // pre-compute the sum of squares
  srcsumsq = fMRIsumSquare(srctmp, 0, NULL);

  // Loop thru all voxels
  nhits = 0;
  for (c = 0; c < srctmp->width; c++) {
    for (r = 0; r < srctmp->height; r++) {
      for (s = 0; s < srctmp->depth; s++) {
        // skip voxel if it's on the edge
        if (c == 0 || r == 0 || s == 0 || c == (srctmp->width - 1) || r == (srctmp->height - 1) ||
            s == (srctmp->depth - 1)) {
          for (f = 0; f < 6; f++) MRIsetVoxVal(ar1, c, r, s, f, 0);
          continue;
        }
        // sum-of-sqares at center voxel
        sumsq0 = MRIgetVoxVal(srcsumsq, c, r, s, 0);
        if (sumsq0 < 1e-6) continue;
        // skip if BOTH voxel and all it's neighbors are
        // not in the mask
        if (mask) {
          skip = 0;
          for (dc = -1; dc < 2; dc++) {
            for (dr = -1; dr < 2; dr++) {
              for (ds = -1; ds < 2; ds++) {
                m = MRIgetVoxVal(mask, c + dc, r + dr, s + ds, 0);
                if (m < 0.5) {
                  for (f = 0; f < 6; f++) MRIsetVoxVal(ar1, c, r, s, f, 0);
                  skip = 1;
                }
              }
            }
          }
          if (skip) continue;
        }
        nhits++;

        // Loop thru all frames
        c1sum = 0;
        c2sum = 0;
        r1sum = 0;
        r2sum = 0;
        s1sum = 0;
        s2sum = 0;
        for (f = 0; f < srctmp->nframes; f++) {
          v0 = MRIgetVoxVal(srctmp, c, r, s, f);  // value at center voxel

          // temporal correlation with vox one col to left
          vc1 = MRIgetVoxVal(srctmp, c - 1, r, s, f);
          c1sum += (v0 * vc1);

          // temporal correlation with vox one col to right
          vc2 = MRIgetVoxVal(srctmp, c + 1, r, s, f);
          c2sum += (v0 * vc2);

          // temporal correlation with vox one row to up
          vr1 = MRIgetVoxVal(srctmp, c, r - 1, s, f);
          r1sum += (v0 * vr1);

          // temporal correlation with vox one row to down
          vr2 = MRIgetVoxVal(srctmp, c, r + 1, s, f);
          r2sum += (v0 * vr2);

          // temporal correlation with vox one slice in
          vs1 = MRIgetVoxVal(srctmp, c, r, s - 1, f);
          s1sum += (v0 * vs1);

          // temporal correlation with vox one slice out
          vs2 = MRIgetVoxVal(srctmp, c, r, s + 1, f);
          s2sum += (v0 * vs2);
        }

        // sum-of-squares at center voxel
        sumsq0 = MRIgetVoxVal(srcsumsq, c, r, s, 0);

        // column AR1
        vc1 = MRIgetVoxVal(srcsumsq, c - 1, r, s, 0);  // variance
        if (vc1 > 1e-6)
          car1 = c1sum / (sqrt(sumsq0 * vc1));
        else
          car1 = 0;
        MRIsetVoxVal(ar1, c, r, s, 0, car1);  // frame 0

        vc2 = MRIgetVoxVal(srcsumsq, c + 1, r, s, 0);
        if (vc2 > 1e-6)
          car1 = c2sum / (sqrt(sumsq0 * vc2));
        else
          car1 = 0;
        MRIsetVoxVal(ar1, c, r, s, 1, car1);  // frame 1

        // rows
        vr1 = MRIgetVoxVal(srcsumsq, c, r - 1, s, 0);
        if (vr1 > 1e-6)
          rar1 = r1sum / (sqrt(sumsq0 * vr1));
        else
          rar1 = 0;
        MRIsetVoxVal(ar1, c, r, s, 2, rar1);  // frame 2

        vr2 = MRIgetVoxVal(srcsumsq, c, r + 1, s, 0);
        if (vr2 > 1e-6)
          rar1 = r2sum / (sqrt(sumsq0 * vr2));
        else
          rar1 = 0;
        MRIsetVoxVal(ar1, c, r, s, 3, rar1);  // frame 3

        // slices
        vs1 = MRIgetVoxVal(srcsumsq, c, r, s - 1, 0);
        if (vs1 > 1e-6)
          sar1 = s1sum / (sqrt(sumsq0 * vs1));
        else
          sar1 = 0;
        MRIsetVoxVal(ar1, c, r, s, 4, sar1);  // frame 4

        vs2 = MRIgetVoxVal(srcsumsq, c, r, s + 1, 0);
        if (vs2 > 1e-6)
          sar1 = s2sum / (sqrt(sumsq0 * vs2));
        else
          sar1 = 0;
        MRIsetVoxVal(ar1, c, r, s, 5, sar1);  // frame 5

      }  // s
    }    // r
  }      // c

  printf("fMRIspatialAR1(): hit %d voxels\n", nhits);
  if (nhits == 0) printf("WARNING: no voxels in AR1 computation\n");

  MRIfree(&srcsumsq);
  if (freetmp) MRIfree(&srctmp);
  return (ar1);
}

/*---------------------------------------------------------------
  fMRIspatialAR2() - computes spatial AR2, ie, the correlation between
  the time course at one voxel and that at a voxel two voxels
  over. There will be six frame, 2 for each dim. If a mask is used,
  then a voxel and all it's neighbors must be in the mask in order to
  compute the AR2 at that point. BUG: this function is sensitive
  to the time series DOF (assumes it is nframes-1).
  --------------------------------------------------------------*/
MRI *fMRIspatialAR2(MRI *src, MRI *mask, MRI *ar2)
{
  int c, r, s, f, nframes, dc, dr, ds, skip;
  MRI *srcvar, *srctmp;
  double m, c1sum, c2sum, r1sum, r2sum, s1sum, s2sum;
  double v0, vc1, vc2, vr1, vr2, vs1, vs2;
  double car2, rar2, sar2;
  int freetmp;

  freetmp = 0;
  if (src->type != MRI_FLOAT) {
    srctmp = MRISeqchangeType(src, MRI_FLOAT, 0, 0, 0);
    freetmp = 1;
  }
  else {
    srctmp = src;
    freetmp = 0;
  }

  // alloc vol with 3 frames
  if (ar2 == NULL) {
    ar2 = MRIcloneBySpace(src, MRI_FLOAT, 6);
    if (ar2 == NULL) {
      printf("ERROR: could not alloc\n");
      return (NULL);
    }
  }

  // pre-compute the variance
  // srcvar = fMRIvariance(srctmp, -1, 0, NULL);
  srcvar = fMRIcovariance(srctmp, 0, -1, mask, NULL);

  nframes = srctmp->nframes - 1;  // assume DOF = nframes - 1

  // Loop thru all voxels
  for (c = 0; c < srctmp->width; c++) {
    for (r = 0; r < srctmp->height; r++) {
      for (s = 0; s < srctmp->depth; s++) {
        // skip voxel if it's on the edge
        if (c < 2 || r < 2 || s < 2 || c >= (srctmp->width - 2) || r >= (srctmp->height - 2) ||
            s >= (srctmp->depth - 2)) {
          for (f = 0; f < 6; f++) MRIsetVoxVal(ar2, c, r, s, f, 0);
          continue;
        }
        // variance at center voxel
        v0 = MRIgetVoxVal(srcvar, c, r, s, 0);
        if (v0 < 1e-6) continue;
        // skip if BOTH voxel and all it's neighbors are
        // not in the mask
        if (mask) {
          skip = 0;
          for (dc = -2; dc < 3; dc += 2) {
            for (dr = -2; dr < 3; dr += 2) {
              for (ds = -2; ds < 3; ds += 2) {
                m = MRIgetVoxVal(mask, c + dc, r + dr, s + ds, 0);
                if (m < 0.5) {
                  for (f = 0; f < 6; f++) MRIsetVoxVal(ar2, c, r, s, f, 0);
                  skip = 1;
                }
              }
            }
          }
          if (skip) continue;
        }

        // Loop thru all frames
        c1sum = 0;
        c2sum = 0;
        r1sum = 0;
        r2sum = 0;
        s1sum = 0;
        s2sum = 0;
        for (f = 0; f < srctmp->nframes; f++) {
          v0 = MRIgetVoxVal(srctmp, c, r, s, f);  // value at center voxel

          // temporal correlation with vox two col to left
          vc1 = MRIgetVoxVal(srctmp, c - 2, r, s, f);
          c1sum += (v0 * vc1);

          // temporal correlation with vox two col to right
          vc2 = MRIgetVoxVal(srctmp, c + 2, r, s, f);
          c2sum += (v0 * vc2);

          // temporal correlation with vox two row to up
          vr1 = MRIgetVoxVal(srctmp, c, r - 2, s, f);
          r1sum += (v0 * vr1);

          // temporal correlation with vox two row to down
          vr2 = MRIgetVoxVal(srctmp, c, r + 2, s, f);
          r2sum += (v0 * vr2);

          // temporal correlation with vox two slice in
          vs1 = MRIgetVoxVal(srctmp, c, r, s - 2, f);
          s1sum += (v0 * vs1);

          // temporal correlation with vox two slice out
          vs2 = MRIgetVoxVal(srctmp, c, r, s + 2, f);
          s2sum += (v0 * vs2);
        }

        // variance at center voxel
        v0 = MRIgetVoxVal(srcvar, c, r, s, 0);

        // column AR1
        vc1 = MRIgetVoxVal(srcvar, c - 2, r, s, 0);  // variance
        if (vc1 > 1e-6)
          car2 = c1sum / (nframes * sqrt(v0 * vc1));
        else
          car2 = 0;
        MRIsetVoxVal(ar2, c, r, s, 0, car2);  // frame 0

        vc2 = MRIgetVoxVal(srcvar, c + 2, r, s, 0);
        if (vc2 > 1e-6)
          car2 = c2sum / (nframes * sqrt(v0 * vc2));
        else
          car2 = 0;
        MRIsetVoxVal(ar2, c, r, s, 1, car2);  // frame 1

        // rows
        vr1 = MRIgetVoxVal(srcvar, c, r - 2, s, 0);
        if (vr1 > 1e-6)
          rar2 = r1sum / (nframes * sqrt(v0 * vr1));
        else
          rar2 = 0;
        MRIsetVoxVal(ar2, c, r, s, 2, rar2);  // frame 2

        vr2 = MRIgetVoxVal(srcvar, c, r + 2, s, 0);
        if (vr2 > 1e-6)
          rar2 = r2sum / (nframes * sqrt(v0 * vr2));
        else
          rar2 = 0;
        MRIsetVoxVal(ar2, c, r, s, 3, rar2);  // frame 3

        // slices
        vs1 = MRIgetVoxVal(srcvar, c, r, s - 2, 0);
        if (vs1 > 1e-6)
          sar2 = s1sum / (nframes * sqrt(v0 * vs1));
        else
          sar2 = 0;
        MRIsetVoxVal(ar2, c, r, s, 4, sar2);  // frame 4

        vs2 = MRIgetVoxVal(srcvar, c, r, s + 2, 0);
        if (vs2 > 1e-6)
          sar2 = s2sum / (nframes * sqrt(v0 * vs2));
        else
          sar2 = 0;
        MRIsetVoxVal(ar2, c, r, s, 5, sar2);  // frame 5

      }  // s
    }    // r
  }      // c

  MRIfree(&srcvar);
  if (freetmp) MRIfree(&srctmp);
  return (ar2);
}

/*----------------------------------------------------------
  fMRIspatialAR1Mean() - computes gobal mean of spatial AR1
  for col, row, and slice separately.
  ----------------------------------------------------------*/
int fMRIspatialAR1Mean(MRI *ar1, MRI *mask, double *car1mn, double *rar1mn, double *sar1mn)
{
  int c, r, s;
  long nhits;
  double m, car1sum, rar1sum, sar1sum;

  car1sum = 0.0;
  rar1sum = 0.0;
  sar1sum = 0.0;
  nhits = 0;
  for (c = 1; c < ar1->width - 1; c++) {
    for (r = 1; r < ar1->height - 1; r++) {
      for (s = 1; s < ar1->depth - 1; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        if (MRIgetVoxVal(ar1, c, r, s, 0) == 0) continue;

        car1sum += MRIgetVoxVal(ar1, c, r, s, 0);
        car1sum += MRIgetVoxVal(ar1, c, r, s, 1);

        rar1sum += MRIgetVoxVal(ar1, c, r, s, 2);
        rar1sum += MRIgetVoxVal(ar1, c, r, s, 3);

        sar1sum += MRIgetVoxVal(ar1, c, r, s, 4);
        sar1sum += MRIgetVoxVal(ar1, c, r, s, 5);

        nhits++;
      }
    }
  }

  *car1mn = (car1sum / (2 * nhits));
  *rar1mn = (rar1sum / (2 * nhits));
  *sar1mn = (sar1sum / (2 * nhits));

  return (0);
}

/*----------------------------------------------------------
  fMRIspatialAR2Mean() - computes gobal mean of spatial AR2
  for col, row, and slice separately.
  ----------------------------------------------------------*/
int fMRIspatialAR2Mean(MRI *src, MRI *mask, double *car2mn, double *rar2mn, double *sar2mn)
{
  int c, r, s;
  long nhits;
  double m, car2sum, rar2sum, sar2sum;
  MRI *ar2;

  ar2 = fMRIspatialAR2(src, mask, NULL);
  if (ar2 == NULL) return (1);

  car2sum = 0.0;
  rar2sum = 0.0;
  sar2sum = 0.0;
  nhits = 0;
  for (c = 1; c < src->width - 1; c++) {
    for (r = 1; r < src->height - 1; r++) {
      for (s = 1; s < src->depth - 1; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        if (MRIgetVoxVal(ar2, c, r, s, 0) == 0) continue;

        car2sum += MRIgetVoxVal(ar2, c, r, s, 0);
        car2sum += MRIgetVoxVal(ar2, c, r, s, 1);

        rar2sum += MRIgetVoxVal(ar2, c, r, s, 2);
        rar2sum += MRIgetVoxVal(ar2, c, r, s, 3);

        sar2sum += MRIgetVoxVal(ar2, c, r, s, 4);
        sar2sum += MRIgetVoxVal(ar2, c, r, s, 5);

        nhits++;
      }
    }
  }

  *car2mn = (car2sum / (2 * nhits));
  *rar2mn = (rar2sum / (2 * nhits));
  *sar2mn = (sar2sum / (2 * nhits));

  return (0);
}

/*!
  \fn MRI *fMRIaddOffset(MRI *in, MRI *offset, MRI *mask, MRI *out)
  \brief Add voxel-wise offset to input.
  \param in - input volume
  \param offset - offset volume (first frame used)
  \param mask - add offset only in mask (otherwise set to input)
  \param out - prealloced output (or NULL)
  \return out
*/
MRI *fMRIaddOffset(MRI *in, MRI *offset, MRI *mask, MRI *out)
{
  int c, r, s, f;
  double val0, val, m;

  if (MRIdimMismatch(in, offset, 0)) {
    printf("ERROR: fMRIaddOffset: input/offset dim mismatch\n");
    return (NULL);
  }

  if (out == NULL) {
    out = MRIcloneBySpace(in, MRI_FLOAT, -1);
    if (out == NULL) return (NULL);
  }
  else {
    if (MRIdimMismatch(in, out, 1)) {
      printf("ERROR: fMRIaddOffset: input/output dim mismatch\n");
      return (NULL);
    }
  }
  if (in != out) MRIcopy(in, out);

  for (c = 0; c < in->width; c++) {
    for (r = 0; r < in->height; r++) {
      for (s = 0; s < in->depth; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        val0 = MRIgetVoxVal(offset, c, r, s, 0);
        for (f = 0; f < in->nframes; f++) {
          val = MRIgetVoxVal(in, c, r, s, f);
          MRIsetVoxVal(out, c, r, s, f, val + val0);
        }
      }
    }
  }
  return (out);
}
/*-----------------------------------------------------------------*/
MRI *fMRIsubSample(MRI *f, int Start, int Delta, int Stop, MRI *fsub)
{
  int nframessub;
  int frame, subframe;
  int r, c, s;
  double v;

  if (Stop < 0) Stop = f->nframes;

  if (Start > Stop) {
    printf("ERROR: fMRIsubSample(): Start > Stop\n");
    return (NULL);
  }
  if (Delta == 0) {
    printf("ERROR: fMRIsubSample(): Delta == 0\n");
    return (NULL);
  }

  nframessub = ceil(((double)Stop - Start) / Delta);

  if (!fsub) fsub = MRIcloneBySpace(f, MRI_FLOAT, nframessub);
  if (nframessub != fsub->nframes) {
    printf("ERROR: fMRIsubSample(): frame mismatch (%d, %d)\n", nframessub, fsub->nframes);
    return (NULL);
  }

  for (c = 0; c < f->width; c++) {
    for (r = 0; r < f->height; r++) {
      for (s = 0; s < f->depth; s++) {
        subframe = 0;
        for (frame = Start; frame < Stop; frame += Delta) {
          v = MRIgetVoxVal(f, c, r, s, frame);
          MRIsetVoxVal(fsub, c, r, s, subframe, v);
          subframe++;
        }
      }
    }
  }

  return (fsub);
}
/*!
  \fn MRI *fMRIexcludeFrames(MRI *f, int *ExcludeFrames, int nExclude, MRI *fex)
  \brief Creates a new MRI by excluding the given set of rows.
*/
MRI *fMRIexcludeFrames(MRI *f, int *ExcludeFrames, int nExclude, MRI *fex)
{
  int skip, m, nframesNew, c, r, s, frame, subframe;
  double v;

  nframesNew = f->nframes - nExclude;
  if (!fex) fex = MRIcloneBySpace(f, MRI_FLOAT, nframesNew);
  if (nframesNew != fex->nframes) {
    printf("ERROR: fMRIexcludeFrames(): frame mismatch (%d, %d)\n", nframesNew, fex->nframes);
    return (NULL);
  }

  for (c = 0; c < f->width; c++) {
    for (r = 0; r < f->height; r++) {
      for (s = 0; s < f->depth; s++) {
        subframe = 0;
        for (frame = 0; frame < f->nframes; frame++) {
          skip = 0;
          for (m = 0; m < nExclude; m++)
            if (frame == ExcludeFrames[m]) skip = 1;
          if (skip) continue;
          v = MRIgetVoxVal(f, c, r, s, frame);
          MRIsetVoxVal(fex, c, r, s, subframe, v);
          subframe++;
        }
      }
    }
  }
  return (fex);
}
/*!
  \fn MRI *fMRItemporalGaussian(MRI *src, double gstdmsec, MRI *targ)
  \brief Temporal gaussian smoothing.
  \param src - source volume
  \param gstdmsec - gaussian stddev in milisec (divided by src->tr)
  \param targ - output
*/
MRI *fMRItemporalGaussian(MRI *src, double gstdmsec, MRI *targ)
{
  MATRIX *G, *v;
  int c, r, s, f;
  double sum;

  if (targ == NULL) {
    targ = MRIallocSequence(src->width, src->height, src->depth, MRI_FLOAT, src->nframes);
    if (targ == NULL) {
      printf("ERROR: MRItemporalGaussian: could not alloc\n");
      return (NULL);
    }
    MRIcopy(src, targ);
  }
  else {
    if (src->width != targ->width) {
      printf("ERROR: MRItemporalGaussian: width dimension mismatch\n");
      return (NULL);
    }
    if (src->height != targ->height) {
      printf("ERROR: MRItemporalGaussian: height dimension mismatch\n");
      return (NULL);
    }
    if (src->depth != targ->depth) {
      printf("ERROR: MRItemporalGaussian: depth dimension mismatch\n");
      return (NULL);
    }
    if (src->nframes != targ->nframes) {
      printf("ERROR: MRItemporalGaussian: frames dimension mismatch\n");
      return (NULL);
    }
    if (src != targ) MRIcopy(src, targ);
  }

  // Normalized at the edges
  G = GaussianMatrix(src->nframes, gstdmsec / src->tr, 1, NULL);
  for (r = 0; r < src->nframes; r++) {
    sum = 0.0;
    for (c = 0; c < src->nframes; c++) sum += G->rptr[r + 1][c + 1];
    for (c = 0; c < src->nframes; c++) G->rptr[r + 1][c + 1] /= sum;
  }
  // MatrixWriteTxt("g.txt",G);

  v = MatrixAlloc(src->nframes, 1, MATRIX_REAL);
  for (c = 0; c < src->width; c++) {
    for (r = 0; r < src->height; r++) {
      for (s = 0; s < src->depth; s++) {
        for (f = 0; f < src->nframes; f++) v->rptr[f + 1][1] = MRIgetVoxVal(src, c, r, s, f);
        MatrixMultiplyD(G, v, v);
        for (f = 0; f < src->nframes; f++) MRIsetVoxVal(targ, c, r, s, f, v->rptr[f + 1][1]);
      }
    }
  }

  MatrixFree(&G);
  MatrixFree(&v);

  return (targ);
}

MRI *fMRIkurtosis(MRI *y, MRI *mask)
{
  MRI *k;
  int c, r, s, f;
  double v, mn, m4 = 0, m2 = 0, g2, delta, b1, b2, n;
  k = MRIallocSequence(y->width, y->height, y->depth, MRI_FLOAT, 1);
  MRIcopyHeader(y, k);

  n = y->nframes;
  b1 = (n + 1) * (n - 1) / ((n - 2) * (n - 3));
  b2 = ((n - 1) * (n - 1)) / ((n - 2) * (n - 3));

  for (c = 0; c < y->width; c++) {
    for (r = 0; r < y->height; r++) {
      for (s = 0; s < y->depth; s++) {
        if (mask) {
          v = MRIgetVoxVal(mask, c, r, s, 0);
          if (v < 0.0001) {
            MRIsetVoxVal(k, c, r, s, 0, 0.0);
            continue;
          }
        }
        mn = 0;
        for (f = 0; f < y->nframes; f++) mn += MRIgetVoxVal(y, c, r, s, f);
        mn /= y->nframes;
        m2 = 0;
        m4 = 0;
        for (f = 0; f < y->nframes; f++) {
          delta = mn - MRIgetVoxVal(y, c, r, s, f);
          m2 += pow(delta, 2.0);  // sum of squares
          m4 += pow(delta, 4.0);  // sum of quads
        }
        m4 *= y->nframes;
        if (m2 != 0)
          g2 = b1 * (m4 / (m2 * m2)) - 3 * b2;  // 0 mean
        else
          g2 = 0;
        MRIsetVoxVal(k, c, r, s, 0, g2);
      }
    }
  }
  return (k);
}

/*!
  \fn MRI *MRIpkurtosis(MRI *kvals, int dof, MRI *mask, int nsamples)
  \brief Computes the p-value for a kurtosis map. Uses simulation.
  \param kvals - source volume with kurtosis values.
  \param dof - dof that the kurtosis was computed from
  \param mask - skip voxels where mask < 0.0001
  \param nsamples - samples to use in the simulation (eg, 10000)
*/
MRI *MRIpkurtosis(MRI *kvals, int dof, MRI *mask, int nsamples)
{
  MRI *nmri, *kmri, *pkmri;
  double *ksynth, pk, kvox, v;
  int m, c, r, s, f, ind;

  nmri = MRIrandn(nsamples, 1, 1, dof, 0, 1, NULL);
  kmri = fMRIkurtosis(nmri, NULL);

  ksynth = (double *)calloc(sizeof(double), nsamples);
  for (m = 0; m < nsamples; m++) ksynth[m] = MRIgetVoxVal(kmri, m, 0, 0, 0);

  qsort(ksynth, nsamples, sizeof(double), CompareDoubles);

  pkmri = MRIclone(kvals, NULL);
  for (c = 0; c < pkmri->width; c++) {
    for (r = 0; r < pkmri->height; r++) {
      for (s = 0; s < pkmri->depth; s++) {
        if (mask) {
          v = MRIgetVoxVal(mask, c, r, s, 0);
          if (v < 0.0001) {
            for (f = 0; f < pkmri->nframes; f++) MRIsetVoxVal(pkmri, c, r, s, f, 0.0);
            continue;
          }
        }
        for (f = 0; f < pkmri->nframes; f++) {
          kvox = MRIgetVoxVal(kvals, c, r, s, f);
          ind = PDFsearchOrderedTable(kvox, ksynth, nsamples);
          pk = 1.0 - (double)ind / nsamples;
          MRIsetVoxVal(pkmri, c, r, s, f, -log10(pk));
        }
      }
    }
  }
  MRIfree(&nmri);
  MRIfree(&kmri);
  free(ksynth);
  return (pkmri);
}

/*!
  \fn MATRIX *ASLinterpMatrix(). Creates a matrix that will linearly.
  \brief Interpolate the controls to esimate their values during the
   tag and then subtract the tag.  Assumes tag is first tp, and then
   every other.
  \param ntp - number of time points
*/
MATRIX *ASLinterpMatrix(int ntp)
{
  int r, c0, nrows, IsOdd;
  MATRIX *M;

  IsOdd = ntp % 2;

  nrows = ceil(ntp / 2.0);

  M = MatrixAlloc(nrows, ntp, MATRIX_REAL);
  for (r = 1; r <= nrows; r++) {
    c0 = 2 * (r - 1) + 1;
    if (r == 1) {
      M->rptr[r][c0] = +1.0;
      M->rptr[r][c0 + 1] = -1.0;
      continue;
    }
    if (r == nrows && IsOdd) {
      M->rptr[r][c0] = +1.0;
      M->rptr[r][c0 - 1] = -1.0;
      continue;
    }
    M->rptr[r][c0 - 1] = -0.5;
    M->rptr[r][c0] = +1.0;
    M->rptr[r][c0 + 1] = -0.5;
  }
  return (M);
}

/*!
  \fn MATRIX *fMRItoMatrix(MRI *fmri, MATRIX *M)
  \brief Stuffs an MRI into a Matrix. Frames goto rows,
  spatial dims go to cols. Note that cols are fastest, etc.
  This is especially important when comparing to matlab.
  Make sure to use fMRIfromMatrix() to undo it.
*/
MATRIX *fMRItoMatrix(MRI *fmri, MATRIX *M)
{
  int nthcol, nvox, c, r, s, f;
  double v;

  nvox = fmri->width * fmri->height * fmri->depth;

  if (M == NULL) {
    printf("fMRItoMatrix: allocating %d %d\n", fmri->nframes, nvox);
    M = MatrixAlloc(fmri->nframes, nvox, MATRIX_REAL);
    if (M == NULL) {
      printf("fMRItoMatrix: could not alloc\n");
      return (NULL);
    }
  }

  printf("fMRItoMatrix: filling matrix %d %d\n", fmri->nframes, nvox);
  nthcol = 0;
  for (s = 0; s < fmri->depth; s++) {
    for (r = 0; r < fmri->height; r++) {
      for (c = 0; c < fmri->width; c++) {
        for (f = 0; f < fmri->nframes; f++) {
          v = MRIgetVoxVal(fmri, c, r, s, f);
          M->rptr[f + 1][nthcol + 1] = v;
        }
        nthcol++;
      }
    }
  }
  return (M);
}
/*!
  \fn fMRIfromMatrix(MATRIX *M, MRI *fmri)
  \brief Stuffs a Matrix into an MRI. Rows goto frames.
  Cols go to spatial dims. Note that cols are fastest, etc.
  This is especially important when comparing to matlab.
  Make sure to use fMRItoMatrix() to undo it. fmri cannot
  be NULL!
*/
int fMRIfromMatrix(MATRIX *M, MRI *fmri)
{
  int nthcol, nvox, c, r, s, f;
  double v;

  nvox = fmri->width * fmri->height * fmri->depth;

  printf("fMRIfromMatrix: filling fMRI %d %d\n", fmri->nframes, nvox);
  nthcol = 0;
  for (s = 0; s < fmri->depth; s++) {
    for (r = 0; r < fmri->height; r++) {
      for (c = 0; c < fmri->width; c++) {
        for (f = 0; f < fmri->nframes; f++) {
          v = M->rptr[f + 1][nthcol + 1];
          MRIsetVoxVal(fmri, c, r, s, f, v);
        }
        nthcol++;
      }
    }
  }
  return (0);
}

/*!
  \fn MRI *fMRIspatialCorMatrix(MRI *fmri)
  \brief Computes the spatial correlation matrix. The result has the
  same spatial dim as fmri and number of frames equal to the total
  number of voxels. This can be very big!  The time courses are
  de-meaned and normalized to unit variance making this a matrix of
  Pearson correlations.
*/
MRI *fMRIspatialCorMatrix(MRI *fmri)
{
  int nvox;
  MRI *scm;
  MATRIX *M, *Mt, *MtM;

  printf("fMRIspatialCorMatrix: creating matrix\n");
  M = fMRItoMatrix(fmri, NULL);
  if (M == NULL) return (NULL);
  // MatrixWrite(M, "M.mat", "M");

  printf("fMRIspatialCorMatrix: demeaing\n");
  M = MatrixDemean(M, M);
  if (M == NULL) return (NULL);
  // MatrixWrite(M, "Mdm.mat", "Mdm");

  printf("fMRIspatialCorMatrix: normalizing\n");
  M = MatrixNormalizeCol(M, M, NULL);
  if (M == NULL) return (NULL);
  // MatrixWrite(M, "Mnorm.mat", "Mnorm");

  printf("fMRIspatialCorMatrix: creating matrix transpose\n");
  Mt = MatrixTranspose(M, NULL);
  if (Mt == NULL) return (NULL);
  // MatrixWrite(Mt, "Mt.mat", "Mt");

  printf("fMRIspatialCorMatrix: multiplying\n");
  MtM = MatrixMultiplyD(Mt, M, NULL);
  // MatrixWrite(MtM, "MtM.mat", "MtM");

  MatrixFree(&M);
  MatrixFree(&Mt);

  printf("fMRIspatialCorMatrix: allocating SCM\n");
  nvox = fmri->width * fmri->height * fmri->depth;
  scm = MRIallocSequence(fmri->width, fmri->height, fmri->depth, MRI_FLOAT, nvox);
  if (scm == NULL) return (NULL);

  printf("fMRIspatialCorMatrix: filling SCM\n");
  fMRIfromMatrix(MtM, scm);

  MatrixFree(&MtM);
  return (scm);
}

/*!
  \fn MRI *fMRIdistance(MRI *mri, MRI *mask)
  \brief Treats each frame triple as an xyz to compute distance
  \param mri - source volume with triples
  \param mask - skip voxels where mask < 0.0001
*/
MRI *fMRIdistance(MRI *mri, MRI *mask)
{
  MRI *d;
  double dx, dy, dz, v;
  int c, r, s, f, fd;

  // should check nframes/3
  d = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, mri->nframes / 3);
  if (d == NULL) {
    printf("ERROR: fMRIdistance: could not alloc\n");
    return (NULL);
  }
  MRIcopyHeader(mri, d);

  printf("fMRIdistance(): Computing Distance\n");

  for (c = 0; c < mri->width; c++) {
    for (r = 0; r < mri->height; r++) {
      for (s = 0; s < mri->depth; s++) {
        if (mask) {
          if (MRIgetVoxVal(mask, c, r, s, 0) < 0.5) {
            MRIFseq_vox(d, c, r, s, 0) = 0;
            continue;
          }
        }
        fd = 0;
        for (f = 0; f < mri->nframes; f += 3) {
          dx = MRIgetVoxVal(mri, c, r, s, f + 0);
          dy = MRIgetVoxVal(mri, c, r, s, f + 1);
          dz = MRIgetVoxVal(mri, c, r, s, f + 2);
          v = sqrt(dx * dx + dy * dy + dz * dz);
          // printf("%5d %2d %2d %3d   %3d   %g %g %g  %g\n",c,r,s,f,fd,dx,dy,dz,v);
          MRIsetVoxVal(d, c, r, s, fd, v);
          fd++;
        }
      }
    }
  }
  return (d);
}

/*!
  \fn MRI *fMRIcumSum(MRI *inmri, MRI *mask, MRI *outmri)
  \brief Computes cumulative sum over frames.
*/
MRI *fMRIcumSum(MRI *inmri, MRI *mask, MRI *outmri)
{
  int c, r, s, f;
  double val;

  if (outmri == NULL) {
    outmri = MRIallocSequence(inmri->width, inmri->height, inmri->depth, MRI_FLOAT, inmri->nframes);
    if (outmri == NULL) {
      printf("ERROR: fMRIcumSum(): could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(inmri, outmri);
  }
  else {
    if (outmri->width != inmri->width || outmri->height != inmri->height || outmri->depth != inmri->depth ||
        outmri->nframes != inmri->nframes) {
      printf("ERROR: fMRIcumSum(): output dimension mismatch\n");
      return (NULL);
    }
    // Should check mask here too
  }

  for (s = 0; s < outmri->depth; s++) {
    for (r = 0; r < outmri->height; r++) {
      for (c = 0; c < outmri->width; c++) {
        if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
        val = 0;
        for (f = 0; f < outmri->nframes; f++) {
          val += MRIgetVoxVal(inmri, c, r, s, f);
          MRIsetVoxVal(outmri, c, r, s, f, val);
        }
      }
    }
  }
  return (outmri);
}

/*!
  \fn MRI *fMRIcumTrapZ(MRI *y, MATRIX *t, MRI *mask, MRI *yz)
  \brief Computes trapezoidal integration (like matlab cumtrapz)
*/
MRI *fMRIcumTrapZ(MRI *y, MATRIX *t, MRI *mask, MRI *yz)
{
  int c, r, s, f;
  double v, vprev, vsum, dt;

  if (yz == NULL) {
    yz = MRIallocSequence(y->width, y->height, y->depth, MRI_FLOAT, y->nframes);
    if (yz == NULL) {
      printf("ERROR: fMRIcumtrapz: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(y, yz);
  }

  for (c = 0; c < y->width; c++) {
    for (r = 0; r < y->height; r++) {
      for (s = 0; s < y->depth; s++) {
        if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) {
          for (f = 0; f < y->nframes; f++) MRIFseq_vox(yz, c, r, s, f) = 0;
          continue;
        }
        vsum = 0;
        vprev = MRIgetVoxVal(y, c, r, s, 0);
        for (f = 1; f < y->nframes; f++) {
          dt = t->rptr[f + 1][1] - t->rptr[f][1];
          v = MRIgetVoxVal(y, c, r, s, f);
          vsum += (dt * ((v + vprev) / 2));
          MRIsetVoxVal(yz, c, r, s, f, vsum);
          vprev = v;
        }
      }
    }
  }
  return (yz);
}
/*!
  \fn MATRIX *HalfLife2Weight(double HalfLifeMin, MATRIX *tSec)
  \brief Computes a weight vector based on an exponential decay
  with the given half life. The weight vector is normalized
  so that the sum of the squares of the weights is 1. This
  vector is ready to be used in a GLM.
*/
MATRIX *HalfLife2Weight(double HalfLifeMin, MATRIX *tSec)
{
  MATRIX *w;
  int n;
  double TDecayMin, TDecaySec, v, wsum, *wd;

  // Convert the half-life to a decay constant
  TDecayMin = -HalfLifeMin / log(.5);
  TDecaySec = 60 * TDecayMin;

  wd = (double *)calloc(tSec->rows, sizeof(double));
  wsum = 0;
  for (n = 0; n < tSec->rows; n++) {
    v = tSec->rptr[n + 1][1] * exp(-tSec->rptr[n + 1][1] / TDecaySec);
    wd[n] = v;
    wsum += v;
  }

  w = MatrixAlloc(tSec->rows, tSec->cols, MATRIX_REAL);
  for (n = 0; n < tSec->rows; n++) w->rptr[n + 1][1] = sqrt(wd[n] / wsum);

  free(wd);
  return (w);
}

/*!
  \fn MRI *MRIframeNorm(MRI *src, MRI *mask, MRI *fnorm)
  \brief normalizes (ie, removes mean, divides by sqrt
   of the sum of squares) the time course at each voxel.
*/
MRI *MRIframeNorm(MRI *src, MRI *mask, MRI *fnorm)
{
  MRI *mn, *var;
  int c, r, s, f;
  double v, mnv, sss;

  if (fnorm == NULL) {
    fnorm = MRIallocSequence(src->width, src->height, src->depth, MRI_FLOAT, src->nframes);
    if (fnorm == NULL) {
      printf("ERROR: fMRIframeNorm: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(src, fnorm);
  }

  mn = MRIframeMean(src, NULL);
  var = fMRIcovariance(src, 0, -1, mask, NULL);

  for (c = 0; c < src->width; c++) {
    for (r = 0; r < src->height; r++) {
      for (s = 0; s < src->depth; s++) {
        if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
        mnv = MRIgetVoxVal(mn, c, r, s, 0);
        sss = sqrt(MRIgetVoxVal(var, c, r, s, 0) * (src->nframes - 1));
        for (f = 0; f < src->nframes; f++) {
          v = MRIgetVoxVal(src, c, r, s, f);
          v = (v - mnv) / (sss + DBL_MIN);
          MRIsetVoxVal(fnorm, c, r, s, f, v);
        }
      }
    }
  }

  MRIfree(&mn);
  MRIfree(&var);

  return (fnorm);
}

/*!
  \fn MRI *fMRIxcorr(MRI *v1, MRI *v2, MRI *mask, MRI *xcorr)
  \brief Computes voxelwise temporal correlation coefficient between
  two volumes. If v2==NULL, then v2 is formed from v1 by left-right
  reversing the columns.
*/
MRI *fMRIxcorr(MRI *v1, MRI *v2, MRI *mask, MRI *xcorr)
{
  MRI *fnorm1, *fnorm2 = NULL;
  int c, c2, r, s, f;
  double val1, val2, sum;

  if (Gdiag_no > 0) printf("MRIxcorr(): starting\n");
  fflush(stdout);

  if (xcorr == NULL) {
    xcorr = MRIallocSequence(v1->width, v1->height, v1->depth, MRI_FLOAT, 1);
    if (xcorr == NULL) {
      printf("ERROR: fMRIxcorr: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(v1, xcorr);
  }

  if (Gdiag_no > 0) printf("MRIxcorr(): normalizing vol 1\n");
  fflush(stdout);
  fnorm1 = MRIframeNorm(v1, mask, NULL);

  if (v2 != NULL) {
    if (Gdiag_no > 0) printf("MRIxcorr(): normalizing vol 2\n");
    fflush(stdout);
    fnorm2 = MRIframeNorm(v2, mask, NULL);
  }
  else
    printf("MRIxcorr(): reversing columns\n");
  fflush(stdout);

  if (Gdiag_no > 0) printf("MRIxcorr(): computing xcorr\n");
  fflush(stdout);
  for (c = 0; c < v1->width; c++) {
    c2 = v1->width - c - 1;  // for column reverse when v2=NULL
    for (r = 0; r < v1->height; r++) {
      for (s = 0; s < v1->depth; s++) {
        if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
        sum = 0;
        for (f = 0; f < v1->nframes; f++) {
          val1 = MRIgetVoxVal(fnorm1, c, r, s, f);
          if (v2 != NULL)
            val2 = MRIgetVoxVal(fnorm2, c, r, s, f);
          else
            val2 = MRIgetVoxVal(fnorm1, c2, r, s, f);
          sum += val1 * val2;
        }
        MRIsetVoxVal(xcorr, c, r, s, 0, sum);
      }
    }
  }

  MRIfree(&fnorm1);
  if (fnorm2) MRIfree(&fnorm2);
  if (Gdiag_no > 0) printf("MRIxcorr(): done\n");
  fflush(stdout);
  return (xcorr);
}

/*---------------------------------------------------------------------
  MRI *SpatialINorm(MRI *vol, MRI *mask, MRI *outvol) performs spatial
  intensity normalization by subtracting the spatial mean and dividing
  by the spatial stddev. If a mask is supplied, then only computes
  the mean and stddev within the mask.
  *--------------------------------------------------------------*/
MRI *SpatialINorm(MRI *vol, MRI *mask, MRI *outvol)
{
  int c, r, s, f, m;
  double gmean, gstddev, gmax, v;

  outvol = MRIclone(vol, outvol);

  RFglobalStats(vol, mask, &gmean, &gstddev, &gmax);
  printf("gmean = %lf, gstddev = %lf\n", gmean, gstddev);
  for (c = 0; c < vol->width; c++) {
    for (r = 0; r < vol->height; r++) {
      for (s = 0; s < vol->depth; s++) {
        if (mask != NULL) {
          m = (int)MRIgetVoxVal(mask, c, r, s, 0);
          if (!m) continue;
        }
        for (f = 0; f < vol->nframes; f++) {
          v = MRIgetVoxVal(vol, c, r, s, f);
          v = (v - gmean) / gstddev;
          MRIsetVoxVal(outvol, c, r, s, f, v);
        }
      }
    }
  }

  return (outvol);
}

/*---------------------------------------------------------------------
  fMRIspatialARN() - computes spatial ARN, ie, the correlation between
  the time course at one voxel and that at voxel N voxels away. There
  will be N+1 frames, where frame 0 is always 1.0.  This function
  assumess that the mean and any other trends have been removed. It
  works regardless of the DOF of the time series.
  --------------------------------------------------------------------*/
MRI *fMRIspatialARN(MRI *src, MRI *mask, int N, MRI *arN)
{
  int c, r, s, f, dc, dr, ds, nhits, lag;
  MRI *srcsumsq, *srctmp;
  double m, sumsq0, sumsqnbr, ar, arsum, sum, v0, vnbr;
  int freetmp;

  freetmp = 0;
  if (src->type != MRI_FLOAT) {
    srctmp = MRISeqchangeType(src, MRI_FLOAT, 0, 0, 0);
    freetmp = 1;
  }
  else {
    srctmp = src;
    freetmp = 0;
  }

  // alloc vol with N frames
  if (arN == NULL) {
    printf("N=%d\n", N);
    arN = MRIcloneBySpace(src, MRI_FLOAT, N + 1);
    if (arN == NULL) {
      printf("ERROR: could not alloc\n");
      return (NULL);
    }
  }

  // pre-compute the sum of squares
  srcsumsq = fMRIsumSquare(srctmp, 0, NULL);

  // Loop thru all voxels
  nhits = 0;
  for (c = 0; c < srctmp->width; c++) {
    for (r = 0; r < srctmp->height; r++) {
      for (s = 0; s < srctmp->depth; s++) {
        // skip voxel if it's on the edge
        if (c == 0 || r == 0 || s == 0 || c == (srctmp->width - 1) || r == (srctmp->height - 1) ||
            s == (srctmp->depth - 1)) {
          for (f = 0; f <= N; f++) MRIsetVoxVal(arN, c, r, s, f, 0);
          continue;
        }
        // skip voxel if it's not in the mask
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }

        // sum-of-sqares at center voxel
        sumsq0 = MRIgetVoxVal(srcsumsq, c, r, s, 0);
        if (sumsq0 < 1e-6) continue;

        MRIsetVoxVal(arN, c, r, s, 0, 1.0);
        for (lag = 1; lag <= N; lag++) {
          arsum = 0;
          nhits = 0;
          for (dc = -lag; dc < lag + 1; dc += (2 * lag)) {
            if (c + dc <= 0 || c + dc >= (srctmp->width - 1)) continue;
            if (mask) {
              m = MRIgetVoxVal(mask, c + dc, r, s, 0);
              if (m < 0.5) continue;
            }
            sumsqnbr = MRIgetVoxVal(srcsumsq, c + dc, r, s, 0);
            if (sumsqnbr < 1e-6) continue;
            nhits++;
            sum = 0;
            for (f = 0; f < srctmp->nframes; f++) {
              v0 = MRIgetVoxVal(srctmp, c, r, s, f);  // value at center voxel
              vnbr = MRIgetVoxVal(srctmp, c + dc, r, s, f);
              sum += v0 * vnbr;
            }
            ar = sum / sqrt(sumsq0 * sumsqnbr);
            arsum += ar;
          }  // dc
          for (dr = -lag; dr < lag + 1; dr += (2 * lag)) {
            if (r + dr <= 0 || r + dr >= (srctmp->height - 1)) continue;
            if (mask) {
              m = MRIgetVoxVal(mask, c, r + dr, s, 0);
              if (m < 0.5) continue;
            }
            sumsqnbr = MRIgetVoxVal(srcsumsq, c, r + dr, s, 0);
            if (sumsqnbr < 1e-6) continue;
            nhits++;
            sum = 0;
            for (f = 0; f < srctmp->nframes; f++) {
              v0 = MRIgetVoxVal(srctmp, c, r, s, f);  // value at center voxel
              vnbr = MRIgetVoxVal(srctmp, c, r + dr, s, f);
              sum += v0 * vnbr;
            }
            ar = sum / sqrt(sumsq0 * sumsqnbr);
            arsum += ar;
          }  // dr
          for (ds = -lag; ds < lag + 1; ds += (2 * lag)) {
            if (s + ds <= 0 || s + ds >= (srctmp->depth - 1)) continue;
            if (mask) {
              m = MRIgetVoxVal(mask, c, r, s + ds, 0);
              if (m < 0.5) continue;
            }
            sumsqnbr = MRIgetVoxVal(srcsumsq, c, r, s + ds, 0);
            if (sumsqnbr < 1e-6) continue;
            nhits++;
            sum = 0;
            for (f = 0; f < srctmp->nframes; f++) {
              v0 = MRIgetVoxVal(srctmp, c, r, s, f);  // value at center voxel
              vnbr = MRIgetVoxVal(srctmp, c, r, s + ds, f);
              sum += v0 * vnbr;
            }
            ar = sum / sqrt(sumsq0 * sumsqnbr);
            arsum += ar;
          }  // ds
          if (nhits > 0) MRIsetVoxVal(arN, c, r, s, lag, arsum / nhits);
        }  // lag

      }  // s
    }    // r
  }      // c

  MRIfree(&srcsumsq);
  if (freetmp) MRIfree(&srctmp);
  return (arN);
}
