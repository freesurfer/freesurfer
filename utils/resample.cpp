/**
 * @brief code to perform resapling from one space to another
 *
 * Purpose: code to perform resapling from one space to another,
 * including: volume-to-volume, volume-to-surface, and surface-to-surface.
 * Notes:
 * 1. Values on the surface are stored as MRI structures. In general,
 * they cannot be stored as part of the MRIS structure because this
 * does not accomodate multiple values per vertex. When stored as an
 * MRI, width = nvertices, height=1, depth=1, nframes=nframes.
 * 2. Interpolation - the code implies that trilinear and sinc
 * interpolation are available. They are not as of 2/4/01.
 * 3. Float2Int - the software supports three ways to convert floating
 * to integer when computing the index of a voxel: round, floor, tkreg.
 * Ideally, round should be used because it makes the mapping invertible.
 * However, tkregister uses it's own funky scheme (replicated with the
 * tkreg option), and paint uses floor. In the end, the conversion must
 * be compatible with the registration program.
 * 4. Volume-to-Volume - V2V is a necessary step when converting functional
 * data to talairach or painting onto the surface. The model as of 2/4/01
 * assumes that the transformation is completely linear from one space to
 * another, though there are variables for intermediate transformations
 * built in. The four variables are: (1) Quantization matrix Q, (2) Field-
 * of-view matrix F, (3) Warping matrix W, and (4), Registration matrix D.
 * D - Registration matrix. Converts from Anatomical Space to Unwarpded
 *     Scanner Space.
 * W - Warping matrix. Converts from Unwarped Scanner Space to Warped
 *     Scanner Space.
 * F - FOV matrix. Converts from Warpded Scanner Space to Field-of-View
 *     Space (the space created by the axes centered in the FOV and
 *     parallel with the edges of the FOV).
 * Q - Quantization matrix. Converts from FOV space to ColRowSlice Index
 *     Space.
 * The matrix in register.dat = D*W*F
 * In theory, the warping matrix can be replaced by a non-linear warping
 * function to account for warping of the scanner space, but this is yet
 * to be implemented.
 **/
/*
 * Original Author: Douglas N. Greve
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

#define RESAMPLE_SOURCE_CODE_FILE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "timer.h"

#include "romp_support.h"

#include "bfileio.h"
#include "corio.h"
#include "diag.h"
#include "label.h"
#include "matrix.h"
#include "mri.h"
#include "mri2.h"
#include "mrimorph.h"
#include "mrishash.h"
#include "mrisurf.h"
#include "proto.h"  // nint

#include "resample.h"

/*-------------------------------------------------------------------*/
double round(double);  // why is this never defined?!?
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------
  ASEGVOLINDEX - this structure is used with MRIaseg2vol() to help
  perform the mapping. Mainly used for sorting with qsort.
  -------------------------------------------------------------------*/
typedef struct
{
  int asegindex;   // index into the seg volume (instead of col,row,slice)
  int segid;       // segmentation code
  int volindex;    // index into the output volume
  int cv, rv, sv;  // corresponding col, row, slice in the output volume
} ASEGVOLINDEX;

static int CompareAVIndices(const void *i1, const void *i2);
static int MostHitsInVolVox(ASEGVOLINDEX *avindsorted, int N, int *segidmost, COLOR_TABLE *ct);

/*---------------------------------------------------------
  interpolation_code(): gets a code that controls the
  interpolation from a string. Returns -1 if the string is
  unrecoginzed.
  ---------------------------------------------------------*/
int interpolation_code(const char *interpolation_string) { return (MRIinterpCode(interpolation_string)); }
/*---------------------------------------------------------
  float2int_code(): gets a code that controls the floating
  point to integer conversion method from a string. Returns
  -1 if the string is unrecoginzed.
  ---------------------------------------------------------*/
int float2int_code(const char *float2int_string)
{
  if (!strcasecmp(float2int_string, "round") || !strcasecmp(float2int_string, "rint")) return (FLT2INT_ROUND);

  if (!strcasecmp(float2int_string, "floor")) return (FLT2INT_FLOOR);

  if (!strncasecmp(float2int_string, "tkregister", 5)) return (FLT2INT_TKREG);

  return (-1);
}

/*-----------------------------------------------------------
  XYZAnat2CRSFunc_TkReg() - computes the col, row, and slc in
  a functional volume given the x,y,z in the anatomical volume
  along with a registration matrix (generated by tkregister) and
  the dimensions of the functional volume:
    npixels  = number of rows = number of columns
    pixsize  = size of in-plane pixels (mm)
    nslcs    = number of slices
    slcthick = distance between slices (mm)
  When volumes are coronally sliced, then col is LR (sagital),
  row is SI (axial), and slc is AP (cor)
  -----------------------------------------------------------*/
int XYZAnat2CRSFunc_TkReg(int *col,
                          int *row,
                          int *slc,
                          int npixels,
                          float pixsize,
                          int nslcs,
                          float slcthick,
                          float xanat,
                          float yanat,
                          float zanat,
                          MATRIX *Reg)
{
  MATRIX *Qf;
  MATRIX *QfR;
  MATRIX *xyz, *crs;

  /* compute the tkregister quantization matrix */
  Qf = FOVQuantMtx_TkReg(npixels, pixsize, nslcs, slcthick);

  /* QfR = Qf* R*/
  QfR = MatrixMultiply(Qf, Reg, NULL);

  /* load xyz into a vector */
  xyz = MatrixAlloc(4, 1, MATRIX_REAL);
  xyz->rptr[0 + 1][0 + 1] = xanat;
  xyz->rptr[1 + 1][0 + 1] = yanat;
  xyz->rptr[2 + 1][0 + 1] = zanat;
  xyz->rptr[3 + 1][0 + 1] = 1.0;

  /* crs = Qf*R*xyz */
  crs = MatrixMultiply(QfR, xyz, NULL);

  /* extract col, row, slice from the crs vector */
  /* convert floats to ints as only tkregister can */
  float2int_TkReg(col, row, slc, crs->rptr[0 + 1][0 + 1], crs->rptr[1 + 1][0 + 1], crs->rptr[2 + 1][0 + 1]);

  MatrixFree(&Qf);
  MatrixFree(&QfR);
  MatrixFree(&xyz);
  MatrixFree(&crs);

  return (0);
}
/*----------------------------------------------------------
  float2int_TkReg() - converts floats to int as only
  tkregsiter can. When volumes are coronally sliced, then
  col is LR (sagital), row is SI (axial), and slc is AP (cor)
  ----------------------------------------------------------*/
int float2int_TkReg(int *col, int *row, int *slc, float fltcol, float fltrow, float fltslc)
{
  *col = (int)(fltcol);
  *row = (int)(ceil(fltrow));
  *slc = (int)(fltslc);
  return (0);
}
/*-----------------------------------------------------------
  FOVQuantMtx_TkReg() -- computes the FOV quantization matrix
  which converts an x,y,z into a col,row,slc : crs = Q*xyz. This
  particular formulation is the one used by TkRegister. Note:
  Q is a 4x4 matrix in order to accomodate translation.
    npixels  = number of rows = number of columns
    pixsize  = size of in-plane pixels (mm)
    nslcs    = number of slices
    slcthick = distance between slices (mm)
    -----------------------------------------------------------*/
MATRIX *FOVQuantMtx_TkReg(int npixels, float pixsize, int nslcs, float slcthick)
{
  MATRIX *Q;

  Q = MatrixAlloc(4, 4, MATRIX_REAL);
  MatrixClear(Q);

  Q->rptr[0 + 1][0 + 1] = -1 / pixsize;
  Q->rptr[0 + 1][3 + 1] = npixels / 2;

  Q->rptr[1 + 1][2 + 1] = 1 / slcthick;
  Q->rptr[1 + 1][3 + 1] = -nslcs / 2;

  Q->rptr[2 + 1][1 + 1] = -1 / pixsize;
  Q->rptr[2 + 1][3 + 1] = npixels / 2;

  Q->rptr[3 + 1][3 + 1] = 1.0;

  return (Q);
}
/*-----------------------------------------------------------
  FOVDeQuantMatrix() -- computes the volume dequantization matrix which
  converts a col,row,slc into an x,y,z: xyz = deQ*crs
-----------------------------------------------------------*/
MATRIX *FOVDeQuantMatrix(int ncols, int nrows, int nslcs, float colres, float rowres, float slcres)
{
  MATRIX *deQ;

  deQ = MatrixAlloc(4, 4, MATRIX_REAL);
  MatrixClear(deQ);

  deQ->rptr[0 + 1][0 + 1] = -colres;
  deQ->rptr[0 + 1][3 + 1] = colres * (ncols) / 2;

  deQ->rptr[1 + 1][2 + 1] = slcres;
  deQ->rptr[1 + 1][3 + 1] = -slcres * (nslcs) / 2;

  deQ->rptr[2 + 1][1 + 1] = -rowres;
  deQ->rptr[2 + 1][3 + 1] = rowres * (nrows) / 2;

  deQ->rptr[3 + 1][3 + 1] = 1.0;

  return (deQ);
}
/*-----------------------------------------------------------
  FOVQuantMatrix() -- computes the volume quantization matrix which
  converts a x,y,z into col,row,slc : crs = Q*xyz
  -----------------------------------------------------------*/
MATRIX *FOVQuantMatrix(int ncols, int nrows, int nslcs, float colres, float rowres, float slcres)
{
  MATRIX *deQ, *Q;

  deQ = FOVDeQuantMatrix(ncols, nrows, nslcs, colres, rowres, slcres);
  Q = MatrixInverse(deQ, NULL);
  MatrixFree(&deQ);

  return (Q);
}
/*------------------------------------------------------------
  ComputeQFWD() - computes the matrix product of Q, F, W, D.
  If any matrix is NULL, then it is treated as the idenity.
  ------------------------------------------------------------*/
MATRIX *ComputeQFWD(MATRIX *Q, MATRIX *F, MATRIX *W, MATRIX *D, MATRIX *QFWD)
{
  MATRIX *QFWDtmp;

  if (QFWD == NULL)
    QFWDtmp = MatrixAlloc(4, 4, MATRIX_REAL);
  else
    QFWDtmp = QFWD;

  MatrixIdentity(4, QFWDtmp);

  if (Q != NULL) MatrixMultiply(QFWDtmp, Q, QFWDtmp);
  if (F != NULL) MatrixMultiply(QFWDtmp, F, QFWDtmp);
  if (W != NULL) MatrixMultiply(QFWDtmp, W, QFWDtmp);
  if (D != NULL) MatrixMultiply(QFWDtmp, D, QFWDtmp);

  return (QFWDtmp);
}

/*------------------------------------------------------------
  vol2vol_linear() - DO NOT USE THIS FUNCTION. Use MRIvol2Vol()
  as found in mri2.c.

  Resample one volume into another assuming
  that the transformation is completely linear.  Msrc2trg is the
  linear mapping from the Anatomical Space of the Source Subject
  to the Anatomical Space of the Target Subject.
------------------------------------------------------------*/
MRI *vol2vol_linear(MRI *SrcVol,
                    MATRIX *Qsrc,
                    MATRIX *Fsrc,
                    MATRIX *Wsrc,
                    MATRIX *Dsrc,
                    MATRIX *Qtrg,
                    MATRIX *Ftrg,
                    MATRIX *Wtrg,
                    MATRIX *Dtrg,
                    int nrows_trg,
                    int ncols_trg,
                    int nslcs_trg,
                    MATRIX *Msrc2trg,
                    int InterpMethod,
                    int float2int)
{
  MATRIX *QFWDsrc, *QFWDtrg, *invQFWDtrg;
  MATRIX *Tcrs2Scrs, *invMsrc2trg;
  MATRIX *Scrs, *Tcrs, *Txyz;
  MRI *TrgVol;
  int irow_trg, icol_trg, islc_trg; /* integer row, col, slc in target */
  int irow_src, icol_src, islc_src; /* integer row, col, slc in source */
  float srcval;
  int frm;

  if (InterpMethod != INTERP_NEAREST) {
    fprintf(stderr, "vol2vol_linear(): only support for nearest interpolation\n");
    return (NULL);
  }

  /* compute the transforms */
  QFWDsrc = ComputeQFWD(Qsrc, Fsrc, Wsrc, Dsrc, NULL);
  QFWDtrg = ComputeQFWD(Qtrg, Ftrg, Wtrg, Dtrg, NULL);
  invQFWDtrg = MatrixInverse(QFWDtrg, NULL);
  if (Msrc2trg != NULL)
    invMsrc2trg = MatrixInverse(Msrc2trg, NULL);
  else
    invMsrc2trg = MatrixIdentity(4, NULL);
  Tcrs2Scrs = MatrixMultiply(QFWDsrc, invMsrc2trg, NULL);
  MatrixMultiply(Tcrs2Scrs, invQFWDtrg, Tcrs2Scrs);

  /* Tcrs2Scrs - maps Target ColRowSlice to that of the Source */

  /* preallocate the row-col-slc vectors */
  Tcrs = MatrixAlloc(4, 1, MATRIX_REAL);
  Tcrs->rptr[3 + 1][0 + 1] = 1.0;
  Scrs = MatrixAlloc(4, 1, MATRIX_REAL);
  Txyz = MatrixAlloc(4, 1, MATRIX_REAL);
  Txyz->rptr[3 + 1][0 + 1] = 1.0;

  /* allocate a volume to hold the output */
  TrgVol = MRIallocSequence(ncols_trg, nrows_trg, nslcs_trg, MRI_FLOAT, SrcVol->nframes);
  if (TrgVol == NULL) return (NULL);

  /* Go through each target voxel and compute the location of
     the closest source voxel */
  for (islc_trg = 0; islc_trg < nslcs_trg; islc_trg++) {
    for (irow_trg = 0; irow_trg < nrows_trg; irow_trg++) {
      for (icol_trg = 0; icol_trg < ncols_trg; icol_trg++) {
        /* Load the Target col-row-slc vector */
        Tcrs->rptr[0 + 1][0 + 1] = icol_trg;
        Tcrs->rptr[1 + 1][0 + 1] = irow_trg;
        Tcrs->rptr[2 + 1][0 + 1] = islc_trg;

        /* Compute the corresponding Source col-row-slc vector */
        MatrixMultiply(Tcrs2Scrs, Tcrs, Scrs);

        /* nearest neighbor */
        switch (float2int) {
          case FLT2INT_ROUND:
            icol_src = (int)rint(Scrs->rptr[0 + 1][0 + 1]);
            irow_src = (int)rint(Scrs->rptr[1 + 1][0 + 1]);
            islc_src = (int)rint(Scrs->rptr[2 + 1][0 + 1]);
            break;
          case FLT2INT_FLOOR:
            icol_src = (int)floor(Scrs->rptr[0 + 1][0 + 1]);
            irow_src = (int)floor(Scrs->rptr[1 + 1][0 + 1]);
            islc_src = (int)floor(Scrs->rptr[2 + 1][0 + 1]);
            break;
          case FLT2INT_TKREG:
            icol_src = (int)floor(Scrs->rptr[0 + 1][0 + 1]);
            irow_src = (int)ceil(Scrs->rptr[1 + 1][0 + 1]);
            islc_src = (int)floor(Scrs->rptr[2 + 1][0 + 1]);
            break;
          default:
            fprintf(stderr, "vol2vol_linear(): unrecoginized float2int code %d\n", float2int);
            MRIfree(&TrgVol);
            return (NULL);
            break;
        }

        /* make sure the Source Voxel is in the volume */
        if (irow_src < 0 || irow_src >= SrcVol->height || icol_src < 0 || icol_src >= SrcVol->width || islc_src < 0 ||
            islc_src >= SrcVol->depth)
          continue;

        /* map each frame */
        for (frm = 0; frm < SrcVol->nframes; frm++) {
          srcval = MRIFseq_vox(SrcVol, icol_src, irow_src, islc_src, frm);
          MRIFseq_vox(TrgVol, icol_trg, irow_trg, islc_trg, frm) = srcval;
        }
      }
    }
  }

  MatrixFree(&QFWDsrc);
  MatrixFree(&QFWDtrg);
  MatrixFree(&invQFWDtrg);
  MatrixFree(&Tcrs2Scrs);
  MatrixFree(&invMsrc2trg);
  MatrixFree(&Scrs);
  MatrixFree(&Tcrs);

  printf("mri_vol2vol_linear: done\n");

  return (TrgVol);
}
/*-----------------------------------------------------------------------
  label2mask_linear() - converts a label into a masking volume. The
  masking volume is the same FOV as the source volume.  A mask voxel
  has a value of 1 where ever the number of label points that fall
  into that voxel exceeds a threshold and the SrcMskVol is > 0.5; the
  mask is zero everywhere else.  The SrcMskVol is a mask of the same
  FOV as the source; it is ignored if NULL. The threshold is computed
  based on rszthresh and the voxel size. Each label point is assumed
  to occupy 1 cubic mm, so the volume of a voxel filled by the label
  points just equals the number of label points that fall into that
  voxel. This must exceed rszthresh*voxelsize (ie, rszthresh is the
  fraction of the voxel that must be filled by the label points inorder
  for the voxel to be included in the mask). To make a mask with as
  few as one hit in a voxel, just set rszthresh to a very small number.
  ------------------------------------------------------------------------*/
MRI *label2mask_linear(MRI *SrcVol,
                       MATRIX *Qsrc,
                       MATRIX *Fsrc,
                       MATRIX *Wsrc,
                       MATRIX *Dsrc,
                       MRI *SrcMskVol,
                       MATRIX *Msrc2lbl,
                       LABEL *Label,
                       float rszthresh,
                       int float2int,
                       int *nlabelhits,
                       int *nfinalhits)
{
  MATRIX *QFWDsrc;
  MATRIX *Lxyz2Scrs;
  MATRIX *Scrs, *Lxyz, *Mlbl2src;
  MRI *FinalMskVol;
  int irow_src, icol_src, islc_src; /* integer row, col, slc in source */
  float mskval, voxsize;
  int vlbl;
  int c, r, s, nfinalmask;

  if (SrcMskVol) {
    if (MRIdimMismatch(SrcVol, SrcMskVol, 0)) {
      printf("ERROR: label2mask_linear: dimension mismatch\n");
      return (NULL);
    }
  }

  /* compute the transforms */
  QFWDsrc = ComputeQFWD(Qsrc, Fsrc, Wsrc, Dsrc, NULL);
  if (Msrc2lbl != NULL)
    Mlbl2src = MatrixInverse(Msrc2lbl, NULL);
  else
    Mlbl2src = NULL;
  if (Mlbl2src != NULL)
    Lxyz2Scrs = MatrixMultiply(QFWDsrc, Mlbl2src, NULL);
  else
    Lxyz2Scrs = MatrixCopy(QFWDsrc, NULL);

  printf("\n");
  printf("Lxyz2Scrs:\n");
  MatrixPrint(stdout, Lxyz2Scrs);
  printf("\n");

  /* preallocate the row-col-slc vectors */
  Lxyz = MatrixAlloc(4, 1, MATRIX_REAL);
  Lxyz->rptr[3 + 1][0 + 1] = 1.0;
  Scrs = MatrixAlloc(4, 1, MATRIX_REAL);

  /* allocate an output volume -- same size as source*/
  FinalMskVol = MRIallocSequence(SrcVol->width, SrcVol->height, SrcVol->depth, MRI_FLOAT, 1);
  if (FinalMskVol == NULL) return (NULL);
  MRIcopyHeader(SrcVol, FinalMskVol);

  *nlabelhits = 0;
  *nfinalhits = 0;

  /* Go through each point in the label */
  for (vlbl = 0; vlbl < Label->n_points; vlbl++) {
    /* load the label xyz into a vector */
    Lxyz->rptr[0 + 1][0 + 1] = Label->lv[vlbl].x;
    Lxyz->rptr[1 + 1][0 + 1] = Label->lv[vlbl].y;
    Lxyz->rptr[2 + 1][0 + 1] = Label->lv[vlbl].z;

    /* compute the corresponding col, row, and slice in the source vol */
    MatrixMultiply(Lxyz2Scrs, Lxyz, Scrs);

    /* Convert the analog col, row, and slice to integer */
    switch (float2int) {
      case FLT2INT_ROUND:
        icol_src = (int)rint(Scrs->rptr[0 + 1][0 + 1]);
        irow_src = (int)rint(Scrs->rptr[1 + 1][0 + 1]);
        islc_src = (int)rint(Scrs->rptr[2 + 1][0 + 1]);
        break;
      case FLT2INT_FLOOR:
        icol_src = (int)floor(Scrs->rptr[0 + 1][0 + 1]);
        irow_src = (int)floor(Scrs->rptr[1 + 1][0 + 1]);
        islc_src = (int)floor(Scrs->rptr[2 + 1][0 + 1]);
        break;
      case FLT2INT_TKREG:
        icol_src = (int)floor(Scrs->rptr[0 + 1][0 + 1]);
        irow_src = (int)ceil(Scrs->rptr[1 + 1][0 + 1]);
        islc_src = (int)floor(Scrs->rptr[2 + 1][0 + 1]);
        break;
      default:
        fprintf(stderr, "label2mask_linear(): unrecoginized float2int code %d\n", float2int);
        MRIfree(&FinalMskVol);
        return (NULL);
        break;
    }

    /* check that the point is within the source volume */
    if (irow_src < 0 || irow_src >= SrcVol->height || icol_src < 0 || icol_src >= SrcVol->width || islc_src < 0 ||
        islc_src >= SrcVol->depth)
      continue;
    (*nlabelhits)++;

    /* check that the point is within the input mask */
    if (SrcMskVol != NULL) {
      mskval = MRIFseq_vox(SrcMskVol, icol_src, irow_src, islc_src, 0);
      if (mskval < 0.5) continue;
    }

    /* keep count; binarization is done below */
    MRIFseq_vox(FinalMskVol, icol_src, irow_src, islc_src, 0)++;

    /* increment the number of hits */
    (*nfinalhits)++;
  }

  printf("INFO: label2mask_linear: there were %d label hits\n", *nfinalhits);

  /* binarize based on the number of hits in each voxel */
  voxsize = SrcVol->xsize * SrcVol->ysize * SrcVol->zsize;
  printf("voxsize = %g, rszthresh = %g, thresh = %g\n", voxsize, rszthresh, voxsize * rszthresh);
  nfinalmask = 0;
  for (c = 0; c < FinalMskVol->width; c++) {
    for (r = 0; r < FinalMskVol->height; r++) {
      for (s = 0; s < FinalMskVol->depth; s++) {
        mskval = MRIFseq_vox(FinalMskVol, c, r, s, 0);
        if (mskval < rszthresh * voxsize)
          MRIFseq_vox(FinalMskVol, c, r, s, 0) = 0;
        else {
          MRIFseq_vox(FinalMskVol, c, r, s, 0) = 1;
          nfinalmask++;
        }
      }
    }
  }

  printf("INFO: label2mask_linear: there were %d hits in final mask\n", nfinalmask);

  MatrixFree(&QFWDsrc);
  MatrixFree(&Lxyz2Scrs);
  MatrixFree(&Scrs);
  MatrixFree(&Lxyz);
  if (Mlbl2src != NULL) MatrixFree(&Mlbl2src);

  printf("label2mask_linear: done\n");

  return (FinalMskVol);
}

/*-----------------------------------------------------------------------
  vol2maskavg() -- averages all the voxels within a mask.  SrcVol and
  SrcMskVol must be the same size.  Searchs SrcMskVol for voxels whose
  value is > 0.5. For all such voxels, averages the corresponding values
  found in SrcVol.  Returns an MRI "volume" of dimension 1X1X1Xnframes.
  -----------------------------------------------------------------------*/
MRI *vol2maskavg(MRI *SrcVol, MRI *SrcMskVol, int *nhits)
{
  int r, c, s, f;
  MRI *MskAvg;
  float mskval, val;

  /* make sure that SrcVol and SrcMskVol are the same dimension */
  if (SrcVol->width != SrcMskVol->width || SrcVol->height != SrcMskVol->height || SrcVol->depth != SrcMskVol->depth) {
    fprintf(stderr,
            "ERROR: vol2maskavg: SrcVol and SrcMskVol do "
            "not have the same dimension\n");
    return (NULL);
  }

  /* allocate an output "volume" */
  MskAvg = MRIallocSequence(1, 1, 1, MRI_FLOAT, SrcVol->nframes);
  if (MskAvg == NULL) return (NULL);

  /* go through each voxel */
  *nhits = 0;
  for (c = 0; c < SrcVol->width; c++) {
    for (r = 0; r < SrcVol->height; r++) {
      for (s = 0; s < SrcVol->depth; s++) {
        /* get mask value at r,c,s */
        mskval = MRIgetVoxVal(SrcMskVol, c, r, s, 0);
        if (mskval > 0.5) {
          /* accumulate sum over suprathreshold points */
          (*nhits)++;
          for (f = 0; f < SrcVol->nframes; f++) {
            val = MRIgetVoxVal(SrcVol, c, r, s, f);
            MRIFseq_vox(MskAvg, 0, 0, 0, f) += val;
          }
        }
      }
    }
  }

  printf("INFO: vol2maskavg: nhits = %d\n", *nhits);
  if (*nhits != 0) {
    /* divide by the number of hits to get the average*/
    for (f = 0; f < SrcVol->nframes; f++) {
      val = MRIFseq_vox(MskAvg, 0, 0, 0, f);
      MRIFseq_vox(MskAvg, 0, 0, 0, f) = val / (*nhits);
      // printf("%2d %g %g\n",f,val,MRIFseq_vox(MskAvg,0,0,0,f));
    }
  }
  else {
    printf("WARNING: there were no voxels in the input mask > 0.5\n");
  }

  /* return the Masked Average */
  return (MskAvg);
}
/*----------------------------------------------------------------
  ProjNormFracThick() - projects along the surface normal a given
  fraction of the thickness at that point.
  ----------------------------------------------------------------*/
int ProjNormFracThick(float *x, float *y, float *z, const MRI_SURFACE *surf, int vtx, float frac)
{
  float r;
  r = frac * surf->vertices[vtx].curv;
  *x = surf->vertices[vtx].x + r * surf->vertices[vtx].nx;
  *y = surf->vertices[vtx].y + r * surf->vertices[vtx].ny;
  *z = surf->vertices[vtx].z + r * surf->vertices[vtx].nz;
  return (0);
}
/*----------------------------------------------------------------
  ProjNormFracThickNbr() - projects along the direction of the average of
  the surface normals of a vertex and its neighbor.  The distance
  along this direction is a fraction of the thickness at that point.
  ----------------------------------------------------------------*/
int ProjNormFracThickNbr(float *x, float *y, float *z, MRI_SURFACE *surf, int vtxno, float frac, int nthNbr)
{
  float r, nx, ny, nz;
  int nbrvtxno;

  nbrvtxno = surf->vertices_topology[vtxno].v[nthNbr];
  nx = (surf->vertices[vtxno].nx + surf->vertices[nbrvtxno].nx) / 2.0;
  ny = (surf->vertices[vtxno].ny + surf->vertices[nbrvtxno].ny) / 2.0;
  nz = (surf->vertices[vtxno].nz + surf->vertices[nbrvtxno].nz) / 2.0;
  r = frac * surf->vertices[vtxno].curv;
  *x = surf->vertices[vtxno].x + r * nx;
  *y = surf->vertices[vtxno].y + r * ny;
  *z = surf->vertices[vtxno].z + r * nz;
  return (0);
}
/*----------------------------------------------------------------
  ProjNormDist() - projects along the surface normal a given
  distance.
  ----------------------------------------------------------------*/
int ProjNormDist(float *x, float *y, float *z, const MRI_SURFACE *surf, int vtx, float dist)
{
  *x = surf->vertices[vtx].x + dist * surf->vertices[vtx].nx;
  *y = surf->vertices[vtx].y + dist * surf->vertices[vtx].ny;
  *z = surf->vertices[vtx].z + dist * surf->vertices[vtx].nz;
  return (0);
}
/*----------------------------------------------------------------
  ProjNormDistNbr() - projects along the direction of the average of
  the surface normals of a vertex and its neighbor.  The distance
  along this direction is dist.
  ----------------------------------------------------------------*/
int ProjNormDistNbr(float *x, float *y, float *z, MRI_SURFACE *surf, int vtxno, float dist, int nthNbr)
{
  float nx, ny, nz;
  int nbrvtxno;

  nbrvtxno = surf->vertices_topology[vtxno].v[nthNbr];
  nx = (surf->vertices[vtxno].nx + surf->vertices[nbrvtxno].nx) / 2.0;
  ny = (surf->vertices[vtxno].ny + surf->vertices[nbrvtxno].ny) / 2.0;
  nz = (surf->vertices[vtxno].nz + surf->vertices[nbrvtxno].nz) / 2.0;
  *x = surf->vertices[vtxno].x + dist * nx;
  *y = surf->vertices[vtxno].y + dist * ny;
  *z = surf->vertices[vtxno].z + dist * nz;
  return (0);
}
/*------------------------------------------------------------
  vol2surf_linear() - resamples data from a volume onto surface
  vertices assuming the the transformation from the volume into
  anatomical space is fully linear. The voxels actually
  sampled are stored in the SrcHitVol. The value is the number
  of times the voxel was sampled. If ProjFrac is non-zero, then
  it will project along the surface normal a distance equal to
  a fraction ProjFrac of the surface thickness at that point.
  If ProjDistFlag is set then ProjFrac is interpreted as an
  absolute distance. Qsrc is the tkreg ras2vox (can be NULL),
  Dsrc is the register.dat, just set Fsrc and Wsrc to NULL.
  ------------------------------------------------------------*/
MRI *vol2surf_linear(MRI *SrcVol,
                     MATRIX *Qsrc,
                     MATRIX *Fsrc,
                     MATRIX *Wsrc,
                     MATRIX *Dsrc,
                     MRI_SURFACE *TrgSurf,
                     float ProjFrac,
                     int InterpMethod,
                     int float2int,
                     MRI *SrcHitVol,
                     int ProjDistFlag,
                     int nskip)
{
  MATRIX *QFWDsrc;
  MATRIX *Scrs, *Txyz;
  MRI *TrgVol;
  int irow_src, icol_src, islc_src;   /* integer row, col, slc in source */
  float frow_src, fcol_src, fslc_src; /* float row, col, slc in source */
  float srcval, Tx, Ty, Tz;
  int frm, FreeQsrc = 0;
  int vtx, nhits;
  float *valvect;
  double rval;

  if (Qsrc == NULL) {
    Qsrc = MRIxfmCRS2XYZtkreg(SrcVol);
    Qsrc = MatrixInverse(Qsrc, Qsrc);
    FreeQsrc = 1;
  }

  /* compute the transforms */
  QFWDsrc = ComputeQFWD(Qsrc, Fsrc, Wsrc, Dsrc, NULL);
  if (Gdiag_no >= 0) {
    printf("QFWDsrc: vol2surf: ------------------------------\n");
    MatrixPrint(stdout, QFWDsrc);
    printf("--------------------------------------------------\n");
  }

  /* preallocate the row-col-slc vectors */
  Scrs = MatrixAlloc(4, 1, MATRIX_REAL);
  Txyz = MatrixAlloc(4, 1, MATRIX_REAL);
  Txyz->rptr[3 + 1][0 + 1] = 1.0;

  /* allocate a "volume" to hold the output */
  TrgVol = MRIallocSequence(TrgSurf->nvertices, 1, 1, MRI_FLOAT, SrcVol->nframes);
  if (TrgVol == NULL) return (NULL);
  MRIcopyHeader(SrcVol, TrgVol);
  // Dims here are meaningless, but setting to 1 means "volume" will be
  // number of vertices.
  TrgVol->xsize = 1;
  TrgVol->ysize = 1;
  TrgVol->zsize = 1;

  /* Zero the source hit volume */
  if (SrcHitVol != NULL) {
    MRIconst(SrcHitVol->width, SrcHitVol->height, SrcHitVol->depth, 1, 0, SrcHitVol);
  }

  srcval = 0;
  valvect = (float *)calloc(sizeof(float), SrcVol->nframes);
  nhits = 0;
  /*--- loop through each vertex ---*/
  for (vtx = 0; vtx < TrgSurf->nvertices; vtx += nskip) {
    if (ProjFrac != 0.0)
      if (ProjDistFlag)
        ProjNormDist(&Tx, &Ty, &Tz, TrgSurf, vtx, ProjFrac);
      else
        ProjNormFracThick(&Tx, &Ty, &Tz, TrgSurf, vtx, ProjFrac);
    else {
      Tx = TrgSurf->vertices[vtx].x;
      Ty = TrgSurf->vertices[vtx].y;
      Tz = TrgSurf->vertices[vtx].z;
    }

    /* Load the Target xyz vector */
    Txyz->rptr[0 + 1][0 + 1] = Tx;
    Txyz->rptr[1 + 1][0 + 1] = Ty;
    Txyz->rptr[2 + 1][0 + 1] = Tz;

    /* Compute the corresponding Source col-row-slc vector */
    MatrixMultiply(QFWDsrc, Txyz, Scrs);
    fcol_src = Scrs->rptr[1][1];
    frow_src = Scrs->rptr[2][1];
    fslc_src = Scrs->rptr[3][1];

    /* nearest neighbor */
    switch (float2int) {
      case FLT2INT_ROUND:
        icol_src = nint(fcol_src);
        irow_src = nint(frow_src);
        islc_src = nint(fslc_src);
        break;
      case FLT2INT_FLOOR:
        icol_src = (int)floor(fcol_src);
        irow_src = (int)floor(frow_src);
        islc_src = (int)floor(fslc_src);
        break;
      case FLT2INT_TKREG:
        icol_src = (int)floor(fcol_src);
        irow_src = (int)ceil(frow_src);
        islc_src = (int)floor(fslc_src);
        break;
      default:
        fprintf(stderr, "vol2surf_linear(): unrecoginized float2int code %d\n", float2int);
        MRIfree(&TrgVol);
        return (NULL);
        break;
    }

    /* check that the point is in the bounds of the volume */
    if (irow_src < 0 || irow_src >= SrcVol->height || icol_src < 0 || icol_src >= SrcVol->width || islc_src < 0 ||
        islc_src >= SrcVol->depth)
      continue;

    if (Gdiag_no == vtx) {
      printf("diag -----------------------------\n");
      printf("vtx = %d  %g %g %g\n", vtx, Tx, Ty, Tz);
      printf("fCRS  %g %g %g\n", Scrs->rptr[1][1], Scrs->rptr[2][1], Scrs->rptr[3][1]);
      printf("CRS  %d %d %d\n", icol_src, irow_src, islc_src);
    }

    /* only gets here if it is in bounds */
    nhits++;

    /* Assign output volume values */
    if (InterpMethod == SAMPLE_TRILINEAR) {
      MRIsampleSeqVolume(SrcVol, fcol_src, frow_src, fslc_src, valvect, 0, SrcVol->nframes - 1);
      if (Gdiag_no == vtx) printf("val = %f\n", valvect[0]);
      for (frm = 0; frm < SrcVol->nframes; frm++) MRIFseq_vox(TrgVol, vtx, 0, 0, frm) = valvect[frm];
    }
    else {
      for (frm = 0; frm < SrcVol->nframes; frm++) {
        switch (InterpMethod) {
          case SAMPLE_NEAREST:
            srcval = MRIgetVoxVal(SrcVol, icol_src, irow_src, islc_src, frm);
            break;
            srcval = MRIgetVoxVal(SrcVol, icol_src, irow_src, islc_src, frm);
            break;
          case SAMPLE_SINC: /* no multi-frame */
            MRIsincSampleVolume(SrcVol, fcol_src, frow_src, fslc_src, 5, &rval);
            srcval = rval;
            break;
        }  // switch
        MRIFseq_vox(TrgVol, vtx, 0, 0, frm) = srcval;
        if (Gdiag_no == vtx) printf("val[%d] = %f\n", frm, srcval);
      }  // for
    }    // else
    if (SrcHitVol != NULL) MRIFseq_vox(SrcHitVol, icol_src, irow_src, islc_src, 0)++;
  }

  MatrixFree(&QFWDsrc);
  MatrixFree(&Scrs);
  MatrixFree(&Txyz);
  free(valvect);
  if (FreeQsrc) MatrixFree(&Qsrc);

  // printf("vol2surf_linear: nhits = %d/%d\n",nhits,TrgSurf->nvertices);

  return (TrgVol);
}

/*!
\fn MRI *MRISapplyReg(MRI *SrcSurfVals, MRI_SURFACE **SurfReg, int nsurfs,
                  int ReverseMapFlag, int DoJac, int UseHash)
\brief Applies one or more surface registrations with or without jacobian correction.
This should be used as a replacement for surf2surf_nnfr and surf2surf_nnfr_jac
(it gives identical results).
\param MRI *SrcSurfVals - Inputs
\param MRIS **SurfReg - array of surface reg pairs, src1-trg1:src2-trg2:... where
trg1 and src2 are from the same anatomy.
\param int nsurfs - total number of surfs in SurfReg
\param int ReverseMapFlag - perform reverse mapping
\param int DoJac - perform jacobian correction (conserves sum(SrcVals))
\param int UseHash - use hash table (no reason not to, much faster).
*/
MRI *MRISapplyReg(MRI *SrcSurfVals, MRI_SURFACE **SurfReg, int nsurfs, int ReverseMapFlag, int DoJac, int UseHash)
{
  MRI *TrgSurfVals = NULL;
  MRI_SURFACE *SrcSurfReg, *TrgSurfReg;
  int svtx = 0, tvtx, tvtxN, svtxN = 0, f, n, nrevhits, nSrcLost;
  int npairs, kS, kT, nhits;
  // int nunmapped;
  VERTEX *v;
  float dmin;
  MHT **Hash = NULL;
  MRI *SrcHits, *TrgHits;

  npairs = nsurfs / 2;
  printf("MRISapplyReg(): nsurfs = %d, revmap=%d, jac=%d,  hash=%d\n", nsurfs, ReverseMapFlag, DoJac, UseHash);
  printf("  Skipping ripped vertices\n");

  SrcSurfReg = SurfReg[0];
  TrgSurfReg = SurfReg[nsurfs - 1];

  /* check dimension consistency */
  if (SrcSurfVals->width != SrcSurfReg->nvertices) {
    printf("MRISapplyReg: Vals and Reg dimension mismatch\n");
    printf("nVals = %d, nReg %d\n", SrcSurfVals->width, SrcSurfReg->nvertices);
    return (NULL);
  }
  for (n = 0; n < npairs - 1; n++) {
    kS = 2 * n + 1;
    kT = kS + 1;
    if (SurfReg[kT]->nvertices != SurfReg[kS]->nvertices) {
      printf("MRISapplyReg: Reg dimension mismatch %d, %d\n", kT, kS);
      printf("targ = %d, next source = %d\n", SurfReg[kT]->nvertices, SurfReg[kS]->nvertices);
      return (NULL);
    }
  }

  /* allocate a "volume" to hold the output */
  TrgSurfVals = MRIallocSequence(TrgSurfReg->nvertices, 1, 1, MRI_FLOAT, SrcSurfVals->nframes);
  if (TrgSurfVals == NULL) return (NULL);
  MRIcopyHeader(SrcSurfVals, TrgSurfVals);

  /* number of source vertices mapped to each target vertex */
  TrgHits = MRIallocSequence(TrgSurfReg->nvertices, 1, 1, MRI_FLOAT, 1);
  if (TrgHits == NULL) return (NULL);
  MRIcopyHeader(SrcSurfVals, TrgHits);

  /* number of target vertices mapped to by each source vertex */
  SrcHits = MRIallocSequence(SrcSurfReg->nvertices, 1, 1, MRI_FLOAT, 1);
  if (SrcHits == NULL) return (NULL);
  MRIcopyHeader(SrcSurfVals, SrcHits);

  if (UseHash) {
    printf("MRISapplyReg: building hash tables (res=16).\n");
    Hash = (MHT **)calloc(sizeof(MHT *), nsurfs);
    for (n = 0; n < nsurfs; n++) {
      Hash[n] = MHTcreateVertexTable_Resolution(SurfReg[n], CURRENT_VERTICES, 16);
    }
  }

  if (DoJac) {
    // If using jacobian correction, get a list of the number of times
    // that a give source vertex gets sampled.
    for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++) {
      // Compute the source vertex that corresponds to this target vertex
      tvtxN = tvtx;
      for (n = npairs - 1; n >= 0; n--) {
        kS = 2 * n;
        kT = kS + 1;
        v = &(SurfReg[kT]->vertices[tvtxN]);
        /* find closest source vertex */
        if(UseHash) svtx = MHTfindClosestVertexNo2(Hash[kS], SurfReg[kS], SurfReg[kT], v, &dmin);
	if(!UseHash || svtx < 0){
	  if(svtx < 0) printf("Target vertex %d of pair %d unmapped in hash, using brute force\n", tvtxN, n);
	  svtx = MRISfindClosestVertex(SurfReg[kS], v->x, v->y, v->z, &dmin, CURRENT_VERTICES);
	}
        tvtxN = svtx;
      }
      /* update the number of hits and distance */
      MRIFseq_vox((SrcHits), svtx, 0, 0, 0)++;
      MRIFseq_vox((TrgHits), tvtx, 0, 0, 0)++;
    }
  }

  /* Set up to create a text file with source-target vertex pairs (STVP) where the source
     is the fist surface and the target is the last surface. The format will be
         srcvtxno srcx srcy srcz trgvtxno trgx trgy trgz 
     The actual coordinates will come from the TMP_VERTEX v->{tx,ty,tz}, 
     so make sure those are set. This functionality is mostly for debugging purposes.  */
  FILE *stvpairfp = NULL;
  std::string stvpairfile = getenv("FS_MRISAPPLYREG_STVPAIR");
  if(stvpairfile.length() > 0){
    stvpairfp = fopen(stvpairfile.c_str(),"w");
    if(stvpairfp == NULL){
      printf("ERROR: could not open stvpairfile %s\n",stvpairfile.c_str());
    }
  }

  /* Go through the forwad loop (finding closest srcvtx to each trgvtx).
  This maps each target vertex to a source vertex */
  printf("MRISapplyReg: Forward Loop (%d)\n", TrgSurfReg->nvertices);
  // nunmapped = 0;
  for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++) {
    if(TrgSurfReg->vertices[tvtx].ripflag) continue;
    if (!UseHash) {
      if (tvtx % 100 == 0) {
        printf("%5d ", tvtx);
        fflush(stdout);
      }
      if (tvtx % 1000 == 999) {
        printf("\n");
        fflush(stdout);
      }
    }

    // Compute the source vertex that corresponds to this target vertex
    tvtxN = tvtx;
    int skip = 0;
    int bf = 0;
    for (n = npairs - 1; n >= 0; n--) {
      kS = 2 * n;
      kT = kS + 1;
      // printf("%5d %5d %d %d %d\n",tvtx,tvtxN,n,kS,kT);
      v = &(SurfReg[kT]->vertices[tvtxN]);
      if(v->ripflag){
	skip = 1;
	break;
      }
      /* find closest source vertex */
      bf = 0;
      if (UseHash) svtx = MHTfindClosestVertexNo2(Hash[kS], SurfReg[kS], SurfReg[kT], v, &dmin);
      if (!UseHash || svtx < 0) {
        if (svtx < 0) {
	  printf("Target vertex %d (%g,%g,%g) of pair %d unmapped in hash, using brute force\n", 
		 tvtxN, v->x, v->y, v->z, n);
	  bf = 1;
	}
        svtx = MRISfindClosestVertex(SurfReg[kS], v->x, v->y, v->z, &dmin, CURRENT_VERTICES);
	if(bf){
	  VERTEX *vs = &(SurfReg[kS]->vertices[svtx]);
	  printf("  Source vertex %d (%g,%g,%g) of pair %d mapped using brute force\n", 
		 svtx, vs->x, vs->y, vs->z, n);
	  fflush(stdout);
	}
      }
      if(SurfReg[kS]->vertices[svtx].ripflag){
	skip = 1;
	break;
      }
      // DNG added these lines on Dec 9, 2020 (without checking it in), 
      // but now can't remember why, so commented them out
      //if(dmin > 2.0){
      //skip = 1;
      //break;
      //}
      tvtxN = svtx;
    }
    if(skip) continue;

    if(!bf && stvpairfp){
      // Good for debugging
      v = &(SurfReg[0]->vertices[svtx]);
      fprintf(stvpairfp,"%d %8.4f %8.4f %8.4f    ",svtx,v->tx, v->ty, v->tz);
      v = &(SurfReg[nsurfs-1]->vertices[tvtx]);
      fprintf(stvpairfp,"%d %8.4f %8.4f %8.4f\n",tvtx,v->tx, v->ty, v->tz);
    }

    if (!DoJac) {
      /* update the number of hits */
      MRIFseq_vox(SrcHits, svtx, 0, 0, 0)++;
      MRIFseq_vox(TrgHits, tvtx, 0, 0, 0)++;
      nhits = 1;
    }
    else
      nhits = MRIgetVoxVal(SrcHits, svtx, 0, 0, 0);

    /* accumulate mapped values for each frame */
    for (f = 0; f < SrcSurfVals->nframes; f++)
      MRIFseq_vox(TrgSurfVals, tvtx, 0, 0, f) += (MRIFseq_vox(SrcSurfVals, svtx, 0, 0, f) / nhits);
  }
  if(stvpairfp) fclose(stvpairfp);

  /*---------------------------------------------------------------
  Go through the reverse loop (finding closest trgvtx to each srcvtx
  unmapped by the forward loop). This assures that each source vertex
  is represented in the map */
  if (ReverseMapFlag) {
    printf("MRISapplyReg: Reverse Loop (%d)\n", SrcSurfReg->nvertices);
    nrevhits = 0;
    for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++) {
      if (MRIFseq_vox((SrcHits), svtx, 0, 0, 0) != 0) continue;
      nrevhits++;

      // Compute the target vertex that corresponds to this source vertex
      svtxN = svtx;
      for (n = 0; n < npairs; n++) {
        kS = 2 * n;
        kT = kS + 1;
        // printf("%5d %5d %d %d %d\n",svtx,svtxN,n,kS,kT);
        v = &(SurfReg[kS]->vertices[svtxN]);
        /* find closest target vertex */
        if (UseHash) tvtx = MHTfindClosestVertexNo2(Hash[kT], SurfReg[kT], SurfReg[kS], v, &dmin);
        if (!UseHash || tvtx < 0) {
          if (tvtx < 0) printf("Source vertex %d of pair %d unmapped in hash, using brute force\n", svtxN, n);
          tvtx = MRISfindClosestVertex(SurfReg[kT], v->x, v->y, v->z, &dmin, CURRENT_VERTICES);
        }
        svtxN = tvtx;
      }

      /* update the number of hits */
      MRIFseq_vox((SrcHits), svtx, 0, 0, 0)++;
      MRIFseq_vox((TrgHits), tvtx, 0, 0, 0)++;
      /* accumulate mapped values for each frame */
      for (f = 0; f < SrcSurfVals->nframes; f++)
        MRIFseq_vox(TrgSurfVals, tvtx, 0, 0, f) += MRIFseq_vox(SrcSurfVals, svtx, 0, 0, f);
    }
    printf("  Reverse Loop had %d hits\n", nrevhits);
  }

  /*---------------------------------------------------------------
  Finally, divide the value at each target vertex by the number
  of source vertices mapping into it */
  if (!DoJac) {
    printf("MRISapplyReg: Dividing by number of hits (%d)\n", TrgSurfReg->nvertices);
    for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++) {
      n = MRIFseq_vox((TrgHits), tvtx, 0, 0, 0);
      if (n > 1) {
        for (f = 0; f < SrcSurfVals->nframes; f++) MRIFseq_vox(TrgSurfVals, tvtx, 0, 0, f) /= n;
      }
    }
  }

  /* Count lost sources */
  nSrcLost = 0;
  for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++) {
    n = MRIFseq_vox((SrcHits), svtx, 0, 0, 0);
    if (n == 0) nSrcLost++;
  }
  printf("MRISapplyReg: nSrcLost = %d\n", nSrcLost);

  MRIfree(&SrcHits);
  MRIfree(&TrgHits);
  if (UseHash)
    for (n = 0; n < nsurfs; n++) MHTfree(&Hash[n]);
  return (TrgSurfVals);
}

/*----------------------------------------------------------------
  MRI *surf2surf_nnfr() - NOTE: use MRISapplyReg instead!

  resample values on one surface to another
  using nearest-neighbor forward/reverse (nnfr) method. If the
  ReverseMapFlag=0, the reverse loop is not implemented. The forward
  loop assures that each target vertex is mapped to a source vertex.
  A source vertex may map to multiple target vertices (fan-out),
  however, there may be some source vertices that do not map to any
  target vertex ("lost vertices"). If the reverse loop is chosen,
  then the lost vertices will be mapped to their closest target
  vertex; this will also generate fan-in.

  Several diagnostic results are returned:
     SrcHits - number of target vertices mapped to by each source
               vertex (measures fan-out).
     SrcDist - average distance of a source vertex to it's mapping
               target vertices.
     TrgHits - number of source vertices mapped to each target
               vertex (measures fan-in).
     TrgDist - average distance of a target vertex to it's mapping
               source vertices.

  See also: surf2surf_nnfr_jac()
  ----------------------------------------------------------------*/
MRI *surf2surf_nnfr(MRI *SrcSurfVals,
                    MRI_SURFACE *SrcSurfReg,
                    MRI_SURFACE *TrgSurfReg,
                    MRI **SrcHits,
                    MRI **SrcDist,
                    MRI **TrgHits,
                    MRI **TrgDist,
                    int ReverseMapFlag,
                    int UseHash)
{
  MRI *TrgSurfVals = NULL;
  int svtx, tvtx, f, n, nrevhits, nSrcLost;
  VERTEX *v;
  MHT *SrcHash, *TrgHash;
  float dmin;
  extern char *ResampleVtxMapFile;
  FILE *fp = NULL;

  /* check dimension consistency */
  if (SrcSurfVals->width != SrcSurfReg->nvertices) {
    fprintf(stderr, "surf2surf_nnfr(): Vals and Reg dimension mismatch\n");
    fprintf(stderr, "nVals = %d, nReg %d\n", SrcSurfVals->width, SrcSurfReg->nvertices);
    return (NULL);
  }

  /* allocate a "volume" to hold the output */
  TrgSurfVals = MRIallocSequence(TrgSurfReg->nvertices, 1, 1, MRI_FLOAT, SrcSurfVals->nframes);
  if (TrgSurfVals == NULL) return (NULL);
  MRIcopyHeader(SrcSurfVals, TrgSurfVals);

  /* number of source vertices mapped to each target vertex */
  *TrgHits = MRIallocSequence(TrgSurfReg->nvertices, 1, 1, MRI_FLOAT, 1);
  if (*TrgHits == NULL) return (NULL);
  MRIcopyHeader(SrcSurfVals, *TrgHits);

  /* Average distance of a target vertex from its source vertices  */
  *TrgDist = MRIallocSequence(TrgSurfReg->nvertices, 1, 1, MRI_FLOAT, 1);
  if (*TrgDist == NULL) return (NULL);
  MRIcopyHeader(SrcSurfVals, *TrgDist);

  /* number of target vertices mapped to by each source vertex */
  *SrcHits = MRIallocSequence(SrcSurfReg->nvertices, 1, 1, MRI_FLOAT, 1);
  if (*SrcHits == NULL) return (NULL);
  MRIcopyHeader(SrcSurfVals, *SrcHits);

  /* Average distance of a source vertex from its target vertices  */
  *SrcDist = MRIallocSequence(SrcSurfReg->nvertices, 1, 1, MRI_FLOAT, 1);
  if (*SrcDist == NULL) return (NULL);
  MRIcopyHeader(SrcSurfVals, *SrcDist);

  /* build hash tables */
  if (UseHash) {
    printf("surf2surf_nnfr: building source hash (res=16).\n");
    SrcHash = MHTcreateVertexTable_Resolution(SrcSurfReg, CURRENT_VERTICES, 16);
  }

  /* Open vertex map file */
  if (ResampleVtxMapFile != NULL) {
    fp = fopen(ResampleVtxMapFile, "w");
    if (fp == NULL) {
      printf("ERROR: could not open %s\n", ResampleVtxMapFile);
      exit(1);
    }
  }

  /*---------------------------------------------------------------
    Go through the forwad loop (finding closest srcvtx to each trgvtx).
    This maps each target vertex to a source vertex */
  printf("Surf2Surf: Forward Loop (%d)\n", TrgSurfReg->nvertices);
  for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++) {
    if (!UseHash) {
      if (tvtx % 100 == 0) {
        printf("%5d ", tvtx);
        fflush(stdout);
      }
      if (tvtx % 1000 == 999) {
        printf("\n");
        fflush(stdout);
      }
    }
    /* find closest source vertex */
    v = &(TrgSurfReg->vertices[tvtx]);
    if (UseHash)
      svtx = MHTfindClosestVertexNo2(SrcHash, SrcSurfReg, TrgSurfReg, v, &dmin);
    else
      svtx = MRISfindClosestVertex(SrcSurfReg, v->x, v->y, v->z, &dmin, CURRENT_VERTICES);

    /* hash table failed, so use brute force */
    if (svtx < 0) svtx = MRISfindClosestVertex(SrcSurfReg, v->x, v->y, v->z, &dmin, CURRENT_VERTICES);

    /* update the number of hits and distance */
    MRIFseq_vox((*SrcHits), svtx, 0, 0, 0)++;
    MRIFseq_vox((*TrgHits), tvtx, 0, 0, 0)++;
    MRIFseq_vox((*SrcDist), svtx, 0, 0, 0) += dmin;
    MRIFseq_vox((*TrgDist), tvtx, 0, 0, 0) += dmin;

    /* accumulate mapped values for each frame */
    for (f = 0; f < SrcSurfVals->nframes; f++)
      MRIFseq_vox(TrgSurfVals, tvtx, 0, 0, f) += MRIFseq_vox(SrcSurfVals, svtx, 0, 0, f);

    if (ResampleVtxMapFile != NULL) {
      fprintf(fp, "%6d  (%6.1f,%6.1f,%6.1f)   ", tvtx, v->x, v->y, v->z);
      v = &(SrcSurfReg->vertices[svtx]);
      fprintf(fp, "%6d  (%6.1f,%6.1f,%6.1f)    %5.4f\n", svtx, v->x, v->y, v->z, dmin);
      fflush(fp);
    }
  }
  printf("\n");
  if (UseHash) MHTfree(&SrcHash);

  if (ResampleVtxMapFile != NULL) fclose(fp);

  /*---------------------------------------------------------------
    Go through the reverse loop (finding closest trgvtx to each srcvtx
    unmapped by the forward loop). This assures that each source vertex
    is represented in the map */
  if (ReverseMapFlag) {
    if (UseHash) {
      MHTfree(&SrcHash);
      printf("surf2surf_nnfr: building target hash (res=16).\n");
      TrgHash = MHTcreateVertexTable_Resolution(TrgSurfReg, CURRENT_VERTICES, 16);
    }
    printf("Surf2Surf: Reverse Loop (%d)\n", SrcSurfReg->nvertices);
    nrevhits = 0;
    for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++) {
      if (MRIFseq_vox((*SrcHits), svtx, 0, 0, 0) == 0) {
        nrevhits++;
        /* find closest target vertex */
        v = &(SrcSurfReg->vertices[svtx]);
        if (UseHash)
          tvtx = MHTfindClosestVertexNo2(TrgHash, TrgSurfReg, SrcSurfReg, v, &dmin);
        else
          tvtx = MRISfindClosestVertex(TrgSurfReg, v->x, v->y, v->z, &dmin, CURRENT_VERTICES);
        /* Hash table failed, so use brute force */
        if (tvtx < 0) tvtx = MRISfindClosestVertex(TrgSurfReg, v->x, v->y, v->z, &dmin, CURRENT_VERTICES);

        /* update the number of hits and distance */
        MRIFseq_vox((*SrcHits), svtx, 0, 0, 0)++;
        MRIFseq_vox((*TrgHits), tvtx, 0, 0, 0)++;
        MRIFseq_vox((*SrcDist), svtx, 0, 0, 0) += dmin;
        MRIFseq_vox((*TrgDist), tvtx, 0, 0, 0) += dmin;
        /* accumulate mapped values for each frame */
        for (f = 0; f < SrcSurfVals->nframes; f++)
          MRIFseq_vox(TrgSurfVals, tvtx, 0, 0, f) += MRIFseq_vox(SrcSurfVals, svtx, 0, 0, f);
      }
    }
    if (UseHash) MHTfree(&TrgHash);
    printf("Reverse Loop had %d hits\n", nrevhits);
  }

  /*---------------------------------------------------------------
    Finally, divide the value at each target vertex by the number
    of source vertices mapping into it */
  printf("Surf2Surf: Dividing by number of hits (%d)\n", TrgSurfReg->nvertices);
  for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++) {
    n = MRIFseq_vox((*TrgHits), tvtx, 0, 0, 0);
    if (n > 1) {
      MRIFseq_vox((*TrgDist), tvtx, 0, 0, 0) /= n; /* average distances */
      for (f = 0; f < SrcSurfVals->nframes; f++) MRIFseq_vox(TrgSurfVals, tvtx, 0, 0, f) /= n;
    }
  }
  /* go through the source loop to average the distance */
  nSrcLost = 0;
  for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++) {
    n = MRIFseq_vox((*SrcHits), svtx, 0, 0, 0);
    if (n == 1) continue;
    if (n > 1)
      MRIFseq_vox((*SrcDist), svtx, 0, 0, 0) /= n;
    else { /* unmapped */
      MRIFseq_vox((*SrcDist), svtx, 0, 0, 0) = -10;
      nSrcLost++;
    }
  }

  printf("INFO: nSrcLost = %d\n", nSrcLost);

  return (TrgSurfVals);
}

/*----------------------------------------------------------------
  surf2surf_nnfr_jac() - NOTE: use MRISapplyReg instead!
  this operates basically the same as
  surf2surf_nnfr with the exception that the sum of values
  across the target vertices should equal that of the source
  (thus it is a kind of jacobian correction). This is good
  for mapping quantities like vertex area, the total of which
  should be preserved across mapping.
  ----------------------------------------------------------------*/
MRI *surf2surf_nnfr_jac(MRI *SrcSurfVals,
                        MRI_SURFACE *SrcSurfReg,
                        MRI_SURFACE *TrgSurfReg,
                        MRI **SrcHits,
                        MRI **SrcDist,
                        MRI **TrgHits,
                        MRI **TrgDist,
                        int ReverseMapFlag,
                        int UseHash)
{
  MRI *TrgSurfVals = NULL;
  int svtx, tvtx, f, n, nrevhits, nSrcLost, nhits;
  // int nunmapped;
  VERTEX *v;
  MHT *SrcHash, *TrgHash;
  float dmin, srcval;

  /* check dimension consistency */
  if (SrcSurfVals->width != SrcSurfReg->nvertices) {
    fprintf(stderr, "surf2surf_nnfr(): Vals and Reg dimension mismatch\n");
    fprintf(stderr, "nVals = %d, nReg %d\n", SrcSurfVals->width, SrcSurfReg->nvertices);
    return (NULL);
  }

  /* allocate a "volume" to hold the output */
  TrgSurfVals = MRIallocSequence(TrgSurfReg->nvertices, 1, 1, MRI_FLOAT, SrcSurfVals->nframes);
  if (TrgSurfVals == NULL) return (NULL);
  MRIcopyHeader(SrcSurfVals, TrgSurfVals);

  /* number of source vertices mapped to each target vertex */
  *TrgHits = MRIallocSequence(TrgSurfReg->nvertices, 1, 1, MRI_FLOAT, 1);
  if (*TrgHits == NULL) return (NULL);

  /* Average distance of a target vertex from its source vertices  */
  *TrgDist = MRIallocSequence(TrgSurfReg->nvertices, 1, 1, MRI_FLOAT, 1);
  if (*TrgDist == NULL) return (NULL);

  /* number of target vertices mapped to by each source vertex */
  *SrcHits = MRIallocSequence(SrcSurfReg->nvertices, 1, 1, MRI_FLOAT, 1);
  if (*SrcHits == NULL) return (NULL);

  /* Average distance of a source vertex from its target vertices  */
  *SrcDist = MRIallocSequence(SrcSurfReg->nvertices, 1, 1, MRI_FLOAT, 1);
  if (*SrcDist == NULL) return (NULL);

  /* build hash tables */
  if (UseHash) {
    printf("surf2surf_nnfr_jac: building source hash (res=16).\n");
    SrcHash = MHTcreateVertexTable_Resolution(SrcSurfReg, CURRENT_VERTICES, 16);
  }

  // First forward loop just counts the number of hits for each src
  printf("Surf2SurfJac: 1st Forward Loop (%d)\n", TrgSurfReg->nvertices);
  // nunmapped = 0;
  for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++) {
    /* find closest source vertex */
    v = &(TrgSurfReg->vertices[tvtx]);
    if (UseHash)
      svtx = MHTfindClosestVertexNo2(SrcHash, SrcSurfReg, TrgSurfReg, v, &dmin);
    else
      svtx = MRISfindClosestVertex(SrcSurfReg, v->x, v->y, v->z, &dmin, CURRENT_VERTICES);
    /* hash table failed, so use brute force */
    if (svtx < 0) svtx = MRISfindClosestVertex(SrcSurfReg, v->x, v->y, v->z, &dmin, CURRENT_VERTICES);

    /* update the number of hits and distance */
    MRIFseq_vox((*SrcHits), svtx, 0, 0, 0)++;  // This is what this loop is for
    MRIFseq_vox((*TrgHits), tvtx, 0, 0, 0)++;
    MRIFseq_vox((*SrcDist), svtx, 0, 0, 0) += dmin;
    MRIFseq_vox((*TrgDist), tvtx, 0, 0, 0) += dmin;
  }

  // Second forward loop accumulates
  printf("Surf2SurfJac: 2nd Forward Loop (%d)\n", TrgSurfReg->nvertices);
  for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++) {
    /* find closest source vertex */
    v = &(TrgSurfReg->vertices[tvtx]);
    if (UseHash)
      svtx = MHTfindClosestVertexNo2(SrcHash, SrcSurfReg, TrgSurfReg, v, &dmin);
    else
      svtx = MRISfindClosestVertex(SrcSurfReg, v->x, v->y, v->z, &dmin, CURRENT_VERTICES);
    /* hash table failed, so use bruce force */
    if (svtx < 0) svtx = MRISfindClosestVertex(SrcSurfReg, v->x, v->y, v->z, &dmin, CURRENT_VERTICES);

    nhits = MRIFseq_vox((*SrcHits), svtx, 0, 0, 0);
    /* Now accumulate mapped values for each frame */
    for (f = 0; f < SrcSurfVals->nframes; f++) {
      srcval = MRIFseq_vox(SrcSurfVals, svtx, 0, 0, f);
      srcval /= nhits;  // divide by number of hits
      MRIFseq_vox(TrgSurfVals, tvtx, 0, 0, f) += srcval;
    }
  }

  /*---------------------------------------------------------------
    Go through the reverse loop (finding closest trgvtx to each srcvtx
    unmapped by the forward loop). This assures that each source vertex
    is represented in the map */
  if (ReverseMapFlag) {
    if (UseHash) {
      MHTfree(&SrcHash);
      printf("surf2surf_nnfr: building target hash (res=16).\n");
      TrgHash = MHTcreateVertexTable_Resolution(TrgSurfReg, CURRENT_VERTICES, 16);
    }
    printf("Surf2SurfJac: Reverse Loop (%d)\n", SrcSurfReg->nvertices);
    nrevhits = 0;
    for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++) {
      if (MRIFseq_vox((*SrcHits), svtx, 0, 0, 0) == 0) {
        nrevhits++;
        /* find closest target vertex */
        v = &(SrcSurfReg->vertices[svtx]);
        if (UseHash)
          tvtx = MHTfindClosestVertexNo2(TrgHash, TrgSurfReg, SrcSurfReg, v, &dmin);
        else
          tvtx = MRISfindClosestVertex(TrgSurfReg, v->x, v->y, v->z, &dmin, CURRENT_VERTICES);
        /* Hash table failed, so use brute force */
        if (tvtx < 0) tvtx = MRISfindClosestVertex(TrgSurfReg, v->x, v->y, v->z, &dmin, CURRENT_VERTICES);
        /* update the number of hits and distance */
        MRIFseq_vox((*SrcHits), svtx, 0, 0, 0)++;
        MRIFseq_vox((*TrgHits), tvtx, 0, 0, 0)++;
        MRIFseq_vox((*SrcDist), svtx, 0, 0, 0) += dmin;
        MRIFseq_vox((*TrgDist), tvtx, 0, 0, 0) += dmin;
        /* accumulate mapped values for each frame */
        for (f = 0; f < SrcSurfVals->nframes; f++)
          MRIFseq_vox(TrgSurfVals, tvtx, 0, 0, f) += MRIFseq_vox(SrcSurfVals, svtx, 0, 0, f);
      }
    }
    if (UseHash) MHTfree(&TrgHash);
    printf("Reverse Loop had %d hits\n", nrevhits);
  }

  // Do NOT normalize target vertices with multiple src vertices

  /* go through the source loop to average the distance */
  nSrcLost = 0;
  for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++) {
    n = MRIFseq_vox((*SrcHits), svtx, 0, 0, 0);
    if (n == 1) continue;
    if (n > 1)
      MRIFseq_vox((*SrcDist), svtx, 0, 0, 0) /= n;
    else { /* unmapped */
      MRIFseq_vox((*SrcDist), svtx, 0, 0, 0) = -10;
      nSrcLost++;
    }
  }
  printf("INFO: nSrcLost = %d\n", nSrcLost);
  printf("surf2surf_nnfr_jac() done\n");
  fflush(stdout);

  return (TrgSurfVals);
}

/*-------------------------------------------------------------
  crs2ind() -- returns linear index into a volume stored by column,
  row, slice.
  --------------------------------------------------------------------*/
int crs2ind(int *ind, int c, int r, int s, int ncols, int nrows, int nslcs)
{
  if (c < 0 || c >= ncols) {
    fprintf(stderr, "crs2ind: col %d out of bounds (0,%d)\n", c, ncols - 1);
    return (1);
  }
  if (r < 0 || r >= nrows) {
    fprintf(stderr, "crs2ind: row %d out of bounds (0,%d)\n", r, nrows - 1);
    return (1);
  }
  if (s < 0 || s >= nslcs) {
    fprintf(stderr, "crs2ind: slc %d out of bounds (0,%d)\n", c, nslcs - 1);
    return (1);
  }

  *ind = c + r * ncols + s * nrows * ncols;

  return (0);
}
/*-------------------------------------------------------------
  ind2crs() -- returns the column, row, and slice of an element in a
  volume given the index of the element in the volume assuming that
  the elements are stored by column, row, then slice.
  --------------------------------------------------------------------*/
int ind2crs(int *c, int *r, int *s, int ind, int ncols, int nrows, int nslcs)
{
  int i = ind, ntot, nrowcols;

  nrowcols = nrows * ncols;
  ntot = nrowcols * nslcs;
  if (i < 0 || i >= ntot) {
    fprintf(stderr, "ind2crs: index %d out of bounds (0,%d)\n", i, ntot - 1);
    return (1);
  }

  *s = (int)(i / nrowcols);
  i = i - *s * (nrowcols);

  *r = (int)(i / ncols);
  i = i - *r * ncols;

  *c = (int)i;

  return (0);
}

/*
  \fn MRI *MRIsurf2VolOpt(MRI *ribbon, MRIS **surfs, MRI **overlays,
                    int nsurfs, LTA *Q, MRI *volsurf)

  \brief Maps surface vertex values to the volume filling the entire
  ribbon. Loops through all the voxels in the volume but only fills
  values that are within the ribbon (3,42). Multiple (nsurfs) surfaces
  can be passed; this number must match the number of overlays.  For a
  given ribbon voxel, the nearest vertex to each of the nsurf surfs is
  determined. The voxel gets the value of this vertex in the overlay
  corresponding to the closest surface. Ribbon is assumed to be in the
  conformed space. If Q is null, the output volume matches the
  ribbon. If Q is non-null, then the output volume matches either the
  src or dst depending on which one is NOT the ribbon (ie, the
  direction of Q does not matter). This function differs from
  MRIsurf2Vol() in that this one will always fill the ribbon (by
  construction).

  This function can be used in several ways. If both lh and rh are to
  be filled, then surfs[0] = lh, surfs[1] = rh, overlays[0] =
  lhoverlay, overlays[1] = rhoverlay. Alternatively, if there are
  multiple overlay layers, eg, at white and pial, then surfs[0] =
  lh.white, surfs[1] = lh.pial, overlays[0] = lh.white.overlay,
  overlays[1] = lh.pial.overlay.

 */
MRI *MRIsurf2VolOpt(MRI *ribbon, MRIS **surfs, MRI **overlays, int nsurfs, LTA *Q, MRI *volsurf)
{
  int n, c, r, s, f, nmin, vtxno, vtxnomin = 0, nframes, ribval, cR, rR, sR;
  MHT **hash = NULL;
  int UseHash = 1;
  MATRIX *T, *invR, *M, *surfRAS = NULL, *crs, *R, *crsRibbon;

  float dmin, d, val;
  LTA *V = NULL, *Q2;

  // Make sure that each overlay has the same number of frames
  for (n = 0; n < nsurfs; n++) {
    if (overlays[n]->nframes != overlays[0]->nframes) {
      printf("ERROR: MRIsurf2VolOpt(): overlay dim mismatch %d\n", n);
      return (NULL);
    }
  }

  if (Q && !LTAmriIsSource(Q, ribbon) && !LTAmriIsTarget(Q, ribbon)) {
    printf("ERROR: MRIsurf2VolOpt(): ribbon is neither source nor target\n");
    return (NULL);
  }

  if (Q == NULL)
    Q2 = TransformRegDat2LTA(ribbon, ribbon, NULL);  // regdat unrelated, just a way to get LTA
  else {
    if (LTAmriIsSource(Q, ribbon)) {
      printf("MRIsurf2VolOpt(): ribbon is source\n");
      Q2 = LTAcopy(Q, NULL);
    }
    else {
      printf("MRIsurf2VolOpt(): inverting LTA\n");
      Q2 = LTAinvert(Q, NULL);
    }
  }
  // Now Q2 goes from output volume to ribbon

  nframes = overlays[0]->nframes;
  if (volsurf == NULL) {
    volsurf = MRIallocFromVolGeom(&(Q2->xforms[0].dst), MRI_FLOAT, nframes, 0);
    if (volsurf == NULL) {
      printf("ERROR: MRIsurf2VolOpt(): could not alloc\n");
      return (NULL);
    }
    MRIcopyPulseParameters(overlays[0], volsurf);
  }
  // else should make sure that Q2->xforms[0]->src is consistent with volsurf

  // Get outputCRS to ribbonCRS
  V = LTAcopy(Q2, NULL);
  V = LTAchangeType(V, LINEAR_VOX_TO_VOX);
  LTAfillInverse(V);

  // M converts output CRS to surface RAS
  R = TransformLTA2RegDat(Q2);  // As a regdat, it maps from ribbonRAS to outRAS
  R = MatrixInverse(R, NULL);
  T = MRIxfmCRS2XYZtkreg(volsurf);
  invR = MatrixInverse(R, NULL);
  M = MatrixMultiply(invR, T, NULL);

  if (UseHash) {
    hash = (MHT **)calloc(sizeof(MHT *), nsurfs);
    for (n = 0; n < nsurfs; n++) {
      hash[n] = MHTcreateVertexTable_Resolution(surfs[n], CURRENT_VERTICES, 16);
    }
  }

  crs = MatrixAlloc(4, 1, MATRIX_REAL);
  crs->rptr[4][1] = 1;
  crsRibbon = MatrixAlloc(4, 1, MATRIX_REAL);
  surfRAS = MatrixAlloc(4, 1, MATRIX_REAL);
  for (c = 0; c < volsurf->width; c++) {
    for (r = 0; r < volsurf->height; r++) {
      for (s = 0; s < volsurf->depth; s++) {
        // Make sure the output voxel is in the ribbon
        crs->rptr[1][1] = c;
        crs->rptr[2][1] = r;
        crs->rptr[3][1] = s;
        crsRibbon = MatrixMultiply(V->inv_xforms[0].m_L, crs, crsRibbon);
        cR = nint(crsRibbon->rptr[1][1]);
        rR = nint(crsRibbon->rptr[1][2]);
        sR = nint(crsRibbon->rptr[1][3]);
        if (cR < 0 || cR >= ribbon->width) continue;
        if (rR < 0 || rR >= ribbon->height) continue;
        if (sR < 0 || sR >= ribbon->depth) continue;
        ribval = MRIgetVoxVal(ribbon, cR, rR, sR, 0);
        if (ribval != 3 && ribval != 42) {
          for (f = 0; f < nframes; f++) MRIsetVoxVal(volsurf, c, r, s, f, 0);
          continue;
        }
        // Compute the surface location of this point
        surfRAS = MatrixMultiply(M, crs, surfRAS);
        float const x = surfRAS->rptr[1][1];
        float const y = surfRAS->rptr[2][1];
        float const z = surfRAS->rptr[3][1];
        // Find surface vertex with closest location
        dmin = 1000;
        nmin = -1;
        d = 999;
        vtxnomin = -1;
        for (n = 0; n < nsurfs; n++) {
          if (surfs[n]->hemisphere == LEFT_HEMISPHERE && ribval != 3) continue;
          if (surfs[n]->hemisphere == RIGHT_HEMISPHERE && ribval != 42) continue;

          if (UseHash)
            vtxno = MHTfindClosestVertexNoXYZ(hash[n], surfs[n], x,y,z, &d);
          else
            vtxno = MRISfindClosestVertex(surfs[n], x, y, z, &d, CURRENT_VERTICES);
            
          if (vtxno < 0) {
            printf("ERROR: MRIsurf2VolOpt(): No Match: %3d %3d %3d    %6.2f %6.2f %6.2f\n", c, r, s, x, y, z);
            vtxno = MRISfindClosestVertex(surfs[n], x, y, z, &d, CURRENT_VERTICES);
            printf("%d    %d %d %d     %g %g %g  %5d\n", n, c, r, s, x, y, z, vtxno);
            printf("V -------------------------\n");
            MatrixPrint(stdout, V->inv_xforms[0].m_L);
            fflush(stdout);
            printf("T = [\n");
            MatrixPrint(stdout, T);
            printf("R = [\n");
            MatrixPrint(stdout, R);
            printf("M = [\n");
            MatrixPrint(stdout, M);
            fflush(stdout);
            return (NULL);
          }
          if (d < dmin) {
            dmin = d;
            nmin = n;  // index of surface with closest vertex
            vtxnomin = vtxno;
          }
        }  // surfs
        // This can happen if only one hemi is specified but the voxel is in the other
        if (nmin == -1) continue;
        // Assign value from vertex to voxel
        for (f = 0; f < nframes; f++) {
          val = MRIgetVoxVal(overlays[nmin], vtxnomin, 0, 0, f);
          MRIsetVoxVal(volsurf, c, r, s, f, val);
        }
      }  // slice
    }    // row
  }      // col

  if (UseHash)
    for (n = 0; n < nsurfs; n++)
      if (UseHash) MHTfree(&hash[n]);

  fflush(stdout);
  MatrixFree(&T);
  MatrixFree(&invR);
  MatrixFree(&R);
  MatrixFree(&surfRAS);
  MatrixFree(&M);
  LTAfree(&V);
  LTAfree(&Q2);
  return (volsurf);
}

/*-------------------------------------------------------------------
  \fn int MRIsurf2Vol(MRI *surfvals, MRI *vol, MRI *map)
  \brief the purpose of this function is to resample values
  from a surface to a volume. Most of the real work is done in the
  creation of map, which has the same size as vol; each voxel in map
  has the vertex number of the corresponding point on the surface. If
  there is no corresponding point, then the value will be -1. See
  MRImapSurf2VolClosest().  MRI vol must have already been allocated.
  Returns -1 on error, otherwise returns the number of voxels in the
  volume that have been hit by the surface. This function differs from
  MRIsurf2VolOpt() in that Opt will always fill the ribbon (by
  construction). Filling by projection can leave holes.
  ----------------------------------------------------------------*/
int MRIsurf2Vol(MRI *surfvals, MRI *vol, MRI *map)
{
  int vtx, c, r, s, f, nhits;
  float val;

  if (vol->width != map->width || vol->height != map->height || vol->depth != map->depth) {
    printf("ERROR: vol and map dimensions differ\n");
    return (-1);
  }

  if (surfvals->nframes != vol->nframes) {
    printf("ERROR: surfvals and vol have different number of frames\n");
    return (-1);
  }

  nhits = 0;
  for (c = 0; c < vol->width; c++) {
    for (r = 0; r < vol->height; r++) {
      for (s = 0; s < vol->depth; s++) {
        vtx = MRIIseq_vox(map, c, r, s, 0);
        for (f = 0; f < vol->nframes; f++) {
          if (vtx < 0) {
#if 0
            MRIsetVoxVal(vol,c,r,s,f,0.0);
#endif
            continue;
          }
          val = MRIgetVoxVal(surfvals, vtx, 0, 0, f);
          MRIsetVoxVal(vol, c, r, s, f, val);
          nhits++;
        }
      }
    }
  }

  /*
  WARNING! nhits will not be correct in the case where
  this function is called multiple times (e.g. for filling
  the cortical ribbon) as it won't take into account the fact that a voxel
  is mapped more than one times (which will happen a lot).
  */
  return (nhits);
}
/*-------------------------------------------------------------------
  MRImapSurf2VolClosest() - the purpose of this function is to create
  a map of a volume in which the value at each voxel is the integer
  value of the vertex that is closest to it. If no vertex is within a
  voxel, its value is set to -1.

  Qa2v is a matrix that maps XYZ in anatomical/surface space to CRS in
  volume space. This matrix can be computed as follows:

          Qa2v = inv(Tvol)*Ma2v,

  where Tvol maps from CRS of the vol to XYZ of the volume, and Ma2v
  maps from XYZ of the reference anatomical to XYZ of the volume.

  If Ma2v is a tkregister-compatible matrix, then

          Tvol = MRIxfmCRS2XYZtkreg(vol);

  If Ma2v maps between the native (ie, scanner) coordinates , then

          Tvol = MRIxfmCRS2XYZ(vol,0);

  If projfrac is non-zero, the XYZ of a vertex is computed based on
  its native XYZ offset in the direction of the surface normal by
  projfrac fraction of the thickness at that point. Obviously, the
  thickness must have been loaded into the surface at this point.
  ------------------------------------------------------------------*/
MRI *MRImapSurf2VolClosest(MRIS *surf, MRI *vol, MATRIX *Qa2v, float projfrac)
{
  MRI *map, *dist2;
  int vtx, vtxmin, c, r, s;
  float xvtx, yvtx, zvtx;
  float d2, d2min;
  MATRIX *xyzvtx, *icrs, *fcrs, *xyzcrs;
  // MATRIX *Tvol;

  /* Alloc map - for a given voxel, holds the number of the
     closest vertex */
  map = MRIalloc(vol->width, vol->height, vol->depth, MRI_INT);
  if (map == NULL) {
    printf("ERROR: MRImapSurf2VolClosest: could not alloc vtx map\n");
    return (NULL);
  }
  /* Alloc dist - for a given voxel, holds the distance of the
     closest vertex to the voxel. */
  dist2 = MRIalloc(vol->width, vol->height, vol->depth, MRI_FLOAT);
  if (dist2 == NULL) {
    printf("ERROR: MRImapSurf2VolClosest: could not alloc dist map\n");
    return (NULL);
  }

  /* Initially set all the voxels in the map to -1 to mark that
     they have not been hit by a surface vertex */
  MRIvalueFill(map, -1);

  /* Intialize matrix stuff */
  // Tvol = MRIxfmCRS2XYZ(vol, 0); /* converts crs to xyz in vol */
  xyzvtx = MatrixAlloc(4, 1, MATRIX_REAL);
  fcrs = MatrixAlloc(4, 1, MATRIX_REAL);   /* vertex row, col, slice */
  icrs = MatrixAlloc(4, 1, MATRIX_REAL);   /* voxel row, col, slice */
  xyzcrs = MatrixAlloc(4, 1, MATRIX_REAL); /* voxel xyz */

  /* Set the 4th item to 1 */
  xyzvtx->rptr[4][1] = 1.0;
  fcrs->rptr[4][1] = 1.0;
  icrs->rptr[4][1] = 1.0;
  xyzcrs->rptr[4][1] = 1.0;

  /*----------- Loop through vertices --------------*/
  for (vtx = 0; vtx < surf->nvertices; vtx++) {
    if (projfrac == 0) {
      xvtx = surf->vertices[vtx].x;
      yvtx = surf->vertices[vtx].y;
      zvtx = surf->vertices[vtx].z;
    }
    else {
      /* Get the xyz of the vertex as projected along the normal a
      distance equal to a fraction of the cortical thickness at
      that point. */
      ProjNormFracThick(&xvtx, &yvtx, &zvtx, surf, vtx, projfrac);
    }

    /* load vertex xyz values into a vector */
    xyzvtx->rptr[1][1] = xvtx;
    xyzvtx->rptr[2][1] = yvtx;
    xyzvtx->rptr[3][1] = zvtx;
    /* [4] is already 1 */

    /* Compute the CRS in volume space */
    /* fcrs = Qa2v * xyzvtx */
    MatrixMultiply(Qa2v, xyzvtx, fcrs);

    /* Round CRS to nearest integer */
    c = nint(fcrs->rptr[1][1]);
    r = nint(fcrs->rptr[2][1]);
    s = nint(fcrs->rptr[3][1]);

    if (Gdiag_no > 0 && vtx == 30756) {
      printf("diag -----------------------------\n");
      printf("vtx = %d  %g %g %g\n", vtx, xvtx, yvtx, zvtx);
      printf("fCRS  %g %g %g\n", fcrs->rptr[1][1], fcrs->rptr[2][1], fcrs->rptr[3][1]);
      printf("CRS  %d %d %d\n", c, r, s);
    }

    /* Check that it is in the volume */
    if (c < 0 || c >= vol->width) continue;
    if (r < 0 || r >= vol->height) continue;
    if (s < 0 || s >= vol->depth) continue;

    /* Load rounded CRS into vector */
    icrs->rptr[1][1] = c;
    icrs->rptr[2][1] = r;
    icrs->rptr[3][1] = s;
    /* [4] is already 1 */

    /* Compute XYZ of rounded CRS in volume space. This
    is the XYZ of the center of the voxel*/
    /* xyzcrs = Tvol * icrs */
    // MatrixMultiply(Tvol,icrs,xyzcrs);

    /* Compute the distance between the voxel and the
       vertex (actually distance squared) */
    d2 = SQR(vol->xsize * (c - icrs->rptr[1][1])) + SQR(vol->ysize * (r - icrs->rptr[2][1])) +
         SQR(vol->zsize * (s - icrs->rptr[3][1]));

    /* Check whether this voxel has been hit. If not
       just set its map vertex and dist to current. */
    vtxmin = MRIIseq_vox(map, c, r, s, 0);
    if (vtxmin < 0) {
      MRIIseq_vox(map, c, r, s, 0) = vtx;
      MRIFseq_vox(dist2, c, r, s, 0) = d2;
      continue;
    }

    /* Check whether this vertex is closer */
    d2min = MRIFseq_vox(dist2, c, r, s, 0);
    if (d2min > d2) {
      MRIIseq_vox(map, c, r, s, 0) = vtx;
      MRIFseq_vox(dist2, c, r, s, 0) = d2;
    }

  } /* end loop over vertices */

  // MatrixFree(&Tvol);
  MatrixFree(&xyzvtx);
  MatrixFree(&icrs);
  MatrixFree(&fcrs);
  MatrixFree(&xyzcrs);
  MRIfree(&dist2);
  return (map);
}
/*
  \fn MRI *MRIseg2SegPVF(MRI *seg, LTA *seg2vol, double resmm, int *segidlist, int nsegs, MRI *mask, int ReInit, MRI
  *out)
  \brief Computes the partial volume fraction (PVF:0->1) for each
  segmentation. The output geometry is derived from the dst in seg2vol
  (seg2vol can go the other direction, it figures it out). The number
  of frames in the output equals the number of segmentations
  (excluding 0).  Each frame is a map where the voxel value is the
  fraction of that voxel that contains the seg corresonding to that
  frame. The algorithm works by sampling each output voxel over a
  uniform grid. For each sample, the location in seg is computed to
  get the segmentation. The number of hits for this segmentation/frame
  is then updated for that voxel. The distance in mm between the
  points of the grid is controlled by resmm. The actual distance will
  be changed so as to assure that the entire output volume is
  uniformly sampled. If resmm is negative, then it is interpreted as
  an upsampling factor (USF=round(abs(resmm))).  If resmm==0, then it
  is set to the minimum voxel size in the segmentation, which is
  probably reasonable.  segidlist is a complete list of segmentation
  ids (excluding 0) and nsegs is the number (equal to the number of
  frames in the output).  This can be obtained with segidlist =
  MRIsegIdListNot0(seg, &nsegs, 0).  The sum of PVFs across seg will
  not always equal 1 because the background is not treated as a
  segmentation.  MRIseg2SegPVF() is built for speed and has several
  optimizations, one of which is that some data are precomputed and
  cached. If the mask or segmentation changes, then call with
  ReInit=1. Run with ReInit=-1 to free the memory in the cache.  This
  function uses OMP.  See also MRIsegPVF2Seg(),
  MRIsegPVF2TissueTypePVF(), MRIseg2TissueType().  If ct is non-NULL,
  then the segmentation is computed based on maximum PVF.  See
  VOXsegPVF2Seg().  There is no specific advantage to upsampling or
  masking seg.  mask is a mask in the output space; anything outside
  of the mask is set to 0.
*/
MRI *MRIseg2SegPVF(
    MRI *seg, LTA *seg2vol, double resmm, int *segidlist, int nsegs, MRI *mask, int ReInit, COLOR_TABLE *ct, MRI *out)
{
  LTA *lta;
  VOL_GEOM *vg;
  double m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34;
  double dc, dr, ds, dcstart, drstart, dsstart;
  float **nperfth, npervox;  // need float for nperf for consistency
  int nhits, segidmax, nthvox, ndc, ndr, nds, nthreads;
  int c, r, s, f, nframesout = 0, outtype = 0, DoSeg;
  static int InitNeeded = 1, nvox = 0, *clist = NULL, *rlist = NULL, *slist = NULL, *seg2frame, nthcall = 0;

  if (ReInit < 0) {
    // This allows the static memory to be freed
    if (clist) free(clist);
    if (rlist) free(rlist);
    if (slist) free(slist);
    if (seg2frame) free(seg2frame);
    InitNeeded = 1;
    nthcall = 0;
    return (NULL);
  }

  Timer timer;

  if (ct)
    DoSeg = 1;
  else
    DoSeg = 0;

  // Get number of threads
  nthreads = 1;
#ifdef HAVE_OPENMP
  nthreads = omp_get_max_threads();  // using max should be ok
#endif

  if (LTAmriIsSource(seg2vol, seg))
    lta = LTAinvert(seg2vol, NULL);
  else
    lta = LTAcopy(seg2vol, NULL);
  if (lta->type != LINEAR_VOX_TO_VOX) LTAchangeType(lta, LINEAR_VOX_TO_VOX);
  vg = &(lta->xforms[0].src);

  if (resmm == 0) resmm = MIN(MIN(seg->xsize, seg->ysize), seg->zsize);

  if (resmm > 0) {
    // compute number of subsamples based on passed resolution
    ndc = round(vg->xsize / resmm);
    if (ndc < 1) ndc = 1;
    ndr = round(vg->ysize / resmm);
    if (ndr < 1) ndr = 1;
    nds = round(vg->zsize / resmm);
    if (nds < 1) nds = 1;
  }
  // otherwise, use the value in resmm as the upsample factor
  else
    ndc = ndr = nds = round(abs(resmm));

  // this assures that the volume is sampled unformly in each dim
  dc = 1.0 / ndc;
  dcstart = -0.5 + dc / 2.0;
  dr = 1.0 / ndr;
  drstart = -0.5 + dr / 2.0;
  ds = 1.0 / nds;
  dsstart = -0.5 + ds / 2.0;
  npervox = ndc * ndr * nds;  // number of samples in each voxel

  /* For speed, pre-compute some variables. If the mask or the
     segs that compose the segmentation change, then the cache
     should be re-computed by passing ReInit=1. The caller is
     responsible for determining whether the cache needs to be
     re-initialized. */
  if (InitNeeded) ReInit = 1;
  if (ReInit) {
    /* Create an array to map from segid to frame (for speed). The array
       may be very large depending on segidmax, but still it should be
       small relative to images. If the segs that compose the segmentation
       change, the cache should be re-initialized.*/
    if (seg2frame) free(seg2frame);
    if (clist) free(clist);
    if (rlist) free(rlist);
    if (slist) free(slist);

    segidmax = 0;
    for (f = 0; f < nsegs; f++)
      if (segidmax < segidlist[f]) segidmax = segidlist[f];
    seg2frame = (int *)calloc(sizeof(int), segidmax + 1);
    for (f = 0; f < nsegs; f++) seg2frame[segidlist[f]] = f;

    /* Precompute the number of voxels in the mask. If the mask
       changes, then the cache should be re-initialized.*/
    if (mask)
      nvox = MRIcountAboveThreshold(mask, 0.5);
    else
      nvox = vg->width * vg->height * vg->depth;
    printf("MRIseg2SegPVF(): Initilizing cache nsegs = %d, nvox = %d, nthreads=%d, DoSeg=%d\n",
           nsegs,
           nvox,
           nthreads,
           DoSeg);
    printf("resmm=%g c: %g %d %g %g, r: %g %d %g %g, s: %g %d %g %g\n",
           resmm,
           vg->xsize,
           ndc,
           dcstart,
           dc,
           vg->ysize,
           ndr,
           drstart,
           dr,
           vg->zsize,
           nds,
           dsstart,
           ds);
    printf("delta %g %g %g, npervox = %g\n", vg->xsize / ndc, vg->ysize / ndr, vg->zsize / nds, npervox);
    fflush(stdout);
    /* Precompute arrays that map voxel index to col,row,slice in
       output. If the mask changes, then the cache should be
       re-initialized.*/
    clist = (int *)calloc(sizeof(int), nvox);
    rlist = (int *)calloc(sizeof(int), nvox);
    slist = (int *)calloc(sizeof(int), nvox);
    nthvox = 0;
    for (c = 0; c < vg->width; c++) {
      for (r = 0; r < vg->height; r++) {
        for (s = 0; s < vg->depth; s++) {
          if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
          clist[nthvox] = c;
          rlist[nthvox] = r;
          slist[nthvox] = s;
          nthvox++;
        }
      }
    }
    InitNeeded = 0;
  }

  // Alloc output and check dimensions
  if (DoSeg) {
    nframesout = 1;
    outtype = MRI_INT;
  }
  else {
    nframesout = nsegs;
    outtype = MRI_FLOAT;
  }
  if (out == NULL) {
    out = MRIallocSequence(vg->width, vg->height, vg->depth, outtype, nframesout);
    if (out == NULL) {
      printf("ERROR: MRIseg2SegPVF(): could not alloc\n");
      return (NULL);
    }
    useVolGeomToMRI(vg, out);
    MRIcopyPulseParameters(seg, out);
  }
  if (out->width != vg->width || out->height != vg->height || out->depth != vg->depth || out->nframes != nframesout) {
    printf("ERROR: MRIseg2SegPVF(): dimension mismatch\n");
    return (NULL);
  }
  if (mask && MRIdimMismatch(out, mask, 0)) {
    printf("ERROR: MRIseg2SegPVF() output/mask dim mismatch\n");
    return (NULL);
  }

  // For speed, put matrix values in their own variables
  m11 = lta->xforms[0].m_L->rptr[1][1];
  m12 = lta->xforms[0].m_L->rptr[1][2];
  m13 = lta->xforms[0].m_L->rptr[1][3];
  m14 = lta->xforms[0].m_L->rptr[1][4];
  m21 = lta->xforms[0].m_L->rptr[2][1];
  m22 = lta->xforms[0].m_L->rptr[2][2];
  m23 = lta->xforms[0].m_L->rptr[2][3];
  m24 = lta->xforms[0].m_L->rptr[2][4];
  m31 = lta->xforms[0].m_L->rptr[3][1];
  m32 = lta->xforms[0].m_L->rptr[3][2];
  m33 = lta->xforms[0].m_L->rptr[3][3];
  m34 = lta->xforms[0].m_L->rptr[3][4];

  if (Gdiag_no > 0) {
    printf("resmm=%g c: %g %d %g %g, r: %g %d %g %g, s: %g %d %g %g\n",
           resmm,
           vg->xsize,
           ndc,
           dcstart,
           dc,
           vg->ysize,
           ndr,
           drstart,
           dr,
           vg->zsize,
           nds,
           dsstart,
           ds);
    printf("delta %g %g %g, npervox = %g\n", vg->xsize / ndc, vg->ysize / ndr, vg->zsize / nds, npervox);
    printf("MRIseg2SegPVF(): starting fill t = %g, nvox=%d\n", timer.seconds(), nvox);
    fflush(stdout);
  }

  // Set up array to hold the number per frame in each thread (n per f th)
  nperfth = (float **)calloc(sizeof(float *), nthreads);
  for (f = 0; f < nthreads; f++) nperfth[f] = (float *)calloc(sizeof(float), nsegs);

  /* The main loop goes over each voxel in the output/mask. This is
     thread-safe because each voxel is handled separately. */
  nhits = 0;  // keep track of the total number of hits
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  // note: removing reduction(+:nhits) slows the speed to that of 1 thread
  #pragma omp parallel for if_ROMP(experimental) shared(nperfth, m13, m23, m33) reduction(+ : nhits)
#endif
  for (nthvox = 0; nthvox < nvox; nthvox++) {
    ROMP_PFLB_begin
    
    int c, r, s, i, j, k, ca, ra, sa, segid, f, threadno;
    double cf, rf, sf;
    double m11cf, m21cf, m31cf, m12rf, m22rf, m32rf;
    float *nperf;

    // Get the thread number
    threadno = 0;
#ifdef HAVE_OPENMP
    threadno = omp_get_thread_num();
#endif
    nperf = nperfth[threadno];

    for (f = 0; f < nsegs; f++) nperf[f] = 0;  // zero n-per-frame

    // Get the col, row, slice of the current voxel
    c = clist[nthvox];
    r = rlist[nthvox];
    s = slist[nthvox];
    // Uniformly sample ndc-by-ndr-by-nds points in each voxel
    for (i = 0; i < ndc; i++) {
      cf = c + dcstart + dc * i;
      m11cf = m11 * cf + m14;
      m21cf = m21 * cf + m24;
      m31cf = m31 * cf + m34;
      for (j = 0; j < ndr; j++) {
        rf = r + drstart + dr * j;
        m12rf = m12 * rf + m11cf;
        m22rf = m22 * rf + m21cf;
        m32rf = m32 * rf + m31cf;
        for (k = 0; k < nds; k++) {
          sf = s + dsstart + ds * k;

          // ca = nint(m11*cf + m12*rf + m13*sf + m14);
          ca = nint(m12rf + m13 * sf);
          if (ca < 0 || ca >= seg->width) continue;

          // ra = nint(m21*cf + m22*rf + m23*sf + m24);
          ra = nint(m22rf + m23 * sf);
          if (ra < 0 || ra >= seg->height) continue;

          // sa = nint(m31*cf + m32*rf + m33*sf + m34);
          sa = nint(m32rf + m33 * sf);
          if (sa < 0 || sa >= seg->depth) continue;

          segid = MRIgetVoxVal(seg, ca, ra, sa, 0);
          if (segid == 0) continue;

          f = seg2frame[segid];  // convert the segid to a frame
          nperf[f]++;            // increment the number of hits for this frame
          nhits++;               // total number of hits
        }                        // ds
      }                          // dr
    }                            // dc
    // Done with this voxel
    // Divide by npervox to compute PVF
    if (DoSeg == 0) {
      for (f = 0; f < nsegs; f++)  // dont compute pvf outside, becomes very slow
        if (nperf[f] > 0) MRIsetVoxVal(out, c, r, s, f, (float)nperf[f] / npervox);
    }
    else {
      // Get seg with highest PVF
      for (f = 0; f < nsegs; f++) nperf[f] = (float)nperf[f] / npervox;
      segid = VOXsegPVF2Seg(nperf, segidlist, nsegs, ct);
      MRIsetVoxVal(out, c, r, s, 0, segid);
    }
    
    ROMP_PFLB_end
  }  // nthvox
  ROMP_PF_end
  
  for (f = 0; f < nthreads; f++) free(nperfth[f]);
  free(nperfth);

  if (nthcall < 4 || Gdiag_no > 0) {
    printf("MRIseg2SegPVF(): done t = %g, nhits=%d, nthreads=%d\n", timer.seconds(), nhits, nthreads);
    fflush(stdout);
  }
  nthcall++;
  return (out);
}
/*
  \fn MRI *MRIsegPVF2Seg(MRI *segpvf, int *segidlist, int nsegs, COLOR_TABLE *ct, MRI *mask, MRI *seg)
  \brief Converts a segmentation-wise PVF created by MRIseg2SegPVF()
  into a segmentation by selecting the segmentation that has the
  largest PVF. If the sum of PVFs at a voxel is less than 0.5, then
  the seg is set to 0.  If no seg has PVF>.5 and ct is non-NULL, then
  the PVFs within each tissue type are computed, and the seg with the
  max PVF within the tissue type with the max PVF is selected. See
  also VoxsegPVF2Seg().
 */
MRI *MRIsegPVF2Seg(MRI *segpvf, int *segidlist, int nsegs, COLOR_TABLE *ct, MRI *mask, MRI *seg)
{
  int c, r, s, f, segid;
  float *vlist;

  if (seg == NULL) {
    seg = MRIallocSequence(segpvf->width, segpvf->height, segpvf->depth, MRI_INT, 1);
    MRIcopyHeader(segpvf, seg);
    MRIcopyPulseParameters(segpvf, seg);
  }
  else
    MRIclear(seg);

  vlist = (float *)calloc(sizeof(float), nsegs);

  for (c = 0; c < segpvf->width; c++) {
    for (r = 0; r < segpvf->height; r++) {
      for (s = 0; s < segpvf->depth; s++) {
        MRIsetVoxVal(seg, c, r, s, 0, 0);
        if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
        // Go through each frame/seg
        for (f = 0; f < nsegs; f++) vlist[f] = MRIgetVoxVal(segpvf, c, r, s, f);
        segid = VOXsegPVF2Seg(vlist, segidlist, nsegs, ct);
        MRIsetVoxVal(seg, c, r, s, 0, segid);
      }
    }
  }
  free(vlist);
  return (seg);
}

/*
  \fn int VOXsegPVF2Seg(double *segpvfvox, int *segidlist, int nsegs, COLOR_TABLE *ct)
  \brief Selects the seg with the greatest PVF.  If the sum of PVFs at
  a voxel is less than 0.5, then the seg is set to 0.  If no seg has
  PVF>.5 and ct is non-NULL, then the PVFs within each tissue type are
  computed, and the seg with the max PVF within the tissue type with
  the max PVF is selected. segpvfvox is a vector of length nsegs with
  the PVF for each of the segs. segidlist is the list of segmentation
  ids (also length nsegs). See also MRIsegPVF2Seg().
*/
int VOXsegPVF2Seg(float *segpvfvox, int *segidlist, int nsegs, COLOR_TABLE *ct)
{
  int segid, f, fmax, tt, ttmax, nTT;
  float v, vsum, vmax, vtt[100], vttmax;

  vsum = 0;
  vmax = 0;
  fmax = 0;
  for (f = 0; f < nsegs; f++) {
    v = segpvfvox[f];
    vsum += v;
    if (v > 0 && vmax < v) {
      vmax = v;
      fmax = f;
    }
  }
  if (vsum < 0.5) return (0);  // mostly background

  if (vmax < 0.5 && ct) {
    // no clear winner, resolve with tissue type
    nTT = ct->ctabTissueType->nentries - 1;
    // compute pvf for each tissue type
    for (tt = 0; tt < nTT; tt++) vtt[tt] = 0;
    for (f = 0; f < nsegs; f++) {
      segid = segidlist[f];
      tt = ct->entries[segid]->TissueType - 1;
      vtt[tt] += segpvfvox[f];
    }
    // find tissue type that has largest pvf
    vttmax = 0;
    ttmax = 0;
    for (tt = 0; tt < nTT; tt++) {
      if (vttmax < vtt[tt]) {
        vttmax = vtt[tt];
        ttmax = tt;
      }
    }
    // select the seg with the largest pvf from within the best TT
    // Note: this may have some dependency on order of the TTs
    vmax = 0;
    fmax = 0;
    for (f = 0; f < nsegs; f++) {
      segid = segidlist[f];
      tt = ct->entries[segid]->TissueType - 1;
      if (tt != ttmax) continue;
      v = segpvfvox[f];
      if (vmax < v) {
        vmax = v;
        fmax = f;
      }
    }
  }

  segid = segidlist[fmax];
  return (segid);
}

/*
  \fn MRI *MRIsegPVF2TissueTypePVF(MRI *segpvf, int *segidlist, int nsegs, COLOR_TABLE *ct, MRI *mask, MRI *pvf)
  \brief Converts a SegPVF (as computed by MRIseg2SegPVF()) to tissue
  type PVF, where the conversion from seg to TT is given in the ct. It
  may be more efficient to convert a highres seg to a tissue type seg
  using MRIseg2TissueType(), then run MRIseg2SegPVF() to generate the
  TT PVF directly.  */
MRI *MRIsegPVF2TissueTypePVF(MRI *segpvf, int *segidlist, int nsegs, COLOR_TABLE *ct, MRI *mask, MRI *pvf)
{
  int nTT, c, r, s, f, tt, segid;
  double vtt[100];
  nTT = ct->ctabTissueType->nentries - 1;

  if (pvf == NULL) {
    pvf = MRIallocSequence(segpvf->width, segpvf->height, segpvf->depth, MRI_FLOAT, nTT);
    if (pvf == NULL) {
      printf("ERROR: MRIsegPVF2TissueTypePVF(): could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(segpvf, pvf);
    MRIcopyPulseParameters(segpvf, pvf);
  }
  else
    MRIclear(pvf);

  for (c = 0; c < segpvf->width; c++) {
    for (r = 0; r < segpvf->height; r++) {
      for (s = 0; s < segpvf->depth; s++) {
        if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
        for (tt = 0; tt < nTT; tt++) vtt[tt] = 0;
        for (f = 0; f < nsegs; f++) {
          segid = segidlist[f];
          tt = ct->entries[segid]->TissueType - 1;
          vtt[tt] += MRIgetVoxVal(segpvf, c, r, s, f);
        }
        for (tt = 0; tt < nTT; tt++) MRIsetVoxVal(pvf, c, r, s, tt, vtt[tt]);
      }
    }
  }
  return (pvf);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------

  MRIaseg2vol() - maps a segmentation volume to another volume thru a
  registration file. voltemp is a geometry template for the new
  volume. The new volume is returned. Also returned (thru the arglist)
  is a pointer to a "hit" volume. This is the same size as the output
  volume. The value at each voxel is the number of aseg voxels mapped
  into that voxel from the seg id assigned to that voxel. If multiple
  seg ids map to a given voxel, then the seg id with the most number
  of hits is chosen.  It is possible that some output voxels will not
  be closest to any aseg voxel. So finds the output voxels that were
  not initially mapped and assigns them the seg id of the closest in
  the seg volume.  The number of hits for the winning seg id must
  exceed fthresh*nhitstot where nhitstot is the total number of aseg
  voxels that land in the given voxel. This cannot necessarily be
  computed from the voxel volumes because the voxel boundaries will
  not necessarily line up. If the number of hits does not exceed this
  value, then the segid is set to 0. 0 < fthresh < 1.

  The registration should work with:
     tkregister2 --targ  aseg --mov voltemp --reg tkR

  -------------------------------------------------------------*/
MRI *MRIaseg2vol(MRI *aseg, MATRIX *tkR, MRI *voltemp, double fthresh, MRI **pvolhit, COLOR_TABLE *ct)
{
  int Na, inda, sa, ra, ca, cv, rv, sv, indv, n, segid, nhits, nhitsmost;
  int nmisses, nfilled;
  MATRIX *Va2v, *Ka, *Kv, *invKv, *Vv2a, *Pa, *Pv;
  ASEGVOLINDEX *avind;
  MRI *volaseg, *volhit;

  if (ct != NULL) {
    if (ct->ctabTissueType == NULL) {
      printf("ERROR: MRIaeg2vol: if passing a ct it must have tissue type info\n");
      return (NULL);
    }
  }
  segid = 0;

  // Compute matrix to map from aseg CRS to vol CRS
  // Va2v = inv(Kv) * tkR * Ka
  Ka = MRIxfmCRS2XYZtkreg(aseg);
  Kv = MRIxfmCRS2XYZtkreg(voltemp);
  invKv = MatrixInverse(Kv, NULL);
  Va2v = MatrixMultiply(invKv, tkR, NULL);
  Va2v = MatrixMultiply(Va2v, Ka, Va2v);
  Vv2a = MatrixInverse(Va2v, NULL);
  MatrixFree(&Ka);
  MatrixFree(&Kv);
  MatrixFree(&invKv);

  Na = aseg->width * aseg->height * aseg->depth;
  avind = (ASEGVOLINDEX *)calloc(Na, sizeof(ASEGVOLINDEX));

  // Build an LUT that maps from aseg voxel to the closest
  // output volume voxel (reverse is done below).
  if (Gdiag_no > 0) printf("ASeg2Vol: Building LUT\n");
  fflush(stdout);
  Pa = MatrixConstVal(0, 4, 1, NULL);
  Pa->rptr[4][1] = 1;
  Pv = MatrixConstVal(0, 4, 1, NULL);
  inda = 0;
  for (sa = 0; sa < aseg->depth; sa++) {
    for (ra = 0; ra < aseg->height; ra++) {
      for (ca = 0; ca < aseg->width; ca++) {
        avind[inda].asegindex = inda;
        avind[inda].segid = MRIgetVoxVal(aseg, ca, ra, sa, 0);
        // map each aseg vox into the volume
        Pa->rptr[1][1] = ca;
        Pa->rptr[2][1] = ra;
        Pa->rptr[3][1] = sa;
        MatrixMultiply(Va2v, Pa, Pv);
        cv = (int)round(Pv->rptr[1][1]);
        rv = (int)round(Pv->rptr[2][1]);
        sv = (int)round(Pv->rptr[3][1]);

        if (cv < 0 || cv >= voltemp->width || rv < 0 || rv >= voltemp->height || sv < 0 || sv >= voltemp->depth) {
          // check for out-of-FOV
          avind[inda].volindex = -1;
          avind[inda].cv = -1;
          avind[inda].rv = -1;
          avind[inda].sv = -1;
        }
        else {
          // save the volume crs and index
          indv = cv + (rv * voltemp->width) + (sv * voltemp->width * voltemp->height);
          avind[inda].cv = cv;
          avind[inda].rv = rv;
          avind[inda].sv = sv;
          avind[inda].volindex = indv;
        }
        inda++;
      }  // col
    }    // row
  }      // slice

  if (Gdiag_no > 0) {
    printf("ASeg2Vol: Sorting \n");
    fflush(stdout);
  }
  qsort(avind, Na, sizeof(ASEGVOLINDEX), CompareAVIndices);

  // Alloc output volume
  volaseg = MRIalloc(voltemp->width, voltemp->height, voltemp->depth, MRI_INT);
  MRIcopyHeader(voltemp, volaseg);

  // Alloc a hit volume (all values should be 0)
  volhit = MRIalloc(voltemp->width, voltemp->height, voltemp->depth, MRI_INT);
  MRIcopyHeader(voltemp, volhit);
  *pvolhit = volhit;

  // Go through each volume voxel and determine which seg id has the
  // most representation.
  if (Gdiag_no > 0) {
    printf("ASeg2Vol: Mapping\n");
    fflush(stdout);
  }
  n = 0;
  while (n < Na) {
    indv = avind[n].volindex;
    if (indv < 0) {
      // aseg vox not in fov of outvol
      n++;
      continue;
    }

    // Count total number of aseg voxels falling into this volume vox
    nhits = 0;
    while ((n + nhits < Na) && (indv == avind[n + nhits].volindex)) nhits++;

    // Determine which segid had the most hits
    nhitsmost = MostHitsInVolVox(&avind[n], nhits, &segid, ct);
    MRIsetVoxVal(volhit, avind[n].cv, avind[n].rv, avind[n].sv, 0, nhitsmost);

    // Threshold (changed method on Jan 4, 2011)
    // New way based on number of segid hits relative to number of voxel hits
    if ((double)nhitsmost / nhits < fthresh) segid = 0;
    // Old way based on voxel volume
    // if (nhitsmost < nhitsthresh) segid=0;

    // Set the output
    MRIsetVoxVal(volaseg, avind[n].cv, avind[n].rv, avind[n].sv, 0, segid);

    if (Gdiag_no > 1) {
      printf("%6d %6d %5d    %6d    %3d %3d %3d     %3d %3d\n",
             n,
             avind[n].asegindex,
             segid,
             avind[n].volindex,
             avind[n].cv,
             avind[n].rv,
             avind[n].sv,
             nhits,
             nhitsmost);
      fflush(stdout);
    }

    n += nhits;
  }

  // Now do the reverse map. It is possible that some output voxels
  // will not be closest to any aseg voxel. So find the output voxels
  // that were not initially mapped and assign them the seg id of the
  // closest in the seg volume.
  if (Gdiag_no > 0) printf("ASeg2Vol: Reverse Map\n");
  nmisses = 0;
  nfilled = 0;
  for (sv = 0; sv < volaseg->depth; sv++) {
    for (rv = 0; rv < volaseg->height; rv++) {
      for (cv = 0; cv < volaseg->width; cv++) {
        nhits = MRIgetVoxVal(volhit, cv, rv, sv, 0);
        if (nhits != 0) continue;
        nmisses++;
        Pv->rptr[1][1] = cv;
        Pv->rptr[2][1] = rv;
        Pv->rptr[3][1] = sv;
        MatrixMultiply(Vv2a, Pv, Pa);
        ca = (int)round(Pa->rptr[1][1]);
        ra = (int)round(Pa->rptr[2][1]);
        sa = (int)round(Pa->rptr[3][1]);
        if (ca < 0 || ca >= aseg->width || ra < 0 || ra >= aseg->height || sa < 0 || sa >= aseg->depth) {
          // out-of-bounds
          MRIsetVoxVal(volaseg, cv, rv, sv, 0, 0);
        }
        else {
          // in-of-bounds
          segid = MRIgetVoxVal(aseg, ca, ra, sa, 0);
          MRIsetVoxVal(volaseg, cv, rv, sv, 0, segid);
          MRIsetVoxVal(volhit, cv, rv, sv, 0, 1);
          nfilled++;
        }
      }
    }
  }
  if (Gdiag_no > 0) printf("nmisses = %d (%d filled)\n", nmisses, nfilled);

  // MRIwrite(volaseg,"volaseg.mgh");
  // MRIwrite(volhit,"volhits.mgh");
  MatrixFree(&Va2v);
  MatrixFree(&Vv2a);
  free(avind);

  if (Gdiag_no > 0) printf("ASeg2Vol: done\n");
  return (volaseg);
}
/*
  \fn MRI *MRIaseg2volMU(MRI *aseg, LTA *aseg2vol, double fthresh, MRI **pvolhit, int USF, COLOR_TABLE *ct)
  \brief Frontend for MRIaseg2vol(). This version takes a LTA and
  USF. It will reduce the FoV of the aseg to the smallest bounding
  box. This will then be upsampled by USF. The result (and a proper
  matrix) is handed off to MRIaseg2vol(). The aseg2vol LTA can point
  in either direction.
*/
MRI *MRIaseg2volMU(MRI *aseg, LTA *aseg2vol, double fthresh, MRI **pvolhit, int USF, COLOR_TABLE *ct)
{
  MRI *asegmu, *OutVol, *TempVol;
  LTA *aseg2asegmu, *asegmu2aseg, *asegmu2vol, *ltaArray[2];
  int nPad = 2;
  VOL_GEOM *vg;

  if (USF != -1) {
    asegmu = MRImaskAndUpsample(aseg, NULL, USF, nPad, 0, &aseg2asegmu);
    asegmu2aseg = LTAinvert(aseg2asegmu, NULL);
    ltaArray[0] = asegmu2aseg;
    ltaArray[1] = aseg2vol;
    asegmu2vol = LTAconcat(ltaArray, 2, 1);  // figures out inversions
    LTAfree(&aseg2asegmu);
  }
  else {
    asegmu = aseg;
    if (LTAmriIsSource(aseg2vol, aseg))
      asegmu2vol = LTAcopy(aseg2vol, NULL);
    else
      asegmu2vol = LTAinvert(aseg2vol, NULL);
  }

  // MRIaseg2vol() requires a register.dat matrix
  if (LTAmriIsSource(asegmu2vol, asegmu)) {
    /* The definitions of mov=src and ref=dst are consistent with
       tkregister2, LTAchangeType() and ltaReadRegisterDat(). This is an
       unfortunate definition because the registration matrix actually
       does from ref to mov. But this was an error introduced a long
       time ago and the rest of the code base has built up around it. */
    LTAinvert(asegmu2vol, asegmu2vol);
  }
  LTAchangeType(asegmu2vol, REGISTER_DAT);

  if (LTAmriIsSource(aseg2vol, aseg))
    vg = &aseg2vol->xforms[0].dst;
  else
    vg = &aseg2vol->xforms[0].src;
  TempVol = MRIallocHeader(vg->width, vg->height, vg->depth, MRI_INT, 1);
  useVolGeomToMRI(vg, TempVol);
  MRIcopyPulseParameters(aseg, TempVol);

  OutVol = MRIaseg2vol(asegmu, asegmu2vol->xforms[0].m_L, TempVol, fthresh, pvolhit, ct);

  if (asegmu != aseg) MRIfree(&asegmu);
  LTAfree(&asegmu2vol);
  MRIfree(&TempVol);
  return (OutVol);
}
/*
  \fn MRI *MRIchangeSegRes(MRI *seg, double xsize, double ysize, double zsize, COLOR_TABLE *ct, LTA **seg2new)
  \brief Changes the resolution of the segmentation to
  {x,y,z}size. I'm not sure how the results will differ from
  MRIaseg2volMU() or MRIaseg2vol(). This uses MRIseg2SegPVF() and MRIaseg2vol()
  uses a voting scheme. I would guess that they are close.
 */
MRI *MRIchangeSegRes(MRI *seg, double xsize, double ysize, double zsize, COLOR_TABLE *ct, LTA **seg2new)
{
  int *segidlist, nsegs, ReGridFactor = -2, ReInitCache = 1;
  MRI *newseg, *mritmp;

  segidlist = MRIsegIdListNot0(seg, &nsegs, 0);

  mritmp = MRIresize(seg, xsize, ysize, zsize, 0);
  *seg2new = TransformRegDat2LTA(mritmp, seg, NULL);
  newseg = MRIseg2SegPVF(seg, *seg2new, ReGridFactor, segidlist, nsegs, NULL, ReInitCache, ct, NULL);
  if (newseg == NULL) return (NULL);

  MRIfree(&mritmp);
  free(segidlist);
  return (newseg);
}

/*-----------------------------------------------------------------------
  MostHitsInVolVox() - determines the segid with the most hits in the
  given volume voxel. In this case, avindsorted points to the first
  entry of the given voxel, and there are N entries for that voxel.
  -----------------------------------------------------------------------*/
static int MostHitsInVolVox(ASEGVOLINDEX *avindsorted, int N, int *segidmost, COLOR_TABLE *ct)
{
  int n, nhits = 0, nmost = 0, segid, nsegs;
  static int segidlist[1000], ttypelist[1000], nhitslist[1000];
  static int nPerTType[100], nPerTTypeMax, TTypeMax;

  bzero(nPerTType, sizeof(int) * 100);

  n = 0;
  nsegs = 0;
  while (n < N) {
    segid = avindsorted[n].segid;
    // count number of hits for this segid
    nhits = 0;
    while ((n + nhits < N) && (segid == avindsorted[n + nhits].segid)) nhits++;
    if (n == 0) {
      // first segid
      *segidmost = segid;
      nmost = nhits;
    }
    if (nmost < nhits) {
      // update if needed
      *segidmost = segid;
      nmost = nhits;
    }
    n += nhits;
    if (ct) {
      segidlist[nsegs] = segid;
      nhitslist[nsegs] = nhits;
      ttypelist[nsegs] = ct->entries[segid]->TissueType;
      if (ct->entries[segid]->TissueType >= 0) nPerTType[ct->entries[segid]->TissueType] += nhits;
      if (ct->entries[segid]->TissueType == -2)
        printf(
            "WARNING: color table entry %d %s has not been assigned a tissue type\n", segid, ct->entries[segid]->name);
    }
    nsegs++;  // keep count of number of segs
  }

  if (ct && (float)nhits / N < 0.5) {
    // no clear winner, choose one from tissue type with most representation
    nPerTTypeMax = 0;
    for (n = 0; n < ct->ctabTissueType->nentries; n++) {
      if (nPerTTypeMax < nPerTType[n]) {
        nPerTTypeMax = nPerTType[n];
        TTypeMax = n;
      }
    }
    *segidmost = segidlist[0];
    nmost = 0;
    for (n = 0; n < nsegs; n++) {
      if (ttypelist[n] != TTypeMax) continue;
      if (nmost > nhitslist[n]) continue;
      nmost = nhitslist[n];
      *segidmost = segidlist[n];
    }
  }

  return (nmost);
}

/*-----------------------------------------------------------------------
  CompareAVIndices() - this compares ASEGVOLINDEX structures in a way
  compatible with qsort. The resulting sorted array will be first sorted
  by volume index. Entries with the same volume index will be subsorted
  by segmentation id. In the end, all entries with the same volume index
  will be contiguous. For a set of entries with the same volume index,
  all entries with the same seg id will be contiguous.
  -----------------------------------------------------------------------*/
static int CompareAVIndices(const void *i1, const void *i2)
{
  ASEGVOLINDEX *avind1 = NULL, *avind2 = NULL;

  avind1 = (ASEGVOLINDEX *)i1;
  avind2 = (ASEGVOLINDEX *)i2;

  // Sort by volume index
  if (avind1->volindex > avind2->volindex) return (+1);
  if (avind1->volindex < avind2->volindex) return (-1);

  // Only gets here if vol indices are the same, so
  // sort by seg id
  if (avind1->segid > avind2->segid) return (+1);
  if (avind1->segid < avind2->segid) return (-1);

  // Same volume index and seg id.
  return (0);
}


/*!
  \fn MRI *MRIapplySpmWarp(MRI *vol, LTA *srclta, MRI *warp, int LRRev, int interp, MRI *out)
  \brief Applies a warp field computed from SPM, eg, with DNG's
  run-vbm (which uses DARTEL).  The input vol is anything that shares
  the scanner space with the input to vbm.  The warp field is
  y_rinput.nii. LRRev=0,1 indicates that the pixels in the anatomical
  were left-right reversed before the warp was computed. I thought
  thiwas good for asym studies, but it turns out it is not
  needed. This will only work when the column of the warp input is in
  the left-right direction. Interp: SAMPLE_NEAREST=0 or
  SAMPLE_TRILINEAR=1.  Handles multiple frames.  The output will be in
  the warp space.
 */
MRI *MRIapplySpmWarp(MRI *vol, LTA *srclta, MRI *warp, int LRRev, int interp, MRI *out)
{
  int c,r,s,k,f,ncols;
  double cc,rr,ss,val;
  MATRIX *Q=NULL;

  if(interp != SAMPLE_NEAREST && interp != SAMPLE_TRILINEAR){
    printf("ERROR: MRIapplySpmWarp():  sample type = %d, must be %d or %d\n",
	   interp,SAMPLE_NEAREST,SAMPLE_TRILINEAR);
    return(NULL);
  }

  if(out==NULL){
    out = MRIallocSequence(warp->width,warp->height,warp->depth,MRI_FLOAT,vol->nframes);
    if(out==NULL) return(NULL);
    MRIcopyHeader(warp,out);
    MRIcopyPulseParameters(vol, out);
  }
  
  // Compute vbminput-to-vol vox2vox. Note: SPM assumes the CRS of the
  // vbm input is is 1-based. Everything else, including the LTA, is
  // 0-based.
  MATRIX *vox2ras1;  // for vbm input
  MATRIX *Vsrc;
  if(srclta){
    if(srclta->type != LINEAR_VOX_TO_VOX){
      printf("MRIapplySpmWarp(): changing Source LTA type to vox2vox\n");
      LTAchangeType(srclta, LINEAR_VOX_TO_VOX) ;
    }
    // If a srclta is specified, then the vol is not the same geom as
    // the vbm input. srclta must map between the vol and the vbm
    // input space. Set the vbm input vox2ras1 from the approp vol geom
    if(LTAmriIsSource(srclta, vol)){
      vox2ras1 = VGgetVoxelToRasXform(&(srclta->xforms[0].dst),NULL,1); // spm crs base=1
      printf("MRIapplySpmWarp(): inverting Source LTA\n");
      Vsrc = MatrixInverse(srclta->xforms[0].m_L,NULL);
      ncols = srclta->xforms[0].dst.width;
    }
    else {
      Vsrc = MatrixCopy(srclta->xforms[0].m_L,NULL);
      vox2ras1 = VGgetVoxelToRasXform(&(srclta->xforms[0].src),NULL,1); // spm crs base=1
      ncols = srclta->xforms[0].src.width;
    }
  }
  else {
    vox2ras1 = MRIxfmCRS2XYZ(vol, 1); // spm crs base=1
    Vsrc = MatrixIdentity(4,NULL);
    ncols = vol->width;
  }
  MATRIX *ras2vox1 = MatrixInverse(vox2ras1,NULL);

  if(LRRev){
    // This is the matrix that realizes the left-right pixel reversal
    // used in mri_convert --left-right-reverse-pix. Only works properly 
    // if the column is in the LR direction 
    Q = MatrixIdentity(4,NULL);
    Q->rptr[1][1] = -1;
    Q->rptr[1][4] = ncols-1;
  }

  MATRIX *vbmiRAS = MatrixAlloc(4,1,MATRIX_REAL);
  vbmiRAS->rptr[4][1] = 1;
  MATRIX *CRS  = MatrixAlloc(4,1,MATRIX_REAL);
  CRS->rptr[4][1] = 1;
  // use 1 to N-1 instead of 1 to N because edge voxels are invalid in the warp
  for(c=1; c < warp->width-1; c++){
    for(r=1; r < warp->height-1; r++){
      for(s=1; s < warp->depth-1; s++){
	// Get the RAS in the vbm input space
	for(k=0; k<3; k++) vbmiRAS->rptr[k+1][1] = MRIgetVoxVal(warp,c,r,s,k);
	// Get the 1-based CRS in the vbm-input space (usually a conformed space)
	CRS = MatrixMultiplyD(ras2vox1,vbmiRAS,CRS);
	// Subtract 1 to make 0-based (could do this in ras2vox1)
	CRS->rptr[1][1] -= 1;
	CRS->rptr[2][1] -= 1;
	CRS->rptr[3][1] -= 1;
	// If left-right rev is needed, do it here
	if(LRRev) CRS = MatrixMultiplyD(Q,CRS,CRS);
	if(srclta != NULL){
	  // Now compute the CRS in the vol space
	  // Could be combined with ras2vox1 if the conversion from 1-to-0-based
	  // is figured out.
	  CRS = MatrixMultiplyD(Vsrc,CRS,CRS);
	}
	cc = CRS->rptr[1][1];
	rr = CRS->rptr[2][1];
	ss = CRS->rptr[3][1];
	if(cc < 0 || cc >= vol->width)  continue;
	if(rr < 0 || rr >= vol->height) continue;
	if(ss < 0 || ss >= vol->depth)  continue;
	for(f=0; f < vol->nframes; f++){
	  MRIsampleVolumeFrameType(vol, cc, rr, ss, f, interp, &val);
	  MRIsetVoxVal(out,c,r,s,f,val);
	}
      }
    }
  }
  MatrixFree(&vox2ras1);
  MatrixFree(&ras2vox1);
  MatrixFree(&vbmiRAS);
  MatrixFree(&CRS);
  MatrixFree(&Vsrc);

  return(out);
}
