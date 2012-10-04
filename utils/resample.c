/**
 * @file  resample.c
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
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2012/10/04 17:51:00 $
 *    $Revision: 1.42 $
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

#define RESAMPLE_SOURCE_CODE_FILE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "diag.h"
#include "matrix.h"
#include "mri.h"
#include "mri2.h"
#include "mrisurf.h"
#include "mrishash.h"
#include "label.h"
#include "resample.h"
#include "bfileio.h"
#include "corio.h"
#include "proto.h" // nint

/*-------------------------------------------------------------------*/
double round(double); // why is this never defined?!?
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------
  ASEGVOLINDEX - this structure is used with MRIaseg2vol() to help
  perform the mapping. Mainly used for sorting with qsort.
  -------------------------------------------------------------------*/
typedef struct
{
  int asegindex; // index into the seg volume (instead of col,row,slice)
  int segid;     // segmentation code
  int volindex;  // index into the output volume
  int cv,rv,sv;  // corresponding col, row, slice in the output volume
}
ASEGVOLINDEX;

static int CompareAVIndices(const void *i1, const void *i2);
static int MostHitsInVolVox(ASEGVOLINDEX *avindsorted, int N, int *segidmost);

/*---------------------------------------------------------
  interpolation_code(): gets a code that controls the
  interpolation from a string. Returns -1 if the string is
  unrecoginzed.
  ---------------------------------------------------------*/
int interpolation_code(char *interpolation_string)
{
  return(MRIinterpCode(interpolation_string));
}
/*---------------------------------------------------------
  float2int_code(): gets a code that controls the floating
  point to integer conversion method from a string. Returns
  -1 if the string is unrecoginzed.
  ---------------------------------------------------------*/
int float2int_code(char *float2int_string)
{
  if (!strcasecmp(float2int_string,"round") ||
      !strcasecmp(float2int_string,"rint"))
    return(FLT2INT_ROUND);

  if (!strcasecmp(float2int_string,"floor")) return(FLT2INT_FLOOR);

  if (!strncasecmp(float2int_string,"tkregister",5)) return(FLT2INT_TKREG);

  return(-1);
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
int XYZAnat2CRSFunc_TkReg(int *col, int *row, int *slc,
                          int npixels, float pixsize,
                          int nslcs,   float slcthick,
                          float xanat, float yanat, float zanat,
                          MATRIX *Reg)
{
  MATRIX *Qf;
  MATRIX *QfR;
  MATRIX *xyz, *crs;

  /* compute the tkregister quantization matrix */
  Qf = FOVQuantMtx_TkReg(npixels, pixsize, nslcs, slcthick);

  /* QfR = Qf* R*/
  QfR = MatrixMultiply(Qf,Reg,NULL);

  /* load xyz into a vector */
  xyz = MatrixAlloc(4,1,MATRIX_REAL);
  xyz->rptr[0+1][0+1] = xanat;
  xyz->rptr[1+1][0+1] = yanat;
  xyz->rptr[2+1][0+1] = zanat;
  xyz->rptr[3+1][0+1] = 1.0;

  /* crs = Qf*R*xyz */
  crs = MatrixMultiply(QfR,xyz,NULL);

  /* extract col, row, slice from the crs vector */
  /* convert floats to ints as only tkregister can */
  float2int_TkReg(col,row,slc,
                  crs->rptr[0+1][0+1],
                  crs->rptr[1+1][0+1],
                  crs->rptr[2+1][0+1]);

  MatrixFree(&Qf);
  MatrixFree(&QfR);
  MatrixFree(&xyz);
  MatrixFree(&crs);

  return(0);
}
/*----------------------------------------------------------
  float2int_TkReg() - converts floats to int as only
  tkregsiter can. When volumes are coronally sliced, then
  col is LR (sagital), row is SI (axial), and slc is AP (cor)
  ----------------------------------------------------------*/
int float2int_TkReg(int *col, int *row, int *slc,
                    float fltcol, float fltrow, float fltslc)
{
  *col = (int)(fltcol);
  *row = (int)(ceil(fltrow));
  *slc = (int)(fltslc);
  return(0);
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
MATRIX * FOVQuantMtx_TkReg(int npixels, float pixsize,
                           int nslcs,   float slcthick)
{
  MATRIX *Q;

  Q = MatrixAlloc(4, 4, MATRIX_REAL );
  MatrixClear(Q);

  Q->rptr[0+1][0+1] = -1/pixsize;
  Q->rptr[0+1][3+1] =  npixels/2;

  Q->rptr[1+1][2+1] =  1/slcthick;
  Q->rptr[1+1][3+1] = -nslcs/2;

  Q->rptr[2+1][1+1] = -1/pixsize;
  Q->rptr[2+1][3+1] =  npixels/2;

  Q->rptr[3+1][3+1] = 1.0;

  return(Q);
}
/*-----------------------------------------------------------
  FOVDeQuantMatrix() -- computes the volume dequantization matrix which
  converts a col,row,slc into an x,y,z: xyz = deQ*crs
-----------------------------------------------------------*/
MATRIX * FOVDeQuantMatrix(int ncols, int nrows, int nslcs,
                          float colres, float rowres, float slcres  )
{
  MATRIX *deQ;

  deQ = MatrixAlloc(4, 4, MATRIX_REAL );
  MatrixClear(deQ);

  deQ->rptr[0+1][0+1] = -colres;
  deQ->rptr[0+1][3+1] =  colres*(ncols)/2;

  deQ->rptr[1+1][2+1] =  slcres;
  deQ->rptr[1+1][3+1] = -slcres*(nslcs)/2;

  deQ->rptr[2+1][1+1] = -rowres;
  deQ->rptr[2+1][3+1] =  rowres*(nrows)/2;

  deQ->rptr[3+1][3+1] = 1.0;

  return(deQ);
}
/*-----------------------------------------------------------
  FOVQuantMatrix() -- computes the volume quantization matrix which
  converts a x,y,z into col,row,slc : crs = Q*xyz
  -----------------------------------------------------------*/
MATRIX * FOVQuantMatrix(int ncols, int nrows, int nslcs,
                        float colres, float rowres, float slcres  )
{
  MATRIX *deQ, *Q;

  deQ = FOVDeQuantMatrix(ncols,nrows,nslcs,colres,rowres,slcres);
  Q   = MatrixInverse(deQ,NULL);
  MatrixFree(&deQ);

  return(Q);
}
/*------------------------------------------------------------
  ComputeQFWD() - computes the matrix product of Q, F, W, D.
  If any matrix is NULL, then it is treated as the idenity.
  ------------------------------------------------------------*/
MATRIX *ComputeQFWD(MATRIX *Q, MATRIX *F, MATRIX *W, MATRIX *D, MATRIX *QFWD)
{
  MATRIX *QFWDtmp;

  if (QFWD==NULL) QFWDtmp = MatrixAlloc(4,4,MATRIX_REAL);
  else           QFWDtmp = QFWD;

  MatrixIdentity(4,QFWDtmp);

  if (Q != NULL)  MatrixMultiply(QFWDtmp, Q, QFWDtmp);
  if (F != NULL)  MatrixMultiply(QFWDtmp, F, QFWDtmp);
  if (W != NULL)  MatrixMultiply(QFWDtmp, W, QFWDtmp);
  if (D != NULL)  MatrixMultiply(QFWDtmp, D, QFWDtmp);

  return(QFWDtmp);
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
                    MATRIX *Qsrc, MATRIX *Fsrc, MATRIX *Wsrc, MATRIX *Dsrc,
                    MATRIX *Qtrg, MATRIX *Ftrg, MATRIX *Wtrg, MATRIX *Dtrg,
                    int   nrows_trg, int   ncols_trg, int   nslcs_trg,
                    MATRIX *Msrc2trg, int InterpMethod, int float2int)
{
  MATRIX *QFWDsrc, *QFWDtrg, *invQFWDtrg;
  MATRIX *Tcrs2Scrs, *invMsrc2trg;
  MATRIX *Scrs, *Tcrs, *Txyz;
  MRI *TrgVol;
  int   irow_trg, icol_trg, islc_trg; /* integer row, col, slc in target */
  int   irow_src, icol_src, islc_src; /* integer row, col, slc in source */
  float srcval;
  int frm;

  if (InterpMethod != INTERP_NEAREST)
  {
    fprintf(stderr,"vol2vol_linear(): only support for nearest interpolation\n");
    return(NULL);
  }

  /* compute the transforms */
  QFWDsrc = ComputeQFWD(Qsrc,Fsrc,Wsrc,Dsrc,NULL);
  QFWDtrg = ComputeQFWD(Qtrg,Ftrg,Wtrg,Dtrg,NULL);
  invQFWDtrg = MatrixInverse(QFWDtrg,NULL);
  if (Msrc2trg != NULL) invMsrc2trg = MatrixInverse(Msrc2trg,NULL);
  else                 invMsrc2trg = MatrixIdentity(4,NULL);
  Tcrs2Scrs = MatrixMultiply(QFWDsrc,invMsrc2trg,NULL);
  MatrixMultiply(Tcrs2Scrs , invQFWDtrg, Tcrs2Scrs);

  /* Tcrs2Scrs - maps Target ColRowSlice to that of the Source */

  /* preallocate the row-col-slc vectors */
  Tcrs = MatrixAlloc(4,1,MATRIX_REAL);
  Tcrs->rptr[3+1][0+1] = 1.0;
  Scrs = MatrixAlloc(4,1,MATRIX_REAL);
  Txyz = MatrixAlloc(4,1,MATRIX_REAL);
  Txyz->rptr[3+1][0+1] = 1.0;

  /* allocate a volume to hold the output */
  TrgVol = MRIallocSequence(ncols_trg, nrows_trg, nslcs_trg,
                            MRI_FLOAT,SrcVol->nframes);
  if (TrgVol == NULL) return(NULL);

  /* Go through each target voxel and compute the location of
     the closest source voxel */
  for (islc_trg = 0; islc_trg < nslcs_trg; islc_trg++)
  {
    for (irow_trg = 0; irow_trg < nrows_trg; irow_trg++)
    {
      for (icol_trg = 0; icol_trg < ncols_trg; icol_trg++)
      {

        /* Load the Target col-row-slc vector */
        Tcrs->rptr[0+1][0+1] = icol_trg;
        Tcrs->rptr[1+1][0+1] = irow_trg;
        Tcrs->rptr[2+1][0+1] = islc_trg;

        /* Compute the corresponding Source col-row-slc vector */
        MatrixMultiply(Tcrs2Scrs,Tcrs,Scrs);

        /* nearest neighbor */
        switch (float2int)
        {
        case FLT2INT_ROUND:
          icol_src = (int)rint(Scrs->rptr[0+1][0+1]);
          irow_src = (int)rint(Scrs->rptr[1+1][0+1]);
          islc_src = (int)rint(Scrs->rptr[2+1][0+1]);
          break;
        case FLT2INT_FLOOR:
          icol_src = (int)floor(Scrs->rptr[0+1][0+1]);
          irow_src = (int)floor(Scrs->rptr[1+1][0+1]);
          islc_src = (int)floor(Scrs->rptr[2+1][0+1]);
          break;
        case FLT2INT_TKREG:
          icol_src = (int)floor(Scrs->rptr[0+1][0+1]);
          irow_src = (int) ceil(Scrs->rptr[1+1][0+1]);
          islc_src = (int)floor(Scrs->rptr[2+1][0+1]);
          break;
        default:
          fprintf(stderr,"vol2vol_linear(): unrecoginized float2int code %d\n",
                  float2int);
          MRIfree(&TrgVol);
          return(NULL);
          break;
        }

        /* make sure the Source Voxel is in the volume */
        if (irow_src < 0 || irow_src >= SrcVol->height ||
            icol_src < 0 || icol_src >= SrcVol->width  ||
            islc_src < 0 || islc_src >= SrcVol->depth ) continue;

        /* map each frame */
        for (frm = 0; frm < SrcVol->nframes; frm++)
        {
          srcval = MRIFseq_vox(SrcVol,icol_src,irow_src,islc_src,frm);
          MRIFseq_vox(TrgVol,icol_trg,irow_trg,islc_trg,frm) = srcval;
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

  return(TrgVol);
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
                       MATRIX *Qsrc, MATRIX *Fsrc, MATRIX *Wsrc, MATRIX *Dsrc,
                       MRI *SrcMskVol, MATRIX *Msrc2lbl, LABEL *Label,
                       float rszthresh, int float2int, int *nlabelhits, int *nfinalhits)
{
  MATRIX *QFWDsrc;
  MATRIX *Lxyz2Scrs;
  MATRIX *Scrs, *Lxyz, *Mlbl2src;
  MRI *FinalMskVol;
  int   irow_src, icol_src, islc_src; /* integer row, col, slc in source */
  float mskval, voxsize;
  int vlbl;
  int c,r,s,nfinalmask;

  if(SrcMskVol){
    if(MRIdimMismatch(SrcVol, SrcMskVol, 0)){
      printf("ERROR: label2mask_linear: dimension mismatch\n");
      return(NULL);
    }
  }

  /* compute the transforms */
  QFWDsrc = ComputeQFWD(Qsrc,Fsrc,Wsrc,Dsrc,NULL);
  if (Msrc2lbl != NULL) Mlbl2src = MatrixInverse(Msrc2lbl,NULL);
  else                 Mlbl2src = NULL;
  if (Mlbl2src != NULL) Lxyz2Scrs = MatrixMultiply(QFWDsrc,Mlbl2src,NULL);
  else                 Lxyz2Scrs = MatrixCopy(QFWDsrc,NULL);

  printf("\n");
  printf("Lxyz2Scrs:\n");
  MatrixPrint(stdout,Lxyz2Scrs);
  printf("\n");

  /* preallocate the row-col-slc vectors */
  Lxyz = MatrixAlloc(4,1,MATRIX_REAL);
  Lxyz->rptr[3+1][0+1] = 1.0;
  Scrs = MatrixAlloc(4,1,MATRIX_REAL);

  /* allocate an output volume -- same size as source*/
  FinalMskVol = MRIallocSequence(SrcVol->width,SrcVol->height,SrcVol->depth,
                                 MRI_FLOAT,1);
  if (FinalMskVol == NULL) return(NULL);
  MRIcopyHeader(SrcVol,FinalMskVol);

  *nlabelhits = 0;
  *nfinalhits = 0;

  /* Go through each point in the label */
  for (vlbl = 0; vlbl < Label->n_points; vlbl++)
  {

    /* load the label xyz into a vector */
    Lxyz->rptr[0+1][0+1] = Label->lv[vlbl].x;
    Lxyz->rptr[1+1][0+1] = Label->lv[vlbl].y;
    Lxyz->rptr[2+1][0+1] = Label->lv[vlbl].z;

    /* compute the corresponding col, row, and slice in the source vol */
    MatrixMultiply(Lxyz2Scrs,Lxyz,Scrs);

    /* Convert the analog col, row, and slice to integer */
    switch (float2int)
    {
    case FLT2INT_ROUND:
      icol_src = (int)rint(Scrs->rptr[0+1][0+1]);
      irow_src = (int)rint(Scrs->rptr[1+1][0+1]);
      islc_src = (int)rint(Scrs->rptr[2+1][0+1]);
      break;
    case FLT2INT_FLOOR:
      icol_src = (int)floor(Scrs->rptr[0+1][0+1]);
      irow_src = (int)floor(Scrs->rptr[1+1][0+1]);
      islc_src = (int)floor(Scrs->rptr[2+1][0+1]);
      break;
    case FLT2INT_TKREG:
      icol_src = (int)floor(Scrs->rptr[0+1][0+1]);
      irow_src = (int) ceil(Scrs->rptr[1+1][0+1]);
      islc_src = (int)floor(Scrs->rptr[2+1][0+1]);
      break;
    default:
      fprintf(stderr,"label2mask_linear(): unrecoginized float2int code %d\n",
              float2int);
      MRIfree(&FinalMskVol);
      return(NULL);
      break;
    }

    /* check that the point is within the source volume */
    if (irow_src < 0 || irow_src >= SrcVol->height ||
        icol_src < 0 || icol_src >= SrcVol->width  ||
        islc_src < 0 || islc_src >= SrcVol->depth ) continue;
    (*nlabelhits)++;

    /* check that the point is within the input mask */
    if (SrcMskVol != NULL)
    {
      mskval = MRIFseq_vox(SrcMskVol,icol_src,irow_src,islc_src,0);
      if (mskval < 0.5) continue;
    }

    /* keep count; binarization is done below */
    MRIFseq_vox(FinalMskVol,icol_src,irow_src,islc_src,0) ++;

    /* increment the number of hits */
    (*nfinalhits)++;
  }

  printf("INFO: label2mask_linear: there were %d label hits\n", *nfinalhits);

  /* binarize based on the number of hits in each voxel */
  voxsize = SrcVol->xsize * SrcVol->ysize * SrcVol->zsize;
  printf("voxsize = %g, rszthresh = %g, thresh = %g\n",
         voxsize,rszthresh,voxsize*rszthresh);
  nfinalmask = 0;
  for (c=0; c < FinalMskVol->width; c++)
  {
    for (r=0; r < FinalMskVol->height; r++)
    {
      for (s=0; s < FinalMskVol->depth;  s++)
      {
        mskval = MRIFseq_vox(FinalMskVol,c,r,s,0);
        if (mskval < rszthresh*voxsize)
          MRIFseq_vox(FinalMskVol,c,r,s,0) = 0;
        else
        {
          MRIFseq_vox(FinalMskVol,c,r,s,0) = 1;
          nfinalmask ++;
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

  return(FinalMskVol);
}

/*-----------------------------------------------------------------------
  vol2maskavg() -- averages all the voxels within a mask.  SrcVol and
  SrcMskVol must be the same size.  Searchs SrcMskVol for voxels whose
  value is > 0.5. For all such voxels, averages the corresponding values
  found in SrcVol.  Returns an MRI "volume" of dimension 1X1X1Xnframes.
  -----------------------------------------------------------------------*/
MRI * vol2maskavg(MRI *SrcVol, MRI *SrcMskVol, int *nhits)
{
  int r,c,s,f;
  MRI *MskAvg;
  float mskval, val;

  /* make sure that SrcVol and SrcMskVol are the same dimension */
  if (SrcVol->width  != SrcMskVol->width ||
      SrcVol->height != SrcMskVol->height ||
      SrcVol->depth  != SrcMskVol->depth)
  {
    fprintf(stderr,"ERROR: vol2maskavg: SrcVol and SrcMskVol do "
            "not have the same dimension\n");
    return(NULL);
  }

  /* allocate an output "volume" */
  MskAvg = MRIallocSequence(1,1,1,MRI_FLOAT,SrcVol->nframes);
  if (MskAvg == NULL) return(NULL);

  /* go through each voxel */
  *nhits = 0;
  for (c=0; c < SrcVol->width; c++)
  {
    for (r=0; r < SrcVol->height; r++)
    {
      for (s=0; s < SrcVol->depth; s++)
      {
        /* get mask value at r,c,s */
        mskval = MRIgetVoxVal(SrcMskVol,c,r,s,0);
        if (mskval > 0.5)
        {
          /* accumulate sum over suprathreshold points */
          (*nhits)++;
          for (f=0; f < SrcVol->nframes; f++)
          {
            val = MRIgetVoxVal(SrcVol,c,r,s,f);
            MRIFseq_vox(MskAvg,0,0,0,f) += val;
          }
        }
      }
    }
  }

  printf("INFO: vol2maskavg: nhits = %d\n",*nhits);
  if (*nhits != 0)
  {
    /* divide by the number of hits to get the average*/
    for (f=0; f < SrcVol->nframes; f++)
    {
      val = MRIFseq_vox(MskAvg,0,0,0,f);
      MRIFseq_vox(MskAvg,0,0,0,f) = val/(*nhits);
      //printf("%2d %g %g\n",f,val,MRIFseq_vox(MskAvg,0,0,0,f));
    }
  }
  else
  {
    printf("WARNING: there were no voxels in the input mask > 0.5\n");
  }

  /* return the Masked Average */
  return(MskAvg);
}
/*----------------------------------------------------------------
  ProjNormFracThick() - projects along the surface normal a given
  fraction of the thickness at that point.
  ----------------------------------------------------------------*/
int ProjNormFracThick( float *x, float *y, float *z,
                       const MRI_SURFACE *surf, int vtx, float frac )
{
  float r;
  r = frac * surf->vertices[vtx].curv;
  *x = surf->vertices[vtx].x + r*surf->vertices[vtx].nx;
  *y = surf->vertices[vtx].y + r*surf->vertices[vtx].ny;
  *z = surf->vertices[vtx].z + r*surf->vertices[vtx].nz;
  return(0);
}
/*----------------------------------------------------------------
  ProjNormFracThickNbr() - projects along the direction of the average of
  the surface normals of a vertex and its neighbor.  The distance
  along this direction is a fraction of the thickness at that point.
  ----------------------------------------------------------------*/
int ProjNormFracThickNbr(float *x, float *y, float *z, MRI_SURFACE *surf, 
			 int vtxno, float frac, int nthNbr)
{
  float r, nx, ny, nz;
  int nbrvtxno;

  nbrvtxno = surf->vertices[vtxno].v[nthNbr];
  nx = (surf->vertices[vtxno].nx + surf->vertices[nbrvtxno].nx)/2.0;
  ny = (surf->vertices[vtxno].ny + surf->vertices[nbrvtxno].ny)/2.0;
  nz = (surf->vertices[vtxno].nz + surf->vertices[nbrvtxno].nz)/2.0;
  r = frac * surf->vertices[vtxno].curv;
  *x = surf->vertices[vtxno].x + r*nx;
  *y = surf->vertices[vtxno].y + r*ny;
  *z = surf->vertices[vtxno].z + r*nz;
  return(0);
}
/*----------------------------------------------------------------
  ProjNormDist() - projects along the surface normal a given
  distance.
  ----------------------------------------------------------------*/
int ProjNormDist( float *x, float *y, float *z,
                  const MRI_SURFACE *surf, int vtx, float dist)
{
  *x = surf->vertices[vtx].x + dist*surf->vertices[vtx].nx;
  *y = surf->vertices[vtx].y + dist*surf->vertices[vtx].ny;
  *z = surf->vertices[vtx].z + dist*surf->vertices[vtx].nz;
  return(0);
}
/*----------------------------------------------------------------
  ProjNormDistNbr() - projects along the direction of the average of
  the surface normals of a vertex and its neighbor.  The distance
  along this direction is dist.
  ----------------------------------------------------------------*/
int ProjNormDistNbr(float *x, float *y, float *z, MRI_SURFACE *surf, 
		    int vtxno, float dist, int nthNbr)
{
  float nx, ny, nz;
  int nbrvtxno;

  nbrvtxno = surf->vertices[vtxno].v[nthNbr];
  nx = (surf->vertices[vtxno].nx + surf->vertices[nbrvtxno].nx)/2.0;
  ny = (surf->vertices[vtxno].ny + surf->vertices[nbrvtxno].ny)/2.0;
  nz = (surf->vertices[vtxno].nz + surf->vertices[nbrvtxno].nz)/2.0;
  *x = surf->vertices[vtxno].x + dist*nx;
  *y = surf->vertices[vtxno].y + dist*ny;
  *z = surf->vertices[vtxno].z + dist*nz;
  return(0);
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
                     MATRIX *Qsrc, MATRIX *Fsrc, MATRIX *Wsrc, MATRIX *Dsrc,
                     MRI_SURFACE *TrgSurf, float ProjFrac,
                     int InterpMethod, int float2int, MRI *SrcHitVol,
                     int ProjDistFlag, int nskip)
{
  MATRIX *QFWDsrc;
  MATRIX *Scrs, *Txyz;
  MRI *TrgVol;
  int   irow_src, icol_src, islc_src; /* integer row, col, slc in source */
  float frow_src, fcol_src, fslc_src; /* float row, col, slc in source */
  float srcval, Tx, Ty, Tz;
  int frm, FreeQsrc=0;
  int vtx,nhits;
  float *valvect;
  double rval;

  if(Qsrc == NULL){
    Qsrc = MRIxfmCRS2XYZtkreg(SrcVol);
    Qsrc = MatrixInverse(Qsrc,Qsrc);
    FreeQsrc = 1;
  }

  /* compute the transforms */
  QFWDsrc = ComputeQFWD(Qsrc,Fsrc,Wsrc,Dsrc,NULL);
  if (Gdiag_no >= 0)
  {
    printf("QFWDsrc: vol2surf: ------------------------------\n");
    MatrixPrint(stdout,QFWDsrc);
    printf("--------------------------------------------------\n");
  }

  /* preallocate the row-col-slc vectors */
  Scrs = MatrixAlloc(4,1,MATRIX_REAL);
  Txyz = MatrixAlloc(4,1,MATRIX_REAL);
  Txyz->rptr[3+1][0+1] = 1.0;

  /* allocate a "volume" to hold the output */
  TrgVol = MRIallocSequence(TrgSurf->nvertices,1,1,MRI_FLOAT,SrcVol->nframes);
  if (TrgVol == NULL) return(NULL);
  MRIcopyHeader(SrcVol,TrgVol);
  // Dims here are meaningless, but setting to 1 means "volume" will be 
  // number of vertices.
  TrgVol->xsize = 1; 
  TrgVol->ysize = 1;
  TrgVol->zsize = 1;

  /* Zero the source hit volume */
  if (SrcHitVol != NULL)
  {
    MRIconst(SrcHitVol->width,SrcHitVol->height,SrcHitVol->depth,
             1,0,SrcHitVol);
  }

  srcval = 0;
  valvect = (float *) calloc(sizeof(float),SrcVol->nframes);
  nhits = 0;
  /*--- loop through each vertex ---*/
  for (vtx = 0; vtx < TrgSurf->nvertices; vtx+=nskip)
  {

    if (ProjFrac != 0.0)
      if (ProjDistFlag)
        ProjNormDist(&Tx,&Ty,&Tz,TrgSurf,vtx,ProjFrac);
      else
        ProjNormFracThick(&Tx,&Ty,&Tz,TrgSurf,vtx,ProjFrac);
    else
    {
      Tx = TrgSurf->vertices[vtx].x;
      Ty = TrgSurf->vertices[vtx].y;
      Tz = TrgSurf->vertices[vtx].z;
    }

    /* Load the Target xyz vector */
    Txyz->rptr[0+1][0+1] = Tx;
    Txyz->rptr[1+1][0+1] = Ty;
    Txyz->rptr[2+1][0+1] = Tz;

    /* Compute the corresponding Source col-row-slc vector */
    MatrixMultiply(QFWDsrc,Txyz,Scrs);
    fcol_src = Scrs->rptr[1][1];
    frow_src = Scrs->rptr[2][1];
    fslc_src = Scrs->rptr[3][1];

    /* nearest neighbor */
    switch (float2int)
    {
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
      irow_src = (int) ceil(frow_src);
      islc_src = (int)floor(fslc_src);
      break;
    default:
      fprintf(stderr,"vol2surf_linear(): unrecoginized float2int code %d\n",
              float2int);
      MRIfree(&TrgVol);
      return(NULL);
      break;
    }

    /* check that the point is in the bounds of the volume */
    if (irow_src < 0 || irow_src >= SrcVol->height ||
        icol_src < 0 || icol_src >= SrcVol->width  ||
        islc_src < 0 || islc_src >= SrcVol->depth ) continue;

    if (Gdiag_no == vtx)
    {
      printf("diag -----------------------------\n");
      printf("vtx = %d  %g %g %g\n",vtx,Tx,Ty,Tz);
      printf("fCRS  %g %g %g\n",Scrs->rptr[1][1],
             Scrs->rptr[2][1],Scrs->rptr[3][1]);
      printf("CRS  %d %d %d\n",icol_src,irow_src,islc_src);
    }


    /* only gets here if it is in bounds */
    nhits ++;

    /* Assign output volume values */
    if (InterpMethod == SAMPLE_TRILINEAR)
    {
      MRIsampleSeqVolume(SrcVol, fcol_src, frow_src, fslc_src,
                         valvect, 0, SrcVol->nframes-1) ;
      if (Gdiag_no == vtx)
        printf("val = %f\n", valvect[0]) ;
      for (frm = 0; frm < SrcVol->nframes; frm++)
        MRIFseq_vox(TrgVol,vtx,0,0,frm) = valvect[frm];
    }
    else
    {
      for (frm = 0; frm < SrcVol->nframes; frm++)
      {
        switch (InterpMethod)
        {
        case SAMPLE_NEAREST:
          srcval = MRIgetVoxVal(SrcVol,icol_src,irow_src,islc_src,frm);
          break ;
          srcval = MRIgetVoxVal(SrcVol,icol_src,irow_src,islc_src,frm);
          break ;
        case SAMPLE_SINC:      /* no multi-frame */
          MRIsincSampleVolume(SrcVol, fcol_src, frow_src, fslc_src, 5, &rval) ;
          srcval = rval;
          break ;
        } //switch
        MRIFseq_vox(TrgVol,vtx,0,0,frm) = srcval;
        if (Gdiag_no == vtx)
          printf("val[%d] = %f\n", frm, srcval) ;
      } // for
    }// else
    if (SrcHitVol != NULL)
      MRIFseq_vox(SrcHitVol,icol_src,irow_src,islc_src,0)++;
  }

  MatrixFree(&QFWDsrc);
  MatrixFree(&Scrs);
  MatrixFree(&Txyz);
  free(valvect);
  if(FreeQsrc) MatrixFree(&Qsrc);

  //printf("vol2surf_linear: nhits = %d/%d\n",nhits,TrgSurf->nvertices);

  return(TrgVol);
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
MRI *MRISapplyReg(MRI *SrcSurfVals, MRI_SURFACE **SurfReg, int nsurfs,
		  int ReverseMapFlag, int DoJac, int UseHash)
{
  MRI *TrgSurfVals = NULL;
  MRI_SURFACE *SrcSurfReg, *TrgSurfReg;
  int svtx=0, tvtx, tvtxN, svtxN=0, f, n, nrevhits,nSrcLost;
  int npairs, kS, kT, nhits, nunmapped;
  VERTEX *v;
  float dmin;
  MHT **Hash=NULL;
  MRI *SrcHits, *TrgHits;

  npairs = nsurfs/2;
  printf("MRISapplyReg: nsurfs = %d, revmap=%d, jac=%d,  hash=%d\n",
	 nsurfs,ReverseMapFlag,DoJac,UseHash);

  SrcSurfReg = SurfReg[0];
  TrgSurfReg = SurfReg[nsurfs-1];

  /* check dimension consistency */
  if (SrcSurfVals->width != SrcSurfReg->nvertices){
    printf("MRISapplyReg: Vals and Reg dimension mismatch\n");
    printf("nVals = %d, nReg %d\n",SrcSurfVals->width,
            SrcSurfReg->nvertices);
    return(NULL);
  }
  for(n=0; n < npairs-1; n++){
    kS = 2*n+1;
    kT = kS+1;
    if(SurfReg[kT]->nvertices != SurfReg[kS]->nvertices){
      printf("MRISapplyReg: Reg dimension mismatch %d, %d\n",kT,kS);
      printf("targ = %d, next source = %d\n",SurfReg[kT]->nvertices,SurfReg[kS]->nvertices);
    return(NULL);
    }
  }

  /* allocate a "volume" to hold the output */
  TrgSurfVals = MRIallocSequence(TrgSurfReg->nvertices,1,1,
                                 MRI_FLOAT,SrcSurfVals->nframes);
  if(TrgSurfVals == NULL) return(NULL);
  MRIcopyHeader(SrcSurfVals,TrgSurfVals);

  /* number of source vertices mapped to each target vertex */
  TrgHits = MRIallocSequence(TrgSurfReg->nvertices,1,1,MRI_FLOAT,1);
  if(TrgHits == NULL) return(NULL);
  MRIcopyHeader(SrcSurfVals,TrgHits);

  /* number of target vertices mapped to by each source vertex */
  SrcHits = MRIallocSequence(SrcSurfReg->nvertices,1,1,MRI_FLOAT,1);
  if (SrcHits == NULL) return(NULL);
  MRIcopyHeader(SrcSurfVals,SrcHits);

  if(UseHash){
    printf("MRISapplyReg: building hash tables (res=16).\n");
    Hash = (MHT **)calloc(sizeof(MHT*),nsurfs);
    for(n=0; n < nsurfs; n++){
      Hash[n] = (MHT *)calloc(sizeof(MHT),1);
      Hash[n] = MHTfillVertexTableRes(SurfReg[n], NULL,CURRENT_VERTICES,16);
    }
  }

  if(DoJac){
    // If using jacobian correction, get a list of the number of times
    // that a give source vertex gets sampled.
    for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++){
      // Compute the source vertex that corresponds to this target vertex
      tvtxN = tvtx;
      for(n=npairs-1; n >= 0; n--){
	kS = 2*n;
	kT = kS + 1;
	v = &(SurfReg[kT]->vertices[tvtxN]);
	/* find closest source vertex */
	if (UseHash) svtx = MHTfindClosestVertexNo(Hash[kS],SurfReg[kS],v,&dmin);
	else         svtx = MRISfindClosestVertex(SurfReg[kS],v->x,v->y,v->z,&dmin);
	tvtxN = svtx;
      }
      /* update the number of hits and distance */
      MRIFseq_vox((SrcHits),svtx,0,0,0) ++;
      MRIFseq_vox((TrgHits),tvtx,0,0,0) ++;
    }
  }

  /* Go through the forwad loop (finding closest srcvtx to each trgvtx).
  This maps each target vertex to a source vertex */
  printf("MRISapplyReg: Forward Loop (%d)\n",TrgSurfReg->nvertices);
  nunmapped = 0;
  for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++){
    if (!UseHash){
      if (tvtx%100 == 0){printf("%5d ",tvtx);fflush(stdout);}
      if(tvtx%1000 == 999){printf("\n");fflush(stdout);}
    }

    // Compute the source vertex that corresponds to this target vertex
    tvtxN = tvtx;
    for(n=npairs-1; n >= 0; n--){
      kS = 2*n;
      kT = kS + 1;
      //printf("%5d %5d %d %d %d\n",tvtx,tvtxN,n,kS,kT);
      v = &(SurfReg[kT]->vertices[tvtxN]);
      /* find closest source vertex */
      if(UseHash) svtx = MHTfindClosestVertexNo(Hash[kS],SurfReg[kS],v,&dmin);
      if(!UseHash || svtx < 0){
	if(svtx < 0) printf("Target vertex %d of pair %d unmapped in hash, using brute force\n",tvtxN,n);
	svtx = MRISfindClosestVertex(SurfReg[kS],v->x,v->y,v->z,&dmin);
      }
      tvtxN = svtx;
    }

    if(! DoJac){
      /* update the number of hits */
      MRIFseq_vox(SrcHits,svtx,0,0,0) ++;
      MRIFseq_vox(TrgHits,tvtx,0,0,0) ++;
      nhits = 1;
    } 
    else nhits = MRIgetVoxVal(SrcHits,svtx,0,0,0);

    /* accumulate mapped values for each frame */
    for (f=0; f < SrcSurfVals->nframes; f++)
      MRIFseq_vox(TrgSurfVals,tvtx,0,0,f) += (MRIFseq_vox(SrcSurfVals,svtx,0,0,f)/nhits);
  }

  /*---------------------------------------------------------------
  Go through the reverse loop (finding closest trgvtx to each srcvtx
  unmapped by the forward loop). This assures that each source vertex
  is represented in the map */
  if(ReverseMapFlag){
    printf("MRISapplyReg: Reverse Loop (%d)\n",SrcSurfReg->nvertices);
    nrevhits = 0;
    for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++) {
      if(MRIFseq_vox((SrcHits),svtx,0,0,0) != 0) continue;
      nrevhits ++;

      // Compute the target vertex that corresponds to this source vertex
      svtxN = svtx;
      for(n=0; n < npairs; n++){
	kS = 2*n;
	kT = kS + 1;
	//printf("%5d %5d %d %d %d\n",svtx,svtxN,n,kS,kT);
	v = &(SurfReg[kS]->vertices[svtxN]);
	/* find closest target vertex */
	if (UseHash) tvtx = MHTfindClosestVertexNo(Hash[kT],SurfReg[kT],v,&dmin);
	if(!UseHash || tvtx < 0){
	  if(tvtx < 0) printf("Source vertex %d of pair %d unmapped in hash, using brute force\n",svtxN,n);
	  tvtx = MRISfindClosestVertex(SurfReg[kT],v->x,v->y,v->z,&dmin);
	}
	svtxN = tvtx;
      }
      
      /* update the number of hits */
      MRIFseq_vox((SrcHits),svtx,0,0,0) ++;
      MRIFseq_vox((TrgHits),tvtx,0,0,0) ++;
      /* accumulate mapped values for each frame */
      for (f=0; f < SrcSurfVals->nframes; f++)
	MRIFseq_vox(TrgSurfVals,tvtx,0,0,f) += MRIFseq_vox(SrcSurfVals,svtx,0,0,f);
    }
    printf("  Reverse Loop had %d hits\n",nrevhits);
  }
  
  /*---------------------------------------------------------------
  Finally, divide the value at each target vertex by the number
  of source vertices mapping into it */
  if(! DoJac){
    printf("MRISapplyReg: Dividing by number of hits (%d)\n",TrgSurfReg->nvertices);
    for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++) {
      n = MRIFseq_vox((TrgHits),tvtx,0,0,0);
      if (n > 1){
	for (f=0; f < SrcSurfVals->nframes; f++)
	  MRIFseq_vox(TrgSurfVals,tvtx,0,0,f) /= n;
      }
    }
  }

  /* Count lost sources */
  nSrcLost = 0;
  for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++){
    n = MRIFseq_vox((SrcHits),svtx,0,0,0);
    if (n == 0)nSrcLost ++;
  }
  printf("MRISapplyReg: nSrcLost = %d\n",nSrcLost);

  MRIfree(&SrcHits);
  MRIfree(&TrgHits);
  if (UseHash) for(n=0; n < nsurfs; n++) MHTfree(&Hash[n]);
  return(TrgSurfVals);
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
MRI *surf2surf_nnfr(MRI *SrcSurfVals, MRI_SURFACE *SrcSurfReg,
                    MRI_SURFACE *TrgSurfReg, MRI **SrcHits,
                    MRI **SrcDist, MRI **TrgHits, MRI **TrgDist,
                    int ReverseMapFlag, int UseHash)
{
  MRI *TrgSurfVals = NULL;
  int svtx, tvtx, f,n, nrevhits,nSrcLost;
  VERTEX *v;
  MHT *SrcHash, *TrgHash;
  float dmin;
  extern char *ResampleVtxMapFile;
  FILE *fp = NULL;

  /* check dimension consistency */
  if (SrcSurfVals->width != SrcSurfReg->nvertices)
  {
    fprintf(stderr,"surf2surf_nnfr(): Vals and Reg dimension mismatch\n");
    fprintf(stderr,"nVals = %d, nReg %d\n",SrcSurfVals->width,
            SrcSurfReg->nvertices);
    return(NULL);
  }

  /* allocate a "volume" to hold the output */
  TrgSurfVals = MRIallocSequence(TrgSurfReg->nvertices,1,1,
                                 MRI_FLOAT,SrcSurfVals->nframes);
  if(TrgSurfVals == NULL) return(NULL);
  MRIcopyHeader(SrcSurfVals,TrgSurfVals);

  /* number of source vertices mapped to each target vertex */
  *TrgHits = MRIallocSequence(TrgSurfReg->nvertices,1,1,MRI_FLOAT,1);
  if(*TrgHits == NULL) return(NULL);
  MRIcopyHeader(SrcSurfVals,*TrgHits);

  /* Average distance of a target vertex from its source vertices  */
  *TrgDist = MRIallocSequence(TrgSurfReg->nvertices,1,1,MRI_FLOAT,1);
  if (*TrgDist == NULL) return(NULL);
  MRIcopyHeader(SrcSurfVals,*TrgDist);

  /* number of target vertices mapped to by each source vertex */
  *SrcHits = MRIallocSequence(SrcSurfReg->nvertices,1,1,MRI_FLOAT,1);
  if (*SrcHits == NULL) return(NULL);
  MRIcopyHeader(SrcSurfVals,*SrcHits);

  /* Average distance of a source vertex from its target vertices  */
  *SrcDist = MRIallocSequence(SrcSurfReg->nvertices,1,1,MRI_FLOAT,1);
  if (*SrcDist == NULL) return(NULL);
  MRIcopyHeader(SrcSurfVals,*SrcDist);

  /* build hash tables */
  if (UseHash)
  {
    printf("surf2surf_nnfr: building source hash (res=16).\n");
    SrcHash = MHTfillVertexTableRes(SrcSurfReg, NULL,CURRENT_VERTICES,16);
  }

  /* Open vertex map file */
  if (ResampleVtxMapFile != NULL)
  {
    fp = fopen(ResampleVtxMapFile,"w");
    if (fp == NULL)
    {
      printf("ERROR: could not open %s\n",ResampleVtxMapFile);
      exit(1);
    }
  }

  /*---------------------------------------------------------------
    Go through the forwad loop (finding closest srcvtx to each trgvtx).
    This maps each target vertex to a source vertex */
  printf("Surf2Surf: Forward Loop (%d)\n",TrgSurfReg->nvertices);
  for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++)
  {
    if (!UseHash){
      if (tvtx%100 == 0) {
        printf("%5d ",tvtx);
        fflush(stdout);
      }
      if (tvtx%1000 == 999){
        printf("\n");
        fflush(stdout);
      }
    }
    /* find closest source vertex */
    v = &(TrgSurfReg->vertices[tvtx]);
    if (UseHash) svtx = MHTfindClosestVertexNo(SrcHash,SrcSurfReg,v,&dmin);
    else         svtx = MRISfindClosestVertex(SrcSurfReg,v->x,v->y,v->z,&dmin);

    /* hash table failed, so use brute force */
    if(svtx < 0) svtx = MRISfindClosestVertex(SrcSurfReg,v->x,v->y,v->z,&dmin);

    /* update the number of hits and distance */
    MRIFseq_vox((*SrcHits),svtx,0,0,0) ++;
    MRIFseq_vox((*TrgHits),tvtx,0,0,0) ++;
    MRIFseq_vox((*SrcDist),svtx,0,0,0) += dmin;
    MRIFseq_vox((*TrgDist),tvtx,0,0,0) += dmin;

    /* accumulate mapped values for each frame */
    for (f=0; f < SrcSurfVals->nframes; f++)
      MRIFseq_vox(TrgSurfVals,tvtx,0,0,f) +=
        MRIFseq_vox(SrcSurfVals,svtx,0,0,f);

    if (ResampleVtxMapFile != NULL)
    {
      fprintf(fp,"%6d  (%6.1f,%6.1f,%6.1f)   ",tvtx,v->x,v->y,v->z);
      v = &(SrcSurfReg->vertices[svtx]);
      fprintf(fp,"%6d  (%6.1f,%6.1f,%6.1f)    %5.4f\n",svtx,v->x,v->y,v->z,dmin);
      fflush(fp);
    }
  }
  printf("\n");
  if(UseHash) MHTfree(&SrcHash);

  if (ResampleVtxMapFile != NULL) fclose(fp);

  /*---------------------------------------------------------------
    Go through the reverse loop (finding closest trgvtx to each srcvtx
    unmapped by the forward loop). This assures that each source vertex
    is represented in the map */
  if (ReverseMapFlag)
  {
    if (UseHash)
    {
      MHTfree(&SrcHash);
      printf("surf2surf_nnfr: building target hash (res=16).\n");
      TrgHash = MHTfillVertexTableRes(TrgSurfReg, NULL,CURRENT_VERTICES,16);
    }
    printf("Surf2Surf: Reverse Loop (%d)\n",SrcSurfReg->nvertices);
    nrevhits = 0;
    for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++)
    {
      if (MRIFseq_vox((*SrcHits),svtx,0,0,0) == 0)
      {
        nrevhits ++;
        /* find closest target vertex */
        v = &(SrcSurfReg->vertices[svtx]);
        if (UseHash) tvtx = MHTfindClosestVertexNo(TrgHash,TrgSurfReg,v,&dmin);
        else         tvtx = MRISfindClosestVertex(TrgSurfReg,v->x,v->y,v->z,&dmin);
	/* Hash table failed, so use brute force */
        if(tvtx < 0) tvtx = MRISfindClosestVertex(TrgSurfReg,v->x,v->y,v->z,&dmin);

        /* update the number of hits and distance */
        MRIFseq_vox((*SrcHits),svtx,0,0,0) ++;
        MRIFseq_vox((*TrgHits),tvtx,0,0,0) ++;
        MRIFseq_vox((*SrcDist),svtx,0,0,0) += dmin;
        MRIFseq_vox((*TrgDist),tvtx,0,0,0) += dmin;
        /* accumulate mapped values for each frame */
        for (f=0; f < SrcSurfVals->nframes; f++)
          MRIFseq_vox(TrgSurfVals,tvtx,0,0,f) +=
            MRIFseq_vox(SrcSurfVals,svtx,0,0,f);
      }
    }
    if(UseHash) MHTfree(&TrgHash);
    printf("Reverse Loop had %d hits\n",nrevhits);
  }

  /*---------------------------------------------------------------
    Finally, divide the value at each target vertex by the number
    of source vertices mapping into it */
  printf("Surf2Surf: Dividing by number of hits (%d)\n",TrgSurfReg->nvertices);
  for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++)
  {
    n = MRIFseq_vox((*TrgHits),tvtx,0,0,0);
    if (n > 1)
    {
      MRIFseq_vox((*TrgDist),tvtx,0,0,0) /= n; /* average distances */
      for (f=0; f < SrcSurfVals->nframes; f++)
        MRIFseq_vox(TrgSurfVals,tvtx,0,0,f) /= n;
    }
  }
  /* go through the source loop to average the distance */
  nSrcLost = 0;
  for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++)
  {
    n = MRIFseq_vox((*SrcHits),svtx,0,0,0);
    if (n == 1) continue;
    if (n > 1) MRIFseq_vox((*SrcDist),svtx,0,0,0) /= n;
    else
    {/* unmapped */
      MRIFseq_vox((*SrcDist),svtx,0,0,0) = -10;
      nSrcLost ++;
    }
  }

  printf("INFO: nSrcLost = %d\n",nSrcLost);

  return(TrgSurfVals);
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
MRI *surf2surf_nnfr_jac(MRI *SrcSurfVals, MRI_SURFACE *SrcSurfReg,
                        MRI_SURFACE *TrgSurfReg, MRI **SrcHits,
                        MRI **SrcDist, MRI **TrgHits, MRI **TrgDist,
                        int ReverseMapFlag, int UseHash)
{
  MRI *TrgSurfVals = NULL;
  int svtx, tvtx, f,n, nunmapped, nrevhits,nSrcLost,nhits;
  VERTEX *v;
  MHT *SrcHash, *TrgHash;
  float dmin,srcval;

  /* check dimension consistency */
  if (SrcSurfVals->width != SrcSurfReg->nvertices)
  {
    fprintf(stderr,"surf2surf_nnfr(): Vals and Reg dimension mismatch\n");
    fprintf(stderr,"nVals = %d, nReg %d\n",SrcSurfVals->width,
            SrcSurfReg->nvertices);
    return(NULL);
  }

  /* allocate a "volume" to hold the output */
  TrgSurfVals = MRIallocSequence(TrgSurfReg->nvertices,1,1,
                                 MRI_FLOAT,SrcSurfVals->nframes);
  if (TrgSurfVals == NULL) return(NULL);
  MRIcopyHeader(SrcSurfVals,TrgSurfVals);

  /* number of source vertices mapped to each target vertex */
  *TrgHits = MRIallocSequence(TrgSurfReg->nvertices,1,1,MRI_FLOAT,1);
  if (*TrgHits == NULL) return(NULL);

  /* Average distance of a target vertex from its source vertices  */
  *TrgDist = MRIallocSequence(TrgSurfReg->nvertices,1,1,MRI_FLOAT,1);
  if (*TrgDist == NULL) return(NULL);

  /* number of target vertices mapped to by each source vertex */
  *SrcHits = MRIallocSequence(SrcSurfReg->nvertices,1,1,MRI_FLOAT,1);
  if (*SrcHits == NULL) return(NULL);

  /* Average distance of a source vertex from its target vertices  */
  *SrcDist = MRIallocSequence(SrcSurfReg->nvertices,1,1,MRI_FLOAT,1);
  if (*SrcDist == NULL) return(NULL);

  /* build hash tables */
  if (UseHash){
    printf("surf2surf_nnfr_jac: building source hash (res=16).\n");
    SrcHash = MHTfillVertexTableRes(SrcSurfReg, NULL,CURRENT_VERTICES,16);
  }

  // First forward loop just counts the number of hits for each src
  printf("Surf2SurfJac: 1st Forward Loop (%d)\n",TrgSurfReg->nvertices);
  nunmapped = 0;
  for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++){
    /* find closest source vertex */
    v = &(TrgSurfReg->vertices[tvtx]);
    if(UseHash) svtx = MHTfindClosestVertexNo(SrcHash,SrcSurfReg,v,&dmin);
    else        svtx = MRISfindClosestVertex(SrcSurfReg,v->x,v->y,v->z,&dmin);
    /* hash table failed, so use brute force */
    if(svtx < 0) svtx = MRISfindClosestVertex(SrcSurfReg,v->x,v->y,v->z,&dmin);

    /* update the number of hits and distance */
    MRIFseq_vox((*SrcHits),svtx,0,0,0) ++;  // This is what this loop is for
    MRIFseq_vox((*TrgHits),tvtx,0,0,0) ++;
    MRIFseq_vox((*SrcDist),svtx,0,0,0) += dmin;
    MRIFseq_vox((*TrgDist),tvtx,0,0,0) += dmin;
  }

  // Second forward loop accumulates
  printf("Surf2SurfJac: 2nd Forward Loop (%d)\n",TrgSurfReg->nvertices);
  for (tvtx = 0; tvtx < TrgSurfReg->nvertices; tvtx++) {
    /* find closest source vertex */
    v = &(TrgSurfReg->vertices[tvtx]);
    if(UseHash) svtx = MHTfindClosestVertexNo(SrcHash,SrcSurfReg,v,&dmin);
    else        svtx = MRISfindClosestVertex(SrcSurfReg,v->x,v->y,v->z,&dmin);
    /* hash table failed, so use bruce force */
    if(svtx < 0) svtx = MRISfindClosestVertex(SrcSurfReg,v->x,v->y,v->z,&dmin);

    nhits = MRIFseq_vox((*SrcHits),svtx,0,0,0);
    /* Now accumulate mapped values for each frame */
    for (f=0; f < SrcSurfVals->nframes; f++){
      srcval = MRIFseq_vox(SrcSurfVals,svtx,0,0,f);
      srcval /= nhits;// divide by number of hits
      MRIFseq_vox(TrgSurfVals,tvtx,0,0,f) += srcval;
    }
  }

  /*---------------------------------------------------------------
    Go through the reverse loop (finding closest trgvtx to each srcvtx
    unmapped by the forward loop). This assures that each source vertex
    is represented in the map */
  if (ReverseMapFlag)
  {
    if (UseHash)
    {
      MHTfree(&SrcHash);
      printf("surf2surf_nnfr: building target hash (res=16).\n");
      TrgHash = MHTfillVertexTableRes(TrgSurfReg, NULL,CURRENT_VERTICES,16);
    }
    printf("Surf2SurfJac: Reverse Loop (%d)\n",SrcSurfReg->nvertices);
    nrevhits = 0;
    for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++)
    {
      if (MRIFseq_vox((*SrcHits),svtx,0,0,0) == 0)
      {
        nrevhits ++;
        /* find closest target vertex */
        v = &(SrcSurfReg->vertices[svtx]);
        if (UseHash) tvtx = MHTfindClosestVertexNo(TrgHash,TrgSurfReg,v,&dmin);
        else        tvtx = MRISfindClosestVertex(TrgSurfReg,v->x,v->y,v->z,&dmin);
	/* Hash table failed, so use brute force */
        if(tvtx < 0) tvtx = MRISfindClosestVertex(TrgSurfReg,v->x,v->y,v->z,&dmin);
        /* update the number of hits and distance */
        MRIFseq_vox((*SrcHits),svtx,0,0,0) ++;
        MRIFseq_vox((*TrgHits),tvtx,0,0,0) ++;
        MRIFseq_vox((*SrcDist),svtx,0,0,0) += dmin;
        MRIFseq_vox((*TrgDist),tvtx,0,0,0) += dmin;
        /* accumulate mapped values for each frame */
        for (f=0; f < SrcSurfVals->nframes; f++)
          MRIFseq_vox(TrgSurfVals,tvtx,0,0,f) +=
            MRIFseq_vox(SrcSurfVals,svtx,0,0,f);
      }
    }
    if (UseHash)MHTfree(&TrgHash);
    printf("Reverse Loop had %d hits\n",nrevhits);
  }

  // Do NOT normalize target vertices with multiple src vertices

  /* go through the source loop to average the distance */
  nSrcLost = 0;
  for (svtx = 0; svtx < SrcSurfReg->nvertices; svtx++)  {
    n = MRIFseq_vox((*SrcHits),svtx,0,0,0);
    if (n == 1) continue;
    if (n > 1) MRIFseq_vox((*SrcDist),svtx,0,0,0) /= n;
    else {/* unmapped */
      MRIFseq_vox((*SrcDist),svtx,0,0,0) = -10;
      nSrcLost ++;
    }
  }
  printf("INFO: nSrcLost = %d\n",nSrcLost);
  printf("surf2surf_nnfr_jac() done\n");
  fflush(stdout);

  return(TrgSurfVals);
}

/*-------------------------------------------------------------
  crs2ind() -- returns linear index into a volume stored by column,
  row, slice.
  --------------------------------------------------------------------*/
int crs2ind(int *ind, int c, int r, int s, int ncols, int nrows, int nslcs )
{
  if (c < 0 || c >= ncols)
  {
    fprintf(stderr,"crs2ind: col %d out of bounds (0,%d)\n",c,ncols-1);
    return(1);
  }
  if (r < 0 || r >= nrows)
  {
    fprintf(stderr,"crs2ind: row %d out of bounds (0,%d)\n",r,nrows-1);
    return(1);
  }
  if (s < 0 || s >= nslcs)
  {
    fprintf(stderr,"crs2ind: slc %d out of bounds (0,%d)\n",c,nslcs-1);
    return(1);
  }

  *ind = c + r * ncols + s * nrows*ncols ;

  return(0);
}
/*-------------------------------------------------------------
  ind2crs() -- returns the column, row, and slice of an element in a
  volume given the index of the element in the volume assuming that
  the elements are stored by column, row, then slice.
  --------------------------------------------------------------------*/
int ind2crs(int *c, int *r, int *s, int ind, int ncols, int nrows, int nslcs)
{
  int i = ind, ntot, nrowcols;

  nrowcols = nrows*ncols;
  ntot = nrowcols*nslcs;
  if (i < 0 || i >= ntot)
  {
    fprintf(stderr,"ind2crs: index %d out of bounds (0,%d)\n",i,ntot-1);
    return(1);
  }

  *s  = (int ) (i/nrowcols);
  i  = i - *s * (nrowcols);

  *r  = (int ) (i/ncols);
  i  = i - *r * ncols;

  *c  = (int ) i;

  return(0);
}
/*-------------------------------------------------------------*/
int MatrixFill(MATRIX *M, float val)
{
  int r,c;

  for (r=0;r < M->rows; r++)
  {
    for (c=0; c < M->cols; c++)
    {
      M->rptr[r+1][c+1] = val;
    }
  }
  return(0);
}
/*-------------------------------------------------------------*/
MATRIX * RandMatrix(int rows, int cols)
{
  int r,c;
  MATRIX *M;

  M = MatrixAlloc(rows, cols, MATRIX_REAL);

  for (r=0;r < M->rows; r++)
  {
    for (c=0; c < M->cols; c++)
    {
      M->rptr[r+1][c+1] = drand48();
    }
  }
  return(M);
  ;
}
/*-------------------------------------------------------------*/
MATRIX * ConstMatrix(int rows, int cols, float val)
{
  int r,c;
  MATRIX *M;

  M = MatrixAlloc(rows, cols, MATRIX_REAL);

  for (r=0;r < M->rows; r++)
  {
    for (c=0; c < M->cols; c++)
    {
      M->rptr[r+1][c+1] = val;
    }
  }
  return(M);
  ;
}


/*-------------------------------------------------------------------
  MRIsurf2Vol() - the purpose of this function is to resample values
  from a surface to a volume. Most of the real work is done in the
  creation of map, which has the same size as vol; each voxel in map
  has the vertex number of the corresponding point on the surface. If
  there is no corresponding point, then the value will be -1. See
  MRImapSurf2VolClosest().  MRI vol must have already been allocated.
  Returns -1 on error, otherwise returns the number of voxels in
  the volume that have been hit by the surface.
  ----------------------------------------------------------------*/
int MRIsurf2Vol(MRI *surfvals, MRI *vol, MRI *map)
{
  int vtx, c, r, s, f, nhits;
  float val ;

  if (vol->width  != map->width ||
      vol->height != map->height||
      vol->depth  != map->depth )
  {
    printf("ERROR: vol and map dimensions differ\n");
    return(-1);
  }

  if (surfvals->nframes != vol->nframes)
  {
    printf("ERROR: surfvals and vol have different number of frames\n");
    return(-1);
  }

  nhits = 0;
  for (c=0; c < vol->width; c++)
  {
    for (r=0; r < vol->height; r++)
    {
      for (s=0; s < vol->depth; s++)
      {
        vtx = MRIIseq_vox(map,c,r,s,0);
        for (f = 0; f < vol->nframes; f++)
        {
          if (vtx < 0)
          {
#if 0
            MRIsetVoxVal(vol,c,r,s,f,0.0);
#endif
            continue;
          }
          val = MRIgetVoxVal(surfvals,vtx,0,0,f);
          MRIsetVoxVal(vol,c,r,s,f,val);
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
  return(nhits);
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
  float xvtx,yvtx,zvtx;
  float d2, d2min;
  MATRIX *xyzvtx, *icrs, *fcrs, *xyzcrs;
  //MATRIX *Tvol;

  /* Alloc map - for a given voxel, holds the number of the
     closest vertex */
  map = MRIalloc(vol->width, vol->height, vol->depth, MRI_INT);
  if (map == NULL)
  {
    printf("ERROR: MRImapSurf2VolClosest: could not alloc vtx map\n");
    return(NULL);
  }
  /* Alloc dist - for a given voxel, holds the distance of the
     closest vertex to the voxel. */
  dist2 = MRIalloc(vol->width, vol->height, vol->depth, MRI_FLOAT);
  if (dist2 == NULL)
  {
    printf("ERROR: MRImapSurf2VolClosest: could not alloc dist map\n");
    return(NULL);
  }

  /* Initially set all the voxels in the map to -1 to mark that
     they have not been hit by a surface vertex */
  MRIvalueFill(map,-1);

  /* Intialize matrix stuff */
  //Tvol = MRIxfmCRS2XYZ(vol, 0); /* converts crs to xyz in vol */
  xyzvtx = MatrixAlloc(4,1,MATRIX_REAL);
  fcrs = MatrixAlloc(4,1,MATRIX_REAL);   /* vertex row, col, slice */
  icrs = MatrixAlloc(4,1,MATRIX_REAL);   /* voxel row, col, slice */
  xyzcrs = MatrixAlloc(4,1,MATRIX_REAL); /* voxel xyz */

  /* Set the 4th item to 1 */
  xyzvtx->rptr[4][1] = 1.0;
  fcrs->rptr[4][1]   = 1.0;
  icrs->rptr[4][1]   = 1.0;
  xyzcrs->rptr[4][1] = 1.0;

  /*----------- Loop through vertices --------------*/
  for (vtx = 0; vtx < surf->nvertices; vtx++)
  {

    if (projfrac == 0)
    {
      xvtx = surf->vertices[vtx].x;
      yvtx = surf->vertices[vtx].y;
      zvtx = surf->vertices[vtx].z;
    }
    else
    {
      /* Get the xyz of the vertex as projected along the normal a
      distance equal to a fraction of the cortical thickness at
      that point. */
      ProjNormFracThick(&xvtx,&yvtx,&zvtx,surf,vtx,projfrac);
    }

    /* load vertex xyz values into a vector */
    xyzvtx->rptr[1][1] = xvtx;
    xyzvtx->rptr[2][1] = yvtx;
    xyzvtx->rptr[3][1] = zvtx;
    /* [4] is already 1 */

    /* Compute the CRS in volume space */
    /* fcrs = Qa2v * xyzvtx */
    MatrixMultiply(Qa2v,xyzvtx,fcrs);

    /* Round CRS to nearest integer */
    c = nint(fcrs->rptr[1][1]);
    r = nint(fcrs->rptr[2][1]);
    s = nint(fcrs->rptr[3][1]);

    if (Gdiag_no > 0 && vtx == 30756)
    {
      printf("diag -----------------------------\n");
      printf("vtx = %d  %g %g %g\n",vtx,xvtx,yvtx,zvtx);
      printf("fCRS  %g %g %g\n",fcrs->rptr[1][1],
             fcrs->rptr[2][1],fcrs->rptr[3][1]);
      printf("CRS  %d %d %d\n",c,r,s);
    }

    /* Check that it is in the volume */
    if (c < 0 || c >= vol->width)  continue;
    if (r < 0 || r >= vol->height) continue;
    if (s < 0 || s >= vol->depth)  continue;

    /* Load rounded CRS into vector */
    icrs->rptr[1][1] = c;
    icrs->rptr[2][1] = r;
    icrs->rptr[3][1] = s;
    /* [4] is already 1 */

    /* Compute XYZ of rounded CRS in volume space. This
    is the XYZ of the center of the voxel*/
    /* xyzcrs = Tvol * icrs */
    //MatrixMultiply(Tvol,icrs,xyzcrs);

    /* Compute the distance between the voxel and the
       vertex (actually distance squared) */
    d2 = SQR( vol->xsize*(c - icrs->rptr[1][1]) ) +
         SQR( vol->ysize*(r - icrs->rptr[2][1]) ) +
         SQR( vol->zsize*(s - icrs->rptr[3][1]) );

    /* Check whether this voxel has been hit. If not
       just set its map vertex and dist to current. */
    vtxmin = MRIIseq_vox(map,c,r,s,0);
    if (vtxmin < 0)
    {
      MRIIseq_vox(map,c,r,s,0)   = vtx;
      MRIFseq_vox(dist2,c,r,s,0) = d2;
      continue;
    }

    /* Check whether this vertex is closer */
    d2min = MRIFseq_vox(dist2,c,r,s,0);
    if (d2min > d2)
    {
      MRIIseq_vox(map,c,r,s,0)   = vtx;
      MRIFseq_vox(dist2,c,r,s,0) = d2;
    }

  } /* end loop over vertices */

  //MatrixFree(&Tvol);
  MatrixFree(&xyzvtx);
  MatrixFree(&icrs);
  MatrixFree(&fcrs);
  MatrixFree(&xyzcrs);
  MRIfree(&dist2);
  return(map);
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
MRI *MRIaseg2vol(MRI *aseg, MATRIX *tkR, MRI *voltemp,
                 double fthresh, MRI **pvolhit)
{
  int Na,inda,sa,ra,ca,cv,rv,sv,indv,n,segid,nhits,nhitsmost;
  int nmisses, nfilled;
  MATRIX *Va2v, *Ka, *Kv, *invKv, *Vv2a, *Pa, *Pv;
  ASEGVOLINDEX *avind;
  MRI *volaseg, *volhit;

  segid=0;

  // Compute matrix to map from aseg CRS to vol CRS
  // Va2v = inv(Kv) * tkR * Ka
  Ka = MRIxfmCRS2XYZtkreg(aseg);
  Kv = MRIxfmCRS2XYZtkreg(voltemp);
  invKv = MatrixInverse(Kv,NULL);
  Va2v = MatrixMultiply(invKv,tkR,NULL);
  Va2v = MatrixMultiply(Va2v,Ka,Va2v);
  Vv2a = MatrixInverse(Va2v,NULL);
  MatrixFree(&Ka);
  MatrixFree(&Kv);
  MatrixFree(&invKv);

  Na = aseg->width * aseg->height * aseg->depth;
  avind = (ASEGVOLINDEX *) calloc(Na,sizeof(ASEGVOLINDEX));

  // Build an LUT that maps from aseg voxel to the closest
  // output volume voxel (reverse is done below).
  printf("ASeg2Vol: Building LUT\n");
  fflush(stdout);
  Pa = MatrixConstVal(0,4,1,NULL);
  Pa->rptr[4][1] = 1;
  Pv = MatrixConstVal(0,4,1,NULL);
  inda = 0;
  for (sa=0; sa < aseg->depth; sa++)
  {
    for (ra=0; ra < aseg->height; ra++)
    {
      for (ca=0; ca < aseg->width; ca++)
      {
        avind[inda].asegindex = inda;
        avind[inda].segid = MRIgetVoxVal(aseg,ca,ra,sa,0);
        // map each aseg vox into the volume
        Pa->rptr[1][1] = ca;
        Pa->rptr[2][1] = ra;
        Pa->rptr[3][1] = sa;
        MatrixMultiply(Va2v,Pa,Pv);
        cv = (int)round(Pv->rptr[1][1]);
        rv = (int)round(Pv->rptr[2][1]);
        sv = (int)round(Pv->rptr[3][1]);

        if (cv < 0 || cv >= voltemp->width ||
            rv < 0 || rv >= voltemp->height ||
            sv < 0 || sv >= voltemp->depth)
        {
          // check for out-of-FOV
          avind[inda].volindex = -1;
          avind[inda].cv = -1;
          avind[inda].rv = -1;
          avind[inda].sv = -1;
        }
        else
        {
          // save the volume crs and index
          indv = cv + (rv * voltemp->width) + (sv * voltemp->width*voltemp->height);
          avind[inda].cv = cv;
          avind[inda].rv = rv;
          avind[inda].sv = sv;
          avind[inda].volindex = indv;
        }
        inda++;
      } // col
    } // row
  } // slice

  printf("ASeg2Vol: Sorting \n");
  fflush(stdout);
  qsort(avind,Na,sizeof(ASEGVOLINDEX),CompareAVIndices);

  // Alloc output volume
  volaseg = MRIalloc(voltemp->width,voltemp->height,voltemp->depth,MRI_INT);
  MRIcopyHeader(voltemp,volaseg);

  // Alloc a hit volume (all values should be 0)
  volhit = MRIalloc(voltemp->width,voltemp->height,voltemp->depth,MRI_INT);
  MRIcopyHeader(voltemp,volhit);
  *pvolhit = volhit;

  // Go through each volume voxel and determine which seg id has the
  // most representation.
  printf("ASeg2Vol: Mapping\n");
  fflush(stdout);
  n = 0;
  while (n < Na)
  {
    indv = avind[n].volindex;
    if (indv < 0)
    {
      // aseg vox not in fov of outvol
      n++;
      continue;
    }

    // Count total number of aseg voxels falling into this volume vox
    nhits = 0;
    while ((n+nhits < Na) && (indv == avind[n+nhits].volindex)) nhits++;

    // Determine which segid had the most hits
    nhitsmost = MostHitsInVolVox(&avind[n], nhits, &segid);
    MRIsetVoxVal(volhit,avind[n].cv,avind[n].rv,avind[n].sv,0,nhitsmost);

    // Threshold (changed method on Jan 4, 2011)
    // New way based on number of segid hits relative to number of voxel hits
    if((double)nhitsmost/nhits < fthresh) segid=0;
    // Old way based on voxel volume
    //if (nhitsmost < nhitsthresh) segid=0; 

    // Set the output
    MRIsetVoxVal(volaseg,avind[n].cv,avind[n].rv,avind[n].sv,0,segid);

    if (Gdiag_no > 1)
    {
      printf("%6d %6d %5d    %6d    %3d %3d %3d     %3d %3d\n",
             n, avind[n].asegindex, segid, avind[n].volindex,
             avind[n].cv, avind[n].rv, avind[n].sv,  nhits, nhitsmost);
      fflush(stdout);
    }

    n += nhits;
  }

  // Now do the reverse map. It is possible that some output voxels
  // will not be closest to any aseg voxel. So find the output voxels
  // that were not initially mapped and assign them the seg id of the
  // closest in the seg volume.
  printf("ASeg2Vol: Reverse Map\n");
  nmisses = 0;
  nfilled = 0;
  for (sv=0; sv < volaseg->depth; sv++)
  {
    for (rv=0; rv < volaseg->height; rv++)
    {
      for (cv=0; cv < volaseg->width; cv++)
      {
        nhits = MRIgetVoxVal(volhit,cv,rv,sv,0);
        if (nhits != 0) continue;
        nmisses++;
        Pv->rptr[1][1] = cv;
        Pv->rptr[2][1] = rv;
        Pv->rptr[3][1] = sv;
        MatrixMultiply(Vv2a,Pv,Pa);
        ca = (int)round(Pa->rptr[1][1]);
        ra = (int)round(Pa->rptr[2][1]);
        sa = (int)round(Pa->rptr[3][1]);
        if (ca < 0 || ca >= aseg->width ||
            ra < 0 || ra >= aseg->height ||
            sa < 0 || sa >= aseg->depth)
        {
          // out-of-bounds
          MRIsetVoxVal(volaseg,cv,rv,sv,0,0);
        }
        else
        {
          // in-of-bounds
          segid = MRIgetVoxVal(aseg,ca,ra,sa,0);
          MRIsetVoxVal(volaseg,cv,rv,sv,0,segid);
          MRIsetVoxVal(volhit,cv,rv,sv,0,1);
          nfilled++;
        }
      }
    }
  }
  printf("nmisses = %d (%d filled)\n",nmisses,nfilled);


  //MRIwrite(volaseg,"volaseg.mgh");
  //MRIwrite(volhit,"volhits.mgh");
  MatrixFree(&Va2v);
  MatrixFree(&Vv2a);
  free(avind);

  printf("ASeg2Vol: done\n");
  return(volaseg);
}
/*-----------------------------------------------------------------------
  MostHitsInVolVox() - determines the segid with the most hits in the
  given volume voxel. In this case, avindsorted points to the first
  entry of the given voxel, and there are N entries for that voxel.
  -----------------------------------------------------------------------*/
static int MostHitsInVolVox(ASEGVOLINDEX *avindsorted, int N, int *segidmost)
{
  int n,nhits,nmost=0,segid;

  n=0;
  while (n < N)
  {
    segid = avindsorted[n].segid;
    nhits = 0;
    while ( (n+nhits < N) && (segid == avindsorted[n+nhits].segid)) nhits++;
    if (n==0)
    {
      *segidmost = segid;
      nmost = nhits;
    }
    if (nmost < nhits)
    {
      *segidmost = segid;
      nmost = nhits;
    }
    n += nhits;
  }
  return(nmost);
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
  ASEGVOLINDEX *avind1 = NULL, *avind2=NULL;

  avind1 = (ASEGVOLINDEX *) i1;
  avind2 = (ASEGVOLINDEX *) i2;

  // Sort by volume index
  if (avind1->volindex > avind2->volindex) return(+1);
  if (avind1->volindex < avind2->volindex) return(-1);

  // Only gets here if vol indices are the same, so
  // sort by seg id
  if (avind1->segid > avind2->segid) return(+1);
  if (avind1->segid < avind2->segid) return(-1);

  // Same volume index and seg id.
  return(0);
}
