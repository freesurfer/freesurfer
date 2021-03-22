/**
 * @brief resample volume to isotropic 1mm^3
 *
 */
/*
 * Original Author: Christian Haselgrove
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

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "error.h"
#include "histo.h"
#include "mri.h"

#include "mri_conform.h"

/*-------------------------------------------------------------------*/
double round(double);  // why is this never defined?!?
/*-------------------------------------------------------------------*/

extern int errno;

MATRIX *MRIgetConformMatrix(MRI *mri)
{
  MRI *templ;
  MATRIX *m_resample;

  if (mri->ras_good_flag == 0) {
    setDirectionCosine(mri, MRI_CORONAL);
  }

  templ = MRIallocHeader(256, 256, 256, MRI_UCHAR, mri->nframes);

  templ->imnr0 = 1;
  templ->imnr1 = 256;
  templ->thick = 1.0;
  templ->ps = 1.0;
  templ->xsize = templ->ysize = templ->zsize = 1.0;
  templ->xstart = templ->ystart = templ->zstart = -128.0;
  templ->xend = templ->yend = templ->zend = 128.0;
  setDirectionCosine(templ, MRI_CORONAL);  // sets c_(r,a,s) = 0
  // retain the src c_(r,a,s)
  templ->c_r = mri->c_r;
  templ->c_a = mri->c_a;
  templ->c_s = mri->c_s;
  templ->ras_good_flag = 1;  // use c_(r,a,s)
  templ->tr = mri->tr;
  templ->te = mri->te;
  templ->flip_angle = mri->flip_angle;
  templ->ti = mri->ti;

  m_resample = MRIgetResampleMatrix(mri, templ);

  MRIfree(&templ);

  return (m_resample);
}

MRI *MRIconform(MRI *mri)
{
  MRI *templ, *mri2, *res;
  // int conform_width;
  // double conform_size;
  // int KeepDC;

  res = MRIcopy(mri, NULL); /* don't mess with the input */

  if (res->ras_good_flag == 0) setDirectionCosine(res, MRI_CORONAL);

// conform_width = 256;
// conform_size = 1;
// KeepDC = 0;

#if 0
  templ = MRIconformedTemplate(mri, conform_width, conform_size, KeepDC);
#else
  templ = MRIallocHeader(256, 256, 256, MRI_UCHAR, mri->nframes);
  templ->imnr0 = 1;
  templ->imnr1 = 256;
  templ->thick = 1.0;
  templ->ps = 1.0;
  templ->xsize = templ->ysize = templ->zsize = 1.0;
  templ->xstart = templ->ystart = templ->zstart = -128.0;
  templ->xend = templ->yend = templ->zend = 128.0;
  setDirectionCosine(templ, MRI_CORONAL);  // sets c_(r,a,s) = 0
  // retain the c_(r,a,s)
  templ->c_r = res->c_r;
  templ->c_a = res->c_a;
  templ->c_s = res->c_s;
  templ->ras_good_flag = 1;  // use c_(r,a,s)
  templ->tr = mri->tr;
  templ->te = mri->te;
  templ->flip_angle = mri->flip_angle;
  templ->ti = mri->ti;
#endif

  /* ----- change type if necessary ----- */
  if (res->type != templ->type) {
    // FALSE means don't rescale
    mri2 = MRIchangeType(res, templ->type, 0.0, 0.999, FALSE);
    MRIfree(&res);
    if (mri2 == NULL) return (NULL);
    res = mri2;
  }

  /* ----- reslice if necessary ----- */
  if (res->xsize != templ->xsize || res->ysize != templ->ysize || res->zsize != templ->zsize ||
      res->width != templ->width || res->height != templ->height || res->depth != templ->depth ||
      res->x_r != templ->x_r || res->x_a != templ->x_a || res->x_s != templ->x_s || res->y_r != templ->y_r ||
      res->y_a != templ->y_a || res->y_s != templ->y_s || res->z_r != templ->z_r || res->z_a != templ->z_a ||
      res->z_s != templ->z_s) {
    // Need to be able to spec CUBIC here
    mri2 = MRIresample(res, templ, SAMPLE_TRILINEAR);
    MRIfree(&res);
    if (mri2 == NULL) return (NULL);
    res = mri2;
  }

  MRIfree(&templ);
  return (res);

} /*  end MRIconform()  */

MRI *MRIconformedTemplate(MRI *mri, int conform_width, double conform_size, int KeepDC)
{
  MRI *mri_template;
  char ostr[4];
  int iLR, iIS, iAP, Nvox[3], FoV[3], conform_FoV, c;
  double delta[3], pad, step;
  MATRIX *K, *invK, *Smri, *Stemp;

  mri_template = MRIallocHeader(conform_width, conform_width, conform_width, MRI_UCHAR, mri->nframes);
  MRIcopyHeader(mri, mri_template);
  MRIcopyPulseParameters(mri, mri_template);
  mri_template->imnr0 = 1;
  mri_template->imnr1 = conform_width;
  mri_template->thick = conform_size;
  mri_template->ps = conform_size;
  mri_template->xsize = mri_template->ysize = mri_template->zsize = conform_size;
  mri_template->xstart = mri_template->ystart = mri_template->zstart = -conform_width / 2;
  mri_template->xend = mri_template->yend = mri_template->zend = conform_width / 2;

  if (KeepDC) {
    conform_FoV = conform_width * conform_size;
    MRIdircosToOrientationString(mri, ostr);
    for (iLR = 0; iLR < 3; iLR++)
      if (ostr[iLR] == 'L' || ostr[iLR] == 'R') break;
    for (iIS = 0; iIS < 3; iIS++)
      if (ostr[iIS] == 'I' || ostr[iIS] == 'S') break;
    for (iAP = 0; iAP < 3; iAP++)
      if (ostr[iAP] == 'A' || ostr[iAP] == 'P') break;
    printf("keeping DC %d %d %d\n", iLR, iIS, iAP);
    printf("ostr %s, width %d, size %g\n", ostr, conform_width, conform_size);

    Nvox[0] = mri->width;
    Nvox[1] = mri->height;
    Nvox[2] = mri->depth;
    delta[0] = mri->xsize;
    delta[1] = mri->ysize;
    delta[2] = mri->zsize;
    for (c = 0; c < 3; c++) FoV[c] = Nvox[c] * delta[c];

    // K maps voxels in mri to voxels in mri_template
    K = MatrixAlloc(4, 4, MATRIX_REAL);
    K->rptr[4][4] = 1;

    // If the delta=conform_size, then no interpolation will result
    // Otherwise, there will be interpolation that depends on voxel size
    // Using round() forces no interpolation at the edge of the FoV
    // pad is the number of conformed voxels of padding when Nvox != conform_width
    // set pad this way makes the C_RASs be about the same under general conditions

    step = delta[iLR] / conform_size;
    pad = round(((conform_FoV - FoV[iLR]) / 2.0) / conform_size);
    if (ostr[iLR] == 'L') {
      K->rptr[1][iLR + 1] = step;
      K->rptr[1][4] = pad;
    }
    else {
      K->rptr[1][iLR + 1] = -step;
      K->rptr[1][4] = conform_width - pad;
    }

    step = delta[iIS] / conform_size;
    pad = round(((conform_FoV - FoV[iIS]) / 2.0) / conform_size);
    if (ostr[iIS] == 'I') {
      K->rptr[2][iIS + 1] = step;
      K->rptr[2][4] = pad;
    }
    else {
      K->rptr[2][iIS + 1] = -step;
      K->rptr[2][4] = conform_width - pad;
    }

    step = delta[iAP] / conform_size;
    pad = round(((conform_FoV - FoV[iAP]) / 2.0) / conform_size);
    if (ostr[iAP] == 'A') {
      K->rptr[3][iAP + 1] = step;
      K->rptr[3][4] = pad;
    }
    else {
      K->rptr[3][iAP + 1] = -step;
      K->rptr[3][4] = conform_width - pad;
    }

    invK = MatrixInverse(K, NULL);
    Smri = MRIxfmCRS2XYZ(mri, 0);
    Stemp = MatrixMultiplyD(Smri, invK, NULL);
    MRIsetVox2RASFromMatrix(mri_template, Stemp);

    printf("K ---------------\n");
    MatrixPrint(stdout, K);
    printf("Kinv ---------------\n");
    MatrixPrint(stdout, invK);
    printf("Smri ---------------\n");
    MatrixPrint(stdout, Smri);
    printf("Stemp ---------------\n");
    MatrixPrint(stdout, Stemp);
    printf("----------------------\n");

    MatrixFree(&K);
    MatrixFree(&invK);
    MatrixFree(&Smri);
    MatrixFree(&Stemp);
  }
  else {
    // replicates old method exactly
    mri_template->x_r = -1.0;
    mri_template->x_a = 0.0;
    mri_template->x_s = 0.0;
    mri_template->y_r = 0.0;
    mri_template->y_a = 0.0;
    mri_template->y_s = -1.0;
    mri_template->z_r = 0.0;
    mri_template->z_a = 1.0;
    mri_template->z_s = 0.0;
  }

  return (mri_template);
}
