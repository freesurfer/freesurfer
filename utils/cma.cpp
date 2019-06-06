/**
 * @file  cma.c
 * @brief constants for neuroanatomical structures.
 *
 * constants and macros for neuroanatomical and some vascular structures.
 * Originally it was just cma labels, but it has been generalized.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2015/10/05 23:59:03 $
 *    $Revision: 1.28 $
 *
 * Copyright © 2011-2014 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "error.h"
#include "fio.h"
#include "gtm.h"
#include "mrisutils.h"

/* see ch notebook 2 */

#include "cma.h"

#ifndef TRUE
#define TRUE (1 == 1)
#endif

#ifndef FALSE
#define FALSE (0 == 1)
#endif

extern int errno;

int CMAfreeOutlineField(CMAoutlineField **of)
{
  CMAoutlineField *ofp;
  int i;

  ofp = *of;

  if (ofp->claim_field != NULL) {
    for (i = 0; i < ofp->height; i++) {
      if (ofp->claim_field[i] != NULL) free(ofp->claim_field[i]);
    }
    free(ofp->claim_field);
  }

  if (ofp->fill_field != NULL) {
    for (i = 0; i < ofp->height; i++) {
      if (ofp->fill_field[i] != NULL) free(ofp->fill_field[i]);
    }
    free(ofp->fill_field);
  }

  if (ofp->outline_points_field != NULL) {
    for (i = 0; i < ofp->height; i++) {
      if (ofp->outline_points_field[i] != NULL) free(ofp->outline_points_field[i]);
    }
    free(ofp->outline_points_field);
  }

  *of = NULL;

  return (NO_ERROR);

} /* end CMAfreeOutlineField() */

CMAoutlineField *CMAoutlineFieldAlloc(int width, int height)
{
  CMAoutlineField *of;
  int i;

  of = (CMAoutlineField *)malloc(sizeof(CMAoutlineField));
  if (of == NULL) ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating structure"));

  of->claim_field = NULL;
  of->fill_field = NULL;
  of->outline_points_field = NULL;
  of->width = width;
  of->height = height;

  of->claim_field = (CMAoutlineClaim **)malloc(height * sizeof(CMAoutlineClaim *));
  if (of->claim_field == NULL) {
    CMAfreeOutlineField(&of);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating claim field"));
  }
  memset(of->claim_field, 0x00, height * sizeof(CMAoutlineClaim *));

  // of->fill_field = (unsigned char **)malloc(height * sizeof(unsigned char *));
  of->fill_field = (short **)malloc(height * sizeof(short *));
  if (of->fill_field == NULL) {
    CMAfreeOutlineField(&of);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating fill field"));
  }
  // memset(of->fill_field, 0x00, height * sizeof(unsigned char *));
  memset(of->fill_field, 0x00, height * sizeof(short *));

  of->outline_points_field = (unsigned char **)malloc(height * sizeof(unsigned char *));
  if (of->outline_points_field == NULL) {
    CMAfreeOutlineField(&of);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating outline points field"));
  }
  memset(of->outline_points_field, 0x00, height * sizeof(unsigned char *));

  for (i = 0; i < height; i++) {
    of->claim_field[i] = (CMAoutlineClaim *)malloc(width * sizeof(CMAoutlineClaim));
    if (of->claim_field[i] == NULL) {
      CMAfreeOutlineField(&of);
      ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating claim field"));
    }
    memset(of->claim_field[i], 0x00, width * sizeof(CMAoutlineClaim));

    // of->fill_field[i] = (unsigned char *)malloc(width * sizeof(unsigned char));
    of->fill_field[i] = (short *)malloc(width * sizeof(short));
    if (of->fill_field[i] == NULL) {
      CMAfreeOutlineField(&of);
      ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating fill field"));
    }
    // memset(of->fill_field[i], 0x00, width * sizeof(unsigned char));
    memset(of->fill_field[i], 0x00, width * sizeof(short));

    of->outline_points_field[i] = (unsigned char *)malloc(width * sizeof(unsigned char));
    if (of->outline_points_field[i] == NULL) {
      CMAfreeOutlineField(&of);
      ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating outline points field"));
    }
    memset(of->outline_points_field[i], 0x00, width * sizeof(unsigned char));
  }

  return (of);

} /* end CMAoutlineFieldAlloc() */

int CMAclearFillField(CMAoutlineField *field)
{
  int i;

  for (i = 0; i < field->height; i++)
    // memset(field->fill_field[i], 0x00, field->width * sizeof(unsigned char));
    memset(field->fill_field[i], 0x00, field->width * sizeof(short));

  return (NO_ERROR);

} /* end CMAclearFillField() */

/* fills with CMA_FILL_INTERIOR to non-zero values */
int CMAfill(CMAoutlineField *field, short seed_x, short seed_y)
{
  if (seed_x < 0 || seed_x >= field->width) return (NO_ERROR);

  if (seed_y < 0 || seed_y >= field->height) return (NO_ERROR);

  if (field->fill_field[seed_y][seed_x] != 0) return (NO_ERROR);

  field->fill_field[seed_y][seed_x] = CMA_FILL_INTERIOR;

  CMAfill(field, seed_x + 1, seed_y);
  CMAfill(field, seed_x - 1, seed_y);
  CMAfill(field, seed_x, seed_y + 1);
  CMAfill(field, seed_x, seed_y - 1);

  return (NO_ERROR);

} /* end CMAfill() */

int CMAclaimPoints(CMAoutlineField *field, short label, short *points, int n_points, short seed_x, short seed_y)
{
  int i, j;
  short x, y;

  if (label < 0 || label > MAX_CMA_LABEL)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "CMAclaimPoints(): label out of range (label is %d, MAX_CMA_LABEL is %d)",
                 label,
                 MAX_CMA_LABEL));

  if (seed_x < 0 || seed_x >= field->width)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "CMAclaimPoints(): seed point out of range (seed_x = %d, field width = %d)",
                 seed_x,
                 field->width));
  if (seed_y < 0 || seed_y >= field->height)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "CMAclaimPoints(): seed point out of range (seed_y = %d, field height = %d)",
                 seed_y,
                 field->height));

  CMAclearFillField(field);

  for (i = 0; i < n_points; i++) {
    x = points[2 * i];
    y = points[2 * i + 1];
    if (x < 0 || x >= field->width)
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): outline point out of range (x)"));
    if (y < 0 || y >= field->height)
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): outline point out of range (y)"));
    field->fill_field[y][x] = CMA_FILL_OUTLINE;
    field->outline_points_field[y][x] = 1;
  }

  CMAfill(field, seed_x, seed_y);

  for (i = 0; i < field->width; i++) {
    for (j = 0; j < field->height; j++) {
      if (field->fill_field[j][i] == CMA_FILL_INTERIOR) {
        field->claim_field[j][i].n_claims = 1;
        field->claim_field[j][i].interior_claim_flag = TRUE;
        field->claim_field[j][i].claim_labels[0] = label;
      }

      if (field->fill_field[j][i] == CMA_FILL_OUTLINE) {
        if (field->claim_field[j][i].n_claims < MAX_OUTLINE_CLAIMS) {
          field->claim_field[j][i].claim_labels[field->claim_field[j][i].n_claims] = label;
          field->claim_field[j][i].n_claims++;
          field->claim_field[j][i].interior_claim_flag = FALSE;
        }
      }
    }
  }

  return (NO_ERROR);

} /* end CMAclaimPoints() */

int CMAvalueClaims(CMAoutlineClaim *claim)
{
  int i;

  for (i = 0; i < MAX_OUTLINE_CLAIMS; i++) claim->claim_values[i] = 0.0;

  if (claim->n_claims == 0) {
    claim->no_label_claim = 0.5;
  }
  else if (claim->n_claims == 1) {
    if (claim->interior_claim_flag == 1)
      claim->claim_values[0] = 1.0;
    else {
      claim->claim_values[0] = 0.5;
      claim->no_label_claim = 0.5;
    }
  }
  else {
    float ct = 1.0 / (float)claim->n_claims;
    for (i = 0; i < claim->n_claims; i++) claim->claim_values[i] = ct;
  }

  return (NO_ERROR);

} /* end CMAvalueClaims() */

int CMAvalueAllClaims(CMAoutlineField *field)
{
  int i, j;

  for (i = 0; i < field->width; i++)
    for (j = 0; j < field->height; j++) CMAvalueClaims(&(field->claim_field[j][i]));

  return (NO_ERROR);

} /* end CMAvalueAllClaims() */

int CMAaddWeightedTotals(CMAoutlineClaim *claim, float weight, float *claim_totals)
{
  int i;

  /* we've checked the label range in CMAclaimPoints() */

  for (i = 0; i < claim->n_claims; i++) claim_totals[claim->claim_labels[i]] += weight * claim->claim_values[i];

  claim_totals[0] += weight * claim->no_label_claim;

  return (NO_ERROR);

} /* end CMAaddWeightedTotals() */

/* returns the label with the greatest claim to a given pixel */
short CMAtotalClaims(CMAoutlineField *field, int x, int y)
{
  float claim_totals[MAX_CMA_LABEL + 1];
  float best_claim;
  short best_index;
  int i;

  if (x < 0 || x >= field->width)
    ErrorReturn(-1, (ERROR_BADPARM, "CMAtotalClaims(): x out of range (x = %d, field->width = %d)", x, field->width));

  if (y < 0 || y >= field->height)
    ErrorReturn(-1, (ERROR_BADPARM, "CMAtotalClaims(): y out of range (y = %d, field->height = %d)", y, field->height));

  for (i = 0; i < MAX_CMA_LABEL; i++) claim_totals[i] = 0.0;

  /* center point (x, y checks done above) */
  CMAaddWeightedTotals(&(field->claim_field[y][x]), OTL_CLAIM_WEIGHT_CENTER, claim_totals);

  /* immediately adjoining points */
  if (x - 1 >= 0) CMAaddWeightedTotals(&(field->claim_field[y][x - 1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x + 1 < field->width) CMAaddWeightedTotals(&(field->claim_field[y][x + 1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (y - 1 >= 0) CMAaddWeightedTotals(&(field->claim_field[y - 1][x]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (y + 1 < field->height) CMAaddWeightedTotals(&(field->claim_field[y + 1][x]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);

  /* diagonally adjoining points */
  if (x - 1 >= 0 && y - 1 >= 0)
    CMAaddWeightedTotals(&(field->claim_field[y - 1][x - 1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x - 1 >= 0 && y + 1 < field->height)
    CMAaddWeightedTotals(&(field->claim_field[y + 1][x - 1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x + 1 < field->width && y - 1 >= 0)
    CMAaddWeightedTotals(&(field->claim_field[y - 1][x + 1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x + 1 < field->width && y + 1 < field->height)
    CMAaddWeightedTotals(&(field->claim_field[y + 1][x + 1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);

  /* find the highest claim and its index */
  best_claim = claim_totals[0];
  best_index = 0;
  for (i = 1; i <= MAX_CMA_LABEL; i++) {
    if (claim_totals[i] > best_claim) {
      best_claim = claim_totals[i];
      best_index = i;
    }
  }

  return (best_index);

} /* end CMAtotalClaims() */

int CMAassignLabels(CMAoutlineField *field)
{
  int i, j;

  CMAclearFillField(field);

  CMAvalueAllClaims(field);

  for (i = 0; i < field->width; i++)
    for (j = 0; j < field->height; j++) field->fill_field[j][i] = CMAtotalClaims(field, i, j);

  return (NO_ERROR);

} /* end CMAassignLabels() */

int CMAzeroOutlines(CMAoutlineField *field)
{
  int i, j;

  for (i = 0; i < field->width; i++)
    for (j = 0; j < field->height; j++)
      if (field->outline_points_field[j][i]) field->fill_field[j][i] = 0;

  return (NO_ERROR);

} /* end CMAzeroOutlines() */

#include "diag.h"
#include "error.h"
#include "macros.h"
#include "mrisurf.h"
int insert_ribbon_into_aseg(MRI *mri_src_aseg, MRI *mri_aseg, MRI_SURFACE *mris_white, MRI_SURFACE *mris_pial, int hemi)
{
  MRI *mri_ribbon, *mri_white;
  int x, y, z, gm_label, wm_label, label, nbr_label, dont_change;

  if (mri_src_aseg != mri_aseg) mri_aseg = MRIcopy(mri_src_aseg, mri_aseg);

  gm_label = hemi == LEFT_HEMISPHERE ? Left_Cerebral_Cortex : Right_Cerebral_Cortex;
  wm_label = hemi == LEFT_HEMISPHERE ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;

  mri_white = MRIclone(mri_aseg, NULL);
  mri_ribbon = MRISribbon(mris_white, mris_pial, mri_aseg, NULL);
  MRISfillInterior(mris_white, mri_aseg->xsize, mri_white);

  for (x = 0; x < mri_aseg->width; x++)
    for (y = 0; y < mri_aseg->height; y++)
      for (z = 0; z < mri_aseg->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        label = nint(MRIgetVoxVal(mri_aseg, x, y, z, 0));
        if (nint(MRIgetVoxVal(mri_ribbon, x, y, z, 0)) == 255)  // in ribbon, set to GM
        {
          if (IS_CORTEX(label) || IS_WHITE_CLASS(label)) {
            int xi, yi, zi, xk, yk, zk;
            // check to make sure we are really in cortex
            for (dont_change = 0, xk = -1; xk <= 1; xk++) {
              xi = mri_aseg->xi[xk + x];
              for (yk = -1; yk <= 1; yk++) {
                yi = mri_aseg->yi[yk + y];
                for (zk = -1; zk <= 1; zk++) {
                  zi = mri_aseg->zi[zk + z];
                  nbr_label = (int)MRIgetVoxVal(mri_aseg, xi, yi, zi, 0);
                  switch (nbr_label) {
                    default:
                      break;
                    case Left_Hippocampus:
                    case Right_Hippocampus:
                    case Left_Amygdala:
                    case Right_Amygdala:
                    case Left_VentralDC:
                    case Right_VentralDC:
                    case Brain_Stem:
                    case Left_Lateral_Ventricle:
                    case Right_Lateral_Ventricle:
                    case Left_Inf_Lat_Vent:
                    case Right_Inf_Lat_Vent:
                    case Left_Thalamus:
                    case Right_Thalamus:
                    case Left_choroid_plexus:
                    case Right_choroid_plexus:
                    case CC_Posterior:
                    case CC_Mid_Posterior:
                    case CC_Central:
                    case CC_Mid_Anterior:
                    case CC_Anterior:
                      dont_change = 1;
                      break;
                  }
                }
              }
            }
            if (dont_change == 0) MRIsetVoxVal(mri_aseg, x, y, z, 0, gm_label);
          }
        }
        else  // not in ribbon
        {
          if (MRIgetVoxVal(mri_white, x, y, z, 0) > 0)  // inside white surface - disambiguate
          {
            if (label == gm_label)  // gm inside white surface should be wm
              MRIsetVoxVal(mri_aseg, x, y, z, 0, wm_label);
          }
          else if (label == gm_label)  // gm outside ribbon should be unknown
            MRIsetVoxVal(mri_aseg, x, y, z, 0, Unknown);
        }
      }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRIwrite(mri_ribbon, "r.mgz");
    MRIwrite(mri_white, "w.mgz");
  }
  MRIfree(&mri_ribbon);
  MRIfree(&mri_white);
  return (NO_ERROR);
}

/*!
\fn int IsSubCorticalGray(int SegId)
\brief Returns a 1 if the given seg is subcortical gray. Note:
subcortical gray does not include brainstem or cerebellum cortex.
\param SegId - segmentation id number
*/
int IsSubCorticalGray(int SegId)
{
  if (SegId == Left_Thalamus) return (1);
  if (SegId == Right_Thalamus) return (1);
  if (SegId == Left_Thalamus) return (1);
  if (SegId == Right_Thalamus) return (1);
  if (SegId == Left_Caudate) return (1);
  if (SegId == Right_Caudate) return (1);
  if (SegId == Left_Putamen) return (1);
  if (SegId == Right_Putamen) return (1);
  if (SegId == Left_Pallidum) return (1);
  if (SegId == Right_Pallidum) return (1);
  // if(SegId == Brain_Stem) return(1);
  if (SegId == Left_Hippocampus) return (1);
  if (SegId == Right_Hippocampus) return (1);
  if (SegId == Left_Amygdala) return (1);
  if (SegId == Right_Amygdala) return (1);
  if (SegId == Left_Accumbens_area) return (1);
  if (SegId == Right_Accumbens_area) return (1);
  if (SegId == Left_VentralDC) return (1);
  if (SegId == Right_VentralDC) return (1);
  if (SegId == Left_Substancia_Nigra) return (1);
  if (SegId == Right_Substancia_Nigra) return (1);
  // if(SegId == Left_Cerebellum_Cortex) return(1);
  // if(SegId == Right_Cerebellum_Cortex) return(1);
  return (0);
}

/*!
\fn double SupraTentorialVolCorrection(MRI *aseg, MRI *ribbon)
\brief Returns the volume of supratentorial structures that do not
fall inside the pial surface or are cut by the pial surface.  The idea
is that the volume of everything in the pial surface can be computed
using the surface and that this function can be used to compute
everything else. Note that there is no partial volume correction.
Note 4/2/2019: there was a bug here in that the voxel volume was
defined as an int instead of double which means that the used voxel
volume would have been truncated and affecting the total.
\param aseg - aseg.mgz or aparc+aseg.mgz
\param ribbon is the ribbon.mgz, which has non-zero values for
everything inside the pial surf.
*/
double SupraTentorialVolCorrection(MRI *aseg, MRI *ribbon)
{
  int c, r, s, SegId;
  double vol = 0;
  int RibbonVal;
  double VoxSize; // was int

  VoxSize = aseg->xsize * aseg->ysize * aseg->zsize;
  for (c = 0; c < aseg->width; c++) {
    for (r = 0; r < aseg->height; r++) {
      for (s = 0; s < aseg->depth; s++) {
        // If this voxel is inside the pial, then skip it because it
        // will be part of the surface-based volume measure
        // This was the location of a bug found by Mike Harms.
        RibbonVal = MRIgetVoxVal(ribbon, c, r, s, 0);
        if (RibbonVal > 0) continue;

        // If it gets here, it means that the voxel was not within
        // the pial surface. It could be in a structure that should
        // be part of the supratentorium.

        // These are midline, medial wall, or unknown structures
        // that the pial could cut through.
        SegId = MRIgetVoxVal(aseg, c, r, s, 0);
        if (SegId == Left_Lateral_Ventricle) vol += VoxSize;
        if (SegId == Right_Lateral_Ventricle) vol += VoxSize;
        if (SegId == Left_choroid_plexus) vol += VoxSize;
        if (SegId == Right_choroid_plexus) vol += VoxSize;
        if (SegId == Left_Inf_Lat_Vent) vol += VoxSize;
        if (SegId == Right_Inf_Lat_Vent) vol += VoxSize;
        if (SegId == WM_hypointensities) vol += VoxSize;
        if (SegId == Left_WM_hypointensities) vol += VoxSize;
        if (SegId == Right_WM_hypointensities) vol += VoxSize;
        if (SegId == Left_Thalamus) vol += VoxSize;
        if (SegId == Right_Thalamus) vol += VoxSize;
        if (SegId == Left_Thalamus) vol += VoxSize;
        if (SegId == Right_Thalamus) vol += VoxSize;
        if (SegId == CC_Posterior) vol += VoxSize;
        if (SegId == CC_Mid_Posterior) vol += VoxSize;
        if (SegId == CC_Central) vol += VoxSize;
        if (SegId == CC_Mid_Anterior) vol += VoxSize;
        if (SegId == CC_Anterior) vol += VoxSize;
        if (SegId == Left_VentralDC) vol += VoxSize;
        if (SegId == Right_VentralDC) vol += VoxSize;
        if (SegId == Left_Hippocampus) vol += VoxSize;
        if (SegId == Right_Hippocampus) vol += VoxSize;
        if (SegId == Left_Amygdala) vol += VoxSize;
        if (SegId == Right_Amygdala) vol += VoxSize;

        // These are unlikely to have the pial surface cut through
        // them, but no harm to include them
        if (SegId == Left_Caudate) vol += VoxSize;
        if (SegId == Right_Caudate) vol += VoxSize;
        if (SegId == Left_Putamen) vol += VoxSize;
        if (SegId == Right_Putamen) vol += VoxSize;
        if (SegId == Left_Pallidum) vol += VoxSize;
        if (SegId == Right_Pallidum) vol += VoxSize;
        if (SegId == Left_Accumbens_area) vol += VoxSize;
        if (SegId == Right_Accumbens_area) vol += VoxSize;
      }
    }
  }
  return (vol);
}

/*!
\fn double CorticalGMVolCorrection(MRI *aseg, MRI *ribbon)
\brief Computes the volume of non-cortical structures within the
ribbon. This should be subtracted from the cortical GM volume computed
by subtracting the volume inside the white from the volume inside the
pial. Note 4/2/2019: there was a bug here in that the voxel volume was
defined as an int instead of double which means that the used voxel
volume would have been truncated and affecting the total.
\param aseg - segmentation
\param ribbon - ribbon
\param hemi - 1=left, 2=right
*/
double CorticalGMVolCorrection(MRI *aseg, MRI *ribbon, int hemi)
{
  int c, r, s, SegId;
  double vol = 0, vol2 = 0;
  int RibbonVal;
  double VoxSize;

  VoxSize = aseg->xsize * aseg->ysize * aseg->zsize;
  for (c = 0; c < aseg->width; c++) {
    for (r = 0; r < aseg->height; r++) {
      for (s = 0; s < aseg->depth; s++) {
        // If this voxel is not inside the ribbon, then skip it
        RibbonVal = MRIgetVoxVal(ribbon, c, r, s, 0);
        if (hemi == 1 && RibbonVal != 3) continue;
        if (hemi == 2 && RibbonVal != 42) continue;

        // If it gets here, it means that the voxel was within the
        // ribbon. It could be in a structure that is not part of the
        // cortical GM

        SegId = MRIgetVoxVal(aseg, c, r, s, 0);
        // This could be done in two ways. (1) If the aseg is not
        // cortex (SegID=3 or 42), then make the correction. or (2) If
        // the SegId is in a list of non-cortical structures, then
        // make the correction. #2 is better because it assumes that
        // the aseg is correct for non-cortex (but requires listing
        // all the structures that could be affected). #1 is easier
        // but assumes that the aseg cortex label always correctly
        // declares cortex to be cortex.

        // This uses method 1 (for testing) - gives very similar value as method 2
        if (SegId != 3 && SegId != 42 && SegId != 2 && SegId != 41 && SegId != 0) vol2 += VoxSize;

        // Method 2
        // By far, most of the volume is going to be Hip and Amyg
        if (SegId == Left_Hippocampus) vol += VoxSize;
        if (SegId == Right_Hippocampus) vol += VoxSize;
        if (SegId == Left_Amygdala) vol += VoxSize;
        if (SegId == Right_Amygdala) vol += VoxSize;

        if (SegId == Left_Lateral_Ventricle) vol += VoxSize;
        if (SegId == Right_Lateral_Ventricle) vol += VoxSize;
        if (SegId == Left_choroid_plexus) vol += VoxSize;
        if (SegId == Right_choroid_plexus) vol += VoxSize;
        if (SegId == Left_Inf_Lat_Vent) vol += VoxSize;
        if (SegId == Right_Inf_Lat_Vent) vol += VoxSize;
        if (SegId == WM_hypointensities) vol += VoxSize;
        if (SegId == Left_WM_hypointensities) vol += VoxSize;
        if (SegId == Right_WM_hypointensities) vol += VoxSize;
        if (SegId == Left_Thalamus) vol += VoxSize;
        if (SegId == Right_Thalamus) vol += VoxSize;
        if (SegId == Left_Thalamus) vol += VoxSize;
        if (SegId == Right_Thalamus) vol += VoxSize;
        if (SegId == CC_Posterior) vol += VoxSize;
        if (SegId == CC_Mid_Posterior) vol += VoxSize;
        if (SegId == CC_Central) vol += VoxSize;
        if (SegId == CC_Mid_Anterior) vol += VoxSize;
        if (SegId == CC_Anterior) vol += VoxSize;
        if (SegId == Left_VentralDC) vol += VoxSize;
        if (SegId == Right_VentralDC) vol += VoxSize;
        if (SegId == Left_Caudate) vol += VoxSize;
        if (SegId == Right_Caudate) vol += VoxSize;
        if (SegId == Left_Putamen) vol += VoxSize;
        if (SegId == Right_Putamen) vol += VoxSize;
        if (SegId == Left_Pallidum) vol += VoxSize;
        if (SegId == Right_Pallidum) vol += VoxSize;
        if (SegId == Left_Accumbens_area) vol += VoxSize;
        if (SegId == Right_Accumbens_area) vol += VoxSize;
        if (SegId == Left_Cerebellum_Cortex) vol += VoxSize;
        if (SegId == Right_Cerebellum_Cortex) vol += VoxSize;
        if (SegId == Left_vessel) vol += VoxSize;
      }
    }
  }
  printf("CorticalGMVolCorrection(): hemi=%d, vol1=%g, vol2=%g\n", hemi, vol, vol2);
  return (vol);
}

/*!
  \fn int MRIasegContraLatLabel(int id)
  \brief returns the seg id of the contralateral segmentation of the
   passed seg id.  Should handle, aseg, aparc+aseg, aparc.a2009s+aseg,
   and wmparc, and samseg. If these segmentations change, then this
   program (written on 12/19/2011) may fail. If the seg id is unlateralized,
   then just returns id. If it is unrecognized, then it prints a message
   and returns -1.
 */
int MRIasegContraLatLabel(int id)
{ 
  int id2 = -1;

  if(id==0) return(0);

  if (id > 999) {
    // This should handle the cortical and wmparc labels in aparc+aseg, etc
    if (id >= 1000 && id < 1999) id2 = id + 1000;
    if (id >= 2000 && id < 2999) id2 = id - 1000;
    if (id >= 3000 && id < 3999) id2 = id + 1000;
    if (id >= 4000 && id < 4999) id2 = id - 1000;
    if (id >= 11100 && id < 11999) id2 = id + 1000;
    if (id >= 12100 && id < 12999) id2 = id - 1000;
    if (id >= 13100 && id < 13999) id2 = id + 1000;
    if (id >= 14100 && id < 14999) id2 = id - 1000;
    if (id == 5001) id2 = 5002;
    if (id == 5002) id2 = 5001;
    return(id2);
  }
  // This should handle all the subcortical stuff
  switch (id) {
  case Left_Cerebral_White_Matter:
    id2 = Right_Cerebral_White_Matter;
    break;
  case Left_Cerebral_Exterior:
    id2 = Right_Cerebral_Exterior;
    break;
  case Left_Cerebral_Cortex:
    id2 = Right_Cerebral_Cortex;
    break;
  case Left_Lateral_Ventricle:
    id2 = Right_Lateral_Ventricle;
    break;
  case Left_Inf_Lat_Vent:
    id2 = Right_Inf_Lat_Vent;
    break;
  case Left_Cerebellum_Exterior:
    id2 = Right_Cerebellum_Exterior;
    break;
  case Left_Cerebellum_White_Matter:
    id2 = Right_Cerebellum_White_Matter;
    break;
  case Left_Cerebellum_Cortex:
    id2 = Right_Cerebellum_Cortex;
    break;
  case Left_Thalamus:
    id2 = Right_Thalamus;
    break;
  case Left_Caudate:
    id2 = Right_Caudate;
    break;
  case Left_Putamen:
    id2 = Right_Putamen;
    break;
  case Left_Pallidum:
    id2 = Right_Pallidum;
    break;
  case Left_Hippocampus:
    id2 = Right_Hippocampus;
    break;
  case Left_Amygdala:
    id2 = Right_Amygdala;
    break;
  case Left_Insula:
    id2 = Right_Insula;
    break;
  case Left_Operculum:
    id2 = Right_Operculum;
    break;
  case Left_Lesion:
    id2 = Right_Lesion;
    break;
  case Left_Accumbens_area:
    id2 = Right_Accumbens_area;
    break;
  case Left_Substancia_Nigra:
    id2 = Right_Substancia_Nigra;
    break;
  case Left_VentralDC:
    id2 = Right_VentralDC;
    break;
  case Left_undetermined:
    id2 = Right_undetermined;
    break;
  case Left_vessel:
    id2 = Right_vessel;
    break;
  case Left_choroid_plexus:
    id2 = Right_choroid_plexus;
    break;
  case Left_F3orb:
    id2 = Right_F3orb;
    break;
  case Left_lOg:
    id2 = Right_lOg;
    break;
  case Left_aOg:
    id2 = Right_aOg;
    break;
  case Left_mOg:
    id2 = Right_mOg;
    break;
  case Left_pOg:
    id2 = Right_pOg;
    break;
  case Left_Stellate:
    id2 = Right_Stellate;
    break;
  case Left_Porg:
    id2 = Right_Porg;
    break;
  case Left_Aorg:
    id2 = Right_Aorg;
    break;
  case Left_Interior:
    id2 = Right_Interior;
    break;
  case Left_WM_hypointensities:
    id2 = Right_WM_hypointensities;
    break;
  case Left_non_WM_hypointensities:
    id2 = Right_non_WM_hypointensities;
    break;
  case Left_F1:
    id2 = Right_F1;
    break;
  case Left_Amygdala_Anterior:
    id2 = Right_Amygdala_Anterior;
    break;
  case Right_Cerebral_Exterior:
    id2 = Left_Cerebral_Exterior;
    break;
  case Right_Cerebral_White_Matter:
    id2 = Left_Cerebral_White_Matter;
    break;
  case Right_Cerebral_Cortex:
    id2 = Left_Cerebral_Cortex;
    break;
  case Right_Inf_Lat_Vent:
    id2 = Left_Inf_Lat_Vent;
    break;
  case Right_Cerebellum_Exterior:
    id2 = Left_Cerebellum_Exterior;
    break;
  case Right_Cerebellum_White_Matter:
    id2 = Left_Cerebellum_White_Matter;
    break;
  case Right_Cerebellum_Cortex:
    id2 = Left_Cerebellum_Cortex;
    break;
  case Right_Thalamus:
    id2 = Left_Thalamus;
    break;
  case Right_Caudate:
    id2 = Left_Caudate;
    break;
  case Right_Putamen:
    id2 = Left_Putamen;
    break;
  case Right_Pallidum:
    id2 = Left_Pallidum;
    break;
  case Right_Hippocampus:
    id2 = Left_Hippocampus;
    break;
  case Right_Amygdala:
    id2 = Left_Amygdala;
    break;
  case Right_Insula:
    id2 = Left_Insula;
    break;
  case Right_Operculum:
    id2 = Left_Operculum;
    break;
  case Right_Lesion:
    id2 = Left_Lesion;
    break;
  case Right_Accumbens_area:
    id2 = Left_Accumbens_area;
    break;
  case Right_Substancia_Nigra:
    id2 = Left_Substancia_Nigra;
    break;
  case Right_VentralDC:
    id2 = Left_VentralDC;
    break;
  case Right_undetermined:
    id2 = Left_undetermined;
    break;
  case Right_vessel:
    id2 = Left_vessel;
    break;
  case Right_choroid_plexus:
    id2 = Left_choroid_plexus;
    break;
  case Right_F3orb:
    id2 = Left_F3orb;
    break;
  case Right_lOg:
    id2 = Left_lOg;
    break;
  case Right_aOg:
    id2 = Left_aOg;
    break;
  case Right_mOg:
    id2 = Left_mOg;
    break;
  case Right_pOg:
    id2 = Left_pOg;
    break;
  case Right_Stellate:
    id2 = Left_Stellate;
    break;
  case Right_Porg:
    id2 = Left_Porg;
    break;
  case Right_Aorg:
    id2 = Left_Aorg;
    break;
  case Right_Interior:
    id2 = Left_Interior;
    break;
  case Right_Lateral_Ventricle:
    id2 = Left_Lateral_Ventricle;
    break;
  case Right_WM_hypointensities:
    id2 = Left_WM_hypointensities;
    break;
  case Right_non_WM_hypointensities:
    id2 = Left_non_WM_hypointensities;
    break;
  case Right_F1:
    id2 = Left_F1;
    break;
  case Right_Amygdala_Anterior:
    id2 = Left_Amygdala_Anterior;
    break;
  case 265: // Left-Eyeball
    id2 = 266; 
    break;
  case 266: // Right-Eyeball
    id2 = 265; 
    break;
  // These are unlateralized
  case 72: // 5th vent
  case Optic_Chiasm:
  case Third_Ventricle:
  case Fourth_Ventricle:
  case Brain_Stem:
  case CC_Posterior:
  case CC_Mid_Posterior:
  case CC_Central:
  case CC_Mid_Anterior:
  case CC_Anterior:
  case CSF_ExtraCerebral:
  case Head_ExtraCerebral:
  case 165: // Skull
  case 259: // SkullApprox
  case 260: // BoneOrAir
  case 262: // Sinus
  case WM_hypointensities:
  case non_WM_hypointensities:
  case CSF:
    id2 = id;
    break;
  default:
    printf("MRIasegContraLatLabel() id=%d unrecognized\n",id);
    break;
  }
  return(id2);
}

/*!
\fn MRI *MRIlrswapAseg(MRI *aseg)
\brief Performs a left-right swap of segmentation labels. Should handle
aseg, aparc+aseg, aparc.a2009s+aseg, and wmparc. If these
segmentations change, then this program (written on 12/19/2011) may
fail. It only changes the voxel values and has no effect on the
volume geometry. See MRIasegContraLatLabel().
\param aseg - segmentation
*/
MRI *MRIlrswapAseg(MRI *aseg)
{
  MRI *asegswap;
  int c, r, s, id, id2;

  asegswap = MRIclone(aseg, NULL);

  for (c = 0; c < aseg->width; c++) {
    for (r = 0; r < aseg->height; r++) {
      for (s = 0; s < aseg->depth; s++) {
        id = MRIgetVoxVal(aseg, c, r, s, 0);
	id2 = MRIasegContraLatLabel(id);
	if(id2 == -1) id2 = id;
        MRIsetVoxVal(asegswap, c, r, s, 0, id2);
      }
    }
  }
  return (asegswap);
}
/*!
\fn MRI *MRIfixAsegWithRibbon(MRI *aseg, MRI *ribbon, MRI *asegfixed)
\brief Compares aseg.mgz to ribbon.mgz and replaces the aseg values
with ribbon values if the aseg is CtxGM or CtxWM or unknown.
\param aseg - aseg.mgz segmentation
\param ribbon - ribbon.mgz segmentation
*/
MRI *MRIfixAsegWithRibbon(MRI *aseg, MRI *ribbon, MRI *asegfixed)
{
  int c, r, s, asegid, ribbonid;

  asegfixed = MRIcopy(aseg, asegfixed);

  for (c = 0; c < aseg->width; c++) {
    for (r = 0; r < aseg->height; r++) {
      for (s = 0; s < aseg->depth; s++) {
        asegid = MRIgetVoxVal(aseg, c, r, s, 0);
        if (asegid == 2 || asegid == 41 || asegid == 3 || asegid == 42 || asegid == 0) {
          ribbonid = MRIgetVoxVal(ribbon, c, r, s, 0);
          MRIsetVoxVal(asegfixed, c, r, s, 0, ribbonid);
        }
      }
    }
  }
  return (asegfixed);
}

/*!
\fn double *ComputeBrainVolumeStats(char *subject)
\brief computes various brain volume statistics and returns them as a
double array.  These include BrainSegVol, BrainSegVolNotVent,
SupraTentVol, SubCortGM, CtxGM, CtxWM, etc.  The hope is that this one
function will be able to consistently define all of these
parameters. Where possible, this function returns values based on
surface-based analysis. It also computes the same values based on
volume-based analysis to check against the surface-based results.
\param subject
*/
double *ComputeBrainVolumeStats(char *subject, char *suff, char *sdir)
{
  char tmpstr[2000];
  char *SUBJECTS_DIR;
  MRI *aseg, *ribbon, *asegfixed, *brainmask;
  MRIS *mris;
  int c, r, s, asegid, ribbonid, asegfixedid, ribbonRead;
  double lhwhitevolTot, rhwhitevolTot, lhpialvolTot, rhpialvolTot;
  double lhCtxGM, rhCtxGM, lhCtxWM, rhCtxWM;
  double lhCtxGMCor, rhCtxGMCor, lhCtxWMCor, rhCtxWMCor;
  double lhCtxGMCount, rhCtxGMCount, lhCtxWMCount, rhCtxWMCount;
  double CCVol, TFFC, eBSVnvSurf;
  double SupraTentVol, SupraTentVolCor, SupraTentVolNotVent, eSTV, eSTVnv;
  double SubCortGMVol, CerebellumVol, CerebellumGMVol, VentChorVol;
  double BrainSegVol, eBSV, BrainSegVolNotVent, MaskVol, VesselVol;
  double OptChiasmVol, CSFVol;
  double *stats = NULL;
  double VoxelVol;


  const char *suffix;
  if (suff == NULL) {
    suffix = "";
  } else {
    suffix = suff;
  }

  if (sdir)
    SUBJECTS_DIR = sdir;
  else
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  sprintf(tmpstr, "%s/%s/surf/lh.white%s", SUBJECTS_DIR, subject, suffix);
  mris = MRISread(tmpstr);
  if (mris == NULL) return (NULL);
  lhwhitevolTot = MRISvolumeInSurf(mris);
  MRISfree(&mris);

  sprintf(tmpstr, "%s/%s/surf/rh.white%s", SUBJECTS_DIR, subject, suffix);
  mris = MRISread(tmpstr);
  if (mris == NULL) return (NULL);
  rhwhitevolTot = MRISvolumeInSurf(mris);
  MRISfree(&mris);

  sprintf(tmpstr, "%s/%s/surf/lh.pial%s", SUBJECTS_DIR, subject, suffix);
  mris = MRISread(tmpstr);
  if (mris == NULL) return (NULL);
  lhpialvolTot = MRISvolumeInSurf(mris);
  MRISfree(&mris);

  sprintf(tmpstr, "%s/%s/surf/rh.pial%s", SUBJECTS_DIR, subject, suffix);
  mris = MRISread(tmpstr);
  if (mris == NULL) return (NULL);
  rhpialvolTot = MRISvolumeInSurf(mris);
  MRISfree(&mris);

  sprintf(tmpstr, "%s/%s/mri/brainmask.mgz", SUBJECTS_DIR, subject);
  brainmask = MRIread(tmpstr);
  if (brainmask == NULL) return (NULL);

  sprintf(tmpstr, "%s/%s/mri/aseg.presurf%s.mgz", SUBJECTS_DIR, subject, suffix);
  if (FileExists(tmpstr) == 0) {
    printf("%s doesn't exist, using old naming conventions\n", tmpstr);
    sprintf(tmpstr, "%s/%s/mri/aseg%s.mgz", SUBJECTS_DIR, subject, suffix);
    if (FileExists(tmpstr) == 0) sprintf(tmpstr, "%s/%s/mri/aseg.long%s.mgz", SUBJECTS_DIR, subject, suffix);
  }
  aseg = MRIread(tmpstr);
  if (aseg == NULL) return (NULL);

  sprintf(tmpstr, "%s/%s/mri/ribbon.mgz", SUBJECTS_DIR, subject);
  if (fio_FileExistsReadable(tmpstr)) {
    ribbon = MRIread(tmpstr);
    if (ribbon == NULL) return (NULL);
    asegfixed = MRIfixAsegWithRibbon(aseg, ribbon, NULL);
    ribbonRead = 1;
  }
  else {
    printf("WARNING: %s does not exist, ribbon based measurements will be inaccurate\n", tmpstr);
    ribbon = aseg;
    asegfixed = aseg;
    ribbonRead = 0;
  }
  VoxelVol = aseg->xsize * aseg->ysize * aseg->zsize;
  printf("ComputeBrainVolumeStats() using version with fixed volume, VoxelVol=%g\n",VoxelVol);

  lhCtxGMCor = 0;
  rhCtxGMCor = 0;
  lhCtxWMCor = 0;
  rhCtxWMCor = 0;
  lhCtxGMCount = 0;
  rhCtxGMCount = 0;
  lhCtxWMCount = 0;
  rhCtxWMCount = 0;
  CCVol = 0;
  SubCortGMVol = 0;
  CerebellumVol = 0;
  CerebellumGMVol = 0;
  VentChorVol = 0;
  BrainSegVol = 0;
  MaskVol = 0;
  VesselVol = 0;
  OptChiasmVol = 0;
  CSFVol = 0;
  TFFC = 0;
  for (c = 0; c < aseg->width; c++) {
    for (r = 0; r < aseg->height; r++) {
      for (s = 0; s < aseg->depth; s++) {
        asegid = MRIgetVoxVal(aseg, c, r, s, 0);
        asegfixedid = MRIgetVoxVal(asegfixed, c, r, s, 0);
        ribbonid = MRIgetVoxVal(ribbon, c, r, s, 0);
        // Corpus Callosum
        if (asegid == 251 || asegid == 252 || asegid == 253 || asegid == 254 || asegid == 255) CCVol += VoxelVol;
        // Correct CtxGM by anything in the ribbon that is not GM, WM, or Unkown in the aseg
        if (ribbonid == 3 && asegid != 3 && asegid != 2 && asegid != 0) lhCtxGMCor += VoxelVol;
        if (ribbonid == 42 && asegid != 42 && asegid != 41 && asegid != 0) rhCtxGMCor += VoxelVol;
        // Correct CtxWM by anything in the WMribbon that is not WM, eg, GM structures.
        // Does not use PVC for subcort GM. Make sure to include hypointensities (77)
        if (ribbonid == 2 && asegfixedid != 2 && asegfixedid != 77 && asegfixedid != 251 && asegfixedid != 252 &&
            asegfixedid != 253 && asegfixedid != 254 && asegfixedid != 255)
          lhCtxWMCor += VoxelVol;
        if (ribbonid == 41 && asegfixedid != 41 && asegfixedid != 77 && asegfixedid != 251 && asegfixedid != 252 &&
            asegfixedid != 253 && asegfixedid != 254 && asegfixedid != 255)
          rhCtxWMCor += VoxelVol;
        // Subcortical GM structures (does not use PVC)
        if (IsSubCorticalGray(asegfixedid)) SubCortGMVol += VoxelVol;
        // Cerebellum GM volume
        if (asegid == Left_Cerebellum_Cortex || asegid == Right_Cerebellum_Cortex) CerebellumGMVol += VoxelVol;
        // Cerebellum (GM+WM) volume
        if (asegid == Left_Cerebellum_Cortex || asegid == Right_Cerebellum_Cortex ||
            asegid == Right_Cerebellum_White_Matter || asegid == Left_Cerebellum_White_Matter)
          CerebellumVol += VoxelVol;
        // Ventricle Volume
        if (asegid == Left_choroid_plexus || asegid == Right_choroid_plexus || asegid == Left_Lateral_Ventricle ||
            asegid == Right_Lateral_Ventricle || asegid == Left_Inf_Lat_Vent || asegid == Right_Inf_Lat_Vent)
          VentChorVol += VoxelVol;
        // 3rd, 4th, 5th, CSF
        if (asegid == Third_Ventricle || asegid == Fourth_Ventricle || asegid == Fifth_Ventricle || asegid == CSF)
          TFFC += VoxelVol;
        // Other
        if (asegid == 30 || asegid == 62) VesselVol += VoxelVol;
        if (asegid == 85) OptChiasmVol += VoxelVol;
        if (asegid == 24) CSFVol += VoxelVol;
        // Total number of voxels in the segmentation. Use fixed to exclude
        // stuff outside of the brain (eg, dura)
        if (asegfixedid != 0 && asegfixedid != Brain_Stem) BrainSegVol += VoxelVol;
        // Total number of voxels in the brainmask
        if (MRIgetVoxVal(brainmask, c, r, s, 0) > 0) MaskVol += VoxelVol;
        // Total number of voxels in the CtxGM, CtxWM of asegfixed
        // This is to check against the surface-based measures
        if (asegfixedid == 3) lhCtxGMCount += VoxelVol;
        if (asegfixedid == 42) rhCtxGMCount += VoxelVol;
        // For CtxWM, include hypointensities. The hypos are not lateralized,
        // so just lateralize them based on column (not perfect, but it is only a check)
        if (asegfixedid == 2 || asegfixedid == 78 || (asegfixedid == 77 && c < 128)) lhCtxWMCount += VoxelVol;
        if (asegfixedid == 41 || asegfixedid == 79 || (asegfixedid == 77 && c >= 128)) rhCtxWMCount += VoxelVol;
      }  // c
    }    // r
  }      // s

  // CtxGM = everything inside pial surface minus everything in white surface
  // minus stuff in the ribbon that is not cortex
  lhCtxGM = lhpialvolTot - lhwhitevolTot - lhCtxGMCor;
  rhCtxGM = rhpialvolTot - rhwhitevolTot - rhCtxGMCor;

  // CtxWM = everything inside of white surface minus stuff that is not WM
  lhCtxWM = lhwhitevolTot - lhCtxWMCor;
  rhCtxWM = rhwhitevolTot - rhCtxWMCor;

  // Add half of CC to each hemi for counting,
  lhCtxWMCount += CCVol / 2.0;
  rhCtxWMCount += CCVol / 2.0;

  // Supratentorial volume is everything inside the pial surface plus
  // stuff that is ouside the surface but still in the ST (eg, hippo, amyg)
  SupraTentVolCor = SupraTentorialVolCorrection(aseg, ribbon);
  SupraTentVol = lhpialvolTot + rhpialvolTot + SupraTentVolCor;
  SupraTentVolNotVent = SupraTentVol - VentChorVol;
  // Estimated STV based - should these be exactly the same? Might depend on how
  // much of CSF and OptChiasm are in or out of the surface.
  // eSTV = lhCtxGM + rhCtxGM + lhCtxWM + rhCtxWM + SubCortGMVol + VentChorVol + VesselVol;
  eSTV = lhCtxGMCount + rhCtxGMCount + lhCtxWMCount + rhCtxWMCount + SubCortGMVol + VentChorVol + VesselVol;
  eSTVnv = lhCtxGMCount + rhCtxGMCount + lhCtxWMCount + rhCtxWMCount + SubCortGMVol + VesselVol;

  // Estimated BrainSegVolume - mostly a comparison between surface- and
  // volume-based
  eBSV = eSTV + CerebellumVol + CSFVol + OptChiasmVol;

  // Surface-based est of brainseg not vent
  eBSVnvSurf = lhCtxGM + rhCtxGM + lhCtxWM + rhCtxWM + SubCortGMVol + CerebellumVol + VesselVol;

  BrainSegVolNotVent = BrainSegVol - VentChorVol - TFFC;

  printf("lhCtxGM: %9.3f %9.3f  diff=%7.1f  pctdiff=%6.3f\n",
         lhCtxGM,
         lhCtxGMCount,
         lhCtxGM - lhCtxGMCount,
         100 * (lhCtxGM - lhCtxGMCount) / lhCtxGM);
  printf("rhCtxGM: %9.3f %9.3f  diff=%7.1f  pctdiff=%6.3f\n",
         rhCtxGM,
         rhCtxGMCount,
         rhCtxGM - rhCtxGMCount,
         100 * (rhCtxGM - rhCtxGMCount) / rhCtxGM);
  printf("lhCtxWM: %9.3f %9.3f  diff=%7.1f  pctdiff=%6.3f\n",
         lhCtxWM,
         lhCtxWMCount,
         lhCtxWM - lhCtxWMCount,
         100 * (lhCtxWM - lhCtxWMCount) / lhCtxWM);
  printf("rhCtxWM: %9.3f %9.3f  diff=%7.1f  pctdiff=%6.3f\n",
         rhCtxWM,
         rhCtxWMCount,
         rhCtxWM - rhCtxWMCount,
         100 * (rhCtxWM - rhCtxWMCount) / rhCtxWM);
  printf("SubCortGMVol  %9.3f\n", SubCortGMVol);
  printf("SupraTentVol  %9.3f (%9.3f) diff=%6.3f pctdiff=%4.3f\n",
         SupraTentVol,
         eSTV,
         SupraTentVol - eSTV,
         100 * (SupraTentVol - eSTV) / SupraTentVol);
  printf("SupraTentVolNotVent  %9.3f (%9.3f) diff=%6.3f pctdiff=%4.3f\n",
         SupraTentVolNotVent,
         eSTVnv,
         SupraTentVolNotVent - eSTVnv,
         100 * (SupraTentVolNotVent - eSTVnv) / SupraTentVolNotVent);
  printf("BrainSegVol  %9.3f (%9.3f) diff=%6.3f pctdiff=%4.3f\n",
         BrainSegVol,
         eBSV,
         BrainSegVol - eBSV,
         100 * (BrainSegVol - eBSV) / BrainSegVol);
  printf("BrainSegVolNotVent  %9.3f (%9.3f) diff=%6.3f pctdiff=%4.3f\n",
         BrainSegVolNotVent,
         eBSVnvSurf,
         BrainSegVolNotVent - eBSVnvSurf,
         100 * (BrainSegVolNotVent - eBSVnvSurf) / BrainSegVolNotVent);

  printf("BrainSegVolNotVent  %9.3f\n", BrainSegVolNotVent);
  printf("CerebellumVol %9.3f\n", CerebellumVol);
  printf("VentChorVol   %9.3f\n", VentChorVol);
  printf("3rd4th5thCSF  %9.3f\n", TFFC);
  printf("CSFVol %9.3f, OptChiasmVol %9.3f\n", CSFVol, OptChiasmVol);
  printf("MaskVol %9.3f\n", MaskVol);

  if (ribbonRead) {
    MRIfree(&ribbon);
    MRIfree(&asegfixed);
  }
  MRIfree(&aseg);
  MRIfree(&brainmask);

  if (!ribbonRead) {
    lhCtxGM = 0;
    rhCtxGM = 0;
  }

  stats = (double *)calloc(16, sizeof(double));
  stats[0] = BrainSegVol;
  stats[1] = BrainSegVolNotVent;
  stats[2] = SupraTentVol;
  stats[3] = SupraTentVolNotVent;
  stats[4] = SubCortGMVol;
  stats[5] = lhCtxGM;
  stats[6] = rhCtxGM;
  stats[7] = lhCtxGM + rhCtxGM;
  stats[8] = SubCortGMVol + lhCtxGM + rhCtxGM + CerebellumGMVol;  // total GM Vol
  stats[9] = lhCtxWM;
  stats[10] = rhCtxWM;
  stats[11] = lhCtxWM + rhCtxWM;
  stats[12] = MaskVol;
  stats[13] = eSTVnv;       // voxel-based supratentorial not vent volume
  stats[14] = eBSVnvSurf;   // surface-based brain  not vent volume
  stats[15] = VentChorVol;  // volume of ventricles + choroid

  return (stats);
}

/*!
  \fn MRI *MRIseg2TissueType(MRI *seg, COLOR_TABLE *ct, MRI *tt)
  \brief Creates a segmentation volume where the segmentation is
  that of tissue type (tissue type info in the ctab).
*/
MRI *MRIseg2TissueType(MRI *seg, COLOR_TABLE *ct, MRI *tt)
{
  int c, r, s, segid;

  if (ct->ctabTissueType == NULL) {
    printf("ERROR: MRIseg2TissueType() ctab tissue type not set\n");
    return (NULL);
  }

  tt = MRIcopy(seg, tt);
  MRIclear(tt);

  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        segid = MRIgetVoxVal(seg, c, r, s, 0);
        if (ct->entries[segid] == NULL) {
          printf("ERROR: MRIseg2TTypeMap() no entry for seg %d\n", segid);
          return (NULL);
        }
        if (ct->entries[segid]->TissueType < 0) {
          printf("ERROR: MRIseg2TTypeMap() tissue type for seg %d %s not set\n", segid, ct->entries[segid]->name);
          return (NULL);
        }
        MRIsetVoxVal(tt, c, r, s, 0, ct->entries[segid]->TissueType);
      }
    }
  }
  return (tt);
}

/*!
  \fn MRI *MRIextractTissueTypeSeg(MRI *seg, COLOR_TABLE *ct, int tt, MRI
  *ttseg)
  \brief Creates a volume with only the segmentations in the given tissue type.
  Tissue type info in the ctab.
*/
MRI *MRIextractTissueTypeSeg(MRI *seg, COLOR_TABLE *ct, int tt, MRI *ttseg)
{
  int c, r, s, segid;

  if (ct->ctabTissueType == NULL) {
    printf("ERROR: MRIextractTissueTypeSeg() ctab tissue type not set\n");
    return (NULL);
  }
  if (tt >= ct->ctabTissueType->nentries) {
    printf("ERROR: MRIextractTissueTypeSeg() tissue type %d exceeds or equals number of tissue types %d\n",
           tt,
           ct->ctabTissueType->nentries);
    return (NULL);
  }
  ttseg = MRIcopy(seg, ttseg);
  MRIclear(ttseg);

  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        segid = MRIgetVoxVal(seg, c, r, s, 0);
        if (ct->entries[segid] == NULL) {
          printf("ERROR: MRIseg2TTypeMap() no entry for seg %d\n", segid);
          return (NULL);
        }
        if (ct->entries[segid]->TissueType < 0) {
          printf("ERROR: MRIseg2TTypeMap() tissue type for seg %d %s not set\n", segid, ct->entries[segid]->name);
          return (NULL);
        }
        if (ct->entries[segid]->TissueType == tt) MRIsetVoxVal(ttseg, c, r, s, 0, segid);
      }
    }
  }
  return (ttseg);
}

/*!
  \fn int CheckSegTissueType(MRI *seg, COLOR_TABLE *ct)
  \brief Make sure that the each segmentation has a tissue type.
  Tissue type info in the ctab. Counts number of voxels for
  each entry.
  \return 0 if no error, 1 if error
*/
int CheckSegTissueType(MRI *seg, COLOR_TABLE *ct)
{
  int c, r, s, n, segid, err;

  err = 1;
  if (ct->ctabTissueType == NULL) {
    printf("ERROR: CheckSegTissueType() ctab tissue type not set\n");
    return (err);
  }

  for (n = 0; n < ct->nentries; n++)
    if (ct->entries[n]) ct->entries[n]->count = 0;

  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        segid = MRIgetVoxVal(seg, c, r, s, 0);
        if (ct->entries[segid] == NULL) {
          printf("ERROR: CheckSegTissueType() no entry for seg %d\n", segid);
          return (err);
        }
        if (ct->entries[segid]->TissueType < 0) {
          printf("ERROR: CheckSegTissueType() tissue type for seg %d %s not set\n", segid, ct->entries[segid]->name);
          return (err);
        }
        ct->entries[segid]->count++;
      }
    }
  }
  err = 0;
  return (err);
}

/*!
  \fn MRI **MRIdilateSegWithinTT(MRI *seg, int nDils, COLOR_TABLE *ct)
  \brief Returns an array of MRI structures of length equal to the number
  of tissue types. Each MRI contains only the segmentation structures
  of the given tissue type. The segmentations have been dilated to
  fill in voxels from the masked out tissue types. This technique
  is used for GTM partial volume correction. Tissue type info in the ctab.
*/
MRI **MRIdilateSegWithinTT(MRI *seg, int nDils, COLOR_TABLE *ct, MRI **r)
{
  MRI *segtt = NULL;
  int nc, tt;
  // char tmpstr[1000];

  if (ct->ctabTissueType == NULL) {
    printf("ERROR: MRIdilateSegWithinTT() ctab tissue type not set\n");
    return (NULL);
  }

  if (r == NULL) r = (MRI **)calloc(sizeof(MRI *), ct->ctabTissueType->nentries - 1);
  for (tt = 1; tt < ct->ctabTissueType->nentries; tt++) {
    // printf("tt = %d\n",tt);
    segtt = MRIextractTissueTypeSeg(seg, ct, tt, segtt);
    if (segtt == NULL) return (NULL);
    r[tt - 1] = MRIdilateSegmentation(segtt, NULL, nDils, NULL, 0, 0, &nc);
    if (r[tt - 1] == NULL) return (NULL);
    // sprintf(tmpstr,"seg.dil%d.tt%d.mgh",nDils,tt);
    // MRIwrite(r[tt-1],tmpstr);
  }

  MRIfree(&segtt);
  return (r);
}

/*!
  \fn int Seg2NbrNonBrainWrapper(char *subject, char *segname, char *ctab, char
  *statname, double threshmm)
  \brief Wrapper around Seg2NbrNonBrain() that operates on the FS directory
  structure and performs
  GTMdefaultSegReplacmentList(). Writes out stats.
*/
int Seg2NbrNonBrainWrapper(char *subject, char *segname, COLOR_TABLE *ctab, char *statname, double threshmm)
{
  char *SUBJECTS_DIR, tmpstr[2000];
  MRI *seg, *mritmp;
  int nReplace, SrcReplace[1000], TrgReplace[1000];
  SEGSTAT *segstat;
  FILE *fp;

  printf("Seg2NbrNonBrainWrapper()  %s %s %s %g\n", subject, segname, statname, threshmm);
  fflush(stdout);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(tmpstr, "%s/%s/mri/%s", SUBJECTS_DIR, subject, segname);
  seg = MRIread(tmpstr);
  if (seg == NULL) exit(1);

  printf("Replacing\n");
  fflush(stdout);
  GTMdefaultSegReplacmentList(&nReplace, &SrcReplace[0], &TrgReplace[0]);
  mritmp = MRIreplaceList(seg, SrcReplace, TrgReplace, nReplace, NULL, NULL);
  MRIfree(&seg);
  seg = mritmp;

  printf("Computing stats\n");
  fflush(stdout);
  segstat = Seg2NbrNonBrain(seg, ctab, 3);

  sprintf(tmpstr, "%s/%s/stats/%s", SUBJECTS_DIR, subject, statname);
  fp = fopen(tmpstr, "w");
  fprintf(fp, "# Non-brain tissue volume near each segmentation\n");
  fprintf(fp, "# Seg2NbrNonBrainWrapper()  %s %s %s %g\n", subject, segname, statname, threshmm);
  fprintf(fp, "# seg voxsize  %g %g %g\n", seg->xsize, seg->ysize, seg->zsize);
  PrintSegStat(fp, segstat);
  fclose(fp);

  MRIfree(&seg);

  printf("Seg2NbrNonBrainWrapper() done\n");
  fflush(stdout);

  return (0);
}

/*!
  \fn SEGSTAT *Seg2NbrNonBrain(MRI *seg, COLOR_TABLE *ctab, double threshmm)
  \brief Computes the volume of non-brain voxels with in threshmm of a
  given brain tissue segmentation. ctab is used to determine which seg
  numbers are brain and non-brain. If ctab is NULL, then the default
  is used. A SEGSTAT table is returned. The mean intensity field
  contains the volume of the actual segmentation for a given structure
  so as to compare against other measures as a double check. Can use
  PrintSegStat(fp, segstat). Note that a single non-brain voxel may
  be counted multiple times.
*/
SEGSTAT *Seg2NbrNonBrain(MRI *seg, COLOR_TABLE *ctab, double threshmm)
{
  int c, r, s, cB, rB, sB, nthseg, segno, segnoB, FreeCTab;
  int *segnolist, *count, *segcount, nsegs;
  double threshvox, d2, dc2, dr2, dc, dr, ds, voxsize, threshmm2;
  int cBmin, cBmax, rBmin, rBmax, sBmin, sBmax;
  SEGSTAT *segstat;
  MRI *hitmap;

  FreeCTab = 0;
  if (ctab == NULL) {
    ctab = TissueTypeSchema(NULL, "default-jan-2014+head");
    FreeCTab = 1;
  }

  segnolist = MRIsegIdListNot0(seg, &nsegs, 0);
  count = (int *)calloc(nsegs, sizeof(int));
  segcount = (int *)calloc(nsegs, sizeof(int));

  threshmm2 = threshmm * threshmm;  // distance threshold squared

  // for the cube bounding box
  threshvox = ceil(threshmm / (MIN(MIN(seg->xsize, seg->ysize), seg->zsize)));

  voxsize = seg->xsize * seg->ysize * seg->zsize;
  printf("threshmm = %g, threshvox = %lf voxsize = %lf\n", threshmm, threshvox, voxsize);

  printf("Allocating hitmap\n");
  fflush(stdout);
  hitmap = MRIallocSequence(seg->width, seg->height, seg->depth, MRI_UCHAR, nsegs);
  if (hitmap == NULL) exit(1);
  printf("Copying header\n");
  fflush(stdout);
  MRIcopyHeader(seg, hitmap);
  printf("Copying pulse params\n");
  fflush(stdout);
  MRIcopyPulseParameters(seg, hitmap);

  printf("Starting loop\n");
  fflush(stdout);
  // Go through each source segmentation voxel
  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        segno = MRIgetVoxVal(seg, c, r, s, 0);
        if (segno == 0) continue;

        // If the source voxel is not Ctx, SubCtxGM, or WM, skip
        if (ctab->entries[segno]->TissueType > 3) continue;

        // Get the nthseg for this segno
        for (nthseg = 0; nthseg < nsegs; nthseg++)
          if (segno == segnolist[nthseg]) break;
        segcount[nthseg]++;

        // cube around the CRS source voxel
        cBmin = nint(MAX(0, c - threshvox));
        cBmax = nint(MIN(c + threshvox, seg->width - 1));
        rBmin = nint(MAX(0, r - threshvox));
        rBmax = nint(MIN(r + threshvox, seg->height - 1));
        sBmin = nint(MAX(0, s - threshvox));
        sBmax = nint(MIN(s + threshvox, seg->depth - 1));

        // Loop through the cube
        for (cB = cBmin; cB < cBmax; cB++) {
          dc = seg->xsize * (c - cB);
          dc2 = (dc * dc);
          for (rB = rBmin; rB < rBmax; rB++) {
            dr = seg->ysize * (r - rB);
            dr2 = (dr * dr);
            for (sB = sBmin; sB < sBmax; sB++) {
              segnoB = MRIgetVoxVal(seg, cB, rB, sB, 0);
              // Only consider target voxels that are NOT Ctx, SubCtxGM, or WM
              if (ctab->entries[segnoB]->TissueType <= 3) continue;
              // Skip if this voxel has already been counted for this segno
              if (MRIgetVoxVal(hitmap, cB, rB, sB, nthseg) > 0.5) continue;
              // Compute distance-squared from source to target
              ds = seg->zsize * (s - sB);
              d2 = dc2 + dr2 + ds * ds;
              if (d2 > threshmm2) continue;
              count[nthseg]++;
              MRIsetVoxVal(hitmap, cB, rB, sB, nthseg, 1);  // mark as hit
            }                                               // sB
          }                                                 // rB
        }                                                   // cB

      }  // s
    }    // r
  }      // c

  // This can be used to verify
  // MRIwrite(hitmap,"hitmap.mgh");
  printf("Freeing hitmap\n");
  fflush(stdout);
  MRIfree(&hitmap);

  segstat = (SEGSTAT *)calloc(sizeof(SEGSTAT), 1);
  segstat->nentries = nsegs;
  segstat->entry = (STATSUMENTRY *)calloc(sizeof(STATSUMENTRY), nsegs);
  segstat->UseName = 1;
  segstat->IsSurf = 0;
  segstat->DoIntensity = 1;
  segstat->InIntensityName = "SegVolInBrain";
  segstat->InIntensityUnits = "mm3";
  for (nthseg = 0; nthseg < nsegs; nthseg++) {
    segno = segnolist[nthseg];
    // printf("%3d %4d %d %-25s
    // %5d\n",nthseg,segno,ctab->entries[segno]->TissueType,ctab->entries[segno]->name,count[nthseg]);
    sprintf(segstat->entry[nthseg].name, "%s", ctab->entries[segno]->name);
    segstat->entry[nthseg].id = segno;
    segstat->entry[nthseg].nhits = count[nthseg];
    segstat->entry[nthseg].vol = count[nthseg] * voxsize;
    segstat->entry[nthseg].mean = segcount[nthseg] * voxsize;
  }

  free(segnolist);
  if (FreeCTab) CTABfree(&ctab);
  return (segstat);
}

/* eof */
