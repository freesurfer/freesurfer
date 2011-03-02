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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:42 $
 *    $Revision: 1.9 $
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


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "error.h"

/* see ch notebook 2 */

#include "cma.h"

#ifndef TRUE
#define TRUE (1==1)
#endif

#ifndef FALSE
#define FALSE (0==1)
#endif

extern int errno;

int CMAfreeOutlineField(CMAoutlineField **of)
{

  CMAoutlineField *ofp;
  int i;

  ofp = *of;

  if (ofp->claim_field != NULL)
  {
    for (i = 0;i < ofp->height;i++)
    {
      if (ofp->claim_field[i] != NULL)
        free(ofp->claim_field[i]);
    }
    free(ofp->claim_field);
  }

  if (ofp->fill_field != NULL)
  {
    for (i = 0;i < ofp->height;i++)
    {
      if (ofp->fill_field[i] != NULL)
        free(ofp->fill_field[i]);
    }
    free(ofp->fill_field);
  }

  if (ofp->outline_points_field != NULL)
  {
    for (i = 0;i < ofp->height;i++)
    {
      if (ofp->outline_points_field[i] != NULL)
        free(ofp->outline_points_field[i]);
    }
    free(ofp->outline_points_field);
  }

  *of = NULL;

  return(NO_ERROR);

} /* end CMAfreeOutlineField() */

CMAoutlineField *CMAoutlineFieldAlloc(int width, int height)
{

  CMAoutlineField *of;
  int i;

  of = (CMAoutlineField *)malloc(sizeof(CMAoutlineField));
  if (of == NULL)
    ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating structure"));

  of->claim_field = NULL;
  of->fill_field = NULL;
  of->outline_points_field = NULL;
  of->width = width;
  of->height = height;

  of->claim_field = (CMAoutlineClaim **)malloc(height * sizeof(CMAoutlineClaim *));
  if (of->claim_field == NULL)
  {
    CMAfreeOutlineField(&of);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating claim field"));
  }
  memset(of->claim_field, 0x00, height * sizeof(CMAoutlineClaim *));

  //of->fill_field = (unsigned char **)malloc(height * sizeof(unsigned char *));
  of->fill_field = (short **)malloc(height * sizeof(short *));
  if (of->fill_field == NULL)
  {
    CMAfreeOutlineField(&of);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating fill field"));
  }
  //memset(of->fill_field, 0x00, height * sizeof(unsigned char *));
  memset(of->fill_field, 0x00, height * sizeof(short *));

  of->outline_points_field = (unsigned char **)malloc(height * sizeof(unsigned char *));
  if (of->outline_points_field == NULL)
  {
    CMAfreeOutlineField(&of);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating outline points field"));
  }
  memset(of->outline_points_field, 0x00, height * sizeof(unsigned char *));

  for (i = 0;i < height;i++)
  {

    of->claim_field[i] = (CMAoutlineClaim *)malloc(width * sizeof(CMAoutlineClaim));
    if (of->claim_field[i] == NULL)
    {
      CMAfreeOutlineField(&of);
      ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating claim field"));
    }
    memset(of->claim_field[i], 0x00, width * sizeof(CMAoutlineClaim));

    //of->fill_field[i] = (unsigned char *)malloc(width * sizeof(unsigned char));
    of->fill_field[i] = (short *)malloc(width * sizeof(short));
    if (of->fill_field[i] == NULL)
    {
      CMAfreeOutlineField(&of);
      ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating fill field"));
    }
    //memset(of->fill_field[i], 0x00, width * sizeof(unsigned char));
    memset(of->fill_field[i], 0x00, width * sizeof(short));

    of->outline_points_field[i] = (unsigned char *)malloc(width * sizeof(unsigned char));
    if (of->outline_points_field[i] == NULL)
    {
      CMAfreeOutlineField(&of);
      ErrorReturn(NULL, (ERROR_NOMEMORY, "CMAoutlineFieldAlloc(): error allocating outline points field"));
    }
    memset(of->outline_points_field[i], 0x00, width * sizeof(unsigned char));

  }

  return(of);

} /* end CMAoutlineFieldAlloc() */

int CMAclearFillField(CMAoutlineField *field)
{

  int i;

  for (i = 0;i < field->height;i++)
    //memset(field->fill_field[i], 0x00, field->width * sizeof(unsigned char));
    memset(field->fill_field[i], 0x00, field->width * sizeof(short));

  return(NO_ERROR);

} /* end CMAclearFillField() */

/* fills with CMA_FILL_INTERIOR to non-zero values */
int CMAfill(CMAoutlineField *field, short seed_x, short seed_y)
{

  if (seed_x < 0 || seed_x >= field->width)
    return(NO_ERROR);

  if (seed_y < 0 || seed_y >= field->height)
    return(NO_ERROR);

  if (field->fill_field[seed_y][seed_x] != 0)
    return(NO_ERROR);

  field->fill_field[seed_y][seed_x] = CMA_FILL_INTERIOR;

  CMAfill(field, seed_x + 1, seed_y    );
  CMAfill(field, seed_x - 1, seed_y    );
  CMAfill(field, seed_x    , seed_y + 1);
  CMAfill(field, seed_x    , seed_y - 1);

  return(NO_ERROR);

} /* end CMAfill() */

int CMAclaimPoints(CMAoutlineField *field, short label, short *points, int n_points, short seed_x, short seed_y)
{

  int i, j;
  short x, y;

  if (label < 0 || label > MAX_CMA_LABEL)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): label out of range (label is %d, MAX_CMA_LABEL is %d)", label, MAX_CMA_LABEL));

  if (seed_x < 0 || seed_x >= field->width)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): seed point out of range (seed_x = %d, field width = %d)", seed_x, field->width));
  if (seed_y < 0 || seed_y >= field->height)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): seed point out of range (seed_y = %d, field height = %d)", seed_y, field->height));

  CMAclearFillField(field);

  for (i = 0;i < n_points;i++)
  {
    x = points[2*i];
    y = points[2*i+1];
    if (x < 0 || x >= field->width)
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): outline point out of range (x)"));
    if (y < 0 || y >= field->height)
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CMAclaimPoints(): outline point out of range (y)"));
    field->fill_field[y][x] = CMA_FILL_OUTLINE;
    field->outline_points_field[y][x] = 1;
  }

  CMAfill(field, seed_x, seed_y);

  for (i = 0;i < field->width;i++)
  {
    for (j = 0;j < field->height;j++)
    {

      if (field->fill_field[j][i] == CMA_FILL_INTERIOR)
      {
        field->claim_field[j][i].n_claims = 1;
        field->claim_field[j][i].interior_claim_flag = TRUE;
        field->claim_field[j][i].claim_labels[0] = label;
      }

      if (field->fill_field[j][i] == CMA_FILL_OUTLINE)
      {
        if (field->claim_field[j][i].n_claims < MAX_OUTLINE_CLAIMS)
        {
          field->claim_field[j][i].claim_labels[field->claim_field[j][i].n_claims] = label;
          field->claim_field[j][i].n_claims++;
          field->claim_field[j][i].interior_claim_flag = FALSE;
        }
      }

    }
  }

  return(NO_ERROR);

} /* end CMAclaimPoints() */

int CMAvalueClaims(CMAoutlineClaim *claim)
{

  int i;

  for (i = 0;i < MAX_OUTLINE_CLAIMS;i++)
    claim->claim_values[i] = 0.0;

  if (claim->n_claims == 0)
  {
    claim->no_label_claim = 0.5;
  }
  else if (claim->n_claims == 1)
  {
    if (claim->interior_claim_flag == 1)
      claim->claim_values[0] = 1.0;
    else
    {
      claim->claim_values[0] = 0.5;
      claim->no_label_claim = 0.5;
    }
  }
  else
  {
    float ct = 1.0 / (float)claim->n_claims;
    for (i = 0;i < claim->n_claims;i++)
      claim->claim_values[i] = ct;
  }

  return(NO_ERROR);

} /* end CMAvalueClaims() */

int CMAvalueAllClaims(CMAoutlineField *field)
{

  int i, j;

  for (i = 0;i < field->width;i++)
    for (j = 0;j < field->height;j++)
      CMAvalueClaims(&(field->claim_field[j][i]));

  return(NO_ERROR);

} /* end CMAvalueAllClaims() */

int CMAaddWeightedTotals(CMAoutlineClaim *claim, float weight, float *claim_totals)
{

  int i;

  /* we've checked the label range in CMAclaimPoints() */

  for (i = 0;i < claim->n_claims;i++)
    claim_totals[claim->claim_labels[i]] += weight * claim->claim_values[i];

  claim_totals[0] += weight * claim->no_label_claim;

  return(NO_ERROR);

} /* end CMAaddWeightedTotals() */

/* returns the label with the greatest claim to a given pixel */
short CMAtotalClaims(CMAoutlineField *field, int x, int y)
{

  float claim_totals[MAX_CMA_LABEL+1];
  float best_claim;
  short best_index;
  int i;

  if (x < 0 || x >= field->width)
    ErrorReturn(-1, (ERROR_BADPARM, "CMAtotalClaims(): x out of range (x = %d, field->width = %d)", x, field->width));

  if (y < 0 || y >= field->height)
    ErrorReturn(-1, (ERROR_BADPARM, "CMAtotalClaims(): y out of range (y = %d, field->height = %d)", y, field->height));

  for (i = 0;i < MAX_CMA_LABEL;i++)
    claim_totals[i] = 0.0;

  /* center point (x, y checks done above) */
  CMAaddWeightedTotals(&(field->claim_field[y][x]), OTL_CLAIM_WEIGHT_CENTER, claim_totals);

  /* immediately adjoining points */
  if (x-1 >= 0)
    CMAaddWeightedTotals(&(field->claim_field[y  ][x-1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x+1 < field->width)
    CMAaddWeightedTotals(&(field->claim_field[y  ][x+1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (y-1 >= 0)
    CMAaddWeightedTotals(&(field->claim_field[y-1][x  ]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (y+1 < field->height)
    CMAaddWeightedTotals(&(field->claim_field[y+1][x  ]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);

  /* diagonally adjoining points */
  if (x-1 >= 0 && y-1 >= 0)
    CMAaddWeightedTotals(&(field->claim_field[y-1][x-1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x-1 >= 0 && y+1 < field->height)
    CMAaddWeightedTotals(&(field->claim_field[y+1][x-1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x+1 < field->width && y-1 >= 0)
    CMAaddWeightedTotals(&(field->claim_field[y-1][x+1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);
  if (x+1 < field->width && y+1 < field->height)
    CMAaddWeightedTotals(&(field->claim_field[y+1][x+1]), OTL_CLAIM_WEIGHT_SIDE, claim_totals);

  /* find the highest claim and its index */
  best_claim = claim_totals[0];
  best_index = 0;
  for (i = 1;i <= MAX_CMA_LABEL;i++)
  {
    if (claim_totals[i] > best_claim)
    {
      best_claim = claim_totals[i];
      best_index = i;
    }
  }

  return(best_index);

} /* end CMAtotalClaims() */

int CMAassignLabels(CMAoutlineField *field)
{

  int i, j;

  CMAclearFillField(field);

  CMAvalueAllClaims(field);

  for (i = 0;i < field->width;i++)
    for (j = 0;j < field->height;j++)
      field->fill_field[j][i] = CMAtotalClaims(field, i, j);

  return(NO_ERROR);

} /* end CMAassignLabels() */

int CMAzeroOutlines(CMAoutlineField *field)
{

  int i, j;

  for (i = 0;i < field->width;i++)
    for (j = 0;j < field->height;j++)
      if (field->outline_points_field[j][i])
        field->fill_field[j][i] = 0;

  return(NO_ERROR);

} /* end CMAzeroOutlines() */

#include "mrisurf.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
int
insert_ribbon_into_aseg(MRI *mri_src_aseg, MRI *mri_aseg, 
                        MRI_SURFACE *mris_white, MRI_SURFACE *mris_pial, 
                        int hemi)
{
  MRI  *mri_ribbon, *mri_white ;
  int  x, y, z, gm_label, wm_label, label, nbr_label, dont_change ;

  if (mri_src_aseg != mri_aseg)
    mri_aseg = MRIcopy(mri_src_aseg, mri_aseg) ;

  gm_label = hemi == LEFT_HEMISPHERE ? 
    Left_Cerebral_Cortex : Right_Cerebral_Cortex ;
  wm_label = hemi == LEFT_HEMISPHERE ? 
    Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;

  mri_white = MRIclone(mri_aseg, NULL) ;
  mri_ribbon = MRISribbon(mris_white, mris_pial, mri_aseg, NULL) ;
  MRISfillInterior(mris_white, mri_aseg->xsize, mri_white) ;

  for (x = 0 ; x < mri_aseg->width ; x++)
    for (y = 0 ; y < mri_aseg->height ; y++)
      for (z = 0 ; z < mri_aseg->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = nint(MRIgetVoxVal(mri_aseg, x, y, z, 0)) ;
        if (nint(MRIgetVoxVal(mri_ribbon, x, y, z, 0)) == 255)  // in ribbon, set to GM
        {
          if (IS_CORTEX(label) || IS_WHITE_CLASS(label))
          {
            int  xi, yi, zi, xk, yk, zk ;
            // check to make sure we are really in cortex
            for (dont_change = 0, xk = -1 ; xk <= 1 ; xk++)
            {
              xi = mri_aseg->xi[xk+x] ;
              for (yk = -1 ; yk <= 1 ; yk++)
              {
                yi = mri_aseg->yi[yk+y] ;
                for (zk = -1 ; zk <= 1 ; zk++)
                {
                  zi = mri_aseg->zi[zk+z] ;
                  nbr_label = (int)MRIgetVoxVal(mri_aseg, xi, yi, zi, 0) ;
                  switch (nbr_label)
                  {
                  default:
                    break ;
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
                  case Left_Thalamus_Proper:
                  case Right_Thalamus_Proper:
                  case Left_choroid_plexus:
                  case Right_choroid_plexus:
                  case CC_Posterior:
                  case CC_Mid_Posterior:
                  case CC_Central:
                  case CC_Mid_Anterior:
                  case CC_Anterior:
                    dont_change = 1 ;
                    break ;
                  }
                }
              }
            }
            if (dont_change == 0)
              MRIsetVoxVal(mri_aseg, x, y, z, 0, gm_label) ;
          }
        }
        else  // not in ribbon
        {
          if (MRIgetVoxVal(mri_white, x, y, z, 0) > 0)  // inside white surface - disambiguate
          {
            if (label == gm_label)  // gm inside white surface should be wm
              MRIsetVoxVal(mri_aseg, x, y, z, 0, wm_label) ;
          }
          else if (label == gm_label)  // gm outside ribbon should be unknown
            MRIsetVoxVal(mri_aseg, x, y, z, 0, Unknown) ;
        }
      }


  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_ribbon, "r.mgz") ;
    MRIwrite(mri_white, "w.mgz");
  }
  MRIfree(&mri_ribbon) ; MRIfree(&mri_white) ;
  return(NO_ERROR) ;
}

/*!
\fn int IsSubCorticalGray(int SegId)
\brief Returns a 1 if the given seg is subcortical gray
\param SegId - segmentation id number
*/
int IsSubCorticalGray(int SegId)
{
  if(SegId == Left_Thalamus) return(1);
  if(SegId == Right_Thalamus) return(1);
  if(SegId == Left_Thalamus_Proper) return(1);
  if(SegId == Right_Thalamus_Proper) return(1);
  if(SegId == Left_Caudate) return(1);
  if(SegId == Right_Caudate) return(1);
  if(SegId == Left_Putamen ) return(1);
  if(SegId == Right_Putamen ) return(1);
  if(SegId == Left_Pallidum) return(1);
  if(SegId == Right_Pallidum) return(1);
  if(SegId == Brain_Stem) return(1);
  if(SegId == Left_Hippocampus) return(1);
  if(SegId == Right_Hippocampus) return(1);
  if(SegId == Left_Amygdala) return(1);
  if(SegId == Right_Amygdala) return(1);
  if(SegId == Left_Accumbens_area) return(1);
  if(SegId == Right_Accumbens_area) return(1);
  if(SegId == Left_VentralDC) return(1);
  if(SegId == Right_VentralDC) return(1);
  if(SegId == Left_Substancia_Nigra) return(1);
  if(SegId == Right_Substancia_Nigra) return(1);
  if(SegId == Left_Cerebellum_Cortex) return(1);
  if(SegId == Right_Cerebellum_Cortex) return(1);
  return(0);
}

/*!
\fn double SupraTentorialVolCorrection(MRI *aseg, MRI *ribbon)
\brief Returns the volume of supratentorial structures that do not
fall inside the pial surface or are cut by the pial surface.  The idea
is that the volume of everything in the pial surface can be computed
using the surface and that this function can be used to compute
everything else. Note that there is no partial volume correction.
\param aseg - aseg.mgz or aparc+aseg.mgz
\param ribbon is the ribbon.mgz, which has non-zero values for
everything inside the pial surf.
*/
double SupraTentorialVolCorrection(MRI *aseg, MRI *ribbon)
{
  int c,r,s,SegId;
  double vol = 0;
  int RibbonVal, VoxSize;

  VoxSize = aseg->xsize * aseg->ysize * aseg->zsize;
  for(c=0; c < aseg->width; c++){
    for(r=0; r < aseg->height; r++){
      for(s=0; s < aseg->depth; s++){

	// If this voxel is inside the pial, then skip it because it
	// will be part of the surface-based volume measure
	RibbonVal = MRIgetVoxVal(ribbon,c,r,s,0);
	if(RibbonVal == 0) continue;

	// If it gets here, it means that the voxel was not within
	// the pial surface. It could be in a structure that should
	// be part of the supratentorium.

	// These are midline, medial wall, or unknown structures
	// that the pial could cut through.
	SegId = MRIgetVoxVal(aseg,c,r,s,0);
	if(SegId == Left_Lateral_Ventricles) vol += VoxSize;
	if(SegId == Right_Lateral_Ventricles) vol += VoxSize;
	if(SegId == Left_choroid_plexus) vol += VoxSize;
	if(SegId == Right_choroid_plexus) vol += VoxSize;
	if(SegId == Left_Inf_Lat_Vent) vol += VoxSize;
	if(SegId == Right_Inf_Lat_Vent) vol += VoxSize;
	if(SegId == WM_hypointensities) vol += VoxSize;
	if(SegId == Left_WM_hypointensities) vol += VoxSize;
	if(SegId == Right_WM_hypointensities) vol += VoxSize;
	if(SegId == Left_Thalamus_Proper) vol += VoxSize;
	if(SegId == Right_Thalamus_Proper) vol += VoxSize;
	if(SegId == Left_Thalamus) vol += VoxSize;
	if(SegId == Right_Thalamus) vol += VoxSize;
	if(SegId == CC_Posterior) vol += VoxSize;
	if(SegId == CC_Mid_Posterior) vol += VoxSize;
	if(SegId == CC_Central) vol += VoxSize;
	if(SegId == CC_Mid_Anterior) vol += VoxSize;
	if(SegId == CC_Anterior) vol += VoxSize;
	if(SegId == Left_VentralDC) vol += VoxSize;
	if(SegId == Right_VentralDC) vol += VoxSize;

	// These are unlikely to have the pial surface cut through
	// them, but no harm to include them
	if(SegId == Left_Caudate) vol += VoxSize;
	if(SegId == Right_Caudate) vol += VoxSize;
	if(SegId == Left_Putamen ) vol += VoxSize;
	if(SegId == Right_Putamen ) vol += VoxSize;
	if(SegId == Left_Pallidum) vol += VoxSize;
	if(SegId == Right_Pallidum) vol += VoxSize;
	if(SegId == Left_Hippocampus) vol += VoxSize;
	if(SegId == Right_Hippocampus) vol += VoxSize;
	if(SegId == Left_Amygdala) vol += VoxSize;
	if(SegId == Right_Amygdala) vol += VoxSize;
	if(SegId == Left_Accumbens_area) vol += VoxSize;
	if(SegId == Right_Accumbens_area) vol += VoxSize;
      }
    }
  }
  return(vol);
}






/* eof */
