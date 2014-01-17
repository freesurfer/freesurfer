/**
 * @file  fcd.c
 * @brief I/O and analysis algorithms for FCDs (focal cortical dysplasias)
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2014/01/17 15:44:00 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include "fcd.h"
#include "error.h"
#include "diag.h"
#include "cma.h"
#include "const.h"
#include "macros.h"
#include "mrisegment.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

static int
most_frequent_label(MRI *mri_seg, MRI_SEGMENT *mseg)
{
  int  label_counts[MAX_CMA_LABELS], i, max_count, max_label, label ;

  memset(label_counts, 0, sizeof(label_counts)) ;
  for (i = max_count = max_label = 0 ; i < mseg->nvoxels ; i++)
  {
    label = MRIgetVoxVal(mri_seg, mseg->voxels[i].x, mseg->voxels[i].y, mseg->voxels[i].z, 0) ;
    label_counts[label]++ ;
  }
  for (i = max_count = max_label = 0 ; i < MAX_CMA_LABELS ; i++)
  {
    if (label_counts[i] > max_count)
    {
      max_count = label_counts[i] ;
      max_label = i ;
    }
  }
  return(max_label) ;
}


FCD_DATA   *
FCDloadData(char *sdir, char *subject)
{
  FCD_DATA *fcd ;
  char     fname[STRLEN] ;

  fcd = (FCD_DATA *)calloc(1, sizeof(FCD_DATA)) ;

  sprintf(fname, "%s/%s/surf/lh.white", sdir, subject) ;
  fcd->mris_lh = MRISread(fname) ;
  if (fcd->mris_lh == NULL)
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname) ;
  MRISsaveVertexPositions(fcd->mris_lh, WHITE_VERTICES) ;

  sprintf(fname, "%s/%s/surf/rh.white", sdir, subject) ;
  fcd->mris_rh = MRISread(fname) ;
  if (fcd->mris_rh == NULL)
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname) ;
  MRISsaveVertexPositions(fcd->mris_rh, WHITE_VERTICES) ;

  sprintf(fname, "%s/%s/mri/aparc+aseg.mgz", sdir, subject) ;
  fcd->mri_aseg = MRIread(fname) ;
  if (fcd->mri_aseg == NULL)
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname) ;

  sprintf(fname, "%s/%s/mri/flair.reg.norm.mgz", sdir, subject) ;
  fcd->mri_flair = MRIread(fname) ;
  if (fcd->mri_flair == NULL)
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname) ;

  sprintf(fname, "%s/%s/mri/norm.mgz", sdir, subject) ;
  fcd->mri_norm = MRIread(fname) ;
  if (fcd->mri_norm == NULL)
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname) ;

  fcd->mri_thickness_increase = MRIcloneDifferentType(fcd->mri_aseg, MRI_FLOAT) ;
  fcd->mri_thickness_decrease = MRIcloneDifferentType(fcd->mri_aseg, MRI_FLOAT) ;

  sprintf(fname, "%s/%s/surf/lh.rh.thickness.smooth0.mgz", sdir, subject) ;
  fcd->rh_thickness_on_lh = MRIread(fname) ;
  if (fcd->mri_aseg == NULL)
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname) ;

  sprintf(fname, "%s/%s/surf/rh.thickness.mgz", sdir, subject) ;
  fcd->rh_thickness_on_rh = MRIread(fname) ;
  if (fcd->rh_thickness_on_rh == NULL)
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname) ;


  sprintf(fname, "%s/%s/surf/lh.thickness.mgz", sdir, subject) ;
  fcd->lh_thickness_on_lh = MRIread(fname) ;
  if (fcd->lh_thickness_on_lh == NULL)
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname) ;

  sprintf(fname, "%s/%s/surf/rh.lh.thickness.mgz", sdir, subject) ;
  fcd->lh_thickness_on_rh = MRIread(fname) ;
  if (fcd->lh_thickness_on_rh == NULL)
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname) ;

  
  return(fcd) ;
}

int
FCDcomputeThicknessLabels(FCD_DATA *fcd, double thickness_thresh, double sigma, int size_thresh) 
{
  MRI    *mri_lh, *mri_rh, *mri_lh_diff, *mri_rh_diff ;
  int    niter, vno, s ;
  MRI_SEGMENTATION *mriseg ;

  niter = SIGMA_TO_SURFACE_SMOOTH_STEPS(sigma) ;

  // do LH
  mri_lh = MRIclone(fcd->lh_thickness_on_lh, NULL) ;
  mri_rh = MRIclone(fcd->lh_thickness_on_lh, NULL) ;

  MRISwriteFrameToValues(fcd->mris_lh, fcd->lh_thickness_on_lh, 0) ;
  MRISaverageVals(fcd->mris_lh, niter) ;
  MRISreadFrameFromValues(fcd->mris_lh, mri_lh, 0) ;

  MRISwriteFrameToValues(fcd->mris_lh, fcd->rh_thickness_on_lh, 0) ;
  MRISaverageVals(fcd->mris_lh, niter) ;
  MRISreadFrameFromValues(fcd->mris_lh, mri_rh, 0) ;
  mri_lh_diff = MRIsubtract(mri_lh, mri_rh, NULL) ;  // lh minus rh on lh
  MRIfree(&mri_lh); MRIfree(&mri_rh) ;

  // do RH
  mri_lh = MRIclone(fcd->lh_thickness_on_rh, NULL) ;
  mri_rh = MRIclone(fcd->lh_thickness_on_rh, NULL) ;

  MRISwriteFrameToValues(fcd->mris_rh, fcd->lh_thickness_on_rh, 0) ;
  MRISaverageVals(fcd->mris_rh, niter) ;
  MRISreadFrameFromValues(fcd->mris_rh, mri_lh, 0) ;

  MRISwriteFrameToValues(fcd->mris_rh, fcd->rh_thickness_on_rh, 0) ;
  MRISaverageVals(fcd->mris_rh, niter) ;
  MRISreadFrameFromValues(fcd->mris_rh, mri_rh, 0) ;
  mri_rh_diff = MRIsubtract(mri_lh, mri_rh, NULL) ;  // lh minus rh on rh
  MRIfree(&mri_lh); MRIfree(&mri_rh) ;

  MRIclear(fcd->mri_thickness_increase) ;
  MRIclear(fcd->mri_thickness_decrease) ;
#if 1
#ifdef HAVE_OPENMP
#pragma omp parallel for shared(fcd, mri_lh_diff, Gdiag_no, thickness_thresh) schedule(static,1)
#endif
#endif
  for (vno = 0 ; vno < fcd->mris_lh->nvertices ; vno++)
  {
    double xs, ys, zs, d, xv, yv, zv  ;
    float  val, val2;
    int    label, xvi, yvi, zvi ;
    VERTEX *v ;

    v = &fcd->mris_lh->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    val = MRIgetVoxVal(mri_lh_diff, vno, 0, 0, 0) ;
    if (fabs(val) < thickness_thresh)
      continue ;

    for (d = 0 ; d < 4 ; d += 0.1)
    {
      xs = v->x+d*v->nx ; 
      ys = v->y+d*v->ny ; 
      zs = v->z+d*v->nz ; 
      MRISsurfaceRASToVoxel(fcd->mris_lh, fcd->mri_thickness_increase, xs, ys, zs, &xv, &yv, &zv) ;
      xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
      label = MRIgetVoxVal(fcd->mri_aseg, xvi, yvi, zvi, 0) ;
      if (IS_WM(label) == 0)
	break ;
    }

    if (val >= 0)
    {
      val2 = MRIgetVoxVal(fcd->mri_thickness_increase, xvi, yvi, zvi, 0) ;
      if (val > val2)  // another thread already populated this voxel
	MRIsetVoxVal(fcd->mri_thickness_increase, xvi, yvi, zvi, 0, val) ;
    }
    else 
    {
      val2 = MRIgetVoxVal(fcd->mri_thickness_decrease, xvi, yvi, zvi, 0) ;
      if (val < val2)  // another thread already populated this voxel
	MRIsetVoxVal(fcd->mri_thickness_decrease, xvi, yvi, zvi, 0, val) ;
    }
  }

#if 1
#ifdef HAVE_OPENMP
#pragma omp parallel for shared(fcd, mri_rh_diff, Gdiag_no, thickness_thresh) schedule(static,1)
#endif
#endif
  for (vno = 0 ; vno < fcd->mris_rh->nvertices ; vno++)
  {
    double xv, yv, zv, xs, ys, zs, d  ;
    float  val, val2;
    int   label, xvi, yvi, zvi ;
    VERTEX *v ;

    v = &fcd->mris_rh->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    val = MRIgetVoxVal(mri_rh_diff, vno, 0, 0, 0) ;
    if (fabs(val) < thickness_thresh)
      continue ;

    for (d = 0 ; d < 4 ; d += 0.1)
    {
      xs = v->x+d*v->nx ; 
      ys = v->y+d*v->ny ; 
      zs = v->z+d*v->nz ; 
      MRISsurfaceRASToVoxel(fcd->mris_lh, fcd->mri_thickness_increase, xs, ys, zs, &xv, &yv, &zv) ;
      xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
      label = MRIgetVoxVal(fcd->mri_aseg, xvi, yvi, zvi, 0) ;
      if (IS_WM(label) == 0)
	break ;
    }
    if (val >= 0)
    {
      val2 = MRIgetVoxVal(fcd->mri_thickness_increase, xvi, yvi, zvi, 0) ;
      if (val > val2)
	MRIsetVoxVal(fcd->mri_thickness_increase, xvi, yvi, zvi, 0, val) ;
    }
    else 
    {
      val2 = MRIgetVoxVal(fcd->mri_thickness_decrease, xvi, yvi, zvi, 0) ;
      if (val < val2)
	MRIsetVoxVal(fcd->mri_thickness_decrease, xvi, yvi, zvi, 0, val) ;
    }
  }

  mriseg = MRIsegment(fcd->mri_thickness_increase, thickness_thresh, 1e10) ;
  MRIremoveSmallSegments(mriseg, size_thresh) ;
  printf("segmenting volume at threshold %2.1f yields %d segments\n", thickness_thresh, mriseg->nsegments) ;
  fflush(stdout) ;

  fcd->nlabels = mriseg->nsegments ;
  for (s = 0 ; s < mriseg->nsegments ; s++)
  {
    int label ;

    fcd->labels[s] = MRIsegmentToLabel(mriseg, fcd->mri_thickness_increase, s) ;
    label = most_frequent_label(fcd->mri_aseg, &mriseg->segments[s]) ;
    sprintf(fcd->label_names[s], cma_label_to_name(label)) ;
    printf("%s\n", fcd->label_names[s]) ;
    fflush(stdout) ;
  }
  MRIfree(&mri_lh_diff) ; MRIfree(&mri_rh_diff) ;
  MRIsegmentFree(&mriseg) ;

  return(fcd->nlabels) ;
}


