/**
 * @file  mris_interpolate_warp.c
 * @brief interpolate a surface warp into the volume
 *
 * take two surfaces and interpret the spatial difference in their vertex locations as a warp field, 
 * then interpolate that into a volume warp.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/03/18 13:10:50 $
 *    $Revision: 1.2 $
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "mrisurf.h"
#include "gcamorph.h"
#include "mri.h"
#include "mrinorm.h"
#include "mri_circulars.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

static char *like_vol_name ;

char *Progname ;
static void usage_exit(int code) ;
static int niter = 50 ;

static int
write_surface_warp_into_volume(MRI_SURFACE *mris, MRI *mri, int niter) 
{
  int    vno, xvi, yvi, zvi, frame ;
  VERTEX *v ;
  double  dx, dy, dz, xv, yv, zv, xv1, yv1, zv1 ;
  MRI     *mri_weights, *mri_ctrl, *mri_frame ;
  float    wt ;

  mri_weights = MRIallocSequence(mri->width, mri->height, mri->height, MRI_FLOAT, 1) ;
  mri_ctrl = MRIallocSequence(mri->width, mri->height, mri->height, MRI_UCHAR, 1) ;
  MRIcopyHeader(mri, mri_weights) ;
  MRIcopyHeader(mri, mri_ctrl) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    MRISsurfaceRASToVoxelCached(mris, mri, v->origx, v->origy, v->origz, &xv, &yv, &zv) ;
    MRISsurfaceRASToVoxelCached(mris, mri, v->x, v->y, v->z, &xv1, &yv1, &zv1) ;
    dx = xv1-xv ; dy = yv1-yv ; dz = zv1-zv ;
    //    dx = xv1 ; dy = yv1 ; dz = zv1 ;
    xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
    if (xvi < 0 || xvi >= mri->width || yv < 0 || yv >= mri->height || zv < 0 || zv >= mri->depth)
      continue ;
    MRIinterpolateIntoVolumeFrame(mri, xv, yv, zv, 0, dx) ;
    MRIinterpolateIntoVolumeFrame(mri, xv, yv, zv, 1, dy) ;
    MRIinterpolateIntoVolumeFrame(mri, xv, yv, zv, 2, dz) ;
    MRIinterpolateIntoVolume(mri_weights, xv, yv, zv, 1.0) ;
  }

  for (xvi = 0 ; xvi < mri->width ; xvi++)
    for (yvi = 0 ; yvi < mri->height ; yvi++)
    {
      MRIsetVoxVal(mri_ctrl, xvi, yvi, 0, 0, CONTROL_MARKED) ;
      MRIsetVoxVal(mri, xvi, yvi, 0, 0, 0) ;
      MRIsetVoxVal(mri, xvi, yvi, 0, 1, 0) ;
      MRIsetVoxVal(mri, xvi, yvi, 0, 2, 0) ;

      MRIsetVoxVal(mri_ctrl, xvi, yvi, mri->depth-1, 0, CONTROL_MARKED) ;
      MRIsetVoxVal(mri, xvi, yvi, mri->depth-1, 0, 0) ;
      MRIsetVoxVal(mri, xvi, yvi, mri->depth-1, 1, 0) ;
      MRIsetVoxVal(mri, xvi, yvi, mri->depth-1, 2, 0) ;
    }

  for (xvi = 0 ; xvi < mri->width ; xvi++)
    for (zvi = 0 ; zvi < mri->depth ; zvi++)
    {
      MRIsetVoxVal(mri_ctrl, xvi, 0, zvi, 0, CONTROL_MARKED) ;
      MRIsetVoxVal(mri, xvi, 0, zvi, 0, 0) ;
      MRIsetVoxVal(mri, xvi, 0, zvi, 1, 0) ;
      MRIsetVoxVal(mri, xvi, 0, zvi, 2, 0) ;

      MRIsetVoxVal(mri_ctrl, xvi, mri->height-1, zvi, 0, CONTROL_MARKED) ;
      MRIsetVoxVal(mri, xvi, mri->height-1, zvi, 0, 0) ;
      MRIsetVoxVal(mri, xvi, mri->height-1, zvi, 1, 0) ;
      MRIsetVoxVal(mri, xvi, mri->height-1, zvi, 2, 0) ;
    }

  for (yvi = 0 ; yvi < mri->width ; yvi++)
    for (zvi = 0 ; zvi < mri->depth ; zvi++)
    {
      MRIsetVoxVal(mri_ctrl, 0, yvi, zvi, 0, CONTROL_MARKED) ;
      MRIsetVoxVal(mri, 0, yvi, zvi, 0, 0) ;
      MRIsetVoxVal(mri, 0, yvi, zvi, 1, 0) ;
      MRIsetVoxVal(mri, 0, yvi, zvi, 2, 0) ;

      MRIsetVoxVal(mri_ctrl, mri->width-1, yvi, zvi, 0, CONTROL_MARKED) ;
      MRIsetVoxVal(mri, mri->width-1, yvi, zvi, 0, 0) ;
      MRIsetVoxVal(mri, mri->width-1, yvi, zvi, 1, 0) ;
      MRIsetVoxVal(mri, mri->width-1, yvi, zvi, 2, 0) ;

    }
  for (xvi = 0 ; xvi < mri->width ; xvi++)
    for (yvi = 0 ; yvi < mri->height ; yvi++)
      for (zvi = 0 ; zvi < mri->depth ; zvi++)
      {
        if (xvi == Gx && yvi == Gy && zvi == Gz)
          DiagBreak() ;
        wt = MRIgetVoxVal(mri_weights, xvi, yvi, zvi, 0) ;
        dx = MRIgetVoxVal(mri, xvi, yvi, zvi, 0) ;
        dy = MRIgetVoxVal(mri, xvi, yvi, zvi, 1) ;
        dz = MRIgetVoxVal(mri, xvi, yvi, zvi, 2) ;
        if (FZERO(wt))
          continue ;
        dx /= wt ; dy /= wt ; dz /= wt ;
        MRIsetVoxVal(mri, xvi, yvi, zvi, 0, dx) ;
        MRIsetVoxVal(mri, xvi, yvi, zvi, 1, dy) ;
        MRIsetVoxVal(mri, xvi, yvi, zvi, 2, dz) ;
        MRIsetVoxVal(mri_ctrl, xvi, yvi, zvi, 0, CONTROL_MARKED) ;  // it is a control point
      }

  MRIwrite(mri, "warp0.mgz") ;
  MRIwrite(mri_ctrl, "ctrl.mgz") ;
  for (frame = 0 ; frame < 3 ; frame++)
  {
    printf("interpolating frame %d\n", frame+1) ;
    mri_frame = MRIcopyFrame(mri, NULL, frame, 0) ;
    MRIbuildVoronoiDiagram(mri_frame, mri_ctrl, mri_frame) ;
    MRIsoapBubble(mri_frame, mri_ctrl, mri_frame, niter) ;
    
    {
      int x, y, z ;
      float val ;
      for (x  = 0 ; x < mri_frame->width ; x++)
        for (y  = 0 ; y < mri_frame->height ; y++)
          for (z  = 0 ; z < mri_frame->depth ; z++)
          {
            val = MRIgetVoxVal(mri_frame, x, y, z, 0) ;
            switch (frame)
            {
            default:
            case 0: val += x ; break ;
            case 1: val += y ; break ;
            case 2: val += z ; break ;
            }
            
            MRIsetVoxVal(mri_frame, x, y, z, 0, val) ;
          }
    }
    MRIcopyFrame(mri_frame, mri, 0, frame) ;
    MRIfree(&mri_frame) ;
  }
  MRIfree(&mri_weights) ; MRIfree(&mri_ctrl) ;
  return(NO_ERROR) ;
}

int
main(int argc, char *argv[]) {
  char        **av, *out_name ;
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  MRI_SURFACE  *mris ;
  GCA_MORPH    *gcam ;
  MRI          *mri = NULL ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_interpolate_warp.c,v 1.2 2011/03/18 13:10:50 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;

  mris = MRISread(argv[2]) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read source surface %s\n", Progname,argv[2]) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  if (MRISreadVertexPositions(mris, argv[1]) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read target surface %s\n", Progname,argv[1]) ;
  if (like_vol_name == NULL)
  {
    mri = MRIallocSequence(mris->vg.width, mris->vg.height, mris->vg.depth, MRI_FLOAT, 3) ;
    MRIcopyVolGeomToMRI(mri, &mris->vg) ;
  }
  else
  {
    MRI *mri_tmp ;
    mri_tmp = MRIread(like_vol_name) ;
    if (mri_tmp == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not like volume %s\n", like_vol_name) ;
    mri = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, 3) ;
    MRIcopyHeader(mri_tmp, mri) ;
    MRIfree(&mri_tmp) ;
  }
  write_surface_warp_into_volume(mris, mri, niter) ;

  MRIwrite(mri, "warp.mgz") ;
  gcam = GCAMalloc(mri->width, mri->height, mri->depth) ;
  GCAMinitVolGeom(gcam, mri, mri) ;
  GCAMreadWarpFromMRI(gcam, mri) ;
  //  GCAsetVolGeom(gca, &gcam->atlas);
  out_name = argv[3] ;
  GCAMwrite(gcam, out_name) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "warp field calculation took  %d minutes and %d seconds.\n", minutes, seconds) ;
  exit(0) ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  switch (toupper(*option)) 
  {
  case 'I':
    niter = atoi(argv[2]) ;
    printf("using %d soap bubble iterations\n", niter) ;
    nargs = 1 ;
    break ;
  case 'L':
    like_vol_name = argv[2] ;
    printf("using volume %s to specify geometry of volumetric warp field\n", like_vol_name) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s [options] <start surface> <end surface> <warp field>.m3z\n", Progname) ;
  exit(code) ;
}





