/**
 * @brief interpolate a surface warp into the volume
 *
 * take two surfaces and interpret the spatial difference in their
 * vertex locations as a warp field, then interpolate that into a volume warp.
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

const char *Progname ;
static void usage_exit(int code) ;
static int niter = 500 ;
static MRI *mri_in = NULL ;
static char *out_fname ;

static int no_write = 0 ;

static int pad = 20 ;

static int
write_surface_warp_into_volume(MRI_SURFACE *mris, MRI *mri, int niter)
{
  int    vno, xvi, yvi, zvi, frame ;
  VERTEX *v ;
  double  dx, dy, dz, xv, yv, zv, xv1, yv1, zv1 ;
  MRI     *mri_weights, *mri_ctrl, *mri_frame ;
  float    wt ;

  mri_weights = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, 1) ;
  mri_ctrl = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_UCHAR, 1) ;
  MRIcopyHeader(mri, mri_weights) ; MRIcopyHeader(mri, mri_ctrl) ;
  //   build a 3 frame volume with the voxel-coords warp (dx, dy, dz) in frames 0, 1 and 2 respectively
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    MRISsurfaceRASToVoxelCached(mris, mri, v->origx, v->origy, v->origz, &xv, &yv, &zv) ;
    MRISsurfaceRASToVoxelCached(mris, mri, v->x, v->y, v->z, &xv1, &yv1, &zv1) ;
    dx = xv1-xv ; dy = yv1-yv ; dz = zv1-zv ;
    xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
    if (vno == Gdiag_no)
    {
      printf("surface vertex %d: inflated (%2.0f, %2.0f, %2.0f), orig (%2.0f, %2.0f, %2.0f), "
	     "dx=(%2.0f, %2.0f, %2.0f)\n", vno, xv1, yv1, zv1, xv, yv, zv, dx, dy, dz) ;
      DiagBreak() ;
    }
    if (xvi < 0 || xvi >= mri->width || yv < 0 || yv >= mri->height || zv < 0 || zv >= mri->depth)
    {
      continue ;
    }
    MRIinterpolateIntoVolumeFrame(mri, xv, yv, zv, 0, dx) ;
    MRIinterpolateIntoVolumeFrame(mri, xv, yv, zv, 1, dy) ;
    MRIinterpolateIntoVolumeFrame(mri, xv, yv, zv, 2, dz) ;
    MRIinterpolateIntoVolume(mri_weights, xv, yv, zv, 1.0) ;
  }

#if 0
  // set boundary conditions in the edge planes to be 0 warping
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
#endif

  // normalize the warp field using a weighted average of all vertices that map to every voxel
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

  if (Gdiag & DIAG_WRITE)
  {
    MRIwrite(mri, "warp0.mgz") ;
    MRIwrite(mri_ctrl, "ctrl.mgz") ;
  }

  for (frame = 0 ; frame < mri->nframes ; frame++)
  {
    printf("interpolating frame %d\n", frame+1) ;
    mri_frame = MRIcopyFrame(mri, NULL, frame, 0) ;
    MRIbuildVoronoiDiagram(mri_frame, mri_ctrl, mri_frame) ;
    MRIsoapBubble(mri_frame, mri_ctrl, mri_frame, niter, mri_frame->xsize*.05) ;
#if 0
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
            case 0:
              val += x ;
              break ;
            case 1:
              val += y ;
              break ;
            case 2:
              val += z ;
              break ;
            }

            MRIsetVoxVal(mri_frame, x, y, z, 0, val) ;
          }
    }
#endif
    MRIcopyFrame(mri_frame, mri, 0, frame) ;
    MRIfree(&mri_frame) ;
  }
  MRIfree(&mri_weights) ;
  MRIfree(&mri_ctrl) ;
  return(NO_ERROR) ;
}

static MRI *
expand_mri_to_fit_surface(MRI_SURFACE *mris, MRI *mri)
{
  int      vno ;
  VERTEX   *v ;
  double   xv, yv, zv, xmin, ymin, zmin, xmax, ymax, zmax, x0, y0, z0, x1, y1, z1, 
    dx_dx, dx_dy, dx_dz,
    dy_dx, dy_dy, dy_dz,
    dz_dx, dz_dy, dz_dz ;
  MRI        *mri_dst ;

  
  MRISsurfaceRASToVoxelCached(mris, mri, 0, 0, 0, &x0, &y0, &z0) ;
  MRISsurfaceRASToVoxelCached(mris, mri, 1, 0, 0, &x1, &y1, &z1) ;
  dx_dx = x1-x0 ; dx_dy = y1-y0 ; dx_dz = z1-z0 ;
  MRISsurfaceRASToVoxelCached(mris, mri, 0, 1, 0, &x1, &y1, &z1) ;
  dy_dx = x1-x0 ; dy_dy = y1-y0 ; dy_dz = z1-z0 ;
  MRISsurfaceRASToVoxelCached(mris, mri, 0, 0, 1, &x1, &y1, &z1) ;
  dz_dx = x1-x0 ; dz_dy = y1-y0 ; dz_dz = z1-z0 ;
  xmax=ymax=zmax=xmin=ymin=zmin=0 ;  // for compiler warnings
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    MRISsurfaceRASToVoxelCached(mris, mri, v->origx, v->origy, v->origz, &xv, &yv, &zv) ;
    if (vno == 0)
    {
      xmin = xmax = xv ; ymin = ymax = yv ; zmin = zmax = zv ;
    }
    if (xv < xmin)
      xmin = xv ;
    if (yv < ymin)
      ymin = yv ;
    if (zv < zmin)
      zmin = zv ;
    if (xv > xmax)
      xmax = xv ;
    if (yv > ymax)
      ymax = yv ;
    if (zv > zmax)
      zmax = zv ;
#if 0
    v->origy += 30 ;
    v->y += 30 ;
#endif
  }
  if (xmin > 0)
    xmin = 0 ;
  if (ymin > 0)
    ymin = 0 ;
  if (zmin > 0)
    zmin = 0 ;
  if (xmax < mri->width-1)
    xmax = mri->width-1 ;
  if (ymax < mri->height-1)
    ymax = mri->height-1 ;
  if (zmax < mri->depth-1)
    zmax = mri->depth-1 ;

  mri_dst = MRIextractRegionAndPad(mri, NULL, NULL, pad) ;
  return(mri_dst) ;
}

int
main(int argc, char *argv[])
{
  char        **av, *out_name ;
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  Timer start ;
  MRI_SURFACE  *mris ;
  GCA_MORPH    *gcam ;
  MRI          *mri = NULL ;

  nargs = handleVersionOption(argc, argv, "mris_interpolate_warp");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
  {
    usage_exit(1) ;
  }


  /*
    note that a "forward" morph means a retraction, so we reverse the order of the argvs here.
    This means that for every voxel in the inflated image we have a vector that points to where in
    the original image it came from, and *NOT* the reverse.
  */
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
    {
      ErrorExit(ERROR_NOFILE, "%s: could not like volume %s\n", like_vol_name) ;
    }
    mri = MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth, MRI_FLOAT, 3) ;
    MRIcopyHeader(mri_tmp, mri) ;
    MRIfree(&mri_tmp) ;
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    double xv, yv, zv ;
    VERTEX *v = &mris->vertices[0] ;
    MRISsurfaceRASToVoxel(mris, mri, v->x, v->y, v->z, &xv, &yv, &zv) ;
    printf("v 0: sras (%f, %f, %f) --> vox (%f, %f, %f)\n", v->x,v->y,v->z,xv,yv,zv);
    MRISsurfaceRASToVoxelCached(mris, mri, v->x, v->y, v->z, &xv, &yv, &zv) ;
    printf("v 0: sras (%f, %f, %f) --> vox (%f, %f, %f)\n", v->x,v->y,v->z,xv,yv,zv);
    DiagBreak() ;
  }
  {
    MRI *mri_tmp ;
    mri_tmp = expand_mri_to_fit_surface(mris, mri) ;
    MRIfree(&mri) ; mri = mri_tmp ;
  }
  write_surface_warp_into_volume(mris, mri, niter) ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri, "warp.mgz") ;
  gcam = GCAMalloc(mri->width, mri->height, mri->depth) ;
  GCAMinitVolGeom(gcam, mri, mri) ;
  GCAMremoveSingularitiesAndReadWarpFromMRI(gcam, mri) ;
//  GCAMreadWarpFromMRI(gcam, mri) ;
  //  GCAsetVolGeom(gca, &gcam->atlas);
#if 0
  gcam->gca = gcaAllocMax(1, 1, 1,
			  mri->width, mri->height,
			  mri->depth,
			  0, 0) ;
 GCAMinit(gcam, mri, NULL, NULL, 0) ;
#endif
#if 0
  GCAMinvert(gcam, mri) ;
  GCAMwriteInverseWarpToMRI(gcam, mri) ;
  GCAMremoveSingularitiesAndReadWarpFromMRI(gcam, mri) ;  // should be inverse now
#endif
  if (mri_in)
  {
    MRI *mri_warped, *mri_tmp ;
    printf("applying warp to %s and writing to %s\n", mri_in->fname, out_fname) ;
    mri_tmp = MRIextractRegionAndPad(mri_in, NULL, NULL, pad) ; MRIfree(&mri_in) ; mri_in = mri_tmp ;
    mri_warped = GCAMmorphToAtlas(mri_in, gcam, NULL, -1, SAMPLE_TRILINEAR) ;
    MRIwrite(mri_warped, out_fname) ;
    if (Gdiag_no >= 0)
    {
      double  xi, yi, zi, xo, yo, zo, val;
      int     xp, yp, zp ;
      GCA_MORPH_NODE *gcamn ;

      VERTEX *v = &mris->vertices[Gdiag_no] ;
      MRISsurfaceRASToVoxelCached(mris, mri, v->origx, v->origy, v->origz, &xi, &yi, &zi) ;
      MRISsurfaceRASToVoxelCached(mris, mri, v->x, v->y, v->z, &xo, &yo, &zo) ;
      printf("surface vertex %d: inflated (%2.0f, %2.0f, %2.0f), orig (%2.0f, %2.0f, %2.0f)\n", Gdiag_no, xi, yi, zi, xo, yo, zo) ;
      MRIsampleVolume(mri_in, xo, yo, zo, &val) ;
      xp = nint(xi) ; yp = nint(yi) ; zp = nint(zi) ;
      gcamn = &gcam->nodes[xp][yp][zp] ;
      printf("warp = (%2.1f, %2.1f, %2.1f), orig (%2.1f %2.1f %2.1f) = %2.1f \n", 
	     gcamn->x, gcamn->y, gcamn->z,
	     gcamn->origx, gcamn->origy, gcamn->origz,val) ;
      DiagBreak() ;
    }
  }
  if (no_write == 0)
  {
    out_name = argv[3] ;
    GCAMwrite(gcam, out_name) ;
  }
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "warp field calculation took %d minutes and %d seconds.\n", minutes, seconds) ;
  exit(0) ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  switch (toupper(*option))
  {
  case 'N':
    no_write = 1 ;
    printf("not saving final warp (specify -a to save warped volume)\n") ;
    break ;
  case 'A':
    mri_in = MRIread(argv[2]) ;
    if (mri_in == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read volume %s", Progname, argv[2]) ;
    out_fname = argv[3] ;
    nargs = 2 ;
    break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
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
usage_exit(int code)
{
  printf("usage: %s [options] <start surface> <end surface> <warp field>.m3z\n", Progname) ;
  exit(code) ;
}





