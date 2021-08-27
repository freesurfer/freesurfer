/**
 * @brief update a conformed segmentation with hires data
 *
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

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

static MRI *MRIupdateSegmentation(MRI *mri_hires, MRI *mri_wm_in,  MRI *mri_wm_out, MRI *mri_mask) ;

static float mag_thresh = 3 ;
static float intensity_thresh = 4 ;
static float dist_thresh = .5 ;
static int nclose = 1 ;
static MRI *mri_mask = NULL ;

int
main(int argc, char *argv[]) {
  char         **av ;
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  Timer start ;
  MRI          *mri_wm, *mri_hires ;
  TRANSFORM    *transform ;
  MRI_SURFACE  *mris ;

  nargs = handleVersionOption(argc, argv, "mri_multiscale_segment");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
    usage_exit(1) ;

  mris = MRISread(argv[1]) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s\n",
              argv[1]) ;
  transform = TransformRead(argv[2]) ;
  if (transform == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s\n",
              argv[2]) ;
  mri_hires = MRIread(argv[3]) ;
  if (mri_hires == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not hires volume from %s\n",
              argv[3]) ;

  printf("filling interior at %2.3fmm resolution...\n", mri_hires->xsize) ;
  mri_wm = MRIclone(mri_hires, NULL) ;
  MRISfillInterior(mris, mri_hires->xsize, mri_wm) ;
  printf("updating segmentation...\n") ;

  MRIwrite(mri_wm, "wmb.mgz") ;
  MRIupdateSegmentation(mri_hires, mri_wm, mri_wm, mri_mask) ;
  
  printf("writing output segmentation to %s\n", argv[4]) ;
  if (nclose > 1)
  {
    int i ;
    MRI *mri_tmp = MRIcopy(mri_wm, NULL) ;
    for (i = 0; i < nclose ; i++)
    {
      MRIdilate(mri_wm, mri_tmp) ;
      MRIcopy(mri_tmp, mri_wm) ;
    }
    for (i = 0; i < nclose ; i++)
    {
      MRIerode(mri_wm, mri_tmp) ;
      MRIcopy(mri_tmp, mri_wm) ;
    }
    MRIfree(&mri_tmp) ;
  }
  else if (nclose == 1)
    MRIclose(mri_wm, mri_wm) ;
      
  MRIwrite(mri_wm, argv[4]) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "multiscale segmentation took %d minutes and %d seconds\n", 
          minutes, seconds) ;
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
  if (!stricmp(option, "debug_voxel"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "dthresh"))
  {
    printf("using distance threshold %2.2f (default = %2.2f)\n", atof(argv[2]), dist_thresh) ;
    dist_thresh = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "ithresh"))
  {
    printf("using intensity threshold %2.2f (default = %2.2f)\n", atof(argv[2]), intensity_thresh) ;
    intensity_thresh = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "mthresh"))
  {
    printf("using gradient magnitude threshold %2.2f (default = %2.2f)\n", atof(argv[2]), mag_thresh) ;
    mag_thresh = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "mask"))
  {
    nargs = 2 ;
    printf("reading %s as mask volume thresholded at %2.1f\n", argv[2], atof(argv[3])) ;
    mri_mask = MRIread(argv[2]) ;
    if (mri_mask == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load mask volume from %s", Progname, argv[2]) ;
    MRIthreshold(mri_mask, mri_mask, atof(argv[2])) ;
    MRIwrite(mri_mask, "m.mgz") ;
  }
  else switch (toupper(*option)) {
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
  printf("usage: %s [options] <inverse operator> <EEG/MEG data file>",
         Progname) ;
  printf(
    "\tf <f low> <f hi> - apply specified filter (not implemented yet)\n"
  );
  printf("\tn - noise-sensitivity normalize inverse (default=1)") ;
  exit(code) ;
}




float
MRIclosestOnVal(MRI *mri_src, MRI *mri_mask, int x, int y, int z, int whalf) 
{
  float     val, closest_val, val0 ; 
  int       xi, yi, zi, xk, yk, zk ;

  closest_val = -1 ;
  val0 = MRIgetVoxVal(mri_src, x, y, z, 0) ;
  for (xk = -1 ; xk <= 1 ; xk++)
  {
    xi = mri_src->xi[x+xk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri_src->yi[y+yk] ;
      for (zk = -1 ; zk <= 1 ; zk++)
      {
        zi = mri_src->zi[z+zk] ;
        if (MRIgetVoxVal(mri_mask, xi, yi, zi, 0) == 0)
          continue ;
        val = MRIgetVoxVal(mri_src, xi, yi, zi, 0) ;
        if (fabs(val-val0) < fabs(val-closest_val))
          closest_val = val ;
      }
    }
  }
        
  return(closest_val) ;
}

static MRI *
MRIupdateSegmentation(MRI *mri_hires, MRI *mri_wm_in,  MRI *mri_wm_out, MRI *mri_mask)
{
  MRI   *mri_means, *mri_stds, *mri_tmp, *mri_eroded  ;
  int   x, y, z, i, removed, added, dxi, dyi, dzi, xi, yi, zi ;
  float mean, std, val, dist, val0, theta, dx, dy, dz, norm ;
#if 1
  MRI   *mri_mag, *mri_grad, *mri_smooth ;

  mri_smooth = MRIgaussianSmooth(mri_hires, 0.5, 1, NULL) ;
  mri_mag = MRIcloneDifferentType(mri_hires, MRI_FLOAT) ;
  mri_grad = MRIsobel(mri_smooth, NULL, mri_mag) ;
  MRIwrite(mri_grad, "grad.mgz") ;
  MRIwrite(mri_mag, "mag.mgz") ;
#endif
  mri_means = MRImeanInMask(mri_hires, NULL, mri_wm_in, 9) ;
  MRIwrite(mri_means, "m.mgz") ;
  mri_stds = MRIstdInMask(mri_hires, NULL, mri_means, mri_wm_in, 9) ;
  MRIwrite(mri_stds, "s.mgz") ;

  mri_wm_out = MRIcopy(mri_wm_in, mri_wm_out) ;
  mri_tmp = MRIcopy(mri_wm_in, NULL) ;
  mri_eroded = MRIcopy(mri_wm_in, NULL) ;

  for (i = 0 ; i < 100 ; i++)
  {
    added = removed = 0 ;
    mean = MRImeanInLabel(mri_hires, mri_wm_out, 1) ;
    std = MRImeanInLabel(mri_stds, mri_wm_out, 1) ;
    MRIcopy(mri_wm_out, mri_tmp) ;
    printf("iter %d: wm mean = %2.1f +- %2.1f\n", i, mean, std) ;

#if 0
    // first remove voxels
    for (x = 0 ; x < mri_means->width ; x++)
      for (y = 0 ; y < mri_means->height ; y++)
        for (z = 0 ; z < mri_means->depth ; z++)
        {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          val = MRIgetVoxVal(mri_hires, x, y,z, 0) ;
          dist = (val - mean) / std ;
          if (MRIgetVoxVal(mri_wm_out, x,y, z,0) == 0)
            continue ;
          //        mean = MRIgetVoxVal(mri_means, x, y, z, 0) ;
          //        std = MRIgetVoxVal(mri_stds, x, y, z, 0) ;
          if (fabs(dist) > 2)
          {
            removed++ ;
            MRIsetVoxVal(mri_tmp, x, y, z, 0, 0) ;
          }
        }
    printf("%d voxels removed from wm\n", removed) ;
    MRIcopy(mri_tmp, mri_wm_out) ;
    mean = MRImeanInLabel(mri_hires, mri_wm_out, 1) ;
    std = MRImeanInLabel(mri_stds, mri_wm_out, 1) ;
    printf("        wm mean = %2.1f +- %2.1f\n", i, mean, std) ;
    {
      char fname[STRLEN] ;
      sprintf(fname, "wr%d.mgz", i+1) ;
      printf("writing segmentation to %s\n", fname) ;
      MRIwrite(mri_wm_out, fname) ;
    }
#endif
    MRIerode(mri_tmp, mri_eroded) ;

    for (x = 0 ; x < mri_means->width ; x++)
      for (y = 0 ; y < mri_means->height ; y++)
        for (z = 0 ; z < mri_means->depth ; z++)
        {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          if (mri_mask && MRIgetVoxVal(mri_mask, x, y, z, 0) == 0)
            continue ;

          // find voxels in the mask that are adjacent to ones that are not
          dx = MRIgetVoxVal(mri_grad, x, y, z, 0) ;
          dy = MRIgetVoxVal(mri_grad, x, y, z, 1) ;
          dz = MRIgetVoxVal(mri_grad, x, y, z, 2) ;
          norm = MRIgetVoxVal(mri_mag, x, y, z, 0) ;
          if (norm > mag_thresh)
            continue ;
          if (FZERO(norm) == 0)
          {
            dx /= norm ; dy /= norm ; dz /= norm ;
          }

          if (MRIgetVoxVal(mri_wm_out, x, y, z, 0) > 0 &&
              MRIgetVoxVal(mri_eroded, x, y, z, 0) == 0)
          {
            val0 = MRIgetVoxVal(mri_hires, x, y, z, 0) ;
            for (dxi = -1 ; dxi <= 1 ; dxi++)
            {
              xi = mri_wm_out->xi[x+dxi] ;
              for (dyi = -1 ; dyi <= 1 ; dyi++)
              {
                yi = mri_wm_out->yi[y+dyi] ;
                for (dzi = -1 ; dzi <= 1 ; dzi++)
                {
                  zi = mri_wm_out->zi[z+dzi] ;
                  if (abs(dxi) + abs(dyi) + abs(dzi) != 1)
                    continue ;   // only consider 6-connected nbrs
                  if (xi == Gx && yi == Gy && zi == Gz)
                    DiagBreak() ;
                  if (MRIgetVoxVal(mri_wm_out,xi,yi,zi,0) > 0)
                    continue ; // look for voxels that aren't on
                  theta = acos(dx*dxi + dy*dyi + dz*dzi) ;
                  if (fabs(theta) > M_PI/6)  // not in grad dir
                  {
                    val = MRIgetVoxVal(mri_hires, xi, yi, zi, 0) ;
                    dist = (val - mean) / std ;
                    if ((fabs(val0 - val) <= intensity_thresh) && (dist < dist_thresh))
                    {
                      added++ ;
                      MRIsetVoxVal(mri_tmp, xi, yi, zi, 0, 1) ;
                    }
                  }
                }
              }
            }
          }
        }
    printf("%d voxels added to wm\n", added) ;
    MRIcopy(mri_tmp, mri_wm_out) ;
    {
      char fname[STRLEN] ;
      sprintf(fname, "w%d.mgz", i+1) ;
      printf("writing segmentation to %s\n", fname) ;
      MRIwrite(mri_wm_out, fname) ;
    }
    if (removed == 0 && added == 0)
      break ;
  }

  MRIfree(&mri_means) ; MRIfree(&mri_stds) ; 
  MRIfree(&mri_tmp) ; MRIfree(&mri_eroded) ;
#if 1
  MRIfree(&mri_grad) ; MRIfree(&mri_mag) ; 
#endif
  return(mri_wm_out) ;
}

