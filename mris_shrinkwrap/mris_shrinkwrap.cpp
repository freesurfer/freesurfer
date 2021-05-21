/**
 * @brief program for shrinkwrapping BEM surfaces onto a segmentation volume
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
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"

#include "mri.h"
#include "mrimorph.h"
#include "mrinorm.h"
#include "mrisegment.h"
#include "mrisurf.h"
#include "mrisurf_project.h"
#include "icosahedron.h"

#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"
#include "version.h"
#include "cma.h"


int main(int argc, char *argv[]) ;

#define INNER_SKULL_OUTER_SKULL_SEPARATION 4
#define BORDER_VAL           128
#define OUTSIDE_BORDER_STEP  16
#define TARGET_VAL           (BORDER_VAL-OUTSIDE_BORDER_STEP/2)

static MRI *pad_volume(MRI_SURFACE *mris, MRI *mri_src, MRI *mri_dst) ;

static int pad = 40 ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static MRI *create_brain_volume(MRI *mri_labeled, 
                                MRI *mri_brain, 
                                int target_label) ;
static MRI *create_skull_volume(MRI *mri_labeled, 
                                MRI *mri_brain) ;
static MRI *create_skin_volume(MRI *mri_labeled, 
                               MRI *mri_brain) ;
static int initialize_surface_position(MRI_SURFACE *mris, 
                                       MRI *mri_masked, 
                                       int outside) ;
static MRI *create_distance_map(MRI *mri_masked, 
                                MRI *mri_distance, 
                                int border_val, 
                                int outside_border_step) ;
static MRI *remove_small_segments(MRI  *mri_src, MRI *mri_dst) ;

const char *Progname ;

static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE    1.0
static float base_dt_scale = BASE_DT_SCALE ;
static int target_label = -1 ;
static const char *suffix = "" ;
static const char *output_suffix = "" ;
static double l_tsmooth = 0.0 ;
static double l_surf_repulse = 5.0 ;
static int smooth = 5 ;
static int nbrs = 2 ;
static int ic = 5 ;
static int inner_skull_only = 0;
static int embed = 0 ;

static float threshold = 0.0 ;

int
main(int argc, char *argv[])
{
  char          **av, fname[STRLEN], *vol_name, *output_dir, *mdir ;
  int           ac, nargs, msec ;
  MRI_SURFACE   *mris ;
  MRI           *mri_labeled, *mri_masked, *mri_masked_smooth, *mri_tmp;
  MRI           *mri_kernel, *mri_dist ;
  Timer then ;
  double        l_spring ;

  nargs = handleVersionOption(argc, argv, "mris_shrinkwrap");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  // memset(&parms, 0, sizeof(parms)) ; Have proper constructor now
  parms.projection = NO_PROJECTION ;
  parms.tol = 0.05 ;
  parms.dt = 0.5f ;
  parms.base_dt = BASE_DT_SCALE*parms.dt ;

  parms.l_intensity = 0.1 ;
  parms.l_spring = 1.0f ;
  parms.l_shrinkwrap = 0.05 ;

  parms.niterations = 0 ;
  parms.write_iterations = 0 /*WRITE_ITERATIONS */;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = 0.0 /*0.8*/ ;
  parms.dt_increase = 1.0 /* DT_INCREASE */;
  parms.dt_decrease = 0.50 /* DT_DECREASE*/ ;
  parms.error_ratio = 50.0 /*ERROR_RATIO */;
  /*  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;*/
  parms.l_surf_repulse = 0.0 ;
  parms.l_repulse = 1 ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  mdir = getenv("FREESURFER_HOME") ;
  if (!mdir)
    ErrorExit(ERROR_BADPARM, "FREESURFER_HOME not defined in environment") ;

  if (argc < 3)
    usage_exit() ;

  /* set default parameters for white and gray matter surfaces */
  parms.niterations = 1000 ;
  if (parms.momentum < 0.0)
    parms.momentum = 0.0 /*0.75*/ ;

  then.reset() ;
  vol_name = argv[1] ;
  output_dir = argv[2] ;
  fprintf(stderr, "reading volume %s...\n", vol_name) ;
  mri_labeled = MRIread(vol_name) ;
  if (!mri_labeled)
    ErrorExit(ERROR_NOFILE, 
              "%s: could not read input volume %s", 
              Progname, vol_name) ;
  if (embed)
  {
    MRI       *mri_tmp ;

    mri_tmp = MRIextractRegionAndPad(mri_labeled, NULL, NULL, pad) ;
    MRIfree(&mri_labeled) ; mri_labeled = mri_tmp ;
  }

  if (threshold > 0)
    MRIbinarize(mri_labeled, mri_labeled, threshold, 0, target_label) ;

  ////////////////////////////// we can handle only conformed volumes
//  setMRIforSurface(mri_labeled);
  sprintf(fname, "%s/lib/bem/ic%d.tri", mdir, ic) ;
  fprintf(stderr, "reading %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, 
              "%s: could not read icosahedron %s", 
              Progname, fname) ;

  /* 
   * first create brain volume 
   */
  fprintf(stderr, "creating brain volume...\n") ;
  mri_masked = create_brain_volume(mri_labeled, NULL, target_label) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_masked, "brain.mgz") ;
  parms.sigma = 8 ;

  initialize_surface_position(mris, mri_masked, 1) ;
  mri_tmp = pad_volume(mris, mri_masked, NULL) ;
  MRIfree(&mri_masked) ; mri_masked = mri_tmp ;
  MRIScopyVolGeomFromMRI(mris, mri_masked) ;
  if (target_label < 0)
  {
    mri_dist = 
      create_distance_map(mri_masked, NULL, BORDER_VAL, OUTSIDE_BORDER_STEP) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_dist, "d.mgz") ;
    MRISsetVals(mris,parms.sigma) ;
    MRIScopyValToVal2(mris) ;
    MRISsetVals(mris, TARGET_VAL) ;

    mri_kernel = MRIgaussian1d(parms.sigma, 0) ;
    mri_masked_smooth = MRIconvolveGaussian(mri_dist, NULL, mri_kernel) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_masked_smooth, "brain_smooth.mgh") ;
    sprintf(parms.base_name, "%s_inner_skull%s%s", 
            vol_name, output_suffix, suffix) ;
    parms.l_intensity = 0 ;   // use shrinkwrap term
  }
  else    // shrink wrapping to a labeled or threshold volume - build distance transform
  {
    mri_dist = 
      MRIdistanceTransform(mri_masked, 
                           NULL, BORDER_VAL, 10000, DTRANS_MODE_SIGNED,NULL) ;
    parms.l_shrinkwrap = 0 ;  // use intensity term
    parms.l_spring_norm = parms.l_spring ;
    parms.l_spring = 0 ;
    sprintf(parms.base_name, "%s%s%s", 
            FileNameOnly(output_dir, fname), output_suffix, suffix) ;
    MRISsetVals(mris,parms.sigma) ;
    MRIScopyValToVal2(mris) ;
    MRISsetVals(mris, 0) ;  // 0 distance is the target

    mri_kernel = MRIgaussian1d(parms.sigma, 0) ;
    mri_masked_smooth = MRIconvolveGaussian(mri_dist, NULL, mri_kernel) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_masked_smooth, "brain_smooth.mgh") ;
  }
  parms.mri_brain = mri_dist ;
  l_spring = parms.l_spring ;
  MRISpositionSurface(mris, mri_dist, mri_masked_smooth, &parms) ;

  if (target_label>=0)
  {
    char fname[STRLEN] ;

    sprintf(fname, "label%d.tri", target_label) ;
    strcpy(fname, output_dir) ;  // last argv is name of file to write to
    printf("saving surface to %s\n", fname) ;
    MRISwrite(mris, fname) ;
    exit(0) ;
  }
  else
    MRISwrite(mris, "inner_skull.tri") ;

  if (inner_skull_only) goto done;

  /* 
   * now create outer skull surface 
   */
  fprintf(stderr, "creating outer skull surface ...\n") ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  create_skull_volume(mri_labeled, mri_masked) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_masked, "brain.mgh") ;
  parms.l_surf_repulse = 5 ;  /* don't let outer_skull come close 
                                 to inner skull */
  parms.sigma = 8 ;
  parms.start_t = 0 ;

  create_distance_map(mri_masked, mri_dist, BORDER_VAL, OUTSIDE_BORDER_STEP) ;
  MRISsetVals(mris,parms.sigma) ;
  MRIScopyValToVal2(mris) ;
  MRISsetVals(mris, TARGET_VAL) ;

  mri_kernel = MRIgaussian1d(parms.sigma, 0) ;
  mri_masked_smooth = MRIconvolveGaussian(mri_dist, NULL, mri_kernel) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_masked_smooth, "brain_smooth.mgh") ;

  sprintf(parms.base_name, "%s_outer_skull%s%s", 
          vol_name, output_suffix, suffix) ;
  parms.mri_brain = mri_dist ;
  MRISpositionSurface(mris, mri_dist, mri_masked_smooth, &parms) ;
  parms.l_spring *= 0.5 ;
  MRISpositionSurface(mris, mri_dist, mri_masked_smooth, &parms) ;
  parms.l_spring *= 0.5 ;
  MRISpositionSurface(mris, mri_dist, mri_masked_smooth, &parms) ;
  parms.l_spring *= 0.5 ;
  MRISpositionSurface(mris, mri_dist, mri_masked_smooth, &parms) ;
  parms.l_spring = l_spring ;

  MRISwrite(mris, "outer_skull.tri") ;

  /* 
   * now build outer skin surface 
   */
  fprintf(stderr, "building outer skin surface...\n") ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  create_skin_volume(mri_labeled, mri_masked) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_masked, "brain.mgh") ;
  parms.l_surf_repulse = 5 ;  /* don't let outer_skin come close 
                                 to outer skull */
  parms.sigma = 8 ;
  parms.start_t = 0 ;

  create_distance_map(mri_masked, mri_dist, BORDER_VAL, OUTSIDE_BORDER_STEP) ;
  MRISsetVals(mris,parms.sigma) ;
  MRIScopyValToVal2(mris) ;
  MRISsetVals(mris, TARGET_VAL) ;

  mri_kernel = MRIgaussian1d(parms.sigma, 0) ;
  mri_masked_smooth = MRIconvolveGaussian(mri_dist, NULL, mri_kernel) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_masked_smooth, "brain_smooth.mgh") ;

  sprintf(parms.base_name, "%s_outer_skin%s%s", 
          vol_name, output_suffix, suffix) ;
  parms.mri_brain = mri_dist ;
  initialize_surface_position(mris, mri_masked, 1) ;
  MRISpositionSurface(mris, mri_dist, mri_masked_smooth, &parms) ;
  parms.l_spring *= .5 ;
  MRISpositionSurface(mris, mri_dist, mri_masked_smooth, &parms) ;
#if 0
  parms.l_spring *= 0.5 ;
  MRISpositionSurface(mris, mri_dist, mri_masked_smooth, &parms) ;
  parms.l_spring *= 0.5 ;
  MRISpositionSurface(mris, mri_dist, mri_masked_smooth, &parms) ;
#endif

  MRISwrite(mris, "outer_skin.tri") ;

 done:
  msec = then.milliseconds() ;
  fprintf(stderr,"positioning took %2.1f minutes\n", (float)msec/(60*1000.0f));

  exit(0) ;
  return(0) ;  /* for ansi */
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "fine"))
    ic = 5 ;
  else if (!stricmp(option, "coarse"))
    ic = 4 ;
  else if (!stricmp(option, "debug_voxel"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "ic"))
  {
    ic = atoi(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    fprintf(stderr,  "using neighborhood size = %d\n", nbrs) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "shrink"))
  {
    parms.l_shrinkwrap = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "label"))
  {
    target_label = atoi(argv[2]) ;
    printf("shrinkwrapping label %s (%d)\n", 
           cma_label_to_name(target_label), target_label) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "name"))
  {
    strcpy(parms.base_name, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "base name = %s\n", parms.base_name) ;
  }
  else if (!stricmp(option, "dt"))
  {
    parms.dt = atof(argv[2]) ;
    parms.base_dt = base_dt_scale*parms.dt ;
    parms.integration_type = INTEGRATE_MOMENTUM ;
    fprintf(stderr,  "using dt = %2.1e\n", parms.dt) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "spring"))
  {
    parms.l_spring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_spring = %2.3f\n", parms.l_spring) ;
  }
  else if (!stricmp(option, "tsmooth"))
  {
    l_tsmooth = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tsmooth = %2.3f\n", l_tsmooth) ;
  }
  else if (!stricmp(option, "grad"))
  {
    parms.l_grad = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_grad = %2.3f\n", parms.l_grad) ;
  }
  else if (!stricmp(option, "tspring"))
  {
    parms.l_tspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tspring = %2.3f\n", parms.l_tspring) ;
  }
  else if (!stricmp(option, "nspring"))
  {
    parms.l_nspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nspring = %2.3f\n", parms.l_nspring) ;
  }
  else if (!stricmp(option, "curv"))
  {
    parms.l_curv = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_curv = %2.3f\n", parms.l_curv) ;
  }
  else if (!stricmp(option, "smooth"))
  {
    smooth = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "smoothing for %d iterations\n", smooth) ;
  }
  else if (!stricmp(option, "output"))
  {
    output_suffix = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "appending %s to output names...\n", output_suffix) ;
  }
  else if (!stricmp(option, "intensity"))
  {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_intensity = %2.3f\n", parms.l_intensity) ;
  }
  else if (!stricmp(option, "embed"))
  {
    embed = 1 ;
    printf("embedding input in 2x blank volume\n") ;
  }
  else if (!stricmp(option, "lm"))
  {
    parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
    fprintf(stderr, "integrating with line minimization\n") ;
  }
  else if (!stricmp(option, "inner_skull_only"))
  {
    inner_skull_only = 1; // exit after writing inner_skull.tri
    fprintf(stderr, "will create only inner_skull.tri\n") ;
  }
  else switch (toupper(*option))
    {
    case 'S':
      suffix = argv[2] ;
      fprintf(stderr, "using %s as suffix\n", suffix) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'Q':
      parms.flags |= IPFLAG_NO_SELF_INT_TEST ;
      fprintf(stderr,
              "doing quick (no self-intersection) surface positioning.\n") ;
      break ;
    case 'M':
      parms.integration_type = INTEGRATE_MOMENTUM ;
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "momentum = %2.2f\n", parms.momentum) ;
      break ;
    case 'R':
      l_surf_repulse = atof(argv[2]) ;
      fprintf(stderr, "l_surf_repulse = %2.3f\n", l_surf_repulse) ;
      nargs = 1 ;
      break ;
    case 'B':
      base_dt_scale = atof(argv[2]) ;
      parms.base_dt = base_dt_scale*parms.dt ;
      nargs = 1;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'T':
      threshold = atof(argv[2]) ;
      target_label = BORDER_VAL ;
      nargs = 1 ;
      printf("thresholding input volume at %2.1f\n", threshold) ;
      break ;
    case 'W':
      sscanf(argv[2], "%d", &parms.write_iterations) ;
      nargs = 1 ;
      fprintf(stderr, "write iterations = %d\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'N':
      sscanf(argv[2], "%d", &parms.niterations) ;
      nargs = 1 ;
      fprintf(stderr, "niterations = %d\n", parms.niterations) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      print_help() ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, "usage: %s [options] <volume> <output name>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr,
          "\nThis program produces three surface files which are\n"
          "shrink-wrapped tesselations of the input volume:\n"
          "  inner_skull.tri\n"
          "  outer_skull.tri\n"
          "  outer_skin.tri\n");
  fprintf(stderr,"\t-t <threshold>  apply threshold to image then deform on distance transform\n") ;
          
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr,
          "(see the source code!)\n") ;
  exit(1) ;
}


static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}


static MRI *
create_brain_volume(MRI *mri_labeled, MRI *mri_brain, int target_label)
{
  int x, y, z, label ;

  if (!mri_brain)
  {
    mri_brain = MRIclone(mri_labeled, NULL) ;
    MRIcopyHeader(mri_labeled, mri_brain) ;
  }
  else
    MRIclear(mri_brain) ;

  if (target_label >= 0)
  {
    MRIcopyLabel(mri_labeled, mri_brain, target_label) ;
    MRIbinarize(mri_brain, mri_brain, 1, 0, BORDER_VAL) ;
  }
  else
  {
    for (x = 0 ; x < mri_labeled->width ; x++)
    {
      for (y = 0 ; y < mri_labeled->height ; y++)
      {
        for (z = 0 ; z < mri_labeled->depth ; z++)
        {
          label = MRIgetVoxVal(mri_labeled,x,y,z, 0) ;
          if (IS_BRAIN(label) || label == CSF_SA || label == Dura)
            MRIsetVoxVal(mri_brain, x, y, z, 0, BORDER_VAL) ;
          else
            MRIsetVoxVal(mri_brain, x, y, z, 0, 0) ;
        }
      }
    }
  }

#if 0
  MRIclose(mri_brain, mri_brain) ;  /* remove small holes */
  MRIopen(mri_brain, mri_brain) ;   /* remove small islands */
#else
  /* 2nd level close  - remove small holes */
  if (target_label <0)
  {
    MRIdilate(mri_brain, mri_brain) ;
    MRIdilate(mri_brain, mri_brain) ;
    MRIerode(mri_brain, mri_brain) ;
    MRIerode(mri_brain, mri_brain) ;

    /* 2nd level open - remove small islands  */
    MRIerode(mri_brain, mri_brain) ;
    MRIerode(mri_brain, mri_brain) ;
    MRIdilate(mri_brain, mri_brain) ;
    MRIdilate(mri_brain, mri_brain) ;
    remove_small_segments(mri_brain, mri_brain) ;
  }
#endif

  return(mri_brain) ;
}

static int
initialize_surface_position(MRI_SURFACE *mris, MRI *mri_masked, int outside)
{
  MRI    *mri_dilated ;
  int    x, y, z ;
  double x0, y0, z0, radius, dist, num, max_r ;
  double xs, ys, zs ;

  if (outside)
  {
    mri_dilated = MRIdilate(mri_masked, NULL) ;

    MRIsubtract(mri_dilated, mri_masked, mri_dilated) ;

    num = x0 = y0 = z0 = 0 ;
    for (x = 0 ; x < mri_dilated->width ; x++)
    {
      for (y = 0 ; y < mri_dilated->height ; y++)
      {
        for (z = 0 ; z < mri_dilated->depth ; z++)
        {
          if (MRIgetVoxVal(mri_dilated, x, y, z, 0) > 0)
          {
            MRIvoxelToSurfaceRAS(mri_dilated, x, y, z, &xs, &ys, &zs) ;
            x0 += xs ;
            y0 += ys ;
            z0 += zs ;
            num++ ;
          }
        }
      }
    }
    x0 /= num ;
    y0 /= num ;
    z0 /= num ;
    printf("centroid at (%2.1f, %2.1f, %2.1f)\n", x0,  y0, z0) ;
    max_r = 0 ;

    num = radius = 0 ;
    for (x = 0 ; x < mri_dilated->width ; x++)
    {
      for (y = 0 ; y < mri_dilated->height ; y++)
      {
        for (z = 0 ; z < mri_dilated->depth ; z++)
        {
          if (MRIgetVoxVal(mri_dilated, x, y, z, 0) > 0)
          {
            MRIvoxelToSurfaceRAS(mri_dilated, x, y, z, &xs, &ys, &zs) ;
            dist = sqrt(SQR(xs-x0)+SQR(ys-y0)+SQR(zs-z0)) ;
            radius += dist ;
            if (dist > max_r)
              max_r = dist ;
            num++ ;
          }
        }
      }
    }

    radius /= num ;
    printf("average radius = %2.3f, max = %2.3f\n", radius, max_r) ;

    MRIfree(&mri_dilated) ;

    MRISprojectOntoSphere(mris, mris, max_r*1.1) ;
    MRIStranslate(mris,x0,y0,z0);
    MRIScomputeMetricProperties(mris) ;
  }

  return(NO_ERROR) ;
}

#define MAX_DIST 40
static MRI *
create_distance_map(MRI *mri_masked, 
                    MRI *mri_distance, 
                    int border_val, 
                    int outside_border_step)
{
  int    target_val ;
  int    i ;
  MRI    *mri_dilated = NULL, *mri_tmp = NULL, *mri_tmp2 = NULL ;

  mri_distance = MRIcopy(mri_masked, mri_distance) ;
  MRIreplaceValues(mri_distance, mri_distance, border_val, 255) ;

  /* build outward distances */
  mri_dilated = MRIcopy(mri_masked, NULL) ;
  target_val = border_val-outside_border_step  ;  /* outside border */
  for (i = 0 ; i < MAX_DIST ; i++, target_val -= 4)
  {
    if (target_val <= 0)
      break ;
    mri_tmp = MRIdilate(mri_dilated, mri_tmp) ;
    mri_tmp2 = MRIcopy(mri_tmp, mri_tmp2) ;
    mri_tmp = MRIsubtract(mri_tmp, mri_dilated, mri_tmp) ;  /* mri_tmp should 
                                                               be the next 
                                                               ring */
    MRIreplaceValues(mri_tmp, mri_tmp, border_val, target_val) ;
    MRIcopyLabel(mri_tmp, mri_distance, target_val) ;
    MRIcopy(mri_tmp2, mri_dilated) ;
  }

  /* build inward distances */
  MRIcopy(mri_masked, mri_dilated) ;
  target_val = border_val+outside_border_step  ;  /* outside border */
  for (i = 0 ; i < MAX_DIST ; i++, target_val += 4)
  {
    if (target_val > 255)
      break ;
    mri_tmp = MRIerode(mri_dilated, mri_tmp) ;
    mri_tmp2 = MRIcopy(mri_tmp, mri_tmp2) ;
    mri_tmp = MRIsubtract(mri_tmp, mri_dilated, mri_tmp) ;  /* mri_tmp should 
                                                               be the next 
                                                               ring */

    if (i == 0)
    {
      MRIreplaceValues(mri_tmp, mri_tmp, border_val, border_val) ;
      MRIcopyLabel(mri_tmp, mri_distance, border_val) ;
    }
    else
    {
      MRIreplaceValues(mri_tmp, mri_tmp, border_val, target_val) ;
      MRIcopyLabel(mri_tmp, mri_distance, target_val) ;
    }
    MRIcopy(mri_tmp2, mri_dilated) ;
  }

  MRIfree(&mri_dilated) ;
  MRIfree(&mri_tmp) ;
  MRIfree(&mri_tmp2) ;
  return(mri_distance) ;
}

static MRI *
create_skull_volume(MRI *mri_labeled, MRI *mri_brain)
{
  int x, y, z, label, i ;
  MRI *mri_tmp = NULL ;

  /* don't let outer skull be within 4 mm of brain */
  for (i = 0 ; i < INNER_SKULL_OUTER_SKULL_SEPARATION ; i++)
  {
    mri_tmp = MRIdilate(mri_brain, mri_tmp) ;
    MRIcopy(mri_tmp, mri_brain) ;
  }

  for (x = 0 ; x < mri_labeled->width ; x++)
  {
    for (y = 0 ; y < mri_labeled->height ; y++)
    {
      for (z = 0 ; z < mri_labeled->depth ; z++)
      {
        label = MRIgetVoxVal(mri_labeled,x,y,z, 0) ;
        if (label == Bone || label == Cranium)
          MRIsetVoxVal(mri_brain, x, y, z, 0, BORDER_VAL) ;
      }
    }
  }

  /* 2nd level close  - remove small holes */
  MRIdilate(mri_brain, mri_brain) ;
  MRIdilate(mri_brain, mri_brain) ;
  MRIerode(mri_brain, mri_brain) ;
  MRIerode(mri_brain, mri_brain) ;

  /* 2nd level open - remove small islands  */
  MRIerode(mri_brain, mri_brain) ;
  MRIerode(mri_brain, mri_brain) ;
  MRIdilate(mri_brain, mri_brain) ;
  MRIdilate(mri_brain, mri_brain) ;

  MRIfree(&mri_tmp) ;
  remove_small_segments(mri_brain, mri_brain) ;
  return(mri_brain) ;
}


static MRI *
create_skin_volume(MRI *mri_labeled, MRI *mri_brain)
{
  int x, y, z, label, i ;
  MRI *mri_tmp = NULL ;

  /* don't let outer skull be within 4 mm of brain */
  for (i = 0 ; i < INNER_SKULL_OUTER_SKULL_SEPARATION ; i++)
  {
    mri_tmp = MRIdilate(mri_brain, mri_tmp) ;
    MRIcopy(mri_tmp, mri_brain) ;
  }

  for (x = 0 ; x < mri_labeled->width ; x++)
  {
    for (y = 0 ; y < mri_labeled->height ; y++)
    {
      for (z = 0 ; z < mri_labeled->depth ; z++)
      {
        label = MRIgetVoxVal(mri_labeled,x,y,z,0) ;
        if (!IS_UNKNOWN(label))
          MRIsetVoxVal(mri_brain, x, y, z, 0, BORDER_VAL) ;
      }
    }
  }

  /* 2nd level close  - remove small holes */
  MRIdilate(mri_brain, mri_brain) ;
  MRIdilate(mri_brain, mri_brain) ;
  MRIerode(mri_brain, mri_brain) ;
  MRIerode(mri_brain, mri_brain) ;

  /* 2nd level open - remove small islands  */
  MRIerode(mri_brain, mri_brain) ;
  MRIerode(mri_brain, mri_brain) ;
  MRIdilate(mri_brain, mri_brain) ;
  MRIdilate(mri_brain, mri_brain) ;

  MRIfree(&mri_tmp) ;
  remove_small_segments(mri_brain, mri_brain) ;
  return(mri_brain) ;
}

static MRI *
remove_small_segments(MRI  *mri_src, MRI *mri_dst)
{
  MRI_SEGMENTATION  *mriseg ;
  int               i, max_vox, same, max_i ;

  same = mri_dst == mri_src ;
  if (same || !mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  mriseg = MRIsegment(mri_src, 1, 255) ;

  max_vox = mriseg->segments[0].nvoxels ;
  max_i = 0 ;
  for (i = 0 ; i < mriseg->nsegments ; i++)
    if (mriseg->segments[i].nvoxels > max_vox)
    {
      max_vox = mriseg->segments[i].nvoxels ;
      max_i = i ;
    }

  MRIsegmentToImage(mri_src, mri_dst, mriseg, max_i) ;

  if (same) /* mri_dst initially same, this proc allocated it for tmp space */
  {
    MRIcopy(mri_dst, mri_src) ;
    MRIfree(&mri_dst) ;
    mri_dst = mri_src ;
  }

  MRIsegmentFree(&mriseg) ;
  return(mri_dst) ;
}

static MRI *
pad_volume(MRI_SURFACE *mris, MRI *mri_src, MRI *mri_dst)
{
  int   vno ;
  double  x, y, z, pad ;
  

  pad = 0.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    MRISvertexToVoxel(mris, &mris->vertices[vno], mri_src, &x, &y, &z) ;
    pad = MAX(pad, -x) ;  pad = MAX(pad, -y) ;  pad = MAX(pad, -z) ; 
    pad = MAX(pad, x-(mri_src->width-1)) ;
    pad = MAX(pad, y-(mri_src->height-1)) ;
    pad = MAX(pad, z-(mri_src->depth-1)) ;
  }

  pad = ceil(pad) ;
  printf("padding volume by %d\n", (int)pad) ;
  mri_dst = MRIextractRegionAndPad(mri_src, NULL, NULL, (int)pad) ;
  {
    double xs, ys, zs, xd, yd, zd ;
    MRISvertexToVoxel(mris, &mris->vertices[0], mri_dst, &xd, &yd, &zd) ;
    MRISvertexToVoxel(mris, &mris->vertices[0], mri_src, &xs, &ys, &zs) ;
    DiagBreak() ;
  }
  return(mri_dst) ;
}

