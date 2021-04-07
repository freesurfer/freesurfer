/*
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
#include "mrisegment.h"
#include "mrisegment.h"
#include "mrimorph.h"
#include "mrinorm.h"
#include "mrisurf.h"
#include "mrisurf_project.h"
#include "icosahedron.h"

#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"
#include "version.h"
#include "cma.h"

static double compute_surface_dist_sse(MRI_SURFACE *mris, MRI *mri_dist) ;
static int MRISrepositionToInnerSkull(MRI_SURFACE *mris, MRI *mri_smooth, INTEGRATION_PARMS *parms) ;


int main(int argc, char *argv[]) ;


#define INNER_SKULL_OUTER_SKULL_SEPARATION 4
#define BORDER_VAL           128
#define OUTSIDE_BORDER_STEP  16
#define TARGET_VAL           (BORDER_VAL-OUTSIDE_BORDER_STEP/2)

static int compute_rigid_gradient(MRI_SURFACE *mris, MRI *mri_dist, double *pdx, double *pdy, double *pdz) ;
static void apply_rigid_gradient(MRI_SURFACE *mris, double dx, double dy, double dz) ;
static int MRISfindOptimalRigidPosition(MRI_SURFACE *mris, MRI *mri_dist, INTEGRATION_PARMS *parms) ;
MRI_SURFACE  *MRISprojectOntoTranslatedSphere(MRI_SURFACE *mris_src,
    MRI_SURFACE *mris_dst, double r,
    double x0, double y0, double z0) ;
int MRISpositionOptimalSphere(MRI_SURFACE *mris, MRI *mri_inner, float sample_dist) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int initialize_surface_position(MRI_SURFACE *mris, MRI *mri_masked, int outside,
                                       INTEGRATION_PARMS *parms) ;

static MRI *MRIfindInnerBoundary(MRI *mri_src, MRI *mri_grad, MRI *mri_dst, float dist) ;
static double compute_surface_sse(MRI_SURFACE *mris, MRI *mri, float sample_dist) ;
const char *Progname ;


static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE    1.0
static float base_dt_scale = BASE_DT_SCALE ;


static const char *suffix = "" ;
static const char *output_suffix = "" ;
static double l_tsmooth = 0.0 ;
static double l_surf_repulse = 5.0 ;

static int smooth = 5 ;

static int nbrs = 2 ;
static int ic_init = 3 ;
static int   finitep(float f) ;
static int
finitep(float f) {
#if 0
  if (!finite(f))
    return(0) ;
  if (fabs(f) > 1e5)
    return(1) ;
  return(1) ;
#else
  return(!devIsnan(f)) ;
#endif
}

int
main(int argc, char *argv[]) {
  char          **av, fname[STRLEN], *T1_fname, *PD_fname, *output_dir, *mdir ;
  int           ac, nargs, msec, s ;
  MRI_SURFACE   *mris ;
  MRI           *mri_flash1, *mri_flash2, *mri_masked, *mri_masked_smooth, *mri_kernel,
  *mri_mean, *mri_dif, *mri_binary, *mri_distance ;
  MRI *mri_smooth, *mri_grad, *mri_inner ;
  Timer then ;
  double        l_spring ;
  MRI_SEGMENTATION *mriseg ;


  nargs = handleVersionOption(argc, argv, "mris_AA_shrinkwrap");
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
  parms.check_tol = 1 ;
  parms.ignore_energy = 1 ;
  parms.dt = 0.5f ;
  parms.base_dt = BASE_DT_SCALE*parms.dt ;

  parms.l_spring_norm = 1 ;
  parms.l_shrinkwrap = 0 ;
  parms.l_intensity = 1 ;

  parms.niterations = 0 ;
  parms.write_iterations = 0 /*WRITE_ITERATIONS */;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = 0.0 /*0.8*/ ;
  parms.l_intensity = 1 ;
  parms.dt_increase = 1.0 /* DT_INCREASE */;
  parms.dt_decrease = 0.50 /* DT_DECREASE*/ ;
  parms.error_ratio = 50.0 /*ERROR_RATIO */;
  /*  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;*/
  parms.l_surf_repulse = 0.0 ;
  parms.l_repulse = 0 /*1*/ ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  mdir = getenv("FREESURFER_HOME") ;
  if (!mdir)
    ErrorExit(ERROR_BADPARM, "FREESURFER_HOME not defined in environment") ;

  if (argc < 4)
    usage_exit() ;

  /* set default parameters for white and gray matter surfaces */
  parms.niterations = 1000 ;
  if (parms.momentum < 0.0)
    parms.momentum = 0.0 /*0.75*/ ;

  then.reset() ;
  T1_fname = argv[1] ;
  PD_fname = argv[2] ;
  output_dir = argv[3] ;
  fprintf(stderr, "reading volume %s...\n", T1_fname) ;
  mri_flash1 = MRIread(T1_fname) ;
  if (!mri_flash1)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s", Progname, T1_fname) ;

  mri_flash2 = MRIread(PD_fname) ;
  if (!mri_flash2)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s", Progname, T1_fname) ;

  //  setMRIforSurface(mri_flash1);
  sprintf(fname, "%s/lib/bem/ic%d.tri", mdir, ic_init) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read icosahedron %s", Progname, fname) ;

  mri_mean = MRImean(mri_flash1, NULL, 5) ;
  MRIwrite(mri_mean, "mean.mgz") ;

  mri_dif = MRIabsdiff(mri_flash1, mri_flash2, NULL) ;
  MRIwrite(mri_dif, "dif.mgz") ;

  mriseg = MRIsegment(mri_mean, 30, 100000) ;
  s = MRIsegmentMax(mriseg) ;
  mri_masked = MRIsegmentToImage(mri_flash1, NULL, mriseg, s) ;
  MRIwrite(mri_masked, "mask.mgz") ;
  MRIsegmentFree(&mriseg) ;

  // MRIthresholdMask(mri_dif, mri_masked, mri_dif, 1, 0) ;
  // MRIwrite(mri_dif, "dif_masked.mgz") ;


  mri_kernel = MRIgaussian1d(2, 0) ;
  mri_smooth = MRIconvolveGaussian(mri_dif, NULL, mri_kernel) ;
  MRIwrite(mri_smooth, "smooth.mgz") ;
  MRIScopyVolGeomFromMRI(mris, mri_smooth) ;
  mris->useRealRAS = 1 ;

  initialize_surface_position(mris, mri_dif, 1, &parms) ;
  MRISwrite(mris, "init") ;
  MRISrepositionToInnerSkull(mris, mri_smooth, &parms) ;

  exit(0) ;
  mri_grad = MRIsobel(mri_smooth, NULL, NULL) ;
  MRIwrite(mri_grad, "grad.mgz") ;
  mri_inner = MRIfindInnerBoundary(mri_dif, mri_grad, NULL, 5.0) ;
  MRIwrite(mri_inner, "inner.mgz") ;
  MRIbinarize(mri_inner, mri_inner, 10, 0, 128) ;

  MRISpositionOptimalSphere(mris, mri_inner, 6) ;
  MRISwrite(mris, "optimal") ;
  exit(0) ;
  parms.sigma = 4 / mri_flash1->xsize ;
  // mri_dist = create_distance_map(mri_masked, NULL, BORDER_VAL, OUTSIDE_BORDER_STEP) ;
  MRISsetVals(mris,parms.sigma) ;
  MRIScopyValToVal2(mris) ;
  MRISsetVals(mris, 0) ;
  sprintf(parms.base_name, "%s_inner_skull%s%s", "test", output_suffix, suffix) ;
  parms.mri_brain = mri_masked ;
  l_spring = parms.l_spring_norm ;
  mri_kernel = MRIgaussian1d(parms.sigma, 0) ;


  mri_binary = MRIbinarize(mri_dif, mri_binary, 40, 0, 128) ;
  MRIwrite(mri_binary, "bin.mgz") ;
  mri_distance = MRIdistanceTransform(mri_binary, NULL, 128, 100, DTRANS_MODE_SIGNED, NULL) ;
  MRIwrite(mri_distance, "dist.mgz") ;
  mri_masked_smooth = MRIconvolveGaussian(mri_distance, NULL, mri_kernel) ;
  MRIfree(&mri_kernel) ;
  MRIwrite(mri_masked_smooth, "dif_smooth.mgz") ;

  MRISwrite(mris, "inner_skull.tri") ;

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
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "fine"))
    ic_init = 5 ;
  else if (!stricmp(option, "coarse"))
    ic_init = 4 ;
  else if (!stricmp(option, "ic")) {
    ic_init = atoi(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "nbrs")) {
    nbrs = atoi(argv[2]) ;
    fprintf(stderr,  "using neighborhood size = %d\n", nbrs) ;
    nargs = 1 ;
  } else if (!stricmp(option, "shrink")) {
    parms.l_shrinkwrap = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "name")) {
    strcpy(parms.base_name, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "base name = %s\n", parms.base_name) ;
  } else if (!stricmp(option, "tol")) {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "tol = %f\n", parms.tol) ;
  } else if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else if (!stricmp(option, "dt")) {
    parms.dt = atof(argv[2]) ;
    parms.base_dt = base_dt_scale*parms.dt ;
    parms.integration_type = INTEGRATE_MOMENTUM ;
    fprintf(stderr,  "using dt = %2.1e\n", parms.dt) ;
    nargs = 1 ;
  } else if (!stricmp(option, "spring")) {
    parms.l_spring_norm = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_spring = %2.3f\n", parms.l_spring_norm) ;
  } else if (!stricmp(option, "tsmooth")) {
    l_tsmooth = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tsmooth = %2.3f\n", l_tsmooth) ;
  } else if (!stricmp(option, "grad")) {
    parms.l_grad = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_grad = %2.3f\n", parms.l_grad) ;
  } else if (!stricmp(option, "tspring")) {
    parms.l_tspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tspring = %2.3f\n", parms.l_tspring) ;
  } else if (!stricmp(option, "nspring")) {
    parms.l_nspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nspring = %2.3f\n", parms.l_nspring) ;
  } else if (!stricmp(option, "curv")) {
    parms.l_curv = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_curv = %2.3f\n", parms.l_curv) ;
  } else if (!stricmp(option, "smooth")) {
    smooth = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "smoothing for %d iterations\n", smooth) ;
  } else if (!stricmp(option, "output")) {
    output_suffix = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "appending %s to output names...\n", output_suffix) ;
  } else if (!stricmp(option, "intensity")) {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_intensity = %2.3f\n", parms.l_intensity) ;
  } else if (!stricmp(option, "lm")) {
    parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
    fprintf(stderr, "integrating with line minimization\n") ;
  } else switch (toupper(*option)) {
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
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr, "usage: %s [options] <T1 vol> <PD vol> <output dir>\n",Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program positions the tessellation of the cortical surface\n"
          "at the white matter surface, then the gray matter surface\n"
          "and generate surface files for these surfaces as well as a\n"
          "'curvature' file for the cortical thickness, and a surface file\n"
          "which approximates layer IV of the cortical sheet.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr,
          "-q    omit self-intersection and only generate "
          "gray/white surface.\n") ;
  fprintf(stderr,
          "-c    create curvature and area files from white matter surface\n"
         );
  fprintf(stderr,
          "-a <avgs>   average curvature values <avgs> times (default=10)\n");
  fprintf(stderr,
          "-whiteonly  only generate white matter surface\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}


static int
initialize_surface_position(MRI_SURFACE *mris, MRI *mri_masked, int outside, INTEGRATION_PARMS *parms) {
  double radius = 0;
  
  if (outside) {

    MRI    *mri_dilated ;
    mri_dilated = MRIdilate(mri_masked, NULL) ;

    MRIsubtract(mri_dilated, mri_masked, mri_dilated) ;
    MRIwrite(mri_dilated, "outside.mgz") ;

    int    x, y, z;
    double x0, y0, z0, dist, num ;
    double xs, ys, zs ;

    num = x0 = y0 = z0 = 0 ;
    for (x = 0 ; x < mri_dilated->width ; x++) {
      for (y = 0 ; y < mri_dilated->height ; y++) {
        for (z = 0 ; z < mri_dilated->depth ; z++) {
          if (MRIgetVoxVal(mri_dilated, x, y, z,0) > 0) {
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

    num = radius = 0 ;
    for (x = 0 ; x < mri_dilated->width ; x++) {
      for (y = 0 ; y < mri_dilated->height ; y++) {
        for (z = 0 ; z < mri_dilated->depth ; z++) {
          if (MRIgetVoxVal(mri_dilated, x, y, z,0) > 0) {
            MRIvoxelToSurfaceRAS(mri_dilated, x, y, z, &xs, &ys, &zs) ;
            dist = sqrt(SQR(xs-x0)+SQR(ys-y0)+SQR(zs-z0)) ;
            radius += dist ;
            num++ ;
          }
        }
      }
    }

    radius /= num ;
    printf("average radius = %2.3f\n", radius) ;

    MRIfree(&mri_dilated) ;
    
    MRISprojectOntoSphere(mris, mris, radius*1.25) ;
    MRIStranslate(mris, x0,y0,z0);
    MRIScomputeMetricProperties(mris) ;
  }

  parms->target_radius = radius ;

  MRISsetOriginalXYZfromXYZ(mris);

  return(NO_ERROR) ;
}

static MRI *
MRIfindInnerBoundary(MRI *mri_src, MRI *mri_grad, MRI *mri_dst, float dist) {
  int    x, y, z, nsamples  ;
  double x1, y1, z1, dx, dy, dz, dot, d, inside, outside, norm, val, xc, yc, zc, wt ;

  if (mri_dst == NULL)
    mri_dst = MRIclone(mri_src, NULL) ;

  xc = yc = zc = 0.0 ;
  for (wt = 0.0, x = 0 ; x < mri_src->width ; x++) {
    for (y = 0 ; y < mri_src->height ; y++) {
      for (z = 0 ; z < mri_src->depth ; z++) {
        val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
        xc += val*x ;
        yc += val*y ;
        zc += val*z ;
        wt += val ;
      }
    }
  }
  xc /= wt ;
  yc /= wt ;
  zc /= wt ;

  for (x = 0 ; x < mri_src->width ; x++) {
    for (y = 0 ; y < mri_src->height ; y++) {
      for (z = 0 ; z < mri_src->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        dx = MRIgetVoxVal(mri_grad, x, y, z, 0) ;
        dy = MRIgetVoxVal(mri_grad, x, y, z, 1) ;
        dz = MRIgetVoxVal(mri_grad, x, y, z, 2) ;
        norm = sqrt(dx*dx + dy*dy + dz*dz) ;
        if (DZERO(norm))
          continue ;
        dx /= norm ;
        dy /= norm ;
        dz /= norm ;
        dot = dx * (x-xc) + dy * (y-yc) + dz * (z-zc);
        if (dot < 0)
          continue ;
        for (nsamples = 0, outside = 0, d = .1 ; d < dist ; d += 0.1, nsamples++) {
          x1 = x+d*dx ;
          y1 = y+d*dy ;
          z1 = z+d*dz ;
          MRIsampleVolume(mri_src, x1, y1, z1, &val) ;
          outside += val ;
        }
        for (nsamples = 0, inside = 0, d = .1 ; d < dist ; d += 0.1, nsamples++) {
          x1 = x-d*dx ;
          y1 = y-d*dy ;
          z1 = z-d*dz ;
          MRIsampleVolume(mri_src, x1, y1, z1, &val) ;
          inside += val ;
        }
        outside /= nsamples ;
        inside /= nsamples ;
        MRIsetVoxVal(mri_dst, x, y, z, 0, outside-inside) ;
      }
    }
  }

  return(mri_dst) ;
}

#define DELTA_R 0.5
int
MRISpositionOptimalSphere(MRI_SURFACE *mris, MRI *mri_inner, float sample_dist) {
  double  r, rmin, rmax, min_sse, min_r, min_x0, min_y0, min_z0, x0, y0, z0, sse,
  xmin, xmax, ymin, ymax, zmin, zmax, delta_r, delta_x, delta_y, delta_z ;
  int     scale ;

  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  rmin = mris->radius*.5 ;
  rmax = mris->radius*1.5 ;
  min_x0 = mris->xctr ;
  min_y0 = mris->yctr ;
  min_z0 = mris->zctr ;
  min_r = mris->radius ;

  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISprojectOntoTranslatedSphere(mris, mris, min_r, min_x0, min_y0, min_z0) ;
  min_sse = compute_surface_sse(mris, mri_inner, sample_dist) ;
  delta_r = DELTA_R ;
  for (r = rmin ; r <= rmax ; r += delta_r) {
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISprojectOntoTranslatedSphere(mris, mris, r, min_x0, min_y0, min_z0) ;
    sse = compute_surface_sse(mris, mri_inner, sample_dist) ;
    if (sse < min_sse) {
      min_sse = sse ;
      min_r = r ;
      printf("new min sse %2.0f, found at r=%2.1f mm, c = (%2.1f, %2.1f, %2.1f)\n",
             min_sse, min_r, min_x0, min_y0, min_z0) ;
    }
  }
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISprojectOntoTranslatedSphere(mris, mris, min_r, min_x0, min_y0, min_z0) ;
  MRISwrite(mris, "lh.minr") ;
  for (scale = 8 ; scale >= 1 ; scale /=2) {
    printf("scale = %d\n", scale) ;
    xmin = min_x0-scale ;
    xmax = min_x0+scale ;
    ymin = min_y0-scale ;
    ymax = min_y0+scale ;
    zmin = min_z0-scale ;
    zmax = min_z0+scale ;
    rmin = min_r*(1.0-scale/128.0) ;
    rmax = min_r*(1.0+scale/128.0) ;
    delta_x = mri_inner->xsize*scale/16 ;
    delta_y = mri_inner->ysize*scale/16 ;
    delta_z = mri_inner->zsize*scale/16 ;
    delta_r = DELTA_R*scale ;
    if (delta_x > (xmax-xmin)/3)
      delta_x = (xmax-xmin)/3 ;
    if (delta_y > (ymax-ymin)/3)
      delta_y = (ymax-ymin)/3 ;
    if (delta_z > (zmax-zmin)/3)
      delta_z = (zmax-zmin)/3 ;
    if (delta_r > (rmax-rmin)/3)
      delta_r = (rmax-rmin)/3 ;

    for (x0 = xmin ; x0 <= xmax ; x0 += delta_x) {
      for (y0 = ymin ; y0 <= ymax ; y0 += delta_y) {
        for (z0 = zmin ; z0 <= zmax ; z0 += delta_z) {
          for (r = rmin ; r <= rmax ; r += delta_r) {
            MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
            MRISprojectOntoTranslatedSphere(mris, mris, r, x0, y0, z0) ;
            sse = compute_surface_sse(mris, mri_inner, sample_dist) ;
            if (sse < min_sse) {
              min_sse = sse ;
              min_r = r ;
              min_x0 = x0 ;
              min_y0 = y0 ;
              min_z0 = z0 ;
              printf("new min sse %2.0f, found at r=%2.1f mm, c = (%2.1f, %2.1f, %2.1f)\n",
                     min_sse, min_r, min_x0, min_y0, min_z0) ;
            }
          }
        }
      }
    }
    {
      char fname[STRLEN] ;
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
      MRISprojectOntoTranslatedSphere(mris, mris, min_r, min_x0, min_y0, min_z0) ;
      sprintf(fname, "lh.min_scale%d", scale) ;
      MRISwrite(mris, fname) ;
    }
  }

  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISprojectOntoTranslatedSphere(mris, mris, min_r, min_x0, min_y0, min_z0) ;
  return(NO_ERROR) ;
}


/*
 compute sse that wants bright stuff inside and dark stuff outside
*/
#define DELTA 1 //mm
static double
compute_surface_sse(MRI_SURFACE *mris, MRI *mri, float sample_dist) {
  int       vno, nsamples ;
  VERTEX    *v ;
  float     sse, dist ;
  double    val, xw, yw, zw, x, y, z, nx, ny, nz, error ;

  for (nsamples = 0, sse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    if (!std::isfinite(v->x) || !std::isfinite(v->y) || !std::isfinite(v->z))
      DiagBreak() ;
    // sample outside - want bright stuff out here
    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;
    for (error = 0.0, dist = 0 ; dist <= sample_dist ; dist += DELTA) {
      x = v->x + dist*nx ;
      y = v->y + dist*ny ;
      z = v->z + dist*nz ;
      MRISvertexToVoxel(mris, v, mri, &xw, &yw, &zw);
      MRIsampleVolume(mri, xw, yw, zw, &val) ;
      val = 1000-val ;
      error += (val*val) ;
      nsamples++ ;
    }
    // sample inwards - want dark stuff
    for (dist = DELTA ; dist <= sample_dist ; dist += DELTA) {
      x = v->x - dist*nx ;
      y = v->y - dist*ny ;
      z = v->z - dist*nz ;
      MRISvertexToVoxel(mris, v, mri, &xw, &yw, &zw);
      MRIsampleVolume(mri, xw, yw, zw, &val) ;
      val = 0-val ;
      error += 100*(val*val) ;
      nsamples++ ;
    }
    sse += error / (float)nsamples ;
  }

  return(sse) ;
}

#ifdef TARGET_VAL
#undef TARGET_VAL
#endif
#define TARGET_VAL 128
static int
MRISrepositionToInnerSkull(MRI_SURFACE *mris, MRI *mri_smooth, INTEGRATION_PARMS *parms) {
  MRI   *mri_dist, *mri_bin, *mri_kernel, *mri_bin_smooth, *mri_dist_smooth ;
  float l_spring, sigma ;
  int   i, ic_order, avgs ;

  parms->niterations = 1000 ;
  if (parms->momentum < 0.0)
    parms->momentum = 0.0 /*0.75*/ ;

  mri_bin = MRIbinarize(mri_smooth, NULL, 15, 0, TARGET_VAL) ;
  mri_dist = MRIdistanceTransform(mri_bin, NULL, TARGET_VAL, 10*mri_bin->width, DTRANS_MODE_SIGNED, NULL) ;
  MRIwrite(mri_bin, "bin.mgz") ;
  MRIwrite(mri_dist, "dist.mgz") ;
  MRISscaleBrain(mris, mris, 0.5) ;  // start inside

  mri_kernel = MRIgaussian1d(2, 0) ;
  mri_bin_smooth = MRIconvolveGaussian(mri_bin, NULL, mri_kernel) ;
  MRIwrite(mri_bin_smooth, "bin_smooth.mgz") ;
  MRISfindOptimalRigidPosition(mris, mri_bin_smooth, parms) ;
  MRIfree(&mri_kernel) ;
  MRIfree(&mri_bin_smooth) ;
  avgs = parms->n_averages = 32 ;
  l_spring = parms->l_spring_norm ;
  for (ic_order = 3 ; ic_order <= 3 ; ic_order++) {
    if (ic_order != ic_init) {
      MRI_SURFACE *mris_new ;
      char fname[STRLEN], *mdir ;

      mdir = getenv("FREESURFER_HOME") ;
      if (!mdir)
        ErrorExit(ERROR_BADPARM, "FREESURFER_HOME not defined in environment") ;

      sprintf(fname, "%s/lib/bem/ic%d.tri", mdir, ic_order) ;
      mris_new = MRISread(fname) ;
      MRISupsampleIco(mris, mris_new) ;
      MRISfree(&mris) ;
      mris = mris_new ;
    }

    printf("********************** using ICO order %d *********************\n", ic_order) ;
    parms->n_averages = avgs ;
    parms->l_spring_norm = l_spring ;
    for (sigma = 16.0, i = 0 ; i < 7 ; i++, sigma /= 2) {
      printf("******************** pass %d, sigma = %2.2f, avgs = %d ******************\n",
             i+1, sigma, parms->n_averages) ;
      parms->sigma = sigma ;
      MRISsetVals(mris,parms->sigma) ;
      MRIScopyValToVal2(mris) ;
      MRISsetVals(mris, 0) ;  // 0 mm from fat
      parms->mri_brain = mri_dist ;


      mri_kernel = MRIgaussian1d(sigma, 0) ;
      mri_dist_smooth = MRIconvolveGaussian(mri_dist, NULL, mri_kernel) ;
      MRIfree(&mri_kernel) ;
      if (i == 0) {
        MRIwrite(mri_dist_smooth, "dist_smooth.mgz") ;

        MRISwrite(mris, "lh.0000") ;
      }
      MRISsetVals(mris, 0) ;
      MRISpositionSurface(mris, mri_dist, mri_dist_smooth, parms) ;
      parms->l_spring_norm /= 2;
      parms->n_averages /= 2 ;
    }
  }

  MRIfree(&mri_bin) ;
  MRIfree(&mri_dist) ;
  return(NO_ERROR) ;
}

static int
MRISfindOptimalRigidPosition(MRI_SURFACE *mris, MRI *mri, INTEGRATION_PARMS *parms) {
  double    dx, dy, dz, old_sse, sse, pct_change, dt, dx_total, dy_total, dz_total ;
  int       i ;

  i = 0 ;
  sse = compute_surface_dist_sse(mris, mri) ;
  dt = .1 ;

  dx_total = dy_total = dz_total = 0;
  do {
    old_sse = sse ;
    compute_rigid_gradient(mris, mri, &dx, &dy, &dz) ;
    apply_rigid_gradient(mris, dx*dt, dy*dt, dz*dt) ;
    sse = compute_surface_dist_sse(mris, mri) ;
    dx_total += dt*dx ;
    dy_total += dt*dy ;
    dz_total += dt*dz ;
    pct_change = 100.0 * ((old_sse - sse) / old_sse) ;
    if (pct_change < 0) {
      dt *=-1 ;
      apply_rigid_gradient(mris, dx*dt, dy*dt, dz*dt) ;
      dx_total += dt*dx ;
      dy_total += dt*dy ;
      dz_total += dt*dz ;
    }

    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN] ;
      sprintf(fname, "lh.%4.4d", ++parms->start_t) ;
      MRISwrite(mris, fname) ;
    }
    if (FZERO(sse))
      break ;
  } while (pct_change > .01);

  printf("after rigid positioning, delta = (%2.0f, %2.0f, %2.0f)\n", dx_total, dy_total, dz_total) ;
  return(NO_ERROR) ;
}

static double
compute_surface_dist_sse(MRI_SURFACE *mris, MRI *mri) {
  double sse = 0.0 ;
  int    vno ;
  VERTEX *v ;
  double val, error, xw, yw, zw ;

  for (vno = 0, sse = 0.0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    MRISvertexToVoxel(mris, v, mri, &xw, &yw, &zw);
    MRIsampleVolume(mri, xw, yw, zw, &val) ;
    error = v->val - val ;
    sse += error*error ;
  }

  return(sse) ;
}

static int
compute_rigid_gradient(MRI_SURFACE *mris, MRI *mri, double *pdx, double *pdy, double *pdz) {
  int    vno ;
  VERTEX *v ;
  double val, xw, yw, zw, dx, dy, dz, delV, x, y, z, Ix, Iy, Iz, xv, yv, zv ;

  dx = dy = dz = 0.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    x = v->x ;
    y = v->y ;
    z = v->z ;

    MRISvertexToVoxel(mris, v, mri, &xw, &yw, &zw);
    MRIsampleVolume(mri, xw, yw, zw, &val) ;

    MRIsampleVolumeGradient(mri, xw, yw, zw,  &Ix, &Iy, &Iz) ;
    // convert back to surface coords
    xw += Ix ;
    yw += Iy ;
    zw += Iz ;
    if (mris->useRealRAS)
      MRIvoxelToWorld(mri, xw, yw, zw, &xv, &yv, &zv) ;
    else
      MRIvoxelToSurfaceRAS(mri, xw, yw, zw, &xv, &yv, &zv) ;
    Ix = xv-v->x ;
    Iy = yv-v->y;
    Iz = zv-v->z ;

    delV = v->val - val ;

    dx += delV * Ix ;
    dy += delV * Iy ;
    dz += delV * Iz ;
    if (!finitep((float)dx))
      DiagBreak() ;
    if (!finitep((float)dy))
      DiagBreak() ;
    if (!finitep((float)dz))
      DiagBreak() ;
  }
  dx /= mris->nvertices ;
  dy /= mris->nvertices ;
  dz /= mris->nvertices ;
  *pdx = dx ;
  *pdy = dy ;
  *pdz = dz ;
  return(NO_ERROR) ;
}

static void apply_rigid_gradient(MRIS* mris, double dx, double dy, double dz) {
  MRIStranslate(mris, dx,dy,dz);
}

