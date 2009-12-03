/**
 * @file  mris_exvivo_deform.c
 * @brief program for deforming a surface to lie at the gray/white or pial boundary from
 *  ultra-high res ex vivo data
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2009/12/03 22:25:45 $
 *    $Revision: 1.9 $
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "transform.h"
#include "mrisurf.h"
#include "label.h"

typedef struct
{
  int    which_border ;
  double inside_val ;
  double infra_granular_val;
  double supra_granular_val ;
  int    error_type ;
  char   base_name[STRLEN] ;
  double ncorr_weight ;   // for linear combinations
  int    use_max_grad ;
  double max_intensity_offset ;
} DEFORMATION_PARMS, DP ;

#define ERROR_FUNC_L1         0
#define ERROR_FUNC_L2         1
#define ERROR_FUNC_NORM_CORR  2
#define ERROR_FUNC_NCORR_L1   3

#define MAX_DIST              3.5
#define MAX_INTENSITY_STEPS   3

static int is_outlier(MRI_SURFACE *mris, int vno) ;
static int recompute_target_locations(MRI_SURFACE *mris, MRI *mri) ;
static double find_optimal_distance(MRI_SURFACE *mris, MRI *mri, int vno,
                                    double max_inwards, double max_outwards, DP *dp);
static double compute_linear_combination(MRI_SURFACE *mris,MRI*mri, 
                                         double x0,double y0,double z0, 
                                         double nx, double ny, double nz, 
                                         double *kernel, double dist, double step, 
                                         int nsamples, char *plot_fname);
static double compute_profile_L1(MRI_SURFACE *mris,MRI*mri, double x0,double y0,double z0, 
                                 double nx, double ny, double nz, 
                                 double *kernel, double dist, double step, 
                                  int nsamples, char *plot_fname);
static double compute_profile_L2(MRI_SURFACE *mris,MRI*mri, double x0,double y0,double z0, 
                                 double nx, double ny, double nz, 
                                 double *kernel, double dist, double step, 
                                  int nsamples, char *plot_fname);
static double compute_profile_norm_corr(MRI_SURFACE *mris,MRI* mri, 
                                        double x0,double y0,double z0, 
                                        double nx, double ny, double nz, 
                                        double *kernel, double dist, double step, 
                                        int nsamples, char *plot_fname);
static int construct_profile(double *kernel,double inside_val, double infra_granular_val, 
                             double supra_granular_val,
                             int inside_len, int infra_granular_len, int supra_granular_len);
double MRISnormalIntensityGradient(MRI_SURFACE *mris, MRI *mri, int vno, double xr, double yr, double zr, double nx, double ny, double zn, double sigma) ;
static double compute_targets(MRI_SURFACE *mris, MRI *mri, int contrast_type, double sigma, DP *dp);
static INTEGRATION_PARMS parms ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;

#define T1 0
#define T2 1

#define GRAY_WHITE_BORDER    0
#define LAYER_IV_BORDER      1
#define PIAL_BORDER          2

static int invert = 0 ;
static int contrast_type = T2 ;
static double sigma = 1.0 ;
static double min_inside_val = 70; 
static double max_outside_val = 180 ;
static double max_inside_val = 150 ;
static double min_outside_val = 120 ;

static double inside_val = 110 ;
static double infra_granular_val = 150;
static double supra_granular_val = 180 ;
static int vavgs = 5 ;
static int min_averages = 2 ;
static  LABEL *label = NULL;
static double max_dist = 2 ;
static DP dp ;

int
main(int argc, char *argv[]) {
  char         **av, *cp ;
  int          ac, nargs, max_averages ;
  double       current_sigma, mean_dist ;
  int          msec, minutes, seconds, n_averages, i, j ;
  struct timeb start ;
  MRI          *mri ;
  MRI_SURFACE  *mris ;
  TRANSFORM    *transform ;
  LTA          *lta ;

  dp.max_intensity_offset = 10 ;
  dp.use_max_grad = 0 ;
  dp.error_type = ERROR_FUNC_NORM_CORR ;
  dp.which_border = GRAY_WHITE_BORDER;
  dp.inside_val = inside_val ;
  dp.infra_granular_val = infra_granular_val;
  dp.supra_granular_val = supra_granular_val ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_deform.c,v 1.9 2009/12/03 22:25:45 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  parms.l_spring = 1.0f ;
  parms.niterations = 100 ;
  parms.n_averages = 16 ;
  parms.l_tspring = 1.0f ;
  parms.l_nspring = 0.5f ;
  parms.l_curv = 1.0 ;
  parms.tol = 1e-6 ;
  parms.dt = 0.5f ;
  parms.l_intensity = 0.1 ;
  parms.l_intensity = 1e-5;
  parms.l_location = 10 ;

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

  if (argc < 5)
    usage_exit(1) ;

  mris = MRISread(argv[1]) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s", Progname, argv[1]) ;
  MRISresetNeighborhoodSize(mris, 3) ; // to allow calculation of nbhd stats

  mri = MRIread(argv[2]) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s", Progname, argv[2]) ;

  transform = TransformRead(argv[3]) ;
  if (transform == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", Progname, argv[3]) ;

  cp = strstr(argv[4], "lh.") ;
  if (cp == NULL)
    cp = strstr(argv[4], "rh.") ;
  if (cp == NULL)
    FileNameExtension(argv[4],parms.base_name) ;  // remove hemi (e.g. lh.)
  else
    strcpy(parms.base_name, cp+3) ;
  printf("using %s as base name\n", parms.base_name) ;
  strcpy(dp.base_name, parms.base_name) ;

  TransformInvert(transform, NULL) ;
  if (transform->type == MNI_TRANSFORM_TYPE ||
      transform->type == TRANSFORM_ARRAY_TYPE ||
      transform->type  == REGISTER_DAT) {
    lta = (LTA *)(transform->xform) ;
    
    if (invert) {
      VOL_GEOM vgtmp;
      LT *lt;
      MATRIX *m_tmp = lta->xforms[0].m_L ;
      lta->xforms[0].m_L = MatrixInverse(lta->xforms[0].m_L, NULL) ;
      MatrixFree(&m_tmp) ;
      lt = &lta->xforms[0];
      if (lt->dst.valid == 0 || lt->src.valid == 0) {
        fprintf(stderr, "WARNING:***************************************************************\n");
        fprintf(stderr, "WARNING:dst volume infor is invalid.  Most likely produce wrong inverse.\n");
        fprintf(stderr, "WARNING:***************************************************************\n");
      }
      copyVolGeom(&lt->dst, &vgtmp);
      copyVolGeom(&lt->src, &lt->dst);
      copyVolGeom(&vgtmp, &lt->src);
    }
  }

  if (stricmp(argv[3], "identity.nofile") != 0)
    MRIStransform(mris, NULL, transform, NULL) ;
  {
    double xv, yv, zv ;
    MRISvertexToVoxel(mris, &mris->vertices[0], mri, &xv, &yv, &zv) ;
    DiagBreak() ;
  }
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISsaveVertexPositions(mris, TARGET_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  if (label == NULL)
    label = LabelInFOV(mris, mri, 2*mri->xsize) ;
  LabelRipRestOfSurface(label, mris) ;

  if (dp.which_border == LAYER_IV_BORDER)
  {
    parms.l_surf_repulse = 5.0 ;
  }
  max_averages = parms.n_averages ;
  for (i = 0, current_sigma = sigma, n_averages = max_averages ; 
       n_averages >= min_averages ; 
       n_averages /= 2, current_sigma /= 2, i++)
  {
    char fname[STRLEN];

    parms.n_averages = n_averages ;
    printf("----------- deforming surfaces with %d smoothing iterations, sigma=%2.3f -------\n", 
           n_averages, current_sigma) ;
    for (j = 0 ; j < 3 ; j++)
    {
      mean_dist = compute_targets(mris, mri, contrast_type, current_sigma, &dp) ;

      printf("----------- positioning surface, j=%d, mean_dist=%2.5f ----------\n", j,mean_dist);

      if (Gdiag & DIAG_WRITE)
      {
        static int ino = 0 ;
        sprintf(fname, "%s.%s.targets.%d.mgz",
                mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", parms.base_name,ino);
        MRISwriteValues(mris, fname) ;

        MRISstoreRipFlags(mris) ;
        MRISunrip(mris) ;
        sprintf(fname, "%s.%s.target_distance.%d.mgz",
                mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", parms.base_name, ino);
        MRISwriteD(mris, fname) ;

        MRISrestoreVertexPositions(mris, TARGET_VERTICES) ;
        sprintf(fname, "%s.%s.target.%d",
                mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", parms.base_name, ino);
        printf("writing target locations to %s\n", fname) ;
        MRISwrite(mris, fname) ;
        MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
        MRISrestoreRipFlags(mris) ;
        ino++ ;
      }
      
      if (vavgs) {
        fprintf
          (stderr,
           "averaging target %s for %d iterations...\n",
           dp.use_max_grad ? "vals" : "D", vavgs) ;
        if (dp.use_max_grad)
        {
          MRISaverageMarkedVals(mris, vavgs) ;
          if (Gdiag_no > 0) {
            VERTEX *v ;
            v = &mris->vertices[Gdiag_no] ;
            fprintf
              (stderr,
               "v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f\n",
               Gdiag_no, v->val, v->mean, v->d) ;
          }
        }
        else
        {
          MRISaverageD(mris, vavgs) ;
          recompute_target_locations(mris, NULL) ;
          if (Gdiag_no > 0) {
            VERTEX *v ;
            double xv, yv, zv ;
            v = &mris->vertices[Gdiag_no] ;
            MRISsurfaceRASToVoxelCached(mris, mri, v->targx, v->targy, v->targz, &xv, &yv, &zv);
            fprintf
              (stderr,
               "v %d, target location = (%2.1f, %2.1f, %2.1f), vox=(%2.0f, %2.0f, %2.0f), "
               "mag = %2.1f, dist=%2.2f\n",
               Gdiag_no, v->targx, v->targy, v->targz, xv, yv, zv, v->mean, v->d) ;
        }
        }
      }
      MRISpositionSurface(mris, mri, mri, &parms) ;
#define MIN_MEAN_DIST 0.25
      if (mean_dist < MIN_MEAN_DIST && j > 0)
        break ;
    }
  }
  if (dp.which_border == GRAY_WHITE_BORDER &&
      dp.use_max_grad == 0 && 0)
  {
    printf("*** using gradient information to deform surface to lie at exact g/w boundary ***\n") ;
    dp.use_max_grad = 1 ;
    parms.n_averages = 4 ;
    parms.l_intensity = 0.1 ;
    parms.l_location = 0 ;
    compute_targets(mris, mri, contrast_type, sigma, &dp) ;
    MRISpositionSurface(mris, mri, mri, &parms) ;
  }

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "surface deformation took %d minutes"
          " and %d seconds.\n", minutes, seconds) ;
  MRISunrip(mris) ;
  printf("writing final surface position to %s\n", argv[4]); 
  MRISwrite(mris, argv[4]) ;
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
  if (!stricmp(option, "vavgs")) {
    vavgs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "smoothing values for %d iterations\n", vavgs) ;
  } else if (!stricmp(option, "loc")) {
    parms.l_location = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using location coefficient = %2.3f\n", parms.l_location) ;
  } else if (!stricmp(option, "layerIv")) {
    dp.which_border = LAYER_IV_BORDER ;
    fprintf(stderr, "repositioning surface to layer IV border\n") ;
  } else if (!stricmp(option, "ncorr")) {
    dp.error_type = ERROR_FUNC_NORM_CORR ;
    fprintf(stderr, "using normalized correlation for error functional\n") ;
  } else if (!stricmp(option, "L1")) {
    dp.error_type = ERROR_FUNC_L1 ;
    fprintf(stderr, "using L1 norm for error functional\n") ;
  } else if (!stricmp(option, "L2")) {
    dp.error_type = ERROR_FUNC_L2 ;
    fprintf(stderr, "using L2 norm for error functional\n") ;
  } else if (!stricmp(option, "L1ncorr")) {
    dp.error_type = ERROR_FUNC_NCORR_L1 ;
    dp.ncorr_weight = atof(argv[2]) ;
    fprintf(stderr, "using weighted combination of L1 and ncorr (%2.3f)\n", dp.ncorr_weight) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
  case 'S':
    sigma = atof(argv[2]) ;

    nargs = 1 ;
    printf("using sigma = %2.3f\n", sigma) ;
    break;
  case 'A':
    parms.n_averages = atoi(argv[2]) ;
    nargs = 1 ;
    printf("smoothing gradient %d times\n", parms.n_averages) ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    printf("debugging vertex %d\n", Gdiag_no) ;
    nargs = 1 ;
    break ;
  case 'L':
    label = LabelRead(NULL, argv[2]) ;
    if (label == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read label file %s", Progname,argv[2]) ;
    nargs = 1 ;
    break ;
  case 'W':
    parms.write_iterations = atoi(argv[2]) ;
    nargs = 1 ;
    printf("writing snapshots of deformation every %d iterations\n",parms.write_iterations);
    Gdiag |= DIAG_WRITE ;
    break ;
  case 'I':
    invert = 1 ;
    printf("inverting xform before applying\n") ;
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
  printf("usage: %s [options] <input surface> <input volume> <xform> <output surface>\n",
         Progname) ;
  printf(
         "\t\n") ;
  exit(code) ;
}

static double
compute_targets(MRI_SURFACE *mris, MRI *mri, int contrast_type, double sigma, DP *dp)
{
  VERTEX  *v ;
  int     vno, nfound, nmissed ;
  MRI     *mri_tmp = MRIclone(mri, NULL) ;
  Real    xv, yv, zv, d, target_val, xr, yr, zr, grad, max_grad, 
    target_dist,mean_border, mean_dist, inward_dist, outward_dist, mean_abs_dist ;

  inward_dist = outward_dist = MAX_DIST ;
  switch (dp->which_border)
  {
  default:
  case GRAY_WHITE_BORDER: break ;
    // assume wm surface is about right
  case LAYER_IV_BORDER: inward_dist = 1.0; outward_dist = 1.0; break ;
  case PIAL_BORDER:       inward_dist += MAX_DIST ;   break ;
  }
  MRISclearMarks(mris) ;
  mean_abs_dist = 0.0 ;
  for (mean_dist=mean_border = 0.0, nfound = nmissed = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if ((vno % 5000) == 0)
      printf("%d of %d vertices processed (%2.1f%%)\n",
             vno, mris->nvertices, 100.0*vno/mris->nvertices) ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    MRISvertexToVoxel(mris, v, mri_tmp, &xv, &yv, &zv) ;
    if (MRIindexNotInVolume(mri, xv, yv, zv))
    {
      v->ripflag = 1 ;
      continue ;
    }
    if (dp->use_max_grad)
    {
      max_grad = 0 ; target_val = -1 ; target_dist = -1000 ;
      for (d = -max_dist ; d <= max_dist ; d += .1)
      {
        xr = v->x + d*v->nx ;
        yr = v->y + d*v->ny ;
        zr = v->z + d*v->nz ;
        MRISsurfaceRASToVoxelCached(mris, mri, xr, yr, zr, &xv, &yv, &zv);
        grad = MRISnormalIntensityGradient(mris, mri, vno,
                                           xr, yr, zr, v->nx, v->ny, v->nz, sigma) ;
        if (contrast_type == T2 && grad < 0)
          grad = 0 ;
        else if (contrast_type == T1 && grad > 0)
          grad = 0 ;
        if (grad > max_grad)
        {
          Real val_inside, val_outside, xs, ys, zs ;
          xs = xr-0.5*v->nx ; ys = yr-0.5*v->ny ; zs = zr-0.5*v->nz ;
          MRIsampleVolume(mri, xv, yv, zv, &val_inside) ;
          xs = xr+0.5*v->nx ; ys = yr+0.5*v->ny ; zs = zr+0.5*v->nz ;
          MRIsampleVolume(mri, xv, yv, zv, &val_outside) ;
          if (val_outside > min_outside_val && val_inside < max_inside_val &&
              val_outside < max_outside_val && val_inside > min_inside_val)
          {
            max_grad = grad;
            MRIsampleVolume(mri, xv, yv, zv, &target_val) ;
            target_dist = d ;
          }
        }
      }
    }
    else  // compute correlation with predicted profile
    {
      target_dist = 
        find_optimal_distance(mris, mri, vno, inward_dist, outward_dist, dp) ;
      xr = v->x + target_dist*v->nx ;
      yr = v->y + target_dist*v->ny ;
      zr = v->z + target_dist*v->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri, xr, yr, zr, &xv, &yv, &zv);
      if (MRIindexNotInVolume(mri, xv, yv, zv))
        target_val = -1 ;
      else
        MRIsampleVolume(mri, xv, yv, zv, &target_val) ;
    }

    if (target_val < 0)  // couldn't find a reasonable guess at border
    {
      target_dist = 0.0 ;
      nmissed++ ;
    }
    else
    {
      mean_border += target_val ;
      mean_dist += target_dist ;
      mean_abs_dist += fabs(target_dist) ;
      nfound++ ;
      v->marked = 1 ;
    }

    v->targx = v->x + target_dist*v->nx ;
    v->targy = v->y + target_dist*v->ny ;
    v->targz = v->z + target_dist*v->nz ;
    MRISsurfaceRASToVoxelCached(mris, mri, v->targx, v->targy, v->targz, &xv, &yv, &zv);
    MRIsetVoxVal(mri_tmp, nint(xv), nint(yv), nint(zv), 0, target_val) ;
    v->val2 = sigma ;
    v->val = target_val ;
    v->d = target_dist ;
  }

  if (nfound > 0)
  {
    mean_border /= (double)nfound ;
    mean_dist /= (double)nfound ;
    mean_abs_dist /= (double)nfound ;
  }

  printf("%d vertices found (%2.2f%%, %d missed), mean dist %2.2f (%2.3f), mean val %2.2f\n",
         nfound, 100*nfound / (double)(nfound+nmissed), nmissed, mean_dist,
         mean_abs_dist, mean_border) ;
  if (nfound == 0)
    DiagBreak() ;
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    static int i = 0 ;
    sprintf(fname, "%s.%s.tvals.%3.3d.mgz",
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", dp->base_name,i++) ;
    printf("writing target vals to %s\n", fname) ;
    MRIwrite(mri_tmp, fname) ;
  }
  if (dp->use_max_grad)
    MRISsoapBubbleVals(mris, 100) ;
  else
  {
    int    n, num, outliers ;
    double mn, std ;
    VERTEX *vn ;

    //    MRISsoapBubbleD(mris, 100) ;
    // now do consistency check on distances
    for (outliers = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->ripflag)
        continue ;
      for (mn = std = 0.0, num = n = 0 ; n < v->vtotal ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked == 0 || vn->ripflag)
          continue ;
        num++ ;
        mn += vn->d ;
        std += vn->d * vn->d ;
      }
      if (num <= 1)
        continue ;
      std = sqrt((std - mn*mn/num)/(num-1)) ;
      mn /= num ;
      if (vno == Gdiag_no)
      {
        FILE *fp  ;
        fp = fopen("out.dat", "w") ;
        for (n = 0 ; n < v->vtotal ; n++)
        {
          vn = &mris->vertices[v->v[n]] ;
          if (vn->marked == 0 || vn->ripflag)
            continue ;
          fprintf(fp, "%d %f\n", v->v[n], vn->d) ;
        }
        fclose(fp) ;
      }
      if ((fabs(v->d-mn) > 4*std) || (is_outlier(mris, vno)))
      {
        if (vno == Gdiag_no)
          printf("replacing vno %d distance %2.2f with %2.2f (std=%2.2f)\n",
                 vno, v->d, mn, std) ;

        v->d = mn ;
        outliers++ ;
      }
    }
    printf("%d outliers replaced\n", outliers) ;
    MRIclear(mri_tmp) ;
    recompute_target_locations(mris, mri_tmp) ;
  }
  MRIfree(&mri_tmp) ;
  return(mean_abs_dist) ;
}


double
MRISnormalIntensityGradient(MRI_SURFACE *mris, MRI *mri, int vno, double xr, double yr, double zr, double nx, double ny, double nz, double sigma_mm)
{
  int     n;
  double  dist, step_size, val_outside, val_inside, k, xw, yw, zw, val,
    ktotal_inside, ktotal_outside, two_sigma_sq  ;
  VERTEX  *v ;

  v = &mris->vertices[vno] ;
  step_size = 0.25 ;
  if (step_size > 2*sigma_mm)
    step_size = sigma_mm ;
  two_sigma_sq = 2.0 * sigma_mm * sigma_mm ;
  ktotal_outside = ktotal_inside = 0 ;
  for (n = 0, val_outside = val_inside = 0.0, dist = step_size ;
       dist <= 2*sigma_mm;
       dist += step_size, n++)
  {
    k = exp(-dist*dist/two_sigma_sq) ;
    xw = xr + dist*nx ; yw = yr + dist*ny ; zw = zr + dist*nz ;
    ktotal_outside += k ;
    MRISsurfaceRASToVoxelCached(mris, mri, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsampleVolume(mri, xw, yw, zw, &val) ;
    val_outside += k*val ;

    xw = xr - dist*nx ; yw = yr - dist*ny ; zw = zr - dist*nz ;
    ktotal_inside += k ;
    MRISsurfaceRASToVoxelCached(mris, mri, xw, yw, zw, &xw, &yw, &zw) ;
    MRIsampleVolume(mri, xw, yw, zw, &val) ;
    val_inside += k*val ;
  }
  if (ktotal_inside> 0)
    val_inside /= (double)ktotal_inside ;
  if (ktotal_outside> 0)
    val_outside /= (double)ktotal_outside ;
  return((val_outside-val_inside) / (2*sigma_mm)) ;
}

static double
find_optimal_distance(MRI_SURFACE *mris, MRI *mri, int vno, 
                      double max_inwards, double max_outwards, DP *dp)

{
  double infra_granular_val, supra_granular_val, inside_val ;
  static double *kernel = NULL ;
  double  min_ig_width = .75, max_ig_width = 2, min_sg_width=.75, max_sg_width=1.5, inside_width=1.5 ;
  double  max_dist, dist, step, sse, min_sse, ig_width, sg_width, offset ;
  int     ig_len, sg_len, in_len, ksize, best_ig_len, best_sg_len, min_ig_len, max_ig_len,
    min_sg_len, max_sg_len;
  static int kmax ;
  double  best_intensity_offset=0.0,min_intensity_offset,max_intensity_offset,intensity_step,
    intensity_offset ; 
  VERTEX  *v = &mris->vertices[vno] ;
  double  xr, yr, zr, xv, yv, zv ;
  int error_type = dp->error_type ;
  double (*profile_error_func)(MRI_SURFACE *mris, MRI *mri,double x0,double y0, double z0, 
                               double nx, double ny, double nz, 
                               double *kernel, double dist, double step, int nsamples, 
                               char *fname);

  infra_granular_val = dp->infra_granular_val ;
  supra_granular_val = dp->supra_granular_val ;
  inside_val  = dp->inside_val ;

  switch (error_type)
  {
  default:
  case ERROR_FUNC_L1:        profile_error_func = compute_profile_L1 ;        break ;
  case ERROR_FUNC_L2:        profile_error_func = compute_profile_L2 ;        break ;
  case ERROR_FUNC_NORM_CORR: profile_error_func = compute_profile_norm_corr ; break ;
  case ERROR_FUNC_NCORR_L1:  profile_error_func = compute_linear_combination; break ;
  }
  switch (error_type)
  {
  default:
  case ERROR_FUNC_L1:        
  case ERROR_FUNC_L2:        
    min_intensity_offset = -dp->max_intensity_offset ;
    max_intensity_offset = dp->max_intensity_offset ;
    intensity_step = 1+(max_intensity_offset - min_intensity_offset) / (MAX_INTENSITY_STEPS) ;
    break ;
  case ERROR_FUNC_NORM_CORR: 
  case ERROR_FUNC_NCORR_L1:  
    min_intensity_offset = max_intensity_offset = 0.0 ;
    intensity_step = 1.0 ;
    break ;
  }
  step = mri->xsize/2 ;

  ig_width = (max_ig_width + min_ig_width) / 2.0 ;
  sg_width = (max_sg_width + min_sg_width) / 2.0 ;

  in_len = nint(inside_width / step) ;
  if (kernel == NULL)   // allocate max profile buffer
  {
    ig_len = nint(max_ig_width/step) ; sg_len = nint(max_sg_width/step) ;
    kmax = ig_len + sg_len + in_len ; 
    kernel = (double *)calloc(kmax, sizeof(double)) ;
  }

  /// pick means and initialize min error based on them and 0 distance
  ig_len = nint(ig_width/step) ; sg_len = nint(sg_width/step) ;
  in_len = nint(inside_width / step) ;
  ksize = ig_len + sg_len + in_len ; 
  best_ig_len = ig_len ; best_sg_len = sg_len ;
  construct_profile(kernel, inside_val, infra_granular_val, supra_granular_val,
                    in_len, ig_len, sg_len) ;
  min_sse = (*profile_error_func)(mris, mri, v->x, v->y, v->z, v->nx, v->ny, v->nz,kernel,
                                  0-inside_width, step, ig_len+sg_len+in_len, NULL) ;
  max_dist = 0.0 ;
  switch (dp->which_border)
  {
  default:
  case GRAY_WHITE_BORDER: offset = 0.0 ;             break ;
  case LAYER_IV_BORDER:   offset = ig_width ;          break ;
  case PIAL_BORDER:       offset = ig_width + sg_width ; break ;
  }
  min_ig_len = nint(min_ig_width/step) ;
  max_ig_len = nint(max_ig_width/step) ;
  min_sg_len = nint(min_sg_width/step) ;
  max_sg_len = nint(max_sg_width/step) ;

  for (intensity_offset = min_intensity_offset ; 
       intensity_offset <= max_intensity_offset ;
       intensity_offset += intensity_step)
  {
    infra_granular_val = dp->infra_granular_val+intensity_offset ;
    supra_granular_val = dp->supra_granular_val+intensity_offset ;
    inside_val  = dp->inside_val+intensity_offset ;
    for (ig_len = min_ig_len ; ig_len <= max_ig_len ; ig_len += 2)
    {
      for (sg_len = min_sg_len ; sg_len <= max_sg_len ; sg_len += 2)
      {
        construct_profile(kernel, inside_val, infra_granular_val, supra_granular_val,
                          in_len, ig_len, sg_len) ;
        switch (dp->which_border)
        {
        default:
        case GRAY_WHITE_BORDER: offset = 0.0 ;                    break ;
        case LAYER_IV_BORDER:   offset = ig_len*step ;            break ;
        case PIAL_BORDER:       offset = (ig_len + sg_len)*step ; break ;
        }
        for (dist = -max_inwards ; dist <= max_outwards ; dist+= step)
        {
          xr = v->x + (dist+offset)*v->nx ;
          yr = v->y + (dist+offset)*v->ny ;
          zr = v->z + (dist+offset)*v->nz ;
          MRISsurfaceRASToVoxelCached(mris, mri, xr, yr, zr, &xv, &yv, &zv);
          if (MRIindexNotInVolume(mri, xv, yv, zv))
            continue ;
          sse = (*profile_error_func)(mris, mri, v->x, v->y, v->z, v->nx, v->ny, v->nz,kernel,
                                      dist-inside_width, step, ig_len+sg_len+in_len, NULL) ;
          if (sse < min_sse)
          {
            if (vno == Gdiag_no)
            {
              (*profile_error_func)(mris, mri, v->x, v->y, v->z, v->nx, v->ny, v->nz,kernel,
                                    dist-inside_width, step, ig_len+sg_len+in_len, "profile.dat") ;
            }
            best_intensity_offset = intensity_offset ;
            max_dist = dist ;
            min_sse = sse ;
            best_ig_len = ig_len ;
            best_sg_len = sg_len ;
          }
        }
      }
    }
  }
  switch (dp->which_border)
  {
  default:
  case GRAY_WHITE_BORDER: offset = 0.0 ;                              break ;
  case LAYER_IV_BORDER:   offset = best_ig_len*step ;                 break ;
  case PIAL_BORDER:       offset = (best_ig_len + best_sg_len)*step ; break ;
  }
  if (vno == Gdiag_no)
  {
    infra_granular_val = dp->infra_granular_val+best_intensity_offset ;
    supra_granular_val = dp->supra_granular_val+best_intensity_offset ;
    inside_val  = dp->inside_val+best_intensity_offset ;
    construct_profile(kernel, inside_val, infra_granular_val, supra_granular_val,
                      in_len, best_ig_len, best_sg_len) ;
    (*profile_error_func)(mris, mri, v->x, v->y, v->z, v->nx, v->ny, v->nz,kernel,
                          max_dist-inside_width, step, best_ig_len+best_sg_len+in_len, 
                          "profile.dat") ;
    DiagBreak() ;
  }
  return(max_dist+offset) ;
}



static double
compute_profile_L1(MRI_SURFACE *mris, MRI *mri, double x0, double y0, double z0, 
                      double nx, double ny, double nz, 
                      double *kernel, double dist, double step, int nsamples, 
                      char *fname)
{
  double   val, kval, x, y, z, xv, yv, zv, dx, dy, dz, L1 ;
  int      i ;
  FILE     *fp = NULL ;
  if (fname)
    fp = fopen(fname, "w") ;

  x = x0 + dist*nx ; y = y0 + dist*ny ; z = z0 + dist*nz ;
  MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
  dx = nx*step ; dy = ny*step ; dz = nz*step ;
  MRISsurfaceRASToVoxelCached(mris, mri, x+dx, y+dy, z+dz, &dx, &dy, &dz) ;
  dx -= xv ; dy -= yv ; dz -= zv ;
  for (L1 = 0.0, i = 0 ; i < nsamples ; i++, dist += step)
  {
    kval = *kernel++ ;
    MRIsampleVolume(mri, xv, yv, zv, &val) ;
    xv += dx ; yv += dy ; zv += dz ;
    if (fp)
      fprintf(fp, "%f %f %f %f %f %f\n", dist, kval, val, xv, yv,zv) ;
    val -= kval ; 
    L1 += fabs(val) ;
  }
  if (fp)
    fclose(fp) ;
  return(L1/nsamples) ;
}

static double
compute_profile_L2(MRI_SURFACE *mris, MRI *mri, double x0, double y0, double z0, 
                      double nx, double ny, double nz, 
                      double *kernel, double dist, double step, int nsamples, 
                      char *fname)
{
  double   sse, val, kval, x, y, z, xv, yv, zv, dx, dy, dz ;
  int      i ;
  FILE     *fp = NULL ;
  if (fname)
    fp = fopen(fname, "w") ;

  x = x0 + dist*nx ; y = y0 + dist*ny ; z = z0 + dist*nz ;
  MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
  dx = nx*step ; dy = ny*step ; dz = nz*step ;
  MRISsurfaceRASToVoxelCached(mris, mri, x+dx, y+dy, z+dz, &dx, &dy, &dz) ;
  dx -= xv ; dy -= yv ; dz -= zv ;
  for (sse = 0.0, i = 0 ; i < nsamples ; i++, dist += step)
  {
    kval = *kernel++ ;
    MRIsampleVolume(mri, xv, yv, zv, &val) ;
    xv += dx ; yv += dy ; zv += dz ;
    if (fp)
      fprintf(fp, "%f %f %f %f %f %f\n", dist, kval, val, xv, yv,zv) ;
    val -= kval ; 
    sse += val*val ;
  }
  if (fp)
    fclose(fp) ;
  return(sqrt(sse/nsamples)) ;
}

static double
compute_profile_norm_corr(MRI_SURFACE *mris, MRI *mri, double x0, double y0, double z0, 
                          double nx, double ny, double nz, 
                          double *kernel, double dist, double step, int nsamples, 
                          char *fname)
{
  double   *k, ival, kval, x, y, z, xv, yv, zv, dx, dy, dz, kmean, imean, kstd, istd ;
  int      i, nzero ;
  FILE     *fp = NULL ;
  double   ncorr ;

  if (fname)
    fp = fopen(fname, "w") ;

  x = x0 + dist*nx ; y = y0 + dist*ny ; z = z0 + dist*nz ;
  MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
  dx = nx*step ; dy = ny*step ; dz = nz*step ;
  MRISsurfaceRASToVoxelCached(mris, mri, x+dx, y+dy, z+dz, &dx, &dy, &dz) ;
  dx -= xv ; dy -= yv ; dz -= zv ;
  k = kernel ;
  kmean = imean = kstd = istd = 0.0 ;
  for (nzero = i = 0 ; i < nsamples ; i++, dist += step)
  {
    kval = *kernel++ ;
    MRIsampleVolume(mri, xv, yv, zv, &ival) ;
    if (FZERO(ival))
      nzero++ ;
    xv += dx ; yv += dy ; zv += dz ;
    if (fp)
      fprintf(fp, "%f %f %f %f %f %f\n", dist, kval, ival, xv, yv,zv) ;
    kmean += kval ; imean += ival ;
    kstd += kval*kval ; istd += ival*ival ;
  }
  if (nzero > 0)
    return(-1) ;  // don't let profile extend outside the image
  if (fp)
    fclose(fp) ;
  kmean /= nsamples ; imean /= nsamples ;
  kstd = sqrt((kstd - (kmean*kmean*nsamples)) / (nsamples-1)) ;
  istd = sqrt((istd - (imean*imean*nsamples)) / (nsamples-1)) ;

  kernel = k ;
  MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
  for (ncorr = 0.0, i = 0 ; i < nsamples ; i++, dist += step)
  {
    kval = *kernel++ ;
    MRIsampleVolume(mri, xv, yv, zv, &ival) ;
    xv += dx ; yv += dy ; zv += dz ;
    ncorr += (kval-kmean) * (ival-imean);
  }
  if (FZERO(kstd))
    kstd = 1e-4 ;
  if (FZERO(istd))
    istd = 1e-4 ;
  ncorr = ncorr / ((nsamples-1)*kstd*istd) ;
  return(-ncorr) ;
}

static int
construct_profile(double *kernel,
                  double inside_val, double infra_granular_val, double supra_granular_val,
                  int inside_len, int infra_granular_len, int supra_granular_len)
{
  int i, itotal ;

  for (i = 0 ; i < inside_len ; i++)
    kernel[i] = inside_val ;
  for (itotal = i, i = 0 ; i < infra_granular_len ; i++, itotal++)
    kernel[itotal] = infra_granular_val ;
  for (i = 0 ; i < supra_granular_len ; i++, itotal++)
    kernel[itotal] = supra_granular_val ;
  return(NO_ERROR) ;
}

static double
compute_linear_combination(MRI_SURFACE *mris,MRI*mri, 
                           double x0,double y0,double z0, 
                           double nx, double ny, double nz, 
                           double *kernel, double dist, double step, 
                           int nsamples, char *plot_fname)
{
  double  ncorr, L1 ;

  ncorr = compute_profile_norm_corr(mris, mri, x0, y0, z0,  nx, ny, nz, 
                                    kernel, dist, step, nsamples, plot_fname);
  L1 = compute_profile_L1(mris, mri, x0, y0, z0,  nx, ny, nz, 
                          kernel, dist, step, nsamples, plot_fname);
  return(ncorr*dp.ncorr_weight + L1) ;
  
}

static int
recompute_target_locations(MRI_SURFACE *mris, MRI *mri)
{
  int          vno ;
  VERTEX       *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    v->targx = v->x + v->d*v->nx ;
    v->targy = v->y + v->d*v->ny ;
    v->targz = v->z + v->d*v->nz ;
    if (mri)
    {
      double xv, yv, zv ;

      MRISsurfaceRASToVoxelCached(mris, mri, v->targx, v->targy, v->targz, 
                                  &xv, &yv, &zv);
      MRIsetVoxVal(mri, nint(xv), nint(yv), nint(zv), 0, v->val) ;
      
    }
  }
  return(NO_ERROR) ;
}
#if 0
static int compare_doubles(const void *v1, const void *v2)
{
  double  *dp1, *dp2 ;

  dp1 = (double *)v1 ;
  dp2 = (double *)v2 ;

  if (*dp1 > *dp2)
    return(1) ;
  else if (*dp1 < *dp2)
    return(-1) ;

  return(0) ;
}
#endif
static int compare_abs_doubles(const void *v1, const void *v2)
{
  double  d1, d2 ;

  d1 = *(double *)v1 ; d1 = fabs(d1) ;
  d2 = *(double *)v2 ; d2 = fabs(d2) ;

  if (d1 > d2)
    return(1) ;
  else if (d1 < d2)
    return(-1) ;

  return(0) ;
}

#define MAX_DISTS 100
static int
is_outlier(MRI_SURFACE *mris, int vno)
{
  int     n, num, imax ;
  VERTEX  *vn, *v ;
  double  dists[MAX_DISTS], mn, std, md ;

  memset(dists, 0, sizeof(dists)) ;
  v = &mris->vertices[vno] ;

  for (num = n = 0 ; n < v->vtotal ; n++)
  {
    vn = &mris->vertices[v->v[n]];
    if (vn->ripflag || vn->marked == 0)
      continue ;

    if (num >= MAX_DISTS)
      break ;
    dists[num] = vn->d ;
    num++ ;
  }

  if (vno == Gdiag_no)
  {
    FILE *fp  ;
    fp = fopen("out.dat", "w") ;
    for (n = 0 ; n < v->vtotal ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->marked == 0 || vn->ripflag)
        continue ;
      fprintf(fp, "%d %f\n", v->v[n], vn->d) ;
    }
    fclose(fp) ;
  }
  if (num <= 1)
    return(0) ;
  md = dists[num/2] ;
  for (n = 0 ; n < num ; n++)
    dists[n] -= md ; // remove the median
  qsort(dists, num, sizeof(double), compare_abs_doubles) ;

  imax = num-(int)(ceil(num/6)) ;  // exclude the top 10% when computing std
  if (imax <= 2)
    return(0) ;

  for (mn = std = 0.0, n = 0 ; n <= imax ; n++)
  {
    mn += dists[n] ;
    std += dists[n]*dists[n] ;
  }
  num  = imax+1 ;
  std = sqrt((std - mn*mn/num)/(num-1)) ;
  mn /= num ;
  if (fabs(v->d-md) > 5*std)
    return(1) ;
  return(0) ;
}
   
