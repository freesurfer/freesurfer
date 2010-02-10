/**
 * @file  mris_deform.c
 * @brief program for deforming a surface to lie at the gray/white or pial boundary from 
 *  ultra-high res data
 *
 * Fit a generative piecewise constant model to the data to determine
 * target locations and deform the surface to match them.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2010/02/10 15:32:13 $
 *    $Revision: 1.17 $
 *
 * Copyright (C) 2002-2009,
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
#include "tritri.h"

#define MIN_INTENISTY_PCT_OFFSET 0.15
#define STRIA_OFFSET .0
#define PROFILE_GENERIC   0
#define PROFILE_V1        1

typedef struct
{
  double wtx, wty, wtz ;   // white target locations
  double ptx, pty, ptz ;   // pial target locations
  double wx, wy, wz ;      // white vertex locations
  double px, py, pz ;      // pial vertex locations
  double l4x, l4y, l4z ;
  double l4tx, l4ty, l4tz ; // layer IV target locations
  double white_val ;
  double pial_val ;
  double l4_val ;
  int    found ;
  double white_dist ;
  double pial_dist ;
  double l4_dist ;
  double wm_intensity_offset ;
  double sg_intensity_offset ;
  double ig_intensity_offset ;
  double wm_val ;
  double nx ;
  double ny ;
  double nz ;               // direction connecting white with pial surface
  int    which_profile ;    // what was the best fitting areal profile
  double min_rms ;
  double pval_v1 ;          // p of being in v1
} VERTEX_PARMS, VP ;
typedef struct
{
  int    which_border ;
  double wm_val ;
  double infra_granular_val;
  double supra_granular_val ;
  double stria_val ;
  double outside_val ;
  int    error_type ;
  char   base_name[STRLEN] ;
  double ncorr_weight ;   // for linear combinations
  int    use_max_grad ;
  double max_wm_intensity_offset ;
  double max_sg_intensity_offset ;
  double max_ig_intensity_offset ;
  int    contrast_type ;
  double outside_width ;
  double stria_width ;
  double inside_width ;
  double max_dist ;
  double sigma ; // for gradient calculations
  double intensity_sigma ;  // for priors on white matter intensity
  double max_sg_width ;
  double max_ig_width ;
  double min_sg_width ;
  double min_ig_width ;
  int    fix_wm_intensity_offset ;
  int    dark_csf ;  // for ex vivo scans
  double step ;
  int    try_stria_model ;
  int    use_prior ;
  double max_stria_intensity_offset ;
} DEFORMATION_PARMS, DP ;

#define WM_INTENSITY_OFFSET   100
#define SG_INTENSITY_OFFSET   101
#define IG_INTENSITY_OFFSET   102
#define RMS                   103

#define ERROR_FUNC_L1         0
#define ERROR_FUNC_L2         1
#define ERROR_FUNC_NORM_CORR  2
#define ERROR_FUNC_NCORR_L1   3

#define MAX_DIST              5
#define MAX_INTENSITY_STEPS   5

// these shouldn't overlap with the ones like WHITE_VERTICES defined in mrisurf.h
#define PIAL_TARGETS          100
#define WHITE_TARGETS         101
#define LAYERIV_TARGETS       102

static int compute_normal(MRI_SURFACE *mris, int vno, VERTEX_PARMS *vp) ;
static int compare_doubles(const void *v1, const void *v2);
static double compute_optimal_intensity(double *iprofile, int i0, int i1) ;
static int compute_best_neighborhood_profile(MRI_SURFACE *mris, MRI *mri,
                                             int vno, 
                                             DP *dp) ;
static double find_min_rms(double *kernel, double *intensity_profile, int profile_len, 
                           char *fname, DP *dp, double wm_dist, double wm_offset, 
                           double min_gm_offset, double max_gm_offset, double gm_offset_step, 
                           int which_profile, int wm_len, int ig_len, int sg_len,   
                           double (*profile_error_func)(double *kernel, double *intensity_profile, 
                                                        int nsamples, char *fname,double step, 
                                                        double wm_dist, double *errors));
static LABEL *label_v1(MRI_SURFACE *mris, MRI *mri, DP *dp) ;
static LABEL *vp_make_v1_label(MRI_SURFACE *mris, double thresh) ;
static double vp_mean(MRI_SURFACE *mris, int which) ;
static int vp_set(MRI_SURFACE *mris, int which, double val) ;
static int vp_copy_to_surface_vals(MRI_SURFACE *mris, int which) ;
static int vp_copy_dist_to_surface_vals(MRI_SURFACE *mris, int which) ;
static int vp_copy_to_surface_dist(MRI_SURFACE *mris, int which) ;
static int vp_copy_from_surface_dist(MRI_SURFACE *mris, int which) ;
static int vp_copy_from_surface(MRI_SURFACE *mris, int which_src, int which_dst) ;
static int vp_copy_to_surface(MRI_SURFACE *mris, int which_src, int which_dst) ;
static int is_outlier(MRI_SURFACE *mris, int vno, int which) ;
static int copy_surface_vals_to_vp(MRI_SURFACE *mris, int which) ;
static int recompute_target_locations(MRI_SURFACE *mris, MRI *mri_white, MRI *mri_l4, 
                                      MRI *mri_pial, DP *dp) ;
static double find_optimal_locations(MRI_SURFACE *mris, MRI *mri, int vno,
                                     double max_inwards, double max_outwards, DP *dp,
                                     VP *vp, int skip);
static double compute_linear_combination(double *kernel, double *iprofile, int nsamples, 
                                         char *plot_fname, double step, double wm_dist,
                                         double *errors);
static double compute_profile_L1(double *kernel, double *iprofile, int nsamples, 
                                 char *plot_fname, double step, double wm_dist, double *errors);
static double compute_profile_L2(double *kernel, double *iprofile, int nsamples, 
                                 char *plot_fname, double step, double wm_dist, double *errors);
static double compute_profile_norm_corr(double *kernel, double *iprofile, int nsamples, 
                                        char *plot_fname, double step, double wm_dist, double *errors);
static int construct_model_profile(double *kernel, DP *dp, 
                                   double wm_intensity_offset, double ig_intensity_offset, 
                                   double sg_intensity_offset, 
                                   double stria_intensity_offset, int inside_len,
                                   int infra_granular_len,int supra_granular_len,int out_len,
                                   int which);
double MRISnormalIntensityGradient(MRI_SURFACE *mris, MRI *mri, 
                                   double xr, double yr, double zr, double nx, double ny, 
                                   double zn, double sigma) ;
static double compute_targets(MRI_SURFACE *mris, MRI *mri, double sigma, DP *dp, int skip);
static double compute_intensity_offsets(MRI_SURFACE *mris, MRI *mri, MRI *mri_dist, DP *dp);

static INTEGRATION_PARMS parms ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int fill_aseg(MRI *mri, char *aseg_fname, DP *dp) ;


char *Progname ;
static void usage_exit(int code) ;

#define T1 0
#define T2 1

#define GRAY_WHITE_BORDER    0
#define LAYER_IV_BORDER      1
#define PIAL_BORDER          2

static double v1_thresh = 0.01 ;
static int invert = 0, deform_debug = 0 ;
static char *aseg_fname = NULL ;
static int read_flag = 0 ;
static char *label_fname = NULL ;

// these are only used for dp.use_max_grad != 0 (not the default)
static double min_wm_val = 70; 
static double max_outside_val = 180 ;
static double max_wm_val = 150 ;
static double min_outside_val = 120 ;

static int vavgs = 0 ;
static int min_averages = 2 ;
static  LABEL *label = NULL;
static double max_dist = 2 ;
static DP dp ;

int
main(int argc, char *argv[]) {
  char         **av, *cp ;
  int          ac, nargs, max_averages ;
  double       current_sigma, mean_dist, l_location ;
  int          msec, minutes, seconds, n_averages, i, j, piter, witer, l4iter, skip ;
  struct timeb start ;
  MRI          *mri, *mri_dist ;
  MRI_SURFACE  *mris ;
  TRANSFORM    *transform ;
  LTA          *lta ;
  VERTEX_PARMS *vp ;
  char         base_name[STRLEN], fname[STRLEN], *hemi ;

  witer = piter = l4iter = 0 ;
  dp.max_dist = MAX_DIST ;
  dp.max_ig_intensity_offset = 40 ;
  dp.max_sg_intensity_offset = 50 ;
  dp.max_wm_intensity_offset = 10 ;
  dp.dark_csf = 0 ;
  dp.try_stria_model = 1 ;
  dp.outside_width = 0.0 ;
  dp.use_prior = 0 ;
  dp.inside_width = 1.5 ;
  dp.use_max_grad = 0 ;
  dp.error_type = ERROR_FUNC_L1 ;
  dp.which_border = GRAY_WHITE_BORDER;
  dp.wm_val = 110 ;
  dp.infra_granular_val = 150;
  dp.max_ig_intensity_offset = 10 ;
  dp.max_stria_intensity_offset = 20 ;
  dp.stria_val = 135 ;
  dp.stria_width = 0.3 ;
  dp.supra_granular_val = 180 ;  // for T2*-weighted ex vivo dat
  dp.contrast_type = T2 ;
  parms.sigma = dp.sigma ;
  dp.intensity_sigma = 15 ;
  parms.check_tol = 1 ;
  dp.max_sg_width = 3 ;
  dp.max_ig_width = 3 ;
  dp.min_sg_width = .75 ;
  dp.min_ig_width = .75 ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_deform.c,v 1.17 2010/02/10 15:32:13 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  parms.l_external = 0 ;
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
  l_location = parms.l_location = 10 ;

  Gx = Gy = Gz = -1 ;
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
  hemi = mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh" ;
  vp = (VERTEX_PARMS *)calloc(mris->nvertices, sizeof(VERTEX_PARMS)) ;
  if (vp == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate VERTEX_PARMS", Progname) ;
  parms.user_parms = (void *)vp ;
  for (i = 0 ; i < mris->nvertices ; i++)
    mris->vertices[i].vp = (void *)(&(vp[i])) ;

  // initialize all the  vp fields to something reasonable for diagnostics
  vp_copy_from_surface(mris, CURRENT_VERTICES, WHITE_VERTICES);
  vp_copy_from_surface(mris, CURRENT_VERTICES, PIAL_VERTICES);
  vp_copy_from_surface(mris, CURRENT_VERTICES, LAYERIV_VERTICES);
  vp_copy_from_surface(mris, CURRENT_VERTICES, WHITE_TARGETS);
  vp_copy_from_surface(mris, CURRENT_VERTICES, PIAL_TARGETS);
  vp_copy_from_surface(mris, CURRENT_VERTICES, LAYERIV_TARGETS);

  MRISresetNeighborhoodSize(mris, 3) ; // to allow calculation of nbhd stats
  MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
  MRISsaveVertexPositions(mris, PIAL_VERTICES) ;

  mri = MRIread(argv[2]) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s", Progname, argv[2]) ;

  dp.step = mri->xsize/2 ;
  if (aseg_fname)
    fill_aseg(mri, aseg_fname, &dp) ;

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
  strcpy(base_name, parms.base_name) ;

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

  if (read_flag) // initialize  with previously computed surfaces
  {
    LABEL *v1, *v1_prior ;
    
    if (label_fname)
    {
      v1_prior = LabelRead(NULL, label_fname) ;
      if (v1_prior == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not load V1 prior label %s\n",
                  Progname, label_fname) ;
      LabelCopyStatsToSurface(v1_prior, mris, VERTEX_STATS) ;
      dp.use_prior = 1 ;
    }
    else
      v1_prior = NULL;

    sprintf(fname, "%s.white", argv[4]) ;
    if (MRISreadWhiteCoordinates(mris, fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read white coords from %s",
                Progname, fname) ;
    vp_copy_from_surface(mris, WHITE_VERTICES, WHITE_VERTICES);

    sprintf(fname, "%s.layerIV", argv[4]) ;
    if (MRISreadWhiteCoordinates(mris, fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read layer IV coords from %s",
                Progname, fname) ;
    vp_copy_from_surface(mris, WHITE_VERTICES, LAYERIV_VERTICES);

    sprintf(fname, "%s.pial", argv[4]) ;
    if (MRISreadWhiteCoordinates(mris, fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read pial coords from %s",
                Progname, fname) ;
    vp_copy_from_surface(mris, WHITE_VERTICES, PIAL_VERTICES);
    vp_copy_to_surface(mris, PIAL_VERTICES, PIAL_VERTICES) ;

    if (dp.try_stria_model)
    {
      v1 = label_v1(mris, mri, &dp) ;
      sprintf(fname, "%s.%s.V1.label", hemi, parms.base_name);
      printf("writing v1 estimated position to %s\n", fname) ;
#if 0
      if (v1_prior)
      {
        int n ;
        MRIRmstVals(mris, 0) ;
        LabelCopyStatsToSurface(v1_prior, mris, VERTEX_VALS) ;
        for (n= 0 ; n < v1->n_points ; n++)
          v1->lv[n].stat *= mris->vertices[v1->lv[n].vno].val ;
      }
#endif
      if (v1)
      {
        LabelWrite(v1, fname) ;
        LabelFree(&v1) ;
      }
    }
    exit(0) ;
  }
  if (dp.which_border == LAYER_IV_BORDER)
  {
    parms.l_surf_repulse = 5.0 ;
  }
  max_averages = parms.n_averages ;
  mri_dist = MRIcloneDifferentType(mri, MRI_FLOAT) ;
  for (skip = 4, i = 0, current_sigma = dp.sigma, n_averages = max_averages ; 
       n_averages >= min_averages ; 
       n_averages /= 2, current_sigma /= 2, i++, skip--)
  {
    char fname[STRLEN];

    if (skip < 1)
      skip = 1 ;
    if (current_sigma < 0.5)
      current_sigma = 1 ;
    dp.sigma = parms.sigma = current_sigma ;
    dp.max_dist = MIN(4*dp.sigma, dp.max_dist) ;  // no point going out that far
    parms.n_averages = n_averages ;
    printf("----------- deforming surfaces with %d smoothing iterations, sigma=%2.3f, skip=%d, max_dist=%2.1f -------\n", 
           n_averages, current_sigma, skip, dp.max_dist) ;
#define MAX_ITERATIONS_PER_SCALE    3
    for (j = 0 ; j < MAX_ITERATIONS_PER_SCALE ; j++)
    {
      dp.max_dist = MIN(4*dp.sigma, dp.max_dist) ;  // no point going out that far
      if (0)
      {
        printf("computing distance transform\n") ;
        vp_copy_to_surface(mris, WHITE_VERTICES, CURRENT_VERTICES) ;
        MRIScomputeDistanceToSurface(mris, mri_dist, mri_dist->xsize) ;
        MRIscalarMul(mri_dist, mri_dist, -1) ; // make inside negative
        vp_set(mris, WM_INTENSITY_OFFSET, 0) ;
        compute_intensity_offsets(mris, mri, mri_dist, &dp) ;
        vp_copy_to_surface_vals(mris, WM_INTENSITY_OFFSET) ;
        MRISsoapBubbleVals(mris, 100) ; 
        MRISmedianFilterVals(mris, 1) ;
        copy_surface_vals_to_vp(mris, WM_INTENSITY_OFFSET) ;

        dp.fix_wm_intensity_offset = 1 ;
      }
      mean_dist = compute_targets(mris, mri, current_sigma, &dp, skip) ;
      vp_copy_to_surface_vals(mris, WM_INTENSITY_OFFSET) ;
      MRISsoapBubbleVals(mris, 100) ; 
      MRISmedianFilterVals(mris, 1) ;
      copy_surface_vals_to_vp(mris, WM_INTENSITY_OFFSET) ;
      dp.fix_wm_intensity_offset = 0 ;

      printf("--- positioning surface, skip=%d, dist=%2.0f, j=%d, mean_dist=%2.5f, sigma=%2.1f ---\n", 
             skip, dp.max_dist, j,mean_dist, dp.sigma);

      if (Gdiag & DIAG_WRITE)
      {
        static int ino = 0 ;
        LABEL *v1 ;

        if (dp.try_stria_model)
        {
          v1 = vp_make_v1_label(mris, v1_thresh) ;
          
          if (v1)
          {
            sprintf(fname, "%s.%s.V1.%d.label", hemi, parms.base_name,ino);
            printf("writing v1 estimated position to %s\n", fname) ;
            LabelWrite(v1, fname) ;
            LabelFree(&v1) ;
          }
        }

        vp_copy_dist_to_surface_vals(mris, LAYERIV_VERTICES) ;
        sprintf(fname, "%s.%s.layerIV.height.%d.mgz", hemi, parms.base_name,ino);
        MRISwriteValues(mris, fname) ;

        vp_copy_to_surface_vals(mris, WHITE_VERTICES) ;
        sprintf(fname, "%s.%s.wtargets.%d.mgz", hemi, parms.base_name,ino);
        MRISwriteValues(mris, fname) ;

        vp_copy_to_surface_vals(mris, RMS) ;
        sprintf(fname, "%s.%s.rms.%d.mgz", hemi, parms.base_name,ino);
        MRISwriteValues(mris, fname) ;

        vp_copy_to_surface_vals(mris, PIAL_VERTICES) ;
        sprintf(fname, "%s.%s.ptargets.%d.mgz", hemi, parms.base_name,ino);
        MRISwriteValues(mris, fname) ;
        vp_copy_to_surface_vals(mris, LAYERIV_VERTICES) ;
        sprintf(fname, "%s.%s.l4targets.%d.mgz", hemi, parms.base_name,ino);
        MRISwriteValues(mris, fname) ;

        MRISstoreRipFlags(mris) ; MRISunrip(mris) ;

        vp_copy_to_surface(mris, WHITE_TARGETS, CURRENT_VERTICES) ;
        sprintf(fname, "%s.%s.wtarget.%d", hemi, parms.base_name, ino);
        printf("writing white target locations to %s\n", fname) ;
        MRISwrite(mris, fname) ;

        vp_copy_to_surface(mris, PIAL_TARGETS, CURRENT_VERTICES) ;
        sprintf(fname, "%s.%s.ptarget.%d", hemi, parms.base_name, ino);
        printf("writing pial target locations to %s\n", fname) ;
        MRISwrite(mris, fname) ;

        vp_copy_to_surface(mris, LAYERIV_TARGETS, CURRENT_VERTICES) ;
        sprintf(fname, "%s.%s.l4target.%d", hemi, parms.base_name, ino);
        printf("writing layer IV target locations to %s\n", fname) ;
        MRISwrite(mris, fname) ;

        MRISrestoreRipFlags(mris) ;
        ino++ ;
      }
      
      recompute_target_locations(mris, NULL, NULL, NULL, &dp) ;
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

      // do white matter surface
      sprintf(parms.base_name, "%s.white", base_name) ;
      sprintf(dp.base_name, "%s.white", base_name) ;
      vp_copy_to_surface(mris, WHITE_VERTICES, CURRENT_VERTICES) ;
      vp_copy_to_surface(mris, WHITE_TARGETS, TARGET_VERTICES) ;
      parms.l_surf_repulse = 0 ;
      vp_copy_to_surface_vals(mris, WHITE_VERTICES);
      MRISsoapBubbleVals(mris, 100) ;
      parms.start_t = witer ;
      printf("positioning white matter surface\n") ;
      MRISpositionSurface(mris, mri, mri, &parms) ;  // move white matter surface
      witer = parms.start_t ; 
      vp_copy_from_surface(mris, CURRENT_VERTICES, WHITE_VERTICES) ;
      if (i == 0 && j == 0) // first time only
      {
        vp_copy_from_surface(mris, CURRENT_VERTICES, PIAL_VERTICES);
        vp_copy_from_surface(mris, CURRENT_VERTICES, LAYERIV_VERTICES);
      }

      MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ; /* save white-matter */
      // do pial surface
      sprintf(parms.base_name, "%s.pial", base_name) ;
      sprintf(dp.base_name, "%s.pial", base_name) ;
      vp_copy_to_surface(mris, PIAL_VERTICES, CURRENT_VERTICES) ;
      vp_copy_to_surface(mris, PIAL_TARGETS, TARGET_VERTICES) ;
      parms.l_surf_repulse = 5 ;  // repel from white matter surface
      vp_copy_to_surface_vals(mris, PIAL_VERTICES);
      MRISsoapBubbleVals(mris, 100) ;
      parms.start_t = piter ;
      printf("positioning pial surface\n") ;
      MRISpositionSurface(mris, mri, mri, &parms) ;  // move pial surface
      piter = parms.start_t ;
      vp_copy_from_surface(mris, CURRENT_VERTICES, PIAL_VERTICES) ;
      MRISsaveVertexPositions(mris, PIAL_VERTICES) ;

      // do layer IV surface
      sprintf(parms.base_name, "%s.layerIV", base_name) ;
      sprintf(dp.base_name, "%s.layerIV", base_name) ;
      vp_copy_to_surface(mris, LAYERIV_VERTICES, CURRENT_VERTICES) ;
      vp_copy_to_surface(mris, LAYERIV_TARGETS, TARGET_VERTICES) ;
      vp_copy_to_surface(mris, WHITE_VERTICES, ORIGINAL_VERTICES) ; // for surf repulse
      parms.l_surf_repulse = 5 ;  // repel from white matter surface
#if 0
      parms.l_osurf_repulse = 5 ;  // repel from pial surface inwards
#endif
      vp_copy_to_surface_vals(mris, LAYERIV_VERTICES);
      MRISsoapBubbleVals(mris, 100) ;
      parms.start_t = l4iter ;
      printf("positioning layer IV surface\n") ;
      MRISpositionSurface(mris, mri, mri, &parms) ;  // move layer IV surface
      parms.l_osurf_repulse = 0 ;  // turn off repulsion from pial surface 
      l4iter = parms.start_t ;
      vp_copy_from_surface(mris, CURRENT_VERTICES, LAYERIV_VERTICES) ;
      strcpy(parms.base_name, base_name) ; // restore default

#define MIN_MEAN_DIST 0.25
      mean_dist = 
        (vp_mean(mris, WHITE_VERTICES) +
         vp_mean(mris, PIAL_VERTICES) +
         vp_mean(mris, LAYERIV_VERTICES)) /3 ;
      printf("iteration %d completed with mean dist = %2.3f\n", j, mean_dist) ;
      if (mean_dist < mri->xsize/2 && j > 0)
        break ;
    }
    //    dp.max_dist /= 2 ;
    if (dp.max_dist < 1)
      dp.max_dist = 1 ;
  }
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  MRISunrip(mris) ;
  sprintf(fname, "%s.white", argv[4]) ;
  vp_copy_to_surface(mris, WHITE_VERTICES, CURRENT_VERTICES) ;
  printf("writing white final surface position to %s\n", fname) ;
  MRISwrite(mris, fname) ;

  vp_copy_dist_to_surface_vals(mris, LAYERIV_VERTICES) ;
  sprintf(fname, "%s.%s.layerIV.height.mgz", hemi, parms.base_name);
  printf("writing layer IV height to %s\n", fname) ;
  MRISwriteValues(mris, fname) ;

  vp_copy_to_surface_vals(mris, SG_INTENSITY_OFFSET) ;
  sprintf(fname, "%s.%s.gm_intensity.mgz", hemi, parms.base_name);
  printf("writing GM intensity offset to %s\n", fname) ;
  MRISwriteValues(mris, fname) ;

  vp_copy_to_surface_vals(mris, IG_INTENSITY_OFFSET) ;
  sprintf(fname, "%s.%s.ig_intensity.mgz", hemi, parms.base_name);
  printf("writing infragranular intensity offset to %s\n", fname) ;
  MRISwriteValues(mris, fname) ;

  sprintf(fname, "%s.pial", argv[4]) ;
  vp_copy_to_surface(mris, PIAL_VERTICES, CURRENT_VERTICES) ;
  printf("writing pial final surface position to %s\n", fname) ;
  MRISwrite(mris, fname) ;

  sprintf(fname, "%s.layerIV", argv[4]) ;
  vp_copy_to_surface(mris, LAYERIV_VERTICES, CURRENT_VERTICES) ;
  printf("writing pial final surface position to %s\n", fname) ;
  MRISwrite(mris, fname) ;
  fprintf(stderr, "surface deformation took %d minutes"
          " and %d seconds.\n", minutes, seconds) ;
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
  } else if (!stricmp(option, "max_offset") || !stricmp(option, "max_wm_offset")) {
    dp.max_wm_intensity_offset = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using max intensity offset %2.0f\n", dp.max_wm_intensity_offset) ;
  } else if (!stricmp(option, "max_sg_offset")) {
    dp.max_sg_intensity_offset = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using max gm intensity offset %2.0f\n", dp.max_sg_intensity_offset) ;
  } else if (!stricmp(option, "max_dist")) {
    dp.max_dist = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using maximum search range %2.1f\n", dp.max_dist) ;
  } else if (!stricmp(option, "tol")) {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using tol = %e\n", parms.tol) ;
  } else if (!stricmp(option, "dark_csf")) {
    dp.dark_csf = 1 ;
    fprintf(stderr, "assuming liquid outside brain is dark\n") ;
  } else if (!stricmp(option, "T1")) {
    dp.wm_val = 110 ;
    dp.max_wm_intensity_offset = 10 ;
    dp.max_ig_intensity_offset = 5 ;
    dp.max_sg_intensity_offset = 10 ;
    dp.inside_width = 1.0 ;  // important for calcarine to avoid ventricle
    dp.infra_granular_val = 85;
    dp.supra_granular_val = 70 ;  // for T1-weighted mp-rage
    dp.contrast_type = T1 ;
    dp.outside_width = 1.0 ; // can't always see csf outside s make it small
    dp.outside_val = 30 ;
    fprintf(stderr, "using location coefficient = %2.3f\n", parms.l_location) ;
  } else if (!stricmp(option, "sg")) {
    dp.supra_granular_val = atof(argv[2]) ;
    nargs = 1 ;
    printf("using %2.1f as supra-granular mean value\n", dp.supra_granular_val) ;
    nargs = 1 ;
  } else if (!stricmp(option, "ig")) {
    dp.infra_granular_val = atof(argv[2]) ;
    nargs = 1 ;
    printf("using %2.1f as infra-granular mean value\n", dp.infra_granular_val) ;
    nargs = 1 ;
  } else if (!stricmp(option, "stria")) {
    dp.stria_val = atof(argv[2]) ;
    nargs = 1 ;
    printf("using %2.1f as mean stria value\n", dp.stria_val) ;
    nargs = 1 ;
  } else if (!stricmp(option, "T2")) {
    dp.wm_val = 110 ;
    dp.max_wm_intensity_offset = 20 ;
    dp.max_sg_intensity_offset = 50 ;
    dp.max_ig_intensity_offset = 40 ;
    dp.inside_width = 1.0 ;  // important for calcarine to avoid ventricle
    dp.infra_granular_val = 140;
    dp.supra_granular_val = 180 ;  // for T2-weighted
    dp.contrast_type = T2 ;
    dp.outside_width = 1.0 ; // can't always see csf outside s make it small
    dp.outside_val = 60 ;
    fprintf(stderr, "using location coefficient = %2.3f\n", parms.l_location) ;
  } else if (!stricmp(option, "aseg")) {
    aseg_fname = argv[2] ;
    fprintf(stderr, "filling ventricles using aseg %s\n", aseg_fname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "layerIv")) {
    dp.which_border = LAYER_IV_BORDER ;
    fprintf(stderr, "repositioning surface to layer IV border\n") ;
  } else if (!stricmp(option, "pial")) {
    dp.which_border = PIAL_BORDER ;
    dp.outside_width = 1.0 ; // can't always see csf outside s make it small
    fprintf(stderr, "repositioning surface to pial surface\n") ;
  } else if (!stricmp(option, "ncorr")) {
    dp.error_type = ERROR_FUNC_NORM_CORR ;
    fprintf(stderr, "using normalized correlation for error functional\n") ;
  } else if (!stricmp(option, "L1")) {
    dp.error_type = ERROR_FUNC_L1 ;
    fprintf(stderr, "using L1 norm for error functional\n") ;
  } else if (!stricmp(option, "L2")) {
    dp.error_type = ERROR_FUNC_L2 ;
    fprintf(stderr, "using L2 norm for error functional\n") ;
  } else if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else if (!stricmp(option, "v1")) {
    label_fname = argv[2] ;
    nargs = 1 ;
    printf("reading prior probabilities from %s\n", label_fname) ;
  } else if (!stricmp(option, "max_grad")) {
    dp.use_max_grad = 1 ;
    fprintf(stderr, "using gradient max as target value location\n") ;
  } else if (!stricmp(option, "L1ncorr")) {
    dp.error_type = ERROR_FUNC_NCORR_L1 ;
    dp.ncorr_weight = atof(argv[2]) ;
    fprintf(stderr, "using weighted combination of L1 and ncorr (%2.3f)\n", dp.ncorr_weight) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
  case 'S':
    dp.sigma = atof(argv[2]) ;

    nargs = 1 ;
    printf("using sigma = %2.3f\n", dp.sigma) ;
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
  case 'R':
    read_flag = 1 ;
    printf("reading surfaces to initialize\n") ;
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


/* do optimal model fit just to find the intensity offsets so
   that spatial smoothness constraints can be imposed on them.
*/
static double
compute_intensity_offsets(MRI_SURFACE *mris, MRI *mri, MRI *mri_dist, DP *dp)
{
  VERTEX  *v ;
  VERTEX_PARMS *vp = NULL ;
  int     vno, nfound, nmissed, peak, wsize ;
  MRI     *mri_white = MRIclone(mri, NULL) ;
  MRI     *mri_pial = MRIclone(mri, NULL) ;
  HISTOGRAM *h, *hs ;
  MRI     *mri_l4 = MRIclone(mri, NULL) ;
  Real    xv, yv, zv, target_val, 
    target_dist,mean_white_border, mean_pial_border, mean_l4_border,
    mean_dist, inward_dist, outward_dist, mean_abs_dist ;

#define WSIZE_MM  5 // diameter of region to histogram within
  wsize = WSIZE_MM / mri->xsize ;  
  inward_dist = outward_dist = dp->max_dist ;
  outward_dist = MIN(4*dp->sigma, outward_dist) ;
  MRIScomputeSurfaceNormals(mris, CURRENT_VERTICES, 10) ;
  MRIScomputeSurfaceNormals(mris, WHITE_VERTICES, 10) ;
  MRIScomputeSurfaceNormals(mris, PIAL_VERTICES, 10) ;
  MRISclearMarks(mris) ;
  mean_abs_dist = 0.0 ;
  mean_white_border = mean_l4_border = mean_pial_border = 0 ;
  for (mean_dist=0.0, nfound = nmissed = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if ((vno % 5000) == 0)
      printf("%d of %d vertices processed (%2.1f%%)\n",
             vno, mris->nvertices, 100.0*vno/mris->nvertices) ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)  
      DiagBreak() ;          // break here XXX
    if (v->ripflag)
      continue ;
    MRISvertexToVoxel(mris, v, mri, &xv, &yv, &zv) ;
    if (MRIindexNotInVolume(mri, xv, yv, zv))
    {
      v->ripflag = 1 ;
      continue ;
    }
    vp = (VERTEX_PARMS *)(v->vp) ;
    h = MRIhistogramVoxel(mri, 0, NULL, nint(xv), nint(yv), nint(zv), wsize, mri_dist, 0) ;
    hs = HISTOsmooth(h, NULL, 2) ;
    peak = HISTOfindHighestPeakInRegion(hs, 0, hs->nbins) ;
    if (peak < 0)
      DiagBreak() ;
    else
      vp->wm_val = hs->bins[peak] ;
    // compute correlation with predicted profile
    if (peak > 0 && fabs(vp->wm_val-dp->wm_val) < dp->max_wm_intensity_offset)
    {
      vp->wm_intensity_offset = vp->wm_val - dp->wm_val ;
      continue ;
    }
    target_dist = 
      find_optimal_locations(mris, mri, vno, inward_dist, outward_dist, dp, vp, 2) ;
    if (vno == Gdiag_no)
      DiagBreak() ; 
    if (target_dist < -1000)
    {
      target_val = -1 ;
      target_dist = 0 ;
    }
    else
    {
      // make sure both targets are in volume
      MRISsurfaceRASToVoxelCached(mris, mri, vp->wtx, vp->wty, vp->wtz, &xv, &yv, &zv);
      if (MRIindexNotInVolume(mri, xv, yv, zv))
        target_val = -1 ;
      else
      {
        MRISsurfaceRASToVoxelCached(mris, mri, vp->ptx, vp->pty, vp->ptz, &xv, &yv, &zv);
        if (MRIindexNotInVolume(mri, xv, yv, zv))
          target_val = -1 ;
        else
          target_val = 0 ; // mark it as ok - will fill in later
      }
    }
    HISTOfree(&h) ; HISTOfree(&hs) ;
  
    if (target_val < 0)  // couldn't find a reasonable guess at border
    {
      target_dist = 0.0 ;
      nmissed++ ;
      vp->found = 0 ;
    }
    else
    {
      mean_dist += target_dist ;
      mean_abs_dist += fabs(target_dist) ;
      nfound++ ;
      v->marked = 1 ;
      MRISsurfaceRASToVoxelCached(mris, mri, vp->wtx, vp->wty, vp->wtz, &xv, &yv, &zv);
      MRIsampleVolume(mri, xv, yv, zv, &target_val) ;
      MRIsetVoxVal(mri_white, nint(xv), nint(yv), nint(zv), 0, target_val) ;
      vp->white_val = target_val ;
      mean_white_border += target_val ;
      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
        DiagBreak() ;
      
      MRISsurfaceRASToVoxelCached(mris, mri, vp->l4tx, vp->l4ty, vp->l4tz, &xv, &yv, &zv);
      MRIsampleVolume(mri, xv, yv, zv, &target_val) ;
      MRIsetVoxVal(mri_l4, nint(xv), nint(yv), nint(zv), 0, target_val) ;
      mean_l4_border += target_val ;
      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
        DiagBreak() ;
      vp->l4_val = target_val ;
      
      MRISsurfaceRASToVoxelCached(mris, mri, vp->ptx, vp->pty, vp->ptz, &xv, &yv, &zv);
      MRIsampleVolume(mri, xv, yv, zv, &target_val) ;
      mean_pial_border += target_val ;
      MRIsetVoxVal(mri_pial, nint(xv), nint(yv), nint(zv), 0, target_val) ;
      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
        DiagBreak() ;
      vp->pial_val = target_val ;
      vp->found = 1 ;
    }
  }
  
  if (nfound > 0)
  {
    mean_white_border /= (double)nfound ;
    mean_pial_border /= (double)nfound ;
    mean_l4_border /= (double)nfound ;
    mean_dist /= (double)nfound ;
    mean_abs_dist /= (double)nfound ;
  }

  printf("%d vertices found (%2.2f%%, %d missed), mean dist %2.2f (%2.3f), mean val %2.2f : %2.2f : %2.2f\n",
         nfound, 100*nfound / (double)(nfound+nmissed), nmissed, mean_dist,
         mean_abs_dist, mean_white_border, mean_l4_border, mean_pial_border) ;
  if (nfound == 0)
    DiagBreak() ;
  if (Gdiag & DIAG_WRITE && 0)
  {
    char fname[STRLEN] ;
    static int i = 0 ;
    sprintf(fname, "%s.%s.wtvals.%3.3d.mgz",
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", dp->base_name,i) ;
    printf("writing white target vals to %s\n", fname) ;
    MRIwrite(mri_white, fname) ;

    sprintf(fname, "%s.%s.ptvals.%3.3d.mgz",
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", dp->base_name,i) ;
    printf("writing pial target vals to %s\n", fname) ;
    MRIwrite(mri_pial, fname) ;

    i++ ;
  }
  MRIfree(&mri_white) ; MRIfree(&mri_pial) ;
  return(mean_abs_dist) ;
}



static double
compute_targets(MRI_SURFACE *mris, MRI *mri, double sigma, DP *dp, int skip)
{
  VERTEX  *v ;
  VERTEX_PARMS *vp = NULL ;
  int     vno, nfound, nmissed ;
  MRI     *mri_white = MRIclone(mri, NULL) ;
  MRI     *mri_pial = MRIclone(mri, NULL) ;
  MRI     *mri_l4 = MRIclone(mri, NULL) ;
  Real    xv, yv, zv, d, target_val, xr, yr, zr, grad, max_grad, 
    target_dist,mean_white_border, mean_pial_border, mean_l4_border,
    mean_dist, inward_dist, outward_dist, mean_abs_dist, base_sigma = dp->sigma ;
  
  inward_dist = outward_dist = dp->max_dist ;
  MRIScomputeSurfaceNormals(mris, CURRENT_VERTICES, 10) ;
  MRIScomputeSurfaceNormals(mris, WHITE_VERTICES, 10) ;
  MRIScomputeSurfaceNormals(mris, PIAL_VERTICES, 10) ;
  MRISclearMarks(mris) ;
  mean_abs_dist = 0.0 ;
  mean_white_border = mean_l4_border = mean_pial_border = 0 ;
  for (mean_dist=0.0, nfound = nmissed = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if ((vno % 5000) == 0)
      printf("%d of %d vertices processed (%2.1f%%)\n",
             vno, mris->nvertices, 100.0*vno/mris->nvertices) ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)  
      DiagBreak() ;          // break here XXX
    if (v->ripflag)
      continue ;
    MRISvertexToVoxel(mris, v, mri, &xv, &yv, &zv) ;
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
        grad = MRISnormalIntensityGradient(mris, mri,
                                           xr, yr, zr, v->nx, v->ny, v->nz, sigma) ;
        if (dp->contrast_type == T2 && grad < 0)
          grad = 0 ;
        else if (dp->contrast_type == T1 && grad > 0)
          grad = 0 ;
        if (grad > max_grad)
        {
          Real val_inside, val_outside, xs, ys, zs ;
          xs = xr-0.5*v->nx ; ys = yr-0.5*v->ny ; zs = zr-0.5*v->nz ;
          MRIsampleVolume(mri, xv, yv, zv, &val_inside) ;
          xs = xr+0.5*v->nx ; ys = yr+0.5*v->ny ; zs = zr+0.5*v->nz ;
          MRIsampleVolume(mri, xv, yv, zv, &val_outside) ;
          if (val_outside > min_outside_val && val_inside < max_wm_val &&
              val_outside < max_outside_val && val_inside > min_wm_val)
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
#define MAX_SCALES 4
      vp = (VERTEX_PARMS *)(v->vp) ;
      
      compute_normal(mris, vno, vp) ;
      if (1)
      {
        target_dist = 
          find_optimal_locations(mris, mri, vno, inward_dist, outward_dist, dp, vp, skip) ;
      }
      else 
      {
        double        sigma, best_rms, best_target_dist ;
        VERTEX_PARMS vps[MAX_SCALES] ;
        int          best_scale,i, found ;

        best_rms = 1e8 ; best_target_dist = -1000 ;
        for (sigma = base_sigma, found = i = 0 ; i < MAX_SCALES ; sigma *= 2, i++)
        {
          dp->sigma = sigma ;
          vps[i] = *vp ;
          target_dist = 
            find_optimal_locations(mris, mri, vno, inward_dist, outward_dist, dp, &vps[i], 2) ;
          if (vps[i].min_rms < best_rms)
          {
            best_rms = vps[i].min_rms ;
            best_scale = i ;
            best_target_dist = target_dist ;
          }
          if (target_dist > -100)  
          {
            found = 1 ;
            if (best_rms < .1*dp->wm_val)
              break ;   // use the smallest scale that produces a feasible answer
            else
              DiagBreak() ;
          }
          DiagBreak() ;
        }
        if (found)
          *vp = vps[best_scale] ;
        dp->sigma = base_sigma ;
        target_dist = best_target_dist ;
      }
      
      if (vno == Gdiag_no)
        DiagBreak() ; 
      if (target_dist < -100)
      {
        target_val = -1 ;
        target_dist = 0 ;
      }
      else
      {
        // make sure both targets are in volume
        MRISsurfaceRASToVoxelCached(mris, mri, vp->wtx, vp->wty, vp->wtz, &xv, &yv, &zv);
        if (vno == Gdiag_no)
          printf("v %d: wm target = %2.0f %2.0f %2.0f\n", vno, xv, yv, zv) ;
        if (MRIindexNotInVolume(mri, xv, yv, zv))
          target_val = -1 ;
        else
        {
          MRISsurfaceRASToVoxelCached(mris, mri, vp->ptx, vp->pty, vp->ptz, &xv, &yv, &zv);
          if (MRIindexNotInVolume(mri, xv, yv, zv))
            target_val = -1 ;
          else
            target_val = 0 ; // mark it as ok - will fill in later
        }
      }
    }

    if (target_val < 0)  // couldn't find a reasonable guess at border
    {
      target_dist = 0.0 ;
      nmissed++ ;
      vp->found = 0 ;
      vp->white_dist = vp->pial_dist = vp->l4_dist = 0 ;
    }
    else
    {
      mean_dist += target_dist ;
      mean_abs_dist += fabs(target_dist) ;
      nfound++ ;
      v->marked = 1 ;
      MRISsurfaceRASToVoxelCached(mris, mri, vp->wtx, vp->wty, vp->wtz, &xv, &yv, &zv);
      MRIsampleVolume(mri, xv, yv, zv, &target_val) ;
      MRIsetVoxVal(mri_white, nint(xv), nint(yv), nint(zv), 0, target_val) ;
      vp->white_val = target_val ;
      mean_white_border += target_val ;
      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
        DiagBreak() ;
      
      MRISsurfaceRASToVoxelCached(mris, mri, vp->l4tx, vp->l4ty, vp->l4tz, &xv, &yv, &zv);
      MRIsampleVolume(mri, xv, yv, zv, &target_val) ;
      MRIsetVoxVal(mri_l4, nint(xv), nint(yv), nint(zv), 0, target_val) ;
      mean_l4_border += target_val ;
      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
        DiagBreak() ;
      vp->l4_val = target_val ;

      MRISsurfaceRASToVoxelCached(mris, mri, vp->ptx, vp->pty, vp->ptz, &xv, &yv, &zv);
      MRIsampleVolume(mri, xv, yv, zv, &target_val) ;
      mean_pial_border += target_val ;
      MRIsetVoxVal(mri_pial, nint(xv), nint(yv), nint(zv), 0, target_val) ;
      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
        DiagBreak() ;
      vp->pial_val = target_val ;
      vp->found = 1 ;
    }

    v->val2 = dp->sigma ;
    v->val = target_val ;
  }

  if (nfound > 0)
  {
    mean_white_border /= (double)nfound ;
    mean_pial_border /= (double)nfound ;
    mean_l4_border /= (double)nfound ;
    mean_dist /= (double)nfound ;
    mean_abs_dist /= (double)nfound ;
  }

  printf("%d vertices found (%2.2f%%, %d missed), mean dist %2.2f (%2.3f), mean val %2.2f : %2.2f : %2.2f\n",
         nfound, 100*nfound / (double)(nfound+nmissed), nmissed, mean_dist,
         mean_abs_dist, mean_white_border, mean_l4_border, mean_pial_border) ;
  if (nfound == 0)
    DiagBreak() ;
  if (dp->use_max_grad)
    MRISsoapBubbleVals(mris, 100) ;
  else
  {
    int    n, num, outliers ;
    double mn, std ;
    VERTEX *vn ;

    // now do consistency check on white distances
    vp_copy_to_surface_dist(mris, WHITE_VERTICES);
    for (outliers = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->ripflag || v->marked == 0)
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
      if ((fabs(v->d-mn) > 4*std) || (is_outlier(mris, vno, WHITE_VERTICES)))
      {
        if (vno == Gdiag_no)
          printf("replacing vno %d distance %2.2f with %2.2f (std=%2.2f)\n",
                 vno, v->d, mn, std) ;

        v->d = mn ;
        outliers++ ;
      }
    }
    printf("%d white outliers replaced\n", outliers) ;

    vp_copy_to_surface_dist(mris, LAYERIV_VERTICES);

    // now do consistency check on pial distances
    for (outliers = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->ripflag || v->marked != 1)
        continue ;
      for (mn = std = 0.0, num = n = 0 ; n < v->vtotal ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked != 1 || vn->ripflag)
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
        fp = fopen("pout.dat", "w") ;
        for (n = 0 ; n < v->vtotal ; n++)
        {
          vn = &mris->vertices[v->v[n]] ;
          if (vn->marked == 0 || vn->ripflag)
            continue ;
          fprintf(fp, "%d %f\n", v->v[n], vn->d) ;
        }
        fclose(fp) ;
      }
      if ((fabs(v->d-mn) > 4*std) || (is_outlier(mris, vno, PIAL_VERTICES)))
      {
        if (vno == Gdiag_no)
          printf("replacing vno %d distance %2.2f with %2.2f (std=%2.2f)\n",
                 vno, v->d, mn, std) ;

        v->d = mn ;
        outliers++ ;
      }
    }

    printf("%d layer IV outliers replaced\n", outliers) ;
    vp_copy_to_surface_dist(mris, LAYERIV_VERTICES);

    // now do consistency check on pial distances
    vp_copy_to_surface_dist(mris, PIAL_VERTICES);
    for (outliers = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->ripflag || v->marked != 1)
        continue ;
      for (mn = std = 0.0, num = n = 0 ; n < v->vtotal ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked != 1 || vn->ripflag)
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
        fp = fopen("pout.dat", "w") ;
        for (n = 0 ; n < v->vtotal ; n++)
        {
          vn = &mris->vertices[v->v[n]] ;
          if (vn->marked == 0 || vn->ripflag)
            continue ;
          fprintf(fp, "%d %f\n", v->v[n], vn->d) ;
        }
        fclose(fp) ;
      }
      if ((fabs(v->d-mn) > 4*std) || (is_outlier(mris, vno, PIAL_VERTICES)))
      {
        if (vno == Gdiag_no)
          printf("replacing vno %d distance %2.2f with %2.2f (std=%2.2f)\n",
                 vno, v->d, mn, std) ;

        v->d = mn ;
        outliers++ ;
      }
    }

    // now do soap bubble smoothing to fill in the missing values
    printf("%d pial outliers replaced\n", outliers) ;
    vp_copy_to_surface_dist(mris, PIAL_VERTICES);
    MRISsoapBubbleD(mris, 100) ;
    MRISmedianFilterD(mris, 1, 0) ;
    vp_copy_from_surface_dist(mris, PIAL_VERTICES) ;

    vp_copy_to_surface_dist(mris, WHITE_VERTICES);
    MRISsoapBubbleD(mris, 100) ;
    MRISmedianFilterD(mris, 1, 0) ;
    vp_copy_from_surface_dist(mris, WHITE_VERTICES) ;

    vp_copy_to_surface_dist(mris, LAYERIV_VERTICES);
    MRISsoapBubbleD(mris, 100) ;     // fill in missing values
    MRISmedianFilterD(mris, 1, 0) ;  // smooth them a bit
    vp_copy_from_surface_dist(mris, LAYERIV_VERTICES) ;

    MRIclear(mri_white) ; MRIclear(mri_pial) ; MRIclear(mri_l4) ;
    recompute_target_locations(mris, mri_white, mri_l4, mri_pial, dp) ;
  }
  if (Gdiag & DIAG_WRITE && 0)
  {
    char fname[STRLEN] ;
    static int i = 0 ;
    sprintf(fname, "%s.%s.wtvals.%3.3d.mgz",
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", dp->base_name,i) ;
    printf("writing white target vals to %s\n", fname) ;
    MRIwrite(mri_white, fname) ;

    sprintf(fname, "%s.%s.ptvals.%3.3d.mgz",
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", dp->base_name,i) ;
    printf("writing pial target vals to %s\n", fname) ;
    MRIwrite(mri_pial, fname) ;

    i++ ;
  }
  MRIfree(&mri_white) ; MRIfree(&mri_pial) ;
  return(mean_abs_dist) ;
}


double
MRISnormalIntensityGradient(MRI_SURFACE *mris, MRI *mri, double xr, double yr, double zr, double nx, double ny, double nz, double sigma_mm)
{
  int     n;
  double  dist, step_size, val_outside, val_inside, k, xw, yw, zw, val,
    ktotal_inside, ktotal_outside, two_sigma_sq  ;

  step_size = mri->xsize/2 ;
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

static int
extract_image_intensities(MRI_SURFACE *mris, MRI *mri, double in_dist, double out_dist, 
                          double x0, double y0, double z0,
                          double nx, double ny, double nz, double step, int tsteps,
                          double *intensity_profile, char *fname)
{
  double   avg_val, val, x, y, z, xv, yv, zv, dx, dy, dz, dist, e1x, e1y, e1z, e2x, e2y, e2z, xt, yt,zt,
    norm, dxn, dyn, dzn, vstep ;
  int      i, nsamples, u, v ;
  FILE     *fp = NULL ;

  if (fname != NULL)
    fp = fopen(fname, "w") ;

  vstep = step / mri->xsize ;  // voxels/step
  x = x0-in_dist*nx ;  y = y0-in_dist*ny ;  z = z0-in_dist*nz ; 
  MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
  dx = nx*step ; dy = ny*step ; dz = nz*step ;
  MRISsurfaceRASToVoxelCached(mris, mri, x+dx, y+dy, z+dz, &dx, &dy, &dz) ;
  dx -= xv ; dy -= yv ; dz -= zv ;  // (dx, dy, dz) is now (nx, ny, nz) in voxel coords
  nsamples = nint((out_dist + in_dist)/step) ;
  norm = sqrt(dx*dx + dy*dy + dz*dz) ; 
  if (FZERO(norm))
    DiagBreak() ;
  dxn = dx/norm ; dyn = dy/norm ; dzn = dy/norm ; 

  e2x = dyn ; e2y = dzn ; e2z = dxn ;
  CROSS3(e1x, e1y, e1z, dxn, dyn, dzn, e2x, e2y, e2z) ;  // find 2 vectors perpindicular d[xyz]
  norm = sqrt(e1x*e1x + e1y*e1y + e1z*e1z) ; 
  if (FZERO(norm))
    DiagBreak() ;
  e1x /= norm ; e1y /= norm ; e1z /= norm ; 
  CROSS3(e2x, e2y, e2z, dxn, dyn, dzn, e1x, e1y, e1z) ;
  norm = sqrt(e2x*e2x + e2y*e2y + e2z*e2z) ; 
  if (FZERO(norm))
    DiagBreak() ;
  e2x /= norm ; e2y /= norm ; e2z /= norm ; 

  for (dist = -in_dist, i = 0 ; i < nsamples ; i++, dist += step)
  {
    // average over the tangent plane
    for (avg_val = 0.0, u = -(tsteps-1) ; u < tsteps ; u++)
    {
      for (v = -(tsteps-1) ; v < tsteps ; v++)
      {
        xt = xv + e1x*step*u + e2x*vstep*v ; 
        yt = yv + e1y*step*u + e2y*vstep*v ; 
        zt = zv + e1z*step*u + e2z*vstep*v ;
        MRIsampleVolume(mri, xt, yt, zt, &val) ;
        avg_val += val ;
      }
    }
    xv += dx ; yv += dy ; zv += dz ;
    intensity_profile[i] = avg_val/((2*tsteps-1)*(2*tsteps-1));
    if (fp)
      fprintf(fp, "%f %f %f %f %f\n", dist, intensity_profile[i], xv, yv, zv) ;
   }

  if (fp)
    fclose(fp) ;
  return(NO_ERROR) ;
}
#if 0
static double
extend_csf_dist_into_dark_regions(MRI_SURFACE *mris, MRI *mri, double x, double y, double z,
                                  double nx, double ny, double nz, double init_dist,
                                  int grad_sign, double max_dist,
                                  double sigma, double step, double *pval, double thresh)
{
  double x1, y1, z1, grad, dist, val, next_val, xv, yv, zv ;

  for (dist = init_dist ; dist <= init_dist+max_dist ; dist += step)
  {
    x1 = x + dist*nx ; y1 = y + dist*ny ; z1 = z + dist*nz ;
    grad = MRISnormalIntensityGradient(mris, mri,  x1,  y1,  z1, nx,  ny,  nz, sigma);
    if (grad*grad_sign < 0)
      break ;
    MRISsurfaceRASToVoxelCached(mris, mri, x1, y1, z1, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &val) ;
    if (val < thresh)
      break;
  }

  // now compute exact location on by looking for the dist where the intensities decrease
  dist -= step/2 ;
  MRISsurfaceRASToVoxelCached(mris, mri, x+dist*nx, y+dist*ny, z+dist*nz, &xv, &yv, &zv) ;
  dist += step/2 ;
  MRIsampleVolume(mri, xv, yv, zv, &val) ;
  for ( ; dist <= max_dist ; dist += step/2)
  {
    MRISsurfaceRASToVoxelCached(mris, mri, x+dist*nx, y+dist*ny, z+dist*nz, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &next_val) ;
    if ((next_val-val)*grad_sign < 0)
      break ;
    if (val < thresh)
      break;
    val = next_val ;
  }
  if (pval)
    *pval = val ;
  return(dist) ;
}

static double
extend_white_matter_distances(MRI_SURFACE *mris, MRI *mri,  double x, double y, double z, double nx,
                              double ny, double nz, double min_val, double max_val, double init_dist,
                              double max_dist, double step)
{
  double dist, val, xv, yv, zv ;

  for (dist = init_dist+step ; dist < init_dist+max_dist ; dist += step)
  {
    MRISsurfaceRASToVoxelCached(mris, mri, x+dist*nx, y+dist*ny, z+dist*nz, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &val) ;
    if (val < min_val || val > max_val)
    {
      dist -= step ;
      break ;
    }
  }
  return(dist) ;
}
#endif
static double
find_max_distance(MRI_SURFACE *mris, MRI *mri, double x, double y, double z,
                  double nx, double ny, double nz, int grad_sign, double max_dist,
                  double sigma, double step, double *pval, 
                  double low_thresh, double hi_thresh, double max_inside_dist)
{
  double x1, y1, z1, grad, dist, val, next_val, xv, yv, zv, inside_dist = 0 ;
  int    inside = 0 ;  // don't let profile go into the desired region then back out

  for (dist = 0 ; dist <= max_dist ; dist += step)
  {
    x1 = x + dist*nx ; y1 = y + dist*ny ; z1 = z + dist*nz ;
    grad = MRISnormalIntensityGradient(mris, mri,  x1,  y1,  z1, nx,  ny,  nz, sigma);
    if (grad*grad_sign < 0)
      break ;
    MRISsurfaceRASToVoxelCached(mris, mri, x1, y1, z1, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &val) ;
    if (val < low_thresh || val > hi_thresh)  // outside of allowable interior intensity range
    {
      if (inside)   // don't go back outside if were previously in the interior of the region
      {
        dist -= step ;
        x1 = x + dist*nx ; y1 = y + dist*ny ; z1 = z + dist*nz ;
        MRISsurfaceRASToVoxelCached(mris, mri, x1, y1, z1, &xv, &yv, &zv) ;
        MRIsampleVolume(mri, xv, yv, zv, &val) ;
        if (pval)
          *pval = val ;
        return(dist) ;
      }
    }
    else  // in the right range to be in the interior of the region
    {
      if (val < (low_thresh+hi_thresh)/2)
        inside = 1 ;  // found the right intensity range
      inside_dist += step ;
      if (inside_dist >= max_inside_dist)
        return(dist-step) ;
    }
  }

  // now compute exact location on by looking for the dist where the intensities decrease
  dist -= step/2 ;
  MRISsurfaceRASToVoxelCached(mris, mri, x+dist*nx, y+dist*ny, z+dist*nz, &xv, &yv, &zv) ;
  dist += step/2 ;
  MRIsampleVolume(mri, xv, yv, zv, &val) ;
  for ( ; dist <= max_dist ; dist += step/2)
  {
    MRISsurfaceRASToVoxelCached(mris, mri, x+dist*nx, y+dist*ny, z+dist*nz, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &next_val) ;
    if ((next_val-val)*grad_sign < 0 || (val < low_thresh || val > hi_thresh))
      break ;
    val = next_val ;
  }
  if (pval)
    *pval = val ;
  return(dist) ;
}

static int
find_next_peak(double *intensity_profile, int index, double *pmax_val, int len) 
{
  int i, found, j ;

  for (i = index ; i < len-3 ; i++)
  {
    for (found = 1, j = MAX(0, i-3) ;  j <= MIN(len-1,i+3) ; j++)
    {
      if (i == j)
        continue ;
      if (intensity_profile[j] >= intensity_profile[i]) 
      {
        found = 0 ;
        break ;
      }
    }
    if (found)
    {
      *pmax_val = intensity_profile[i] ;
      return(i) ;
    }
  }
  return(-1) ;
}
static int
find_previous_peak(double *intensity_profile, int index, double *pmax_val, int len) 
{
  int i, found, j ;

  for (i = index ; i >= 0 ; i++)
  {
    for (found = 1, j = MAX(0, i-3) ;  j <= MIN(len-1,i+3) ; j++)
    {
      if (i == j)
        continue ;
      if (intensity_profile[j] >= intensity_profile[i]) 
      {
        found = 0 ;
        break ;
      }
    }
    if (found)
    {
      *pmax_val = intensity_profile[i] ;
      return(i) ;
    }
  }
  return(-1) ;
}
static int
find_next_valley(double *intensity_profile, int index, double *pmin_val, int len) 
{
  int i, found, j ;

  for (i = index ; i < len-3 ; i++)
  {
    for (found = 1, j = MAX(0, i-3) ;  j <= MIN(len-1,i+3) ; j++)
    {
      if (i == j)
        continue ;
      if (intensity_profile[j] <= intensity_profile[i]) 
      {
        found = 0 ;
        break ;
      }
    }
    if (found)
    {
      *pmin_val = intensity_profile[i] ;
      return(i) ;
    }
  }
  return(-1) ;
}
static int
mark_sg_index_okay(double *intensity_profile, int *sg_index_okay, double step, int len)
{
  int    index, i, j ;
  double valley_val, peak_val, previous_peak_val ;
  int    valley_index, peak_index, previous_peak ;
  
  for (index = 0 ; index < len ; index++)
  {
    sg_index_okay[index] = 0 ;  // assume it's on an upwards slope unless violated below
    for (i = index ; i <= index+1 && i < len ; i++)
    {
      for (j = i+1 ; j <= index+3 && j < len ; j++)
        if (intensity_profile[j] < intensity_profile[i])
        {
          sg_index_okay[index] = 1 ;
          break ;
        }
      if (sg_index_okay[index])
        break ;  // already found a violation of monotonicity
    }

    if (sg_index_okay[index])  
    {
      valley_index = find_next_valley(intensity_profile, index, &valley_val, len) ;
      if (valley_index < 0)
        continue ;

      if (valley_val > intensity_profile[i]-.15*intensity_profile[i]) // not a deep valley
      {
        previous_peak = find_previous_peak(intensity_profile,valley_index, &previous_peak_val,len);
        peak_index = find_next_peak(intensity_profile, valley_index, &peak_val, len) ;
        if (peak_index >= 0 && previous_peak >= 0 &&
            peak_val > (1+MIN_INTENISTY_PCT_OFFSET)*previous_peak_val)
          sg_index_okay[index] = 0 ;  // if the next peak is a lot higher than this one
      }
    }
    else  // on an upwards slope - mark everthing on the slope no good except peak
    {
      peak_index = find_next_peak(intensity_profile, index, &peak_val, len) ;
      for ( ; index < peak_index ; index++)
        sg_index_okay[index] = 0 ;
    }
  }
  return(NO_ERROR);
} 

static int
mark_wm_index_okay(double *intensity_profile, int *wm_index_okay, 
                   int wm_len, double min_wm_intensity, double max_wm_intensity,
                   int current_index, int len)
{
  int    index, inside, wm_samples, wm_start ;

  memset(wm_index_okay, 0, len*sizeof(int)) ;
  for (inside = wm_samples = 0, index = current_index-(wm_len-1) ; index < len ; index++)
  {
    if (intensity_profile[index] <= max_wm_intensity && intensity_profile[index] >= min_wm_intensity)
    {
      if (wm_samples++ >= wm_len)   // found interior
      {
        if (inside == 0)
          wm_start = index-wm_samples+1 ;
        inside = 1 ;
      }
    }
    else   // not in wm intensity range
    {
      if (inside) // everything else after here is invalid
        break ;
      inside = wm_samples = 0 ;
    }
  }
  if (inside)
  {
    for (index = wm_start+wm_samples ; index < len ; index++)
      wm_index_okay[index] = 0 ;
    for (index = current_index ; index < wm_start + wm_samples ; index++)
      wm_index_okay[index] = 1 ;
  }

  for (inside = wm_samples = 0, index = current_index+wm_len-1 ; index >= 0 ; index--)
  {
    if (intensity_profile[index] <= max_wm_intensity && intensity_profile[index] >= min_wm_intensity)
    {
      if (wm_samples++ >= wm_len)   // found interior
      {
        if (inside == 0)
          wm_start = index+wm_samples-1 ;
        inside = 1 ;
      }
    }
    else   // not in wm intensity range
    {
      if (inside) // everything else after here is invalid
        break ;
      inside = wm_samples = 0 ;
    }
  }
  if (inside)
  {
    for (index = wm_start - wm_samples ; index >= 0 ; index--)
      wm_index_okay[index] = 0 ;
    for (index = current_index ; index >= 0 && index > wm_start - wm_samples ; index--)
      wm_index_okay[index] = 1 ;
  }

  return(NO_ERROR);
} 

#define MAX_PROFILE_LEN 1000
#define DIST_NOT_FOUND  -10000
/*
  find the optimal profile location and fill in the vp->w[xyz] and vp->p[xyz]
  for the predicted pial and white locations.
*/
static double
find_optimal_locations(MRI_SURFACE *mris, MRI *mri, int vno, 
                       double max_inwards, double max_outwards, DP *dp, VP *vp, int skip)

{
  int     start_index, current_index, min_gray_matter_len, max_wm_len,
    max_start_index, csf_len, max_csf_len, max_profile_len, profile_len, wm_len,
    best_wm_len, best_csf_len, best_start_index, sg_index_okay[MAX_PROFILE_LEN], index,
    wm_index_okay[MAX_PROFILE_LEN] ;
  double errors[MAX_PROFILE_LEN], best_wm_dist, prior,
    kernel[MAX_PROFILE_LEN], nx, ny, nz, base_intensity_profile[MAX_PROFILE_LEN], wm_dist ;
  double  wm_val, ig_val, sg_val, *intensity_profile ;
  double  step, rms, min_rms, ig_width, sg_width, offset ;
  int     ig_len, sg_len, best_ig_len=0, best_sg_len=0, min_ig_len, max_ig_len,
    min_sg_len, max_sg_len, best_profile = PROFILE_GENERIC, tsteps, min_start_index ;
  double  best_wm_intensity_offset=0.0,min_wm_intensity_offset,max_wm_intensity_offset,
    wm_intensity_offset, min_wm_intensity, max_wm_intensity,
    best_sg_intensity_offset, sg_intensity_offset, ig_intensity_offset,  best_ig_intensity_offset ;
  
  int error_type = dp->error_type ;
  double (*profile_error_func)(double *kernel, double *intensity_profile, int nsamples, char *fname,
                               double step, double in_dist, double *errors);
  char   *fname1, *fname2 ;
  double  stria_val, stria_intensity_offset, best_stria_intensity_offset = 0.0 ; 

  intensity_profile = base_intensity_profile ;
  fname1 = (vno == Gdiag_no || deform_debug) ? "profile.dat" : NULL;
  fname2 = (vno == Gdiag_no || deform_debug) ? "intensity.dat" : NULL;

  vp->min_rms = min_rms = 1e10 ;
  dp->step = step = mri->xsize/2 ;
  tsteps = (int)ceil(0.25/step) ;
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
    min_wm_intensity_offset = -1.5*dp->max_wm_intensity_offset ;
    max_wm_intensity_offset = dp->max_wm_intensity_offset ;
    break ;
  case ERROR_FUNC_NORM_CORR: 
  case ERROR_FUNC_NCORR_L1:  
    min_wm_intensity_offset = max_wm_intensity_offset = 0.0 ;
    break ;
  }

  nx = vp->nx ; ny = vp->ny ; nz = vp->nz ;

  if (dp->fix_wm_intensity_offset)
  {
    max_wm_intensity_offset = min_wm_intensity_offset = vp->wm_intensity_offset ;
  }

  min_wm_intensity = dp->wm_val + min_wm_intensity_offset ;
  max_wm_intensity = dp->wm_val + max_wm_intensity_offset ;

  //  max_outwards += sqrt(SQR(vp->px-vp->wx) + SQR(vp->py-vp->wy) + SQR(vp->pz-vp->wz)) ;
  extract_image_intensities(mris, mri, max_inwards, max_outwards, vp->wx, vp->wy, vp->wz, nx, ny, nz, 
                            step, tsteps, base_intensity_profile, fname2);
  ig_width = (dp->max_ig_width + dp->min_ig_width) / 2.0 ;
  sg_width = (dp->max_sg_width + dp->min_sg_width) / 2.0 ;

  current_index = nint(max_inwards/step) ;
  max_profile_len = nint((max_inwards + max_outwards)/step) ;
  best_wm_dist = DIST_NOT_FOUND ;
  min_ig_len = nint(dp->min_ig_width/step) ; max_ig_len = nint(dp->max_ig_width/step) ;
  min_sg_len = nint(dp->min_sg_width/step) ; max_sg_len = nint(dp->max_sg_width/step) ;

  min_gray_matter_len = min_ig_len + min_sg_len ;
  max_start_index = max_profile_len - (min_gray_matter_len+1) ; // gm plus 1 CSF sample

  best_wm_len = best_ig_len = best_sg_len = best_start_index = best_csf_len = 0 ;
  best_ig_intensity_offset = best_sg_intensity_offset = best_wm_intensity_offset = 0.0 ;
  for (min_start_index = 0 ; min_start_index <= max_start_index ; min_start_index++)
    if (base_intensity_profile[min_start_index] >= min_wm_intensity &&
        base_intensity_profile[min_start_index] <= max_wm_intensity)
      break ;

  // just to make things run a bit faster
#define MAX_WM_LEN nint(5.0/step)
#define MIN_WM_LEN nint(1.0/step)
#define MAX_CSF_LEN nint(2.0/step)

  mark_sg_index_okay(base_intensity_profile, sg_index_okay, step, max_profile_len) ;
  mark_wm_index_okay(base_intensity_profile, wm_index_okay, nint(.5/step), 
                     min_wm_intensity, max_wm_intensity, current_index, max_profile_len);
  for (start_index = min_start_index ; start_index <= max_start_index; start_index += skip)
  {
    offset = (current_index - start_index)*step ;
    intensity_profile = base_intensity_profile + start_index ;
    max_wm_len = max_profile_len-(min_gray_matter_len+1+start_index) ;
    max_wm_len = MIN(max_wm_len, MAX_WM_LEN) ;
    for (wm_len = MIN_WM_LEN ; wm_len <= max_wm_len; wm_len += skip)
    {
      index = wm_len+start_index ;
      if (wm_index_okay[index] == 0)
        continue ;
      wm_dist = (start_index+wm_len-current_index)*step ;   
      wm_val = compute_optimal_intensity(intensity_profile, 0, wm_len-1) ;
      if (wm_val < min_wm_intensity || wm_val > max_wm_intensity)
        continue ; // not in feasible region
      wm_intensity_offset = wm_val - dp->wm_val ;
      for (ig_len = min_ig_len ; ig_len <= max_ig_len ; ig_len += skip)
      {
        ig_val = compute_optimal_intensity(intensity_profile, wm_len, wm_len+ig_len-1) ;
        if (ig_val < dp->infra_granular_val-dp->max_ig_intensity_offset ||
            ig_val > dp->infra_granular_val+dp->max_ig_intensity_offset)
          continue ; // not in feasible region
        ig_intensity_offset = ig_val - dp->infra_granular_val ;

        if (ig_val <= wm_val+MIN_INTENISTY_PCT_OFFSET*wm_val) // do some sanity checks
          continue ;
            
        for (sg_len = min_sg_len ; sg_len <= max_sg_len ; sg_len += skip)
        {
          index = wm_len+ig_len+sg_len+start_index ;
          if (sg_index_okay[index] == 0)
            continue ;
          if (sg_len > 3*ig_len || ig_len > 3*sg_len)
            continue ; // sanity check
          sg_val = compute_optimal_intensity(intensity_profile, wm_len+ig_len,wm_len+ig_len+sg_len-1);
          if (sg_val < dp->supra_granular_val-dp->max_sg_intensity_offset 
              /*  || sg_val > dp->supra_granular_val+dp->max_sg_intensity_offset*/)
            continue ; // not in feasible region
          sg_intensity_offset = sg_val - dp->supra_granular_val ;
            
          if (sg_val <= ig_val+MIN_INTENISTY_PCT_OFFSET*ig_val) // do some sanity checks
            continue ;
            
          max_csf_len = max_profile_len-(ig_len + sg_len + wm_len + 1) ;
          if (max_csf_len <= 0)
            continue ;  // enforce the presence of at least a bit of csf
          max_csf_len = MIN(max_csf_len, MAX_CSF_LEN) ;

          for (csf_len = 1 ; csf_len <= max_csf_len; csf_len += skip)
          {
            // compute length of subset of total profile being used
            profile_len = wm_len+ig_len+sg_len+csf_len ;
            
            construct_model_profile(kernel, dp, 
                                    wm_intensity_offset, ig_intensity_offset,sg_intensity_offset,
                                    0.0, wm_len, ig_len, sg_len, csf_len, PROFILE_GENERIC) ;
            
            prior = exp(-SQR(wm_dist)/(2*SQR(dp->sigma))) ;
            rms = (*profile_error_func)(kernel, intensity_profile, profile_len, NULL,step,
                                        offset, errors) ;
            if (DZERO(prior))
              rms += 20 ;
            else
              rms -= log(prior) ;
            prior = exp(-SQR(wm_intensity_offset)/(2*SQR(dp->intensity_sigma))) ;
            if (DZERO(prior))
              rms += 20 ;
            else
              rms -= log(prior) ;
            if (rms < min_rms)
            {
              if (vno == Gdiag_no || deform_debug)
              {
                (*profile_error_func)(kernel,intensity_profile,profile_len,"profile.dat",step, offset,
                                      errors);
              }
              best_wm_dist = wm_dist ;
              best_wm_len = wm_len ;
              best_start_index = start_index ;
              best_profile = PROFILE_GENERIC ;
              best_csf_len = csf_len ;
              best_wm_intensity_offset = wm_intensity_offset ;
              best_sg_intensity_offset = sg_intensity_offset ;
              best_ig_intensity_offset = ig_intensity_offset ;
              min_rms = rms ;
              best_ig_len = ig_len ; 
              best_sg_len = sg_len ; 
            }
          }
        }
      }
    }
  }
  intensity_profile = base_intensity_profile ;

  if (dp->try_stria_model)
  {
    int stria_start, stria_end, stria_len = nint(dp->stria_width/dp->step), 
        stria_offset_len = nint(STRIA_OFFSET/dp->step) ;  // amount of layerIV superficial to stria

    for (start_index = min_start_index ; start_index <= max_start_index; start_index += skip)
    {
      offset = (current_index - start_index)*step ;
      intensity_profile = base_intensity_profile + start_index ;
      max_wm_len = max_profile_len-(min_gray_matter_len+1+start_index) ;
      max_wm_len = MIN(max_wm_len, MAX_WM_LEN) ;
      for (wm_len = MIN_WM_LEN ; wm_len <= max_wm_len; wm_len += skip)
      {
        index = wm_len+start_index ;
        if (wm_index_okay[index] == 0)
          continue ;
        wm_dist = (start_index+wm_len-current_index)*step ;   
        wm_val = compute_optimal_intensity(intensity_profile, 0, wm_len-1) ;
        if (wm_val < min_wm_intensity || wm_val > max_wm_intensity)
          continue ; // not in feasible region
        wm_intensity_offset = wm_val - dp->wm_val ;
        for (ig_len = min_ig_len ; ig_len <= max_ig_len ; ig_len += skip)
        {
          ig_val = compute_optimal_intensity(intensity_profile, wm_len, wm_len+ig_len-1) ;
          if (ig_val < dp->infra_granular_val-dp->max_ig_intensity_offset ||
              ig_val > dp->infra_granular_val+dp->max_ig_intensity_offset)
            continue ; // not in feasible region
          ig_intensity_offset = ig_val - dp->infra_granular_val ;
          
          if (ig_val <= wm_val+MIN_INTENISTY_PCT_OFFSET*wm_val) // do some sanity checks
            continue ;
          
          stria_start = wm_len+ig_len-(stria_len+stria_offset_len) ;
          stria_end = stria_start + stria_len - 1  ;
          stria_val = compute_optimal_intensity(intensity_profile, stria_start, stria_end);
          stria_intensity_offset = stria_val - dp->stria_val ;
          
          // do some sanity checks
          if (ig_val <= stria_val)
            continue ;
          if (stria_val < wm_val)
            continue ;
          for (sg_len = min_sg_len ; sg_len <= max_sg_len ; sg_len += skip)
          {
            if (sg_len > 3*ig_len || ig_len > 3*sg_len)
              continue ; // sanity check
            index = wm_len+ig_len+sg_len+start_index ;
            if (sg_index_okay[index] == 0)
              continue ;
            sg_val = 
              compute_optimal_intensity(intensity_profile, wm_len+ig_len,wm_len+ig_len+sg_len-1);
            if (sg_val < dp->supra_granular_val-dp->max_sg_intensity_offset ||
                sg_val > dp->supra_granular_val+dp->max_sg_intensity_offset)
              continue ; // not in feasible region
            sg_intensity_offset = sg_val - dp->supra_granular_val ;
            
            if (sg_val <= ig_val+MIN_INTENISTY_PCT_OFFSET*ig_val) // do some sanity checks
              continue ;
            
            max_csf_len = max_profile_len-(ig_len + sg_len + wm_len + 1) ;
            if (max_csf_len <= 0)
              continue ;  // enforce the presence of at least a bit of csf
            max_csf_len = MIN(max_csf_len, MAX_CSF_LEN) ;
            
            for (csf_len = 1 ; csf_len <= max_csf_len; csf_len += skip)
            {
              // compute length of subset of total profile being used
              profile_len = wm_len+ig_len+sg_len+csf_len ;
              
              construct_model_profile(kernel, dp, 
                                      wm_intensity_offset, ig_intensity_offset,sg_intensity_offset,
                                      stria_intensity_offset, wm_len, ig_len, sg_len, csf_len, 
                                      PROFILE_V1) ;
              
              prior = exp(-SQR(wm_dist)/(2*SQR(dp->sigma))) ;
              rms = (*profile_error_func)(kernel, intensity_profile, profile_len, NULL,step,
                                          offset, errors) ;
              if (DZERO(prior))
                rms += 20 ;
              else
                rms -= log(prior) ;
              prior = exp(-SQR(wm_intensity_offset)/(2*SQR(dp->intensity_sigma))) ;
              if (DZERO(prior))
                rms += 20 ;
              else
                rms -= log(prior) ;
              if (rms < min_rms)
              {
                if (vno == Gdiag_no || deform_debug)
                {
                  (*profile_error_func)(kernel,intensity_profile,profile_len,"profile.dat",step, 
                                        offset, errors);
                }
                best_profile = PROFILE_V1 ;
                best_wm_dist = wm_dist ;
                best_wm_len = wm_len ;
                best_start_index = start_index ;
                best_csf_len = csf_len ;
                best_wm_intensity_offset = wm_intensity_offset ;
                best_stria_intensity_offset = stria_intensity_offset ;
                best_sg_intensity_offset = sg_intensity_offset ;
                best_ig_intensity_offset = ig_intensity_offset ;
                min_rms = rms ;
                best_ig_len = ig_len ; 
                best_sg_len = sg_len ; 
              }
            }
          }
        }
      }
    }
  }
  
  vp->wtx = vp->wx + best_wm_dist * nx ;
  vp->wty = vp->wy + best_wm_dist * ny ;
  vp->wtz = vp->wz + best_wm_dist * nz ;

  offset = (best_ig_len)*step ; 
  vp->l4tx = vp->wx + (best_wm_dist+offset) * nx ;
  vp->l4ty = vp->wy + (best_wm_dist+offset) * ny ;
  vp->l4tz = vp->wz + (best_wm_dist+offset) * nz ;

  offset = (best_ig_len + best_sg_len)*step ; 
  vp->ptx = vp->wx + (best_wm_dist+offset) * nx ;
  vp->pty = vp->wy + (best_wm_dist+offset) * ny ;
  vp->ptz = vp->wz + (best_wm_dist+offset) * nz ;

  vp->l4_dist = (best_wm_dist)+(best_ig_len*step) ;
  vp->white_dist = best_wm_dist ;
  vp->pial_dist = (best_wm_dist)+((best_sg_len+best_ig_len)*step) ;
  
  vp->wm_intensity_offset = best_wm_intensity_offset ;
  vp->sg_intensity_offset = best_sg_intensity_offset ;
  vp->ig_intensity_offset = best_ig_intensity_offset ;
  vp->which_profile = best_profile ;
  vp->min_rms = min_rms ;
  if (min_rms > .2*dp->wm_val && best_wm_dist > -100)
    best_wm_dist = DIST_NOT_FOUND ;
  if ((deform_debug || vno == Gdiag_no) && best_wm_dist > -100)
  {
    double rms ;
    offset = (current_index - best_start_index)*step ;
    profile_len = best_wm_len + best_ig_len + best_sg_len + best_csf_len ;
    intensity_profile = base_intensity_profile + best_start_index ;
    construct_model_profile(kernel, dp, best_wm_intensity_offset, best_ig_intensity_offset, 
                            best_sg_intensity_offset, best_stria_intensity_offset, 
                            best_wm_len, best_ig_len, best_sg_len, best_csf_len, best_profile) ;
    rms = (*profile_error_func)(kernel, intensity_profile,profile_len, "profile.dat",step,offset,
                                errors);
    prior = exp(-SQR(best_wm_dist)/(2*SQR(dp->sigma))) ;
    if (DZERO(prior))
      rms += 20 ;
    else
      rms -= log(prior) ;

    prior = exp(-(SQR(best_wm_intensity_offset)/(2*SQR(dp->intensity_sigma)))) ;
    if (DZERO(prior))
      rms += 20 ;
    else
      rms -= log(prior) ;
    printf("v %d: best_wm_dist = %2.1f, rms = %2.1f\n", vno, best_wm_dist, min_rms) ;
    DiagBreak() ;  // YYY
  }
  return(best_wm_dist) ;
}

static double
compute_optimal_intensity(double *iprofile, int i0, int i1)
{
  double optimal_intensity, intensities[MAX_PROFILE_LEN] ;
  int    num, i ;

  for (num = 0, i = i0 ; i <= i1 ; i++, num++)
    intensities[num] = iprofile[i] ;
  qsort(intensities, num, sizeof(double), compare_doubles) ;
  optimal_intensity = intensities[num/2] ;
  return(optimal_intensity) ;
}

static double
compute_profile_L1(double *kernel, double *iprofile, int nsamples, 
                   char *fname, double step, double in_dist, double *errors)
{
  double   val, kval, L1 ;
  int      i ;
  FILE     *fp = NULL ;

  if (fname)
    fp = fopen(fname, "w") ;

  for (L1 = 0.0, i = 0 ; i < nsamples ; i++)
  {
    kval = *kernel++ ;
    val = *iprofile++ ;
    if (fp)
      fprintf(fp, "%f %f %f\n", (float)i*step - in_dist, kval, val) ;
    val -= kval ; 
    L1 += fabs(val) ;
  }
  if (fp)
    fclose(fp) ;
  return(L1/nsamples) ;
}

static double
compute_profile_L2(double *kernel, double *iprofile, int nsamples, 
                   char *fname, double step, double in_dist, double *errors)
{
  double   rms, val, kval ;
  int      i ;
  FILE     *fp = NULL ;

  if (fname)
    fp = fopen(fname, "w") ;

  for (rms = 0.0, i = 0 ; i < nsamples ; i++)
  {
    kval = *kernel++ ;
    val = *iprofile++ ;
    if (fp)
      fprintf(fp, "%f %f %f\n", (float)i*step - in_dist, kval, val) ;
    val -= kval ; 
    rms += val*val ;
  }
  if (fp)
    fclose(fp) ;
  return(sqrt(rms/nsamples)) ;
}

static double
compute_profile_norm_corr(double *kernel, double *iprofile, int nsamples, char *fname,
                          double step, double in_dist, double *errors)
{
  double   *k, ival, kval, kmean, imean, kstd, istd ;
  int      i, nzero ;
  FILE     *fp = NULL ;
  double   ncorr ;

  if (fname)
    fp = fopen(fname, "w") ;

  k = kernel ; 
  kmean = imean = kstd = istd = 0.0 ;
  for (nzero = i = 0 ; i < nsamples ; i++)
  {
    kval = *kernel++ ;
    ival = iprofile[i] ;
    if (FZERO(ival))
      nzero++ ;
    if (fp)
      fprintf(fp, "%f %f %f\n", (float)i*step - in_dist, kval, ival) ;
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
  for (ncorr = 0.0, i = 0 ; i < nsamples ; i++)
  {
    kval = *kernel++ ;
    ival = *iprofile++ ;
    ncorr += (kval-kmean) * (ival-imean);
  }
  if (FZERO(kstd))
    kstd = 1e-4 ;
  if (FZERO(istd))
    istd = 1e-4 ;
  ncorr = ncorr / ((nsamples-1)*kstd*istd) ;
  return(-ncorr) ;
}

static int use_partial_volume_model = 1 ;
static int
construct_model_profile(double *kernel, DP *dp, double wm_intensity_offset, double ig_intensity_offset, 
                        double sg_intensity_offset,
                        double stria_intensity_offset, int inside_len, int infra_granular_len, 
                        int supra_granular_len, int out_len, int which_profile)
{
  int i, itotal, index, profile_len ;
  double wm_val, ig_val, sg_val, stria_val, out_val ;

  profile_len = inside_len + infra_granular_len + supra_granular_len + out_len ;
  wm_val = dp->wm_val+wm_intensity_offset ;
  ig_val = dp->infra_granular_val+ig_intensity_offset ;
  stria_val = dp->stria_val+stria_intensity_offset ;
  sg_val = dp->supra_granular_val+sg_intensity_offset ;
  out_val = dp->outside_val+wm_intensity_offset ;
  for (i = 0 ; i < inside_len ; i++)
    kernel[i] = wm_val ;
  for (itotal = i, i = 0 ; i < infra_granular_len ; i++, itotal++)
    kernel[itotal] = ig_val ;
  for (i = 0 ; i < supra_granular_len ; i++, itotal++)
    kernel[itotal] = sg_val ;
  for (i = 0 ; i < out_len ; i++, itotal++)
    kernel[itotal] = out_val ;

  for (i = inside_len-1 ; i >= 0 && i>=inside_len-4 ; i--)
    kernel[i] -= 0 ;  // account for dark band (u-fibers?)  was -10
  switch (which_profile)
  {
  default:
  case PROFILE_GENERIC:
  if (use_partial_volume_model)
  {
    index = inside_len+infra_granular_len ;
    if (index > 0)
      kernel[index-1]   = .75*ig_val + .25*sg_val;
    kernel[index]   = .5*ig_val + .5*sg_val;
    if (index < profile_len-1)
      kernel[index+1]   = .25*ig_val + .75*sg_val;
  }
    break ; 
  case PROFILE_V1:    // add the stria in
    {
      int stria_start, stria_len = nint(dp->stria_width/dp->step), 
        stria_offset_len = nint(STRIA_OFFSET/dp->step) ;  // amount of layerIV superficial to stria

      itotal = inside_len+infra_granular_len-(stria_len+stria_offset_len) ;
      for (i = 0 ; i < stria_len ; i++, itotal++)
        kernel[itotal] = stria_val ;
      if (use_partial_volume_model)
      {
        stria_start = inside_len+infra_granular_len-(stria_len+stria_offset_len) ;
        if (stria_start > 0)
          kernel[stria_start-1] = .25*stria_val + .75*ig_val;
        kernel[stria_start] = .5*stria_val + .5*ig_val;
        if (stria_start < profile_len-1)
          kernel[stria_start+1] = .75*stria_val + .25*ig_val;

        index = stria_start+stria_len  ;  // end of the stria
        kernel[index]   = .5*stria_val + .5*sg_val ;
        if (index > 0)
          kernel[index-1]   = .75*stria_val + .25*sg_val ;
        if (index < profile_len-1)
          kernel[index+1]   = .25*stria_val + .75*sg_val ;
      }
    }
    break ;
  }
  itotal = inside_len+infra_granular_len+supra_granular_len ;

  if (use_partial_volume_model)
  {
    kernel[inside_len] = .5*wm_val + .5*ig_val ;
    if (inside_len > 0)
      kernel[inside_len-1] = .75*wm_val + .25*ig_val ;
    if (inside_len < profile_len-1)
      kernel[inside_len+1] = .25*wm_val + .75*ig_val ;
    if (out_len > 0)
    {
      index = inside_len+infra_granular_len+supra_granular_len ;
      if (index >= 0 && index < profile_len)
        kernel[index] =   .5*sg_val + .5*out_val ;
      if (index > 0)
        kernel[index-1] =   .75*sg_val + .25*out_val ;
      if (index < profile_len-1)
        kernel[index+1] =   .25*sg_val + .75*out_val ;
    }
  }
  
  return(NO_ERROR) ;
}

static double
compute_linear_combination(double *kernel, double *iprofile, int nsamples, char *plot_fname,
                           double step, double in_dist, double *errors)
{
  double  ncorr, L1 ;

  ncorr = compute_profile_norm_corr(kernel, iprofile, nsamples, plot_fname, step, in_dist, errors);
  L1 = compute_profile_L1(kernel, iprofile, nsamples, plot_fname, step, in_dist, errors);

  return(ncorr*dp.ncorr_weight + L1) ;
}

static int
recompute_target_locations(MRI_SURFACE *mris, MRI *mri_white, MRI *mri_l4, MRI *mri_pial, DP *dp)
{
  int          vno ;
  VERTEX_PARMS *vp ;
  VERTEX       *v ;
  double       nx, ny, nz, norm ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    vp = (VERTEX_PARMS *)v->vp ;
    nx = (v->wnx + v->pnx)/2 ; ny = (v->wny + v->pny)/2 ; nz = (v->wnz + v->pnz)/2 ;
    norm = sqrt(nx*nx + ny*ny + nz*nz) ; 
    if (!FZERO(norm))
    {
      nx /= norm ; ny /= norm ; nz /= norm ;
    }
    else
      nx = 1.0 ;
    nx = vp->nx ; ny = vp->ny ; nz = vp->nz ;
    vp->wtx = vp->wx + vp->white_dist*nx ;
    vp->wty = vp->wy + vp->white_dist*ny ;
    vp->wtz = vp->wz + vp->white_dist*nz ;

    if (mri_white)
    {
      double xv, yv, zv ;

      MRISsurfaceRASToVoxelCached(mris, mri_white, vp->wtx, vp->wty, vp->wtz, &xv, &yv, &zv);
      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
        DiagBreak() ;
      if (MRIindexNotInVolume(mri_white, xv, yv, zv) == 0)
        MRIsetVoxVal(mri_white, nint(xv), nint(yv), nint(zv), 0, v->val) ;
      else
      {
        v->marked = 0 ;
      }
    }

    vp->ptx = vp->wx + vp->pial_dist*nx ;
    vp->pty = vp->wy + vp->pial_dist*ny ;
    vp->ptz = vp->wz + vp->pial_dist*nz ;
    if (mri_pial)
    {
      double xv, yv, zv ;

      MRISsurfaceRASToVoxelCached(mris, mri_pial, vp->ptx, vp->pty, vp->ptz, &xv, &yv, &zv);
      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
        DiagBreak() ;
      if (MRIindexNotInVolume(mri_pial, xv, yv, zv) == 0)
        MRIsetVoxVal(mri_pial, nint(xv), nint(yv), nint(zv), 0, v->val) ;
      else
      {
        v->marked = 0 ;
      }
    }

    vp->l4tx = vp->wx + vp->l4_dist*nx ;
    vp->l4ty = vp->wy + vp->l4_dist*ny ;
    vp->l4tz = vp->wz + vp->l4_dist*nz ;
    if (mri_l4)
    {
      double xv, yv, zv ;

      MRISsurfaceRASToVoxelCached(mris, mri_l4, vp->l4tx, vp->l4ty, vp->l4tz, &xv, &yv, &zv);
      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
        DiagBreak() ;
      if (MRIindexNotInVolume(mri_l4, xv, yv, zv) == 0)
        MRIsetVoxVal(mri_l4, nint(xv), nint(yv), nint(zv), 0, v->val) ;
      else
      {
        v->marked = 0 ;
      }
    }
  }
  return(NO_ERROR) ;
}
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
#if 0
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
#endif

#include "cma.h"   
static int
fill_aseg(MRI *mri, char *aseg_fname, DP *dp)
{
  MRI  *mri_aseg ;
  int  x, y, z, label ;

  mri_aseg = MRIread(aseg_fname) ;
  if (mri_aseg == NULL)
    ErrorExit(ERROR_NOFILE, "%s: couldn't read aseg volume from %s", Progname, aseg_fname) ;

  for (x = 0 ; x < mri_aseg->width ; x++)
    for (y = 0 ; y < mri_aseg->height ; y++)
      for (z = 0 ; z < mri_aseg->depth ; z++)
      {
        label = nint(MRIgetVoxVal(mri_aseg, x, y, z, 0)) ;
        if (label == Left_Lateral_Ventricle ||
            label == Right_Lateral_Ventricle)
          MRIsetVoxVal(mri, x, y, z, 0, dp->wm_val) ;
      }
  return(NO_ERROR) ;
}
#if 0
static double
sum_L1_abs_distances(double *dists, int num, double d)
{
  double  sum ;
  int     n ;

  for (sum = 0.0, n = 0 ; n < num ; n++)
    sum += fabs(dists[n]-d) ;
  return(sum) ;
}
#endif
#define MAX_DISTS 100
static int
is_outlier(MRI_SURFACE *mris, int vno, int which)
{
  int       n, num ;
  VERTEX    *vn, *v ;
  double    dists[MAX_DISTS], mn, std, min_dist, max_dist, dist ;
  HISTOGRAM *h ;
  VERTEX_PARMS *vp ;
#if 0
  int       imax ;
  double    md ;
#endif

  memset(dists, 0, sizeof(dists)) ;
  v = &mris->vertices[vno] ;

  vp = (VERTEX_PARMS *)(v->vp) ;
  dist = which == WHITE_VERTICES ? vp->white_dist : vp->pial_dist ; ;
  min_dist = max_dist = dist ;
  for (num = n = 0 ; n < v->vtotal ; n++)
  {
    vn = &mris->vertices[v->v[n]];
    if (vn->ripflag || vn->marked != 1)
      continue ;

    if (num >= MAX_DISTS)
      break ;
    dists[num] = which == WHITE_VERTICES ? vp->white_dist : vp->pial_dist ;
    if (dists[num] < min_dist)
      min_dist = dists[num];
    if (dists[num] > max_dist)
      max_dist = dists[num] ;
    num++ ;
  }
  if (FEQUAL(min_dist, max_dist))
    return(0) ;
  h = HISTObins(100, min_dist, max_dist) ;
  for (n = 0 ; n < num ; n++)
    HISTOaddSample(h, (float)dists[n], min_dist, max_dist) ;

  HISTOsoapBubbleZeros(h, h, 100) ;
  HISTOsmooth(h, h, 2) ;
  HISTOrobustGaussianFit(h, .75, &mn, &std) ;
  
  if (vno == Gdiag_no)
  {
    FILE *fp  ;
    VERTEX_PARMS *vnp ;
    fp = fopen("out.dat", "w") ;
    HISTOplot(h, "h.plt") ;
    for (n = 0 ; n < v->vtotal ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->marked != 1 || vn->ripflag)
        continue ;
      vnp = (VERTEX_PARMS *)vn->vp ;
      fprintf(fp, "%d %f\n", v->v[n], which == WHITE_VERTICES ? vnp->white_dist : vnp->pial_dist);
    }
    fclose(fp) ;
  }
  if (num <= 1)
    return(0) ;

#if 0
  qsort(dists, num, sizeof(double), compare_doubles) ;
  md = dists[num/2];
  for (n = 0 ; n < num ; n++)
    dists[n] -= md ; // remove the median
  qsort(dists, num, sizeof(double), compare_abs_doubles) ;
  if (vno == Gdiag_no)
  {
    FILE *fp  ;
    fp = fopen("dist.dat", "w") ;
    for (n = 0 ; n < num ; n++)
      fprintf(fp, "%f\n", dists[n]) ;
    fclose(fp) ;
  }

  imax = num-(int)(ceil(num/4)) ;  // exclude the top 10% when computing std
  if (imax <= 2)
    return(0) ;

  for (mn = std = 0.0, n = 0 ; n <= imax ; n++)
  {
    mn += fabs(dists[n]) ;
    std += dists[n]*dists[n] ;
  }
  num  = imax+1 ;
  std = sqrt((std - mn*mn/num)/(num-1)) ;
  mn /= num ;
  if (fabs(dist-md) > 4*mn)
    return(1) ;
#else
  if (fabs(dist-mn) > 8*std)
    return(1) ;
#endif
  HISTOfree(&h) ;
  return(0) ;
}


static int
vp_copy_to_surface(MRI_SURFACE *mris, int which_src, int which_dst)
{
  int          vno ;
  VERTEX       *v ;
  double       x, y, z; 
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    switch (which_src)
    {
    default:
    case WHITE_VERTICES: x = vp->wx ;  y = vp->wy ;  z = vp->wz ; break ;
    case PIAL_VERTICES:  x = vp->px ;  y = vp->py ;  z = vp->pz ; break ;
    case WHITE_TARGETS:  x = vp->wtx ; y = vp->wty ; z = vp->wtz ; break ;
    case LAYERIV_VERTICES: x = vp->l4x ; y = vp->l4y ; z = vp->l4z ; break;
    case PIAL_TARGETS:   x = vp->ptx ; y = vp->pty ; z = vp->ptz ; break ;
    case LAYERIV_TARGETS:x = vp->l4tx ; y = vp->l4ty ; z = vp->l4tz ; break ;
      break ;
    }

    switch (which_dst)
    {
    case WHITE_VERTICES:   v->whitex = x ; v->whitey = y; v->whitez = z ;break ;
    case PIAL_VERTICES:    v->pialx = x ;  v->pialy = y;  v->pialz = z ; break ;
    case ORIG_VERTICES:    v->origx = x ;  v->origy = y;  v->origz = z ; break ;
    case CURRENT_VERTICES: v->x = x ;      v->y = y;      v->z = z ;     break ;
    case TARGET_VERTICES:  v->targx = x ;  v->targy = y ; v->targz = z ; break ;
    default:
      break ;
    }
  }
  return(NO_ERROR) ;
}

static int
vp_copy_from_surface(MRI_SURFACE *mris, int which_src, int which_dst)
{
  int          vno ;
  VERTEX       *v ;
  double       x, y, z; 
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    switch (which_src)
    {
    default:
    case WHITE_VERTICES:   x = v->whitex  ; y = v->whitey ; z = v->whitez  ;break ;
    case PIAL_VERTICES:    x = v->pialx  ;  y = v->pialy ;  z = v->pialz  ; break ;
    case ORIG_VERTICES:    x = v->origx  ;  y = v->origy ;  z = v->origz  ; break ;
    case CURRENT_VERTICES: x = v->x  ;      y = v->y ;      z = v->z  ;     break ;
    case TARGET_VERTICES:  x = v->targx ;   y = v->targy ;  z = v->targz ;  break ;
      break ;
    }

    switch (which_dst)
    {
    case LAYERIV_VERTICES: vp->l4x = x ; vp->l4y = y ; vp->l4z = z ; break ;
    case WHITE_VERTICES:   vp->wx = x ;  vp->wy = y ;  vp->wz = z ; break ;
    case PIAL_VERTICES:    vp->px = x ;  vp->py = y ;  vp->pz = z ; break ;
    case WHITE_TARGETS:    vp->wtx = x ; vp->wty = y ; vp->wtz = z ; break ;
    case PIAL_TARGETS:     vp->ptx = x ; vp->pty = y ; vp->ptz = z ; break ;
    case LAYERIV_TARGETS:  vp->l4tx = x ; vp->l4ty = y ; vp->l4tz = z ; break ;
    default:
      break ;
    }
  }
  return(NO_ERROR) ;
}
static int
vp_copy_to_surface_vals(MRI_SURFACE *mris, int which)
{
  int          vno ;
  VERTEX       *v ;
  double       val; 
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    switch (which)
    {
    default:
    case WHITE_VERTICES:   val = vp->white_val ; break ;
    case PIAL_VERTICES:    val = vp->pial_val ; break ;
    case LAYERIV_VERTICES: val = vp->l4_val ; break ;
    case WM_INTENSITY_OFFSET: val = vp->wm_intensity_offset ; break ;
    case SG_INTENSITY_OFFSET: val = vp->sg_intensity_offset ; break ;
    case RMS:                 val = vp->min_rms ; break ;
    case IG_INTENSITY_OFFSET: val = vp->ig_intensity_offset ; break ;
    }
    v->marked = vp->found ;
    v->val = val ;
  }
  return(NO_ERROR) ;
}
static int
vp_copy_dist_to_surface_vals(MRI_SURFACE *mris, int which)
{
  int          vno ;
  VERTEX       *v ;
  double       dist; 
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    switch (which)
    {
    default:
    case WHITE_VERTICES:   dist = vp->white_dist ; break ;
    case PIAL_VERTICES:    dist = vp->pial_dist ; break ;
    case LAYERIV_VERTICES: dist = vp->l4_dist ; break ;
    }
    v->val = dist ;
  }
  return(NO_ERROR) ;
}
static int
copy_surface_vals_to_vp(MRI_SURFACE *mris, int which)
{
  int          vno ;
  VERTEX       *v ;
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    switch (which)
    {
    default:
    case WHITE_VERTICES:   vp->white_val = v->val ; break ;
    case PIAL_VERTICES:    vp->pial_val = v->val ; break ;
    case LAYERIV_VERTICES: vp->l4_val = v->val ; break ;
    case WM_INTENSITY_OFFSET: vp->wm_intensity_offset = v->val ; break ;
    case SG_INTENSITY_OFFSET: vp->sg_intensity_offset = v->val ; break ;
    case RMS:                 vp->min_rms = v->val ; break ;
    case IG_INTENSITY_OFFSET: vp->ig_intensity_offset = v->val ; break ;
    }
  }
  return(NO_ERROR) ;
}
static int
vp_copy_to_surface_dist(MRI_SURFACE *mris, int which)
{
  int          vno ;
  VERTEX       *v ;
  double       val; 
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    switch (which)
    {
    default:
    case WHITE_VERTICES:   val = vp->white_dist ; break ;
    case PIAL_VERTICES:    val = vp->pial_dist ; break ;
    case LAYERIV_VERTICES: val = vp->l4_dist ; break ;
    }
    v->marked = vp->found ;
    v->d = val ;
  }
  return(NO_ERROR) ;
}
static int
vp_copy_from_surface_dist(MRI_SURFACE *mris, int which)
{
  int          vno ;
  VERTEX       *v ;
  double       val; 
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    val = v->d ;
    switch (which)
    {
    default:
    case WHITE_VERTICES:   vp->white_dist = val ; break ;
    case PIAL_VERTICES:    vp->pial_dist = val ; break ;
    case LAYERIV_VERTICES: vp->l4_dist = val ; break ;
    }
  }
  return(NO_ERROR) ;
}
static int
vp_set(MRI_SURFACE *mris, int which, double val)
{
  int          vno ;
  VERTEX       *v ;
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    switch (which)
    {
    default:
    case WHITE_VERTICES:   vp->white_val = val ; break ;
    case PIAL_VERTICES:    vp->pial_val = val ; break ;
    case LAYERIV_VERTICES: vp->l4_val = val ; break ;
    case WM_INTENSITY_OFFSET: vp->wm_intensity_offset = val ; break ;
    case SG_INTENSITY_OFFSET: vp->sg_intensity_offset = val ; break ;
    case RMS:                 vp->min_rms = val ; break ;
    case IG_INTENSITY_OFFSET: vp->ig_intensity_offset = val ; break ;
    }
  }
  return(NO_ERROR) ;
}

static double
vp_mean(MRI_SURFACE *mris, int which)
{
  int          vno, num ;
  VERTEX       *v ;
  double       dist, mean; 
  VERTEX_PARMS *vp ;

  for (mean = 0.0, num = 0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    if (vp->found == 0)
      continue ;
    switch (which)
    {
    default:
    case WHITE_VERTICES:   dist = vp->white_dist ; break ;
    case PIAL_VERTICES:    dist = vp->pial_dist ; break ;
    case LAYERIV_VERTICES: dist = vp->l4_dist ; break ;
    }
    mean += fabs(dist) ;
    num++ ;
  }
  return(mean/num) ;
}

static LABEL *
vp_make_v1_label(MRI_SURFACE *mris, double thresh)
{
  int     vno ;
  VERTEX  *v ;
  VP      *vp ;
  LABEL   *v1 ;

  MRIScopyMarkedToMarked2(mris) ;
  MRISclearMarks(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    vp = (VERTEX_PARMS *)v->vp ;
    if (vp->which_profile != PROFILE_V1 && vp->pval_v1 < thresh)
      continue ;
    v->marked = 1 ;
    v->stat = vp->pval_v1 ;
  }

  v1 = LabelFromMarkedSurface(mris) ;
  MRIScopyMarked2ToMarked(mris) ;
  return(v1) ;
}

static LABEL *
label_v1(MRI_SURFACE *mris, MRI *mri, DP *dp)
{
  LABEL  *v1 ;
  int    vno ;
  VERTEX *v ;

  dp->step = mri->xsize/2 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    compute_normal(mris, vno, (VERTEX_PARMS *)(v->vp)) ;
    compute_best_neighborhood_profile(mris, mri, vno, dp) ; 
  }

  v1 = vp_make_v1_label(mris, v1_thresh) ;

  return(v1) ;
}

static int
compute_best_neighborhood_profile(MRI_SURFACE *mris, MRI *mri,int vno, DP *dp)
{
  VP     *vp ;
  int    best_profile, wm_len, ig_len, sg_len, i, profile_len, n, stria_len, num, tsteps ;
  double kernel[MAX_PROFILE_LEN], nx, ny, nz, wm_intensity_offset, sg_intensity_offset, in_dist, wm_val,
         intensity_profile[MAX_PROFILE_LEN],norm, ig_dist, sg_dist, v1_rms_total, 
    generic_rms_total, generic_rms, v1_rms ;
  VERTEX *v, *vn ;
  double (*profile_error_func)(double *kernel, double *intensity_profile, 
                               int nsamples, char *fname,double step, double in_dist,
                               double *errors);
  char   *fname, *fname2 ;

  tsteps = (int)ceil(0.25/dp->step) ;

  if (vno == Gdiag_no || deform_debug)
  {
    fname = "profile.dat" ;
    fname2 = "intensity.dat" ;
  }
  else
    fname = NULL ;

  stria_len = nint((dp->stria_width+STRIA_OFFSET) / dp->step) ;
  v = &mris->vertices[vno] ;
  switch (dp->error_type)
  {
  default:
  case ERROR_FUNC_L1:        profile_error_func = compute_profile_L1 ;        break ;
  case ERROR_FUNC_L2:        profile_error_func = compute_profile_L2 ;        break ;
  case ERROR_FUNC_NORM_CORR: profile_error_func = compute_profile_norm_corr ; break ;
  case ERROR_FUNC_NCORR_L1:  profile_error_func = compute_linear_combination; break ;
  }

  for (num = 0, v1_rms_total = generic_rms_total = 0.0, n = -1 ; n < 0 ; n++)
  {
    if (n < 0)
      vn = v ;
    else
      vn = &mris->vertices[v->v[n]] ;
    if (vn->ripflag)
      continue ;
    
    num++ ;
    vp = (VERTEX_PARMS *)(vn->vp) ;
    ig_dist = sqrt(SQR(vp->l4x-vp->wx) + SQR(vp->l4y-vp->wy) + SQR(vp->l4z-vp->wz)) ;
    sg_dist = sqrt(SQR(vp->px-vp->l4x) + SQR(vp->py-vp->l4y) + SQR(vp->pz-vp->l4z)) ;
    nx = vp->nx ; ny = vp->ny ; nz = vp->nz ;
    
    in_dist = find_max_distance(mris, mri, vp->wx, vp->wy, vp->wz, -nx, -ny, -nz,
                                dp->contrast_type == T1 ? 1 : -1, 2, dp->sigma, dp->step,
                                &wm_val, dp->wm_val-dp->max_wm_intensity_offset,
                                dp->wm_val+dp->max_wm_intensity_offset, dp->max_dist) ;
    wm_len = nint(in_dist / dp->step) ;
    ig_len = nint(ig_dist / dp->step) ;
    sg_len = nint(sg_dist / dp->step) ;
    profile_len = wm_len+ig_len+sg_len ;
    if (ig_len+sg_len <= stria_len)
      continue ;
    extract_image_intensities(mris, mri, in_dist, ig_dist+sg_dist, vp->wx, vp->wy, vp->wz, nx, ny, nz, 
                              dp->step, tsteps, intensity_profile,fname2) ;
    if (wm_len > 0)
    {
      for (wm_intensity_offset = 0.0, i = 0 ; i < wm_len ; i++)
        wm_intensity_offset += intensity_profile[i] ;
      wm_intensity_offset /= wm_len ;
      wm_intensity_offset -= dp->wm_val;
    }
    else
      wm_intensity_offset = wm_val - dp->wm_val ;
    
    for (i = wm_len, sg_intensity_offset = 0.0  ; i < wm_len+(ig_len-stria_len) ; i++)
      sg_intensity_offset += (intensity_profile[i] - (dp->infra_granular_val+wm_intensity_offset)) ;
    for (  ; i < wm_len+ig_len+sg_len ; i++)
      sg_intensity_offset += (intensity_profile[i] - (dp->supra_granular_val+wm_intensity_offset)) ;
    if (stria_len < ig_len)
      sg_intensity_offset /= (ig_len+sg_len-stria_len) ;
    else
      sg_intensity_offset /= (sg_len) ;

    generic_rms = find_min_rms(kernel, intensity_profile, profile_len, fname, dp, in_dist, wm_intensity_offset, 
                               sg_intensity_offset-dp->max_sg_intensity_offset, sg_intensity_offset+dp->max_sg_intensity_offset, 1, 
                               PROFILE_GENERIC, wm_len, ig_len, sg_len, profile_error_func);
    v1_rms = find_min_rms(kernel, intensity_profile, profile_len, fname, dp, in_dist, wm_intensity_offset, 
                               sg_intensity_offset-dp->max_sg_intensity_offset, sg_intensity_offset+dp->max_sg_intensity_offset, 1, 
                      PROFILE_V1, wm_len, ig_len, sg_len, profile_error_func);
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
      if (generic_rms < v1_rms)
        DiagBreak() ;
      else
        DiagBreak() ;
    }
    generic_rms_total += generic_rms ;
    v1_rms_total += v1_rms ;
  }

  if (vno == Gdiag_no)
    DiagBreak() ;
#define PSMALL 30
  if (dp->use_prior)
  {
    if (FZERO(v->stat))
      v1_rms_total += PSMALL ;
    else
      v1_rms_total += -log10(v->stat) ;
    if (FEQUAL(v->stat, 1.0))
      generic_rms_total += PSMALL ;
    else
      generic_rms_total += -log10(1.0-v->stat) ;
  }
  vp = (VERTEX_PARMS *)(v->vp) ;
  if (v1_rms_total < generic_rms_total)
    best_profile = PROFILE_V1 ;
  else
    best_profile = PROFILE_GENERIC ;
  vp->which_profile = best_profile ;
  if (num > 0)
  {
    vp->pval_v1 = exp(-v1_rms_total) ;
    norm = vp->pval_v1 + exp(-generic_rms_total) ;
    vp->pval_v1 = vp->pval_v1 / norm ;
  }
  return(best_profile) ;
}

static double
find_min_rms(double *kernel, double *intensity_profile, int profile_len, 
             char *fname, DP *dp, double in_dist, double wm_intensity_offset, 
             double min_sg_intensity_offset, double max_sg_offset, double sg_intensity_offset_step, 
             int which_profile, int wm_len, int ig_len, int sg_len,   
             double (*profile_error_func)(double *kernel, double *intensity_profile, 
                                          int nsamples, char *fname,double step, 
                                          double in_dist, double *errors))
{
  double    min_rms, sg_intensity_offset, best_sg_intensity_offset = 0, rms, ig_intensity_offset, best_ig_intensity_offset = 0.0,
    stria_intensity_offset, best_stria_intensity_offset = 0.0 ; 
  int       i, stria_len = nint((dp->stria_width+STRIA_OFFSET) / dp->step) ;

  min_rms = 1e8 ;
  for (stria_intensity_offset = -dp->max_stria_intensity_offset ; 
       stria_intensity_offset <= dp->max_stria_intensity_offset ; 
       stria_intensity_offset += (dp->max_stria_intensity_offset/2))
  for (ig_intensity_offset = 
         -dp->max_ig_intensity_offset ; 
       ig_intensity_offset <= dp->max_ig_intensity_offset ; 
       ig_intensity_offset += (dp->max_ig_intensity_offset/2))
  for (sg_intensity_offset = min_sg_intensity_offset ; 
       sg_intensity_offset <= max_sg_offset ; 
       sg_intensity_offset += sg_intensity_offset_step)
  {
    construct_model_profile(kernel, dp, wm_intensity_offset, ig_intensity_offset, 
                            sg_intensity_offset, stria_intensity_offset, wm_len, ig_len, 
                            sg_len, 0,which_profile) ; ;
#if 1
    rms = (*profile_error_func)(kernel, intensity_profile, profile_len, fname,
                                dp->step, in_dist, NULL) ;
    rms /= profile_len ;
#else
    for (rms = 0.0, i = wm_len+ig_len-stria_len ; i < wm_len+ig_len ; i++)
      rms += fabs(intensity_profile[i] - kernel[i]) ;
    rms /= stria_len ;
#endif
    if (rms < min_rms)
    {
      min_rms = rms ;
      best_sg_intensity_offset = sg_intensity_offset ;
      best_ig_intensity_offset = ig_intensity_offset ;
      best_stria_intensity_offset = stria_intensity_offset ;
    }
  }

  construct_model_profile(kernel, dp, wm_intensity_offset, best_ig_intensity_offset, 
                          best_sg_intensity_offset, best_stria_intensity_offset,
                          wm_len, ig_len, sg_len, 0, which_profile) ; ;
  for (rms = 0.0, i = wm_len+ig_len-(stria_len+1) ; i <= wm_len+ig_len ; i++)
    rms += fabs(intensity_profile[i] - kernel[i]) ;
  rms /= (stria_len+2) ;
  if (fname)
  {
    rms = (*profile_error_func)(kernel, intensity_profile, profile_len, fname,
                                dp->step, in_dist, NULL) ;
    DiagBreak()  ;
  }
  return(min_rms) ;
}

static int
compute_normal(MRI_SURFACE *mris, int vno, VERTEX_PARMS *vp)
{
  double nx, ny, nz, norm ;
  VERTEX *v ;

  v = &mris->vertices[vno] ;
  nx = vp->px - vp->wx ; ny = vp->py - vp->wy ; nz = vp->pz - vp->wz ;
  norm = sqrt(nx*nx + ny*ny + nz*nz) ; 
  if (FZERO(norm)) // they start in the same place
  {
    nx = (v->wnx + v->pnx)/2 ; ny = (v->wny + v->pny)/2 ; nz = (v->wnz + v->pnz)/2 ;
    norm = sqrt(nx*nx + ny*ny + nz*nz) ; 
    if (FZERO(norm))
      nx = 1.0 ;
  }
  nx /= norm ; ny /= norm ; nz /= norm ;
  vp->nx = nx ; vp->ny = ny ; vp->nz = nz ;
  return(NO_ERROR) ;
}
