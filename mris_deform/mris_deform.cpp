/**
 * @brief program for deforming a surface to lie at the gray/white or pial boundary from 
 *  ultra-high res data
 *
 * Fit a generative piecewise constant model to the data to determine
 * target locations and deform the surface to match them.
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
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "transform.h"
#include "mrisurf.h"
#include "label.h"
#include "tritri.h"
#include "filter.h"
#include "icosahedron.h"

#define MAX_PROFILE_LEN 1000
#define SMOOTH_LEN      3
#define PV_LEN          3
#define MIN_INTENISTY_PCT_OFFSET 0.15
#define STRIA_OFFSET .0
#define PROFILE_GENERIC   0
#define PROFILE_V1        1
#define PROFILE_AGRANULAR 2

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
  double deep_ratio ;       // ratio of bottom of IG to top of IG (good for stria)
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
  int    use_intensity ;  // position surfaces to intensity targets
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
  int    fix_intensities ;
  int    dark_csf ;  // for ex vivo scans
  double step ;
  int    try_stria_model ;
  int    use_prior ;
  double max_stria_intensity_offset ;
  int    search_multiple_angles ;
} DEFORMATION_PARMS, DP ;

#define WM_INTENSITY_OFFSET   100
#define SG_INTENSITY_OFFSET   101
#define IG_INTENSITY_OFFSET   102
#define RMS                   103
#define WM_INTENSITY          104
#define SG_INTENSITY          105
#define IG_INTENSITY          106
#define IG_WM_RATIO           107
#define SG_WM_RATIO           108
#define SG_IG_RATIO           109
#define DEEP_RATIO            110

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
#define SURFACE_NORMALS       103

static int pad_voxels = 2 ;
static float threshold = 0 ;
static int white_only = 0 ;

static int extract_image_intensities(MRI_SURFACE *mris, MRI *mri, double in_dist, 
                                     double out_dist, double x0, double y0, double z0,
                                     double nx, double ny, double nz, double step, 
                                     int tsteps,double *intensity_profile, const char *fname);
static int rip_dark_vertices(MRI_SURFACE *mris, MRI *mri, float threshold) ;
static int rip_vertices_with_no_faces(MRI_SURFACE *mris) ;
//static int filter_coords(MRI_SURFACE *mris, int filter_type)  ;
static int vp_copy_found_from_marked(MRI_SURFACE *mris) ;
//static int vp_copy_found_to_marked(MRI_SURFACE *mris) ;
static int dilate_not_found(MRI_SURFACE *mris, int niter)  ;
static int compute_laminar_ratios(MRI_SURFACE *mris, MRI *mri, DP *dp) ;
static int compute_distances_from_positions(MRI_SURFACE *mris) ;
static int filter_vals(MRI_SURFACE *mris, int filter_type)  ;
static int filter_distances(MRI_SURFACE *mris, int filter_type)  ;
//static int compute_normal(MRI_SURFACE *mris, int vno, VERTEX_PARMS *vp) ;
static int compute_normals(MRI_SURFACE *mris, VERTEX_PARMS *vp) ;
static int compare_doubles(const void *v1, const void *v2);
static double compute_optimal_intensity(double *iprofile, int i0, int i1, int filter_type) ;
static int compute_best_neighborhood_profile(MRI_SURFACE *mris, MRI *mri,
                                             int vno, 
                                             DP *dp) ;
static double find_min_rms(double *kernel, double *intensity_profile, int profile_len, 
                           const char *fname, DP *dp, double wm_dist, double wm_offset, 
                           double min_gm_offset, double max_gm_offset, double gm_offset_step, 
                           int which_profile, int wm_len, int ig_len, int sg_len,   
                           double (*profile_error_func)(double *kernel, double *intensity_profile, 
                                                        int nsamples, const char *fname,double step, 
                                                        double wm_dist, double *errors));
static LABEL *label_v1(MRI_SURFACE *mris, MRI *mri, DP *dp) ;
static LABEL *vp_make_v1_label(MRI_SURFACE *mris, double thresh) ;
static double vp_mean(MRI_SURFACE *mris, int which) ;
static int vp_set(MRI_SURFACE *mris, int which, double val) ;
static int vp_copy_to_surface_vals(MRI_SURFACE *mris, int which, DP *dp) ;
static int vp_copy_from_surface_vals(MRI_SURFACE *mris, int which, DP *dp) ;
static int vp_copy_dist_to_surface_vals(MRI_SURFACE *mris, int which) ;
static int vp_copy_to_surface_dist(MRI_SURFACE *mris, int which) ;
static int vp_copy_from_surface_dist(MRI_SURFACE *mris, int which) ;
static int vp_copy_from_surface(MRI_SURFACE *mris, int which_src, int which_dst) ;
static int vp_copy_to_surface(MRI_SURFACE *mris, int which_src, int which_dst) ;
static int is_outlier(MRI_SURFACE *mris, int vno, int which) ;
static int recompute_target_locations(MRI_SURFACE *mris, MRI *mri_white, MRI *mri_l4, 
                                      MRI *mri_pial, DP *dp) ;
static double find_optimal_locations(MRI_SURFACE *mris, MRI *mri, int vno,
                                     double max_inwards, double max_outwards, DP *dp,
                                     VP *vp, int skip);
static double compute_linear_combination(double *kernel, double *iprofile, int nsamples, 
                                         const char *plot_fname, double step, double wm_dist,
                                         double *errors);
static double compute_profile_L1(double *kernel, double *iprofile, int nsamples, 
                                 const char *plot_fname, double step, double wm_dist, double *errors);
static double compute_profile_L2(double *kernel, double *iprofile, int nsamples, 
                                 const char *plot_fname, double step, double wm_dist, double *errors);
static double compute_profile_norm_corr(double *kernel, double *iprofile, int nsamples, 
                                        const char *plot_fname, double step, double wm_dist, double *errors);
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
static int compute_intensity_offsets(MRI_SURFACE *mris, MRI *mri, MRI *mri_dist, DP *dp);

static INTEGRATION_PARMS parms ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int fill_aseg(MRI *mri, char *aseg_fname, DP *dp) ;


static int use_partial_volume_model = 1 ;
static int filter_type = FILTER_MEAN ;
static char *read_name = NULL ;
const char *Progname ;
static void usage_exit(int code) ;

#define T1 0
#define T2 1

#define GRAY_WHITE_BORDER    0
#define LAYER_IV_BORDER      1
#define PIAL_BORDER          2

static int label_only = 0 ;
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
  double       current_sigma, mean_dist ;
  int          msec, minutes, seconds, n_averages, i, j, piter, witer, l4iter, skip ;
  Timer start ;
  MRI          *mri, *mri_dist ;
  MRI_SURFACE  *mris ;
  TRANSFORM    *transform ;
  LTA          *lta ;
  VERTEX_PARMS *vp ;
  char         base_name[STRLEN], fname[STRLEN];
  const char *hemi ;

  witer = piter = l4iter = 0 ;
  dp.max_dist = MAX_DIST ;
  dp.max_ig_intensity_offset = 50 ;
  dp.search_multiple_angles = 0 ;

  dp.max_sg_intensity_offset = 50 ;
  dp.max_wm_intensity_offset = 30 ;
  dp.max_stria_intensity_offset = 20 ;
  dp.dark_csf = 0 ;
  dp.try_stria_model = 1 ;
  dp.outside_width = 0.0 ;
  dp.use_prior = 0 ;
  dp.inside_width = 1.5 ;
  dp.use_max_grad = 0 ;
  dp.error_type = ERROR_FUNC_L2 ;
  dp.which_border = GRAY_WHITE_BORDER;
  dp.wm_val = 110 ;
  dp.infra_granular_val = 160;
  dp.stria_val = 135 ;
  dp.stria_width = 0.3 ;
  dp.supra_granular_val = 190 ;  // for T2*-weighted ex vivo dat
  dp.contrast_type = T2 ;
  parms.sigma = dp.sigma = 3;
  dp.intensity_sigma = 15 ;
  parms.check_tol = 1 ;
  parms.tol = 1e-3 ;
  dp.max_sg_width = 3 ;
  dp.max_ig_width = 6 ;
  dp.min_sg_width = .5 ;
  dp.min_ig_width = .75 ;

  nargs = handleVersionOption(argc, argv, "mris_deform");
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
  parms.dt = 0.5 ;
  parms.l_intensity = 0.1 ;
  parms.l_intensity = 1e-5;
  parms.l_location = 10 ;

  Gx = Gy = Gz = -1 ;
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
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s", Progname, argv[1]) ;
  MRIScomputeMetricProperties(mris) ;
  MRISstoreMetricProperties(mris) ;


  if (MRISreadCanonicalCoordinates(mris, "sphere") != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "") ;

  hemi = mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh" ;
#if 0
  {
    MRI_SURFACE  *mris_sphere ;
    char         fname[STRLEN], path[STRLEN] ;
    FileNamePath(argv[1], path) ;
    sprintf(fname, "%s/%s.sphere", path, hemi) ;
    mris_sphere = MRISread(fname) ;
    mris_sphere->status = MRIS_SPHERE ;
    MRIScomputeMetricProperties(mris_sphere) ;

    mris->status = MRIS_SPHERE ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;

    DiagBreak() ;
  }
#endif
  vp = (VERTEX_PARMS *)calloc(mris->nvertices, sizeof(VERTEX_PARMS)) ;
  if (vp == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate VERTEX_PARMS", Progname) ;
  parms.user_parms = (void *)vp ;
  mris->user_parms = (void *)vp ;
  for (i = 0 ; i < mris->nvertices ; i++)
    mris->vertices[i].vp = (void *)(&(vp[i])) ;

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

  if (read_flag) // initialize  with previously computed surfaces
  {
    LABEL *v1, *v1_prior ;
    char  base_read_name[STRLEN] ;
    
    cp = strstr(read_name, "lh.") ;
    if (cp == NULL)
      cp = strstr(read_name, "rh.") ;
    if (cp == NULL)
      FileNameExtension(read_name,base_read_name) ;  // remove hemi (e.g. lh.)
    else
      strcpy(base_read_name, cp+3) ;
    printf("using %s as base read name\n", base_read_name) ;
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

    sprintf(fname, "%s.white", read_name) ;
    if (MRISreadWhiteCoordinates(mris, fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read white coords from %s",
                Progname, fname) ;
    vp_copy_from_surface(mris, WHITE_VERTICES, WHITE_VERTICES);

    sprintf(fname, "%s.layerIV", read_name) ;
    if (MRISreadWhiteCoordinates(mris, fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read layer IV coords from %s",
                Progname, fname) ;
    vp_copy_from_surface(mris, WHITE_VERTICES, LAYERIV_VERTICES);

    sprintf(fname, "%s.pial", read_name) ;
    if (MRISreadWhiteCoordinates(mris, fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read pial coords from %s",
                Progname, fname) ;
    vp_copy_from_surface(mris, WHITE_VERTICES, PIAL_VERTICES);
    vp_copy_to_surface(mris, PIAL_VERTICES, PIAL_VERTICES) ;
    vp_copy_to_surface(mris, WHITE_VERTICES, WHITE_VERTICES) ;
    compute_distances_from_positions(mris) ;

    if (label_only)
    {
      v1 = label_v1(mris, mri, &dp) ;
      sprintf(fname, "%s.%s.V1.label", hemi, parms.base_name);
      printf("writing v1 estimated position to %s\n", fname) ;
      if (v1)
      {
        LabelWrite(v1, fname) ;
        LabelFree(&v1) ;
      }
      vp_copy_to_surface_vals(mris, DEEP_RATIO, &dp) ;
      //      MRISsoapBubbleVals(mris, 100) ; 
      MRISaverageVals(mris, vavgs) ;
      sprintf(fname, "%s.%s.deep_ratio.mgz", hemi, parms.base_name);
      printf("writing ratio of deep IG to layer IV to %s\n", fname) ;
      MRISwriteValues(mris, fname) ;

      // compute midway surface and see how it does
      {
        int vno ;
        VERTEX_PARMS *vp ;
        VERTEX       *v ;
        for (vno = 0 ; vno < mris->nvertices ; vno++)
        {
          v = &mris->vertices[vno] ;
          vp = (VERTEX_PARMS *)(v->vp) ;
          if (v->ripflag)
            continue ;
          vp->l4_dist = vp->pial_dist/2 ;
        }
        label_v1(mris, mri, &dp) ;
        vp_copy_to_surface_vals(mris, DEEP_RATIO, &dp) ;
        MRISaverageVals(mris, vavgs) ;
        sprintf(fname, "%s.%s.mid.deep_ratio.mgz", hemi, parms.base_name);
        printf("writing ratio of deep middle surface to superficial to %s\n", fname) ;
        MRISwriteValues(mris, fname) ;
      }
      exit(0) ;
    }
    sprintf(fname, "%s.%s.marked", hemi, base_read_name);
    if (MRISreadMarked(mris, fname)!= NO_ERROR)
      ErrorPrintf(Gerror, "could not read marked from %s", fname) ;
    MRISaverageVals(mris, vavgs) ;
    vp_copy_found_from_marked(mris) ;

    if (dp.use_intensity) // use precomputed intensity info instead of optimizing it
    {
      sprintf(fname, "%s.%s.ig_intensity.mgz", hemi, base_read_name);
      if (MRISreadValues(mris, fname)!= NO_ERROR)
        ErrorExit(Gerror, "") ;
      MRISaverageVals(mris, vavgs) ;
      vp_copy_from_surface_vals(mris, IG_INTENSITY, &dp) ;

      sprintf(fname, "%s.%s.wm_intensity.mgz", hemi, base_read_name);
      if (MRISreadValues(mris, fname)!= NO_ERROR)
        ErrorExit(Gerror, "") ;
      MRISaverageVals(mris, vavgs) ;
      vp_copy_from_surface_vals(mris, WM_INTENSITY, &dp) ;

      sprintf(fname, "%s.%s.sg_intensity.mgz", hemi, base_read_name);
      if (MRISreadValues(mris, fname)!= NO_ERROR)
        ErrorExit(Gerror, "") ;
      MRISaverageVals(mris, vavgs) ;
      vp_copy_from_surface_vals(mris, SG_INTENSITY, &dp) ;
    }
  }
  else // initialize all the  vp fields to something reasonable 
  {
    INTEGRATION_PARMS eparms ;

    memset(&eparms, 0, sizeof(eparms)) ;
    eparms.l_spring = .1;
    eparms.l_location = 1 ;
    // eparms.l_curv = 1.0 ;
    eparms.n_averages = 16 ;
    eparms.min_averages = 4 ;
    eparms.l_surf_repulse = .0 ;
    eparms.dt = 0.25 ;
    eparms.min_averages = 8 ;

    MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
    vp_copy_from_surface(mris, CURRENT_VERTICES, WHITE_VERTICES);
    if (Gdiag & DIAG_WRITE)
      INTEGRATION_PARMS_setFp(&eparms, stdout) ;
    MRISexpandSurface(mris, MIN(1, 4*mri->xsize), &eparms, 0, 1) ;
    if (Gdiag & DIAG_WRITE)
    {
      char fname[STRLEN] ;
      sprintf(fname, "%s.%s.expanded", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", parms.base_name);
      printf("writing expanded 'pial' surface to %s\n", fname) ;
      MRISwrite(mris, fname) ;
    }
    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
    vp_copy_from_surface(mris, CURRENT_VERTICES, PIAL_VERTICES);
    vp_copy_from_surface(mris, CURRENT_VERTICES, LAYERIV_VERTICES);
    vp_copy_from_surface(mris, CURRENT_VERTICES, WHITE_TARGETS);
    vp_copy_from_surface(mris, CURRENT_VERTICES, PIAL_TARGETS);
    vp_copy_from_surface(mris, CURRENT_VERTICES, LAYERIV_TARGETS);
    vp_copy_from_surface(mris, WHITE_VERTICES, WHITE_VERTICES);
    vp_copy_from_surface(mris, PIAL_VERTICES, PIAL_VERTICES);
  }
  MRISresetNeighborhoodSize(mris, 3) ; // to allow calculation of nbhd stats

  TransformInvert(transform, NULL) ;
  if (transform->type == MNI_TRANSFORM_TYPE ||
      transform->type == TRANSFORM_ARRAY_TYPE ||
      transform->type  == REGISTER_DAT) 
  {
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
    label = LabelInFOV(mris, mri, pad_voxels*mri->xsize) ;
  LabelRipRestOfSurface(label, mris) ;
  if (threshold > 0)
    rip_dark_vertices(mris, mri, threshold) ;
  rip_vertices_with_no_faces(mris) ;

  max_averages = parms.n_averages ;
  mri_dist = MRIcloneDifferentType(mri, MRI_FLOAT) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    char fname[STRLEN] ;
    MRISstoreRipFlags(mris) ; MRISunrip(mris) ;
    sprintf(fname, "%s.%s.wnormals.mgz", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", parms.base_name);
    printf("writing surface normals to %s\n", fname) ;
    MRISwriteNormals(mris, fname) ;
    MRISrestoreRipFlags(mris) ;
  }
  for (skip = 4, i = 0, current_sigma = dp.sigma, n_averages = max_averages ; 
       n_averages >= min_averages ; 
       n_averages /= 2, current_sigma /= 2, i++, skip--)
  {
    char fname[STRLEN];

    if (skip < 1)
      skip = 1 ;
    dp.sigma = parms.sigma = current_sigma ;
    dp.max_dist = MIN(3*dp.sigma, dp.max_dist) ;  // no point going out that far
    parms.n_averages = n_averages ;
    printf("----------- deforming surfaces with %d smoothing iterations, sigma=%2.3f, skip=%d, max_dist=%2.1f -------\n", 
           n_averages, current_sigma, skip, dp.max_dist) ;
#define MAX_ITERATIONS_PER_SCALE    3
    for (j = 0 ; j < MAX_ITERATIONS_PER_SCALE ; j++)
    {
      dp.max_dist = MIN(3*dp.sigma, dp.max_dist) ;  // no point going out that far
      if (dp.use_intensity != 0)
      {
        vp_copy_to_surface(mris, WHITE_VERTICES, CURRENT_VERTICES) ;
        vp_set(mris, WM_INTENSITY_OFFSET, 0) ;
        if (0)
        {
          printf("computing distance transform\n") ;
          MRIScomputeDistanceToSurface(mris, mri_dist, mri_dist->xsize) ;
          MRIscalarMul(mri_dist, mri_dist, -1) ; // make inside negative
          compute_intensity_offsets(mris, mri, mri_dist, &dp) ;
        }
        vp_copy_to_surface_vals(mris, WM_INTENSITY_OFFSET, &dp) ;
        MRISsoapBubbleVals(mris, 100) ; 
        filter_vals(mris, filter_type) ;
        vp_copy_from_surface_vals(mris, WM_INTENSITY_OFFSET, &dp) ;

        dp.fix_intensities = 1 ;
      }
      mean_dist = compute_targets(mris, mri, current_sigma, &dp, skip) ;
      vp_copy_to_surface_vals(mris, WM_INTENSITY_OFFSET, &dp) ;
      MRISsoapBubbleVals(mris, 100) ; 
      MRISaverageVals(mris, vavgs) ;
      filter_vals(mris, filter_type) ;
      vp_copy_from_surface_vals(mris, WM_INTENSITY_OFFSET, &dp) ;
      //      dp.fix_intensities = 0 ;

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
            sprintf(fname, "%s.%s.V1.%d.label", hemi, base_name,ino);
            printf("writing v1 estimated position to %s\n", fname) ;
            LabelWrite(v1, fname) ;
            LabelFree(&v1) ;
          }
        }

        vp_copy_dist_to_surface_vals(mris, LAYERIV_VERTICES) ;
        sprintf(fname, "%s.%s.layerIV.height.%d.mgz", hemi, base_name,ino);
        MRISwriteValues(mris, fname) ;

        vp_copy_to_surface_vals(mris, WM_INTENSITY, &dp) ;
        sprintf(fname, "%s.%s.wtargets.%d.mgz", hemi, base_name,ino);
        MRISwriteValues(mris, fname) ;

        vp_copy_to_surface_vals(mris, RMS, &dp) ;
        sprintf(fname, "%s.%s.rms.%d.mgz", hemi, base_name,ino);
        MRISwriteValues(mris, fname) ;

        vp_copy_to_surface_vals(mris, SG_INTENSITY, &dp) ;
        sprintf(fname, "%s.%s.ptargets.%d.mgz", hemi, base_name,ino);
        MRISwriteValues(mris, fname) ;
        vp_copy_to_surface_vals(mris, IG_INTENSITY, &dp) ;
        sprintf(fname, "%s.%s.l4targets.%d.mgz", hemi, base_name,ino);
        MRISwriteValues(mris, fname) ;

        MRISstoreRipFlags(mris) ; MRISunrip(mris) ;

        vp_copy_to_surface(mris, WHITE_TARGETS, CURRENT_VERTICES) ;
        sprintf(fname, "%s.%s.wtarget.%d", hemi, base_name, ino);
        printf("writing white target locations to %s\n", fname) ;
        MRISwrite(mris, fname) ;

        vp_copy_to_surface(mris, PIAL_TARGETS, CURRENT_VERTICES) ;
        sprintf(fname, "%s.%s.ptarget.%d", hemi, base_name, ino);
        printf("writing pial target locations to %s\n", fname) ;
        MRISwrite(mris, fname) ;

        vp_copy_to_surface(mris, LAYERIV_TARGETS, CURRENT_VERTICES) ;
        sprintf(fname, "%s.%s.l4target.%d", hemi, base_name, ino);
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
      vp_copy_to_surface_vals(mris, WM_INTENSITY, &dp);
      MRISsoapBubbleVals(mris, 100) ;
      MRISaverageVals(mris, vavgs) ;
      parms.start_t = witer ;
      printf("positioning white matter surface\n") ;
      MRISpositionSurface(mris, mri, mri, &parms) ;  // move white matter surface
      witer = parms.start_t ; 
      vp_copy_from_surface(mris, CURRENT_VERTICES, WHITE_VERTICES) ;
#if 0
      if (i == 0 && j == 0 && read_flag == 0) // first time only
      {
        vp_copy_from_surface(mris, CURRENT_VERTICES, PIAL_VERTICES);
        vp_copy_from_surface(mris, CURRENT_VERTICES, LAYERIV_VERTICES);
        continue ; // don't trust surface normal in first pass
      }
#endif

      MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ; /* save white-matter */

      if (white_only == 0)
      {
        // do pial surface
        sprintf(parms.base_name, "%s.pial", base_name) ;
        sprintf(dp.base_name, "%s.pial", base_name) ;
        vp_copy_to_surface(mris, PIAL_VERTICES, CURRENT_VERTICES) ;
        vp_copy_to_surface(mris, PIAL_TARGETS, TARGET_VERTICES) ;
        vp_copy_to_surface(mris, LAYERIV_VERTICES, ORIGINAL_VERTICES) ; // for surf repulse
        parms.l_surf_repulse = 5 ;  // repel from layer IV surface
        vp_copy_to_surface_vals(mris, SG_INTENSITY, &dp);
        MRISsoapBubbleVals(mris, 100) ;
        MRISaverageVals(mris, vavgs) ;
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
        vp_copy_to_surface_vals(mris, IG_INTENSITY, &dp);
        MRISsoapBubbleVals(mris, 100) ;
        MRISaverageVals(mris, vavgs) ;
        parms.start_t = l4iter ;
        printf("positioning layer IV surface\n") ;
        MRISpositionSurface(mris, mri, mri, &parms) ;  // move layer IV surface
        parms.l_osurf_repulse = 0 ;  // turn off repulsion from pial surface 
        l4iter = parms.start_t ;
        vp_copy_from_surface(mris, CURRENT_VERTICES, LAYERIV_VERTICES) ;
        strcpy(parms.base_name, base_name) ; // restore default
      }
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


  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  MRISunrip(mris) ;
  sprintf(fname, "%s.white", argv[4]) ;
  vp_copy_to_surface(mris, WHITE_VERTICES, CURRENT_VERTICES) ;
  printf("writing white final surface position to %s\n", fname) ;
  MRISwrite(mris, fname) ;

  vp_copy_dist_to_surface_vals(mris, LAYERIV_VERTICES) ;
  //  MRISsoapBubbleVals(mris, 100) ; 
  sprintf(fname, "%s.%s.layerIV.height.mgz", hemi, base_name);
  printf("writing layer IV height to %s\n", fname) ;
  MRISwriteValues(mris, fname) ;

  vp_copy_to_surface_vals(mris, SG_INTENSITY, &dp) ;
  //  MRISsoapBubbleVals(mris, 100) ; 
  sprintf(fname, "%s.%s.sg_intensity.mgz", hemi, base_name);
  printf("writing SG intensity to %s\n", fname) ;
  MRISwriteValues(mris, fname) ;

  vp_copy_to_surface_vals(mris, WM_INTENSITY, &dp) ;
  //  MRISsoapBubbleVals(mris, 100) ; 
  sprintf(fname, "%s.%s.wm_intensity.mgz", hemi, base_name);
  printf("writing WM intensity to %s\n", fname) ;
  MRISwriteValues(mris, fname) ;

  vp_copy_to_surface_vals(mris, IG_INTENSITY, &dp) ;
  //  MRISsoapBubbleVals(mris, 100) ; 
  sprintf(fname, "%s.%s.ig_intensity.mgz", hemi, base_name);
  printf("writing infragranular intensity to %s\n", fname) ;
  MRISwriteValues(mris, fname) ;

  vp_copy_to_surface_vals(mris, IG_WM_RATIO, &dp) ;
  //  MRISsoapBubbleVals(mris, 100) ; 
  sprintf(fname, "%s.%s.ig_wm_ratio.mgz", hemi, base_name);
  printf("writing infragranular/wm ratio to %s\n", fname) ;
  MRISwriteValues(mris, fname) ;

  vp_copy_to_surface_vals(mris, SG_WM_RATIO, &dp) ;
  //  MRISsoapBubbleVals(mris, 100) ; 
  sprintf(fname, "%s.%s.sg_wm_ratio.mgz", hemi, base_name);
  printf("writing supragranular/wm intensity ratio to %s\n", fname) ;
  MRISwriteValues(mris, fname) ;

  vp_copy_to_surface_vals(mris, SG_IG_RATIO, &dp) ;
  //  MRISsoapBubbleVals(mris, 100) ; 
  sprintf(fname, "%s.%s.sg_ig_ratio.mgz", hemi, base_name);
  printf("writing supragranular/infrgranular intensity ratio to %s\n", fname) ;
  MRISwriteValues(mris, fname) ;

  sprintf(fname, "%s.pial", argv[4]) ;
  vp_copy_to_surface(mris, PIAL_VERTICES, CURRENT_VERTICES) ;
  printf("writing pial final surface position to %s\n", fname) ;
  MRISwrite(mris, fname) ;

  sprintf(fname, "%s.layerIV", argv[4]) ;
  vp_copy_to_surface(mris, LAYERIV_VERTICES, CURRENT_VERTICES) ;
  printf("writing pial final surface position to %s\n", fname) ;
  MRISwrite(mris, fname) ;

  vp_copy_to_surface_vals(mris, DEEP_RATIO, &dp) ;
  //      MRISsoapBubbleVals(mris, 100) ; 
  MRISaverageVals(mris, vavgs) ;
  sprintf(fname, "%s.%s.deep_ratio.mgz", hemi, parms.base_name);
  printf("writing ratio of deep IG to layer IV to %s\n", fname) ;
  MRISwriteValues(mris, fname) ;

  sprintf(fname, "%s.marked", argv[4]) ;
  printf("writing vertex marks to %s\n", fname) ;
  MRISwriteMarked(mris, fname) ;

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
  } else if (!stricmp(option, "white_only")) {
    white_only = 1 ;
    fprintf(stderr, "only deforming white matter surface locations\n")  ;
  } else if (!stricmp(option, "ico")) {
    dp.search_multiple_angles = 1 ;
    fprintf(stderr, "searching multiple surface normal angles in optimization\n") ;
  } else if (!stricmp(option, "max_offset") || !stricmp(option, "max_wm_offset")) {
    dp.max_wm_intensity_offset = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using max intensity offset %2.0f\n", dp.max_wm_intensity_offset) ;
  } else if (!stricmp(option, "max_ig_offset")) {
    dp.max_ig_intensity_offset = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using max IG intensity offset %2.0f\n", dp.max_ig_intensity_offset) ;
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
  } else if (!stricmp(option, "bright_csf")) {
    dp.dark_csf = 0 ;
    dp.outside_val = 300 ;
    fprintf(stderr, "assuming liquid outside brain is bright\n") ;
  } else if (!stricmp(option, "thresh")) {
    threshold = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "ripping vertices that initially have intensities < %f\n", threshold) ;
  } else if (!stricmp(option, "pad")) {
    pad_voxels = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "ripping %d voxel region at orders\n", pad_voxels) ;
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
    if (dp.dark_csf)
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
    printf("using %2.1f as mean stria value\n", dp.stria_val) ;
    nargs = 1 ;
  } else if (!stricmp(option, "mean")) {
    filter_type = FILTER_MEAN ;
    printf("using mean filtering instead of median\n") ;
  } else if (!stricmp(option, "intensity")) {
    dp.use_intensity = 1 ;
    dp.fix_intensities = 1 ;
    printf("using intensity histograms to nudge white matter surface\n") ;
#if 0
    nargs = 1 ;
    printf("using intensity morphing with lambda = %2.3f\n", parms.l_intensity) ;
    parms.l_intensity = atof(argv[2]) ;
    parms.l_location = 0 ;
#endif
  } else if (!stricmp(option, "gaussian")) {
    filter_type = FILTER_GAUSSIAN ;
    printf("using gaussian filtering instead of median\n") ;
  } else if (!stricmp(option, "nofilter")) {
    filter_type = FILTER_NONE ;
    printf("disabling filtering\n") ;
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
    if (dp.dark_csf)
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
  } else if (!stricmp(option, "label_only")) {
    label_only = 1 ;
    printf("labeling V1 and exiting\n") ;
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
    nargs = 1 ;
    read_name = argv[2] ;
    printf("reading surfaces to initialize with names prefixed by %s\n", read_name) ;
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
static int
compute_intensity_offsets(MRI_SURFACE *mris, MRI *mri, MRI *mri_dist, DP *dp)
{
  VERTEX  *v ;
  VERTEX_PARMS *vp = NULL ;
  int       vno, nfound, nmissed, peak, wsize ;
  HISTOGRAM *h, *hs ;
  double    xv, yv, zv, mean_white_border, inward_dist, outward_dist ;

#define WSIZE_MM  5 // diameter of region to histogram within
  wsize = WSIZE_MM / mri->xsize ;  
  inward_dist = outward_dist = dp->max_dist ;
  outward_dist = MIN(4*dp->sigma, outward_dist) ;
  MRIScomputeSurfaceNormals(mris, CURRENT_VERTICES, 10) ;
  MRIScomputeSurfaceNormals(mris, WHITE_VERTICES, 10) ;
  MRIScomputeSurfaceNormals(mris, PIAL_VERTICES, 10) ;
  MRISclearMarks(mris) ;
  mean_white_border = 0 ;
  for (nfound = nmissed = vno = 0 ; vno < mris->nvertices ; vno++)
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
      mean_white_border += vp->wm_val ;
      vp->found = 1;
      nfound++ ;
    }
    else
    {
      vp->found = 0 ;
      nmissed++ ;
    }
    HISTOfree(&h) ; HISTOfree(&hs) ;
  }
  
  if (nfound > 0)
    mean_white_border /= (double)nfound ;

  if (nmissed > 10000)
    DiagBreak() ;
  printf("%d vertices found (%2.2f%%, %d missed), mean wm border %2.1f\n",
         nfound, 100*nfound / (double)(nfound+nmissed), nmissed, mean_white_border) ;
  if (nfound == 0)
    DiagBreak() ;
  return(NO_ERROR) ;
}



static double
compute_targets(MRI_SURFACE *mris, MRI *mri, double sigma, DP *dp, int skip)
{
  static int i= 0 ;
  VERTEX_PARMS *vp = NULL ;
  int     vno, nfound, nmissed ;
  MRI     *mri_white = MRIclone(mri, NULL) ;
  MRI     *mri_pial = MRIclone(mri, NULL) ;
  MRI     *mri_l4 = MRIclone(mri, NULL) ;
  double  xv, yv, zv, d, target_val, xr, yr, zr, grad, max_grad, 
    target_dist,mean_white_border, mean_pial_border, mean_l4_border,
    mean_dist, inward_dist, outward_dist, mean_abs_dist ;

  max_wm_val = dp->wm_val + dp->max_wm_intensity_offset ;
  min_wm_val = dp->wm_val - dp->max_wm_intensity_offset ;
  inward_dist = outward_dist = dp->max_dist ;
  MRIScomputeSurfaceNormals(mris, CURRENT_VERTICES, 10) ;
  MRIScomputeSurfaceNormals(mris, WHITE_VERTICES, 10) ;
  MRIScomputeSurfaceNormals(mris, PIAL_VERTICES, 10) ;
  MRISclearMarks(mris) ;
  mean_abs_dist = 0.0 ;
  mean_white_border = mean_l4_border = mean_pial_border = 0 ;
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    MRISstoreRipFlags(mris) ; MRISunrip(mris) ;
    sprintf(fname, "%s.%s.normals.%d.init.mgz", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", parms.base_name, i);
    printf("writing surface normals to %s\n", fname) ;
    MRIScomputeNormals(mris) ;
    MRISwriteNormals(mris, fname) ;
    MRISrestoreRipFlags(mris) ;
  }
  compute_normals(mris, (VERTEX_PARMS *)mris->user_parms) ;
  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    MRISstoreRipFlags(mris) ; MRISunrip(mris) ;
    sprintf(fname, "%s.%s.normals.%d.final.mgz", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", parms.base_name,i);
    printf("writing surface normals to %s\n", fname) ;
    vp_copy_to_surface(mris, SURFACE_NORMALS, SURFACE_NORMALS) ;
    MRISwriteNormals(mris, fname) ;
    MRISrestoreRipFlags(mris) ;
    i++ ;
  }
  else if (getenv("READ_NORMALS"))
  {
    static int i = 0 ;
    char fname[STRLEN] ;
    if (i == 0)
    {
      sprintf(fname, "%s.%s.normals.mgz", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", parms.base_name);
      printf("reading surface normals from %s\n", fname) ;
      MRISreadNormals(mris, fname) ;
      vp_copy_from_surface(mris, SURFACE_NORMALS, SURFACE_NORMALS) ;
      i++ ;
    }
  }
  for (mean_dist=0.0, nfound = nmissed = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if ((vno % 5000) == 0)
      printf("%d of %d vertices processed (%2.1f%%)\n",
             vno, mris->nvertices, 100.0*vno/mris->nvertices) ;
    VERTEX * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    MRISvertexToVoxel(mris, v, mri, &xv, &yv, &zv) ;
    if (MRIindexNotInVolume(mri, xv, yv, zv))
    {
      v->ripflag = 1 ;
      continue ;
    }
    vp = (VERTEX_PARMS *)(v->vp) ;
    if (dp->use_max_grad)
    {
      max_grad = 0 ; target_val = -1 ; target_dist = -1000 ;
      for (d = -max_dist ; d <= max_dist ; d += dp->step)
      {
        xr = vp->wx + d*vp->nx ;
        yr = vp->wy + d*vp->ny ;
        zr = vp->wz + d*vp->nz ;
        MRISsurfaceRASToVoxelCached(mris, mri, xr, yr, zr, &xv, &yv, &zv);
        grad = MRISnormalIntensityGradient(mris, mri,
                                           xr, yr, zr, vp->nx, vp->ny, vp->nz, sigma) ;
        if (dp->contrast_type == T2 && grad < 0)
          grad = 0 ;
        else if (dp->contrast_type == T1 && grad > 0)
          grad = 0 ;
        if (grad > max_grad)
        {
          double val_inside, val_outside, xs, ys, zs ;
          xs = xr-0.5*vp->nx ; ys = yr-0.5*vp->ny ; zs = zr-0.5*vp->nz ;
          MRISsurfaceRASToVoxelCached(mris, mri, xs, ys, zs, &xv, &yv, &zv);
          MRIsampleVolume(mri, xv, yv, zv, &val_inside) ;
          xs = xr+0.5*vp->nx ; ys = yr+0.5*vp->ny ; zs = zr+0.5*vp->nz ;
          MRISsurfaceRASToVoxelCached(mris, mri, xs, ys, zs, &xv, &yv, &zv);
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
      if (target_dist > -100) // found wm boundary
      {
        if (vno == Gdiag_no)
          DiagBreak() ;
        vp->white_dist = target_dist ;
        vp->wtx = vp->wx + target_dist * vp->nx ;
        vp->wty = vp->wy + target_dist * vp->ny ;
        vp->wtz = vp->wz + target_dist * vp->nz ;
        vp->l4tx = vp->l4x ; vp->l4ty = vp->l4y ; vp->l4tz = vp->l4z ;
        vp->ptx = vp->px ; vp->pty = vp->py ; vp->ptz = vp->pz ;
        vp->l4_dist = sqrt(SQR(vp->l4x-vp->wx)+SQR(vp->l4y-vp->wy)+SQR(vp->l4z-vp->wz)) ;
        vp->pial_dist = sqrt(SQR(vp->px-vp->wx)+SQR(vp->py-vp->wy)+SQR(vp->pz-vp->wz)) ;
      }
      else
      {
        vp->found = 0 ;
        continue ;
      }
    }
#if 0
    else  // compute correlation with predicted profile
#endif
    {
      if (vno == Gdiag_no)  
        DiagBreak() ;          // break here XXX
      target_dist = 
        find_optimal_locations(mris, mri, vno, inward_dist, outward_dist, dp, vp, skip) ;
      if (vno == Gdiag_no)
        DiagBreak() ;  // zzz

      if (dp->use_intensity && target_dist < -100)
      {
        max_grad = 0 ; target_val = -1 ; target_dist = -1000 ;
        for (d = -max_dist ; d <= max_dist ; d += .1)
        {
          xr = v->x + d*vp->nx ;
          yr = v->y + d*vp->ny ;
          zr = v->z + d*vp->nz ;
          MRISsurfaceRASToVoxelCached(mris, mri, xr, yr, zr, &xv, &yv, &zv);
          grad = MRISnormalIntensityGradient(mris, mri,
                                             xr, yr, zr, vp->nx, vp->ny, vp->nz, sigma) ;
          if (dp->contrast_type == T2 && grad < 0)
            grad = 0 ;
          else if (dp->contrast_type == T1 && grad > 0)
            grad = 0 ;
          if (grad > max_grad)
          {
            double val_inside, val_outside, xs, ys, zs ;
            xs = xr-0.5*vp->nx ; ys = yr-0.5*vp->ny ; zs = zr-0.5*vp->nz ;
            MRISsurfaceRASToVoxelCached(mris, mri, xs, ys, zs, &xv, &yv, &zv);
            MRIsampleVolume(mri, xv, yv, zv, &val_inside) ;
            xs = xr+0.5*vp->nx ; ys = yr+0.5*vp->ny ; zs = zr+0.5*vp->nz ;
            MRISsurfaceRASToVoxelCached(mris, mri, xs, ys, zs, &xv, &yv, &zv);
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
        if (target_dist > -100) // found wm boundary
        {
          // leave l4 and pial in current location
          vp->l4tx = vp->l4x ; vp->l4ty = vp->l4y ; vp->l4tz = vp->l4z ;
          vp->ptx = vp->px ; vp->pty = vp->py ; vp->ptz = vp->pz ;
          vp->l4_dist = sqrt(SQR(vp->l4x-vp->wx)+SQR(vp->l4y-vp->wy)+SQR(vp->l4z-vp->wz)) ;
          vp->pial_dist = sqrt(SQR(vp->px-vp->wx)+SQR(vp->py-vp->wy)+SQR(vp->pz-vp->wz)) ;
        }
      }
      else if (dp->dark_csf == 0 && target_dist < -100)
      {
        double intensity_profile[MAX_PROFILE_LEN], step = dp->step ;
        int    i, wm_found, ig_found, sg_found ;
        extract_image_intensities(mris, mri, 0, max_dist, vp->wx, vp->wy, vp->wz,
                                  vp->nx,  vp->ny, vp->nz, dp->step, 1, 
                                  intensity_profile, NULL) ;
        wm_found = ig_found = sg_found = 0 ;
        for (i = 1 ; i < nint(max_dist / dp->step) ;  i++)
        {
          if (wm_found == 0 && intensity_profile[i] > dp->wm_val+dp->max_wm_intensity_offset)
          {
            wm_found = 1 ;
            vp->wtx = vp->wx + (i-1)*step*vp->nx ;
            vp->wty = vp->wy + (i-1)*step*vp->ny ;
            vp->wtz = vp->wz + (i-1)*step*vp->nz ;
            vp->white_dist = target_dist = (i-1)*step ;
          }

          if (wm_found == 1 && ig_found == 0 && 
              intensity_profile[i] > dp->supra_granular_val)
          {
            ig_found = 1 ;
            vp->l4tx = vp->wx + (i-1)*step*vp->nx ;
            vp->l4ty = vp->wy + (i-1)*step*vp->ny ;
            vp->l4tz = vp->wz + (i-1)*step*vp->nz ;
            vp->l4_dist = (i-1)*step ;
          }
          if (wm_found == 1 && ig_found == 1 && sg_found == 0 &&
              intensity_profile[i] > dp->outside_val)
          {
            sg_found = 1 ;
            vp->ptx = vp->wx + (i-1)*step*vp->nx ;
            vp->pty = vp->wy + (i-1)*step*vp->ny ;
            vp->ptz = vp->wz + (i-1)*step*vp->nz ;
            vp->pial_dist = (i-1)*step ;
          }
        }
        if(sg_found == 0 || wm_found == 0 || ig_found == 0)
          target_dist = -1000;
      }


      if (target_dist > -100) // found wm boundary
      {
        if (vno == Gdiag_no)
          DiagBreak() ;
        vp->white_dist = target_dist ;
        vp->wtx = vp->wx + target_dist * vp->nx ;
        vp->wty = vp->wy + target_dist * vp->ny ;
        vp->wtz = vp->wz + target_dist * vp->nz ;
      }
      if (target_dist < -100)
      {
        target_val = -1 ;
        target_dist = 0 ;
      }
      else
      {
        if (vno == Gdiag_no)
        {
          MRISsurfaceRASToVoxelCached(mris, mri, vp->l4tx, vp->l4ty, vp->l4tz, &xv, &yv, &zv);
          printf("v %d: layer IV target = %2.0f %2.0f %2.0f\n", vno, xv, yv, zv) ;
        }
        // make sure both targets are in volume
        MRISsurfaceRASToVoxelCached(mris, mri, vp->wtx, vp->wty, vp->wtz, &xv, &yv, &zv);
        if (vno == Gdiag_no)
          printf("v %d: wm target = %2.0f %2.0f %2.0f\n", vno, xv, yv, zv) ;
        if (MRIindexNotInVolume(mri, xv, yv, zv))
          target_val = -1 ;
        else
        {
          MRISsurfaceRASToVoxelCached(mris, mri, vp->ptx, vp->pty, vp->ptz, &xv, &yv, &zv);
          if (vno == Gdiag_no)
            printf("v %d: pial target = %2.0f %2.0f %2.0f\n", vno, xv, yv, zv) ;
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
      vp->ptx = vp->px ; vp->pty = vp->py ; vp->ptz = vp->pz ;
      vp->wtx = vp->wx ; vp->wty = vp->wy ; vp->wtz = vp->wz ;
      vp->l4tx = vp->l4x ; vp->l4ty = vp->l4y ; vp->l4tz = vp->l4z ;
    }
    else
    {
      mean_dist += target_dist ;
      mean_abs_dist += fabs(target_dist) ;
      nfound++ ;
      v->marked = 1 ; vp->found = 1 ;
      MRISsurfaceRASToVoxelCached(mris, mri, vp->wtx, vp->wty, vp->wtz, &xv, &yv, &zv);
      if (MRIindexNotInVolume(mri_white, xv, yv, zv) == 0)
      {
        MRIsampleVolume(mri, xv, yv, zv, &target_val) ;
        MRIsetVoxVal(mri_white, nint(xv), nint(yv), nint(zv), 0, target_val) ;
        vp->white_val = target_val ;
        mean_white_border += target_val ;
        if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
          DiagBreak() ;
      }
      
      MRISsurfaceRASToVoxelCached(mris, mri, vp->l4tx, vp->l4ty, vp->l4tz, &xv, &yv, &zv);
      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
        DiagBreak() ;
      if (MRIindexNotInVolume(mri_l4, xv, yv, zv) == 0)
      {
        MRIsampleVolume(mri, xv, yv, zv, &target_val) ;
        MRIsetVoxVal(mri_l4, nint(xv), nint(yv), nint(zv), 0, target_val) ;
        mean_l4_border += target_val ;
        vp->l4_val = target_val ;
      }

      MRISsurfaceRASToVoxelCached(mris, mri, vp->ptx, vp->pty, vp->ptz, &xv, &yv, &zv);
      if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz)
        DiagBreak() ;
      if (MRIindexNotInVolume(mri_l4, xv, yv, zv) == 0)
      {
        MRIsampleVolume(mri, xv, yv, zv, &target_val) ;
        mean_pial_border += target_val ;
        MRIsetVoxVal(mri_pial, nint(xv), nint(yv), nint(zv), 0, target_val) ;
        vp->pial_val = target_val ;
      }

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
#if 0
  if (dp->use_max_grad)
    MRISsoapBubbleVals(mris, 100) ;
  else
#endif
  {
    int    n, num, outliers ;
    double mn, std, mean_vdist, sigma_vdist ;

    // now do consistency check on white distances
      vp_copy_to_surface(mris, WHITE_TARGETS, TARGET_VERTICES) ;
      mean_vdist = MRIScomputeVertexSpacingStats(mris, &sigma_vdist, NULL,NULL,NULL,NULL,CURRENT_VERTICES) ;
    vp_copy_to_surface_dist(mris, WHITE_VERTICES);
    for (outliers = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->ripflag || v->marked == 0)
        continue ;
      for (mn = std = 0.0, num = n = 0 ; n < vt->vtotal ; n++)
      {
        VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
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
        for (n = 0 ; n < vt->vtotal ; n++)
        {
          VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
          if (vn->marked == 0 || vn->ripflag)
            continue ;
          fprintf(fp, "%d %f\n", vt->v[n], vn->d) ;
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
      else   // check target locations of this vertex and distance to neighbors
      {
        double        dx, dy, dz,d  ;
        VERTEX_PARMS  *vnp ;
        FILE          *fp = NULL;

        vp = (VERTEX_PARMS *)(v->vp) ;
        if (vno == Gdiag_no)
          fp = fopen("ldist.dat", "w") ;
        for (mn = std = 0.0, num = n = 0 ; n < vt->vnum ; n++)
        {
          VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
          if (vn->marked != 1 || vn->ripflag)
            continue ;
          vnp = (VERTEX_PARMS *)(vn->vp) ;
          dx = vp->wtx - vnp->wtx ;
          dy = vp->wty - vnp->wty ;
          dz = vp->wtz - vnp->wtz ;
          d = sqrt(dx*dx + dy*dy + dz*dz) ;
          num++ ;
          mn += d ;
          std += d*d ;
          if (vno == Gdiag_no)
            fprintf(fp, "%d %f\n", vt->v[n], d) ;
        }
        std = sqrt((std - mn*mn/num)/(num-1)) ;
        mn /= num ;
        d = fabs(mn - mean_vdist) ;
        if (d > 8*sigma_vdist)
          DiagBreak() ;
        if (vno == Gdiag_no)
          fclose(fp) ;
        if (d/sigma_vdist > 8 && d > 2.5)
        {
          if (vno == Gdiag_no)
            printf("discarding vertex %d white as an outlier (%2.1f - %2.1f) > 4*%2.1f\n",
                   vno, mn, mean_vdist, sigma_vdist) ;
          
          v->marked = 0 ; vp->found = 0 ;
          outliers++ ;
        }
      }
    }
    printf("%d white outliers replaced\n", outliers) ;

    vp_copy_to_surface_dist(mris, LAYERIV_VERTICES);
    vp_copy_to_surface(mris, LAYERIV_TARGETS, TARGET_VERTICES) ;
    mean_vdist = MRIScomputeVertexSpacingStats(mris, &sigma_vdist, NULL,NULL,NULL,NULL,CURRENT_VERTICES) ;

    // now do consistency check on pial distances
    for (outliers = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->ripflag || v->marked != 1)
        continue ;
      for (mn = std = 0.0, num = n = 0 ; n < vt->vtotal ; n++)
      {
        VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
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
        FILE *fp = NULL  ;
        fp = fopen("pout.dat", "w") ;
        for (n = 0 ; n < vt->vtotal ; n++)
        {
          VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
          if (vn->marked == 0 || vn->ripflag)
            continue ;
          fprintf(fp, "%d %f\n", vt->v[n], vn->d) ;
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
      else   // check target locations of this vertex and distance to neighbors
      {
        double        dx, dy, dz,d  ;
        VERTEX_PARMS  *vnp ;
        FILE          *fp = NULL ;

        vp = (VERTEX_PARMS *)(v->vp) ;
        if (vno == Gdiag_no)
          fp = fopen("ldist.dat", "w") ;
        for (mn = std = 0.0, num = n = 0 ; n < vt->vnum ; n++)
        {
          VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
          if (vn->marked != 1 || vn->ripflag)
            continue ;
          vnp = (VERTEX_PARMS *)(vn->vp) ;
          dx = vp->l4tx - vnp->l4tx ;
          dy = vp->l4ty - vnp->l4ty ;
          dz = vp->l4tz - vnp->l4tz ;
          d = sqrt(dx*dx + dy*dy + dz*dz) ;
          num++ ;
          mn += d ;
          std += d*d ;
          if (vno == Gdiag_no)
            fprintf(fp, "%d %f\n", vt->v[n], d) ;
        }
        std = sqrt((std - mn*mn/num)/(num-1)) ;
        mn /= num ;
        d = fabs(mn - mean_vdist) ;
        if (d > 8*sigma_vdist)
          DiagBreak() ;
        if (vno == Gdiag_no)
          fclose(fp) ;
        if (d/sigma_vdist > 8 && d > 2.5)
        {
          if (vno == Gdiag_no)
            printf("discarding vertex %d layer IV as an outlier (%2.1f - %2.1f) > 4*%2.1f\n",
                   vno, mn, mean_vdist, sigma_vdist) ;
          
          v->marked = 0 ; vp->found = 0 ;
          outliers++ ;
        }
      }
    }

    printf("%d layer IV outliers replaced\n", outliers) ;
    vp_copy_to_surface_dist(mris, LAYERIV_VERTICES);

    // now do consistency check on pial distances
    vp_copy_to_surface_dist(mris, PIAL_VERTICES);
    vp_copy_to_surface(mris, PIAL_TARGETS, TARGET_VERTICES) ;
    mean_vdist = MRIScomputeVertexSpacingStats(mris, &sigma_vdist, NULL,NULL,NULL,NULL,CURRENT_VERTICES) ;
    for (outliers = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->ripflag || v->marked != 1)
        continue ;
      for (mn = std = 0.0, num = n = 0 ; n < vt->vtotal ; n++)
      {
        VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
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
        for (n = 0 ; n < vt->vtotal ; n++)
        {
          VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
          if (vn->marked == 0 || vn->ripflag)
            continue ;
          fprintf(fp, "%d %f\n", vt->v[n], vn->d) ;
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
      else   // check target locations of this vertex and distance to neighbors
      {
        double        dx, dy, dz,d  ;
        VERTEX_PARMS  *vnp ;
        FILE          *fp = NULL ;

        vp = (VERTEX_PARMS *)(v->vp) ;
        if (vno == Gdiag_no)
          fp = fopen("pdist.dat", "w") ;
        for (mn = std = 0.0, num = n = 0 ; n < vt->vnum ; n++)
        {
          VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
          if (vn->marked != 1 || vn->ripflag)
            continue ;
          vnp = (VERTEX_PARMS *)(vn->vp) ;
          dx = vp->ptx - vnp->ptx ;
          dy = vp->pty - vnp->pty ;
          dz = vp->ptz - vnp->ptz ;
          d = sqrt(dx*dx + dy*dy + dz*dz) ;
          num++ ;
          mn += d ;
          std += d*d ;
          if (vno == Gdiag_no)
            fprintf(fp, "%d %f\n", vt->v[n], d) ;
        }
        std = sqrt((std - mn*mn/num)/(num-1)) ;
        mn /= num ;
        d = fabs(mn - mean_vdist) ;
        if (d > 8*sigma_vdist)
          DiagBreak() ;
        if (vno == Gdiag_no)
          fclose(fp) ;
        if (d/sigma_vdist > 8 && d > 2.5)
        {
          if (vno == Gdiag_no)
            printf("discarding vertex %d pial as an outlier (%2.1f - %2.1f) > 4*%2.1f\n",
                   vno, mn, mean_vdist, sigma_vdist) ;
          
          v->marked = 0 ; vp->found = 0 ;
          outliers++ ;
        }
      }
    }

    if (Gdiag_no >= 0)
    {
      VERTEX const * const v = &mris->vertices[Gdiag_no] ;
      vp = (VERTEX_PARMS *)(v->vp) ;
    }
    // now do soap bubble smoothing to fill in the missing values
    printf("%d pial outliers replaced\n", outliers) ;
    vp_copy_to_surface_dist(mris, PIAL_VERTICES);
    MRISsoapBubbleD(mris, 100) ;
    filter_distances(mris, filter_type) ;
    vp_copy_from_surface_dist(mris, PIAL_VERTICES) ;

    vp_copy_to_surface_dist(mris, WHITE_VERTICES);
    MRISsoapBubbleD(mris, 100) ;
    filter_distances(mris, filter_type) ;
    vp_copy_from_surface_dist(mris, WHITE_VERTICES) ;

    vp_copy_to_surface_dist(mris, LAYERIV_VERTICES);
    MRISsoapBubbleD(mris, 100) ;
    filter_distances(mris, filter_type) ;
    vp_copy_from_surface_dist(mris, LAYERIV_VERTICES) ;

    MRIclear(mri_white) ; MRIclear(mri_pial) ; MRIclear(mri_l4) ;
    recompute_target_locations(mris, mri_white, mri_l4, mri_pial, dp) ;
  }
  if (Gdiag & DIAG_WRITE)
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
                          double *intensity_profile, const char *fname)
{
  double   avg_val, val, x, y, z, xv, yv, zv, dx, dy, dz, dist, e1x, e1y, e1z, e2x, e2y, e2z, xt, yt,zt,
    norm, dxn, dyn, dzn, vstep ;
  int      i, nsamples, u, v, retval = NO_ERROR ;
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
  dxn = dx/norm ; dyn = dy/norm ; dzn = dz/norm ; 

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
        if (MRIindexNotInVolume(mri, xt, yt, zt))
          retval = ERROR_BADPARM ;
      }
    }
    xv += dx ; yv += dy ; zv += dz ;
    intensity_profile[i] = avg_val/((2*tsteps-1)*(2*tsteps-1));
    if (fp)
      fprintf(fp, "%f %f %f %f %f\n", dist, intensity_profile[i], xv, yv, zv) ;
   }

  if (fp)
    fclose(fp) ;
  return(retval) ;
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
find_next_peak(double *intensity_profile, int index, double *pmax_val, int len, int whalf) 
{
  int i, found, j ;

  for (i = index ; i < len ; i++)
  {
    for (found = 1, j = MAX(0, i-whalf) ;  j <= MIN(len-1,i+whalf) ; j++)
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
find_previous_peak(double *intensity_profile, int index, double *pmax_val, int len, int whalf) 
{
  int i, found, j ;

  for (i = index ; i >= 0 ; i--)
  {
    for (found = 1, j = MAX(0, i-whalf) ;  j <= MIN(len-1,i+whalf) ; j++)
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
find_next_valley(double *intensity_profile, int index, double *pmin_val, int len, int whalf) 
{
  int i, found, j ;

  for (i = index ; i < len-whalf ; i++)
  {
    for (found = 1, j = MAX(0, i-whalf) ;  j <= MIN(len-1,i+whalf) ; j++)
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
smooth_profile(double *src_intensity_profile, double *dst_intensity_profile, int len, double sigma) 
{
  int i ;
  HISTOGRAM *h, *hs ;

  h = HISTOalloc(len) ;
  
  for (i = 0 ; i < len ; i++)
  {
    h->bins[i] = i ;
    h->counts[i] = src_intensity_profile[i] ;
  }
  hs = HISTOsmooth(h, NULL, sigma) ;
  for (i = 0 ; i < len ; i++)
    dst_intensity_profile[i] = hs->counts[i] ;

  HISTOfree(&h) ; HISTOfree(&hs) ;
  return(NO_ERROR) ;
}

static int
mark_sg_index_okay_dark_csf(double *orig_intensity_profile, int *sg_index_okay, double step, int len, const char *fname, int current_index, double outside_val, int min_wm_index, int *sg_indices, int *pmax_len, int current_sg,
                   int max_out)
{
  int    index, i, j, min_index, max_index, pv_len ;
  double valley_val, peak_val, previous_peak_val, intensity_profile[MAX_PROFILE_LEN] ;
  int    valley_index, peak_index, previous_peak, whalf ;
  FILE   *fp = NULL ;

  pv_len = SMOOTH_LEN-1 + PV_LEN ;
  smooth_profile(orig_intensity_profile, intensity_profile, len, SMOOTH_LEN) ;
  
  whalf = nint(0.5/step) ;  // valley must be less than everything within 1/2 mm
  if (fname)
    fp = fopen(fname, "w") ;

  for (index = current_index ; index < len ; index++)
    if (orig_intensity_profile[index] <= outside_val)
      break ;
  max_index = index-1 ;
  for (index = current_index ; index >= 0 ; index--)
    if (orig_intensity_profile[index] <= outside_val)
      break ;
  min_index = index+1 ;
    
  for (index = min_wm_index ; index < len ; index++)
  {
    sg_index_okay[index] = 0 ;  // assume it's on an upwards slope unless violated below
    if (index < min_index || index > max_index)
      continue ;
    for (i = index ; i <= index+1 && i < len ; i++)
    {
      for (j = i+1 ; j <= index+SMOOTH_LEN && j < len ; j++)
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
      valley_index = find_next_valley(intensity_profile, index, &valley_val, len, whalf) ;
      if (valley_index < 0)
        continue ;
      valley_val = orig_intensity_profile[valley_index] ;

      if (valley_val > orig_intensity_profile[i]-.15*orig_intensity_profile[i]) // not a deep valley
      {
        previous_peak = find_previous_peak(intensity_profile,valley_index, &previous_peak_val,len,whalf);
        peak_index = find_next_peak(intensity_profile, valley_index, &peak_val, len, whalf) ;
        if (peak_index >= 0 && previous_peak >= 0)
        {
          peak_val = orig_intensity_profile[peak_index] ;
          previous_peak_val = orig_intensity_profile[previous_peak] ;

          if ((peak_val-previous_peak_val) > previous_peak_val*.1)
            sg_index_okay[index] = 0 ;  // if the next peak is a lot higher than this one
        }
      }
    }
    else  // on an upwards slope - mark everthing on the slope no good except peak
    {
      peak_index = find_next_peak(intensity_profile, index, &peak_val, len, whalf) ;
      for ( ; index < peak_index ; index++)
        sg_index_okay[index] = 0 ;
    }
    if (use_partial_volume_model && sg_index_okay[index])
    {
      int j ;
      /* smoothing and partial volume models mean that the edge location can move around
         some from the location of the peak, so test some additional locations
      */
      for (j = 1 ; index+j < len && j<=pv_len ; j++)
        sg_index_okay[index+j] = 1 ;
      index += pv_len; 
    }
  }
  {
    HISTO  *h ;
    int    peak ;
    double peak_val ;
    h = HISTOalloc(len) ;
    
    for (index = 0 ; index < len ; index++)
    {
      h->bins[index] = index ;
      h->counts[index] = orig_intensity_profile[index] ;
    }
    peak = HISTOfindHighestPeakInRegion(h, min_wm_index, max_index) ;
    if (peak >= 0)
    {
      peak_val = h->counts[peak] ;
      for (index = min_index ; index < max_index ; index++)
      {
        if ((orig_intensity_profile[index]*1.1 < peak_val) &&
            (orig_intensity_profile[index] < 190))
          sg_index_okay[index] = 0 ;
        else if (use_partial_volume_model)
          index += pv_len ; // skip partial volume locations
      }
    }
    HISTOfree(&h) ;
  }
  if (fp)
  {
    for (index = 0 ; index < len ; index++)
      fprintf(fp, "%d %d %f\n", index, sg_index_okay[index], intensity_profile[index]) ;
    fclose(fp) ;
  }

  // only allow indices that are within +- max_out of current_sg index
  for (index = 0 ; index < current_sg-max_out ; index++)
    sg_index_okay[index] = 0 ;
  for (index = current_sg+max_out+1 ; index < len ; index++)
    sg_index_okay[index] = 0 ;
  if (sg_indices)
  {
    int nsg ;
    for (index = nsg = 0 ; index < len ; index++)
      if (sg_index_okay[index])
      {
        sg_indices[nsg] = index ;
        nsg++ ;
      }
    *pmax_len  = nsg ;
  }
  return(NO_ERROR);
} 
static int
mark_sg_index_okay_bright_csf(double *orig_intensity_profile, int *sg_index_okay, 
                              double step, int len, const char *fname, int current_index, 
                              double outside_val, int min_wm_index, int *sg_indices, 
                              int *pmax_len, int current_sg,
                              int max_out)
{
  int    index, i, j, min_index, max_index, pv_len, whalf ;
  //  double valley_val, peak_val, previous_peak_val ;
  double intensity_profile[MAX_PROFILE_LEN] ;
  //  int    valley_index, peak_index, previous_peak ;
  FILE   *fp = NULL ;

  pv_len = SMOOTH_LEN-1 + PV_LEN ;
  smooth_profile(orig_intensity_profile, intensity_profile, len, SMOOTH_LEN) ;
  
  whalf = nint(0.5/step) ;  // valley must be less than everything within 1/2 mm
  if (fname)
    fp = fopen(fname, "w") ;

  // find something that is bright enough to be outside
  for (index = current_index ; index < len ; index++)
    if (orig_intensity_profile[index] >= outside_val)
      break ;
  max_index = index-1 ;
  for (index = current_index ; index >= 0 ; index--)
    if (orig_intensity_profile[index] > outside_val)
      break ;
  min_index = index+1 ;
    
  for (index = min_wm_index ; index < len ; index++)
  {
    sg_index_okay[index] = 0 ;  // assume it's on a downards slope unless violated below
    if (index < min_index || index > max_index)
      continue ;
    for (i = index ; i <= index+1 && i < len ; i++)
    {
      for (j = i+1 ; j <= index+SMOOTH_LEN && j < len ; j++)
        if (intensity_profile[j] > intensity_profile[i])
        {
          sg_index_okay[index] = 1 ;
          break ;
        }
      if (sg_index_okay[index])
        break ;  // already found a violation of monotonicity
    }

#if 0
    if (sg_index_okay[index])  
    {
      valley_index = find_next_valley(intensity_profile, index, &valley_val, len, whalf) ;
      if (valley_index < 0)
        continue ;
      valley_val = orig_intensity_profile[valley_index] ;

      if (valley_val > orig_intensity_profile[i]-.15*orig_intensity_profile[i]) // not a deep valley
      {
        previous_peak = find_previous_peak(intensity_profile,valley_index, &previous_peak_val,len,whalf);
        peak_index = find_next_peak(intensity_profile, valley_index, &peak_val, len, whalf) ;
        if (peak_index >= 0 && previous_peak >= 0)
        {
          peak_val = orig_intensity_profile[peak_index] ;
          previous_peak_val = orig_intensity_profile[previous_peak] ;

          if ((peak_val-previous_peak_val) > previous_peak_val*.1)
            sg_index_okay[index] = 0 ;  // if the next peak is a lot higher than this one
        }
      }
    }
    else  // on an upwards slope - mark everthing on the slope no good except peak
    {
      peak_index = find_next_peak(intensity_profile, index, &peak_val, len, whalf) ;
      for ( ; index < peak_index ; index++)
        sg_index_okay[index] = 0 ;
    }
    if (use_partial_volume_model && sg_index_okay[index])
    {
      int j ;
      /* smoothing and partial volume models mean that the edge location can move around
         some from the location of the peak, so test some additional locations
      */
      for (j = 1 ; index+j < len && j<=pv_len ; j++)
        sg_index_okay[index+j] = 1 ;
      index += pv_len; 
    }
#endif
  }
#if 0
  {
    HISTO  *h ;
    int    peak ;
    double peak_val ;
    h = HISTOalloc(len) ;
    
    for (index = 0 ; index < len ; index++)
    {
      h->bins[index] = index ;
      h->counts[index] = orig_intensity_profile[index] ;
    }
    peak = HISTOfindHighestPeakInRegion(h, min_wm_index, max_index) ;
    if (peak >= 0)
    {
      peak_val = h->counts[peak] ;
      for (index = min_index ; index < max_index ; index++)
      {
        if ((orig_intensity_profile[index]*1.1 < peak_val) &&
            (orig_intensity_profile[index] < 190))
          sg_index_okay[index] = 0 ;
        else if (use_partial_volume_model)
          index += pv_len ; // skip partial volume locations
      }
    }
    HISTOfree(&h) ;
  }
#endif
  if (fp)
  {
    for (index = 0 ; index < len ; index++)
      fprintf(fp, "%d %d %f\n", index, sg_index_okay[index], intensity_profile[index]) ;
    fclose(fp) ;
  }

  // only allow indices that are within +- max_out of current_sg index
  for (index = 0 ; index < current_sg-max_out ; index++)
    sg_index_okay[index] = 0 ;
  for (index = current_sg+max_out+1 ; index < len ; index++)
    sg_index_okay[index] = 0 ;
  if (sg_indices)
  {
    int nsg ;
    for (index = nsg = 0 ; index < len ; index++)
      if (sg_index_okay[index])
      {
        sg_indices[nsg] = index ;
        nsg++ ;
      }
    *pmax_len  = nsg ;
  }
  return(NO_ERROR);
} 

static int
mark_wm_index_okay(double *intensity_profile, int *wm_index_okay, 
                   int wm_len, double min_wm_intensity, double max_wm_intensity,
                   int current_index, int len, const char *fname, int max_in)
{
  int    index, inside, wm_samples, wm_start=-1, max_wm_len ;
  FILE   *fp = NULL ;

  if (fname)
    fp = fopen(fname, "w") ;

  memset(wm_index_okay, 0, len*sizeof(int)) ;
  max_wm_len = 0 ;
  for (inside = wm_samples = 0, index = current_index-(wm_len-1) ; index < len ; index++)
  {
    if ((intensity_profile[index] > max_wm_intensity+.75*(max_wm_intensity-min_wm_intensity)) &&
         (max_wm_len > wm_len/2))
      break ;  // don't go through very high intensity regions after in low ones
    if (intensity_profile[index] <= max_wm_intensity && 
        intensity_profile[index] >= min_wm_intensity)
    {
      if (++wm_samples >= wm_len)   // found interior
      {
        if (inside == 0)
          wm_start = index-wm_samples+1 ;
        inside = 1 ;
      }
      if (wm_samples > max_wm_len)
        max_wm_len = wm_samples ;
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

  max_wm_len = 0 ;
  for (inside = wm_samples = 0, index = current_index+wm_len-1 ; index >= 0 ; index--)
  {
    if ((intensity_profile[index] > max_wm_intensity+.75*(max_wm_intensity-min_wm_intensity)) &&
         (max_wm_len > wm_len/2))
      break ;  // don't go through very high intensity regions after in low ones
    if (intensity_profile[index] <= max_wm_intensity && intensity_profile[index] >= min_wm_intensity)
    {
      if (++wm_samples >= wm_len)   // found interior
      {
        if (inside == 0)
          wm_start = index+wm_samples-1 ;
        inside = 1 ;
      }
      if (wm_samples > max_wm_len)
        max_wm_len = wm_samples ;
    }
    else   // not in wm intensity range
    {
      if (inside) // everything else after here is invalid
        break ;
      inside = wm_samples = 0 ;
    }
  }

  // only allow indices that are within +- max_in of current index
  for (index = 0 ; index < current_index-max_in ; index++)
    wm_index_okay[index] = 0 ;
  for (index = current_index+max_in+1 ; index < len ; index++)
    wm_index_okay[index] = 0 ;
  if (inside)
  {
    for (index = wm_start - wm_samples ; index >= 0 ; index--)
      wm_index_okay[index] = 0 ;
    for (index = current_index ; index >= 0 && index > wm_start - wm_samples ; index--)
      wm_index_okay[index] = 1 ;
  }

  for (wm_start = 0 ; wm_start < len ; wm_start++)
    if (wm_index_okay[wm_start] > 0)
      break ;
  if (fp)
  {
    for (index = 0 ; index < len ; index++)
      fprintf(fp, "%d %d\n", index, wm_index_okay[index]) ;
    fclose(fp) ;
  }
  return(wm_start);
} 

#define DIST_NOT_FOUND  -10000
#define MAX_WM_LEN nint(5.0/step)
#define MIN_WM_LEN nint(1.0/step)
#define MAX_CSF_LEN nint(2.0/step)
#define MIN_CSF_LEN nint(.3/step)
/*
  find the optimal profile location and fill in the vp->w[xyz] and vp->p[xyz]
  for the predicted pial and white locations.
*/
static double
find_optimal_locations(MRI_SURFACE *mris, MRI *mri, int vno, 
                       double max_inwards, double max_outwards, DP *dp, VP *vp, int skip)

{
  int     start_index, current_index, min_gray_matter_len, max_wm_len, filter_type,
    max_wm_index, csf_len, max_csf_len, max_profile_len, profile_len, wm_len,
    best_wm_len, best_csf_len, best_start_index, sg_index_okay[MAX_PROFILE_LEN], index,
    wm_index_okay[MAX_PROFILE_LEN], n, sg_indices[MAX_PROFILE_LEN],nsg,wm_skip = skip ;
  double errors[MAX_PROFILE_LEN], best_wm_dist, prior,
    kernel[MAX_PROFILE_LEN], nx, ny, nz, base_intensity_profile[MAX_PROFILE_LEN], wm_dist ;
  double  wm_val, ig_val, sg_val, *intensity_profile, pdist, wm_profile_val ;
  double  step, rms, min_rms, offset ;
  int     ig_len, sg_len, best_ig_len=0, best_sg_len=0, min_ig_len, max_ig_len, ico_vno,
    min_sg_len, max_sg_len, best_profile = PROFILE_GENERIC, tsteps, min_wm_index,
    min_start_index ;
  double  best_wm_intensity_offset=0.0,min_wm_intensity_offset,max_wm_intensity_offset,
    wm_intensity_offset, min_wm_intensity, max_wm_intensity, best_nx, best_ny, best_nz,
    best_sg_intensity_offset, sg_intensity_offset, ig_intensity_offset=0,  best_ig_intensity_offset, 
    min_wm_len ;
  
  int error_type = dp->error_type ;
  double (*profile_error_func)(double *kernel, double *intensity_profile, int nsamples, const char *fname,
                               double step, double in_dist, double *errors);
  const char   *fname1, *fname2, *fname3, *fname4 ;
  double  stria_val, stria_intensity_offset, best_stria_intensity_offset = 0.0 ; 
  static MRI_SURFACE   *mris_ico = NULL;

  max_profile_len = current_index = min_start_index = min_wm_index = 0 ;
  pdist = 0.0 ;

  if (mris_ico == NULL)
  {
    mris_ico = ic642_make_surface(0, 0) ;
    MRISscaleVertexCoordinates(mris_ico, 1.0/100.0) ;
    MRIScomputeMetricProperties(mris_ico) ;
  }

  intensity_profile = base_intensity_profile ;
  if (vno == Gdiag_no || deform_debug)
  {
    fname1 = "profile.dat" ;
    fname2 = "intensity.dat" ;
    fname3 = "sg_okay.dat" ;
    fname4 = "wm_okay.dat" ;
  }
  else
    fname1 = fname2 = fname3 = fname4 = NULL ;

  vp->min_rms = min_rms = 1e10 ;
  dp->step = step = mri->xsize/2 ;
  min_wm_len = MIN_WM_LEN ;
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
  case ERROR_FUNC_L1:        filter_type = FILTER_MEDIAN ;      break ;
  case ERROR_FUNC_L2:        filter_type = FILTER_MEAN ;        break ;
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


  min_wm_intensity = dp->wm_val + min_wm_intensity_offset ;
  max_wm_intensity = dp->wm_val + max_wm_intensity_offset ;
  best_wm_dist = DIST_NOT_FOUND ;
  min_ig_len = nint(dp->min_ig_width/step) ; max_ig_len = nint(dp->max_ig_width/step) ;
  min_sg_len = nint(dp->min_sg_width/step) ; max_sg_len = nint(dp->max_sg_width/step) ;

  min_gray_matter_len = min_ig_len + min_sg_len ;

  best_wm_len = best_ig_len = best_sg_len = best_start_index = best_csf_len = 0 ;

  ico_vno = MRISfindClosestVertex(mris_ico, vp->nx, vp->ny, vp->nz, NULL, CURRENT_VERTICES);
  VERTEX_TOPOLOGY const * const vico = &mris_ico->vertices_topology[ico_vno] ;
  best_nx = vp->nx ; best_ny = vp->ny ; best_nz = vp->nz ;
  best_ig_intensity_offset = best_sg_intensity_offset = best_wm_intensity_offset = 0.0 ;
  for (n = -1 ; n < vico->vnum ; n++)
  {
    if (n >= 0 && dp->search_multiple_angles == 0)
      break ;
    if (n < 0)
    {
      nx = vp->nx ; ny = vp->ny ; nz = vp->nz ;
    }
    else
    {
      nx = mris_ico->vertices[vico->v[n]].nx ; 
      ny = mris_ico->vertices[vico->v[n]].ny ; 
      nz = mris_ico->vertices[vico->v[n]].nz ; 
    }

    // just to make things run a bit faster

    pdist = sqrt(SQR(vp->px-vp->wx) + SQR(vp->py-vp->wy) + SQR(vp->pz-vp->wz)) ;
    profile_len=max_profile_len = nint((max_inwards + max_outwards+pdist)/step) ;
    current_index = nint(max_inwards/step) ;
    if (dp->use_max_grad)  // white matter location is already found
    {
      int index, max_in_offset, max_out_offset, sg_index ;
      HISTOGRAM *h ;

      wm_skip = 1 ;
      current_index = nint((2*max_inwards)/step) ;
      max_profile_len = nint((2*max_inwards + 2*max_outwards+pdist)/step) ;
      extract_image_intensities(mris, mri, 2*max_inwards, 2*max_outwards+pdist, vp->wx, vp->wy, vp->wz, 
                                nx,ny,nz, step, tsteps, base_intensity_profile, fname2);

      h = HISTOalloc(max_profile_len) ;
      for (index = 0 ; index < max_profile_len ; index++)
      {
        h->bins[index] = index ;
        h->counts[index] = base_intensity_profile[index] ;
      }
#define WHALF  5
      min_wm_index = HISTOfindMaxDerivative(h, min_wm_intensity, max_wm_intensity, WHALF, 1) ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      min_wm_index = MAX(1, min_wm_index-(WHALF-1)) ;
      vp->white_dist = (min_wm_index - current_index) * step ;

      max_in_offset = nint(max_inwards / step) ;
      max_out_offset = nint(max_outwards / step) ;
      memset(wm_index_okay, 0, max_profile_len*sizeof(int)) ;
      memset(sg_index_okay, 0, max_profile_len*sizeof(int)) ;
      min_wm_len = 1 ; max_wm_len = min_wm_index ;
      wm_index_okay[min_wm_index] = 1 ;
      sg_index = current_index + nint(vp->pial_dist/step) ;
      for (nsg = 0, index = MAX(0,sg_index-max_in_offset) ; 
           index <= MIN(max_profile_len-1,sg_index+max_out_offset) ;
           index++)
      {
        sg_indices[nsg++] = index ;
        sg_index_okay[index] = 1 ;
      }
      min_start_index = 0 ;
      if (vno == Gdiag_no)
      {
        FILE *fp ;

        fp = fopen(fname4, "w") ;
        for (index = 0 ; index < max_profile_len ; index++)
          fprintf(fp, "%d %d %f\n", index, sg_index_okay[index], base_intensity_profile[index]);
        fclose(fp) ;

        fp = fopen(fname3, "w") ;
        for (index = 0 ; index < max_profile_len ; index++)
          fprintf(fp, "%d %d %f\n", index, sg_index_okay[index], base_intensity_profile[index]) ;
        fclose(fp) ;
      }
    }
#if 0
    else if (dp->fix_intensities)// mark a region around current locations
    {
      int index, max_in_offset, max_out_offset, sg_index ;

      max_in_offset = nint(max_inwards / step) ;
      max_out_offset = nint(max_outwards / step) ;
      extract_image_intensities(mris, mri, max_inwards, max_outwards+pdist, vp->wx, vp->wy, vp->wz, 
                                nx,ny,nz, step, tsteps, base_intensity_profile, fname2);
      memset(wm_index_okay, 0, max_profile_len*sizeof(int)) ;
      memset(sg_index_okay, 0, max_profile_len*sizeof(int)) ;
      for (min_wm_index = index = MAX(0,current_index-max_in_offset) ; 
           index <= MIN(max_profile_len-1,current_index+max_out_offset) ;
           index++)
        wm_index_okay[index] = 1 ;
      sg_index = current_index + nint(vp->pial_dist/step) ;
      for (nsg = 0, index = MAX(0,sg_index-max_in_offset) ; 
           index <= MIN(max_profile_len-1,sg_index+max_out_offset) ;
           index++)
      {
        sg_indices[nsg++] = index ;
        sg_index_okay[index] = 1 ;
      }
      min_start_index = 0 ;
      if (vno == Gdiag_no)
      {
        FILE *fp ;

        fp = fopen(fname4, "w") ;
        for (index = 0 ; index < max_profile_len ; index++)
          fprintf(fp, "%d %d %f\n", index, sg_index_okay[index], base_intensity_profile[index]);
        fclose(fp) ;

        fp = fopen(fname3, "w") ;
        for (index = 0 ; index < max_profile_len ; index++)
          fprintf(fp, "%d %d %f\n", index, sg_index_okay[index], base_intensity_profile[index]) ;
        fclose(fp) ;
      }
    }
#endif
    else
    {
      extract_image_intensities(mris, mri, max_inwards, max_outwards+pdist, vp->wx, vp->wy, vp->wz, 
                                nx,ny,nz, step, tsteps, base_intensity_profile, fname2);
      min_wm_index = mark_wm_index_okay(base_intensity_profile, wm_index_okay, nint(.4/step), 
                                        min_wm_intensity, max_wm_intensity, current_index, 
                                        max_profile_len, fname4, nint(max_inwards/step));
      if (dp->dark_csf)
        mark_sg_index_okay_dark_csf(base_intensity_profile, sg_index_okay, step, 
                                    max_profile_len, fname3, 
                                    current_index, dp->outside_val+10, min_wm_index,
                                    sg_indices, &nsg, current_index+nint(pdist/step), 
                                    nint(max_outwards/step)) ;
      else
        mark_sg_index_okay_bright_csf(base_intensity_profile, sg_index_okay, step, 
                                      max_profile_len, fname3, 
                                      current_index, dp->outside_val-10, min_wm_index,
                                      sg_indices, &nsg, current_index+nint(pdist/step), 
                                      nint(max_outwards/step)) ;

      min_start_index = min_wm_index ;
    }

    for (max_wm_index = min_wm_index+1 ; max_wm_index < max_profile_len-min_wm_len ; max_wm_index++)
      if (wm_index_okay[max_wm_index] == 0)
        break ;
    
    for (start_index = min_start_index ; start_index <= max_wm_index; start_index++)
    {
      offset = (current_index - start_index)*step ;
      intensity_profile = base_intensity_profile + start_index ;
      max_wm_len = max_profile_len-(min_gray_matter_len+1+start_index) ;
      max_wm_len = MIN(max_wm_len, MAX_WM_LEN) ;
      for (wm_len = min_wm_len ; wm_len <= max_wm_len; wm_len += wm_skip)
      {
        index = wm_len+start_index ;
        if (wm_index_okay[index] == 0)
          continue ;
        wm_dist = (start_index+wm_len-current_index)*step ;   
        if (dp->fix_intensities)
        {
          wm_val = dp->wm_val + vp->wm_intensity_offset ;
          wm_profile_val = compute_optimal_intensity(intensity_profile, 0, wm_len-2, filter_type) ;
        }
        else
        {
          wm_val = compute_optimal_intensity(intensity_profile, 0, wm_len-2, filter_type) ;
          if (wm_val < min_wm_intensity || wm_val > max_wm_intensity)
            continue ; // not in feasible region
          wm_profile_val = wm_val ;
        }
        wm_intensity_offset = wm_val - dp->wm_val ;
        for (ig_len = min_ig_len ; ig_len <= max_ig_len ; ig_len += skip)
        {
#if 0
          if (dp->fix_intensities)
            ig_val = dp->infra_granular_val + vp->ig_intensity_offset ;
          else
#endif
          {
            ig_val = compute_optimal_intensity(intensity_profile, wm_len+1, 
                                               wm_len+ig_len-2,filter_type);
            if (ig_val < dp->infra_granular_val-dp->max_ig_intensity_offset ||
                ig_val > dp->infra_granular_val+dp->max_ig_intensity_offset)
              continue ; // not in feasible region
          }
          ig_intensity_offset = ig_val - dp->infra_granular_val ;
          
          if (ig_val <= wm_profile_val+MIN_INTENISTY_PCT_OFFSET*wm_profile_val) // do some sanity checks
            continue ;

          for (index = 0 ; index < nsg ; index++)
          {
            sg_len = sg_indices[index]-(wm_len+ig_len+start_index) ;
            if (sg_len < min_sg_len || sg_len > max_sg_len)
              continue ;
#if 0
            if (dp->fix_intensities)
              sg_val = dp->supra_granular_val + vp->sg_intensity_offset ;
            else
#endif
            {
              sg_val = compute_optimal_intensity(intensity_profile, 
                                                 wm_len+ig_len+1,wm_len+ig_len+sg_len-2,
                                                 filter_type);
              if (sg_val < dp->supra_granular_val-dp->max_sg_intensity_offset 
                  /*  || sg_val > dp->supra_granular_val+dp->max_sg_intensity_offset*/)
                continue ; // not in feasible region
            }
            sg_intensity_offset = sg_val - dp->supra_granular_val ;
            
            if (sg_val <= 1.1*ig_val) // do some sanity checks
              continue ;
            
            max_csf_len = max_profile_len-(ig_len + sg_len + wm_len + 1) ;
            if (max_csf_len <= 0)
              continue ;  // enforce the presence of at least a bit of csf
            max_csf_len = MIN(max_csf_len, MAX_CSF_LEN) ;
            
            for (csf_len = MIN_CSF_LEN ; csf_len <= max_csf_len; csf_len += skip)
            {
              // compute length of subset of total profile being used
              profile_len = wm_len+ig_len+sg_len+csf_len ;
              if ((sg_val-ig_val > 2*(ig_val-wm_val)) && (ig_val-wm_val < 25) &&
                  dp->use_max_grad == 0)
                continue ;  // enforce comparable contrast
              
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
                  (*profile_error_func)(kernel,intensity_profile,profile_len,fname1,step, 
                                        offset,errors);
                }
                best_wm_dist = wm_dist ;
                best_nx = nx ; best_ny = ny ; best_nz = nz ;
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
      
      for (start_index = min_start_index ; start_index <= max_wm_index; start_index++)
      {
        offset = (current_index - start_index)*step ;
        intensity_profile = base_intensity_profile + start_index ;
        max_wm_len = max_profile_len-(min_gray_matter_len+1+start_index) ;
        max_wm_len = MIN(max_wm_len, MAX_WM_LEN) ;
        for (wm_len = min_wm_len ; wm_len <= max_wm_len; wm_len += skip)
        {
          index = wm_len+start_index ;
          if (wm_index_okay[index] == 0)
            continue ;
          wm_dist = (start_index+wm_len-current_index)*step ;   
          if (dp->fix_intensities)
            wm_val = dp->wm_val + vp->wm_intensity_offset ;
          else
          {
            wm_val = compute_optimal_intensity(intensity_profile, 0, wm_len-2,filter_type) ;
            if (wm_val < min_wm_intensity || wm_val > max_wm_intensity)
              continue ; // not in feasible region
          }
          wm_intensity_offset = wm_val - dp->wm_val ;
          for (ig_len = min_ig_len ; ig_len <= max_ig_len ; ig_len += skip)
          {
            int i1 = MAX(wm_len+2, wm_len+ig_len-(stria_len+stria_offset_len)-2) ;
            ig_val = compute_optimal_intensity(intensity_profile,wm_len+1, i1,filter_type);
            if (ig_val < dp->infra_granular_val-dp->max_ig_intensity_offset ||
                ig_val > dp->infra_granular_val+dp->max_ig_intensity_offset)
              continue ; // not in feasible region
            ig_intensity_offset = ig_val - dp->infra_granular_val ;
            
            if (ig_val <= wm_val+MIN_INTENISTY_PCT_OFFSET*wm_val) // do some sanity checks
              continue ;
            
            stria_start = wm_len+ig_len-(stria_len+stria_offset_len) ;
            stria_end = stria_start + stria_len - 1  ;
            stria_val = compute_optimal_intensity(intensity_profile,stria_start,stria_end,
                                                  filter_type);
            stria_intensity_offset = stria_val - dp->stria_val ;
            
            // do some sanity checks
            if (ig_val <= stria_val)
              continue ;
            if (stria_val < wm_val)
              continue ;
            for (index = 0 ; index < nsg ; index++)
            {
              sg_len = sg_indices[index]-(wm_len+ig_len+start_index) ;
              if (sg_len < min_sg_len || sg_len > max_sg_len)
                continue ;
              sg_val = compute_optimal_intensity(intensity_profile, 
                                                 wm_len+ig_len+1,wm_len+ig_len+sg_len-2, filter_type);
              if (sg_val < dp->supra_granular_val-dp->max_sg_intensity_offset ||
                  sg_val > dp->supra_granular_val+dp->max_sg_intensity_offset)
                continue ; // not in feasible region
              sg_intensity_offset = sg_val - dp->supra_granular_val ;
              
              if (sg_val <= 1.1*ig_val) // do some sanity checks
                continue ;
              if ((sg_val-ig_val > 2*(ig_val-wm_val))  && (ig_val-wm_val < 25))
                continue ;  // enforce comparable contrast
              
              max_csf_len = max_profile_len-(ig_len + sg_len + wm_len + 1) ;
              if (max_csf_len <= 0)
                continue ;  // enforce the presence of at least a bit of csf
              max_csf_len = MIN(max_csf_len, MAX_CSF_LEN) ;
              
              for (csf_len = MIN_CSF_LEN ; csf_len <= max_csf_len; csf_len += skip)
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
                    (*profile_error_func)(kernel,intensity_profile,profile_len,fname1,step, 
                                          offset, errors);
                  }
                  best_profile = PROFILE_V1 ;
                  best_wm_dist = wm_dist ;
                  best_nx = nx ; best_ny = ny ; best_nz = nz ;
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
  }

#if 0
  if (min_rms > 1e4)  // try two compartment model instead of 3
#endif
  {
    nx = vp->nx ; ny = vp->ny ; nz = vp->nz ;
    ig_len = 0 ;
    if (/*dp->fix_intensities == 0 && */dp->use_max_grad == 0)
    {
      min_wm_index = mark_wm_index_okay(base_intensity_profile, wm_index_okay, nint(.4/step), 
                                        min_wm_intensity, max_wm_intensity, current_index, 
                                        max_profile_len, fname4, nint(max_inwards/step));
      if (dp->dark_csf)
        mark_sg_index_okay_dark_csf(base_intensity_profile, sg_index_okay, step, 
                                    max_profile_len, fname3, 
                                    current_index, dp->outside_val+10, min_wm_index, 
                                    sg_indices, &nsg,
                                    current_index+nint(pdist/step), nint(max_outwards/step));
      else
        mark_sg_index_okay_bright_csf(base_intensity_profile, sg_index_okay, step, 
                                      max_profile_len, fname3, 
                                      current_index, dp->outside_val-10, min_wm_index, 
                                      sg_indices, &nsg,
                                      current_index+nint(pdist/step),
                                      nint(max_outwards/step));

      min_start_index = min_wm_index ;
    }
    for (max_wm_index = min_wm_index+1 ; max_wm_index < max_profile_len-min_wm_len ; max_wm_index++)
      if (wm_index_okay[max_wm_index] == 0)
        break ;
    for (start_index = min_start_index ; start_index <= max_wm_index; start_index++)
    {
      offset = (current_index - start_index)*step ;
      intensity_profile = base_intensity_profile + start_index ;
      max_wm_len = max_profile_len-(min_gray_matter_len+1+start_index) ;
      max_wm_len = MIN(max_wm_len, MAX_WM_LEN) ;
      for (wm_len = min_wm_len ; wm_len <= max_wm_len; wm_len += skip)
      {
        index = wm_len+start_index ;
        if (wm_index_okay[index] == 0)
          continue ;
        wm_dist = (start_index+wm_len-current_index)*step ;   
        if (dp->fix_intensities)
          wm_val = dp->wm_val + vp->wm_intensity_offset ;
        else
        {
          wm_val = compute_optimal_intensity(intensity_profile, 0, wm_len-2,filter_type) ;
          if (wm_val < min_wm_intensity || wm_val > max_wm_intensity)
            continue ; // not in feasible region
        }
        wm_intensity_offset = wm_val - dp->wm_val ;
        for (index = 0 ; index < nsg ; index++)
        {
          sg_len = sg_indices[index]-(wm_len+ig_len+start_index) ;
          if (sg_len < min_sg_len || sg_len > max_sg_len)
            continue ;
#if 0
          if (dp->fix_intensities)
            sg_val = dp->supra_granular_val + vp->sg_intensity_offset ;
          else
#endif
          {
            sg_val = compute_optimal_intensity(intensity_profile, wm_len+1,wm_len+sg_len-2,filter_type);
            if (sg_val < dp->infra_granular_val-dp->max_ig_intensity_offset )
              continue ; // not in feasible region
          }
          sg_intensity_offset = sg_val - dp->supra_granular_val ;
          
          max_csf_len = max_profile_len-(sg_len + wm_len + 1) ;
          if (max_csf_len <= 0)
            continue ;  // enforce the presence of at least a bit of csf
          max_csf_len = MIN(max_csf_len, MAX_CSF_LEN) ;
          
          for (csf_len = MIN_CSF_LEN ; csf_len <= max_csf_len; csf_len += skip)
          {
            // compute length of subset of total profile being used
            profile_len = wm_len+sg_len+csf_len ;
            
            construct_model_profile(kernel, dp, 
                                    wm_intensity_offset, ig_intensity_offset,sg_intensity_offset,
                                    0.0, wm_len, 0, sg_len, csf_len, PROFILE_AGRANULAR) ;
            
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
                (*profile_error_func)(kernel,intensity_profile,profile_len,fname1,step, 
                                      offset,errors);
              }
              best_wm_dist = wm_dist ;
              best_nx = nx ; best_ny = ny ; best_nz = nz ;
              best_wm_len = wm_len ;
              best_start_index = start_index ;
              best_profile = PROFILE_AGRANULAR ;
              best_csf_len = csf_len ;
              best_wm_intensity_offset = wm_intensity_offset ;
              best_sg_intensity_offset = sg_intensity_offset ;
              best_ig_len = 0 ;
              best_ig_intensity_offset = 
                (((dp->supra_granular_val+sg_intensity_offset) + 
                  (dp->wm_val + wm_intensity_offset)) / 2) -
                dp->infra_granular_val ;
              min_rms = rms ;
              best_sg_len = sg_len ; 
            }
          }
        }
      }
    }
  }
  if (dp->fix_intensities && best_wm_dist > -1000 && 0)
  {   // do linear fit to intensities in range of wm boundary to get subvoxel
    VECTOR *vY ;
    MATRIX *mX, *mP, *mXinv ;
    int    i, wm_index ;
    double border_intensity ;
#define NPOINTS 2
    wm_index = current_index + nint(best_wm_dist/step) ;
    if (wm_index > -NPOINTS)  // enough points on both sides to do a fit
    {
      mX = MatrixAlloc(2*NPOINTS+1, 2, MATRIX_REAL) ;
      vY = VectorAlloc(2*NPOINTS+1, MATRIX_REAL) ;

      for (i = -NPOINTS ; i <= NPOINTS ; i++)
      {
        *MATRIX_RELT(mX, i+NPOINTS+1, 1) = i*step ;
        *MATRIX_RELT(mX, i+NPOINTS+1, 2) = 1 ;
        *MATRIX_RELT(vY, i+NPOINTS+1, 1) = base_intensity_profile[i+wm_index] ;
      }
      mXinv = MatrixPseudoInverse(mX, NULL) ;
      mP = MatrixMultiply(mXinv, vY, NULL) ;
      border_intensity = (dp->wm_val+best_wm_intensity_offset+dp->infra_granular_val+best_ig_intensity_offset)/2;
      best_wm_dist = (border_intensity - *MATRIX_RELT(mP, 1, 2)) / *MATRIX_RELT(mP, 1, 1) ;
      MatrixFree(&mX) ; VectorFree(&vY) ; MatrixFree(&mXinv) ; MatrixFree(&mP) ;
    }
  }

  vp->nx = nx = best_nx ; vp->ny = ny = best_ny ; vp->nz = nz = best_nz ;
  vp->wtx = vp->wx + best_wm_dist * nx ;
  vp->wty = vp->wy + best_wm_dist * ny ;
  vp->wtz = vp->wz + best_wm_dist * nz ;

  if (best_profile == PROFILE_AGRANULAR)
    offset = (best_sg_len/2)*step ;
  else
    offset = (best_ig_len)*step ; 
  vp->l4tx = vp->wx + (best_wm_dist+offset) * nx ;
  vp->l4ty = vp->wy + (best_wm_dist+offset) * ny ;
  vp->l4tz = vp->wz + (best_wm_dist+offset) * nz ;
  vp->l4_dist = (best_wm_dist)+(best_ig_len*step) ;
  if ((best_profile != PROFILE_AGRANULAR) && best_wm_dist > -100)
  {
    int   i0, i1, best_index ;
    double target_intensity, min_dif  ;
    i0 = MAX(0, best_start_index + best_wm_len + best_ig_len - (PV_LEN-1)) ;
    i1 = MIN(max_profile_len-1, i0 + 2*(PV_LEN-1)) ;
    target_intensity = 
      ((dp->supra_granular_val+best_sg_intensity_offset) + 
       (dp->infra_granular_val+best_ig_intensity_offset))/2;
    best_index = i0 ;
    min_dif = fabs(base_intensity_profile[i0]-target_intensity);
    for (index = i0+1 ; index <= i1 ; index++)
      if (fabs(base_intensity_profile[index]-target_intensity) < min_dif)
      {
        best_index = index ;
        min_dif = fabs(base_intensity_profile[index]-target_intensity) ;
      }
    offset = (best_index - (best_start_index+best_wm_len))*step ;
    vp->l4tx = vp->wx + (best_wm_dist+offset) * nx ;
    vp->l4ty = vp->wy + (best_wm_dist+offset) * ny ;
    vp->l4tz = vp->wz + (best_wm_dist+offset) * nz ;
    vp->l4_dist = best_wm_dist+offset ;
  }

  offset = (best_ig_len + best_sg_len)*step ; 
  vp->ptx = vp->wx + (best_wm_dist+offset) * nx ;
  vp->pty = vp->wy + (best_wm_dist+offset) * ny ;
  vp->ptz = vp->wz + (best_wm_dist+offset) * nz ;

  vp->white_dist = best_wm_dist ;
  vp->pial_dist = (best_wm_dist)+((best_sg_len+best_ig_len)*step) ;
  
  vp->wm_intensity_offset = best_wm_intensity_offset ;
  vp->sg_intensity_offset = best_sg_intensity_offset ;
  vp->ig_intensity_offset = best_ig_intensity_offset ;
  vp->which_profile = best_profile ;
  vp->min_rms = min_rms ;
  if (min_rms > .3*dp->wm_val && best_wm_dist > -100)
    best_wm_dist = DIST_NOT_FOUND ;
  if ((deform_debug || vno == Gdiag_no) && best_wm_dist > -100)
  {
    double rms, wx, wy, wz, px, py, pz, lx, ly, lz ;
    if (/*dp->fix_intensities || */ dp->use_max_grad)
    {
      double pdist ;
      pdist = sqrt(SQR(vp->px-vp->wx) + SQR(vp->py-vp->wy) + SQR(vp->pz-vp->wz)) ;
      extract_image_intensities(mris, mri, 2*max_inwards, 2*max_outwards+pdist, vp->wx, vp->wy, vp->wz, 
                                best_nx,best_ny,best_nz, step, tsteps, base_intensity_profile, fname2);
    }
    else
      extract_image_intensities(mris, mri, max_inwards, max_outwards+pdist, vp->wx, vp->wy, vp->wz, 
                                best_nx,best_ny,best_nz,
                                step, tsteps, base_intensity_profile, fname2);
    min_wm_index = mark_wm_index_okay(base_intensity_profile, wm_index_okay, nint(.4/step), 
                                      min_wm_intensity, max_wm_intensity, current_index, 
                                      max_profile_len, fname4, nint(max_inwards/step));
    if (dp->dark_csf)
      mark_sg_index_okay_dark_csf(base_intensity_profile, sg_index_okay, step, 
                                  max_profile_len, fname3, current_index, dp->outside_val+10,
                                  min_wm_index, sg_indices, &nsg,
                                  current_index+nint(pdist/step), nint(max_outwards/step)) ;
    else
      mark_sg_index_okay_bright_csf(base_intensity_profile, sg_index_okay, step, 
                                    max_profile_len, fname3, current_index, 
                                    dp->outside_val-10,
                                    min_wm_index, sg_indices, &nsg,
                                    current_index+nint(pdist/step), nint(max_outwards/step));
    offset = (current_index - best_start_index)*step ;
    profile_len = best_wm_len + best_ig_len + best_sg_len + best_csf_len ;
    intensity_profile = base_intensity_profile + best_start_index ;
    construct_model_profile(kernel, dp, best_wm_intensity_offset, best_ig_intensity_offset, 
                            best_sg_intensity_offset, best_stria_intensity_offset, 
                            best_wm_len, best_ig_len, best_sg_len, best_csf_len, best_profile) ;
    rms = (*profile_error_func)(kernel, intensity_profile,profile_len, fname1,step,offset,
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

    MRISsurfaceRASToVoxelCached(mris, mri, vp->wtx, vp->wty, vp->wtz, &wx, &wy, &wz);
    MRISsurfaceRASToVoxelCached(mris, mri, vp->ptx, vp->pty, vp->ptz, &px, &py, &pz);
    MRISsurfaceRASToVoxelCached(mris, mri, vp->l4tx, vp->l4ty, vp->l4tz, &lx, &ly, &lz);
    printf("targets (white, IV, pial) at dists (%2.2f %2.2f %2.2f)\n",
           vp->white_dist,vp->l4_dist,vp->pial_dist) ;
    printf("W:  %2.0f %2.0f %2.0f\nL4: %2.0f %2.0f %2.0f\nP: %2.0f %2.0f %2.0f\n",
           wx, wy, wz, lx, ly, lz, px, py, pz) ;

    DiagBreak() ;  // YYY
  }
  return(best_wm_dist) ;
}

static double
compute_optimal_intensity(double *iprofile, int i0, int i1,int filter_type)
{
  double optimal_intensity, intensities[MAX_PROFILE_LEN] ;
  int    num, i ;

  switch (filter_type)
  {
  case FILTER_MEDIAN:
    for (num = 0, i = i0 ; i <= i1 ; i++, num++)
      intensities[num] = iprofile[i] ;
    qsort(intensities, num, sizeof(double), compare_doubles) ;
    optimal_intensity = intensities[num/2] ;
    break ;
  default:
  case FILTER_MEAN:
    for (optimal_intensity = 0.0, num = 0, i = i0 ; i <= i1 ; i++, num++)
      optimal_intensity += iprofile[i] ;
    optimal_intensity /= num ;
    break ;
  }
  return(optimal_intensity) ;
}

static double
compute_profile_L1(double *kernel, double *iprofile, int nsamples, 
                   const char *fname, double step, double in_dist, double *errors)
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
                   const char *fname, double step, double in_dist, double *errors)
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
compute_profile_norm_corr(double *kernel, double *iprofile, int nsamples, const char *fname,
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
  case PROFILE_AGRANULAR: // 2 compartments
    if (use_partial_volume_model)
    {
      index = inside_len ;
      if (index > 0)
        kernel[index-1]   = .75*wm_val + .25*sg_val;
      kernel[index]   = .5*wm_val + .5*sg_val;
      if (index < profile_len-1)
        kernel[index+1]   = .25*wm_val + .75*sg_val;
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
    if (which_profile == PROFILE_AGRANULAR) // 2 compartments
    {
      kernel[inside_len] = .5*wm_val + .5*sg_val ;
      if (inside_len > 0)
        kernel[inside_len-1] = .75*wm_val + .25*sg_val ;
      if (inside_len < profile_len-1)
        kernel[inside_len+1] = .25*wm_val + .75*sg_val ;
      if (out_len > 0)
      {
        index = inside_len+supra_granular_len ;
        if (index >= 0 && index < profile_len)
          kernel[index] =   .5*sg_val + .5*out_val ;
        if (index > 0)
          kernel[index-1] =   .75*sg_val + .25*out_val ;
        if (index < profile_len-1)
          kernel[index+1] =   .25*sg_val + .75*out_val ;
      }
    }
    else
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
  }
  
  return(NO_ERROR) ;
}

static double
compute_linear_combination(double *kernel, double *iprofile, int nsamples, const char *plot_fname,
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
  double       nx, ny, nz ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    vp = (VERTEX_PARMS *)v->vp ;
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
      if (vno == Gdiag_no)
        printf("v %d: wm target = %2.0f %2.0f %2.0f\n", vno, xv, yv, zv) ;
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
      if (vno == Gdiag_no)
        printf("v %d: pial target = %2.0f %2.0f %2.0f\n", vno, xv, yv, zv) ;
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
      if (vno == Gdiag_no)
        printf("v %d: layer IV target = %2.0f %2.0f %2.0f\n", vno, xv, yv, zv) ;
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
  double    dists[MAX_DISTS], mn, std, min_dist, max_dist, dist ;
  HISTOGRAM *h ;
  VERTEX_PARMS *vp ;
#if 0
  int       imax ;
  double    md ;
#endif

  memset(dists, 0, sizeof(dists)) ;
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX          const * const v  = &mris->vertices         [vno];

  vp = (VERTEX_PARMS *)(v->vp) ;
  dist = which == WHITE_VERTICES ? vp->white_dist : vp->pial_dist ; ;
  min_dist = max_dist = dist ;
  for (num = n = 0 ; n < vt->vtotal ; n++)
  {
    VERTEX const * const vn = &mris->vertices[vt->v[n]];
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
    for (n = 0 ; n < vt->vtotal ; n++)
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      if (vn->marked != 1 || vn->ripflag)
        continue ;
      vnp = (VERTEX_PARMS *)vn->vp ;
      fprintf(fp, "%d %f\n", vt->v[n], which == WHITE_VERTICES ? vnp->white_dist : vnp->pial_dist);
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
  int vno ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    double x = 0, y = 0, z = 0; 
    VERTEX* v = &mris->vertices[vno] ;
    VERTEX_PARMS *vp = (VERTEX_PARMS *)(v->vp) ;
    switch (which_src)
    {
    case SURFACE_NORMALS:x = vp->nx ;  y = vp->ny ;  z = vp->nz ; break ;
    case WHITE_VERTICES: x = vp->wx ;  y = vp->wy ;  z = vp->wz ; break ;
    case PIAL_VERTICES:  x = vp->px ;  y = vp->py ;  z = vp->pz ; break ;
    case WHITE_TARGETS:  x = vp->wtx ; y = vp->wty ; z = vp->wtz ; break ;
    case LAYERIV_VERTICES: x = vp->l4x ; y = vp->l4y ; z = vp->l4z ; break;
    case PIAL_TARGETS:   x = vp->ptx ; y = vp->pty ; z = vp->ptz ; break ;
    case LAYERIV_TARGETS:x = vp->l4tx ; y = vp->l4ty ; z = vp->l4tz ; break ;
      break ;
    default:
      cheapAssert(!"vp_copy_to_surface bad which_src");
    }

    switch (which_dst)
    {
    case SURFACE_NORMALS:  v->nx = x ; v->ny = y ; v->nz = z ; break ;
    case WHITE_VERTICES:   v->whitex = x ; v->whitey = y; v->whitez = z ;break ;
    case PIAL_VERTICES:    v->pialx = x ;  v->pialy = y;  v->pialz = z ; break ;
    case CURRENT_VERTICES: MRISsetXYZ(mris,vno, x, y, z) ;     break ;
    case TARGET_VERTICES:  v->targx = x ;  v->targy = y ; v->targz = z ; break ;
    case ORIG_VERTICES:    cheapAssert(!"vp_copy_to_surface should not use ORIG_VERTICES") ; break ;
    default:
      cheapAssert(!"vp_copy_to_surface bad which_dst") ;
   }
    v->marked = vp->found ;
  }
  return(NO_ERROR) ;
}

static int
vp_copy_from_surface(MRI_SURFACE *mris, int which_src, int which_dst)
{
  int          vno ;
  VERTEX       *v ;
  double       x = 0, y = 0, z = 0; 
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    switch (which_src)
    {
    case SURFACE_NORMALS:  x = v->nx ;      y = v->ny ;     z = v->nz ;     break ;
    case WHITE_VERTICES:   x = v->whitex  ; y = v->whitey ; z = v->whitez ; break ;
    case PIAL_VERTICES:    x = v->pialx  ;  y = v->pialy ;  z = v->pialz  ; break ;
    case CURRENT_VERTICES: x = v->x  ;      y = v->y ;      z = v->z  ;     break ;
    case TARGET_VERTICES:  x = v->targx ;   y = v->targy ;  z = v->targz ;  break ;
    case ORIG_VERTICES:    cheapAssert(!"vp_copy_from_surface should not use ORIG_VERTICES") ; break ;
    default:
      cheapAssert(!"vp_copy_from_surface bad which_src");
    }

    switch (which_dst)
    {
    case SURFACE_NORMALS:  vp->nx = x ;  vp->ny = y ;  vp->nz = z ;  break ;
    case LAYERIV_VERTICES: vp->l4x = x ; vp->l4y = y ; vp->l4z = z ; break ;
    case WHITE_VERTICES:   vp->wx = x ;  vp->wy = y ;  vp->wz = z ;  break ;
    case PIAL_VERTICES:    vp->px = x ;  vp->py = y ;  vp->pz = z ;  break ;
    case WHITE_TARGETS:    
      vp->wtx = x ; vp->wty = y ; vp->wtz = z ; 
      vp->pial_dist = sqrt(SQR(x-vp->wx)+SQR(y-vp->wy)+SQR(z-vp->wz)) ;
      break ;
    case PIAL_TARGETS:     
      vp->ptx = x ; vp->pty = y ; vp->ptz = z ; 
      vp->pial_dist = sqrt(SQR(x-vp->wx)+SQR(y-vp->wy)+SQR(z-vp->wz)) ;
      break ;
    case LAYERIV_TARGETS:  
      vp->l4tx = x ; vp->l4ty = y ; vp->l4tz = z ; 
      vp->l4_dist = sqrt(SQR(x-vp->wx)+SQR(y-vp->wy)+SQR(z-vp->wz)) ;
      break ;
    default:
      cheapAssert(!"vp_copy_from_surface bad which_dst");
      break ;
    }
  }
  return(NO_ERROR) ;
}
static int
vp_copy_to_surface_vals(MRI_SURFACE *mris, int which, DP *dp)
{
  int          vno ;
  VERTEX       *v ;
  double       val; 
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
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
    case IG_INTENSITY:        val = vp->ig_intensity_offset+dp->infra_granular_val ; break ;
    case SG_INTENSITY:        val = vp->sg_intensity_offset+dp->supra_granular_val ; break ;
    case WM_INTENSITY:        val = vp->wm_intensity_offset+dp->wm_val ; break ;
    case IG_WM_RATIO:         
      val = (vp->ig_intensity_offset+dp->infra_granular_val) / (vp->wm_intensity_offset+dp->wm_val) ; 
      break ;
    case SG_WM_RATIO:         
      val = (vp->sg_intensity_offset+dp->supra_granular_val) / (vp->wm_intensity_offset+dp->wm_val) ; 
    case SG_IG_RATIO:         
      val = (vp->sg_intensity_offset+dp->supra_granular_val) / (vp->ig_intensity_offset+dp->infra_granular_val) ; 
      break ;
    case DEEP_RATIO:         
      val = vp->deep_ratio ;
      break ;
    }
    v->marked = vp->found ;
    v->val = val ;
  }
  return(NO_ERROR) ;
}
static int
vp_copy_from_surface_vals(MRI_SURFACE *mris, int which, DP *dp)
{
  int          vno ;
  VERTEX       *v ;
  VERTEX_PARMS *vp = NULL;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    if (v->ripflag)
    {
      vp->found = 0 ;
      continue ;
    }
    switch (which)
    {
    case WHITE_VERTICES:      vp->white_val = v->val ; break ;
    case PIAL_VERTICES:       vp->pial_val = v->val ; break ;
    case LAYERIV_VERTICES:    vp->l4_val = v->val ; break ;
    case WM_INTENSITY_OFFSET: vp->wm_intensity_offset = v->val ; break ;
    case SG_INTENSITY_OFFSET: vp->sg_intensity_offset = v->val ; break ;
    case RMS:                 vp->min_rms = v->val ; break ;
    case IG_INTENSITY_OFFSET: vp->ig_intensity_offset = v->val ; break ;
    case IG_INTENSITY:        vp->ig_intensity_offset = v->val-dp->infra_granular_val ; break ;
    case SG_INTENSITY:        vp->sg_intensity_offset= v->val-dp->supra_granular_val ; break ;
    case WM_INTENSITY:        vp->wm_intensity_offset = v->val - dp->wm_val ; break ;
    case DEEP_RATIO:          vp->deep_ratio = v->val ;
    default:
    case IG_WM_RATIO:         
    case SG_WM_RATIO:         
    case SG_IG_RATIO:         
      ErrorExit(ERROR_UNSUPPORTED, "vp_copy_from_surface_vals: unsupported which %d",which);
      break ;
    }
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
    if (v->ripflag)
      continue ;
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
#if 0
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
    if (v->ripflag)
    {
      vp->found = 0 ;
      continue ;
    }
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
#endif
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

  MRIScomputeSurfaceNormals(mris, CURRENT_VERTICES, 10) ;
  MRIScomputeSurfaceNormals(mris, WHITE_VERTICES, 10) ;
  MRIScomputeSurfaceNormals(mris, PIAL_VERTICES, 10) ;
  dp->step = mri->xsize/2 ;
  compute_normals(mris, (VERTEX_PARMS *)(mris->vertices[0].vp)) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    compute_best_neighborhood_profile(mris, mri, vno, dp) ; 
  }

  compute_laminar_ratios(mris, mri, dp) ;
  dilate_not_found(mris, 3) ;
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

  double (*profile_error_func)(double *kernel, double *intensity_profile, 
                               int nsamples, const char *fname,double step, double in_dist,
                               double *errors);
  const char   *fname, *fname2 ;

  wm_val = dp->wm_val ;  // for compiler warning
  tsteps = (int)ceil(0.25/dp->step) ;

  if (vno == Gdiag_no || deform_debug)
  {
    fname = "profile.dat" ;
    fname2 = "intensity.dat" ;
  }
  else
    fname = fname2 = NULL ;

  stria_len = nint((dp->stria_width+STRIA_OFFSET) / dp->step) ;
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX          const * const v  = &mris->vertices         [vno];
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
    VERTEX const * const vn = (n < 0) ? v : &mris->vertices[vt->v[n]] ;
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
             const char *fname, DP *dp, double in_dist, double wm_intensity_offset, 
             double min_sg_intensity_offset, double max_sg_offset, double sg_intensity_offset_step, 
             int which_profile, int wm_len, int ig_len, int sg_len,   
             double (*profile_error_func)(double *kernel, double *intensity_profile, 
                                          int nsamples, const char *fname,double step, 
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

#if 0
static int
compute_normal(MRI_SURFACE *mris, int vno, VERTEX_PARMS *vp)
{
  double nx, ny, nz, norm ;
  VERTEX *v ;

  v = &mris->vertices[vno] ;
  nx = vp->px - vp->wx ; ny = vp->py - vp->wy ; nz = vp->pz - vp->wz ;
  norm = sqrt(nx*nx + ny*ny + nz*nz) ; 
  if (1 || FZERO(norm)) // they start in the same place
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
#endif

static int
filter_vals(MRI_SURFACE *mris, int filter_type) 
{
  switch (filter_type)
  {
  case FILTER_MEDIAN: MRISmedianFilterVals(mris, 1) ; break ;
  case FILTER_MEAN:   MRISaverageVals(mris, 1) ; break ;
  default: break ;
  }
  return(NO_ERROR) ;
}

static int
filter_distances(MRI_SURFACE *mris, int filter_type) 
{
  switch (filter_type)
  {
  case FILTER_MEDIAN:     MRISmedianFilterD(mris, 1, 0) ; break ;
  case FILTER_MEAN:       MRISaverageD(mris, 1) ; break ;
  case FILTER_GAUSSIAN:   MRISgaussianFilterD(mris, .8) ; break ;
  case FILTER_NONE:
  default: break ;
  }
  return(NO_ERROR) ;
}
#if 0
static int
filter_coords(MRI_SURFACE *mris, int filter_type) 
{
  switch (filter_type)
  {
  case FILTER_MEDIAN: MRISmedianFilterD(mris, 1, 0) ; break ;
  case FILTER_MEAN:   MRISaverageVertexPositions(mris, 1) ; break ;
  default: break ;
  }
  return(NO_ERROR) ;
}
#endif
static int
compute_distances_from_positions(MRI_SURFACE *mris)
{
  int           vno ;
  VERTEX        *v ;
  VERTEX_PARMS  *vp ;
  double        dx, dy, dz ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    vp->white_dist = 0 ;
    dx = vp->px - vp->wx ; dy = vp->py - vp->wy ; dz = vp->pz - vp->wz ;
    vp->pial_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    dx = vp->l4x - vp->wx ; dy = vp->l4y - vp->wy ; dz = vp->l4z - vp->wz ;
    vp->l4_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    vp->wtx = vp->wx ; vp->wty = vp->wy ; vp->wtz = vp->wz ;
    vp->ptx = vp->px ; vp->pty = vp->py ; vp->ptz = vp->pz ;
    vp->l4tx = vp->l4x ; vp->l4ty = vp->l4y ; vp->l4tz = vp->l4z ;
  }
  return(NO_ERROR) ;
}
#define WM_DIST  0.5  // how much white matter to include
static int
compute_laminar_ratios(MRI_SURFACE *mris, MRI *mri, DP *dp)
{
  int          vno, tsteps, wm_index, ig_index, sg_index, num_bottom, num_top ;
  VERTEX       *v ;
  VERTEX_PARMS *vp ;
  double       intensity_profile[MAX_PROFILE_LEN], pdist, bottom, top, peak_bottom, valley_top;
  const char         *fname ;
  int          stria_len = nint(dp->stria_width/dp->step), stria_index, index, 
    stria_offset_len = nint(STRIA_OFFSET/dp->step) ;  // amount of layerIV superficial to stria

  tsteps = (int)ceil(0.25/dp->step) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
    {
      fname = "intensity.dat" ;
      DiagBreak() ;
    }
    else
      fname = NULL ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    pdist = sqrt(SQR(vp->px-vp->wx) + SQR(vp->py-vp->wy) + SQR(vp->pz-vp->wz)) ;
    if (extract_image_intensities(mris, mri, WM_DIST, pdist, vp->wx, vp->wy, vp->wz, 
                                  vp->nx, vp->ny, vp->nz, dp->step, tsteps, 
                                  intensity_profile, fname)!= NO_ERROR)
    {
      vp->found = 0;
      vp->deep_ratio = 0 ;
      continue ;
    }
    wm_index = nint(WM_DIST/dp->step) ;
    ig_index = nint((WM_DIST + vp->l4_dist)/dp->step) ;
    sg_index = nint((WM_DIST + vp->pial_dist)/dp->step) ;
    stria_index = ig_index - (stria_len+stria_offset_len) ;
    peak_bottom = intensity_profile[wm_index] ;
    for (num_bottom  = 0, bottom = 0.0, index = wm_index ; index < stria_index ; 
         index++, num_bottom++)
    {
      bottom += intensity_profile[index] ;
      if (intensity_profile[index] > peak_bottom)
        peak_bottom = intensity_profile[index] ;
    }
    if (num_bottom > 0)
      bottom /= num_bottom ;
    valley_top = intensity_profile[stria_index] ;
    for (num_top = 0, top = 0.0, index = stria_index  ; 
         index <= ig_index ; 
         index++, num_top++)
    {
      top += intensity_profile[index] ;
      if (intensity_profile[index] < valley_top)
        valley_top = intensity_profile[index];
    }
    vp->deep_ratio = peak_bottom / valley_top ;
    if (num_top > 0)
      top /= num_top ;
      
#if 0
    if (num_top > 0 && num_bottom >0)
      vp->deep_ratio = bottom / top ;
    if (FZERO(stria_index-wm_index))
    {
      vp->found = 0 ;
      continue ;
    }
    if (FZERO(ig_index-stria_index))
    {
      vp->found = 0 ;
      continue ;
    }
    top /= (ig_index - stria_index) ;
    if (FZERO(top))
    {
      vp->found = 0 ;
      continue ;
    }
    vp->deep_ratio = bottom / top ;
#endif
    vp->found = 1 ;
  }
  return(NO_ERROR) ;
}
static int
dilate_not_found(MRI_SURFACE *mris, int niter) 
{
  int          vno ;
  VERTEX       *v ;
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    v->marked  = vp->found ;
  }
  MRISerodeMarked(mris, niter) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    vp->found = v->marked ;
  }
  return(NO_ERROR) ;
}
static int
vp_copy_found_from_marked(MRI_SURFACE *mris)
{
  int          vno ;
  VERTEX       *v ;
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    vp->found = v->marked ;
  }
  return(NO_ERROR) ;
}

#if 0
static int
vp_copy_found_to_marked(MRI_SURFACE *mris)
{
  int          vno ;
  VERTEX       *v ;
  VERTEX_PARMS *vp ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    vp = (VERTEX_PARMS *)(v->vp) ;
    v->marked = vp->found ;
  }
  return(NO_ERROR) ;
}
#endif
static int
compute_normals(MRI_SURFACE *mris, VERTEX_PARMS *vp)
{
  INTEGRATION_PARMS thick_parms ;
  static MRI_SURFACE       *mris_ico = NULL ;
  char              fname[STRLEN], *cp ;
  int               vno ;
  VERTEX            *v ;
  double            nx, ny, nz, norm ;
  static int cno  = 0 ;

  vp_copy_from_surface(mris, SURFACE_NORMALS, SURFACE_NORMALS) ;
  //  return(NO_ERROR) ;

  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  memset(&thick_parms, 0, sizeof(thick_parms)) ; 
  if (cno++ > 0)
    thick_parms.dt = 0.01 ; 
  else
    thick_parms.dt = 0.01 ; 
  thick_parms.momentum = 0; thick_parms.l_nlarea = 0 ;
  thick_parms.tol = 2e-5 ; 

  thick_parms.l_thick_normal = 1.0 ;
  thick_parms.l_thick_parallel = 1.0 ;
  thick_parms.niterations = 1000 ;
  thick_parms.l_thick_min = .1 ;
  if (mris_ico == NULL)
  {
    cp = getenv("FREESURFER_HOME");
    if (cp == NULL)
      ErrorExit(ERROR_BADPARM, "%s: FREESURFER_HOME not defined in environment", cp) ;
    sprintf(fname,"%s/lib/bem/ic7.tri",cp);
    mris_ico = MRISread(fname) ;
    if (!mris_ico)
      ErrorExit(ERROR_NOFILE, "%s: could not open surface file %s",Progname, fname) ;
    MRISscaleBrain(mris_ico, mris_ico, mris->radius/mris_ico->radius) ;
    MRISsaveVertexPositions(mris_ico, CANONICAL_VERTICES) ;
  }
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISstoreRipFlags(mris) ; MRISunrip(mris) ;
  MRIScomputeNormals(mris) ;
  MRISminimizeThicknessFunctional(mris, &thick_parms, 10.0) ;
  MRISrestoreRipFlags(mris) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nx = v->tx - v->whitex ; ny = v->ty - v->whitey ; nz = v->tz - v->whitez ;
    norm = sqrt(nx*nx + ny*ny + nz*nz) ;
    if (norm > 0)
    {
      nx /= norm ; ny /= norm ; nz /= norm ;
      vp[vno].nx = nx ; vp[vno].ny = ny ; vp[vno].nz = nz ;
      if (vno == Gdiag_no)
        printf("setting normal for vno %d to (%2.2f, %2.2f, %2.2f)\n", vno, nx, ny, nz) ;
    }
  }
  vp_copy_to_surface(mris, WHITE_VERTICES, CURRENT_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  return(NO_ERROR) ;
}

static int
rip_dark_vertices(MRI_SURFACE *mris, MRI *mri, float threshold)
{
  int    vno, ripped ;
  double xv, yv, zv ;
  VERTEX *v ;
  float  val ;

  for (ripped = 0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    MRISvertexToVoxel(mris, v, mri, &xv, &yv, &zv) ;
    val = MRIgetVoxVal(mri, nint(xv), nint(yv), nint(zv), 0) ;
    if (val < threshold)
    {
      ripped++ ;
      v->ripflag = 1 ;
    }
  }

  printf("%d dark vertices ripped, %d remain valid\n", ripped, MRISvalidVertices(mris));
  return(NO_ERROR) ;
}

static int
rip_vertices_with_no_faces(MRI_SURFACE *mris)
{
  int    vno, n, fnum ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    for (n = fnum = 0; n < vt->num ; n++)
      if (mris->faces[vt->f[n]].ripflag == 0)
        fnum++ ;
    if (fnum == 0)
      v->ripflag = 1 ;
  }
  return(NO_ERROR) ;
}


