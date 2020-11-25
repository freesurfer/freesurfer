/*
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
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "mrimorph.h"
#include "mrinorm.h"
#include "mri_conform.h"
#include "cvector.h"
#include "histo.h"
#include "version.h"


int main(int argc, char *argv[]) ;

#define TOO_SMALL(dist)     (dist < 1)
#define MAX_FLASH_VOLUMES   50
#define ORIG_EXPANSION_DIST 1.0  /* average thickness of the cortex in mm (approx) */
#define MAX_SAMPLES         1000
static int MAX_VNO =        5000000 ;
static int MIN_VNO = 0 ;

#define EPSILON             0.25

#define IS_CSF(T1, PD)  (!((PD > MIN_NONBRAIN_PD) && (T1 > MIN_NONBRAIN_T1 && T1 < MIN_CSF_T1)))
#define IS_PV_CSF(T1,PD) (T1 > MIN_PARTIAL_VOLUMED_CSF_T1)
#if 0
#define IS_WM(T1,PD)   ((T1 >= MIN_WM_T1) && (T1 <= MAX_WM_T1) && (PD <= MAX_WM_PD) && (PD >= MIN_WM_PD))
#define IS_GM(T1,PD)   ((T1 >= MIN_GM_T1) && (T1 <= MAX_GM_T1) && (PD >= MIN_GM_PD))
#define IS_BRAIN(T1,PD)  (IS_GM(T1,PD) || IS_WM(T1,PD))

#else

#define IS_WM(T1,PD,vno,ep)   ((T1 >= ep->cv_min_wm_T1[vno]) && (T1 <= ep->cv_max_wm_T1[vno]) && (PD <= ep->cv_max_wm_PD[vno]) && (PD >= ep->cv_min_wm_PD[vno]))
#define IS_GM(T1,PD,vno,ep)   ((T1 >= ep->cv_min_gm_T1[vno]) && (T1 <= ep->cv_max_gm_T1[vno]) && (PD >= ep->cv_min_gm_PD[vno]))
#define IS_BRAIN(T1,PD,vno,ep)  (IS_GM(T1,PD,vno,ep) || IS_WM(T1,PD,vno,ep))

#endif

#define WM_PARM_DIST(T1,PD)   (T1 > MAX_WM_T1 ? (T1-MAX_WM_T1) : T1 < MIN_WM_T1 ? MIN_WM_T1-T1 : 0)
#define CSF_PARM_DIST(T1,PD)   (fabs(MIN_CSF_T1-T1))

typedef struct {
  float *cv_wm_T1 ;
  float *cv_wm_PD ;
  float *cv_gm_T1 ;
  float *cv_gm_PD ;
  float *cv_csf_T1 ;
  float *cv_csf_PD ;
  float *cv_inward_dists ;
  float *cv_outward_dists ;
  float *cv_last_pialx ;
  float *cv_last_pialy ;
  float *cv_last_pialz ;
  float *cv_last_whitex ;
  float *cv_last_whitey ;
  float *cv_last_whitez ;
  float *cv_min_wm_PD ;
  float *cv_max_wm_PD ;
  float *cv_min_gm_PD ;
  float *cv_max_gm_PD ;
  float *cv_min_wm_T1 ;
  float *cv_max_wm_T1 ;
  float *cv_min_gm_T1 ;
  float *cv_max_gm_T1 ;
  MRI   **mri_flash ;
  int   nvolumes ;
  int   *nearest_pial_vertices ;
  int   *nearest_white_vertices ;
  double current_sigma ;
  double dstep ;
  double max_inward_dist ;
  double max_outward_dist ;
  int    nvertices ;   /* # of vertices in surface on previous invocation */
  double scale ;       /* PD scaling */
  MRI    *mri_T1 ;
  MRI    *mri_PD ;
}
EXTRA_PARMS ;


static int   compute_PD_T1_limits(MRI_SURFACE *mris, EXTRA_PARMS *ep, int navgs) ;
static int   plot_stuff = 0 ;
static int min_max_scale = 100 ;
#if 0
static float histo_sigma = 1.0f ;
#endif
static int fix_T1 = 0 ;
static int  build_lookup_table(double tr, double flip_angle,
                               double min_T1, double max_T1, double step) ;

static double    scale_all_images(MRI **mri_flash, int nvolumes, MRI_SURFACE *mris,
                                  float target_pd_wm, EXTRA_PARMS *ep) ;
static int    compute_T1_PD(double *image_vals, MRI **mri_flash, int nvolumes,
                            double *pT1, double *pPD);
#if 0
static int    compute_T1_PD_slow(double *image_vals, MRI **mri_flash, int nvolumes,
                                 double *pT1, double *pPD);
#endif
static double lookup_flash_value(double TR, double flip_angle, double PD, double T1) ;
#if 0
static float subsample_dist = 10.0 ;
#endif

#define PD_STD          125
#define MIN_T1          10
#define MAX_T1          10000
#define T1_STEP_SIZE    5

#define MIN_WM_T1  500
#define MAX_WM_T1  1200
#define MIN_GM_T1  900
#define MAX_GM_T1  1600
#if 1
#define MIN_WM_PD  500
#define MIN_GM_PD  400
#define MAX_GM_PD  1500
#endif
#define MAX_WM_PD  1300
#define MIN_CSF_T1 1900 /*1600*/
#define MIN_PARTIAL_VOLUMED_CSF_T1  1500
#define MIN_CSF_PD 1500
#define MEAN_WM_PD 800

static int sample_type = SAMPLE_NEAREST ;
static double MIN_NONBRAIN_T1 = 500.0 ;
static double MIN_NONBRAIN_PD = 550 ;
static double MIN_RELIABLE_PD = 150 ;

//#define BRIGHT_LABEL         130
//#define BRIGHT_BORDER_LABEL  100
static double  DSTEP = 0.25 ;
static double MAX_DSTEP = 0.5 ;   /* max sampling distance */
#define MAX_WHITE_DIST  8               /* max distance to sample in and out */
#define MAX_PIAL_DIST   8
#define DEFORM_WHITE         0
#define DEFORM_PIAL          1

static int       compute_parameter_maps(MRI **mri_flash, int nvolumes, MRI **pmri_T1,
                                        MRI **pmri_PD) ;
/*static int init_lookup_table(MRI *mri) ;*/
static int       store_current_positions(MRI_SURFACE *mris, EXTRA_PARMS *parms) ;
static int       update_distances(MRI_SURFACE *mris, EXTRA_PARMS *parms) ;
static int       update_parameters(MRI_SURFACE *mris, EXTRA_PARMS *parms) ;
#if 0
static double compute_sse(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double vertex_error(MRI_SURFACE *mris, int vno, EXTRA_PARMS *ep, double *prms) ;
#endif
static double compute_optimal_parameters(MRI_SURFACE *mris, int vno, EXTRA_PARMS *ep,
    double *pwhite_delta, double *pial_delta) ;
static double compute_optimal_vertex_positions(MRI_SURFACE *mris, int vno,
    EXTRA_PARMS *ep, double *pwhite_delta,
    double *ppial_delta,
    int debug_flag) ;
static double image_forward_model(double TR, double flip_angle, double dist,
                                  double thickness, double T1_wm, double PD_wm,
                                  double T1_gm, double PD_gm, double T1_csf, double PD_csf) ;
static int write_map(MRI_SURFACE *mris, float *cv, const char *name, int suffix, const char *output_suffix) ;
static int write_maps(MRI_SURFACE *mris, EXTRA_PARMS *ep, int suffix, const char *output_suffix) ;
static int smooth_maps(MRI_SURFACE *mris, EXTRA_PARMS *ep, int smooth_parms) ;
static int smooth_marked_maps(MRI_SURFACE *mris, EXTRA_PARMS *ep, int smooth_parms) ;
static int smooth_map(MRI_SURFACE *mris, float *cv, int navgs) ;
static int smooth_csf_map(MRI_SURFACE *mris, float *cv_T1, float *cv_PD, int navgs) ;
static int smooth_marked_map(MRI_SURFACE *mris, float *cv, int navgs) ;
static int smooth_marked_csf_map(MRI_SURFACE *mris, float *cv_T1, float *cv_PD, int navgs) ;
static int constrain_parameters(int nvertices, EXTRA_PARMS *ep) ;
#if 0
static int soap_bubble_maps(MRI_SURFACE *mris, EXTRA_PARMS *ep, int smooth_parms) ;
static int soap_bubble_map(MRI_SURFACE *mris, float *cv, int navgs) ;
#endif
static int compute_maximal_distances(MRI_SURFACE *mris, float sigma, MRI **mri_flash,
                                     int nvolumes, float *cv_inward_dists,
                                     float *cv_outward_dists, int *nearest_pial_vertices,
                                     int *nearest_white_vertices, double dstep,
                                     double max_inward_dist, double max_outward_dist,
                                     MRI *mri_T1, MRI *mri_PD, EXTRA_PARMS *ep) ;
static int find_nearest_pial_vertices(MRI_SURFACE *mris, int *nearest_pial_vertices,
                                      int  *nearest_white_vertices) ;
static int find_nearest_white_vertices(MRI_SURFACE *mris, int *nearest_white_vertices) ;
static int sample_parameter_maps(MRI_SURFACE *mris, MRI *mri_T1, MRI *mri_PD,
                                 float *cv_wm_T1, float *cv_wm_PD,
                                 float *cv_gm_T1, float *cv_gm_PD,
                                 float *cv_csf_T1, float *cv_csf_PD,
                                 float *cv_inward_dists, float *cv_outward_dists,
                                 EXTRA_PARMS *ep, int fix_T1, INTEGRATION_PARMS *parms,
                                 int start_vno) ;
#if 0
static int sample_parameter_map(MRI_SURFACE *mris, MRI *mri, MRI *mri_res,
                                float *cv_parm, float *cv_dists, float dir,
                                int which_vertices,
                                char *name, HISTOGRAM *prior_histo,
                                double dstep,
                                double max_dist) ;
#endif

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

MRI *MRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm) ;

static char brain_name[STRLEN] = "brain" ;
#if 0
static double dM_dT1(double flip_angle, double TR, double PD, double T1) ;
static double dM_dPD(double flip_angle, double TR, double PD, double T1) ;
#endif
static double FLASHforwardModel(double flip_angle, double TR, double PD,
                                double T1) ;
static double FLASHforwardModelLookup(double flip_angle, double TR, double PD,
                                      double T1) ;

static double compute_vertex_sse(EXTRA_PARMS *ep, double image_vals[MAX_FLASH_VOLUMES][MAX_SAMPLES],
                                 int max_j,
                                 double white_dist, double cortical_dist, double T1_wm, double PD_wm,
                                 double T1_gm, double PD_gm, double T1_csf, double PD_csf, int debug,
                                 double *T1_vals, double *PD_vals, int vno) ;
const char *Progname ;
static char *gSdir = NULL ;

static int graymid = 0 ;
static int curvature_avgs = 10 ;
static int create = 1 ;
static int smoothwm = 0 ;
static int white_only = 0 ;

static int apply_median_filter = 0 ;

static int nbhd_size = 5 ;

static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE    1.0
static float base_dt_scale = BASE_DT_SCALE ;

static int add = 0 ;

static int orig_flag = 1 ;
static char *start_white_name = NULL ;
static char *start_pial_name = NULL ;
static int smooth_parms = 10 ;
static int smooth = 0 ;
static int vavgs = 0 ;

static int nbrs = 2 ;
static int write_vals = 0 ;

static const char *suffix = "" ;
static const char *output_suffix = "ms" ;
static char *xform_fname = NULL ;

static char pial_name[STRLEN] = "pial" ;
static char white_matter_name[STRLEN] = WHITE_MATTER_NAME ;
static char orig_name[STRLEN] = ORIG_NAME ;

static int lh_label = LH_LABEL ;
static int rh_label = RH_LABEL ;
#if 0
static int max_averages = 0 ;
static int min_averages = 0 ;
static float sigma = 0.0f ;
#else
static int max_averages = 8 ;
static int min_averages = 2 ;
static float sigma = 2.0f ;
#endif
static float max_thickness = 5.0 ;

static const char *map_dir = "parameter_maps" ;

static int ms_errfunc_rip_vertices(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double ms_errfunc_gradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double ms_errfunc_sse(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static int ms_errfunc_timestep(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double ms_errfunc_rms(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;

static float *cv_inward_dists, *cv_outward_dists ;
static float *cv_wm_T1, *cv_wm_PD, *cv_gm_T1, *cv_gm_PD, *cv_csf_T1, *cv_csf_PD ;
static MRI *mri_T1 = NULL, *mri_PD = NULL ;

int
main(int argc, char *argv[]) {
  char          **av, *hemi, *sname, sdir[STRLEN], *cp, fname[STRLEN], mdir[STRLEN],
  *xform_fname ;
  int           ac, nargs, i, /*label_val, replace_val,*/ msec, n_averages, nvolumes,
  index, j, start_t, replace_val, label_val ;
  double        current_sigma ;
  MRI_SURFACE   *mris ;
  MRI           *mri_template = NULL, *mri_filled,
                                /* *mri_labeled ,*/ *mri_flash[MAX_FLASH_VOLUMES] ;
  float         max_len ;
  Timer then ;
  LTA           *lta ;
  EXTRA_PARMS   ep ;

  nargs = handleVersionOption(argc, argv, "mris_ms_refine");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.projection = NO_PROJECTION ;
  parms.tol = 1e-3 ;
  parms.dt = 0.5f ;
  parms.base_dt = BASE_DT_SCALE*parms.dt ;
  parms.l_spring = 0.0f ;
  parms.l_curv = .1 ;
  parms.l_intensity = 0.0 ;
  parms.l_tspring = 1.0f ;
  parms.l_nspring = 0.5 ;
  parms.l_repulse = 1 ;
  parms.l_surf_repulse = 5 ;
  parms.l_external = 1 ;

  parms.niterations = 100 ;
  parms.write_iterations = 0 /*WRITE_ITERATIONS */;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = 0.0 /*0.8*/ ;
  /*  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;*/
  parms.user_parms = (void *)&ep ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    usage_exit() ;

  /* set default parameters for white and gray matter surfaces */
  if (parms.momentum < 0.0)
    parms.momentum = 0.0 /*0.75*/ ;

  gMRISexternalGradient = ms_errfunc_gradient ;
  gMRISexternalSSE = ms_errfunc_sse ;
  gMRISexternalRMS = ms_errfunc_rms ;
  gMRISexternalTimestep = ms_errfunc_timestep ;
  gMRISexternalRipVertices = ms_errfunc_rip_vertices ;

  then.reset() ;
  sname = argv[1] ;
  hemi = argv[2] ;

  if (gSdir)
    strcpy(sdir, gSdir) ;
  else {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  xform_fname = argv[3] ;
  int req = snprintf(fname, STRLEN, "%s/%s/mri/flash/%s/%s", sdir, sname, map_dir, xform_fname) ; 
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  lta = LTAread(fname) ;
  if (!lta)
    ErrorExit(ERROR_NOFILE, "%s: could not read FLASH transform from %s...\n",
              Progname, fname) ;

  cp = getenv("FREESURFER_HOME") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM,
              "%s: FREESURFER_HOME not defined in environment.\n", Progname) ;
  strcpy(mdir, cp) ;


#define FLASH_START 4
  nvolumes = argc - FLASH_START ;
  printf("reading %d flash volumes...\n", nvolumes) ;
  for (i = FLASH_START ; i < argc ; i++) {
    MRI *mri ;

    index = i-FLASH_START ;
    int req = snprintf(fname, STRLEN, "%s/%s/mri/flash/%s/%s", sdir, sname, map_dir, argv[i]) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("reading FLASH volume %s...\n", fname) ;
    mri = mri_flash[index] = MRIread(fname) ;
    if (!mri_flash[index])
      ErrorExit(ERROR_NOFILE, "%s: could not read FLASH volume from %s...\n",
                Progname, fname) ;
    /*    init_lookup_table(mri_flash[index]) ;*/

    if (!mri_template) {
      mri = MRIcopy(mri_flash[index], NULL) ;
      mri_template = MRIconform(mri) ;
      mri = mri_flash[index] ;
    }

    switch (sample_type) {
    case SAMPLE_TRILINEAR:
      mri_flash[index] = MRIresample(mri, mri_template, SAMPLE_TRILINEAR) ;
      break ;
    case SAMPLE_SINC:
      mri_flash[index] = MRIresample(mri, mri_template, SAMPLE_SINC) ;
      break ;
    case SAMPLE_NEAREST:
      mri_flash[index] = MRIresample(mri, mri_template, SAMPLE_NEAREST) ;
      break ;
    default:
      ErrorExit(ERROR_UNSUPPORTED, "%s: unsupported sample type %d",
                Progname, sample_type) ;
      break ;
    }

    mri_flash[index]->tr = mri->tr ;
    mri_flash[index]->flip_angle = mri->flip_angle ;
    mri_flash[index]->te = mri->te ;
    mri_flash[index]->ti = mri->ti ;
    MRIfree(&mri) ;
    mri = mri_flash[index] ;
    printf("TR = %2.1f msec, flip angle = %2.0f degrees, TE = %2.1f msec\n",
           mri->tr, DEGREES(mri->flip_angle), mri->te) ;
    build_lookup_table(mri->tr, mri->flip_angle, MIN_T1, MAX_T1, T1_STEP_SIZE) ;
  }

  index = i-FLASH_START ;

  if (mri_T1 || mri_PD) {
    MRI    *mri_tmp ;

    if (mri_T1) {
      mri_tmp = MRIresample(mri_T1, mri_template, SAMPLE_TRILINEAR) ;
      MRIfree(&mri_T1) ;
      mri_T1 = mri_tmp ;
    }
    if (mri_PD) {
      mri_tmp = MRIresample(mri_PD, mri_template, SAMPLE_TRILINEAR) ;
      MRIfree(&mri_PD) ;
      mri_PD = mri_tmp ;
    }
  }
  MRIfree(&mri_template) ;

  req = snprintf(fname, STRLEN, "%s/%s/mri/filled", sdir, sname) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_filled = MRIread(fname) ;
  if (!mri_filled)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

  if (!stricmp(hemi, "lh")) {
    label_val = lh_label ;
    replace_val = rh_label ;
  } else {
    label_val = rh_label ;
    replace_val = lh_label ;
  }
  /* remove other hemi */
  MRIdilateLabel(mri_filled, mri_filled, replace_val, 1) ;
  if (replace_val == RH_LABEL)
    MRIdilateLabel(mri_filled, mri_filled, RH_LABEL2, 1) ;

  for (i = 0 ; i < nvolumes ; i++) {
    if (replace_val == RH_LABEL)
      MRImask(mri_flash[i], mri_filled, mri_flash[i], RH_LABEL2,0) ;

    MRImask(mri_flash[i], mri_filled, mri_flash[i], replace_val,0) ;
  }
  MRIfree(&mri_filled) ;

  if (orig_flag) {
    int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s", sdir, sname, hemi, orig_name, suffix) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("reading orig surface position from %s...\n", fname) ;
  } else {
    if (start_pial_name) {
      int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s", sdir, sname, hemi,
			 start_pial_name, suffix) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s", sdir, sname, hemi, pial_name,
			 suffix) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
  }

  if (add)
    mris = MRISreadOverAlloc(fname, 1.5) ;
  else
    mris = MRISreadOverAlloc(fname, 1.0) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  if (add)
    printf("surface read - max vertices %d\n", mris->max_vertices) ;

  MRISallocExtraGradients(mris) ;

  sprintf(parms.base_name, "%s%s", output_suffix, suffix) ;
  if (orig_flag) {
    MRISaverageVertexPositions(mris, 2) ;  /* so normals will be reasonable */
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    MRISstoreMetricProperties(mris) ;
    parms.start_t = 0 ;
  } else {
    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
    if (start_white_name) {
      if (MRISreadOriginalProperties(mris, start_white_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE,
                  "%s: could not read gray/white surface from %s",
                  Progname, start_white_name) ;
    } else {
      if (MRISreadOriginalProperties(mris, white_matter_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE,
                  "%s: could not read gray/white surface from %s",
                  Progname, white_matter_name) ;
    }
  }

  ep.mri_flash = mri_flash ;
  ep.nvolumes = nvolumes ;
  ep.cv_inward_dists = cv_inward_dists = cvector_alloc(mris->max_vertices) ;
  ep.cv_outward_dists = cv_outward_dists = cvector_alloc(mris->max_vertices) ;

  ep.cv_last_pialx = cvector_alloc(mris->max_vertices) ;
  ep.cv_last_pialy = cvector_alloc(mris->max_vertices) ;
  ep.cv_last_pialz = cvector_alloc(mris->max_vertices) ;

  ep.cv_last_whitex = cvector_alloc(mris->max_vertices) ;
  ep.cv_last_whitey = cvector_alloc(mris->max_vertices) ;
  ep.cv_last_whitez = cvector_alloc(mris->max_vertices) ;

  ep.cv_min_wm_PD = cvector_alloc(mris->max_vertices) ;
  ep.cv_max_wm_PD = cvector_alloc(mris->max_vertices) ;
  ep.cv_min_gm_PD = cvector_alloc(mris->max_vertices) ;
  ep.cv_max_gm_PD = cvector_alloc(mris->max_vertices) ;

  ep.cv_min_wm_T1 = cvector_alloc(mris->max_vertices) ;
  ep.cv_max_wm_T1 = cvector_alloc(mris->max_vertices) ;
  ep.cv_min_gm_T1 = cvector_alloc(mris->max_vertices) ;
  ep.cv_max_gm_T1 = cvector_alloc(mris->max_vertices) ;

  ep.cv_wm_T1 = cv_wm_T1 = cvector_alloc(mris->max_vertices) ;
  ep.cv_wm_PD = cv_wm_PD = cvector_alloc(mris->max_vertices) ;
  ep.cv_gm_T1 = cv_gm_T1 = cvector_alloc(mris->max_vertices) ;
  ep.cv_gm_PD = cv_gm_PD = cvector_alloc(mris->max_vertices) ;
  ep.cv_csf_T1 = cv_csf_T1 = cvector_alloc(mris->max_vertices) ;
  ep.cv_csf_PD = cv_csf_PD = cvector_alloc(mris->max_vertices) ;
  ep.nearest_pial_vertices = (int *)calloc(mris->max_vertices, sizeof(int)) ;
  ep.nearest_white_vertices = (int *)calloc(mris->max_vertices, sizeof(int)) ;
  if (!cv_inward_dists || !cv_outward_dists || !cv_wm_T1 || !cv_wm_PD
      || !cv_gm_T1 || !cv_gm_PD || !cv_csf_T1 || !cv_csf_PD || !ep.nearest_pial_vertices
      || !ep.nearest_white_vertices || !ep.cv_min_wm_PD  || !ep.cv_max_wm_PD  || !ep.cv_min_gm_PD
      || !ep.cv_max_gm_PD)
    ErrorExit(ERROR_NOMEMORY, "%s: could allocate %d len cvector arrays",
              Progname, mris->nvertices) ;
  ep.nvertices = mris->nvertices ;
  ep.dstep = MAX_DSTEP ;
  if (orig_flag) {
    ep.max_outward_dist = 1.5*MAX_PIAL_DIST ;
    ep.max_inward_dist = 1.5*MAX_WHITE_DIST ;
  } else {
    ep.max_outward_dist = MAX_PIAL_DIST ;
    ep.max_inward_dist = MAX_WHITE_DIST ;
  }
  ep.scale = scale_all_images(mri_flash, nvolumes, mris, MEAN_WM_PD, &ep) ;
  if (smooth) {
    printf("smoothing surface for %d iterations...\n", smooth) ;
    MRISaverageVertexPositions(mris, smooth) ;
  }
  if (Gdiag & DIAG_WRITE)
    write_maps(mris, &ep, 0, output_suffix) ;

  if (nbrs > 1)
    MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;

  MRIScomputeMetricProperties(mris) ;    /* recompute surface normals */

  if (add) {
    printf("adding vertices to initial tessellation...\n") ;
    for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
    while (MRISdivideLongEdges(mris, max_len) > 0) {}
  }

  if (!mri_T1) {
    printf("computing parameter maps...\n") ;
    compute_parameter_maps(mri_flash, nvolumes, &mri_T1, &mri_PD) ;
    printf("done.\n") ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      printf("writing parameter maps to file...\n") ;
      MRIwrite(mri_T1, "T1.mgh") ;
      MRIwrite(mri_PD, "PD.mgh") ;
    }
  }
  ep.mri_T1 = mri_T1 ;
  ep.mri_PD = mri_PD ;

  current_sigma = sigma ;
  for (n_averages = max_averages, i = 0 ;
       n_averages >= min_averages ;
       n_averages /= 2, current_sigma /= 2, i++) {
    ep.current_sigma = current_sigma ;
    compute_PD_T1_limits(mris, &ep, min_max_scale) ;
    printf("computing inner and outer bounds on error functional, avgs=%d,sigma=%2.2f,dstep=%2.3f,"
           "max_dist=%2.2f:%2.2f...\n",
           n_averages, current_sigma, ep.dstep, ep.max_inward_dist, ep.max_outward_dist) ;
    {
      for (j = 0 ; j < mris->nvertices ; j++)
        ep.nearest_white_vertices[j] = -1 ;
    }

    start_t = parms.start_t ;
    find_nearest_pial_vertices(mris, ep.nearest_pial_vertices, ep.nearest_white_vertices);
    find_nearest_white_vertices(mris, ep.nearest_white_vertices) ;
    compute_maximal_distances(mris, current_sigma, mri_flash, nvolumes,
                              cv_inward_dists, cv_outward_dists,
                              ep.nearest_pial_vertices, ep.nearest_white_vertices,
                              ep.dstep, ep.max_inward_dist, ep.max_outward_dist,
                              ep.mri_T1, ep.mri_PD, &ep) ;
#if 0
    if (orig_flag && start_t == 0) {
      for (i = 0 ; i < mris->nvertices ; i++)
        cv_outward_dists[i] = 4.0 ;   /* first time through allow large-scale search */
    }
#endif
    printf("estimating tissue parameters for gray/white/csf...\n") ;
    sample_parameter_maps(mris, mri_T1, mri_PD, cv_wm_T1, cv_wm_PD,
                          cv_gm_T1, cv_gm_PD, cv_csf_T1, cv_csf_PD,
                          cv_inward_dists, cv_outward_dists, &ep, fix_T1, &parms, MIN_VNO) ;
    fix_T1 = 0 ; /* only on 1st time through */
    smooth_maps(mris, &ep, smooth_parms) ;
    constrain_parameters(mris->nvertices, &ep) ;
    if (Gdiag & DIAG_WRITE)
      write_maps(mris, &ep, n_averages, output_suffix) ;
    if (Gdiag_no >= 0 && Gdiag_no < mris->nvertices) {
      printf("v %d: wm = (%2.0f, %2.0f), gm = (%2.0f, %2.0f), csf = (%2.0f, %2.0f)\n",
             Gdiag_no,ep.cv_wm_T1[Gdiag_no], ep.scale*ep.cv_wm_PD[Gdiag_no],
             ep.cv_gm_T1[Gdiag_no], ep.scale*ep.cv_gm_PD[Gdiag_no],
             ep.cv_csf_T1[Gdiag_no], ep.scale*ep.cv_csf_PD[Gdiag_no]);
      DiagBreak() ;
    }

    parms.sigma = current_sigma ;
    parms.n_averages = n_averages ;
    MRISprintTessellationStats(mris, stderr) ;

    if (write_vals) {
      sprintf(fname, "./%s-white%2.2f.w", hemi, current_sigma) ;
      MRISwriteValues(mris, fname) ;
    }
    MRISclearD(mris) ;
    MRISpositionSurfaces(mris, mri_flash, nvolumes, &parms);
    if (add) {
      for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
      while (MRISdivideLongEdges(mris, max_len) > 0) {}
    }
    if (!n_averages)
      break ;
#if 0
    ep.dstep /= 2 ;
    if (ep.dstep < 0.2)
      ep.dstep = 0.2 ;  /* no point in sampling finer (and takes forever) */
#endif
    ep.max_outward_dist = MAX_PIAL_DIST ;
  }

  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s%s",
		 sdir, sname,hemi,white_matter_name,
		 output_suffix,suffix);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("writing white matter surface to %s...\n", fname) ;
  MRISaverageVertexPositions(mris, smoothwm) ;
  MRISwrite(mris, fname) ;

  if (create)   /* write out curvature and area files */
  {
    MRIScomputeMetricProperties(mris) ;
    MRIScomputeSecondFundamentalForm(mris) ;
    MRISuseMeanCurvature(mris) ;
    MRISaverageCurvatures(mris, curvature_avgs) ;
    int req = snprintf(fname, STRLEN, "%s.curv%s%s",
		       mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", output_suffix,
		       suffix);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("writing smoothed curvature to %s\n", fname) ;
    MRISwriteCurvature(mris, fname) ;
    req = snprintf(fname, STRLEN, "%s.area%s",
		   mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", suffix);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRISprintTessellationStats(mris, stderr) ;
  }

  ep.max_outward_dist = ep.max_inward_dist = 1.0 ;
  ep.dstep = 0.25 ;
  sample_parameter_maps(mris, mri_T1, mri_PD, cv_wm_T1, cv_wm_PD,
                        cv_gm_T1, cv_gm_PD, cv_csf_T1, cv_csf_PD,
                        cv_inward_dists, cv_outward_dists, &ep, 0, &parms, 0) ;
  smooth_maps(mris, &ep, smooth_parms) ;

  write_maps(mris, &ep, -1, output_suffix) ;
  msec = then.milliseconds() ;
  fprintf(stderr,"positioning took %2.1f minutes\n", (float)msec/(60*1000.0f));
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s%s", sdir, sname, hemi,
		 white_matter_name, output_suffix, suffix) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  printf("writing final white matter position to %s...\n", fname) ;
  MRISwrite(mris, fname) ;
  MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s%s", sdir, sname, hemi,
		 pial_name, output_suffix, suffix) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("writing final pial surface position to %s...\n", fname) ;
  MRISwrite(mris, fname) ;


  /*  if (!(parms.flags & IPFLAG_NO_SELF_INT_TEST))*/
  {
    printf("measuring cortical thickness...\n") ;
    MRISmeasureCorticalThickness(mris, nbhd_size, max_thickness) ;
    printf(
      "writing cortical thickness estimate to 'thickness' file.\n") ;
    int req = snprintf(fname, STRLEN, "thickness%s%s", output_suffix, suffix) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRISwriteCurvature(mris, fname) ;

    /* at this point, the v->curv slots contain the cortical surface. Now
       move the white matter surface out by 1/2 the thickness as an estimate
       of layer IV.
    */
    if (graymid) {
      MRISsaveVertexPositions(mris, TMP_VERTICES) ;
      mrisFindMiddleOfGray(mris) ;
      int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s",
			 sdir, sname, hemi, GRAYMID_NAME,
			 suffix) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing layer IV surface to %s...\n", fname) ;
      MRISwrite(mris, fname) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    }
  }
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
  else if (!stricmp(option, "nbrs")) {
    nbrs = atoi(argv[2]) ;
    printf( "using neighborhood size = %d\n", nbrs) ;
    nargs = 1 ;
  } else if (!stricmp(option, "sample")) {
    if (stricmp(argv[2], "nearest") == 0) {
      sample_type = SAMPLE_NEAREST ;
      printf("using nearest nbr sampling\n") ;
    } else if (stricmp(argv[2], "trilinear") == 0) {
      sample_type = SAMPLE_TRILINEAR ;
      printf("using trilinear interpolation\n") ;
    } else
      ErrorExit(ERROR_BADPARM,  "unknown sampling type %s\n", argv[2]) ;

    nargs = 1 ;
  } else if (!stricmp(option, "fix_T1")) {
    fix_T1 = 1 ;
    printf( "fixing T1 values of gray and white matter\n");
  } else if (!stricmp(option, "min_max_scale")) {
    min_max_scale = atoi(argv[2]) ;
    nargs = 1 ;
    printf( "using min/max filter size %d (%2.3f mm)\n", min_max_scale,sqrt((float)min_max_scale));
  } else if (!stricmp(option, "orig")) {
    orig_flag = 1 ;
    printf( "using orig surface as initial surface placement\n") ;
  } else if (!stricmp(option, "start")) {
    start_white_name = argv[2] ;
    start_pial_name = argv[3] ;
    printf( "using %s and %s surfaces for initial placement\n",
            start_white_name, start_pial_name) ;
    orig_flag = 0 ;
    nargs = 2 ;
  } else if (!stricmp(option, "maxv")) {
    MAX_VNO = atoi(argv[2]) ;
    printf( "limiting  calculations to 1st %d vertices\n", MAX_VNO) ;
    nargs = 1 ;
  } else if (!stricmp(option, "minv")) {
    MIN_VNO = atoi(argv[2]) ;
    printf( "starting  calculations at vertex # %d\n", MIN_VNO) ;
    nargs = 1 ;
  } else if (!stricmp(option, "dstep")) {
    MAX_DSTEP = (double)atof(argv[2]) ;
    printf( "sampling volume every %2.2f mm\n", DSTEP) ;
    nargs = 1 ;
  } else if (!stricmp(option, "T1")) {
    printf("reading T1 parameter map from %s...\n", argv[2]) ;
    mri_T1 = MRIread(argv[2]) ;
    if (!mri_T1)
      ErrorExit(ERROR_NOFILE, "%s: could not read T1 volume from %s",
                Progname, argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "TOL")) {
    parms.tol = atof(argv[2]) ;
    printf("using integration tol %e\n", parms.tol) ;
    nargs = 1 ;
  } else if (!stricmp(option, "PD")) {
    printf("reading PD parameter map from %s...\n", argv[2]) ;
    mri_PD = MRIread(argv[2]) ;

    if (!mri_PD)
      ErrorExit(ERROR_NOFILE, "%s: could not read PD volume from %s",
                Progname, argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "brain")) {
    strcpy(brain_name, argv[2]) ;
    printf("using %s as brain volume...\n", brain_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "SDIR")) {
    gSdir = argv[2] ;
    printf("using %s as SUBJECTS_DIR...\n", gSdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "median")) {
    apply_median_filter = 1 ;
  } else if (!stricmp(option, "map_dir")) {
    map_dir = argv[2] ;
    nargs = 1 ;
    printf("reading parameter maps and residuals from %s...\n", map_dir) ;
  } else if (!stricmp(option, "graymid")) {
    graymid = 1 ;
    printf("generating graymid surface...\n") ;
  } else if (!strcmp(option, "rval")) {
    rh_label = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"using %d as fill val for right hemisphere.\n", rh_label);
  } else if (!strcmp(option, "nbhd_size")) {
    nbhd_size = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"using %d size nbhd for thickness calculation.\n",
            nbhd_size);
  } else if (!strcmp(option, "lval")) {
    lh_label = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"using %d as fill val for left hemisphere.\n", lh_label);
  } else if (!stricmp(option, "whiteonly")) {
    white_only = 1 ;
    printf( "only generating white matter surface\n") ;
  } else if (!stricmp(option, "pial")) {
    strcpy(pial_name, argv[2]) ;
    printf( "writing pial surface to file named %s\n", pial_name) ;
  } else if (!stricmp(option, "write_vals")) {
    write_vals = 1 ;
    printf( "writing gray and white surface targets to w files\n") ;
  } else if (!stricmp(option, "name")) {
    strcpy(parms.base_name, argv[2]) ;
    nargs = 1 ;
    printf("base name = %s\n", parms.base_name) ;
  } else if (!stricmp(option, "dt")) {
    parms.dt = atof(argv[2]) ;
    parms.base_dt = base_dt_scale*parms.dt ;
    parms.integration_type = INTEGRATE_MOMENTUM ;
    printf( "using dt = %2.1e\n", parms.dt) ;
    nargs = 1 ;
  } else if (!stricmp(option, "spring")) {
    parms.l_spring = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_spring = %2.3f\n", parms.l_spring) ;
  } else if (!stricmp(option, "repulse")) {
    parms.l_repulse = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_repulse = %2.3f\n", parms.l_repulse) ;
  } else if (!stricmp(option, "grad")) {
    parms.l_grad = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_grad = %2.3f\n", parms.l_grad) ;
  } else if (!stricmp(option, "external")) {
    parms.l_external = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_external = %2.3f\n", parms.l_external) ;
  } else if (!stricmp(option, "tspring")) {
    parms.l_tspring = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_tspring = %2.3f\n", parms.l_tspring) ;
  } else if (!stricmp(option, "nspring")) {
    parms.l_nspring = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_nspring = %2.3f\n", parms.l_nspring) ;
  } else if (!stricmp(option, "curv")) {
    parms.l_curv = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_curv = %2.3f\n", parms.l_curv) ;
  } else if (!stricmp(option, "smooth")) {
    smooth = atoi(argv[2]) ;
    nargs = 1 ;
    printf("smoothing for %d iterations\n", smooth) ;
  } else if (!stricmp(option, "smooth_parms")) {
    smooth_parms = atoi(argv[2]) ;
    nargs = 1 ;
    printf("smoothing parameter maps for %d iterations\n", smooth_parms) ;
  } else if (!stricmp(option, "output")) {
    output_suffix = argv[2] ;
    nargs = 1 ;
    printf("appending %s to output names...\n", output_suffix) ;
  } else if (!stricmp(option, "vavgs")) {
    vavgs = atoi(argv[2]) ;
    nargs = 1 ;
    printf("smoothing values for %d iterations\n", vavgs) ;
  } else if (!stricmp(option, "white")) {
    strcpy(white_matter_name, argv[2]) ;
    nargs = 1 ;
    printf("using %s as white matter name...\n", white_matter_name) ;
  } else if (!stricmp(option, "intensity")) {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_intensity = %2.3f\n", parms.l_intensity) ;
  } else if (!stricmp(option, "lm")) {
    parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
    printf("integrating with line minimization\n") ;
  } else if (!stricmp(option, "smoothwm")) {
    smoothwm = atoi(argv[2]) ;
    printf("writing smoothed (%d iterations) wm surface\n",
           smoothwm) ;
    nargs = 1 ;
  } else if (!stricmp(option, "sigma")) {
    sigma = atof(argv[2]) ;
    printf( "smoothing volume with Gaussian sigma = %2.1f\n",
            sigma) ;
    nargs = 1 ;
  } else if (!stricmp(option, "add")) {
    add = 1 ;
    printf("adding vertices to tessellation during deformation.\n");
    parms.flags |= IPFLAG_ADD_VERTICES ;
  } else switch (toupper(*option)) {
    case 'S':
      suffix = argv[2] ;
      printf("using %s as suffix\n", suffix) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'T':
      xform_fname = argv[2] ;
      nargs = 1;
      printf("applying ventricular xform %s\n", xform_fname);
      break ;
    case 'O':
      strcpy(orig_name, argv[2]) ;
      nargs = 1 ;
      printf("reading original vertex positions from %s\n", orig_name);
      break ;
    case 'Q':
      parms.flags |= IPFLAG_NO_SELF_INT_TEST ;
      printf(
        "doing quick (no self-intersection) surface positioning.\n") ;
      break ;
    case 'P':
      plot_stuff = 1 ;
      printf("plotting intensity profiles...\n") ;
      break ;
    case 'A':
      max_averages = atoi(argv[2]) ;
      printf("using max_averages = %d\n", max_averages) ;
      nargs = 1 ;
      if (isdigit(*argv[3])) {
        min_averages = atoi(argv[3]) ;
        printf("using min_averages = %d\n", min_averages) ;
        nargs++ ;
      }
      break ;
    case 'M':
      parms.integration_type = INTEGRATE_MOMENTUM ;
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      printf("momentum = %2.2f\n", parms.momentum) ;
      break ;
    case 'R':
      parms.l_surf_repulse = atof(argv[2]) ;
      printf("l_surf_repulse = %2.3f\n", parms.l_surf_repulse) ;
      nargs = 1 ;
      break ;
    case 'B':
      base_dt_scale = atof(argv[2]) ;
      parms.base_dt = base_dt_scale*parms.dt ;
      nargs = 1;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      printf("printing diagnostic info for vertex %d\n", Gdiag_no) ;
      nargs = 1 ;
      break ;
    case 'C':
      create = !create ;
      printf("%screating area and curvature files for wm surface...\n",
             create ? "" : "not ") ;
      break ;
    case 'W':
      sscanf(argv[2], "%d", &parms.write_iterations) ;
      nargs = 1 ;
      printf("write iterations = %d\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'N':
      sscanf(argv[2], "%d", &parms.niterations) ;
      nargs = 1 ;
      printf("niterations = %d\n", parms.niterations) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
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
  printf("usage: %s [options] <subject name> <hemisphere> <xform> <flash 1> <flash 2> .. <residuals>\n",
         Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  printf(
    "\nThis program positions the tessellation of the cortical surface\n"
    "at the white matter surface, then the gray matter surface\n"
    "and generate surface files for these surfaces as well as a\n"
    "'curvature' file for the cortical thickness, and a surface file\n"
    "which approximates layer IV of the cortical sheet.\n");
  printf("\nvalid options are:\n\n") ;
  printf(
    "-q    omit self-intersection and only generate "
    "gray/white surface.\n") ;
  printf(
    "-c    create curvature and area files from white matter surface\n"
  );
  printf(
    "-a <avgs>   average curvature values <avgs> times (default=10)\n");
  printf(
    "-whiteonly  only generate white matter surface\n") ;
  exit(1) ;
}

static void
print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}

#define MAX_DEFORM_DIST  3
static double
ms_errfunc_gradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) {
  double      nx, ny, nz, lambda, white_delta, pial_delta, len, sse ;
  int         vno, white_vno, pial_vno, ndone ;
  VERTEX      *v, *v_pial, *v_white ;
  EXTRA_PARMS *ep ;


  ep = (EXTRA_PARMS *)parms->user_parms ;
  store_current_positions(mris, ep) ;
  update_parameters(mris, ep) ;
  lambda = parms->l_external ;
  ndone = 0 ;
#if 1
  compute_maximal_distances(mris, ep->current_sigma, ep->mri_flash, ep->nvolumes,
                            ep->cv_inward_dists, ep->cv_outward_dists,
                            ep->nearest_pial_vertices, ep->nearest_white_vertices,
                            ep->dstep, ep->max_inward_dist, ep->max_outward_dist,
                            ep->mri_T1, ep->mri_PD, ep) ;
#endif
  sse = 0.0 ;
  for (vno = MIN_VNO ; vno < mris->nvertices ; vno++) {
#if 0
    if ((vno+1) % (mris->nvertices/10) == 0)
      printf("%d ", vno) ;
#endif
    if (vno > MAX_VNO)
      break ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;

    compute_optimal_vertex_positions(mris, vno, ep, &white_delta,&pial_delta,1);
    sse += (white_delta*white_delta) + (pial_delta*pial_delta) ;
    if (Gdiag_no == vno)
      printf("v %d: sse %2.2f\n", vno, (white_delta*white_delta) + (pial_delta*pial_delta)) ;

    if (fabs(white_delta) < 0.3 && fabs(pial_delta) < 0.3)
      ndone++ ;
    white_vno = vno ;
#if 0
    pial_vno = ep->nearest_pial_vertices[vno] ;
#else
    pial_vno = vno ;
#endif


    if (fabs(white_delta) > 1)
      white_delta /= fabs(white_delta) ;
    if (fabs(pial_delta) > 1)
      pial_delta /= fabs(pial_delta) ;
    v_white = &mris->vertices[white_vno] ;
    v_pial = &mris->vertices[pial_vno] ;
    nx = v_pial->pialx - v_white->origx ;
    ny = v_pial->pialy - v_white->origy ;
    nz = v_pial->pialz - v_white->origz ;
    len = sqrt(nx*nx + ny*ny + nz*nz) ;
    if (TOO_SMALL(len)) {
      nx = v_white->nx ;
      ny = v_white->ny ;
      nz = v_white->nz ;
      len = sqrt(nx*nx + ny*ny + nz*nz) ;
      if (FZERO(len))
        len = 0.01 ;
    }
    nx /= len ;
    ny /= len ;
    nz /= len ;
    if (Gdiag_no == white_vno || Gdiag_no == pial_vno) {
      printf("v %d,%d: delta = %2.2f, %2.2f along n = (%2.2f, %2.2f, %2.2f)\n",
             white_vno, pial_vno, lambda*white_delta, lambda*pial_delta, nx, ny, nz) ;
      printf("      current location   = (%2.2f, %2.2f, %2.2f) --> (%2.2f, %2.2f, %2.2f)\n",
             v_white->origx, v_white->origy, v_white->origz,
             v_pial->pialx, v_pial->pialy, v_pial->pialz);
      printf("      predicted location = (%2.2f, %2.2f, %2.2f) --> (%2.2f, %2.2f, %2.2f)\n",
             v_white->origx+lambda * white_delta * nx,
             v_white->origy+lambda * white_delta * ny,
             v_white->origz+lambda * white_delta * nz,
             v_pial->pialx+lambda * pial_delta * nx,
             v_pial->pialy+lambda * pial_delta * ny,
             v_pial->pialz+lambda * pial_delta * nz) ;

    }
    v_white->dx += lambda * white_delta * nx ;
    v_white->dy += lambda * white_delta * ny ;
    v_white->dz += lambda * white_delta * nz ;

    mris->dx2[pial_vno] += lambda * pial_delta * nx ;
    mris->dy2[pial_vno] += lambda * pial_delta * ny ;
    mris->dz2[pial_vno] += lambda * pial_delta * nz ;
    if (white_vno == Gdiag_no)
      printf("v %d white intensity   :  (%2.3f, %2.3f, %2.3f)\n",
             vno, lambda * white_delta * nx, lambda * white_delta * ny,
             lambda * white_delta * nz) ;
    if (pial_vno == Gdiag_no)
      printf("v %d pial intensity    :  (%2.3f, %2.3f, %2.3f)\n",
             vno, lambda * pial_delta * nx, lambda * pial_delta * ny,
             lambda * pial_delta * nz) ;
  }

  if (Gdiag_no >= 0 && DIAG_VERBOSE_ON) {
    int   n, n_vno ;
    float  dx, dy, dz, dot ;

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[Gdiag_no];
    dx = mris->dx2[Gdiag_no] ;
    dy = mris->dy2[Gdiag_no] ;
    dz = mris->dz2[Gdiag_no] ;
    for (n = 0 ; n < vt->vnum ; n++) {
      n_vno = vt->v[n] ;
      dot = mris->dx2[n_vno]*dx + mris->dy2[n_vno]*dy + mris->dz2[n_vno]*dz ;
      if (dot < 0) {
        printf("vertex %d: dx = (%2.1f, %2.1f, %2.1f), dot = %2.2f\n",
               n_vno, mris->dx2[n_vno], mris->dy2[n_vno], mris->dz2[n_vno], dot) ;
      }
    }
  }

  printf("%d vertices asymptoted (%2.3f%%)\n", ndone,
         100.0f*(float)ndone/mris->nvertices) ;
#if 0
  printf("\n") ;
#endif
  return(sse) ;
}
static int
ms_errfunc_timestep(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) {
  EXTRA_PARMS *ep ;

  ep = (EXTRA_PARMS *)parms->user_parms ;
  update_distances(mris, ep) ;  /* will mark vertices that need parms recomputed */
  store_current_positions(mris, ep) ;
  if (ep->nvertices < mris->nvertices) {
    printf("resampling parameters for %d new vertices...\n",
           mris->nvertices-ep->nvertices) ;
    sample_parameter_maps(mris, mri_T1, mri_PD, ep->cv_wm_T1, ep->cv_wm_PD,
                          ep->cv_gm_T1, ep->cv_gm_PD, ep->cv_csf_T1, ep->cv_csf_PD,
                          ep->cv_inward_dists, ep->cv_outward_dists, ep, fix_T1, parms,
                          ep->nvertices) ;
    ep->nvertices = mris->nvertices ;

  }
  return(NO_ERROR) ;
}

static double last_sse[300000];
static double
ms_errfunc_sse(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) {
  double      sse, total_sse, white_delta, pial_delta, wm_total, pial_total ;
  int         vno ;
  VERTEX      *v ;
  EXTRA_PARMS *ep ;

  ep = (EXTRA_PARMS *)parms->user_parms ;
#if 0
  compute_maximal_distances(mris, ep->current_sigma, ep->mri_flash, ep->nvolumes,
                            ep->cv_inward_dists, ep->cv_outward_dists,
                            ep->nearest_pial_vertices, ep->nearest_white_vertices,ep->dstep,
                            ep->max_inward_dist, ep->max_outward_dist,
                            ep->mri_T1, ep->mri_PD, ep);
#endif
  wm_total = pial_total = total_sse = 0.0 ;
  for (vno = MIN_VNO ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (vno > MAX_VNO)
      break ;
    if (v->ripflag)
      continue ;
#if 0
    sse = vertex_error(mris, vno, ep, NULL) ;
#else
    compute_optimal_vertex_positions(mris, vno, ep, &white_delta,&pial_delta,
                                     plot_stuff);
    sse = (white_delta*white_delta) + (pial_delta*pial_delta) ;
    wm_total += fabs(white_delta) ;
    pial_total += fabs(pial_delta) ;
#endif
    if (!std::isfinite(sse))
      DiagBreak() ;
    total_sse += sse ;
    if (Gdiag_no == vno)
      printf("flash intensity term v %d: sse = %2.2f (%2.2f, %2.2f)\n",
             vno, sse, white_delta, pial_delta) ;
    if (sse > last_sse[vno] && !FZERO(last_sse[vno]))
      DiagBreak() ;
    last_sse[vno] = sse ;
  }

  if (mris->nvertices < MAX_VNO) {
    wm_total /= (mris->nvertices - MIN_VNO) ;
    pial_total /= (mris->nvertices - MIN_VNO) ;
  } else {
    wm_total /= (MAX_VNO - MIN_VNO) ;
    pial_total /= (MAX_VNO - MIN_VNO) ;
  }
  printf("mean distances = (%2.2f, %2.2f)\n", wm_total, pial_total) ;
  return(total_sse) ;
}

static double
ms_errfunc_rms(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) {
  double      rms, rms_total, white_delta, pial_delta ;
  int         vno ;
  VERTEX      *v ;
  EXTRA_PARMS *ep ;

  ep = (EXTRA_PARMS *)parms->user_parms ;
  rms_total = 0.0 ;
  for (vno = MIN_VNO ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (vno > MAX_VNO)
      break ;
    if (v->ripflag)
      continue ;
#if 0
    vertex_error(mris, vno, ep, &rms) ;
#else
    compute_optimal_vertex_positions(mris, vno, ep, &white_delta,&pial_delta,
                                     plot_stuff);
    rms = sqrt((white_delta*white_delta) + (pial_delta*pial_delta)) ;
#endif
    rms_total += rms ;
    if (!std::isfinite(rms))
      DiagBreak() ;
    if (Gdiag_no == vno)
      printf("v %d: rms = %2.3f\n", vno, rms) ;
  }

  return(rms_total/(double)vno) ;
}

#if 0
static double
vertex_error(MRI_SURFACE *mris, int vno, EXTRA_PARMS *ep, double *prms) {
  VERTEX  *v_white, *v_pial ;
  double  dist, dx, dy, dz, cortical_dist, sigma,
  T1_wm, T1_gm, T1_csf, PD_wm, PD_gm, PD_csf, sse, inward_dist, outward_dist ;
  double  xp, yp, zp, xw, yw, zw, x, y, z,
  image_vals[MAX_FLASH_VOLUMES][MAX_SAMPLES] ;
  MRI     *mri ;
  int     i, j, max_j ;

  v_white = &mris->vertices[vno] ;
#if 0
  v_pial = &mris->vertices[ep->nearest_pial_vertices[vno]] ;
#else
  v_pial = v_white ;
#endif

  sigma = ep->current_sigma ;
  sse = 0.0 ;
  T1_wm = ep->cv_wm_T1[vno] ;
  PD_wm = ep->cv_wm_PD[vno] ;
  T1_gm = ep->cv_gm_T1[vno] ;
  PD_gm = ep->cv_gm_PD[vno] ;
  T1_csf = ep->cv_csf_T1[vno] ;
  PD_csf = ep->cv_csf_PD[vno] ;
  inward_dist = ep->cv_inward_dists[vno] ;
  outward_dist = ep->cv_outward_dists[vno] ;
  mri = ep->mri_flash[0] ;
  // MRIworldToVoxel(mri, v_white->origx, v_white->origy, v_white->origz,
  //                 &xw, &yw, &zw) ;
  // MRIworldToVoxel(mri, v_pial->pialx, v_pial->pialy, v_pial->pialz,
  //                &xp, &yp, &zp) ;
  MRIsurfaceRASToVoxel(mri, v_white->origx, v_white->origy, v_white->origz,
                       &xw, &yw, &zw) ;
  MRIsurfaceRASToVoxel(mri, v_pial->pialx, v_pial->pialy, v_pial->pialz,
                       &xp, &yp, &zp) ;

  dx = xp-xw ;
  dy = yp-yw ;
  dz = zp-zw ;
  cortical_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
  if (TOO_SMALL(cortical_dist)) {
    // MRIworldToVoxel(mri,
    //                 v_pial->pialx+v_pial->nx,
    //                 v_pial->pialy+v_pial->nx,
    //                 v_pial->pialz+v_pial->nz,
    //                 &x, &y, &z) ;
    MRIsurfaceRASToVoxel(mri,
                         v_pial->pialx+v_pial->nx,
                         v_pial->pialy+v_pial->nx,
                         v_pial->pialz+v_pial->nz,
                         &x, &y, &z) ;

    dx = x-xw ;
    dy = y-yw ;
    dz = z-zw ;
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (!FZERO(dist)) {
      dx /= dist ;
      dy /= dist ;
      dz /= dist ;
    }
  } else {
    dx /= cortical_dist ;
    dy /= cortical_dist ;
    dz /= cortical_dist ;
  }

  max_j = (int)((inward_dist + cortical_dist + outward_dist) / ep->dstep) ;


  for (i = 0 ; i < ep->nvolumes ; i++) {
    mri = ep->mri_flash[i] ;
    for (j = 0, dist = -inward_dist ; dist <= cortical_dist+outward_dist ; dist += ep->dstep, j++) {
      x = xw + dist * dx ;
      y = yw + dist * dy ;
      z = zw + dist * dz ;
#if 0
      MRIsampleVolumeDirectionScale(mri, x, y, z, dx, dy, dz, &image_vals[i][j], sigma) ;
#else
      MRIsampleVolumeType(mri, x, y, z, &image_vals[i][j], sample_type) ;
#endif
      if (!finite(image_vals[i][j]))
        DiagBreak() ;
    }

    if ((plot_stuff == 1 && Gdiag_no == vno) ||
        plot_stuff > 1) {
      FILE *fp ;
      char fname[STRLEN] ;

      sprintf(fname, "sse_vol%d.plt", i) ;
      fp = fopen(fname, "w") ;
      for (j = 0, dist = -inward_dist ; dist <= cortical_dist+outward_dist ; dist += ep->dstep, j++)
        fprintf(fp, "%d %f %f\n", j, dist, image_vals[i][j]) ;
      fclose(fp) ;
    }
  }
  sse = compute_vertex_sse(ep, image_vals, max_j, inward_dist, cortical_dist,
                           T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf, plot_stuff && vno == Gdiag_no,
                           vno) ;
  if (prms)
    *prms = sqrt(sse) ;

  return(sse) ;
}
#endif


#if 0
#define MAX_TABLES     10
static float  TRs[MAX_TABLES] ;
static double *exp_minus_TR_div_T1[MAX_TABLES] ;

static int
init_lookup_table(MRI *mri) {
  int    i ;
  double TR = mri->tr ;
  double T1 ;

  for (i = 0 ; i < MAX_TABLES ; i++)
    if ((NULL != exp_minus_TR_div_T1[i]) && FEQUAL(TRs[i], TR))
      return(NO_ERROR) ;   /* table alread created */

  for (i = 0 ; i < MAX_TABLES ; i++)
    if (NULL == exp_minus_TR_div_T1[i])
      break ;
  if (i >= MAX_TABLES)
    return(NO_ERROR) ;    /* can't fit it */

  TRs[i] = TR ;
  exp_minus_TR_div_T1[i] = (double *)calloc(MAX_T1, sizeof(double)) ;
  if (!exp_minus_TR_div_T1[i])
    ErrorExit(ERROR_NOMEMORY, "init_lookup_table(%2.1f) - could not alloc table",
              TR) ;

  for (T1 = 0.0 ; T1 < MAX_T1 ; T1 += 1.0)
    exp_minus_TR_div_T1[i][(int)T1] = exp(-TR/T1) ;

  printf("creating lookup table for TR=%2.1f ms\n", (float)TR) ;

  return(NO_ERROR) ;
}

double *
find_lookup_table(double TR) {
  int i ;

  for (i = 0 ; i < MAX_TABLES ; i++)
    if ((NULL != exp_minus_TR_div_T1[i]) && FEQUAL(TRs[i], TR))
      return(exp_minus_TR_div_T1[i]) ;

  return(NULL) ;
}
#endif

static double
FLASHforwardModel(double flip_angle, double TR, double PD, double T1) {
  static double  CFA = 1, SFA = 0 ;
  static double last_flip_angle = -1 ;
  double  E1, FLASH ;

  if (!DZERO(flip_angle-last_flip_angle)) {
    CFA = cos(flip_angle) ;
    SFA = sin(flip_angle) ;
    last_flip_angle = flip_angle ;
  }

  E1 = exp(-TR/T1) ;

  FLASH = PD * SFA ;
  if (!DZERO(T1))
    FLASH *= (1-E1)/(1-CFA*E1);
  return(FLASH) ;
}

#if 1

static double
FLASHforwardModelLookup(double flip_angle, double TR, double PD, double T1) {
  double FLASH ;

  FLASH = lookup_flash_value(TR, flip_angle, PD, T1) ;
  return(FLASH) ;
}
#else
static double
FLASHforwardModelLookup(double flip_angle, double TR, double PD, double T1) {
  double FLASH, E1, *table ;
  static double  CFA = 1, SFA = 0 ;
  static double last_flip_angle = -1 ;

  if (!DZERO(flip_angle-last_flip_angle)) {
    CFA = cos(flip_angle) ;
    SFA = sin(flip_angle) ;
    last_flip_angle = flip_angle ;
  }

  table = find_lookup_table(TR) ;
  if (NULL == table)
    E1 = exp(-TR/T1) ;
  else
    E1 = table[nint(T1)] ;

  FLASH = PD * SFA ;
  if (!DZERO(T1))
    FLASH *= (1-E1)/(1-CFA*E1);
  return(FLASH) ;
}
#endif

#if 0
static double
dM_dT1(double flip_angle, double TR, double PD, double T1) {
  double  dT1, E1 ;
  static double  CFA, SFA2_3, CFA2 ;


  CFA = cos(flip_angle) ;
  SFA2_3 = sin(flip_angle/2) ;
  SFA2_3 = SFA2_3*SFA2_3*SFA2_3 ;
  CFA2 = cos(flip_angle/2) ;
  E1 = exp(TR/T1) ;

  dT1 = -4*E1*PD*TR*CFA2*SFA2_3/ (T1*T1*pow(-E1 + CFA,2)) ;

  return(dT1) ;
}

static double
dM_dPD(double flip_angle, double TR, double PD, double T1) {
  double  dPD, E1 ;
  static double  CFA, SFA ;

  CFA = cos(flip_angle) ;
  SFA = sin(flip_angle) ;
  E1 = exp(TR/T1) ;
  dPD = (-1 + E1)*SFA/(E1 - CFA) ;

  return(dPD) ;
}
#endif

static int
compute_maximal_distances(MRI_SURFACE *mris, float sigma, MRI **mri_flash, int nvolumes,
                          float *cv_inward_dists, float *cv_outward_dists,
                          int *nearest_pial_vertices, int *nearest_white_vertices,
                          double dstep, double max_inward_dist, double max_outward_dist_total,
                          MRI *mri_T1, MRI *mri_PD, EXTRA_PARMS *ep) {
  int    vno, i, found_csf, j, max_j = 0, found_wm, past_wm, wm_dist, csf_dist, min_j, low_pd_csf;
  /*  int    nreversed, max_reversed,done;*/
  double T1, PD, max_outward_dist ;
  VERTEX *v_white, *v_pial ;
  MRI    *mri ;
  float  nx, ny, nz, min_inward_dist, min_outward_dist, dist,
  min_parm_dist, parm_dist ;
  double xw, yw, zw, xp, yp, zp, xo, yo, zo, cortical_dist,
  image_vals[MAX_FLASH_VOLUMES][MAX_SAMPLES] ;
  double T1_vals[MAX_SAMPLES], PD_vals[MAX_SAMPLES] ;
  /*  double dIdN_start[MAX_FLASH_VOLUMES], dIdN ;*/



  nx = ny = nz = 0 ;   /* remove compiler warning */
  for (vno = MIN_VNO ; vno < mris->nvertices ; vno++) {
    if (vno >= MAX_VNO)
      break ;
    v_white = &mris->vertices[vno] ;
#if 0
    v_pial = &mris->vertices[nearest_pial_vertices[vno]] ;
#else
    v_pial = v_white ;
#endif
    if (v_white->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    min_inward_dist = max_inward_dist ;
    max_j = (max_inward_dist-dstep) / dstep ;
    for (i = 0 ; i < nvolumes ; i++) {
      mri = mri_flash[i] ;
      // converting surface vertex to volume voxel
      // MRIworldToVoxel(mri, v_white->origx, v_white->origy, v_white->origz, &xw, &yw, &zw) ;
      // MRIworldToVoxel(mri, v_pial->pialx, v_pial->pialy, v_pial->pialz, &xp, &yp, &zp) ;
      MRIsurfaceRASToVoxel(mri, v_white->origx, v_white->origy, v_white->origz, &xw, &yw, &zw) ;
      MRIsurfaceRASToVoxel(mri, v_pial->pialx, v_pial->pialy, v_pial->pialz, &xp, &yp, &zp) ;
      nx = xp - xw ;
      ny = yp - yw ;
      nz = zp - zw ;
      dist = sqrt(nx*nx + ny*ny + nz*nz) ;
      if (TOO_SMALL(dist)) {
        // MRIworldToVoxel(mri,
        //                 v_pial->pialx+v_white->nx,
        //                 v_pial->pialy+v_white->ny,
        //                 v_pial->pialz+v_white->nz, &xp, &yp, &zp) ;
        MRIsurfaceRASToVoxel(mri,
                             v_pial->pialx+v_white->nx,
                             v_pial->pialy+v_white->ny,
                             v_pial->pialz+v_white->nz, &xp, &yp, &zp) ;
        nx = xp - xw ;
        ny = yp - yw ;
        nz = zp - zw ;
        dist = sqrt(nx*nx + ny*ny + nz*nz) ;
      }
      if (dist < 0.0001) dist = 0.001 ;
      nx /= dist ;
      ny /= dist ;
      nz /= dist ;

      /* first search inwards */
      for (j = 0, dist = dstep ; j <= max_j ; dist += dstep, j++) {
        xo = xw-dist*nx ;
        yo = yw-dist*ny ;
        zo = zw-dist*nz ;
        MRIsampleVolumeType(mri, xo, yo, zo, &image_vals[i][j], sample_type) ;
        MRIsampleVolumeType(mri_T1, xo, yo, zo, &T1_vals[j], sample_type) ;
        MRIsampleVolumeType(mri_PD, xo, yo, zo, &PD_vals[j], sample_type) ;
      }
      dist -= dstep ;
      dist = MAX(dist, dstep) ;
      if (min_inward_dist > dist)
        min_inward_dist = dist ;
    }
    /* search to see if T1/PD pairs are reasonable for WM */
    wm_dist = 0 ;
    for (j = found_wm = 0 ; j <= max_j ; j++) {
      T1 = T1_vals[j] ;
      PD = PD_vals[j] ;
      if (found_wm) {
        if (found_wm > 1 && IS_GM(T1,PD,vno,ep))
          break ;
        else if (!IS_GM(T1,PD,vno,ep))
          found_wm = 2 ;
        if (!IS_WM(T1,PD,vno,ep))
          break ;
#if 0
        if (++wm_dist >= 5)
          break ;
#endif
      } else if (IS_WM(T1,PD,vno,ep)) {
        found_wm = 1 ;
        if (!IS_GM(T1,PD,vno,ep))
          found_wm=2 ;
      }
    }
    if (found_wm == 0) {
      wm_dist = 0 ;
      min_parm_dist = MAX_T1*MAX_T1 ;
      min_j = 0 ;
      for (j = 0 ; j <= max_j ; j++) {
        T1 = T1_vals[j] ;
        PD = PD_vals[j] ;
        parm_dist = WM_PARM_DIST(T1,PD) ;
        if (parm_dist < min_parm_dist) {
          min_parm_dist = parm_dist ;
          min_j = j ;
        }
      }
      j = min_j+1 ;
    }

    min_inward_dist = ((float)j)*dstep ;
    cv_inward_dists[vno] = min_inward_dist ;
    if (vno == Gdiag_no)
      printf("v %d: inward_dist = %2.2f, X=(%2.1f, %2.1f, %2.1f)\n",
             vno, min_inward_dist,v_white->origx,v_white->origy,v_white->origz);
  }

  for (vno = MIN_VNO ; vno < mris->nvertices ; vno++) {
    if (vno >= MAX_VNO)
      break ;
    v_pial = &mris->vertices[vno] ;
    if (v_pial->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
#if 0
    v_white = &mris->vertices[nearest_white_vertices[vno]] ;
#else
    v_white = v_pial ;
#endif

    found_csf = 0 ;
    // MRIworldToVoxel(mri_flash[0],v_white->origx, v_white->origy, v_white->origz, &xw,&yw,&zw);
    // MRIworldToVoxel(mri_flash[0],v_pial->pialx, v_pial->pialy, v_pial->pialz, &xp, &yp, &zp) ;
    MRIsurfaceRASToVoxel(mri_flash[0],v_white->origx, v_white->origy, v_white->origz, &xw,&yw,&zw);
    MRIsurfaceRASToVoxel(mri_flash[0],v_pial->pialx, v_pial->pialy, v_pial->pialz, &xp, &yp, &zp) ;
    nx = xp - xw ;
    ny = yp - yw ;
    nz = zp - zw ;
    cortical_dist = dist = sqrt(nx*nx + ny*ny + nz*nz) ;
    max_outward_dist = max_outward_dist_total /* - cortical_dist*/ ;

    min_outward_dist = max_outward_dist ;
    max_j = (max_outward_dist-dstep) / dstep ;
    for (i = 0 ; i < nvolumes ; i++) {
      mri = mri_flash[i] ;
      // MRIworldToVoxel(mri, v_white->origx, v_white->origy, v_white->origz, &xw, &yw, &zw) ;
      // MRIworldToVoxel(mri, v_pial->pialx, v_pial->pialy, v_pial->pialz, &xp, &yp, &zp) ;
      MRIsurfaceRASToVoxel(mri, v_white->origx, v_white->origy, v_white->origz, &xw, &yw, &zw) ;
      MRIsurfaceRASToVoxel(mri, v_pial->pialx, v_pial->pialy, v_pial->pialz, &xp, &yp, &zp) ;
      nx = xp - xw ;
      ny = yp - yw ;
      nz = zp - zw ;
      cortical_dist = dist = sqrt(nx*nx + ny*ny + nz*nz) ;
      if (TOO_SMALL(dist)) {
        double xpn, ypn, zpn ;

        // MRIworldToVoxel(mri,
        //                 v_pial->pialx+v_white->nx,
        //                 v_pial->pialy+v_white->ny,
        //                 v_pial->pialz+v_white->nz, &xpn, &ypn, &zpn) ;
        MRIsurfaceRASToVoxel(mri,
                             v_pial->pialx+v_white->nx,
                             v_pial->pialy+v_white->ny,
                             v_pial->pialz+v_white->nz, &xpn, &ypn, &zpn) ;
        nx = xpn - xp ;
        ny = ypn - yp ;
        nz = zpn - zp ;
        dist = sqrt(nx*nx + ny*ny + nz*nz) ;
      }

      if (dist < 0.0001) dist = 0.001 ;
      nx /= dist ;
      ny /= dist ;
      nz /= dist ;

      /* search outwards - match sampling of compute_optimal_parameters */
      dist = (int)(cortical_dist/dstep)*dstep + dstep ;
      for (j = 0 ; j <= max_j ; dist += dstep, j++) {
        xo = xw+dist*nx ;
        yo = yw+dist*ny ;
        zo = zw+dist*nz ;
        MRIsampleVolumeType(mri, xo, yo, zo, &image_vals[i][j], sample_type) ;
        MRIsampleVolumeType(mri_T1, xo, yo, zo, &T1_vals[j], sample_type) ;
        MRIsampleVolumeType(mri_PD, xo, yo, zo, &PD_vals[j], sample_type) ;
      }
      dist -= dstep ;
      dist = MAX(dstep, dist) ;
      if (min_outward_dist > dist)
        min_outward_dist = dist ;
    }

    /* search to see if T1/PD pairs are reasonable for CSF */
    csf_dist = 0 ;
    past_wm = 0 ;
    low_pd_csf = 0 ;
    for (j = found_csf = 0 ; j <= max_j ; j++) {
      T1 = T1_vals[j] ;
      PD = PD_vals[j] ;
      if (past_wm && IS_WM(T1,PD,vno,ep))  /* don't go through wm twice */
        break ;
      if (!IS_WM(T1,PD,vno,ep))
        past_wm = 1 ;
      if (found_csf) {
        if (!IS_CSF(T1,PD))
          break ;
        /* don't include multiple types of non-brain stuff (the model doesn't support it) */
        if ((PD < MIN_NONBRAIN_PD && low_pd_csf == 0) ||
            (PD > MIN_NONBRAIN_PD && low_pd_csf != 0)) {
          if (vno == Gdiag_no)
            DiagBreak() ;
          break ;
        }
        if (low_pd_csf != 0 && PD < MIN_RELIABLE_PD) {
          if (vno == Gdiag_no)
            DiagBreak() ;
          break ;  /* don't include very low pd stuff, as it's T1 will be arbitrary */
        }

        if (++csf_dist >= 5)
          break ;
      } else if (!IS_BRAIN(T1,PD,vno,ep))  /* not brain */
      {
        if (PD < MIN_NONBRAIN_PD)
          low_pd_csf = 1 ;
        else
          low_pd_csf = 0 ;
        found_csf = 1 ;
      }
    }
    if (!found_csf)   /* look for partial-volumed csf */
    {
      if (vno == Gdiag_no)
        printf("couldn't find unambiguous csf - searching for partial volumed csf...") ;
      csf_dist = 0 ;
      min_parm_dist = MAX_T1*MAX_T1 ;
      min_j = 0 ;
      past_wm = 0 ;
      for (j = 0 ; j <= max_j ; j++) {
        T1 = T1_vals[j] ;
        PD = PD_vals[j] ;
        if (!IS_WM(T1,PD,vno,ep))
          past_wm = 1 ;

        if (past_wm && (IS_WM(T1,PD,vno,ep) || (!IS_CSF(T1,PD) && !IS_GM(T1,PD,vno,ep))))
          break ;   /* don't cross through wm to get to csf */
        parm_dist = CSF_PARM_DIST(T1,PD) ;
        if (parm_dist < min_parm_dist && IS_PV_CSF(T1,PD)) {
          found_csf = 1 ;
          min_parm_dist = parm_dist ;
          min_j = j ;
        }
      }
      min_outward_dist = ((float)min_j+1)*dstep ;
      if (vno == Gdiag_no) {
        if (found_csf)
          printf("csf found at distance %2.1f\n", min_outward_dist) ;
        else
          printf("couldn't find pv csf.\n") ;
      }
    }


#if 1
    if (!found_csf) {
      for (j = 0 ; j < max_j ; j++)    /* assume T1 and PD are monotonically increasing */
      {
        if ((T1_vals[j+1] < .99*T1_vals[j])
            /*            ||(PD_vals[j+1] < .99*PD_vals[j])*/
           )
          break ;
      }
    }
#else
    if (found_csf) {
      min_outward_dist = ((float)j)*dstep ;
      max_reversed = nvolumes ;
    } else
      max_reversed = 1 ;

    /* compute original directional derivatives */
    for (i = 0; i < nvolumes ; i++)
      MRIsampleVolumeDerivativeScale(mri_flash[i], xp-nx, yp-ny, zp-nz,
                                     nx, ny, nz, &dIdN_start[i], sigma) ;
    /* find distance outward over which directional derivative doesn't
    change sign.
    */
    for (done = 0, j = 0 ; !done && j <= max_j ; j++) {
      nreversed = 0 ;
      dist = j*dstep ;
      xo = xp+dist*nx ;
      yo = yp+dist*ny ;
      zo = zp+dist*nz ;
      for (i = 0 ; i < nvolumes ; i++) {
        MRIsampleVolumeDerivativeScale(mri_flash[i], xo, yo, zo,
                                       nx, ny, nz, &dIdN, sigma) ;
        if (dIdN * dIdN_start[i] < 0) {
          if (++nreversed >= max_reversed) {
            done = 1 ;
            break ;
          }
        }
      }
      if (done)
        break ;
#if 0
      /* ????????????????? */
      T1 = T1_vals[j] ;
      PD = PD_vals[j] ;
      if (!IS_BRAIN(T1,PD,vno,ep))
        break ;
#endif
    }
#endif
    min_outward_dist = ((float)j)*dstep ;

#if 0
    {
      csf_dist = 0 ;
      min_parm_dist = MAX_T1*MAX_T1 ;
      min_j = 0 ;
      past_wm = 0 ;
      for (j = 0 ; j <= max_j ; j++) {
        T1 = T1_vals[j] ;
        PD = PD_vals[j] ;
        if (!IS_WM(T1,PD,vno,ep))
          past_wm = 1 ;

        if (past_wm && (IS_WM(T1,PD,vno,ep) || (!IS_CSF(T1,PD) && !IS_GM(T1,PD,vno,ep))))
          break ;   /* don't cross through wm to get to csf */
        parm_dist = CSF_PARM_DIST(T1,PD) ;
        if (parm_dist < min_parm_dist) {
          min_parm_dist = parm_dist ;
          min_j = j ;
        }
      }
      min_outward_dist = ((float)min_j+1)*dstep ;
    }
#endif

    cv_outward_dists[vno] = min_outward_dist ; /* j=0 is dist=dstep */
    if (vno == Gdiag_no && min_outward_dist >= 2)
      DiagBreak() ;
    if (vno == Gdiag_no)
      printf("v %d: outward_dist = %2.2f, X=(%2.1f, %2.1f, %2.1f), N=(%2.1f, %2.1f, %2.1f)\n",
             vno, min_outward_dist, v_pial->pialx, v_pial->pialy, v_pial->pialz, nx, ny, nz);
  }

  return(NO_ERROR) ;
}
#define MIN_CSF   2000
#define CSF_T1    3000
#define WM_T1     750
#define GM_T1     1050

#define WM_T1_STD 100
#define GM_T1_STD 100

#define CSF_PD    1500
#define WM_PD     700
#define GM_PD     1050

#define BIN_SIZE 20
#define MAX_PARM (BIN_SIZE*MAX_BINS)

static int
sample_parameter_maps(MRI_SURFACE *mris, MRI *mri_T1, MRI *mri_PD,
                      float *cv_wm_T1, float *cv_wm_PD,
                      float *cv_gm_T1, float *cv_gm_PD,
                      float *cv_csf_T1, float *cv_csf_PD,
                      float *cv_inward_dists, float *cv_outward_dists,
                      EXTRA_PARMS *ep, int fix_T1, INTEGRATION_PARMS *parms,
                      int start_vno) {
  VERTEX    *v ;
  int       vno, nholes =0/*, bno*/ ;
  double    /*inward_dist, outward_dist,*/ dstep,
  best_white_delta = 0, best_pial_delta = 0, sse ;

#if 1

  if (fix_T1) {
    int vno ;

    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      cv_gm_T1[vno] = GM_T1 ;
      cv_wm_T1[vno] = WM_T1 ;
    }
  }

#if 0
  /* compute csf parameters */
  sample_parameter_map(mris, mri_T1, mri_res, cv_csf_T1, cv_outward_dists, 1,
                       PIAL_VERTICES, "csf T1", NULL, ep->dstep, ep->max_dist);
  sample_parameter_map(mris, mri_PD, mri_res, cv_csf_PD, cv_outward_dists, 1,
                       PIAL_VERTICES, "csf PD", NULL, ep->dstep, ep->max_dist);
#endif
#if 0
  MRISclearFixedValFlags(mris) ;
  vno = MRISsubsampleDist(mris, subsample_dist) ;
  printf("computing optimal parameters for %d vertices...\n", vno) ;
  dstep = ep->dstep ;
  ep->dstep = 0.25 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (FEQUAL(v->val, 1))
      v->fixedval = TRUE ;
  }
#endif

  for (nholes = 0, vno = start_vno ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (vno > MAX_VNO)
      break ;
    if (((vno+1) % (mris->nvertices/25)) == 0) {
      printf("%2.0f%% ", 100.0f*(float)(vno+1)/(float)(mris->nvertices)) ;
      fflush(stdout) ;
    }
#if 0
    if (v->fixedval != TRUE )  /* not in subsampled surface */
      continue ;
#endif


    dstep = ep->dstep ;
    ep->dstep = 0.5 ;
    sse = compute_optimal_parameters(mris, vno, ep, &best_white_delta,
                                     &best_pial_delta) ;
    ep->dstep = dstep ;
    if (vno == Gdiag_no) {
      printf("v %d: wm = (%2.0f, %2.0f), gm = (%2.0f, %2.0f), "
             "csf = (%2.0f, %2.0f), "
             "deltas = (%2.2f, %2.2f)\n",
             vno,ep->cv_wm_T1[vno], ep->scale*ep->cv_wm_PD[vno], ep->cv_gm_T1[vno],
             ep->scale*ep->cv_gm_PD[vno], ep->cv_csf_T1[vno],
             ep->scale*ep->cv_csf_PD[vno],
             best_white_delta, best_pial_delta);
      DiagBreak() ;
    }
  }

#else
  HISTOGRAM *h_prior, *h_tmp ;

  h_prior = HISTOalloc(MAX_BINS) ;
  h_tmp = HISTOalloc(MAX_BINS) ;
  if (mri_T1) {
    if (fix_T1) {
      int vno ;

      for (vno = 0 ; vno < mris->nvertices ; vno++) {
        cv_gm_T1[vno] = GM_T1 ;
        cv_wm_T1[vno] = WM_T1 ;
      }
    } else {
      for (bno = 0 ; bno < MAX_BINS ; bno++)
        h_tmp->bins[bno] = BIN_SIZE*bno+0.5*BIN_SIZE ;

      for (bno = 0 ; bno < MAX_BINS ; bno++)
        if ((h_tmp->bins[bno] > WM_T1-WM_T1_STD) && (h_tmp->bins[bno] < WM_T1+WM_T1_STD))
          h_tmp->counts[bno] = 1 ;
      h_prior = HISTOsmooth(h_tmp, NULL, WM_T1_STD/BIN_SIZE) ;
      sample_parameter_map(mris, mri_T1, mri_res, cv_wm_T1, cv_inward_dists, -1,
                           ORIGINAL_VERTICES, "wm T1", h_prior, ep->dstep, ep->max_dist);
      HISTOclearCounts(h_tmp, h_tmp) ;
      for (bno = 0 ; bno < MAX_BINS ; bno++)
        if ((h_tmp->bins[bno] > GM_T1-GM_T1_STD) && (h_tmp->bins[bno] < GM_T1+GM_T1_STD))
          h_tmp->counts[bno] = 1 ;
      HISTOsmooth(h_tmp, h_prior, GM_T1_STD/BIN_SIZE) ;
      sample_parameter_map(mris, mri_T1, mri_res, cv_gm_T1, NULL, 1,
                           ORIGINAL_VERTICES, "gm T1", h_prior, ep->dstep, ep->max_dist);
    }
    sample_parameter_map(mris, mri_T1, mri_res, cv_csf_T1, cv_outward_dists, 1,
                         PIAL_VERTICES, "csf T1", NULL, ep->dstep, ep->max_dist);
  } else {
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      cv_gm_T1[vno] = GM_T1 ;
      cv_wm_T1[vno] = WM_T1 ;
      cv_csf_T1[vno] = CSF_T1 ;
    }
  }
  if (mri_PD) {
    sample_parameter_map(mris, mri_PD, mri_res, cv_wm_PD, cv_inward_dists, -1,
                         ORIGINAL_VERTICES, "wm PD", NULL, ep->dstep, ep->max_dist);
    sample_parameter_map(mris, mri_PD, mri_res, cv_gm_PD, NULL, 1,
                         ORIGINAL_VERTICES, "gm PD", NULL, ep->dstep, ep->max_dist);
    sample_parameter_map(mris, mri_PD, mri_res, cv_csf_PD, cv_outward_dists, 1,
                         PIAL_VERTICES, "csf PD", NULL, ep->dstep, ep->max_dist);
  } else {
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      cv_gm_PD[vno] = GM_PD ;
      cv_wm_PD[vno] = WM_PD ;
      cv_csf_PD[vno] = CSF_PD ;
    }
  }


  HISTOfree(&h_tmp) ;
  HISTOfree(&h_prior) ;
#endif
  return(NO_ERROR) ;
}

#if 0
static int
soap_bubble_maps(MRI_SURFACE *mris, EXTRA_PARMS *ep, int smooth_parms) {
  soap_bubble_map(mris, ep->cv_wm_T1, smooth_parms) ;
  soap_bubble_map(mris, ep->cv_wm_PD, smooth_parms) ;
  soap_bubble_map(mris, ep->cv_gm_T1, smooth_parms) ;
  soap_bubble_map(mris, ep->cv_gm_PD, smooth_parms) ;
#if 0
  soap_bubble_map(mris, ep->cv_csf_T1, smooth_parms) ;
  soap_bubble_map(mris, ep->cv_csf_PD, smooth_parms) ;
#endif
  return(NO_ERROR) ;
}


static int
soap_bubble_map(MRI_SURFACE *mris, float *cv, int navgs) {
  MRISimportValVector(mris, cv) ;
  MRISsoapBubbleVals(mris, navgs) ;
  MRISexportValVector(mris, cv) ;
  return(NO_ERROR) ;
}
#endif

static int
smooth_csf_map(MRI_SURFACE *mris, float *cv_T1, float *cv_PD, int navgs) {
  int    i, vno, vnb, vnum, n_vno ;
  float  num, T1, PD, T1_nbr, PD_nbr, T1_avg, PD_avg ;

  for (i = 0 ; i < navgs ; i++) {
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag)
        continue ;
      T1 = cv_T1[vno] ;
      PD = cv_PD[vno] ;
      v->tdx = PD ;
      v->tdy = T1 ;
      if ((PD < MIN_NONBRAIN_PD) || (T1 < MIN_NONBRAIN_T1))
        continue ;
      int const * pnb = vt->v ;
      vnum = vt->vnum ;
      T1_avg = T1 ;
      PD_avg = PD ;
      for (num = 1.0f, vnb = 0 ; vnb < vnum ; vnb++) {
        n_vno = *pnb++ ;
        VERTEX const * const vn = &mris->vertices[n_vno] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        T1_nbr = cv_T1[n_vno] ;
        PD_nbr = cv_PD[n_vno] ;
        if ((PD_nbr < MIN_NONBRAIN_PD) || (T1_nbr < MIN_NONBRAIN_T1))
          continue ;
        T1_avg += T1_nbr ;
        PD_avg += PD_nbr ;
        num++ ;
      }
      v->tdx = PD_avg / num ;
      v->tdy = T1_avg / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      VERTEX const * const v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      cv_T1[vno] = v->tdy ;
      cv_PD[vno] = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}

static int
smooth_map(MRI_SURFACE *mris, float *cv, int navgs) {
  MRISimportCurvatureVector(mris, cv) ;
#if  0
  MRISmedianFilterCurvature(mris,  navgs) ;
#else
  MRISaverageCurvatures(mris, navgs) ;
#endif
  MRISextractCurvatureVector(mris, cv) ;
  return(NO_ERROR) ;
}

static int
smooth_marked_csf_map(MRI_SURFACE *mris, float *cv_T1, float *cv_PD, int navgs) {
  int    i, vno, vnb, vnum, n_vno ;
  float  num, T1, PD, T1_nbr, PD_nbr, T1_avg, PD_avg ;

  for (i = 0 ; i < navgs ; i++) {
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->marked == 0)
        continue ;

      T1 = cv_T1[vno] ;
      PD = cv_PD[vno] ;
      v->tdx = PD ;
      v->tdy = T1 ;
      if ((PD < MIN_NONBRAIN_PD) || (T1 < MIN_NONBRAIN_T1))
        continue ;
      int const * pnb = vt->v ;
      vnum = vt->vnum ;
      T1_avg = T1 ;
      PD_avg = PD ;
      for (num = 1.0f, vnb = 0 ; vnb < vnum ; vnb++) {
        n_vno = *pnb++ ;
        VERTEX const * const vn = &mris->vertices[n_vno] ;    /* neighboring vertex pointer */
        if (vn->ripflag)
          continue ;
        T1_nbr = cv_T1[n_vno] ;
        PD_nbr = cv_PD[n_vno] ;
        if ((PD_nbr < MIN_NONBRAIN_PD) || (T1_nbr < MIN_NONBRAIN_T1))
          continue ;
        T1_avg += T1_nbr ;
        PD_avg += PD_nbr ;
        num++ ;
      }
      v->tdx = PD_avg / num ;
      v->tdy = T1_avg / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      VERTEX const * const v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked == 0)
        continue ;
      cv_T1[vno] = v->tdy ;
      cv_PD[vno] = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}

static int
smooth_marked_map(MRI_SURFACE *mris, float *cv, int navgs) {
  MRISimportCurvatureVector(mris, cv) ;
  MRISaverageMarkedCurvatures(mris, navgs) ;
  MRISextractCurvatureVector(mris, cv) ;
  return(NO_ERROR) ;
}

static int
write_map(MRI_SURFACE *mris, float *cv, const char *name, int suffix, const char *output_suffix) {
  char fname[STRLEN] ;

  if (suffix >= 0)
    sprintf(fname, "%s.%s_%d%s",
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", name, suffix, output_suffix) ;
  else
    sprintf(fname, "%s.%s%s",
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", name, output_suffix) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("writing vector to curvature file %s...\n",fname) ;
  MRISimportCurvatureVector(mris, cv) ;
  MRISwriteCurvature(mris, fname) ;
  return(NO_ERROR) ;
}

static int
smooth_marked_maps(MRI_SURFACE *mris, EXTRA_PARMS *ep, int smooth_parms) {
  smooth_marked_map(mris, ep->cv_wm_T1, smooth_parms) ;
  smooth_marked_map(mris, ep->cv_wm_PD, smooth_parms) ;
  smooth_marked_map(mris, ep->cv_gm_T1, smooth_parms) ;
  smooth_marked_map(mris, ep->cv_gm_PD, smooth_parms) ;
#if 1
  smooth_marked_csf_map(mris, ep->cv_csf_T1, ep->cv_csf_PD, smooth_parms) ;
#endif
  return(NO_ERROR) ;
}
static int
smooth_maps(MRI_SURFACE *mris, EXTRA_PARMS *ep, int smooth_parms) {
  smooth_map(mris, ep->cv_wm_T1, smooth_parms) ;
  smooth_map(mris, ep->cv_wm_PD, smooth_parms) ;
  smooth_map(mris, ep->cv_gm_T1, smooth_parms) ;
  smooth_map(mris, ep->cv_gm_PD, smooth_parms) ;
#if 1
  smooth_csf_map(mris, ep->cv_csf_T1, ep->cv_csf_PD, smooth_parms) ;
#endif
  return(NO_ERROR) ;
}

static int
write_maps(MRI_SURFACE *mris, EXTRA_PARMS *ep, int suffix, const char *output_suffix) {
  write_map(mris, ep->cv_inward_dists, "inward_dists", suffix, output_suffix) ;
  write_map(mris, ep->cv_outward_dists, "outward_dists", suffix, output_suffix) ;
  write_map(mris, ep->cv_wm_T1, "wm_T1", suffix, output_suffix) ;
  write_map(mris, ep->cv_wm_PD, "wm_PD", suffix, output_suffix) ;
  write_map(mris, ep->cv_gm_T1, "gm_T1", suffix, output_suffix) ;
  write_map(mris, ep->cv_gm_PD, "gm_PD", suffix, output_suffix) ;
  write_map(mris, ep->cv_csf_T1, "csf_T1", suffix, output_suffix) ;
  write_map(mris, ep->cv_csf_PD, "csf_PD", suffix, output_suffix) ;
  return(NO_ERROR) ;
}

#define HALF_VOXEL_SIZE 0.0005

static double
image_forward_model(double TR, double flip_angle, double dist, double thickness,
                    double T1_wm, double PD_wm, double T1_gm, double PD_gm,
                    double T1_csf, double PD_csf) {
  static double   predicted_white, predicted_gray, predicted_csf, last_TR=0,
      last_flip_angle=0, last_T1_wm, last_PD_wm, last_T1_gm, last_PD_gm, last_T1_csf,
                      last_PD_csf;
  double predicted_val, wt ;

  if (!FZERO(last_TR - TR) || !FZERO(last_flip_angle - flip_angle) ||
      !FZERO(last_T1_wm-T1_wm) || !FZERO(last_PD_wm-PD_wm) ||
      !FZERO(last_T1_gm-T1_gm) || !FZERO(last_PD_gm-PD_gm) ||
      !FZERO(last_T1_csf-T1_csf) || !FZERO(last_PD_csf-PD_csf)) {
    last_flip_angle = flip_angle ;
    last_TR = TR ;
    last_T1_wm = T1_wm ;
    last_PD_wm = PD_wm ;
    last_T1_gm = T1_gm ;
    last_PD_gm = PD_gm ;
    last_T1_csf = T1_csf ;
    last_PD_csf = PD_csf ;

    predicted_white = FLASHforwardModelLookup(flip_angle, TR, PD_wm, T1_wm) ;
    predicted_gray = FLASHforwardModelLookup(flip_angle, TR, PD_gm, T1_gm) ;
    predicted_csf = FLASHforwardModelLookup(flip_angle, TR, PD_csf, T1_csf) ;
  }

  if (dist < -HALF_VOXEL_SIZE)   /* in white matter */
    return(predicted_white) ;
  else if (dist > thickness+HALF_VOXEL_SIZE)
    return(predicted_csf) ;
  else if (dist > HALF_VOXEL_SIZE && dist < thickness-HALF_VOXEL_SIZE)
    return(predicted_gray) ;
  else if (dist <= HALF_VOXEL_SIZE)   /* mixure of white and gray */
  {
    wt = (HALF_VOXEL_SIZE - dist) / HALF_VOXEL_SIZE ;
    predicted_val = wt*predicted_white + (1.0-wt)*predicted_gray ;
  } else  /* mixture of gray and csf */
  {
    dist -= thickness ;
    wt = (HALF_VOXEL_SIZE - dist) / HALF_VOXEL_SIZE ;
    predicted_val = wt*predicted_gray + (1.0-wt)*predicted_csf ;
  }

  return(predicted_val) ;
}

static int
find_nearest_pial_vertices(MRI_SURFACE *mris, int *nearest_pial_vertices,
                           int *nearest_white_vertices) {
#if 1
  int   vno ;

  for (vno = 0 ; vno < mris->max_vertices ; vno++)
    nearest_pial_vertices[vno] = nearest_white_vertices[vno] = vno ;
#else
  int     vno, n, vlist[100000], vtotal, ns, i, vnum, nbr_count[100], min_n ;
  VERTEX  *v, *vn, *vn2 ;
  float   dx, dy, dz, dist, min_dist, nx, ny, nz, dot ;

  memset(nbr_count, 0, 100*sizeof(int)) ;

  /* pial vertex positions are gray matter, orig are white matter */
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    dx = v->pialx - v->origx ;
    dy = v->pialy - v->origy ;
    dz = v->pialz - v->origz ;
    min_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    min_n = 0 ;
    v->marked = 1 ;
    vtotal = 1 ;
    vlist[0] = vno ;
    nearest_pial_vertices[vno] = vno ;
    min_n = 0 ;
    for (ns = 1 ; ns <= nbhd_size ; ns++) {
      vnum = 0 ;  /* will be # of new neighbors added to list */
      for (i = 0 ; i < vtotal ; i++) {
        vn = &mris->vertices[vlist[i]] ;
        if (vn->ripflag)
          continue ;
        if (vn->marked && vn->marked < ns-1)
          continue ;
        for (n = 0 ; n < vn->vnum ; n++) {
          vn2 = &mris->vertices[vn->v[n]] ;
          if (vn2->ripflag || vn2->marked)  /* already processed */
            continue ;
          vlist[vtotal+vnum++] = vn->v[n] ;
          vn2->marked = ns ;
          dx = vn2->pialx-v->origx ;
          dy = vn2->pialy-v->origy ;
          dz = vn2->pialz-v->origz ;
          dot = dx*nx + dy*ny + dz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dot = vn2->nx*nx + vn2->ny*ny + vn2->nz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dist = sqrt(dx*dx + dy*dy + dz*dz) ;
          if (dist < min_dist) {
            min_n = ns ;
            min_dist = dist ;
            nearest_pial_vertices[vno] = vn->v[n] ;
            if (min_n == nbhd_size && DIAG_VERBOSE_ON)
              fprintf(stdout, "%d --> %d = %2.3f\n",
                      vno,vn->v[n], dist) ;
          }
        }
      }
      vtotal += vnum ;
    }

    nearest_white_vertices[nearest_pial_vertices[vno]] = vno ;
    nbr_count[min_n]++ ;
    for (n = 0 ; n < vtotal ; n++) {
      vn = &mris->vertices[vlist[n]] ;
      if (vn->ripflag)
        continue ;
      vn->marked = 0 ;
    }
  }


  for (n = 0 ; n <= nbhd_size ; n++)
    printf("%d vertices at %d distance\n", nbr_count[n], n) ;
#endif
  return(NO_ERROR) ;
}
static int
find_nearest_white_vertices(MRI_SURFACE *mris, int *nearest_white_vertices) {
  int     vno, n, vlist[100000], vtotal, ns, i, vnum, nbr_count[100], min_n, nfilled ;
  float   dx, dy, dz, dist, min_dist, nx, ny, nz, dot ;

  memset(nbr_count, 0, 100*sizeof(int)) ;

  /* pial vertex positions are gray matter, orig are white matter */
  for (nfilled = 0, vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX * const v  = &mris->vertices[vno];
    if (v->ripflag)
      continue ;
    if (nearest_white_vertices[vno] >= 0)
      continue ;
    nfilled++ ;
    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    dx = v->origx - v->pialx ;
    dy = v->origy - v->pialy ;
    dz = v->origz - v->pialz ;
    min_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    nearest_white_vertices[vno] = vno ;
    min_n = 0 ;
    v->marked = 1 ;
    vtotal = 1 ;
    vlist[0] = vno ;
    min_n = 0 ;
    for (ns = 1 ; ns <= nbhd_size ; ns++) {
      vnum = 0 ;  /* will be # of new neighbors added to list */
      for (i = 0 ; i < vtotal ; i++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vlist[i]];
        VERTEX                * const vn  = &mris->vertices         [vlist[i]] ;
        if (vn->ripflag)
          continue ;
        if (vn->marked && vn->marked < ns-1)
          continue ;
        for (n = 0 ; n < vnt->vnum ; n++) {
          VERTEX * const vn2 = &mris->vertices[vnt->v[n]] ;
          if (vn2->ripflag || vn2->marked)  /* already processed */
            continue ;
          vlist[vtotal+vnum++] = vnt->v[n] ;
          vn2->marked = ns ;
          dx = vn2->origx-v->pialx ;
          dy = vn2->origy-v->pialy ;
          dz = vn2->origz-v->pialz ;
          dot = dx*nx + dy*ny + dz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dot = vn2->nx*nx + vn2->ny*ny + vn2->nz*nz ;
          if (dot < 0) /* must be outwards from surface */
            continue ;
          dist = sqrt(dx*dx + dy*dy + dz*dz) ;
          if (dist < min_dist) {
            min_n = ns ;
            min_dist = dist ;
            nearest_white_vertices[vno] = vnt->v[n] ;
            if (min_n == nbhd_size && DIAG_VERBOSE_ON)
              fprintf(stdout, "%d --> %d = %2.3f\n",
                      vno,vnt->v[n], dist) ;
          }
        }
      }
      vtotal += vnum ;
    }

    nbr_count[min_n]++ ;
    for (n = 0 ; n < vtotal ; n++) {
      VERTEX * const vn = &mris->vertices[vlist[n]] ;
      if (vn->ripflag)
        continue ;
      vn->marked = 0 ;
    }
  }

#if 0
  printf("%d holes filled in nearest wm vertex finding...\n", nfilled) ;
  for (n = 0 ; n <= nbhd_size ; n++)
    printf("%d vertices at %d distance\n", nbr_count[n], n) ;
#endif
  return(NO_ERROR) ;
}

static double
compute_optimal_parameters(MRI_SURFACE *mris, int vno,
                           EXTRA_PARMS *ep, double *pwhite_delta,
                           double *ppial_delta) {
  double       predicted_val, dx, dy, dz,
  dist, inward_dist, outward_dist, cortical_dist, PD_wm, PD_gm, PD_csf,
  T1_wm, T1_gm, T1_csf, white_dist, pial_dist, sse, best_sse,
  best_white_dist, best_pial_dist, orig_cortical_dist, total_dist,
  orig_white_index, orig_pial_index, sigma, best_T1_wm, best_PD_wm,
  best_T1_gm, best_PD_gm, best_T1_csf, best_PD_csf, T1, PD  ;
  int          i, j, white_vno, pial_vno, max_j, best_white_index, csf_len, bad,
  best_pial_index, pial_index, white_index, max_white_index, best_csf_len;
  VERTEX       *v_white, *v_pial ;
  MRI          *mri, *mri_T1, *mri_PD ;
  double       image_vals[MAX_FLASH_VOLUMES][MAX_SAMPLES], xw, yw, zw, xp, yp, zp,
  x, y, z, T1_vals[MAX_SAMPLES], PD_vals[MAX_SAMPLES],
  best_image_vals[MAX_FLASH_VOLUMES][MAX_SAMPLES] ;

  if (vno == Gdiag_no)
    DiagBreak() ;

  mri_T1 = ep->mri_T1 ;
  mri_PD = ep->mri_PD ;
  white_vno = vno ;
#if 0
  pial_vno = ep->nearest_pial_vertices[vno] ;
#else
  pial_vno = vno ;
#endif
  v_white = &mris->vertices[white_vno] ;
  v_pial = &mris->vertices[pial_vno] ;
#if 1
  inward_dist = ep->cv_inward_dists[white_vno] ;
  outward_dist = ep->cv_outward_dists[white_vno] ;
#else
  inward_dist = ep->max_inward_dist ;
  outward_dist = ep->max_outward_dist ;
#endif
  sigma = ep->current_sigma ;

  // MRIworldToVoxel(ep->mri_flash[0], v_white->origx, v_white->origy,v_white->origz,
  //                &xw,&yw,&zw);
  // MRIworldToVoxel(ep->mri_flash[0], v_pial->pialx, v_pial->pialy,v_pial->pialz,
  //                 &xp,&yp,&zp);
  MRIsurfaceRASToVoxel(ep->mri_flash[0], v_white->origx, v_white->origy,v_white->origz,
                       &xw,&yw,&zw);
  MRIsurfaceRASToVoxel(ep->mri_flash[0], v_pial->pialx, v_pial->pialy,v_pial->pialz,
                       &xp,&yp,&zp);

  dx = xp - xw ;
  dy = yp - yw ;
  dz = zp - zw ;
  orig_cortical_dist = cortical_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
  if (TOO_SMALL(cortical_dist)) {
    dx = v_pial->nx ;
    dy = v_pial->ny ;
    dz = v_pial->nz ;
    // MRIworldToVoxel(ep->mri_flash[0],
    //                 v_pial->pialx+v_pial->nx,
    //                v_pial->pialy+v_pial->ny,
    //                v_pial->pialz+v_pial->nz,
    //                &xp,&yp,&zp);
    MRIsurfaceRASToVoxel(ep->mri_flash[0],
                         v_pial->pialx+v_pial->nx,
                         v_pial->pialy+v_pial->ny,
                         v_pial->pialz+v_pial->nz,
                         &xp,&yp,&zp);
    dx = xp - xw ;
    dy = yp - yw ;
    dz = zp - zw ;
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (FZERO(dist))
      dist = 1.0 ;
    dx /= dist ;
    dy /= dist ;
    dz /= dist ;
    xp = xw ;
    yp = yw ;
    zp = zw ;   /* true surface position */
  } else {
    /* surface normal in volume coords */
    dx /= cortical_dist ;
    dy /= cortical_dist ;
    dz /= cortical_dist ;
  }

  /* fill image_vals[][] array */
  orig_white_index = inward_dist / ep->dstep ;
  orig_pial_index = orig_white_index + cortical_dist / ep->dstep ;
  max_j = (int)((inward_dist + cortical_dist + outward_dist) / ep->dstep) ;
  for (j = 0, dist = -inward_dist ; dist <= cortical_dist+outward_dist ;
       dist += ep->dstep, j++) {
    x = xw + dist * dx ;
    y = yw + dist * dy ;
    z = zw + dist * dz ;
    for (i = 0 ; i < ep->nvolumes ; i++) {
      mri = ep->mri_flash[i] ;
#if 0
      MRIsampleVolumeDirectionScale(mri, x, y, z, dx, dy, dz, &image_vals[i][j], sigma) ;
#else
      MRIsampleVolumeType(mri, x, y, z, &image_vals[i][j], sample_type) ;
      if (j > orig_pial_index)
        DiagBreak() ;
      if (!std::isfinite(image_vals[i][j]))
        DiagBreak() ;
#endif
    }
    MRIsampleVolumeType(mri_T1, x, y, z, &T1_vals[j], sample_type) ;
    MRIsampleVolumeType(mri_PD, x, y, z, &PD_vals[j], sample_type) ;
  }
  if ((plot_stuff == 1 && Gdiag_no == vno) || plot_stuff > 1) {
    FILE *fp ;
    char fname[STRLEN] ;

    for (i = 0 ; i < ep->nvolumes ; i++) {
      sprintf(fname, "vol%d.plt", i) ;
      fp = fopen(fname, "w") ;
      for (j = 0, dist = -inward_dist ; dist <= cortical_dist+outward_dist ;
           dist += ep->dstep, j++)
        fprintf(fp, "%d %f %f\n", j, dist, image_vals[i][j]) ;
      fclose(fp) ;
    }
  }
  if ((plot_stuff == 1 && Gdiag_no == vno) ||
      plot_stuff > 1) {
    FILE *fp ;

    fp = fopen("orig_white.plt", "w") ;
    fprintf(fp, "%f %f %f\n%f %f %f\n",
            orig_white_index, 0.0, 0.0, orig_white_index, 0.0, 120.0) ;
    fclose(fp) ;
    fp = fopen("orig_pial.plt", "w") ;
    fprintf(fp, "%f %f %f\n%f %f %f\n",
            orig_pial_index, cortical_dist, 0.0, orig_pial_index, cortical_dist, 120.0) ;
    fclose(fp) ;
  }


  total_dist = outward_dist + orig_cortical_dist + inward_dist ;
  best_white_index = orig_white_index ;
  best_pial_index = orig_pial_index ;
  best_csf_len = outward_dist*ep->dstep ;
  best_white_dist = orig_white_index*ep->dstep;
  best_pial_dist = orig_pial_index*ep->dstep;
#if 0
  max_white_index = nint((inward_dist + orig_cortical_dist/2)/ep->dstep) ;
#else
  max_white_index = max_j - 1/ep->dstep ;
#endif
  cortical_dist = best_pial_dist - best_white_dist ;
  best_sse = -1 ;
  best_T1_wm = WM_T1 ;
  best_PD_wm = WM_PD ;
  best_T1_gm = GM_T1 ;
  best_PD_gm = GM_PD ;
  best_T1_csf = CSF_T1 ;
  best_PD_csf = CSF_PD ;

  for (white_index = 0 ; white_index <= max_white_index ; white_index++) {
    white_dist = white_index * ep->dstep ;
    /*    for (pial_index = white_index + nint(1.0/ep->dstep) ; pial_index <= max_j ;pial_index++)*/
    for (pial_index = white_index+1 ; pial_index <= max_j ;pial_index++) {
      /*
        for this pair of white/pial offsets, compute the mean wm,gm,csf vals
        and the T1/PD pairs that match them. Then compute the sse of the
        total and see if it is better than the best previous set of parms.
      */
      pial_dist = pial_index * ep->dstep ;
      cortical_dist = pial_dist - white_dist ;
      csf_len = 0 ;
      bad = 0 ;

      /* compute mean wm for these positions */
      if (white_index < 1)   /* just sample 1mm inwards */
      {
        T1_wm = T1_vals[white_index] ;
        PD_wm = PD_vals[white_index] ;
        if (T1_wm < ep->cv_min_wm_T1[vno] || T1_wm > ep->cv_max_wm_T1[vno])
          bad = 1 ;
      } else  /* average inwards up to 1/2 mm from white position */
      {
        for (T1_wm = PD_wm = 0.0, j = 0 ; j < white_index ; j++) {
          T1 = T1_vals[j] ;
          PD = PD_vals[j] ;
          if (T1 < ep->cv_min_wm_T1[vno] || T1 > ep->cv_max_wm_T1[vno])
            bad = 1 ;
          if (T1 < ep->cv_min_wm_T1[vno])
            T1 = ep->cv_min_wm_T1[vno] ;
          if (T1 > ep->cv_max_wm_T1[vno])
            T1 = ep->cv_max_wm_T1[vno] ;
          if (PD < ep->cv_min_wm_PD[vno])
            PD = ep->cv_min_wm_PD[vno] ;
          T1_wm += T1 ;
          PD_wm += PD ;
        }

        PD_wm /= (double)white_index ;
        T1_wm /= (double)white_index ;
      }

      /* compute mean gm for these positions */
      if (pial_index - white_index <= 1) {
        for (T1_gm = PD_gm = 0.0, j = white_index ; j <= pial_index ; j++) {
          T1 = T1_vals[j] ;
          PD = PD_vals[j] ;
          T1_gm += T1 ;
          PD_gm += PD ;
        }
        T1_gm /= (double)(pial_index-white_index+1) ;
        PD_gm /= (double)(pial_index-white_index+1) ;
        if (T1_gm < ep->cv_min_gm_T1[vno] || T1_gm > ep->cv_max_gm_T1[vno])
          bad = 1 ;
        if (PD_gm < ep->cv_min_gm_PD[vno])
          PD_gm = ep->cv_min_gm_PD[vno] ;
#if 0
        if (PD_gm > ep->cv_max_gm_PD[vno])
          PD_gm = ep->cv_max_gm_PD[vno] ;
#endif
      } else {
        for (T1_gm = PD_gm = 0.0, j = white_index+1 ; j < pial_index ; j++) {
          T1 = T1_vals[j] ;
          PD = PD_vals[j] ;
          if (T1 < ep->cv_min_gm_T1[vno] || T1 > ep->cv_max_gm_T1[vno])
            bad = 1 ;
          if (T1 < ep->cv_min_gm_T1[vno])
            T1 = ep->cv_min_gm_T1[vno] ;
          if (T1 > ep->cv_max_gm_T1[vno])
            T1 = ep->cv_max_gm_T1[vno] ;
          if (PD < ep->cv_min_gm_PD[vno])
            PD = ep->cv_min_gm_PD[vno] ;
          T1_gm += T1 ;
          PD_gm += PD ;
        }
        T1_gm /= (double)(pial_index-white_index-1) ;
        PD_gm /= (double)(pial_index-white_index-1) ;
      }
      /* compute mean csf for these positions */
      if (max_j - pial_index < 1) {
        T1_csf = T1_vals[max_j] ;
        PD_csf = PD_vals[max_j] ;
      } else {
        for (T1_csf = PD_csf = 0.0, j = pial_index+1 ; j <= max_j ; j++) {
          T1 = T1_vals[j] ;
          PD = PD_vals[j] ;
          T1_csf += T1 ;
          PD_csf += PD ;
        }
        T1_csf /= (double)(max_j - pial_index) ;
        PD_csf /= (double)(max_j - pial_index) ;
      }

      if (bad)   /* disallow certain T1/PD values for gm and wm */
        continue ;

      /* do some bounds checking */
      if (PD_wm < ep->cv_min_wm_PD[vno])
        PD_wm = ep->cv_min_wm_PD[vno] ;
      if (T1_wm < ep->cv_min_wm_T1[vno])
        T1_wm = ep->cv_min_wm_T1[vno] ;
      if (T1_wm > ep->cv_max_wm_T1[vno])
        T1_wm = ep->cv_max_wm_T1[vno];
      if (T1_gm < ep->cv_min_gm_T1[vno])
        T1_gm = ep->cv_min_gm_T1[vno] ;
      if (T1_gm > ep->cv_max_gm_T1[vno])
        T1_gm = ep->cv_max_gm_T1[vno]  ;
      if (PD_csf > MIN_NONBRAIN_PD && T1_csf < MIN_CSF_T1)
        T1_csf = MIN_CSF_T1 ;
      if (PD_gm < ep->cv_min_gm_PD[vno])
        PD_gm = ep->cv_min_gm_PD[vno] ;
      if ((PD_csf > MIN_NONBRAIN_PD) && (T1_gm > T1_csf*.9))
        T1_csf = T1_gm/0.9 ;
      if (T1_wm > T1_gm*.9)
        T1_wm = T1_gm *.9 ;
      sse = compute_vertex_sse(ep, image_vals, max_j, white_dist,
                               cortical_dist,
                               T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf, 0,
                               T1_vals, PD_vals, vno) ;

      /* do some bounds checking */
      if (T1_wm < ep->cv_min_wm_T1[vno] || T1_wm > ep->cv_max_wm_T1[vno] ||
          T1_gm < ep->cv_min_gm_T1[vno] || T1_gm > ep->cv_max_gm_T1[vno])
        continue ;

      if (!std::isfinite(sse))
        ErrorPrintf(ERROR_BADPARM, "sse not finite at v %d", vno) ;

      if (sse < best_sse || best_sse < 0) {
        best_csf_len = csf_len ;
        best_white_index = white_index ;
        best_pial_index = pial_index ;
        best_pial_dist = pial_dist ;
        best_white_dist = white_dist ;
        best_sse = sse ;
        best_T1_wm = T1_wm ;
        best_PD_wm = PD_wm ;
        best_T1_gm = T1_gm ;
        best_PD_gm = PD_gm ;
        best_T1_csf = T1_csf ;
        best_PD_csf = PD_csf ;
        mris->vertices[vno].marked = 1 ;
        for (i = 0 ; i < ep->nvolumes ; i++) {
          for (j = 0 ; j <= max_j ; j++)
            best_image_vals[i][j] = image_vals[i][j] ;
        }
      }
    }
  }


  if ((plot_stuff == 1 && Gdiag_no == vno) ||
      plot_stuff > 1) {
    FILE *fp ;
    char fname[STRLEN] ;

    cortical_dist = best_pial_dist - best_white_dist ;
    for (i = 0 ; i < ep->nvolumes ; i++) {
      sprintf(fname, "vol%d_best.plt", i) ;
      fp = fopen(fname, "w") ;
      mri = ep->mri_flash[i] ;
      for (j = 0 ; j <= max_j ; j++) {
        dist = (double)j*ep->dstep-best_white_dist ;
        predicted_val =
          image_forward_model(mri->tr, mri->flip_angle, dist, cortical_dist,
                              best_T1_wm, best_PD_wm, best_T1_gm, best_PD_gm, best_T1_csf,
                              best_PD_csf) ;
        fprintf(fp, "%d %f %f\n", j, dist, predicted_val) ;
      }
      fclose(fp) ;
    }
  }

  if ((plot_stuff == 1 && Gdiag_no == vno) ||
      plot_stuff > 1) {
    FILE *fp ;

    fp = fopen("best_white.plt", "w") ;
    fprintf(fp, "%d %f %f\n%d %f %f\n",
            best_white_index, 0.0, 0.0, best_white_index, 0.0, 120.0) ;
    fclose(fp) ;
    fp = fopen("best_pial.plt", "w") ;
    fprintf(fp, "%d %f %f\n%d %f %f\n",
            best_pial_index, cortical_dist, 0.0, best_pial_index, cortical_dist, 120.0) ;
    fclose(fp) ;
  }


  *pwhite_delta = best_white_dist - (orig_white_index * ep->dstep) ;
  *ppial_delta =  best_pial_dist  - (orig_pial_index  * ep->dstep) ;
  if (*pwhite_delta < 0.1 || *ppial_delta > 0.1)
    DiagBreak() ;
  if (vno == Gdiag_no) {
    double pwhite[MAX_FLASH_VOLUMES], pgray[MAX_FLASH_VOLUMES], pcsf[MAX_FLASH_VOLUMES];

    printf("\n") ;
    for (i = 0 ; i < ep->nvolumes ; i++) {
      mri = ep->mri_flash[i] ;
      pwhite[i] =
        image_forward_model(mri->tr, mri->flip_angle, -1, cortical_dist,
                            best_T1_wm, best_PD_wm, best_T1_gm, best_PD_gm,
                            best_T1_csf, best_PD_csf) ;
      pgray[i] =
        image_forward_model(mri->tr, mri->flip_angle, cortical_dist/2,
                            cortical_dist,
                            best_T1_wm, best_PD_wm, best_T1_gm, best_PD_gm,
                            best_T1_csf, best_PD_csf) ;
      pcsf[i] =
        image_forward_model(mri->tr, mri->flip_angle, cortical_dist+1,
                            cortical_dist, best_T1_wm, best_PD_wm,
                            best_T1_gm, best_PD_gm, best_T1_csf, best_PD_csf) ;
      printf("\tpredicted image intensities %d: (%2.0f, %2.0f, %2.0f)\n",
             i, ep->scale*pwhite[i], ep->scale*pgray[i], ep->scale*pcsf[i]) ;
    }

    printf("v %d: best_white_delta = %2.2f, best_pial_delta = %2.2f, best_sse = %2.2f\n",
           vno, *pwhite_delta, *ppial_delta, best_sse) ;
    printf("      inward_dist = %2.2f, outward_dist = %2.2f, cortical_dist = %2.2f\n",
           inward_dist, outward_dist, orig_cortical_dist) ;
    cortical_dist = best_pial_dist - best_white_dist ;
    sse = compute_vertex_sse(ep, image_vals, max_j, best_white_dist,
                             cortical_dist, best_T1_wm, best_PD_wm,
                             best_T1_gm, best_PD_gm, best_T1_csf,
                             best_PD_csf, plot_stuff, T1_vals, PD_vals, vno) ;
    if (!std::isfinite(sse))
      ErrorPrintf(ERROR_BADPARM, "sse not finite at v %d", vno) ;

  }
  if (best_T1_wm < ep->cv_min_wm_T1[vno])
    best_T1_wm = ep->cv_min_wm_T1[vno] ;
  if (best_T1_wm > ep->cv_max_wm_T1[vno])
    best_T1_wm = ep->cv_max_wm_T1[vno] ;
  if (best_T1_gm < ep->cv_min_gm_T1[vno])
    best_T1_gm = ep->cv_min_gm_T1[vno]  ;
  if (best_T1_gm > ep->cv_max_gm_T1[vno])
    best_T1_gm = ep->cv_max_gm_T1[vno] ;
#if 0
  if (best_sse > 0)
    ep->cv_outward_dists[vno] =
      ep->dstep*((float)(best_pial_index+csf_len) - orig_pial_index) ;
#endif

  ep->cv_wm_T1[vno] = best_T1_wm ;
  ep->cv_wm_PD[vno] = best_PD_wm ;
  ep->cv_gm_T1[vno] = best_T1_gm ;
  ep->cv_gm_PD[vno] = best_PD_gm ;
  ep->cv_csf_T1[vno] = best_T1_csf ;
  ep->cv_csf_PD[vno] = best_PD_csf ;
  return(best_sse) ;
}


static double
compute_vertex_sse(EXTRA_PARMS *ep, double image_vals[MAX_FLASH_VOLUMES][MAX_SAMPLES], int max_j,
                   double white_dist, double cortical_dist, double T1_wm, double PD_wm,
                   double T1_gm, double PD_gm, double T1_csf, double PD_csf, int debug,
                   double *T1_vals, double *PD_vals, int vno) {
  double sse, dist, predicted_val, error, scale, T1, PD ;
  MRI    *mri, *mri_T1, *mri_PD ;
  int    i, j ;
  static int callno = 0 ;
  FILE   *fp ;

  mri_T1 = ep->mri_T1 ;
  mri_PD = ep->mri_PD ;
  sse = 0.0 ;
  if (debug) {
    char fname[1000] ;
    sprintf(fname, "sse%d.dat", callno) ;
    fp = fopen(fname, "w") ;
    fprintf(fp, "call %d: T1 = (%2.0f, %2.0f, %2.0f), PD = (%2.0f, %2.0f, %2.0f)\n",
            callno, T1_wm, T1_gm, T1_csf, PD_wm, PD_gm, PD_csf) ;
    fprintf(fp, "max_j = %d, white_dist = %2.2f, cortical_dist = %2.2f\n",
            max_j, white_dist, cortical_dist) ;
    printf("call %d: T1 = (%2.0f, %2.0f, %2.0f), PD = (%2.0f, %2.0f, %2.0f)\n",
           callno, T1_wm, T1_gm, T1_csf, PD_wm, PD_gm, PD_csf) ;
    printf("max_j = %d, white_dist = %2.2f, cortical_dist = %2.2f\n",
           max_j, white_dist, cortical_dist) ;
    callno++ ;
  } else
    fp = NULL ;
  for (i = 0 ; i < ep->nvolumes ; i++) {
    mri = ep->mri_flash[i] ;
    for (j = 0 ; j <= max_j ; j++) {
      dist = (double)j*ep->dstep-white_dist ;
      T1 = T1_vals[j] ;
      PD = PD_vals[j] ;
      if (dist < (cortical_dist-.5) && !IS_BRAIN(T1,PD,vno,ep))
        scale = 10 ;
      else
        scale = 1 ;
      predicted_val =
        image_forward_model(mri->tr, mri->flip_angle, dist, cortical_dist,
                            T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf) ;
      error = scale*(image_vals[i][j] - predicted_val) ;
      sse += (error*error) ;
      if (!std::isfinite(sse)) {
        printf("sse not finite predicted_val=%2.1f, "
               "tissue parms=(%2.0f,%2.0f,%2.0f|%2.0f, %2.0f, %2.0f)\n",
               predicted_val, T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf) ;
        DiagBreak() ;
      }
      if (debug)
        fprintf(fp, "%f %f %f %f %d %d %f\n",
                dist, image_vals[i][j], predicted_val, error, i, j, sse) ;
    }
  }
  sse /= ((max_j+1)*ep->nvolumes) ; /* so it matches vertex_error */
  if (debug) {
    fprintf(fp, "sse = %2.3f\n", sse) ;
    printf("sse = %2.3f\n", sse) ;
    fclose(fp) ;
  }
  return(sse) ;
}




#if 0
static int
sample_parameter_map(MRI_SURFACE *mris, MRI *mri, MRI *mri_res,
                     float *cv_parm, float *cv_dists, float dir,
                     int which_vertices, char *name, HISTOGRAM *h_prior,
                     double dstep, double max_sampling_dist) {
  VERTEX    *v ;
  int       vno, bpeak, bsmooth_peak, bno ;
  float     dist, max_dist, len, dist1, dist2 ;
  double    dx, dy, dz, res, parm_sample, total_wt, wt, x0, y0, z0, x1, y1, z1,
  parm, xs, ys, zs, xe, ye, ze, tx1, ty1, tz1, tx2, ty2, tz2 ;
  HISTOGRAM *histo, *hsmooth ;

  histo = HISTOalloc(MAX_BINS) ;
  hsmooth = HISTOalloc(MAX_BINS) ;
  for (bno = 0 ; bno < MAX_BINS ; bno++)
    histo->bins[bno] = BIN_SIZE*bno+0.5*BIN_SIZE ;

  if (!mri)
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    switch (which_vertices) {
    default:
      ErrorReturn(ERROR_BADPARM,
                  (ERROR_BADPARM, "sample_parameter_map: %d bad vertex set", which_vertices));
      break ;
    case PIAL_VERTICES:
      xs = v->pialx ;
      ys = v->pialy ;
      zs = v->pialz ;
      break ;
    case ORIGINAL_VERTICES:
      xs = v->origx ;
      ys = v->origy ;
      zs = v->origz ;
      break ;
    }

    dx = v->pialx - v->origx ;
    dy = v->pialy - v->origy ;
    dz = v->pialz - v->origz ;
    if (cv_dists)
      max_dist = cv_dists[vno] ;
    else  /* no distance specified - assume it is cortical thickness */
      max_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (h_prior)
      max_dist += 2 ;

    if (FZERO(max_dist))
      max_dist = dstep ;

    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (TOO_SMALL(dist)) {
      dx = v->nx ;
      dy = v->ny ;
      dz = v->nz ;
      dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    }
    if (FZERO(dist))
      dist = dstep ;
    dx /= dist ;
    dy /= dist ;
    dz /= dist ;
    xe = xs + max_dist*dir*dx ;
    ye = ys + max_dist*dir*dy ;
    ze = zs + max_dist*dir*dz ;

    // MRIworldToVoxel(mri, xs, ys, zs, &xs, &ys, &zs) ;
    // MRIworldToVoxel(mri, xe, ye, ze, &xe, &ye, &ze) ;
    MRIsurfaceRASToVoxel(mri, xs, ys, zs, &xs, &ys, &zs) ;
    MRIsurfaceRASToVoxel(mri, xe, ye, ze, &xe, &ye, &ze) ;

    dx = xe-xs ;
    dy = ye-ys ;
    dz = ze-zs ;
    max_dist = sqrt(dx*dx + dy*dy + dz*dz) ;  /* in units of voxels now */
    if (FZERO(max_dist))
      max_dist = dstep ;
    dx /= max_dist ;
    dy /= max_dist ;
    dz /= max_dist ;

    /* avoid boundaries of region */
    xs += dstep*dx ;
    ys += dstep*dy ;
    zs += dstep*dz ;
    max_dist = max_dist - 1.5*dstep ;   /* so we will sample the last point */
    if (max_dist <= dstep)
      max_dist = 1.5*dstep ;

    /*  now compute tangent vectors */
    tx2 = dy ;
    ty2 = dz ;
    tz2 = dx ;  /* any independent vector */
    tx1 = dy*tz2 - dz*ty2 ;           /* cross product */
    ty1 = dz*tx2 - dx*tz2 ;
    tz1 = dx*ty2 - dy*tx2 ;
    len = sqrt(tx1*tx1 + ty1*ty1 + tz1*tz1) ;
    if (FZERO(len))
      len = 1 ;
    tx1 /= len ;
    ty1 /= len ;
    tz1 /= len ;

    tx2 = dy*tz1 - dz*ty1 ;           /* cross product */
    ty2 = dz*tx1 - dx*tz1 ;
    tz2 = dx*ty1 - dy*tx1 ;
    len = sqrt(tx2*tx2 + ty2*ty2 + tz2*tz2) ;
    if (FZERO(len))
      len = 1 ;
    tx2 /= len ;
    ty2 /= len ;
    tz2 /= len ;

    HISTOclearCounts(histo, histo) ;
    total_wt = 0 ;
    parm = 0 ;
    for (dist = 0 ; dist <= max_dist ; dist += dstep) /* along surface normal (sort of) */
    {
      x0 = xs + dist*dx ;
      y0 = ys + dist*dy ;
      z0 = zs + dist*dz ;

      /* sample in plane */
      for (dist1 = -.5 ; dist1 <= 0.5 ; dist1 += dstep) {
        for (dist2 = -.5 ; dist2 <= 0.5 ; dist2 += dstep) {
          x1 = x0 + dist1*tx1 + dist2*tx2 ;
          y1 = y0 + dist1*ty1 + dist2*ty2 ;
          z1 = z0 + dist1*tz1 + dist2*tz2 ;

          MRIsampleVolumeType(mri, x1, y1, z1, &parm_sample, sample_type) ;
          MRIsampleVolumeType(mri_res, x1, y1, z1, &res, sample_type) ;
          wt = 1 / (res + EPSILON) ;   /* avoid dividing by 0 */
          total_wt += wt ;
          parm += (parm_sample*wt) ;
          bno = (int)(parm_sample/BIN_SIZE) ;
          if (bno >= MAX_BINS)
            bno = MAX_BINS-1 ;
          else if (bno < 0)
            bno = 0 ;
          histo->counts[bno]++ ;
        }
      }
    }
    HISTOsmooth(histo, hsmooth, histo_sigma) ;
    if (h_prior) {
      if (Gdiag_no == vno) {
        HISTOplot(hsmooth, "hsmooth_noprior.plt") ;
        bsmooth_peak = HISTOfindHighestPeakInRegion(hsmooth, 0, MAX_BINS) ;
        printf("before application of prior, peak at %2.1f (%d)\n",
               hsmooth->bins[bsmooth_peak], bsmooth_peak) ;
      }
      HISTOmul(histo, h_prior, hsmooth) ;
    }

    bsmooth_peak = HISTOfindHighestPeakInRegion(hsmooth, 0, MAX_BINS) ;
    parm /= total_wt ;
#if 0
    cv_parm[vno] = parm ;
#else
    cv_parm[vno] =  hsmooth->bins[bsmooth_peak] ;
#endif
    if (vno == Gdiag_no) {
      bpeak = HISTOfindHighestPeakInRegion(histo, 0, MAX_BINS) ;
      printf("v %d: %s = %2.0f (max_dist = %2.2f) (%2.1f,%2.1f,%2.1f) --> (%2.1f,%2.1f,%2.1f)\n",
             vno, name, cv_parm[vno], max_dist, xs, ys, zs, xe, ye, ze) ;
      printf("    : max peak %2.0f (%d), smooth %2.0f (%d), weighted mean %2.1f\n",
             histo->bins[bpeak], bpeak,
             hsmooth->bins[bsmooth_peak], bsmooth_peak, parm) ;
      HISTOplot(histo, "histo.plt") ;
      HISTOplot(hsmooth, "hsmooth.plt") ;
    }
    if (!std::isfinite(parm))
      ErrorPrintf(ERROR_BADPARM, "sample_parameter_map(%s): vno %d, parm = %2.2f, total_wt = %2.2f\n",
                  name, vno, parm, total_wt) ;
  }

  HISTOfree(&histo) ;
  HISTOfree(&hsmooth) ;
  return(NO_ERROR) ;
}
#endif

#if 0
static double
compute_sse(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) {
  double sse ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;

  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  sse = MRIScomputeSSE(mris, parms) ;

  MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  sse += MRIScomputeSSE(mris, parms) ;

  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  return(sse) ;
}
#endif

typedef struct {
  double *flash ;    /* forward model f(T1,TR,alpha) */
  double TR ;
  double alpha ;
  double step_size ; /* T1 increment between table elements */
  double min_T1 ;
  int    size ;      /* number of elements in flash and norm */
}
FLASH_LOOKUP_TABLE, FLT ;

static FLT *find_lookup_table(double TR, double flip_angle) ;
static FLT lookup_tables[MAX_FLASH_VOLUMES] ;
static double *norms = NULL ;
static int ntables = 0 ;

static int
build_lookup_table(double tr, double flip_angle,
                   double min_T1, double max_T1, double step) {
  FLT    *flt ;
  int    i ;
  double T1 ;

  flt = find_lookup_table(tr, flip_angle) ;
  if (flt != NULL)
    return(NO_ERROR) ;  /* already created one */

  if (ntables >= MAX_FLASH_VOLUMES)
    ErrorExit(ERROR_NOMEMORY, "%s: MAX_FLASH_VOLUMES %d exceeded",
              Progname, MAX_FLASH_VOLUMES) ;


  flt = &lookup_tables[ntables++] ;
  flt->TR = tr ;
  flt->alpha = flip_angle ;
  flt->step_size = step ;

  flt->size = (int)((max_T1 - min_T1) / step)+1 ;
  flt->flash = (double *)calloc(flt->size, sizeof(*flt->flash)) ;
  if (!flt->flash)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocated %dth lookup table",
              Progname, ntables) ;

  for (i = 0, T1 = MIN_T1 ; T1 <= MAX_T1 ; i++, T1 += step)
    flt->flash[i] = FLASHforwardModel(flip_angle, tr, 1.0, T1) ;
  return(NO_ERROR) ;
}

#define T1_TO_INDEX(T1) (nint((T1 - MIN_T1) / T1_STEP_SIZE))

static double
lookup_flash_value(double TR, double flip_angle, double PD, double T1) {
  int   index ;
  FLT   *flt ;
  double FLASH ;

  flt = find_lookup_table(TR, flip_angle) ;
  if (!flt)
    return(FLASHforwardModel(flip_angle, TR, PD, T1)) ;

  index = T1_TO_INDEX(T1) ;
  if (index < 0)
    index = 0 ;
  if (index >= flt->size)
    index = flt->size-1 ;
  FLASH = PD*flt->flash[index] ;
  return(FLASH) ;
}

static FLT *
find_lookup_table(double TR, double flip_angle) {
  int   i ;

  for (i = 0 ; i < ntables ; i++)
    if (FEQUAL(lookup_tables[i].TR, TR) && FEQUAL(lookup_tables[i].alpha, flip_angle))
      break ;

  if (i >= ntables)
    return(NULL) ;
  return(&lookup_tables[i]) ;
}
static double
compute_optimal_vertex_positions(MRI_SURFACE *mris, int vno, EXTRA_PARMS *ep,
                                 double *pwhite_delta, double *ppial_delta,
                                 int debug_flag) {
  double       predicted_val, dx, dy, dz,
  dist, inward_dist, outward_dist, cortical_dist, PD_wm, PD_gm, PD_csf,
  T1_wm, T1_gm, T1_csf, white_dist, pial_dist, sse, best_sse,
  best_white_dist, best_pial_dist, orig_cortical_dist, total_dist,
  orig_white_index, orig_pial_index, sigma, pial_delta, wm_delta,
  gw_border, pial_border ;
  int          i, j, white_vno, pial_vno, max_j, best_white_index,
  best_pial_index, pial_index, white_index, max_white_index, nwm, npial  ;
  VERTEX       *v_white, *v_pial ;
  MRI          *mri ;
  double       image_vals[MAX_FLASH_VOLUMES][MAX_SAMPLES], xw, yw, zw, xp, yp, zp, x, y, z ;
  double       pwhite[MAX_FLASH_VOLUMES], pgray[MAX_FLASH_VOLUMES], pcsf[MAX_FLASH_VOLUMES],
  T1_vals[MAX_SAMPLES], PD_vals[MAX_SAMPLES] ;

  white_vno = vno ;
#if 0
  pial_vno = ep->nearest_pial_vertices[vno] ;
#else
  pial_vno = vno ;
#endif
  v_white = &mris->vertices[white_vno] ;
  v_pial = &mris->vertices[pial_vno] ;
#if 1
  inward_dist = ep->cv_inward_dists[white_vno] ;
  outward_dist = ep->cv_outward_dists[white_vno] ;
#else
  inward_dist = ep->max_dist ;
  outward_dist = ep->max_dist ;
#endif
  sigma = ep->current_sigma ;

  dx = v_pial->pialx - v_white->origx ;
  dy = v_pial->pialy - v_white->origy ;
  dz = v_pial->pialz - v_white->origz ;
  orig_cortical_dist = dist = sqrt(dx*dx + dy*dy + dz*dz) ;
  // MRIworldToVoxel(ep->mri_flash[0], v_white->origx, v_white->origy,v_white->origz,
  //                 &xw,&yw,&zw);
  // MRIworldToVoxel(ep->mri_flash[0], v_pial->pialx, v_pial->pialy,v_pial->pialz,
  //                 &xp,&yp,&zp);
  MRIsurfaceRASToVoxel(ep->mri_flash[0], v_white->origx, v_white->origy,v_white->origz,
                       &xw,&yw,&zw);
  MRIsurfaceRASToVoxel(ep->mri_flash[0], v_pial->pialx, v_pial->pialy,v_pial->pialz,
                       &xp,&yp,&zp);
  if (TOO_SMALL(dist)) {
    dx = v_white->nx ;
    dy = v_white->ny ;
    dz = v_white->nz ;
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    // MRIworldToVoxel(ep->mri_flash[0],
    //                 v_pial->pialx+v_pial->nx,
    //                 v_pial->pialy+v_pial->ny,
    //                 v_pial->pialz+v_pial->nz,
    //                 &xp,&yp,&zp);
    MRIsurfaceRASToVoxel(ep->mri_flash[0],
                         v_pial->pialx+v_pial->nx,
                         v_pial->pialy+v_pial->ny,
                         v_pial->pialz+v_pial->nz,
                         &xp,&yp,&zp);
  }

  if (FZERO(dist))
    dist = ep->dstep ;
  dx /= dist ;
  dy /= dist ;
  dz /= dist ;

  dx = xp - xw ;
  dy = yp - yw ;
  dz = zp - zw ;
  dist = sqrt(dx*dx + dy*dy + dz*dz) ;
  if (FZERO(dist))
    dist = ep->dstep ;
  dx /= dist ;
  dy /= dist ;
  dz /= dist ;  /* surface normal in volume coords */

  T1_wm = ep->cv_wm_T1[white_vno] ;
  PD_wm = ep->cv_wm_PD[white_vno] ;
  T1_gm = ep->cv_gm_T1[white_vno] ;
  PD_gm = ep->cv_gm_PD[white_vno] ;
  T1_csf = ep->cv_csf_T1[white_vno] ;
  PD_csf = ep->cv_csf_PD[white_vno] ;
  if (white_vno == Gdiag_no)
    printf("v %d: wm = (%2.0f, %2.0f), gm = (%2.0f, %2.0f), csf = (%2.0f, %2.0f)\n",
           Gdiag_no,ep->cv_wm_T1[Gdiag_no], ep->scale*ep->cv_wm_PD[Gdiag_no],
           ep->cv_gm_T1[Gdiag_no], ep->scale*ep->cv_gm_PD[Gdiag_no],
           ep->cv_csf_T1[Gdiag_no], ep->scale*ep->cv_csf_PD[Gdiag_no]);

  orig_white_index = inward_dist / ep->dstep ;
  orig_pial_index = orig_white_index + orig_cortical_dist / ep->dstep ;
  max_j = (int)((inward_dist + orig_cortical_dist + outward_dist) / ep->dstep) ;
  for (j = 0, dist = -inward_dist ; dist <= orig_cortical_dist+outward_dist ;
       dist += ep->dstep, j++) {
    x = xw + dist * dx ;
    y = yw + dist * dy ;
    z = zw + dist * dz ;
    MRIsampleVolumeType(ep->mri_T1, x, y, z, &T1_vals[j], sample_type) ;
    MRIsampleVolumeType(ep->mri_PD, x, y, z, &PD_vals[j], sample_type) ;
  }

  for (i = 0 ; i < ep->nvolumes ; i++) {
    mri = ep->mri_flash[i] ;
    for (j = 0, dist = -inward_dist ; dist <= orig_cortical_dist+outward_dist ;
         dist += ep->dstep, j++) {
      x = xw + dist * dx ;
      y = yw + dist * dy ;
      z = zw + dist * dz ;
#if 0
      MRIsampleVolumeDirectionScale(mri, x, y, z, dx, dy, dz, &image_vals[i][j], sigma) ;
#else
      MRIsampleVolumeType(mri, x, y, z, &image_vals[i][j], sample_type) ;
#endif
    }

    if ((plot_stuff == 1 && Gdiag_no == vno) ||
        plot_stuff > 1) {
      FILE *fp ;
      char fname[STRLEN] ;

      sprintf(fname, "vol%d.plt", i) ;
      fp = fopen(fname, "w") ;
      for (j = 0, dist = -inward_dist ; dist <= orig_cortical_dist+outward_dist ;
           dist += ep->dstep, j++)
        fprintf(fp, "%d %f %f\n", j, dist, image_vals[i][j]) ;
      fclose(fp) ;
    }
  }
  if ((plot_stuff == 1 && Gdiag_no == vno) ||
      plot_stuff > 1) {
    FILE *fp ;

    fp = fopen("orig_white.plt", "w") ;
    fprintf(fp, "%f %f %f\n%f %f %f\n",
            orig_white_index, 0.0, 0.0, orig_white_index, 0.0, 120.0) ;
    fclose(fp) ;
    fp = fopen("orig_pial.plt", "w") ;
    fprintf(fp, "%f %f %f\n%f %f %f\n",
            orig_pial_index, orig_cortical_dist, 0.0, orig_pial_index, orig_cortical_dist, 120.0) ;
    fclose(fp) ;
  }


  total_dist = outward_dist + orig_cortical_dist + inward_dist ;
  best_white_index = orig_white_index ;
  best_pial_index = orig_pial_index ;
  best_white_dist = orig_white_index*ep->dstep;
  best_pial_dist = orig_pial_index*ep->dstep;
#if 0
  max_white_index = nint((inward_dist + orig_cortical_dist/2)/ep->dstep) ;
#else
  max_white_index = nint(max_j - 1/ep->dstep) ;
#endif

  cortical_dist = best_pial_dist - best_white_dist ;
  best_sse =
    compute_vertex_sse(ep, image_vals, max_j, best_white_dist, cortical_dist,
                       T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf, 0,
                       T1_vals, PD_vals, vno) ;
  for (white_index = 0 ; white_index <= max_white_index ; white_index++) {
    white_dist = white_index * ep->dstep ;
    for (pial_index = white_index + nint(1.0/ep->dstep) ; pial_index <= max_j ;pial_index++) {
      pial_dist = pial_index * ep->dstep ;
      cortical_dist = pial_dist - white_dist ;
      sse = compute_vertex_sse(ep, image_vals, max_j, white_dist, cortical_dist,
                               T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf, 0,
                               T1_vals, PD_vals, vno) ;
      if (sse < best_sse) {
        best_white_index = white_index ;
        best_pial_index = pial_index ;
        best_pial_dist = pial_dist ;
        best_white_dist = white_dist ;
        best_sse = sse ;
      }
    }
  }

  if ((plot_stuff == 1 && Gdiag_no == vno) ||
      plot_stuff > 1) {
    FILE *fp ;
    char fname[STRLEN] ;

    cortical_dist = best_pial_dist - best_white_dist ;
    for (i = 0 ; i < ep->nvolumes ; i++) {
      sprintf(fname, "vol%d_best.plt", i) ;
      fp = fopen(fname, "w") ;
      mri = ep->mri_flash[i] ;
      for (j = 0 ; j <= max_j ; j++) {
        dist = (double)j*ep->dstep-best_white_dist ;
        predicted_val =
          image_forward_model(mri->tr, mri->flip_angle, dist, cortical_dist,
                              T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf) ;
        fprintf(fp, "%d %f %f\n", j, dist, predicted_val) ;
      }
      fclose(fp) ;
    }
  }

  if ((plot_stuff == 1 && Gdiag_no == vno) ||
      plot_stuff > 1) {
    FILE *fp ;

    fp = fopen("best_white.plt", "w") ;
    fprintf(fp, "%d %f %f\n%d %f %f\n",
            best_white_index, 0.0, 0.0, best_white_index, 0.0, 120.0) ;
    fclose(fp) ;
    fp = fopen("best_pial.plt", "w") ;
    fprintf(fp, "%d %f %f\n%d %f %f\n",
            best_pial_index, cortical_dist, 0.0, best_pial_index, cortical_dist, 120.0) ;
    fclose(fp) ;
  }



  for (nwm = npial = 0, pial_delta = wm_delta = 0.0, i = 0 ; i < ep->nvolumes ; i++) {
    double I0, I1, Idist ;

    mri = ep->mri_flash[i] ;
    pwhite[i] =
      image_forward_model(mri->tr, mri->flip_angle, -1, cortical_dist,
                          T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf) ;
    pgray[i] =
      image_forward_model(mri->tr, mri->flip_angle, cortical_dist/2, cortical_dist,
                          T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf) ;
    pcsf[i] =
      image_forward_model(mri->tr, mri->flip_angle, cortical_dist+1, cortical_dist,
                          T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf) ;
    gw_border = (pwhite[i] + pgray[i])/2 ;
    pial_border = (pcsf[i] + pgray[i])/2 ;

    I0 = image_vals[i][best_white_index] ;
    if (best_white_index > 0) {
      I1 = image_vals[i][best_white_index-1] ;
      Idist = I1-I0 ;
      if (fabs(gw_border-I0) < fabs(Idist) &&
          fabs(gw_border-I1) < fabs(Idist)) {    /* predicted value is between best and previous */
        wm_delta -= (gw_border - I0) / Idist ;
        nwm++ ;
      }
    }

    if (best_white_index < max_j) {
      I1 = image_vals[i][best_white_index+1] ;
      Idist = I1-I0 ;

      if (fabs(gw_border-I0) < fabs(Idist) &&
          fabs(gw_border-I1) < fabs(Idist)) {    /* predicted value is between best and previous */
        wm_delta += (gw_border - I0)  / Idist ;
        nwm++ ;
      }
    }

    I0 = image_vals[i][best_pial_index] ;
    if (best_pial_index > 0) {
      I1 = image_vals[i][best_pial_index-1] ;
      Idist = I1-I0 ;
      if (fabs(pial_border-I0) < fabs(Idist) &&
          fabs(pial_border-I1) < fabs(Idist)) {    /* predicted value is between best and previous */
        pial_delta -= (pial_border - I0) / Idist ;
        npial++ ;
      }
    }

    if (best_pial_index < max_j) {
      I1 = image_vals[i][best_pial_index+1] ;
      Idist = I1-I0 ;

      if (fabs(pial_border-I0) < fabs(Idist) &&
          fabs(pial_border-I1) < fabs(Idist)) {    /* predicted value is between best and previous */
        pial_delta += (pial_border - I0) / Idist ;
        npial++ ;
      }
    }
  }

  if (nwm)
    wm_delta /= (float)nwm ;
  if (npial)
    pial_delta /= (float)npial ;


  *pwhite_delta = (best_white_dist - (orig_white_index * ep->dstep)) ;
  *ppial_delta =  (best_pial_dist  - (orig_pial_index  * ep->dstep)) ;

  if ((abs(vno-Gdiag_no) < 100) && (ep->scale*image_vals[3][nint(orig_pial_index)] < 100) && (*ppial_delta > 0.5)) {
    FILE *fp ;
    fp = fopen("v.log","a") ;
    fprintf(fp, "v %d, T1=(%2.0f,%2.0f,%2.0f), PD=(%2.0f,%2.0f,%2.0f), D=%2.1f,%2.1f\n",
            vno, T1_wm, T1_gm, T1_csf, PD_wm, PD_gm, PD_csf,
            *pwhite_delta, *ppial_delta) ;
    fclose(fp) ;
    DiagBreak() ;
  }

  if (vno == Gdiag_no && debug_flag) {
    double pwhite[MAX_FLASH_VOLUMES], pgray[MAX_FLASH_VOLUMES], pcsf[MAX_FLASH_VOLUMES];

    printf("current white intensities: ") ;
    for (i = 0 ; i < ep->nvolumes ; i++)
      printf("%2.1f  ", ep->scale*image_vals[i][nint(orig_white_index)]) ;
#if 0
    printf("\nbest    white intensities: ") ;
    for (i = 0 ; i < ep->nvolumes ; i++)
      printf("%2.1f  ", ep->scale*image_vals[i][best_white_index]) ;
#endif

    printf("\ncurrent pial  intensities: ") ;
    for (i = 0 ; i < ep->nvolumes ; i++)
      printf("%2.1f  ", ep->scale*image_vals[i][nint(orig_pial_index)]) ;
#if 0
    printf("\nbest    pial  intensities: ") ;
    for (i = 0 ; i < ep->nvolumes ; i++)
      printf("%2.1f  ", ep->scale*image_vals[i][best_pial_index]) ;
#endif
    printf("\n") ;


    for (i = 0 ; i < ep->nvolumes ; i++) {
      mri = ep->mri_flash[i] ;
      pwhite[i] =
        image_forward_model(mri->tr, mri->flip_angle, -1, cortical_dist,
                            T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf) ;
      pgray[i] =
        image_forward_model(mri->tr, mri->flip_angle, cortical_dist/2, cortical_dist,
                            T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf) ;
      pcsf[i] =
        image_forward_model(mri->tr, mri->flip_angle, cortical_dist+1, cortical_dist,
                            T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf) ;
      printf("\tpredicted image intensities %d: (%2.0f, %2.0f, %2.0f)\n",
             i, ep->scale*pwhite[i], ep->scale*pgray[i], ep->scale*pcsf[i]) ;
    }

    printf("v %d: best_white_delta = %2.2f (%2.2f), best_pial_delta = %2.2f (%2.2f), best_sse = %2.2f\n",
           vno, *pwhite_delta,
           (best_white_dist  - (orig_white_index  * ep->dstep)),
           *ppial_delta,
           (best_pial_dist  - (orig_pial_index  * ep->dstep)), best_sse) ;
    printf("      inward_dist = %2.2f, outward_dist = %2.2f, cortical_dist = %2.2f\n",
           inward_dist, outward_dist, orig_cortical_dist) ;
    cortical_dist = best_pial_dist - best_white_dist ;
    sse = compute_vertex_sse(ep, image_vals, max_j, best_white_dist, cortical_dist,
                             T1_wm, PD_wm, T1_gm, PD_gm, T1_csf, PD_csf, plot_stuff,
                             T1_vals, PD_vals, vno) ;

  }

  return(best_sse) ;
}

#if 0
static int
compute_T1_PD_slow(double *image_vals, MRI **mri_flash, int nvolumes,
                   double *pT1, double *pPD) {
  double    sse, best_sse, best_T1, best_PD, pred, PD, T1, error ;
  int       i ;
  MRI       *mri ;

  if (plot_stuff <= 2)   /* disable it for now */
    return(NO_ERROR) ;


  best_T1 = best_PD = 0 ;
  best_sse = -1 ;
  for (PD = 100.0 ; PD < 1700 ; PD += T1_STEP_SIZE) {
    for (T1 = MIN_T1 ; T1 <= MAX_T1 ; T1 += T1_STEP_SIZE) {
      for (sse = 0.0, i = 0 ; i < nvolumes ; i++) {
        mri = mri_flash[i] ;
        pred = FLASHforwardModel(mri->flip_angle, mri->tr, PD, T1) ;
        error = (pred - image_vals[i]) ;
        sse += (error*error) ;
      }
      if (sse < best_sse || best_sse < 0) {
        best_sse = sse ;
        best_T1 = T1 ;
        best_PD = PD ;
      }
    }
  }

  *pT1 = best_T1 ;
  *pPD = best_PD ;
  return(NO_ERROR) ;
}
#endif
static int
compute_T1_PD(double *image_vals, MRI **mri_flash, int nvolumes,
              double *pT1, double *pPD) {
  double    best_T1, best_PD, norm_im, norm_pred, sse, best_sse, T1,
  pred_vals[MAX_FLASH_VOLUMES], error, upper_T1, lower_T1, mid_T1,
  upper_sse, lower_sse, mid_sse, upper_norm, mid_norm, lower_norm, range ;
  int       i, j, upper_j, lower_j, mid_j, niter ;
  MRI       *mri ;

  if (!norms) {
    norms = (double *)calloc(nint((MAX_T1-MIN_T1)/T1_STEP_SIZE)+1, sizeof(double)) ;
    if (!norms)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate norm table", Progname) ;
    for (j = 0, T1 = MIN_T1 ; T1 < MAX_T1 ; T1 += T1_STEP_SIZE, j++) {
      for (norm_pred = 0.0, i = 0 ; i < nvolumes ; i++) {
        mri = mri_flash[i] ;
        pred_vals[i] = FLASHforwardModelLookup(mri->flip_angle, mri->tr, 1.0, T1) ;
        norm_pred += (pred_vals[i]*pred_vals[i]) ;
      }
      norms[j] = sqrt(norm_pred) ;
      if (FZERO(norms[j])) {
        printf("norms[%d] is zero!\n", j) ;
        DiagBreak() ;
        exit(0) ;
      }
    }
  }

  for (norm_im = i = 0 ; i < nvolumes ; i++)
    norm_im += (image_vals[i]*image_vals[i]) ;
  norm_im = sqrt(norm_im) ;   /* length of image vector */
  if (FZERO(norm_im)) {
    *pT1 = MIN_T1 ;
    *pPD = MIN_T1 ;
    return(ERROR_BADPARM) ;
  }
  for (i = 0 ; i < nvolumes ; i++)
    image_vals[i] /= norm_im ;   /* normalize them */

  mid_T1 = (MAX_T1 - MIN_T1)/2 ;
  mid_j = T1_TO_INDEX(mid_T1) ;
  range = (MAX_T1 - MIN_T1) / 2 ;

  /* compute sse for mid T1 */
  mid_norm = norms[mid_j] ;
  for (mid_sse = 0.0, i = 0 ; i < nvolumes ; i++) {
    mri = mri_flash[i] ;
    pred_vals[i] = FLASHforwardModelLookup(mri->flip_angle, mri->tr, 1.0, mid_T1) ;
    pred_vals[i] /= mid_norm ; /* normalize them */
    error = (pred_vals[i] - image_vals[i]) ;
    mid_sse += (error * error) ;
  }

  best_T1 = mid_T1 ;
  best_PD = norm_im / mid_norm ;
  best_sse = mid_sse ;
  niter = 0 ;
  if (FZERO(mid_norm)) {
    printf("mid norm=0 at %d (%2.1f)\n", mid_j, mid_T1) ;
    DiagBreak() ;
    exit(0) ;
  }
  do {
    upper_T1 = mid_T1 + 0.5*range ;
    lower_T1 = mid_T1 - 0.5*range ;
    if (upper_T1 > MAX_T1)
      upper_T1 = MAX_T1 ;
    if (lower_T1 < MIN_T1)
      lower_T1 = MIN_T1 ;
    upper_j = T1_TO_INDEX(upper_T1) ;
    lower_j = T1_TO_INDEX(lower_T1) ;
    upper_norm = norms[upper_j] ;
    lower_norm = norms[lower_j] ;
    if (FZERO(upper_norm)) {
      printf("upper norm=0 at %d (%2.1f)\n", upper_j, upper_T1) ;
      DiagBreak() ;
      exit(0) ;
    }
    if (FZERO(lower_norm)) {
      printf("lower norm=0 at %d (%2.1f)\n", lower_j, lower_T1) ;
      DiagBreak() ;
      exit(0) ;
    }
    for (lower_sse = upper_sse = 0.0, i = 0 ; i < nvolumes ; i++) {
      mri = mri_flash[i] ;

      pred_vals[i] = FLASHforwardModelLookup(mri->flip_angle, mri->tr, 1.0, upper_T1) ;
      pred_vals[i] /= upper_norm ; /* normalize them */
      error = (pred_vals[i] - image_vals[i]) ;
      upper_sse += (error * error) ;

      pred_vals[i] = FLASHforwardModelLookup(mri->flip_angle, mri->tr, 1.0, lower_T1) ;
      pred_vals[i] /= lower_norm ; /* normalize them */
      error = (pred_vals[i] - image_vals[i]) ;
      lower_sse += (error * error) ;
    }

    if (lower_sse <= mid_sse && lower_sse <= upper_sse) /* make lower new mid */
    {
      mid_sse = lower_sse ;
      mid_j = lower_j ;
      mid_norm = lower_norm ;
      mid_T1 = lower_T1 ;
      best_T1 = lower_T1 ;
      best_PD = norm_im / lower_norm ;
      best_sse = lower_sse ;
    } else if (upper_sse < mid_sse)  /* make upper new mid */
    {
      mid_sse = upper_sse ;
      mid_j = upper_j ;
      mid_norm = upper_norm ;
      mid_T1 = upper_T1 ;
      best_T1 = upper_T1 ;
      best_PD = norm_im / upper_norm ;
      best_sse = upper_sse ;
    }
    if (!std::isfinite(best_PD)) {
      printf("best_PD is not finite at %d (%2.1f)\n", mid_j, mid_T1) ;
      DiagBreak() ;
      exit(0) ;
    }
    range /= 2 ;
    niter++ ;
  } while (upper_j - lower_j > 3) ;

  for (i = 0 ; i < nvolumes ; i++)
    image_vals[i] *= norm_im ;   /* restore them */

  for (sse = 0.0, i = 0 ; i < nvolumes ; i++) {
    mri = mri_flash[i] ;
    pred_vals[i] = FLASHforwardModelLookup(mri->flip_angle, mri->tr, best_PD, best_T1) ;
    error = (pred_vals[i] - image_vals[i]) ;
    sse += (error * error) ;
  }

  *pT1 = best_T1 ;
  *pPD = best_PD ;
  return(NO_ERROR) ;
}


static double
scale_all_images(MRI **mri_flash, int nvolumes, MRI_SURFACE *mris, float target_pd_wm,
                 EXTRA_PARMS *ep) {
  int    vno, i ;
  VERTEX *v ;
  double xw, yw, zw, T1, PD, T1_wm_total, PD_wm_total, T1_gm_total, PD_gm_total ;
  double mean_wm[MAX_FLASH_VOLUMES], scale ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;  /* recompute white surface normals */
  MRISsmoothSurfaceNormals(mris, 5) ;
  T1_wm_total = PD_wm_total = 0 ;
  T1_gm_total = PD_gm_total = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    for (i = 0 ; i < nvolumes ; i++) {
      // MRIworldToVoxel(mri_flash[i], v->x-v->nx, v->y-v->ny, v->z-v->nz, &xw, &yw, &zw) ;
      // converting surface vertex to volume voxel
      MRIsurfaceRASToVoxel(mri_flash[i],v->x-v->nx, v->y-v->ny, v->z-v->nz, &xw, &yw, &zw);
      MRIsampleVolumeType(mri_flash[i], xw, yw, zw, &mean_wm[i], sample_type) ;
    }
    compute_T1_PD(mean_wm, mri_flash, nvolumes, &T1, &PD) ;
    T1_wm_total += T1 ;
    PD_wm_total += PD ;
    if (!std::isfinite(T1) || !std::isfinite(PD) ||
        !std::isfinite(T1_wm_total) || !std::isfinite(PD_wm_total))
      DiagBreak() ;
    ep->cv_wm_T1[vno] = T1 ;
    ep->cv_wm_PD[vno] = PD ;

    for (i = 0 ; i < nvolumes ; i++) {
      // MRIworldToVoxel(mri_flash[i], v->x+v->nx, v->y+v->ny, v->z+v->nz, &xw, &yw, &zw) ;
      // converting surface vertex to volume voxel
      MRIsurfaceRASToVoxel(mri_flash[i],v->x+v->nx, v->y+v->ny, v->z+v->nz, &xw, &yw, &zw);
      MRIsampleVolumeType(mri_flash[i], xw, yw, zw, &mean_wm[i], sample_type) ;
    }
    compute_T1_PD(mean_wm, mri_flash, nvolumes, &T1, &PD) ;
    T1_gm_total += T1 ;
    PD_gm_total += PD ;
    ep->cv_gm_T1[vno] = T1 ;
    ep->cv_gm_PD[vno] = PD ;
  }
  T1_wm_total /= (float)mris->nvertices ;
  PD_wm_total /= (float)mris->nvertices ;
  T1_gm_total /= (float)mris->nvertices ;
  PD_gm_total /= (float)mris->nvertices ;
  scale = target_pd_wm / PD_wm_total ;
  printf("mean wm T1 = %2.0f, mean wm PD = %2.0f, scaling images by %2.3f\n",
         T1_wm_total, PD_wm_total, scale) ;
  printf("mean gm T1 = %2.0f, mean gm PD = %2.0f\n",
         T1_gm_total, PD_gm_total) ;
  for (i = 0 ; i < nvolumes ; i++)
    MRIscalarMul(mri_flash[i], mri_flash[i], scale) ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  return(1.0/scale) ;
}

static int
store_current_positions(MRI_SURFACE *mris, EXTRA_PARMS *parms) {
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    parms->cv_last_pialx[vno] = v->pialx ;
    parms->cv_last_pialy[vno] = v->pialy ;
    parms->cv_last_pialz[vno] = v->pialz ;

    parms->cv_last_whitex[vno] = v->origx ;
    parms->cv_last_whitey[vno] = v->origy ;
    parms->cv_last_whitez[vno] = v->origz ;
  }
  return(NO_ERROR) ;
}

static int
update_parameters(MRI_SURFACE *mris, EXTRA_PARMS *ep) {
  int    vno, nmarked = 0 ;
  VERTEX *v ;
  double  best_white_delta, best_pial_delta ;

  MRISclearMarks(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->d < 1)
      continue ;
    nmarked++ ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    compute_optimal_parameters(mris, vno, ep, &best_white_delta,
                               &best_pial_delta) ;
    if (vno == Gdiag_no) {
      printf("v %d: wm = (%2.0f, %2.0f), gm = (%2.0f, %2.0f), "
             "csf = (%2.0f, %2.0f), "
             "deltas = (%2.2f, %2.2f)\n",
             vno,ep->cv_wm_T1[vno], ep->scale*ep->cv_wm_PD[vno], ep->cv_gm_T1[vno],
             ep->scale*ep->cv_gm_PD[vno], ep->cv_csf_T1[vno],
             ep->scale*ep->cv_csf_PD[vno],
             best_white_delta, best_pial_delta);
      DiagBreak() ;
    }
    v->d = 0 ;  /* reset distance moved to 0 */
    v->marked = 1 ;
  }
  smooth_marked_maps(mris, ep, smooth_parms) ;
  constrain_parameters(mris->nvertices, ep) ;
  if (Gdiag_no > 0 && mris->vertices[Gdiag_no].marked > 0) {
    v = &mris->vertices[vno = Gdiag_no] ;

    printf("v %d: wm = (%2.0f, %2.0f), gm = (%2.0f, %2.0f), "
           "csf = (%2.0f, %2.0f)\n",
           vno,ep->cv_wm_T1[vno], ep->scale*ep->cv_wm_PD[vno],
           ep->cv_gm_T1[vno],
           ep->scale*ep->cv_gm_PD[vno], ep->cv_csf_T1[vno],
           ep->scale*ep->cv_csf_PD[vno]) ;
    DiagBreak() ;
  }
  MRISclearMarks(mris) ;
  printf("Tissue parameters recomputed for %d vertices\n", nmarked) ;

  return(NO_ERROR) ;
}
static int
update_distances(MRI_SURFACE *mris, EXTRA_PARMS *parms) {
  int    vno ;
  VERTEX *v ;
  float  dot, nx, ny, nz, norm, dx, dy, dz ;

  for (vno = MIN_VNO ; vno < mris->nvertices ; vno++) {
    if (vno >= MAX_VNO)
      break ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    nx = parms->cv_last_pialx[vno] - parms->cv_last_whitex[vno] ;
    ny = parms->cv_last_pialy[vno] - parms->cv_last_whitey[vno] ;
    nz = parms->cv_last_pialz[vno] - parms->cv_last_whitez[vno] ;
    norm = sqrt(nx*nx + ny*ny + nz*nz) ;
    if (TOO_SMALL(norm)) {
      nx = v->nx ;
      ny = v->ny ;
      nz = v->nz ;
      norm = sqrt(nx*nx + ny*ny + nz*nz) ;
    }
    if (FZERO(norm))
      norm = 1 ;
    nx /= norm ;
    ny /= norm ;
    nz /= norm ;

    dx = v->pialx - parms->cv_last_pialx[vno] ;
    dy = v->pialy - parms->cv_last_pialy[vno] ;
    dz = v->pialz - parms->cv_last_pialz[vno] ;
    dot = dx*nx + dy*ny + dz*nz ;
    dot *= -1 ;
    v->d += fabs(dot) ;
    parms->cv_outward_dists[vno] += dot ;
    if (vno == Gdiag_no)
      printf("v %d: updating outward dist by %2.1f to %2.1f\n",
             vno, dot, parms->cv_outward_dists[vno]) ;

    dx = v->origx - parms->cv_last_whitex[vno] ;
    dy = v->origy - parms->cv_last_whitey[vno] ;
    dz = v->origz - parms->cv_last_whitez[vno] ;
    dot = dx*nx + dy*ny + dz*nz ;
    parms->cv_inward_dists[vno] += dot ;
    v->d += fabs(dot) ;
    if (vno == Gdiag_no)
      printf("v %d: updating inward dist by %2.1f to %2.1f (d=%2.2f)\n", vno, dot,
             parms->cv_inward_dists[vno],v->d) ;
  }
  return(NO_ERROR) ;
}
static int
constrain_parameters(int nvertices, EXTRA_PARMS *ep) {
  int      i ;
  double   T1_wm, T1_gm, T1_csf, PD_csf, delta, PD_gm ;

  for (i = MIN_VNO ; i < nvertices ; i++) {
    if (i == Gdiag_no)
      DiagBreak() ;
    if (i >= MAX_VNO)
      break ;
    T1_wm = ep->cv_wm_T1[i] ;
    T1_gm = ep->cv_gm_T1[i] ;
    T1_csf = ep->cv_csf_T1[i] ;
    PD_csf = ep->cv_csf_PD[i];
    PD_gm = ep->cv_gm_PD[i] ;

    if (PD_csf > MIN_NONBRAIN_PD) {
#if 0
      if (PD_csf < MIN_CSF_PD)
        PD_csf = MIN_CSF_PD ;
#endif
      if (T1_csf < MIN_CSF_T1)
        T1_csf = MIN_CSF_T1 ;
      if (PD_gm > PD_csf*.9) {
        delta = (PD_csf - PD_csf*.9) - (PD_csf - PD_gm) ;
        PD_csf += delta/2 ;
        PD_gm -= delta/2 ;
      }
    }
    if ((PD_csf > MIN_NONBRAIN_PD) && (T1_gm > T1_csf*.9)) {
      delta = (T1_csf - T1_csf*.9) - (T1_csf - T1_gm) ;
      T1_csf += delta/2 ;
      T1_gm -= delta/2 ;
    }
    if (T1_wm > T1_gm*.9) {
      delta = (T1_gm - T1_gm*.9) - (T1_gm - T1_wm) ;
      T1_gm += delta/2 ;
      T1_wm -= delta/2 ;
    }

    ep->cv_wm_T1[i] = T1_wm ;
    ep->cv_gm_T1[i] = T1_gm ;
    ep->cv_gm_PD[i] = PD_gm ;
    ep->cv_csf_T1[i] = T1_csf ;
    ep->cv_csf_PD[i] = PD_csf ;
  }
  return(NO_ERROR) ;
}

static int
compute_parameter_maps(MRI **mri_flash, int nvolumes, MRI **pmri_T1,
                       MRI **pmri_PD) {
  int   i, x, y, z, width, height, depth ;
  double image_vals[MAX_FLASH_VOLUMES] ;
  MRI   *mri_T1, *mri_PD ;
  double T1, PD ;

  width = mri_flash[0]->width ;
  height = mri_flash[0]->height ;
  depth = mri_flash[0]->depth ;
  mri_T1 = MRIclone(mri_flash[0], NULL) ;
  mri_PD = MRIclone(mri_flash[0], NULL) ;

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak()  ;
        for (i = 0 ; i < nvolumes ; i++)
          MRIsampleVolume(mri_flash[i], x, y, z, &image_vals[i]) ;
        compute_T1_PD(image_vals, mri_flash, nvolumes, &T1, &PD) ;
        MRIsetVoxVal(mri_T1, x, y,z,0,T1);
        MRIsetVoxVal(mri_PD, x, y,z,0,PD);
      }
    }
  }
  *pmri_T1 = mri_T1 ;
  *pmri_PD = mri_PD ;
  return(NO_ERROR) ;
}

static int
ms_errfunc_rip_vertices(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) {
  int         vno, ripped = 0, n ;
  EXTRA_PARMS *ep ;
  double      white_delta, pial_delta ;

  ep = (EXTRA_PARMS *)parms->user_parms ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX * const v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    compute_optimal_vertex_positions(mris, vno, ep,
                                     &white_delta, &pial_delta,0) ;
    if ((fabs(white_delta) < 0.01) && (fabs(pial_delta) < 0.01)) {
      v->ripflag = 1 ;
      ripped++ ;
    }
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    for (n = 0 ; n < vt->vnum ; n++) {
      VERTEX * const vn = &mris->vertices[vt->v[n]] ;
      if (vn->ripflag) {
        vn->ripflag = 0 ;   /* allow neighbors of unripped vertices to move */
        ripped-- ;
      }
    }
  }
  printf("%d vertices ripped (%2.3f%% of surface)\n",
         ripped, 100.0f*(float)ripped/mris->nvertices) ;
  if (Gdiag_no > 0)
    printf("v %d: ripflag = %d\n", Gdiag_no, mris->vertices[Gdiag_no].ripflag);
  return(NO_ERROR) ;
}

#if 0
static int
compute_PD_T1_limits(MRI *mri_flash[], int nvolumes, MRI_SURFACE *mris, EXTRA_PARMS *ep, int navgs) {
  int    vno ;

  MRISimportCurvatureVector(mris, ep->cv_wm_PD) ;
  MRISaverageCurvatures(mris, 25) ;   /* reduce effects of outliers */
  MRISextractCurvatureVector(mris, ep->cv_wm_PD) ;
  MRISminFilterCurvatures(mris, navgs) ;
  MRISextractCurvatureVector(mris, ep->cv_min_wm_PD) ;

  MRISimportCurvatureVector(mris, ep->cv_wm_PD) ;
  MRISmaxFilterCurvatures(mris, navgs) ;
  MRISextractCurvatureVector(mris, ep->cv_max_wm_PD) ;

  MRISimportCurvatureVector(mris, ep->cv_gm_PD) ;
  MRISaverageCurvatures(mris, 25) ;   /* reduce effects of outliers */
  MRISextractCurvatureVector(mris, ep->cv_gm_PD) ;
  MRISminFilterCurvatures(mris, navgs) ;
  MRISextractCurvatureVector(mris, ep->cv_min_gm_PD) ;

  MRISimportCurvatureVector(mris, ep->cv_gm_PD) ;
  MRISmaxFilterCurvatures(mris, navgs) ;
  MRISextractCurvatureVector(mris, ep->cv_max_gm_PD) ;

  if (Gdiag_no > 0) {
    printf("%d: WM PD limits %2.0f: (%2.0f --> %2.0f), GM PD limits %2.0f: (%2.0f --> %2.0f)\n",
           Gdiag_no, ep->cv_wm_PD[Gdiag_no], ep->cv_min_wm_PD[Gdiag_no], ep->cv_max_wm_PD[Gdiag_no],
           ep->cv_gm_PD[Gdiag_no], ep->cv_min_gm_PD[Gdiag_no], ep->cv_max_gm_PD[Gdiag_no]) ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    ep->cv_min_gm_PD[vno] = ep->cv_gm_PD[vno]-2*PD_STD ;
    ep->cv_max_gm_PD[vno] = ep->cv_gm_PD[vno]+2*PD_STD ;

    ep->cv_min_wm_PD[vno] = ep->cv_wm_PD[vno]-2*PD_STD ;  /* do something smarter than this */
    ep->cv_max_wm_PD[vno] = ep->cv_wm_PD[vno]+2*PD_STD ;
  }
  return(NO_ERROR) ;
}
#else
static int
compute_PD_T1_limits(MRI_SURFACE *mris, EXTRA_PARMS *ep, int navgs) {
  int    vno, found ;
  VERTEX *v, *v_white,  *v_pial ;
  double x, y, z, T1, PD, PD_min, PD_max, T1_min, T1_max,  n, cortical_dist,
  xp, yp, zp, xw, yw, zw, dx, dy,  dz ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;  /* recompute white surface normals */
  MRISsmoothSurfaceNormals(mris, 5) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    /* sample inwards to find  wm limits */
    T1_min = MIN_WM_T1  ;
    T1_max = MAX_WM_T1 ;
    PD_min = 100000 ;
    PD_max = 0 ;
    for (found = n = 0 ; n < ep->max_inward_dist ; n++) {
      // MRIworldToVoxel(ep->mri_PD, v->x-n*v->nx, v->y-n*v->ny, v->z-n*v->nz, &x, &y, &z) ;
      MRIsurfaceRASToVoxel(ep->mri_PD, v->x-n*v->nx, v->y-n*v->ny, v->z-n*v->nz, &x, &y, &z) ;
      MRIsampleVolumeType(ep->mri_PD, x, y, z, &PD, sample_type) ;
      MRIsampleVolumeType(ep->mri_T1, x, y, z, &T1, sample_type) ;
      if (T1 < MIN_WM_T1 || T1 > MAX_WM_T1)
        continue ;
      if (PD < PD_min)
        PD_min = PD ;
      if (PD > PD_max)
        PD_max = PD ;
      if (T1 < T1_min)
        T1_min = T1 ;
      if (T1 > T1_max)
        T1_max = T1 ;
      found = 1 ;
    }
    if (!found) {
      PD_min = MIN_WM_PD ;
      PD_max = MAX_WM_PD ;
      T1_min = MIN_WM_T1 ;
      T1_max = MAX_WM_T1 ;
    } else {
      T1_min -= WM_T1_STD ;
      T1_max += WM_T1_STD ;
      PD_min -= PD_STD ;
      PD_max += PD_STD ;
    }
    ep->cv_min_wm_PD[vno] = PD_min ;
    ep->cv_max_wm_PD[vno] = PD_max ;
    ep->cv_min_wm_T1[vno] = T1_min ;
    ep->cv_max_wm_T1[vno] = T1_max ;

    T1_min = MIN_GM_T1 ;
    T1_max  = MAX_GM_T1  ;
    PD_min = 100000 ;
    PD_max = 0 ;
    /* try not to go into non-brain regions. If we don't get out of the wm, it doesn't
    really matter, as it has a lower PD than gray in any case, so will be a reasonable lower
    bound. Note: gm_max_PD isn't currently used */
    v_white = &mris->vertices[vno] ;
    v_pial = v_white ;
    MRIsurfaceRASToVoxel(ep->mri_T1, v_white->origx, v_white->origy, v_white->origz, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(ep->mri_T1, v_pial->pialx, v_pial->pialy, v_pial->pialz, &xp, &yp, &zp) ;
    dx = xp-xw ;
    dy = yp-yw ;
    dz = zp-zw ;
    cortical_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (cortical_dist  < 1)
      cortical_dist = 1 ;
    for (found = 0, n = 1 ; n <= cortical_dist ; n++) {
      // MRIworldToVoxel(ep->mri_PD, v->x+n*v->nx, v->y+n*v->ny, v->z+n*v->nz, &x, &y, &z) ;
      MRIsurfaceRASToVoxel(ep->mri_PD, v->x+n*v->nx, v->y+n*v->ny, v->z+n*v->nz, &x, &y, &z) ;
      MRIsampleVolumeType(ep->mri_PD, x, y, z, &PD, sample_type) ;
      MRIsampleVolumeType(ep->mri_T1, x, y, z, &T1, sample_type) ;
#if 1
      if (T1 < MIN_GM_T1  || PD < MIN_GM_PD || T1 > MIN_CSF_T1)
        continue ;
#endif
      if (PD < PD_min)
        PD_min = PD ;
      if (PD > PD_max)
        PD_max = PD ;
      if (T1-2*GM_T1_STD < T1_min)
        T1_min = T1-2*GM_T1_STD ;
      if (T1+2*GM_T1_STD > T1_max)
        T1_max = T1+2*GM_T1_STD ;
      found = 1 ;
    }
    if (found) {
      /*   T1_min -= WM_T1_STD ; T1_max += WM_T1_STD ;*/
      PD_min -= PD_STD ;
      PD_max += PD_STD ;
    } else {
      PD_min = MIN_GM_PD  ;
      PD_max = MAX_GM_PD;
    }

    ep->cv_min_gm_T1[vno] = T1_min ;
    ep->cv_max_gm_T1[vno] = T1_max ;
    ep->cv_min_gm_PD[vno] = PD_min ;
    ep->cv_max_gm_PD[vno] = PD_max ;
  }
#if 0
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    for  (i = 0 ; i < v->vtotal ; i++) {
      vn = &mris->vertices[v->v[i]] ;
      if  (vn->ripflag)
        continue  ;
      if (ep->cv_min_wm_PD[v->v[i]]  < ep->cv_min_wm_PD[vno])
        ep->cv_min_wm_PD[vno] = ep->cv_min_wm_PD[v->v[i]]  ;
      if (ep->cv_min_wm_T1[v->v[i]]  < ep->cv_min_wm_T1[vno])
        ep->cv_min_wm_T1[vno] = ep->cv_min_wm_T1[v->v[i]]  ;
      if (ep->cv_min_gm_PD[v->v[i]]  < ep->cv_min_gm_PD[vno])
        ep->cv_min_gm_PD[vno] = ep->cv_min_gm_PD[v->v[i]]  ;
      if (ep->cv_min_gm_T1[v->v[i]]  < ep->cv_min_gm_T1[vno])
        ep->cv_min_gm_T1[vno] = ep->cv_min_gm_T1[v->v[i]]  ;

      if (ep->cv_max_wm_PD[v->v[i]]  > ep->cv_max_wm_PD[vno])
        ep->cv_max_wm_PD[vno] = ep->cv_max_wm_PD[v->v[i]]  ;
      if (ep->cv_max_wm_T1[v->v[i]]  > ep->cv_max_wm_T1[vno])
        ep->cv_max_wm_T1[vno] = ep->cv_max_wm_T1[v->v[i]]  ;
      if (ep->cv_max_gm_PD[v->v[i]]  > ep->cv_max_gm_PD[vno])
        ep->cv_max_gm_PD[vno] = ep->cv_max_gm_PD[v->v[i]]  ;
      if (ep->cv_max_gm_T1[v->v[i]]  > ep->cv_max_gm_T1[vno])
        ep->cv_max_gm_T1[vno] = ep->cv_max_gm_T1[v->v[i]]  ;
    }

    if (vno == Gdiag_no) {
      printf("wm PD: (%2.1f --> %2.1f), T1: (%2.1f --> %2.1f)\n", ep->cv_min_wm_PD[vno], ep->cv_max_wm_PD[vno],
             ep->cv_min_wm_T1[vno], ep->cv_max_wm_T1[vno])  ;
      printf("gm PD: (%2.1f --> %2.1f), T1: (%2.1f --> %2.1f)\n", ep->cv_min_gm_PD[vno], ep->cv_max_gm_PD[vno],
             ep->cv_min_gm_T1[vno], ep->cv_max_gm_T1[vno])  ;
    }
  }
#endif
  smooth_map(mris, ep->cv_min_wm_PD, navgs) ;
  smooth_map(mris, ep->cv_max_wm_PD, navgs) ;
  smooth_map(mris, ep->cv_min_gm_PD, navgs) ;
  smooth_map(mris, ep->cv_max_gm_PD, navgs) ;

  smooth_map(mris, ep->cv_min_wm_T1, navgs) ;
  smooth_map(mris, ep->cv_max_wm_T1, navgs) ;
  smooth_map(mris, ep->cv_min_gm_T1, navgs) ;
  smooth_map(mris, ep->cv_max_gm_T1, navgs) ;
  if (Gdiag_no >= 0) {
    printf("after smoothing, wm PD (%2.1f --> %2.1f), T1: (%2.1f --> %2.1f)\n",
           ep->cv_min_wm_PD[Gdiag_no], ep->cv_max_wm_PD[Gdiag_no],
           ep->cv_min_wm_T1[Gdiag_no], ep->cv_max_wm_T1[Gdiag_no]) ;
    printf("after smoothing, gm PD (%2.1f --> %2.1f), T1: (%2.1f --> %2.1f)\n",
           ep->cv_min_gm_PD[Gdiag_no], ep->cv_max_gm_PD[Gdiag_no],
           ep->cv_min_gm_T1[Gdiag_no], ep->cv_max_gm_T1[Gdiag_no]) ;
  }

  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  return(NO_ERROR) ;
}
#endif
