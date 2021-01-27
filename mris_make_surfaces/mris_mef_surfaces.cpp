/**
 * @brief repositions cortical surface to gray/white boundary
 *
 * This program positions the tessellation of the cortical surface
 * at the white matter surface, then the gray matter surface
 * and generates surface files for these surfaces as well as a
 * 'curvature' file for the cortical thickness, and a surface file
 * which approximates layer IV of the cortical sheet.
 */
/*
 * Original Author: Bruce Fischl
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
#include "tags.h"
#include "mrinorm.h"
#include "version.h"
#include "label.h"


int main(int argc, char *argv[]) ;

//#define BRIGHT_LABEL         130
//#define BRIGHT_BORDER_LABEL  100

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

MRI *MRIfillVentricle(MRI *mri_inv_lv, MRI *mri_T1, float thresh,
                      int out_label, MRI *mri_dst);

int LocalMRISfindExpansionRegions(MRI_SURFACE *mris) ;
int MRIsmoothBrightWM(MRI *mri_T1, MRI *mri_wm) ;
MRI *LocalMRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm) ;

int MRISaverageMarkedValbaks(MRI_SURFACE *mris, int navgs);
static int  MRIcomputeClassStatistics_mef(MRI *mri_T1_30, 
                                          MRI *mri_T1_5, 
                                          MRI *mri_em_seg, 
                                          float *white_mean, float *white_std, 
                                          float *gray_mean, float *gray_std) ;

int
MRIScomputeBorderValues_MEF_WHITE(MRI_SURFACE *mris, 
                                  MRI *mri_em_combined, MRI *mri_30,
                                  MRI *mri_5, 
                                  float wm_mean[2], float wm_std[2],
                                  float gm_mean[2], float gm_std[2],
                                  double sigma,
                                  float max_thickness, FILE *log_fp);

int
MRIScomputeBorderValues_MEF_PIAL(MRI_SURFACE *mris, 
                                 MRI *mri_em_combined, MRI *mri_30,
                                 MRI *mri_5, float wm_mean[2], float wm_std[2],
                                 float gm_mean[2], float gm_std[2],
                                 double sigma,
                                 float max_thickness, FILE *log_fp);

// int MRISpositionSurface_mef(MRI_SURFACE *mris, MRI *mri_30, MRI *mri_5, INTEGRATION_PARMS *parms, float weight30, float weight5);

// static double mrisRmsValError_mef(MRI_SURFACE *mris, MRI *mri_30, MRI *mri_5, float weight30, float weight5);

// static int mrisComputeIntensityTerm_mef(MRI_SURFACE *mris, double l_intensity, MRI *mri_30, MRI *mri_5, double sigma_global, float weight30, float weight5);

static int  MRInormalizeMEF(MRI *mri, MRI *mri_em_seg);

static LABEL *highres_label = NULL ;

static char T1_30_name[STRLEN] = 
"flash30_T1" ; //INU corrected flash30, can use EM's output
static char T1_5_name[STRLEN] = 
"flash5_T1" ; //INU corrected flash5, can use EM's output
static char em_name[STRLEN] = 
"atlas_EM_combined" ; //synthesized volume from EM segmentation

static char *white_fname = NULL ;

static char *orig_white = NULL ;
static char *orig_pial = NULL ;

const char *Progname ;

static double std_scale = 1.0;

static int graymid = 0 ;
static int curvature_avgs = 10 ;
static int create = 1 ;
static int smoothwm = 0 ;
static int white_only = 0 ;
static int overlay = 0 ;

static int auto_detect_stats = 1 ;

static int in_out_in_flag = 0 ;
static int apply_median_filter = 0 ;

static int nbhd_size = 20 ;

static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE    1.0
static float base_dt_scale = BASE_DT_SCALE ;

static int add = 0 ;

static double l_tsmooth = 0.0 ;
static double l_surf_repulse = 5.0 ;

static int smooth = 5 ;
static int vavgs = 5 ;
static int nwhite = 20 /*5*/ ;
static int ngray = 30 /*45*/ ;

static int nowhite = 0 ;
static int nbrs = 2 ;
static int write_vals = 0 ;

static const char *orig_name = ORIG_NAME ;
static const char *suffix = "" ;
static const char *output_suffix = "" ;
static char *xform_fname = NULL ;

static char pial_name[STRLEN] = "pial" ;
static char white_matter_name[STRLEN] = WHITE_MATTER_NAME ;

static int lh_label = LH_LABEL ;
static int rh_label = RH_LABEL ;

static int max_pial_averages = 16 ;
static int min_pial_averages = 2 ;
static int max_white_averages = 4 ;
static int min_white_averages = 0 ;
static float pial_sigma = 2.0f ;
static float white_sigma = 2.0f ;
static float max_thickness = 5.0 ;


static char sdir[STRLEN] = "" ;

static int MGZ = 0; // for use with MGZ format

static int longitudinal = 0;


int
main(int argc, char *argv[]) {
  char           *hemi, *sname, *cp, fname[STRLEN], mdir[STRLEN];
  int           nargs, i, replace_val, msec, n_averages, j ;
  MRI_SURFACE   *mris ;
  MRI           *mri_wm, *mri_filled, *mri_T1_30, *mri_T1_5; // *mri_labeled;
  MRI           *mri_em_seg;
  float         max_len ;
  //need to be estimated from EM segmentation; 
  // need to mask out cerebellum though
  float         white_mean[2], white_std[2], gray_mean[2], gray_std[2];
  double        current_sigma ;
  Timer then ;

  std::string cmdline = getAllInfo(argc, argv, "mris_mef_surfaces");

  nargs = handleVersionOption(argc, argv, "mris_mef_surfaces");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.projection = NO_PROJECTION ;
  parms.tol = 1e-4 ;
  parms.dt = 0.5f ;
  parms.base_dt = BASE_DT_SCALE*parms.dt ;
  parms.l_spring = 1.0f ;
  parms.l_curv = 1.0 ;
  parms.l_intensity = 0.2 ;
  parms.l_spring = 0.0f ;
  parms.l_curv = 1.0 ;
  parms.l_intensity = 0.2 ;
  parms.l_tspring = 1.0f ;
  parms.l_nspring = 0.5f ;

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

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit() ;

  /* set default parameters for white and gray matter surfaces */
  parms.niterations = nwhite ;
  if (parms.momentum < 0.0)
    parms.momentum = 0.0 /*0.75*/ ;

  then.reset() ;
  sname = argv[1] ;
  hemi = argv[2] ;
  if (!strlen(sdir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  cp = getenv("FREESURFER_HOME") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM,
              "%s: FREESURFER_HOME not defined in environment.\n", Progname) ;
  strcpy(mdir, cp) ;

  int req = snprintf(fname, STRLEN, "%s/%s/mri/filled", sdir, sname) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (MGZ) strcat(fname, ".mgz");
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_filled = MRIread(fname) ;
  if (!mri_filled)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;
  ////////////////////////////// we can handle only conformed volumes
  setMRIforSurface(mri_filled);

  if (!stricmp(hemi, "lh")) {
    replace_val = rh_label ;
  } else {
    replace_val = lh_label ;
  }

  req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, T1_30_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (MGZ) strcat(fname, ".mgz");
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_T1_30 = MRIread(fname) ;

  if (!mri_T1_30)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

  req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, T1_5_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (MGZ) strcat(fname, ".mgz");
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_T1_5 = MRIread(fname) ;

  if (!mri_T1_5)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

  req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, em_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (MGZ) strcat(fname, ".mgz");
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_em_seg = MRIread(fname) ;

  if (!mri_em_seg)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

  /////////////////////////////////////////
  setMRIforSurface(mri_T1_30);
  setMRIforSurface(mri_T1_5);

  if (apply_median_filter) {
    MRI *mri_tmp ;

    fprintf(stderr, "applying median filter to T1 image...\n") ;
    mri_tmp = MRImedian(mri_T1_30, NULL, 3, NULL) ;
    MRIfree(&mri_T1_30) ;
    mri_T1_30 = mri_tmp ;

    mri_tmp = MRImedian(mri_T1_5, NULL, 3, NULL) ;
    MRIfree(&mri_T1_5) ;
    mri_T1_5 = mri_tmp ;
  }

  //normalize each input volume to have same mean and std within WM
  MRInormalizeMEF(mri_T1_30, mri_em_seg);
  MRInormalizeMEF(mri_T1_5, mri_em_seg);

  //  MRIwrite(mri_T1_30, "normalized_30.mgz");
  // MRIwrite(mri_T1_5, "normalized_5.mgz");


  /* remove other hemi */
  MRIdilateLabel(mri_filled, mri_filled, replace_val, 1) ;

  MRImask(mri_T1_30, mri_filled, mri_T1_30, replace_val,0) ;
  MRImask(mri_T1_5, mri_filled, mri_T1_5, replace_val,0) ;

  MRIfree(&mri_filled) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_T1_30, "r30.mgz") ;

  req = snprintf(fname, STRLEN, "%s/%s/mri/wm", sdir, sname) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (MGZ) strcat(fname, ".mgz");
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_wm = MRIread(fname) ;
  if (!mri_wm)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;
  //////////////////////////////////////////
  setMRIforSurface(mri_wm);

  //the following two steps might as well be ignored
  //  MRIsmoothBrightWM(mri_T1_30, mri_wm) ;
  //  mri_labeled = MRIfindBrightNonWM(mri_T1_30, mri_wm) ;

  MRIcomputeClassStatistics_mef(mri_T1_30, mri_T1_5, mri_em_seg, 
                                white_mean, white_std, gray_mean, gray_std) ;

  printf("WM: mean = (%g, %g) std = (%g, %g)\n", 
         white_mean[0], white_mean[1], white_std[0], white_std[1]);
  printf("GM: mean = (%g, %g) std = (%g, %g)\n", 
         gray_mean[0], gray_mean[1], gray_std[0], gray_std[1]);

  MRIfree(&mri_wm) ;
  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s", sdir, sname, hemi, orig_name, suffix) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fprintf(stderr, "reading original surface position from %s...\n", fname) ;
  mris = MRISreadOverAlloc(fname, 1.1) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  MRISaddCommandLine(mris, cmdline) ;

  if (smooth && !nowhite) {
    printf("smoothing surface for %d iterations...\n", smooth) ;
    MRISaverageVertexPositions(mris, smooth) ;
  }

  if (nbrs > 1)
    MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;

  sprintf(parms.base_name, "%s%s%s", 
          white_matter_name, output_suffix, suffix) ;
  if (orig_white) {
    printf("reading initial white vertex positions from %s...\n", orig_white) ;
    if (MRISreadVertexPositions(mris, orig_white) != NO_ERROR)
      ErrorExit(Gerror, "reading of orig white failed...");
  }
  MRIScomputeMetricProperties(mris) ;    /* recompute surface normals */
  MRISstoreMetricProperties(mris) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;

  if (add) {
    fprintf(stderr, "adding vertices to initial tessellation...\n") ;
    for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
    while (MRISdivideLongEdges(mris, max_len) > 0) {}
  }

  //what if I need a vector-values intensity??
  MRISsetVals(mris, -1) ;  /* clear white matter intensities */
  MRIScopyValToValBak(mris);

#if 0
  if (!nowhite) {
    fprintf(stderr, "repositioning cortical surface to gray/white boundary\n");

    MRImask(mri_T1_30, mri_labeled, mri_T1_30, BRIGHT_LABEL, 0) ;
    MRImask(mri_T1_30, mri_labeled, mri_T1_30, BRIGHT_BORDER_LABEL, 0) ;

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_T1_30, "white_masked.mgz") ;
  }
#endif

  current_sigma = white_sigma ;
  for (n_averages = max_white_averages, i = 0 ;
       n_averages >= min_white_averages ;
       n_averages /= 2, current_sigma /= 2, i++) {
    if (nowhite)
      break ;

    parms.sigma = current_sigma ;

    parms.n_averages = n_averages ;
    MRISprintTessellationStats(mris, stderr) ;

    //This following function needs major modification
    MRIScomputeBorderValues_MEF_WHITE(mris, mri_em_seg, mri_T1_30, mri_T1_5,
                                      white_mean, white_std, gray_mean,
                                      gray_std,
                                      current_sigma,
                                      2*max_thickness, parms.fp) ;

    //what does this do?
    LocalMRISfindExpansionRegions(mris) ;

    if (vavgs) {
      fprintf(stderr, "averaging target values for %d iterations...\n",vavgs) ;
      //the following function needs be modified too to use two channels
      MRISaverageMarkedVals(mris, vavgs) ;
      MRISaverageMarkedValbaks(mris, vavgs) ;
      if (Gdiag_no > 0) {
        VERTEX *v ;
        v = &mris->vertices[Gdiag_no] ;
        fprintf(stderr,"v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f\n",
                Gdiag_no, v->val, v->mean, v->d) ;
      }
    }


    if (write_vals) {
      sprintf(fname, "./%s-white%2.2f.w", hemi, current_sigma) ;
      //MRISwriteValues(mris, fname) ;
      MRISwrite(mris, fname);
    }

    //the following function needs major change
    //and may need different versions for white and pial resp.
    MRISpositionSurface_mef(mris, mri_T1_30, mri_T1_5, &parms, 0.75, 0.25);

    if (add) {
      for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
      while (MRISdivideLongEdges(mris, max_len) > 0) {}
    }
    if (!n_averages)
      break ;
  }

  if (!nowhite) {
    req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s%s", sdir, sname,hemi,white_matter_name,
            output_suffix,suffix);
    if (req >= STRLEN) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stderr, "writing white matter surface to %s...\n", fname) ;
    MRISaverageVertexPositions(mris, smoothwm) ;
    MRISwrite(mris, fname) ;

    if (create)   /* write out curvature and area files */
    {
      MRIScomputeMetricProperties(mris) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISuseMeanCurvature(mris) ;
      MRISaverageCurvatures(mris, curvature_avgs) ;
      sprintf(fname, "%s.curv%s%s",
              mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", output_suffix,
              suffix);
      fprintf(stderr, "writing smoothed curvature to %s\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
      sprintf(fname, "%s.area%s%s",
              mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", output_suffix,
              suffix);
      MRISprintTessellationStats(mris, stderr) ;
    }
  } else   /* read in previously generated white matter surface */
  {
    sprintf(fname, "%s%s", white_matter_name, suffix) ;
    if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
      ErrorExit(Gerror, "%s: could not read white matter surfaces.",
                Progname) ;
    MRIScomputeMetricProperties(mris) ;
  }


  if (white_only) {
    msec = then.milliseconds() ;
    fprintf(stderr,
            "refinement took %2.1f minutes\n", (float)msec/(60*1000.0f));
    MRIfree(&mri_T1_30);
    MRIfree(&mri_T1_5);
    MRIfree(&mri_em_seg);
    exit(0) ;
  }

  ////////////////////////////////////////////////////////////////////////////
  // pial surface
  ////////////////////////////////////////////////////////////////////////////
  parms.t = parms.start_t = 0 ;
  sprintf(parms.base_name, "%s%s%s", pial_name, output_suffix, suffix) ;
  parms.niterations = ngray ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ; /* save white-matter */
  parms.l_surf_repulse = l_surf_repulse ;

  MRISsetVals(mris, -1) ;  /* clear target intensities */
  MRIScopyValToValBak(mris);

#if 0
  if (smooth) //smooth again?? seems too much
  {
    printf("smoothing surface for %d iterations...\n", smooth) ;
    MRISaverageVertexPositions(mris, smooth) ;
  }
#endif

  fprintf(stderr, "repositioning cortical surface to gray/csf boundary.\n") ;
  parms.l_repulse = 0 ;

  if (orig_pial) {
    printf("reading initial pial vertex positions from %s...\n", orig_pial) ;

    if (MRISreadVertexPositions(mris, orig_pial) != NO_ERROR)
      ErrorExit(Gerror, "reading orig pial positions failed") ;

    MRIScomputeMetricProperties(mris) ; /*shouldn't this be done whenever 
                                          orig_pial is used??? Maybe that's 
                                          why the cross-intersection 
                                          was caused */
  }

  /*    parms.l_convex = 1000 ;*/

  for (j = 0 ; j <= 0 ; parms.l_intensity *= 2, j++)  /* only once for now */
  {
    current_sigma = pial_sigma ;
    for (n_averages = max_pial_averages, i = 0 ;
         n_averages >= min_pial_averages ;
         n_averages /= 2, current_sigma /= 2, i++) {

      parms.sigma = current_sigma ;
      parms.n_averages = n_averages ;
      parms.l_tsmooth = l_tsmooth ;
      /*
        replace bright stuff such as eye sockets with 255. 
        Simply zeroing it out
        would make the border always go through the sockets, and ignore subtle
        local minima in intensity at the border of the sockets. Will set to 0
        after border values have been computed so that it 
        doesn't mess up gradients.
      */
      //   MRImask(mri_T1_30, mri_labeled, mri_T1_30, BRIGHT_LABEL, 255) ;
      //  MRImask(mri_T1_30, mri_labeled, mri_T1_30, BRIGHT_BORDER_LABEL, MID_GRAY) ;
      //  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      //  MRIwrite(mri_T1_30, "pial_masked.mgz") ;

      //The following function need major modification
      MRIScomputeBorderValues_MEF_PIAL
        (mris, mri_em_seg, mri_T1_30, mri_T1_5,
         white_mean, white_std, gray_mean, gray_std,
         current_sigma, 2*max_thickness, parms.fp) ;

      //MRImask(mri_T1_30, mri_labeled, mri_T1_30, BRIGHT_LABEL, 0) ;

      if (vavgs) {
        fprintf(stderr, 
                "averaging target values for %d iterations...\n",vavgs) ;
        MRISaverageMarkedVals(mris, vavgs) ;
        MRISaverageMarkedValbaks(mris, vavgs) ;
        if (Gdiag_no > 0) {
          VERTEX *v ;
          v = &mris->vertices[Gdiag_no] ;
          fprintf(stderr,
                  "v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f\n",
                  Gdiag_no, v->val, v->mean, v->d) ;
        }
      }

      if (write_vals) {
        sprintf(fname, "./%s-gray%2.2f.w", hemi, current_sigma) ;
        MRISwriteValues(mris, fname) ;
      }

      //The following function need major modification
      MRISpositionSurface_mef(mris, mri_T1_30, mri_T1_5, &parms, 0.95, 0.05);
      /*    parms.l_nspring = 0 ;*/
      if (!n_averages)
        break ;
    }
  }

  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s%s", sdir, sname, hemi, pial_name,
          output_suffix, suffix) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fprintf(stderr, "writing pial surface to %s...\n", fname) ;
  MRISwrite(mris, fname) ;

  MRIfree(&mri_T1_30);
  MRIfree(&mri_T1_5);
  MRIfree(&mri_em_seg);


  {
    fprintf(stderr, "measuring cortical thickness...\n") ;
    MRISmeasureCorticalThickness(mris, nbhd_size, max_thickness) ;
    fprintf(stderr,
            "writing cortical thickness estimate to 'thickness' file.\n") ;
    sprintf(fname, "thickness%s%s", output_suffix, suffix) ;
    MRISwriteCurvature(mris, fname) ;

    /* at this point, the v->curv slots contain the cortical surface. Now
       move the white matter surface out by 1/2 the thickness as an estimate
       of layer IV.
    */
    if (graymid) {
      MRISsaveVertexPositions(mris, TMP_VERTICES) ;
      mrisFindMiddleOfGray(mris) ;
      req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s", sdir, sname, hemi, GRAYMID_NAME,
              suffix) ;
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stderr, "writing layer IV surface to %s...\n", fname) ;
      MRISwrite(mris, fname) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    }
  }
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
  else if (!stricmp(option, "nbrs")) {
    nbrs = atoi(argv[2]) ;
    fprintf(stderr,  "using neighborhood size = %d\n", nbrs) ;
    nargs = 1 ;
  } else if (!stricmp(option, "wvol")) {
    white_fname = argv[2] ;
    printf("using %s as volume for white matter deformation...\n", 
           white_fname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "hires") || !stricmp(option, "highres")) {
    highres_label = LabelRead(NULL, argv[2]) ;
    if (!highres_label)
      ErrorExit(ERROR_NOFILE, 
                "%s: could not read highres label %s", Progname, argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "long")) {
    longitudinal = 1;
    printf("Using longitudinal scheme\n");
  } else if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "orig_white")) {
    orig_white = argv[2] ;
    printf("using %s starting white location...\n", orig_white) ;
    nargs = 1 ;
  } else if (!stricmp(option, "orig_pial")) {
    orig_pial = argv[2] ;
    printf("using %s starting pial locations...\n", orig_pial) ;
    nargs = 1 ;
  } else if (!stricmp(option, "flash30") ||
             !stricmp(option, "T130") ||
             !stricmp(option, "T1_30")
            ) {
    strcpy(T1_30_name, argv[2]) ;
    printf("using %s as (INU-normalized) flash30 volume...\n", T1_30_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "flash5") ||
             !stricmp(option, "T15") ||
             !stricmp(option, "T1_5")
            ) {
    strcpy(T1_5_name, argv[2]) ;
    printf("using %s as (INU-normalized) flash5 volume...\n", T1_5_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "em") ||
             !stricmp(option, "em_seg") ||
             !stricmp(option, "em_combined")
            ) {
    strcpy(em_name, argv[2]) ;
    printf("using %s as the EM tissue segmentation volume...\n", em_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "median")) {
    apply_median_filter = 1 ;
  } else if (!stricmp(option, "scale_std")) {
    std_scale = atof(argv[2]);
    printf("scale the estimated WM and GM std by %g \n", std_scale) ;
    nargs = 1 ;
  } else if (!stricmp(option, "noauto")) {
    auto_detect_stats = 0 ;
    fprintf(stderr, "disabling auto-detection of border ranges...\n") ;
  } else if (!stricmp(option, "inoutin")) {
    in_out_in_flag = 1 ;
    fprintf(stderr, "applying final white matter deformation after pial\n") ;
  } else if (!stricmp(option, "graymid")) {
    graymid = 1 ;
    fprintf(stderr, "generating graymid surface...\n") ;
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
    fprintf(stderr,  "only generating white matter surface\n") ;
  } else if (!stricmp(option, "overlay")) {
    overlay = !overlay ;
    fprintf(stderr,  "%soverlaying T1 volume with edited white matter\n",
            overlay ? "" : "not") ;
  } else if (!stricmp(option, "pial")) {
    strcpy(pial_name, argv[2]) ;
    fprintf(stderr,  "writing pial surface to file named %s\n", pial_name) ;
  } else if (!stricmp(option, "write_vals")) {
    write_vals = 1 ;
    fprintf(stderr,  "writing gray and white surface targets to w files\n") ;
  } else if (!stricmp(option, "name")) {
    strcpy(parms.base_name, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "base name = %s\n", parms.base_name) ;
  } else if (!stricmp(option, "dt")) {
    parms.dt = atof(argv[2]) ;
    parms.base_dt = base_dt_scale*parms.dt ;
    parms.integration_type = INTEGRATE_MOMENTUM ;
    fprintf(stderr,  "using dt = %2.1e\n", parms.dt) ;
    nargs = 1 ;
  } else if (!stricmp(option, "spring")) {
    parms.l_spring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_spring = %2.3f\n", parms.l_spring) ;
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
  } else if (!stricmp(option, "vavgs")) {
    vavgs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "smoothing values for %d iterations\n", vavgs) ;
  } else if (!stricmp(option, "white")) {
    strcpy(white_matter_name, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using %s as white matter name...\n", white_matter_name) ;
  } else if (!stricmp(option, "intensity")) {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_intensity = %2.3f\n", parms.l_intensity) ;
  } else if (!stricmp(option, "lm")) {
    parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
    fprintf(stderr, "integrating with line minimization\n") ;
  } else if (!stricmp(option, "nwhite")) {
    nwhite = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,
            "integrating gray/white surface positioning for %d time steps\n",
            nwhite) ;
  } else if (!stricmp(option, "nowhite")) {
    nowhite = 1 ;
    fprintf(stderr, "reading previously compute gray/white surface\n") ;
  } else if (!stricmp(option, "smoothwm")) {
    smoothwm = atoi(argv[2]) ;
    fprintf(stderr, "writing smoothed (%d iterations) wm surface\n",
            smoothwm) ;
    nargs = 1 ;
  } else if (!stricmp(option, "ngray")) {
    ngray = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"integrating pial surface positioning for %d time steps\n",
            ngray) ;
  } else if (!stricmp(option, "wsigma")) {
    white_sigma = atof(argv[2]) ;
    fprintf(stderr,  "smoothing volume with Gaussian sigma = %2.1f\n",
            white_sigma) ;
    nargs = 1 ;
  } else if (!stricmp(option, "psigma")) {
    pial_sigma = atof(argv[2]) ;
    fprintf(stderr,  "smoothing volume with Gaussian sigma = %2.1f\n",
            pial_sigma) ;
    nargs = 1 ;
  } else if (!stricmp(option, "pa")) {
    max_pial_averages = atoi(argv[2]) ;
    fprintf(stderr, "using max pial averages = %d\n", max_pial_averages) ;
    nargs = 1 ;
    if (isdigit(*argv[3])) {
      min_pial_averages = atoi(argv[3]) ;
      fprintf(stderr, "using min pial averages = %d\n", min_pial_averages) ;
      nargs++ ;
    }
  } else if (!stricmp(option, "wa")) {
    max_white_averages = atoi(argv[2]) ;
    fprintf(stderr, "using max white averages = %d\n", max_white_averages) ;
    nargs = 1 ;
    if (isdigit(*argv[3])) {
      min_white_averages = atoi(argv[3]) ;
      fprintf(stderr, "using min white averages = %d\n", min_white_averages) ;
      nargs++ ;
    }
  } else if (!stricmp(option, "add")) {
    add = 1 ;
    fprintf(stderr, "adding vertices to tessellation during deformation.\n");
  } else if (!stricmp(option, "max")) {
    max_thickness = atof(argv[2]) ;
    nargs = 1 ;
    printf("using max_thickness = %2.1f\n", max_thickness) ;
  } else if (!stricmp(option, "mgz")) {
    MGZ = 1;
    printf("INFO: assuming MGZ format for volumes.\n");
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
    case 'T':
      xform_fname = argv[2] ;
      nargs = 1;
      fprintf(stderr, "applying ventricular xform %s\n", xform_fname);
      break ;
    case 'O':
      orig_name = argv[2] ;
      nargs = 1 ;
      fprintf(stderr, 
              "reading original vertex positions from %s\n", orig_name);
      break ;
    case 'Q':
      parms.flags |= IPFLAG_NO_SELF_INT_TEST ;
      fprintf(stderr,
              "doing quick (no self-intersection) surface positioning.\n") ;
      break ;
#if 0
    case 'A':
      max_averages = atoi(argv[2]) ;
      fprintf(stderr, "using n_averages = %d\n", max_averages) ;
      nargs = 1 ;
      if (isdigit(*argv[3])) {
        min_averages = atoi(argv[3]) ;
        fprintf(stderr, "using min_averages = %d\n", min_averages) ;
        nargs++ ;
      }
      break ;
#endif
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
    case 'C':
      create = !create ;
      fprintf(stderr, 
              "%screating area and curvature files for wm surface...\n",
              create ? "" : "not ") ;
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
  fprintf(stderr, "usage: %s [options] <subject name> <hemisphere>\n",
          Progname) ;
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

MRI *
MRIfillVentricle(MRI *mri_inv_lv, MRI *mri_T1, float thresh,
                 int out_label, MRI *mri_dst) {
  BUFTYPE   *pdst, *pinv_lv, out_val, inv_lv_val, *pT1 ;
  int       width, height, depth, x, y, z,
  ventricle_voxels;

  if (!mri_dst)
    mri_dst = MRIclone(mri_T1, NULL) ;

  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;
  /* now apply the inverse morph to build an average wm representation
     of the input volume
  */


  ventricle_voxels = 0 ;
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      pT1 = &MRIvox(mri_T1, 0, y, z) ;
      pinv_lv = &MRIvox(mri_inv_lv, 0, y, z) ;
      for (x = 0 ; x < width ; x++) {
        pT1++ ;
        inv_lv_val = *pinv_lv++ ;
        out_val = 0 ;
        if (inv_lv_val >= thresh) {
          ventricle_voxels++ ;
          out_val = out_label ;
        }
        *pdst++ = out_val ;
      }
    }
  }

#if 0
  MRIfillRegion(mri_T1, mri_dst, 30, out_label, 2*ventricle_voxels) ;
  MRIdilate(mri_dst, mri_dst) ;
  MRIdilate(mri_dst, mri_dst) ;
#endif
  return(mri_dst) ;
}

int
LocalMRISfindExpansionRegions(MRI_SURFACE *mris) {
  int    vno, num, n, num_long, total ;
  float  d, dsq, mean, std, dist ;

  d = dsq = 0.0f ;
  for (total = num = vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX const * const v = &mris->vertices[vno] ;
    if (v->ripflag || v->val <= 0)
      continue ;
    num++ ;
    dist = fabs(v->d) ;
    d += dist ;
    dsq += (dist*dist) ;
  }

  mean = d / num ;
  std = sqrt(dsq/num - mean*mean) ;
  fprintf(stderr, "mean absolute distance = %2.2f +- %2.2f\n", mean, std) ;

  for (num = vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    v->curv = 0 ;
    if (v->ripflag || v->val <= 0)
      continue ;
    if (fabs(v->d) < mean+2*std)
      continue ;
    for (num_long = num = 1, n = 0 ; n < vt->vnum ; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      if (vn->val <= 0 || v->ripflag)
        continue ;
      if (fabs(vn->d) >= mean+2*std)
        num_long++ ;
      num++ ;
    }

    if ((float)num_long / (float)num > 0.25) {
      v->curv = fabs(v->d) ;
      total++ ;
#if 0
      fprintf(stderr, "v %d long: (%d of %d)\n", vno, num_long, num) ;
#endif
    }
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "%d vertices more than 2 sigmas from mean.\n", total) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRISwriteCurvature(mris, "long") ;
  return(NO_ERROR) ;
}

int
MRIsmoothBrightWM(MRI *mri_T1, MRI *mri_wm) {
  int     width, height, depth, x, y, z, nthresholded ;
  BUFTYPE *pT1, *pwm, val, wm ;

  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;

  nthresholded = 0 ;
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      pT1 = &MRIvox(mri_T1, 0, y, z) ;
      pwm = &MRIvox(mri_wm, 0, y, z) ;
      for (x = 0 ; x < width ; x++) {
        val = *pT1 ;
        wm = *pwm++ ;
        if (wm >= WM_MIN_VAL)  /* labeled as white */
        {
          if (val > DEFAULT_DESIRED_WHITE_MATTER_VALUE) {
            nthresholded++ ;
            val = DEFAULT_DESIRED_WHITE_MATTER_VALUE ;
          }
        }
        *pT1++ = val ;
      }
    }
  }

  fprintf(stderr, "%d bright wm thresholded.\n", nthresholded) ;

  return(NO_ERROR) ;
}
MRI *
LocalMRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm) {
  int     width, height, depth, x, y, z, nlabeled, nwhite,
  xk, yk, zk, xi, yi, zi;
  BUFTYPE *pT1, *pwm, val, wm ;
  MRI     *mri_labeled, *mri_tmp ;

  mri_labeled = MRIclone(mri_T1, NULL) ;
  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      pT1 = &MRIvox(mri_T1, 0, y, z) ;
      pwm = &MRIvox(mri_wm, 0, y, z) ;
      for (x = 0 ; x < width ; x++) {
        val = *pT1++ ;
        wm = *pwm++ ;

        if (x == Gx && y == Gy && z == Gz)  /* T1=127 */
          DiagBreak() ;
        /* not white matter and bright (e.g. eye sockets) */
        if ((wm < WM_MIN_VAL) && (val > 125)) {
          nwhite = 0 ;
          for (xk = -1 ; xk <= 1 ; xk++) {
            xi = mri_T1->xi[x+xk] ;
            for (yk = -1 ; yk <= 1 ; yk++) {
              yi = mri_T1->yi[y+yk] ;
              for (zk = -1 ; zk <= 1 ; zk++) {
                zi = mri_T1->zi[z+zk] ;
                if (MRIvox(mri_wm, xi, yi, zi) >= WM_MIN_VAL)
                  nwhite++ ;
              }
            }
          }
#define MIN_WHITE  ((3*3*3-1)/2)
          if (nwhite < MIN_WHITE)
            MRIvox(mri_labeled, x, y, z) = BRIGHT_LABEL ;
        }
      }
    }
  }

  /* find all connected voxels that are above 115 */
  MRIdilateThreshLabel(mri_labeled, mri_T1, NULL, BRIGHT_LABEL, 10,115);
  MRIclose(mri_labeled, mri_labeled) ;

  /* expand once more to all neighboring voxels that are bright. At
     worst we will erase one voxel of white matter.
  */
  mri_tmp =
    MRIdilateThreshLabel(mri_labeled, mri_T1, NULL, BRIGHT_LABEL,1,100);
  MRIxor(mri_labeled, mri_tmp, mri_tmp, 1, 255) ;
  MRIreplaceValues(mri_tmp, mri_tmp, 1, BRIGHT_BORDER_LABEL) ;
  MRIunion(mri_tmp, mri_labeled, mri_labeled) ;
#if 0
  fprintf(stderr, "selectively smoothing volume....\n") ;
  MRIsoapBubbleLabel(mri_T1, mri_labeled, mri_T1, BRIGHT_LABEL, 200) ;
#endif
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_labeled, "label.mgz") ;
  /*    MRIwrite(mri_tmp, "tmp.mgz") ;*/
  nlabeled = MRIvoxelsInLabel(mri_labeled, BRIGHT_LABEL) ;
  fprintf(stderr, "%d bright non-wm voxels segmented.\n", nlabeled) ;

  /* dilate outwards if exactly 0 */
  MRIdilateInvThreshLabel(mri_labeled, mri_T1, mri_labeled, 
                          BRIGHT_LABEL, 3, 0) ;

  MRIfree(&mri_tmp) ;
  return(mri_labeled) ;
}

#define STEP_SIZE 0.1

int
MRIScomputeBorderValues_MEF_WHITE(MRI_SURFACE *mris, 
                                  MRI *mri_em_combined, MRI *mri_30,
                                  MRI *mri_5, 
                                  float wm_mean[2], float wm_std[2],
                                  float gm_mean[2], float gm_std[2],
                                  double sigma,
                                  float max_thickness, FILE *log_fp) {
  double  val30, val5, x, y, z, max_mag_val30, 
    max_mag_val5, xw, yw, zw,mag,max_mag, max_mag_dist=0.0f,
    min_val30, min_val5, inward_dist,outward_dist,xw1,yw1,zw1,
    min_val_dist, orig_dist, dx, dy, dz;
  double previous_mag, next_mag, previous_mag30, 
    previous_mag5, next_mag30, next_mag5;
  double mag30, mag5, previous_val30, next_val30;
  int     total_vertices, vno, nmissing = 0, nout = 0, nin = 0, nfound = 0,
    nalways_missing = 0, local_max_found, 
    ngrad_max, ngrad, nmin, num_changed=0 ;
  float   mean_border, mean_in, mean_out, 
    dist, nx, ny, nz, mean_dist, step_size ;
  double  current_sigma ;
  VERTEX  *v ;
  FILE    *fp = NULL ;

  step_size = mri_30->xsize/2 ;

  /* first compute intensity of local gray/white boundary */
  mean_dist = mean_in = mean_out = mean_border = 0.0f ;
  ngrad_max = ngrad = nmin = 0 ;
  MRISclearMarks(mris) ;  /* for soap bubble smoothing later */
  for (total_vertices = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ;
    y = v->y ;
    z = v->z ;
    MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
    x = v->x + v->nx ;
    y = v->y + v->ny ;
    z = v->z + v->nz ;
    MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw1, &yw1, &zw1) ;
    nx = xw1 - xw ;
    ny = yw1 - yw ;
    nz = zw1 - zw ;
    dist = sqrt(SQR(nx)+SQR(ny)+SQR(nz)) ;
    if (FZERO(dist))
      dist = 1 ;
    nx /= dist ;
    ny /= dist ;
    nz /= dist ;


    /*
    find the distance in the directions parallel and anti-parallel to
    the surface normal in which the gradient is pointing 'inwards'.
    The border will then be constrained to be within that region.
    */
    inward_dist = 1.0 ;
    outward_dist = -1.0 ;
    for (current_sigma = sigma; 
         current_sigma <= 10*sigma; 
         current_sigma *= 2) {
      for (dist = 0 ; dist > -max_thickness ; dist -= step_size) {
        dx = v->x-v->origx ;
        dy = v->y-v->origy ;
        dz = v->z-v->origz ;
        orig_dist = fabs(dx*v->nx + dy*v->ny + dz*v->nz) ;
        if (fabs(dist)+orig_dist > max_thickness)
          break ;
        x = v->x + v->nx*dist ;
        y = v->y + v->ny*dist ;
        z = v->z + v->nz*dist ;
        MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, 
                                       nx, ny,nz,&mag30,  current_sigma);
        //       MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, 
        // nx, ny,nz,&mag5,  current_sigma);

        //       if (mag30 >= 0.0 || mag5 <= 0)
        if (mag30 >= 0.0) //this condition doesn't make much sense to me
          break ;
        MRIsampleVolume(mri_30, xw, yw, zw, &val30) ;
        // MRIsampleVolume(mri_5, xw, yw, zw, &val5) ;

        if (val30 > (wm_mean[0] + wm_std[0]))
          break ;
      }
      inward_dist = dist+step_size/2 ;
      for (dist = 0 ; dist < max_thickness ; dist += step_size) {
        dx = v->x-v->origx ;
        dy = v->y-v->origy ;
        dz = v->z-v->origz ;
        orig_dist = fabs(dx*v->nx + dy*v->ny + dz*v->nz) ;
        if (fabs(dist)+orig_dist > max_thickness)
          break ;
        x = v->x + v->nx*dist ;
        y = v->y + v->ny*dist ;
        z = v->z + v->nz*dist ;
        MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, 
                                       nx, ny,nz, &mag30, current_sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, 
                                       nx, ny,nz, &mag5, current_sigma);

        if (mag30 >= 0.0)
          break ;
        MRIsampleVolume(mri_30, xw, yw, zw, &val30) ;
        MRIsampleVolume(mri_5, xw, yw, zw, &val5) ;
        if (val30 < gm_mean[0] || val5 > (gm_mean[1] + gm_std[1]))
          break ;
      }
      outward_dist = dist-step_size/2 ;
      if (!std::isfinite(outward_dist))
        DiagBreak() ;
      if (inward_dist <= 0 || outward_dist >= 0)
        break ;
    }

    if (inward_dist > 0 && outward_dist < 0)
      current_sigma = sigma ;  /* couldn't find anything */

    if (vno == Gdiag_no) {
      char fname[STRLEN] ;
      sprintf(fname, "v%d.%2.0f.log", Gdiag_no, sigma*100) ;
      fp = fopen(fname, "w") ;
      fprintf(stdout,
              "v %d: inward dist %2.2f, outward dist %2.2f, sigma %2.1f\n",
              vno, inward_dist, outward_dist, current_sigma) ;
    }

    v->val2 = current_sigma ;

    // if they can't be found, the vertex will be declared as missing
    //      v->val = 0.5*(wm_mean[0] + gm_mean[0]); //initialize
    //   v->valbak = 0.5*(wm_mean[1] + gm_mean[1]); //initialize

    /*
    search outwards and inwards and find the local gradient maximum
    at a location with a reasonable MR intensity value. This will
    be the location of the edge.
    */

    /* search in the normal direction to find the min value */
    max_mag_val30 = -10.0f ;
    mag = 5.0f ; //is 5 too high?
    max_mag = 0.0f ;
    min_val30 = 10000.0 ;
    max_mag_val5 = -10.0f;
    min_val5 = 10000.0;
    min_val_dist = 0.0f ;
    local_max_found = 0 ;
    for (dist = inward_dist ; dist <= outward_dist ; dist += STEP_SIZE) {
      x = v->x + v->nx*(dist-STEP_SIZE) ;
      y = v->y + v->ny*(dist-STEP_SIZE) ;
      z = v->z + v->nz*(dist-STEP_SIZE) ;
      MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
      MRIsampleVolume(mri_30, xw, yw, zw, &previous_val30) ;

      /* the previous point was inside the surface */
      if (previous_val30 < (wm_mean[0] + 3*wm_std[0]) && 
          previous_val30 > gm_mean[0]) {
        /* see if we are at a local maximum in the gradient magnitude */
        x = v->x + v->nx*dist ;
        y = v->y + v->ny*dist ;
        z = v->z + v->nz*dist ;
        MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolume(mri_30, xw, yw, zw, &val30) ;
        MRIsampleVolume(mri_5, xw, yw, zw, &val5) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, 
                                       nx, ny, nz,&mag30,sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, 
                                       nx, ny, nz,&mag5,sigma);

        x = v->x + v->nx*(dist+STEP_SIZE) ;
        y = v->y + v->ny*(dist+STEP_SIZE) ;
        z = v->z + v->nz*(dist+STEP_SIZE) ;
        MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, nx, ny, nz,
                                       &next_mag30, sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, nx, ny, nz,
                                       &next_mag5, sigma);

        x = v->x + v->nx*(dist-STEP_SIZE) ;
        y = v->y + v->ny*(dist-STEP_SIZE) ;
        z = v->z + v->nz*(dist-STEP_SIZE) ;
        MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, nx, ny, nz,
                                       &previous_mag30, sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, nx, ny, nz,
                                       &previous_mag5, sigma);


        //flash30 and flash5 have oppisite contrast at gray/white boundary
        // the weights are arbitrary for now,
        mag = mag30*0.7 - mag5*0.3;
        previous_mag = previous_mag30*0.7 - previous_mag5*0.3;
        next_mag = next_mag30*0.7 - next_mag5*0.3;

        if (val30 < min_val30) {
          min_val30 = val30 ;  /* used if no gradient max is found */
          min_val5 = val5;
          min_val_dist = dist ;
        }


        if (vno == Gdiag_no)
          fprintf(fp, "%2.3f  %2.3f  %2.3f  %2.3f  %2.3f\n",
                  dist, val30, mag, previous_mag, next_mag) ;

        /*
        if no local max has been found, or this one has a greater magnitude,
        and it is in the right intensity range....
        */
        if (
          (fabs(mag) > fabs(previous_mag)) &&
          (fabs(mag) > fabs(next_mag)) &&
          (val30 <= (wm_mean[0]+wm_std[0])) && 
          (val30 >= gm_mean[0]) //is this too restrictive??
        ) {
          x = v->x + v->nx*(dist+1) ;
          y = v->y + v->ny*(dist+1) ;
          z = v->z + v->nz*(dist+1) ;

          MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
          MRIsampleVolume(mri_30, xw, yw, zw, &next_val30) ;
          /*
            if next val is in the right range, and the intensity at
            this local max is less than the one at the previous local
            max, assume it is the correct one.
          */
          if ((next_val30 >= (gm_mean[0] - gm_std[0])) &&
              (next_val30 <= (wm_mean[0] + wm_std[0]))  &&
              (!local_max_found || (max_mag < fabs(mag)))) {
            local_max_found = 1 ;
            max_mag_dist = dist ;
            max_mag = fabs(mag) ;
            max_mag_val30 = val30 ;
            max_mag_val5 = val5 ;
          }
        } else {
          /*
            if no local max found yet, just used largest gradient
            if the intensity is in the right range.
          */
          if ((local_max_found == 0) &&
              (fabs(mag) > max_mag) &&
              (val30 <= (wm_mean[0] + wm_std[0])) &&
              (val30 >= gm_mean[0])
             ) {
            x = v->x + v->nx*(dist+1) ;
            y = v->y + v->ny*(dist+1) ;
            z = v->z + v->nz*(dist+1) ;
            MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
            MRIsampleVolume(mri_30, xw, yw, zw, &next_val30) ;
            // MRIsampleVolume(mri_5, xw, yw, zw, &next_val5) ;
            if ((next_val30 >= (gm_mean[0] - gm_std[0])) &&
                (next_val30 <= (wm_mean[0] + wm_std[0]))
               ) {
              max_mag_dist = dist ;
              max_mag = fabs(mag) ;
              max_mag_val30 = val30 ;
              max_mag_val5 = val5 ;
            }
          }
        }
      }
    }

    if (vno == Gdiag_no)
      fclose(fp) ;

    if (max_mag_val30 > 0)   /* found the border value */
    {
      if (local_max_found)
        ngrad_max++ ;
      else
        ngrad++ ;
      if (max_mag_dist > 0) {
        nout++ ;
        nfound++ ;
        mean_out += max_mag_dist ;
      } else {
        nin++ ;
        nfound++ ;
        mean_in -= max_mag_dist ;
      }

      mean_dist += max_mag_dist ;
      v->val = max_mag_val30;
      v->valbak = max_mag_val5;
      v->mean = max_mag ;
      mean_border += max_mag_val30;
      total_vertices++ ;
      v->d = max_mag_dist ;
      v->marked = 1 ;
    } else         /* couldn't find the border value */
    {
      if (min_val30 < 1000) {
        nmin++ ;
        v->d = min_val_dist ;
        if (min_val30 < (gm_mean[0] - gm_std[0]))
          min_val30 = gm_mean[0] - gm_std[0];

        //note for channel 1 (flash5), wm has lower intensity
        if (min_val5 < (wm_mean[1] - wm_std[1]))
          min_val5 = wm_mean[1] - wm_std[1];

        v->val = min_val30;
        v->valbak = min_val5;
        mean_border += min_val30 ;
        total_vertices++ ;
        v->marked = 1 ;
      } else {
        /* don't overwrite old target intensity if it was there */
        /*        v->val = -1.0f ;*/
        v->d = 0 ;
        if (v->val < 0) {
          nalways_missing++ ;
          v->marked = 0 ;
        } else
          v->marked = 1 ;
        nmissing++ ;
      }
    }
    if (vno == Gdiag_no)
      fprintf
        (stdout,
         "v %d, target value = %2.1f, mag = %2.1f, dist = %2.2f, %s\n",
         Gdiag_no, v->val, v->mean, v->d,
         local_max_found ? "local max" : max_mag_val30 > 0 ? "grad":"min");
  }

  mean_dist /= (float)(total_vertices-nmissing) ;
  mean_border /= (float)total_vertices ;

  if (nin > 0)
    mean_in /= (float)nin ;
  if (nout > 0)
    mean_out /= (float)nout ;


  fprintf(stdout,
          "mean border=%2.1f, %d (%d) missing vertices, mean dist %2.1f "
          "[%2.1f (%%%2.1f)->%2.1f (%%%2.1f))]\n",
          mean_border, nmissing, nalways_missing, mean_dist,
          mean_in, 100.0f*(float)nin/(float)nfound,
          mean_out, 100.0f*(float)nout/(float)nfound) ;
  fprintf(stdout, "%%%2.0f local maxima, %%%2.0f large gradients "
          "and %%%2.0f min vals, %d gradients ignored\n",
          100.0f*(float)ngrad_max/(float)mris->nvertices,
          100.0f*(float)ngrad/(float)mris->nvertices,
          100.0f*(float)nmin/(float)mris->nvertices, num_changed) ;
  if (log_fp) {
    fprintf(log_fp,
            "mean border=%2.1f, %d (%d) missing vertices, mean dist %2.1f "
            "[%2.1f (%%%2.1f)->%2.1f (%%%2.1f))]\n",
            mean_border, nmissing, nalways_missing, mean_dist,
            mean_in, 100.0f*(float)nin/(float)nfound,
            mean_out, 100.0f*(float)nout/(float)nfound) ;
    fprintf(log_fp, "%%%2.0f local maxima, %%%2.0f large gradients "
            "and %%%2.0f min vals, %d gradients ignored\n",
            100.0f*(float)ngrad_max/(float)mris->nvertices,
            100.0f*(float)ngrad/(float)mris->nvertices,
            100.0f*(float)nmin/(float)mris->nvertices, num_changed) ;
  }
  return(NO_ERROR) ;
}


#define min_gray_em_combined 50

int
MRIScomputeBorderValues_MEF_PIAL(MRI_SURFACE *mris, 
                                 MRI *mri_em_combined, MRI *mri_30,
                                 MRI *mri_5, float wm_mean[2], float wm_std[2],
                                 float gm_mean[2], float gm_std[2],
                                 double sigma,
                                 float max_thickness, FILE *log_fp) 
{
  //for pial surface, dura is a problem if I still use 
  // original images, really need to use the membership functions
  //or equivalently, the EM_combined
  //dura may be still low at both flash30 and flash5, but not that low
  double  val,  val30, val5, x, y, z, max_mag_val30, 
    max_mag_val5, xw, yw, zw,mag,max_mag, max_mag_dist=0.0f,
    min_val30, min_val5, inward_dist,outward_dist,xw1,yw1,zw1,
      min_val_dist, orig_dist, dx, dy, dz;
  double previous_mag, next_mag, previous_mag30, 
    previous_mag5, next_mag30, next_mag5;
  double mag30, mag5, previous_val30, previous_val5, next_val30, next_val5;
  int     total_vertices, vno, 
    nmissing = 0, nout = 0, nin = 0, nfound = 0,
    nalways_missing = 0, local_max_found, 
    ngrad_max, ngrad, nmin, num_changed=0 ;
  float   mean_border, mean_in, mean_out, 
    dist, nx, ny, nz, mean_dist, step_size ;
  double  current_sigma ;
  VERTEX  *v ;
  FILE    *fp = NULL ;

  step_size = mri_30->xsize/2 ;

  mean_dist = mean_in = mean_out = mean_border = 0.0f ;
  ngrad_max = ngrad = nmin = 0 ;
  MRISclearMarks(mris) ;  /* for soap bubble smoothing later */
  for (total_vertices = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ;
    y = v->y ;
    z = v->z ;
    MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
    x = v->x + v->nx ;
    y = v->y + v->ny ;
    z = v->z + v->nz ;
    MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw1, &yw1, &zw1) ;
    nx = xw1 - xw ;
    ny = yw1 - yw ;
    nz = zw1 - zw ;
    dist = sqrt(SQR(nx)+SQR(ny)+SQR(nz)) ;
    if (FZERO(dist))
      dist = 1 ;
    nx /= dist ;
    ny /= dist ;
    nz /= dist ;


    /*
    find the distance in the directions parallel and anti-parallel to
    the surface normal in which the gradient is pointing 'inwards'.
    The border will then be constrained to be within that region.
    */
    inward_dist = 1.0 ;
    outward_dist = -1.0 ;
    for (current_sigma = sigma; 
         current_sigma <= 10*sigma; 
         current_sigma *= 2) {
      for (dist = 0 ; dist > -max_thickness ; dist -= step_size) {
        dx = v->x-v->origx ;
        dy = v->y-v->origy ;
        dz = v->z-v->origz ;
        orig_dist = fabs(dx*v->nx + dy*v->ny + dz*v->nz) ;
        if (fabs(dist)+orig_dist > max_thickness)
          break ;
        x = v->x + v->nx*dist ;
        y = v->y + v->ny*dist ;
        z = v->z + v->nz*dist ;
        MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, 
                                       nx, ny,nz,&mag30,  current_sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, 
                                       nx, ny,nz,&mag5,  current_sigma);

        //The inside check can be restrictive, since pial surface have 
        // to move outside anyway. No not necessary for longitudinal method
        if (mag30 >= 0.0 || mag5 > 0.0)
          break ;
        MRIsampleVolume(mri_30, xw, yw, zw, &val30) ;
        MRIsampleVolume(mri_5, xw, yw, zw, &val5) ;

        if (val5 > (gm_mean[1] - 0.5*gm_std[1]) || 
            val30 > (wm_mean[0] - 0.5*wm_std[0]))
          break ;
      }
      inward_dist = dist+step_size/2 ;
      for (dist = 0 ; dist < max_thickness ; dist += step_size) {
        dx = v->x-v->origx ;
        dy = v->y-v->origy ;
        dz = v->z-v->origz ;
        orig_dist = fabs(dx*v->nx + dy*v->ny + dz*v->nz) ;
        if (fabs(dist)+orig_dist > max_thickness)
          break ;
        x = v->x + v->nx*dist ;
        y = v->y + v->ny*dist ;
        z = v->z + v->nz*dist ;
        MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, 
                                       nx, ny,nz, &mag30, current_sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, 
                                       nx, ny,nz, &mag5, current_sigma);

        if (mag30 >= 0.0)
          break ;
        MRIsampleVolume(mri_30, xw, yw, zw, &val30) ;
        MRIsampleVolume(mri_em_combined, xw, yw, zw, &val5) ; //borrowed
        if (val30 < (gm_mean[0] - 3*gm_std[0]) || val5 < 40)
          break ;
      }
      outward_dist = dist;
      if (!std::isfinite(outward_dist))
        DiagBreak() ;
      if (outward_dist >= (0.5*step_size))
        break ;
    }

    if (inward_dist > 0 && outward_dist <= 0)
      current_sigma = sigma ;  /* couldn't find anything */

    if (vno == Gdiag_no) {
      char fname[STRLEN] ;
      sprintf(fname, "v%d.%2.0f.log", Gdiag_no, sigma*100) ;
      fp = fopen(fname, "w") ;
      fprintf(stdout,
              "v %d: inward dist %2.2f, outward dist %2.2f, sigma %2.1f\n",
              vno, inward_dist, outward_dist, current_sigma) ;
    }

    v->val2 = current_sigma ;

    // if they can't be found, the vertex will be declared as missing
    //      v->val = 0.5*(wm_mean[0] + gm_mean[0]); //initialize
    //   v->valbak = 0.5*(wm_mean[1] + gm_mean[1]); //initialize

    /*
    search outwards and inwards and find the local gradient maximum
    at a location with a reasonable MR intensity value. This will
    be the location of the edge.
    */

    /* search in the normal direction to find the min value */
    max_mag_val30 = -10.0f ;
    mag = 0.0f;
    max_mag = 0.0f ;
    min_val30 = 10000.0 ;
    max_mag_val5 = -10.0f;
    min_val5 = 10000.0;
    min_val_dist = 0.0f ;
    local_max_found = 0 ;
    for (dist = inward_dist ; dist <= outward_dist ; dist += STEP_SIZE) {
      x = v->x + v->nx*(dist-STEP_SIZE) ;
      y = v->y + v->ny*(dist-STEP_SIZE) ;
      z = v->z + v->nz*(dist-STEP_SIZE) ;
      MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
      MRIsampleVolume(mri_30, xw, yw, zw, &previous_val30) ;
      MRIsampleVolume(mri_em_combined, xw, yw, zw, &previous_val5) ; //borrowed

      /* the previous point was inside the surface */
      if (previous_val30 < (wm_mean[0] - wm_std[0]) &&  previous_val5 > 50) {
        /* see if we are at a local maximum in the gradient magnitude */
        x = v->x + v->nx*dist ;
        y = v->y + v->ny*dist ;
        z = v->z + v->nz*dist ;
        MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolume(mri_30, xw, yw, zw, &val30) ;
        MRIsampleVolume(mri_5, xw, yw, zw, &val5) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, 
                                       nx, ny, nz,&mag30,sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, 
                                       nx, ny, nz,&mag5,sigma);

        x = v->x + v->nx*(dist+STEP_SIZE) ;
        y = v->y + v->ny*(dist+STEP_SIZE) ;
        z = v->z + v->nz*(dist+STEP_SIZE) ;
        MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, nx, ny, nz,
                                       &next_mag30, sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, nx, ny, nz,
                                       &next_mag5, sigma);

        x = v->x + v->nx*(dist-STEP_SIZE) ;
        y = v->y + v->ny*(dist-STEP_SIZE) ;
        z = v->z + v->nz*(dist-STEP_SIZE) ;
        MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, nx, ny, nz,
                                       &previous_mag30, sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, nx, ny, nz,
                                       &previous_mag5, sigma);


        //flash30 and flash5 have same contrast at gray/csf boundary
        // the weights are arbitrary for now,
        mag = mag30*0.8 + mag5*0.2;
        previous_mag = previous_mag30*0.8 + previous_mag5*0.2;
        next_mag = next_mag30*0.8 + next_mag5*0.2;

        if (val30 < min_val30) {
          min_val30 = val30 ;  /* used if no gradient max is found */
          min_val5 = val5;
          min_val_dist = dist ;
        }


        //       if (which == GRAY_CSF)
        {
          /*
             sample the next val we would process. If it is too low, then we
             have definitely reached the border, and the current gradient
             should be considered a local max.
          */
          x = v->x + v->nx*(dist+STEP_SIZE) ;
          y = v->y + v->ny*(dist+STEP_SIZE) ;
          z = v->z + v->nz*(dist+STEP_SIZE) ;
          MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
          MRIsampleVolume(mri_em_combined, xw, yw, zw, &next_val5) ; //borrowed
          if (next_val5 < 40)
            next_mag = 0 ;
        }

        if (vno == Gdiag_no)
          fprintf(fp, "%2.3f  %2.3f  %2.3f  %2.3f  %2.3f\n",
                  dist, val, mag, previous_mag, next_mag) ;

        /*
        if no local max has been found, or this one has a greater magnitude,
        and it is in the right intensity range....
        */
        if (
          (fabs(mag) > fabs(previous_mag)) &&
          (fabs(mag) > fabs(next_mag)) &&
          (val30 <= (gm_mean[0]-0.5*gm_std[0])) && 
          (val30 >= (gm_mean[0] - 3*gm_std[0])) //is this too restrictive??
        ) {
          x = v->x + v->nx*(dist+1) ;
          y = v->y + v->ny*(dist+1) ;
          z = v->z + v->nz*(dist+1) ;

          MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
          MRIsampleVolume(mri_30, xw, yw, zw, &next_val30) ;

          /*
            if next val is in the right range, and the intensity at
            this local max is less than the one at the previous local
            max, assume it is the correct one.
          */
          if ((next_val30 <= (gm_mean[0] - gm_std[0])) &&
              (!local_max_found || (max_mag < fabs(mag)))) {
            local_max_found = 1 ;
            max_mag_dist = dist ;
            max_mag = fabs(mag) ;
            max_mag_val30 = val30 ;
            max_mag_val5 = val5 ;
          }
        } else {
          /*
            if no local max found yet, just used largest gradient
            if the intensity is in the right range.
          */
          if ((local_max_found == 0) &&
              (fabs(mag) > max_mag) &&
              (val30 <= (gm_mean[0] - 0.5*gm_std[0])) &&
              (val30 >= (gm_mean[0] - 3*gm_std[0]))
             ) {
            x = v->x + v->nx*(dist+1) ;
            y = v->y + v->ny*(dist+1) ;
            z = v->z + v->nz*(dist+1) ;
            MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
            MRIsampleVolume(mri_30, xw, yw, zw, &next_val30) ;
            // MRIsampleVolume(mri_5, xw, yw, zw, &next_val5) ;
            if ((next_val30 <= (gm_mean[0] - gm_std[0]))
               ) {
              max_mag_dist = dist ;
              max_mag = fabs(mag) ;
              max_mag_val30 = val30 ;
              max_mag_val5 = val5 ;
            }
          }
        }

      } //end of if(previous_val30 ...)
    }

    if (vno == Gdiag_no)
      fclose(fp) ;

    //      if (which == GRAY_CSF && local_max_found == 0 && max_mag_dist > 0)
    {
      float outlen ;
      int   allgray = 1 ;

      /* check to make sure it's not ringing near the gray white boundary,
         by seeing if there is uniform stuff outside that could be gray matter.
      */
      for (outlen = max_mag_dist ; 
           outlen < max_mag_dist+2 ; 
           outlen += STEP_SIZE) {
        x = v->x + v->nx*outlen ;
        y = v->y + v->ny*outlen ;
        z = v->z + v->nz*outlen ;
        MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolume(mri_em_combined, xw, yw, zw, &val) ; //borrowed
        if (val < 50) {
          allgray = 0 ;
          break ;
        }
      }
      if (allgray) {
        if (Gdiag_no == vno)
          printf("v %d: exterior gray matter detected, "
                 "ignoring large gradient at %2.3f (I=%2.1f)\n",
                 vno, max_mag_dist, max_mag_val30) ;
        max_mag_val30 = -10 ;   /* don't worry about largest gradient */
        max_mag_dist = 0 ;
        num_changed++ ;
      }
    }

    if (max_mag_val30 > 0)   /* found the border value */
    {
      if (local_max_found)
        ngrad_max++ ;
      else
        ngrad++ ;
      if (max_mag_dist > 0) {
        nout++ ;
        nfound++ ;
        mean_out += max_mag_dist ;
      } else {
        nin++ ;
        nfound++ ;
        mean_in -= max_mag_dist ;
      }

      mean_dist += max_mag_dist ;
      v->val = max_mag_val30;
      v->valbak = max_mag_val5;
      v->mean = max_mag ;
      mean_border += max_mag_val30;
      total_vertices++ ;
      v->d = max_mag_dist ;
      v->marked = 1 ;
    } else         /* couldn't find the border value */
    {
      if (min_val30 < 1000) {
        nmin++ ;
        v->d = min_val_dist ;
        if (min_val30 < (gm_mean[0] - 3*gm_std[0]))
          min_val30 = gm_mean[0] - 3*gm_std[0];

        //note that flash5 may make the surface trying to move inwards!!
        //so for pial surface movement, 
        // may better only use flash30 for intensity term
        if (min_val5 < (gm_mean[1] - 3*gm_std[1]))
          min_val5 = gm_mean[1] - 3*gm_std[1];

        v->val = min_val30;
        v->valbak = min_val5;
        mean_border += min_val30 ;
        total_vertices++ ;
        v->marked = 1 ;
      } else {
        /* don't overwrite old target intensity if it was there */
        /*        v->val = -1.0f ;*/
        v->d = 0 ;
        if (v->val < 0) {
          nalways_missing++ ;
          v->marked = 0 ;
        } else
          v->marked = 1 ;
        nmissing++ ;
      }
    }

  }

  mean_dist /= (float)(total_vertices-nmissing) ;
  mean_border /= (float)total_vertices ;

  if (nin > 0)
    mean_in /= (float)nin ;
  if (nout > 0)
    mean_out /= (float)nout ;


  fprintf(stdout,
          "mean border=%2.1f, %d (%d) missing vertices, mean dist %2.1f "
          "[%2.1f (%%%2.1f)->%2.1f (%%%2.1f))]\n",
          mean_border, nmissing, nalways_missing, mean_dist,
          mean_in, 100.0f*(float)nin/(float)nfound,
          mean_out, 100.0f*(float)nout/(float)nfound) ;
  fprintf(stdout, "%%%2.0f local maxima, %%%2.0f large gradients "
          "and %%%2.0f min vals, %d gradients ignored\n",
          100.0f*(float)ngrad_max/(float)mris->nvertices,
          100.0f*(float)ngrad/(float)mris->nvertices,
          100.0f*(float)nmin/(float)mris->nvertices, num_changed) ;
  if (log_fp) {
    fprintf(log_fp,
            "mean border=%2.1f, %d (%d) missing vertices, mean dist %2.1f "
            "[%2.1f (%%%2.1f)->%2.1f (%%%2.1f))]\n",
            mean_border, nmissing, nalways_missing, mean_dist,
            mean_in, 100.0f*(float)nin/(float)nfound,
            mean_out, 100.0f*(float)nout/(float)nfound) ;
    fprintf(log_fp, "%%%2.0f local maxima, %%%2.0f large gradients "
            "and %%%2.0f min vals, %d gradients ignored\n",
            100.0f*(float)ngrad_max/(float)mris->nvertices,
            100.0f*(float)ngrad/(float)mris->nvertices,
            100.0f*(float)nmin/(float)mris->nvertices, num_changed) ;
  }
  return(NO_ERROR) ;
}

#if 0
#define MAX_ASYNCH_MM       0.3
#define REDUCTION_PCT      0.5
#define MAX_REDUCTIONS     2

int
MRISpositionSurface_mef(MRI_SURFACE *mris, MRI *mri_30, MRI *mri_5,
                        INTEGRATION_PARMS *parms, 
                        float weight30, float weight5) {
  /*  char   *cp ;*/
  int    avgs, niterations, n, write_iterations, nreductions = 0, done ;
  double delta_t = 0.0, rms, dt, l_intensity, base_dt, last_rms, max_mm;
  MHT    *mht = NULL, *mht_v_orig = NULL, *mht_v_current = NULL ;
  Timer then ;
  int msec ;

  max_mm = MIN(MAX_ASYNCH_MM, MIN(mri_30->xsize, 
                                  MIN(mri_30->ysize, mri_30->zsize))/2) ;

  //note that the following is for pial surface avoid intersection with white
  if (!FZERO(parms->l_surf_repulse))
    mht_v_orig = MHTcreateVertexTable(mris, ORIGINAL_VERTICES) ;

  base_dt = parms->dt ;
  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;
  then.reset() ;
  //the following are used in mrisComputeIntensityError() and computeSSE()
  parms->mri_brain = NULL; //mri_30 ;
  parms->mri_smooth = NULL; //mri_5 ;
  niterations = parms->niterations ; //should be different for white and pial; 
                                     // yeah 25 for white and 30 for pial
  write_iterations = parms->write_iterations ;
  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN] ;

    if (!parms->fp) {
      sprintf(fname, "%s.%s.out",
              mris->hemisphere==RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
      if (!parms->start_t)
        parms->fp = fopen(fname, "w") ;
      else
        parms->fp = fopen(fname, "a") ;
      if (!parms->fp)
        ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
                  Progname, fname) ;
    }
    mrisLogIntegrationParms(parms->fp, mris, parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms(stderr, mris, parms) ;

  mrisClearMomentum(mris) ;
  MRIScomputeMetricProperties(mris) ;
  MRISstoreMetricProperties(mris) ;

  MRIScomputeNormals(mris) ;
  MRISclearD(mris) ;

  MRISclearCurvature(mris) ;  /* curvature will be used to calculate sulc */

  /* write out initial surface */
  if ((parms->write_iterations > 0) && (Gdiag&DIAG_WRITE) && !parms->start_t)
    mrisWriteSnapshot(mris, parms, 0) ;

  avgs = parms->n_averages ;
  last_rms = rms = mrisRmsValError_mef(mris, mri_30, mri_5, weight30, weight5);
  // last_sse = sse = MRIScomputeSSE(mris, parms) ; 
//this computation results were never used

  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "%3.3d: dt: %2.4f, rms=%2.2f\n",
            0, 0.0f, (float)rms);

  if (Gdiag & DIAG_WRITE) {
    fprintf(parms->fp, "%3.3d: dt: %2.4f, rms=%2.2f\n",
            0, 0.0f,  (float)rms);
    fflush(parms->fp) ;
  }

  dt = parms->dt ;
  l_intensity = parms->l_intensity ;
  for (n = parms->start_t ; n < parms->start_t+niterations ; n++) {
    if (!FZERO(parms->l_repulse)) {
      MHTfree(&mht_v_current);
      mht_v_current = MHTcreateVertexTable(mris, CURRENT_VERTICES);
    }
    if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST)) {
      MHTfreeTable(&mht);
      mht = MHTcreateFaceTable(mris) ;
    }
    MRISclearGradient(mris) ;
    mrisComputeIntensityTerm_mef(mris, l_intensity, mri_30, mri_5,
                                 parms->sigma, weight30, weight5);
    //the following term is not used for white, but used for pial!
    mrisComputeSurfaceRepulsionTerm(mris, parms->l_surf_repulse, mht_v_orig);
#if 1
    mrisAverageSignedGradients(mris, avgs) ;
#else
    mrisAverageWeightedGradients(mris, avgs) ;
#endif
    /*  mrisUpdateSulcalGradients(mris, parms) ;*/

    /* smoothness terms */
    mrisComputeSpringTerm(mris, parms->l_spring) ;
    mrisComputeNormalizedSpringTerm(mris, parms->l_spring_norm) ;
    mrisComputeRepulsiveTerm(mris, parms->l_repulse, mht_v_current) ;
    mrisComputeThicknessSmoothnessTerm(mris, parms->l_tsmooth) ;
    mrisComputeNormalSpringTerm(mris, parms->l_nspring) ;
    mrisComputeQuadraticCurvatureTerm(mris, parms->l_curv) ;
    mrisComputeTangentialSpringTerm(mris, parms->l_tspring) ;


    do {
      MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
      delta_t = mrisAsynchronousTimeStep(mris, parms->momentum, dt,mht,
                                         max_mm) ;
      if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
        MHTcheckFaces(mris, mht) ;
      MRIScomputeMetricProperties(mris) ;
      rms = mrisRmsValError_mef(mris, mri_30, mri_5, weight30, weight5) ;
      //   sse = MRIScomputeSSE(mris, parms) ;
      done = 1 ;
      if (rms > last_rms-0.05)  /* error increased - reduce step size */
      {
        nreductions++ ;
        parms->dt *= REDUCTION_PCT ;
        dt = parms->dt ;
        fprintf(stdout,
                "rms = %2.2f, time step reduction %d of %d to %2.3f...\n",
                rms, nreductions, MAX_REDUCTIONS+1, dt) ;
        mrisClearMomentum(mris) ;
        if (rms > last_rms)  /* error increased - reject step */
        {
          MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
          MRIScomputeMetricProperties(mris) ;

          /* if error increased and we've only reduced the time
             step a few times, try taking a smaller step (done=0).
          */
          done = (nreductions > MAX_REDUCTIONS) ;
        }
      }
    } while (!done) ;
    //last_sse = sse ;
    last_rms = rms ;

    mrisTrackTotalDistanceNew(mris) ;  /* computes signed deformation amount */
    rms = mrisRmsValError_mef(mris, mri_30, mri_5, weight30, weight5) ;
    //  sse = MRIScomputeSSE(mris, parms) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "%3.3d: dt: %2.4f, rms=%2.2f\n",
              n+1,(float)delta_t,  (float)rms);

    if (Gdiag & DIAG_WRITE) {
      fprintf(parms->fp, "%3.3d: dt: %2.4f, rms=%2.2f\n",
              n+1,(float)delta_t, (float)rms);
      fflush(parms->fp) ;
    }

    if ((parms->write_iterations > 0) &&
        !((n+1)%write_iterations)&&(Gdiag&DIAG_WRITE))
      mrisWriteSnapshot(mris, parms, n+1) ;

    if ((Gdiag & DIAG_SHOW) && !((n+1)%5))
      MRISprintTessellationStats(mris, stderr) ;
    if (nreductions > MAX_REDUCTIONS) {
      n++ ;  /* count this step */
      break ;
    }
  }

  parms->start_t = n ;
  parms->dt = base_dt ;
  if (Gdiag & DIAG_SHOW) {
    msec = then.milliseconds() ;
    fprintf(stdout,"positioning took %2.1f minutes\n",
            (float)msec/(60*1000.0f));
  }
  if (Gdiag & DIAG_WRITE) {
    fclose(parms->fp) ;
    parms->fp = NULL ;
  }

  /*  MHTcheckSurface(mris, mht) ;*/
  if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST))
    MHTfree(&mht) ;
  if (mht_v_current)
    MHTfree(&mht_v_current) ;
  if (mht_v_orig)
    MHTfree(&mht_v_orig) ;
  return(NO_ERROR) ;
}
#endif

#if 0
static double
mrisRmsValError_mef(MRI_SURFACE *mris,
                    MRI *mri_30, MRI *mri_5, float weight30, float weight5) {
  int     vno, n; //, xv, yv, zv ;
  double  val30, val5, total, delta, x, y, z ;
  VERTEX  *v ;

  for (total = 0.0, n = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < 0)
      continue ;
    n++ ;
    MRISvertexToVoxel(v, mri_30, &x, &y, &z) ;
    //      xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    MRIsampleVolume(mri_30, x, y, z, &val30) ;
    MRIsampleVolume(mri_5, x, y, z, &val5) ;
    delta = (val30 - v->val) ;
    total += delta*delta*weight30;
    delta = (val5 - v->valbak) ;
    total += delta*delta*weight5;
  }
  return(sqrt(total / (double)n)) ;
}

#endif

#if 0
static int
mrisComputeIntensityTerm_mef(MRI_SURFACE *mris, 
                             double l_intensity, MRI *mri_30,
                             MRI *mri_5, 
                             double sigma_global, 
                             float weight30, float weight5) 
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z, nx, ny, nz, dx, dy, dz ;
  double  val0, xw, yw, zw, del, val_outside, val_inside, delI, delV, k,
  ktotal ;
  double  sigma ;

  if (FZERO(l_intensity))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->val < 0)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ;
    y = v->y ;
    z = v->z ;

    MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_30, xw, yw, zw, &val0) ;
    sigma = v->val2 ;
    if (FZERO(sigma))
      sigma = sigma_global ;

    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;

    /* compute intensity gradient using smoothed volume */

    {
      double dist, val, step_size ;
      int  n ;

      step_size = MIN(sigma/2, 
                      MIN(mri_30->xsize, 
                          MIN(mri_30->ysize, mri_30->zsize))*0.5) ;
      ktotal = 0.0 ;
      for (n = 0, val_outside = val_inside = 0.0, dist = step_size ;
           dist <= 2*sigma;
           dist += step_size, n++) {
        k = exp(-dist*dist/(2*sigma*sigma)) ;
        ktotal += k ;
        xw = x + dist*nx ;
        yw = y + dist*ny ;
        zw = z + dist*nz ;
        MRISsurfaceRASToVoxel(mris, mri_30, xw, yw, zw, &xw, &yw, &zw);
        MRIsampleVolume(mri_30, xw, yw, zw, &val) ;
        val_outside += k*val ;

        xw = x - dist*nx ;
        yw = y - dist*ny ;
        zw = z - dist*nz ;
        MRISsurfaceRASToVoxel(mris, mri_30, xw, yw, zw, &xw, &yw, &zw);
        MRIsampleVolume(mri_30, xw, yw, zw, &val) ;
        val_inside += k*val ;
      }
      val_inside /= (double)ktotal ;
      val_outside /= (double)ktotal ;
    }

    delV = v->val - val0 ;
    delI = (val_outside - val_inside) / 2.0 ;

    if (!FZERO(delI))
      delI /= fabs(delI) ;
    else
      delI = -1; //intensity tends to increase inwards for flash30


    if (delV > 5)
      delV = 5 ;
    else if (delV < -5)
      delV = -5 ;

    del = l_intensity * delV * delI ;

    dx = nx * del*weight30 ;
    dy = ny * del*weight30 ;
    dz = nz * del*weight30;

    v->dx += dx ;
    v->dy += dy ;
    v->dz += dz ;

    //now compute flash5
    /* compute intensity gradient using smoothed volume */
    MRISsurfaceRASToVoxel(mris, mri_5, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_5, xw, yw, zw, &val0) ;
    {
      double dist, val, step_size ;
      int  n ;

      step_size = MIN(sigma/2, 
                      MIN(mri_5->xsize, 
                          MIN(mri_5->ysize, mri_5->zsize))*0.5) ;
      ktotal = 0.0 ;
      for (n = 0, val_outside = val_inside = 0.0, dist = step_size ;
           dist <= 2*sigma;
           dist += step_size, n++) {
        k = exp(-dist*dist/(2*sigma*sigma)) ;
        ktotal += k ;
        xw = x + dist*nx ;
        yw = y + dist*ny ;
        zw = z + dist*nz ;
        MRISsurfaceRASToVoxel(mris, mri_5, xw, yw, zw, &xw, &yw, &zw);
        MRIsampleVolume(mri_5, xw, yw, zw, &val) ;
        val_outside += k*val ;

        xw = x - dist*nx ;
        yw = y - dist*ny ;
        zw = z - dist*nz ;
        MRISsurfaceRASToVoxel(mris, mri_5, xw, yw, zw, &xw, &yw, &zw);
        MRIsampleVolume(mri_5, xw, yw, zw, &val) ;
        val_inside += k*val ;
      }
      val_inside /= (double)ktotal ;
      val_outside /= (double)ktotal ;
    }

    delV = v->valbak - val0 ;
    delI = (val_outside - val_inside) / 2.0 ;

    if (!FZERO(delI))
      delI /= fabs(delI) ;
    else
      delI = 0 ;

    if (delV > 5)
      delV = 5 ;
    else if (delV < -5)
      delV = -5 ;

    del = l_intensity * delV * delI ;

    dx = nx * del*weight5 ;
    dy = ny * del*weight5;
    dz = nz * del*weight5;

    v->dx += dx ;
    v->dy += dy ;
    v->dz += dz ;
  }

  return(NO_ERROR) ;
}
#endif

#define STD_TARGET 6

static int  MRInormalizeMEF(MRI *mri_src, MRI *mri_em_seg) {
  /* Normalize the source volume to be 110 mean and 10 variance */
  /* mri_dst and mri_src can be the same */
  /* mri_em_seg should be 120 for WM, and 75 for GM, and 30 for CSF */

  int width, height, depth, x, y, z;
  double mean, variance, std, total, tmpval;
  double val1, val2;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  /* compute mean */
  mean = 0.0;
  total = 0.0;
  variance = 0.0;
  for (z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++) {
        val1 = MRIgetVoxVal(mri_em_seg, x, y, z,0);
        val2 = MRIgetVoxVal(mri_src, x, y, z, 0);
        if (val2 <= 1e-10 || val1 < 100) continue;

        mean += val2;
        variance += val2 * val2;
        total += 1;
      }

  if (total <= 1) {
    printf("Too few WM voxels. Something is wrong. Exit. \n");
    exit(0);
  }

  mean /= (double)total;
  variance /= (double)total;
  std = sqrt(variance - mean*mean);

  printf("mean = %g, std = %g\n", mean, std);

  /* normalization: invert std first to save time */
  // normalize to mean = 110 and std = 10
  std = 1.0/(std + 1e-30);
  for (z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++) {

        tmpval = MRIgetVoxVal(mri_src, x, y, z,0);

        tmpval = (tmpval - mean)*STD_TARGET*std + 110;
        MRIsetVoxVal(mri_src,x,y,z,0,tmpval);
      }

  return 0;
}


static int  MRIcomputeClassStatistics_mef(MRI *mri_T1_30, 
                                          MRI *mri_T1_5, 
                                          MRI *mri_em_seg, 
                                          float *white_mean, float *white_std, 
                                          float *gray_mean, float *gray_std) 
{
  //Do I include all WM and GM or just at boundaries?
  //when mri_em_seg is generated, WM mean is set to 120, GM is 75, CSF is 30
  int channel;
  float val30, val5;
  float val;
  double dist1, dist2, dist3, sumdist;
  int x, y,z, width, height, depth;
  double mem_wm, mem_gm;
  double sum_wm, sum_gm;
  int total;

  width = mri_T1_30->width;
  height =  mri_T1_30->height;
  depth = mri_T1_30->depth;

  white_mean[0] = 0;
  white_mean[1] = 0;
  gray_mean[0] = 0;
  gray_mean[1] = 0;
  white_std[0] = 0;
  white_std[1] = 0;
  gray_std[0] = 0;
  gray_std[1] = 0;
  sum_wm = 0;
  sum_gm = 0;
  total = 0;
  for (z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++) {
        val = MRIgetVoxVal(mri_em_seg, x, y, z,0);
        val30 = MRIgetVoxVal(mri_T1_30, x, y, z,0);
        val5 = MRIgetVoxVal(mri_T1_5, x, y, z,0);

        if (val30 <= 1e-10 || val5 <= 1e-10 || val < 40) continue;

        total++;

        dist1 =  val - 120;
        dist2 = val - 75;
        dist3 = val - 30;

        dist1 =  1.0/(dist1*dist1 + 1e-30);
        dist2 =  1.0/(dist2*dist2 + 1e-30);
        dist3 =  1.0/(dist3*dist3 + 1e-30);

        sumdist = dist1 + dist2 + dist3;

        mem_wm = dist1/sumdist;
        mem_gm = dist2/sumdist;

        white_mean[0] += mem_wm*val30;
        white_mean[1] += mem_wm*val5;

        white_std[0] += mem_wm*val30*val30;
        white_std[1] += mem_wm*val5*val5;

        sum_wm += mem_wm;

        gray_mean[0] += mem_gm*val30;
        gray_mean[1] += mem_gm*val5;

        gray_std[0] += mem_gm*val30*val30;
        gray_std[1] += mem_gm*val5*val5;

        sum_gm += mem_gm;
      }

  if (total < 10) {
    printf("Too few brain voxels. Sth is wrong. Exit.\n");
    exit(0);
  }

  for (channel = 0; channel <= 1; channel ++) {
    white_mean[channel] /= (sum_wm + 1e-30);
    white_std[channel] /= (sum_wm + 1e-30);
    gray_mean[channel] /= (sum_gm + 1e-30);
    gray_std[channel] /= (sum_gm + 1e-30);

    white_std[channel] = 
      sqrt(white_std[channel] - white_mean[channel]*white_mean[channel]);
    gray_std[channel] = 
      sqrt(gray_std[channel] - gray_mean[channel]*gray_mean[channel]);
  }

  return 0;
}



int MRISaverageMarkedValbaks(MRI_SURFACE *mris, int navgs) {
  int    i, vno, vnb, vnum ;
  float  val, num ;

  for (i = 0 ; i < navgs ; i++) {
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->marked == 0)
        continue ;
      val = v->valbak ;
      int const * pnb = vt->v ;
      vnum = vt->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag || v->marked == 0 )
          continue ;
        num++ ;
        val += vn->valbak ;
      }
      num++ ;  /*  account for central vertex */
      v->tdx = val / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      VERTEX * const v = &mris->vertices[vno] ;
      if (v->ripflag ||  v->marked == 0)
        continue ;
      v->valbak = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}
