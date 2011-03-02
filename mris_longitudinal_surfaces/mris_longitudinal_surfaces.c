/**
 * @file  mris_longitudinal_surfaces.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:33 $
 *    $Revision: 1.5 $
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
#include <sys/types.h>
#include <sys/stat.h>

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

static char vcid[] = "$Id: mris_longitudinal_surfaces.c,v 1.5 2011/03/02 00:04:33 nicks Exp $";

int main(int argc, char *argv[]) ;

#define BRIGHT_LABEL         130
#define BRIGHT_BORDER_LABEL  100

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  mrisFindMiddleOfGray(MRI_SURFACE *mris) ;
static MRI  *MRIsmoothMasking(MRI *mri_src, MRI *mri_mask, MRI *mri_dst,
                              int mask_val, int wsize) ;
MRI *MRIfillVentricle(MRI *mri_inv_lv, MRI *mri_T1, float thresh,
                      int out_label, MRI *mri_dst);

int MRISfindExpansionRegions(MRI_SURFACE *mris) ;
int MRIsmoothBrightWM(MRI *mri_T1, MRI *mri_wm) ;
MRI *MRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm) ;

static LABEL *highres_label = NULL ;
static char T1_name[STRLEN] = "brain" ;

static char *white_fname = NULL ;
static int use_mode = 1 ;

static char *orig_white = NULL ;
static char *orig_pial = NULL ;

char *Progname ;

static double std_scale = 1.0;

static int graymid = 0 ;
static int curvature_avgs = 10 ;
static int create = 1 ;
static int smoothwm = 0 ;
static int white_only = 0 ;
static int overlay = 0 ;
static int inverted_contrast = 0 ;
static int auto_detect_stats = 1 ;

static int in_out_in_flag = 0 ;  /* for Arthur (as are most things) */
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
static int nwhite = 25 /*5*/ ;
static int ngray = 30 /*45*/ ;

static int nowhite = 0 ;
static int nbrs = 2 ;
static int write_vals = 0 ;

static char *orig_name = ORIG_NAME ; // "orig"
static char *suffix = "" ;
static char *output_suffix = "" ;
static char *xform_fname = NULL ;

static char pial_name[STRLEN] = "pial" ;
static char white_matter_name[STRLEN] = WHITE_MATTER_NAME ; // "white"

static int lh_label = LH_LABEL ;
static int rh_label = RH_LABEL ;

static int max_pial_averages = 16 ;
static int min_pial_averages = 2 ;
static int max_white_averages = 4 ;
static int min_white_averages = 0 ;
static float pial_sigma = 2.0f ;
static float white_sigma = 2.0f ;
static float max_thickness = 5.0 ;


#define MAX_WHITE             120
#define MAX_BORDER_WHITE      105
#define MIN_BORDER_WHITE       85
#define MIN_GRAY_AT_WHITE_BORDER  70

#define MAX_GRAY               95
#define MIN_GRAY_AT_CSF_BORDER    40
#define MID_GRAY               ((max_gray + min_gray_at_csf_border) / 2)
#define MAX_GRAY_AT_CSF_BORDER    75
#define MIN_CSF                10
#define MAX_CSF                40

static  int   max_border_white_set = 0,
                                     min_border_white_set = 0,
                                                            min_gray_at_white_border_set = 0,
                                                                                           max_gray_set = 0,
                                                                                                          max_gray_at_csf_border_set = 0,
                                                                                                                                       min_gray_at_csf_border_set = 0,
                                                                                                                                                                    min_csf_set = 0,
                                                                                                                                                                                  max_csf_set = 0 ;

static  float   max_border_white = MAX_BORDER_WHITE,
                                   min_border_white = MIN_BORDER_WHITE,
                                                      min_gray_at_white_border = MIN_GRAY_AT_WHITE_BORDER,
                                                                                 max_gray = MAX_GRAY,
                                                                                            max_gray_at_csf_border = MAX_GRAY_AT_CSF_BORDER,
                                                                                                                     min_gray_at_csf_border = MIN_GRAY_AT_CSF_BORDER,
                                                                                                                                              min_csf = MIN_CSF,
                                                                                                                                                        max_csf = MAX_CSF ;
static char sdir[STRLEN] = "" ;

static int MGZ = 0; // for use with MGZ format

static int longitudinal = 0;

static float check_contrast_direction(MRI_SURFACE *mris,MRI *mri_T1) ;
int
main(int argc, char *argv[]) {
  char          **av, *hemi, *sname, *cp, fname[STRLEN], mdir[STRLEN];
  int           ac, nargs, i, label_val, replace_val, msec, n_averages, j ;
  MRI_SURFACE   *mris ;
  MRI           *mri_wm, *mri_kernel = NULL, *mri_smooth = NULL,
                                       *mri_filled, *mri_T1, *mri_labeled, *mri_T1_white = NULL, *mri_T1_pial ;
  float         max_len ;
  float         white_mean, white_std, gray_mean, gray_std ;
  double        l_intensity, current_sigma ;
  struct timeb  then ;
  M3D           *m3d ;

  char cmdline[CMD_LINE_LEN] ;

  make_cmd_version_string (argc, argv, "$Id: mris_longitudinal_surfaces.c,v 1.5 2011/03/02 00:04:33 nicks Exp $", "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_longitudinal_surfaces.c,v 1.5 2011/03/02 00:04:33 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  memset(&parms, 0, sizeof(parms)) ;
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

  ac = argc ;
  av = argv ;
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

  TimerStart(&then) ;
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

  // print out version of this program and mrisurf.c
  printf("%s\n",vcid);
  printf("%s\n",MRISurfSrcVersion());
  fflush(stdout);
  sprintf(fname, "%s/%s/surf/mris_make_surfaces.%s.mrisurf.c.version", sdir, sname, hemi) ;

  sprintf(fname, "%s/%s/mri/filled", sdir, sname) ;
  if (MGZ) sprintf(fname, "%s.mgz",fname);
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_filled = MRIread(fname) ;
  if (!mri_filled)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;
  ////////////////////////////// we can handle only conformed volumes
  setMRIforSurface(mri_filled);

  if (!stricmp(hemi, "lh")) {
    label_val = lh_label ;
    replace_val = rh_label ;
  } else {
    label_val = rh_label ;
    replace_val = lh_label ;
  }

  sprintf(fname, "%s/%s/mri/%s", sdir, sname, T1_name) ;
  if (MGZ) sprintf(fname, "%s.mgz",fname);
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_T1 = mri_T1_pial = MRIread(fname) ;

  if (!mri_T1)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;
  /////////////////////////////////////////
  setMRIforSurface(mri_T1);

  if (white_fname != NULL) {
    sprintf(fname, "%s/%s/mri/%s", sdir, sname, white_fname) ;
    if (MGZ) sprintf(fname, "%s.mgz",fname);
    fprintf(stderr, "reading volume %s...\n", fname) ;
    mri_T1_white = MRIread(fname) ;
    if (!mri_T1_white)
      ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
                Progname, fname) ;
    /////////////////////////////////////////
    if (mri_T1_white->type != MRI_UCHAR) {
      MRI *mri_tmp ;

      MRIeraseNegative(mri_T1_white, mri_T1_white) ;
      mri_tmp = MRIchangeType(mri_T1_white, MRI_UCHAR, 0, 255, 1) ;
      MRIfree(&mri_T1_white) ;
      mri_T1_white = mri_tmp ; // swap
    }

    setMRIforSurface(mri_T1_white);
  }
  if (apply_median_filter) {
    MRI *mri_tmp ;

    fprintf(stderr, "applying median filter to T1 image...\n") ;
    mri_tmp = MRImedian(mri_T1, NULL, 3, NULL) ;
    MRIfree(&mri_T1) ;
    mri_T1 = mri_tmp ;

    if (mri_T1_white) {
      mri_tmp = MRImedian(mri_T1_white, NULL, 3, NULL) ;
      MRIfree(&mri_T1_white) ;
      mri_T1_white = mri_tmp ; // swap
    }
  }

  if (xform_fname) {
    char fname[STRLEN], ventricle_fname[STRLEN] ;
    MRI  *mri_lv, *mri_inv_lv ;

    sprintf(fname, "%s/%s/mri/transforms/%s", sdir, sname,xform_fname) ;
    fprintf(stderr, "reading transform %s...\n", fname) ;
    m3d = MRI3DreadSmall(fname) ;
    if (!m3d)
      ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                Progname, fname) ;
    sprintf(ventricle_fname, "%s/average/%s_ventricle.mgz#0@mgh",
            mdir, !stricmp(hemi, "lh") ? "left" : "right") ;
    fprintf(stderr,"reading ventricle representation %s...\n",ventricle_fname);
    mri_lv = MRIread(ventricle_fname) ;
    if (!mri_lv)
      ErrorExit(ERROR_NOFILE,"%s: could not read %s",Progname,ventricle_fname);
    fprintf(stderr, "applying inverse morph to ventricle...\n") ;
    mri_inv_lv = MRIapplyInverse3DMorph(mri_lv,m3d,NULL);
    MRI3DmorphFree(&m3d) ;
    MRIfree(&mri_lv) ;
    fprintf(stderr, "filling in ventricle...\n") ;
    mri_lv = MRIfillVentricle(mri_inv_lv,mri_T1,100,
                              DEFAULT_DESIRED_WHITE_MATTER_VALUE, NULL);
    MRIfree(&mri_inv_lv) ;
    MRIunion(mri_lv, mri_T1, mri_T1) ;
    MRIfree(&mri_lv) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      sprintf(fname, "%s/%s/mri/T1_filled", sdir, sname) ;
      MRIwrite(mri_T1, fname) ;
    }
  }
  /* remove other hemi */
  MRIdilateLabel(mri_filled, mri_filled, replace_val, 1) ;
  if (replace_val == RH_LABEL) {
    MRIdilateLabel(mri_filled, mri_filled, RH_LABEL2, 1) ;
    MRImask(mri_T1, mri_filled, mri_T1, RH_LABEL2,0) ;
    if (mri_T1_white)
      MRImask(mri_T1_white, mri_filled, mri_T1_white, RH_LABEL2,0) ;
  }

  if (mri_T1_white)
    MRImask(mri_T1_white, mri_filled, mri_T1_white, replace_val,0) ;
  MRImask(mri_T1, mri_filled, mri_T1, replace_val,0) ;
  MRIfree(&mri_filled) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_T1, "r.mgz") ;

  sprintf(fname, "%s/%s/mri/wm", sdir, sname) ;
  if (MGZ) sprintf(fname, "%s.mgz",fname);
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_wm = MRIread(fname) ;
  if (!mri_wm)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;
  //////////////////////////////////////////
  setMRIforSurface(mri_wm);

  MRIsmoothBrightWM(mri_T1, mri_wm) ;
  mri_labeled = MRIfindBrightNonWM(mri_T1, mri_wm) ;
  if (mri_T1_white)
    MRIsmoothBrightWM(mri_T1_white, mri_wm) ;
  if (overlay) {
    fprintf(stderr, "overlaying editing into T1 volume...\n") ;
    MRImask(mri_T1, mri_wm, mri_T1,
            WM_EDITED_ON_VAL, DEFAULT_DESIRED_WHITE_MATTER_VALUE);
    MRIsmoothMasking(mri_T1, mri_wm, mri_T1, WM_EDITED_ON_VAL, 15) ;
    sprintf(fname, "%s/%s/mri/T1_overlay", sdir, sname) ;
    MRIwrite(mri_T1, fname) ;
  }


  sprintf(fname, "%s/%s/surf/%s.%s%s", sdir, sname, hemi, orig_name, suffix) ;
  fprintf(stderr, "reading original surface position from %s...\n", fname) ;
  mris = MRISreadOverAlloc(fname, 1.1) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  MRISaddCommandLine(mris, cmdline) ;

  if (auto_detect_stats) {
    MRI *mri_tmp ;
    float white_mode, gray_mode ;

    mri_tmp = MRIbinarize(mri_wm, NULL, WM_MIN_VAL, MRI_NOT_WHITE, MRI_WHITE) ;
    fprintf(stderr, "computing class statistics...\n");
    MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
    MRIScomputeClassModes(mris, mri_T1, &white_mode, &gray_mode, NULL);
    MRIcomputeClassStatistics(mri_T1, mri_tmp, 30, WHITE_MATTER_MEAN,
                              &white_mean, &white_std, &gray_mean,
                              &gray_std) ;
    if (use_mode) {
      printf("using class modes intead of means....\n") ;
      white_mean = white_mode ;
      gray_mean = gray_mode ;
    }

    white_std /= std_scale;
    gray_std /= std_scale;

    if (!min_gray_at_white_border_set)
      min_gray_at_white_border = gray_mean-gray_std ;
    if (!max_border_white_set)
      max_border_white = white_mean+white_std ;
    if (!max_csf_set)
      max_csf = gray_mean-2*gray_std ;
    if (!min_border_white_set)
      min_border_white = gray_mean ;
    fprintf(stderr, "setting MIN_GRAY_AT_WHITE_BORDER to %2.1f (was %d)\n",
            min_gray_at_white_border, MIN_GRAY_AT_WHITE_BORDER) ;
    fprintf(stderr, "setting MAX_BORDER_WHITE to %2.1f (was %d)\n",
            max_border_white, MAX_BORDER_WHITE) ;
    fprintf(stderr, "setting MIN_BORDER_WHITE to %2.1f (was %d)\n",
            min_border_white, MIN_BORDER_WHITE) ;
    fprintf(stderr, "setting MAX_CSF to %2.1f (was %d)\n",
            max_csf, MAX_CSF) ;

    if (!max_gray_set)
      max_gray = white_mean-white_std ;
    if (!max_gray_at_csf_border_set)
      max_gray_at_csf_border = gray_mean-0.5*gray_std ;
    if (!min_gray_at_csf_border_set)
      min_gray_at_csf_border = gray_mean - 3*gray_std ;
    fprintf(stderr, "setting MAX_GRAY to %2.1f (was %d)\n",
            max_gray, MAX_GRAY) ;
    fprintf(stderr, "setting MAX_GRAY_AT_CSF_BORDER to %2.1f (was %d)\n",
            max_gray_at_csf_border, MAX_GRAY_AT_CSF_BORDER) ;
    fprintf(stderr, "setting MIN_GRAY_AT_CSF_BORDER to %2.1f (was %d)\n",
            min_gray_at_csf_border, MIN_GRAY_AT_CSF_BORDER) ;
    MRIfree(&mri_tmp) ;
  }
  MRIfree(&mri_wm) ;
  inverted_contrast = (check_contrast_direction(mris,mri_T1) < 0) ;
  if (inverted_contrast) {
    printf("inverted contrast detected....\n") ;
  }
  if (highres_label) {
    LabelRipRestOfSurface(highres_label, mris) ;
  }
  if (smooth && !nowhite) {
    printf("smoothing surface for %d iterations...\n", smooth) ;
    MRISaverageVertexPositions(mris, smooth) ;
  }

  if (nbrs > 1)
    MRISsetNeighborhoodSize(mris, nbrs) ;

  sprintf(parms.base_name, "%s%s%s", white_matter_name, output_suffix, suffix) ;
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
  l_intensity = parms.l_intensity ;
  MRISsetVals(mris, -1) ;  /* clear white matter intensities */

  if (!nowhite) {
    fprintf(stderr, "repositioning cortical surface to gray/white boundary\n");

    MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_LABEL, 0) ;
    MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_BORDER_LABEL, 0) ;
    if (mri_T1_white) {
      MRImask(mri_T1_white, mri_labeled, mri_T1_white, BRIGHT_LABEL, 0) ;
      MRImask(mri_T1_white, mri_labeled, mri_T1_white, BRIGHT_BORDER_LABEL, 0) ;
    }
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_T1, "white_masked.mgz") ;
  }
  if (mri_T1_white) {
    MRIfree(&mri_T1);
    mri_T1 = mri_T1_white ; // T1 and T1_white is swapped
  }
  current_sigma = white_sigma ;
  for (n_averages = max_white_averages, i = 0 ;
       n_averages >= min_white_averages ;
       n_averages /= 2, current_sigma /= 2, i++) {
    if (nowhite)
      break ;

    parms.sigma = current_sigma ;
    mri_kernel = MRIgaussian1d(current_sigma, 100) ;
    fprintf(stderr, "smoothing T1 volume with sigma = %2.3f\n",
            current_sigma) ;
    if (!mri_smooth)
      mri_smooth = MRIclone(mri_T1, NULL) ;
#if 0
    MRIconvolveGaussian(mri_T1, mri_smooth, mri_kernel) ;
#endif

    MRIfree(&mri_kernel) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      char fname[STRLEN] ;
      sprintf(fname, "sigma%.0f.mgz", current_sigma) ;
      fprintf(stderr, "writing smoothed volume to %s...\n", fname) ;
      MRIwrite(mri_smooth, fname) ;
    }

    parms.n_averages = n_averages ;
    MRISprintTessellationStats(mris, stderr) ;
    MRIScomputeBorderValues(mris, mri_T1, mri_smooth,
                            MAX_WHITE, max_border_white, min_border_white,
                            min_gray_at_white_border,
                            max_border_white /*max_gray*/, current_sigma,
                            2*max_thickness, parms.fp, GRAY_WHITE, NULL, 0) ;
    MRISfindExpansionRegions(mris) ;
    if (vavgs) {
      fprintf(stderr, "averaging target values for %d iterations...\n",vavgs) ;
      MRISaverageMarkedVals(mris, vavgs) ;
      if (Gdiag_no > 0) {
        VERTEX *v ;
        v = &mris->vertices[Gdiag_no] ;
        fprintf(stderr,"v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f\n",
                Gdiag_no, v->val, v->mean, v->d) ;
      }
    }


    /*
      there are frequently regions of gray whose intensity is fairly
      flat. We want to make sure the surface settles at the innermost
      edge of this region, so on the first pass, set the target
      intensities artificially high so that the surface will move
      all the way to white matter before moving outwards to seek the
      border (I know it's a hack, but it improves the surface in
      a few areas. The alternative is to explicitly put a gradient-seeking
      term in the cost functional instead of just using one to find
      the target intensities).
    */
#if 0
    if (!i) {
      parms.l_nspring = 1.0 ;
      MRISscaleVals(mris, 1.05) ;  /* move inwards on first pass */
    } else
      parms.l_nspring = 0.0 ;
#endif

    if (write_vals) {
      sprintf(fname, "./%s-white%2.2f.w", hemi, current_sigma) ;
      //MRISwriteValues(mris, fname) ;
      MRISwrite(mris, fname);
    }
    MRISpositionSurface(mris, mri_T1, mri_smooth,&parms);
    if (add) {
      for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
      while (MRISdivideLongEdges(mris, max_len) > 0) {}
    }
    if (!n_averages)
      break ;
  }

  if (!nowhite) {
    sprintf(fname, "%s/%s/surf/%s.%s%s%s", sdir, sname,hemi,white_matter_name,
            output_suffix,suffix);
    fprintf(stderr, "writing white matter surface to %s...\n", fname) ;
    MRISaverageVertexPositions(mris, smoothwm) ;
    MRISwrite(mris, fname) ;
#if 0
    if (smoothwm > 0) {
      MRISaverageVertexPositions(mris, smoothwm) ;
      sprintf(fname, "%s/%s/surf/%s.%s%s", sdir, sname,hemi,SMOOTH_NAME,
              suffix);
      fprintf(stderr,"writing smoothed white matter surface to %s...\n",fname);
      MRISwrite(mris, fname) ;
    }
#endif

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
#if 1
      fprintf(stderr, "writing smoothed area to %s\n", fname) ;
      MRISwriteArea(mris, fname) ;
#endif
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
    msec = TimerStop(&then) ;
    fprintf(stderr,
            "refinement took %2.1f minutes\n", (float)msec/(60*1000.0f));
    MRIfree(&mri_T1);
    exit(0) ;
  }

  //////////////////////////////////////////////////////////////////////////////
  // pial surface
  /////////////////////////////////////////////////////////////////////////////
  parms.t = parms.start_t = 0 ;
  sprintf(parms.base_name, "%s%s%s", pial_name, output_suffix, suffix) ;
  parms.niterations = ngray ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ; /* save white-matter */
  parms.l_surf_repulse = l_surf_repulse ;

  MRISsetVals(mris, -1) ;  /* clear target intensities */

  if (smooth && !nowhite) {
    printf("smoothing surface for %d iterations...\n", smooth) ;
    MRISaverageVertexPositions(mris, smooth) ;
  }

  fprintf(stderr, "repositioning cortical surface to gray/csf boundary.\n") ;
  parms.l_repulse = 0 ;
  if (orig_pial) {
    printf("reading initial pial vertex positions from %s...\n", orig_pial) ;

    if (longitudinal) {
      //save final white location
      MRISsaveVertexPositions(mris, TMP_VERTICES);
    }

    if (MRISreadVertexPositions(mris, orig_pial) != NO_ERROR)
      ErrorExit(Gerror, "reading orig pial positions failed") ;

    if (longitudinal) {
      //reset starting point to be in the middle of final white and orig pial
      int vno;
      VERTEX *v;
      //reset the starting position to be slightly inside the orig_pial in the longitudinal case
      for (vno = 0; vno < mris->nvertices; vno++) {
        v = &mris->vertices[vno];
        if (v->ripflag)
          continue;
        v->x = 0.75*v->x + 0.25*v->tx;
        v->y = 0.75*v->y + 0.25*v->ty;
        v->z = 0.75*v->z + 0.25*v->tz;
      }
    }
    MRIScomputeMetricProperties(mris) ; //shouldn't this be done whenever orig_pial is used??? Maybe that's why the cross-intersection was caused
  }
  /*    parms.l_convex = 1000 ;*/
  mri_T1 = mri_T1_pial ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_T1, "p.mgz") ;
  for (j = 0 ; j <= 0 ; parms.l_intensity *= 2, j++)  /* only once for now */
  {
    current_sigma = pial_sigma ;
    for (n_averages = max_pial_averages, i = 0 ;
         n_averages >= min_pial_averages ;
         n_averages /= 2, current_sigma /= 2, i++) {

      parms.sigma = current_sigma ;
      mri_kernel = MRIgaussian1d(current_sigma, 100) ;
      fprintf(stderr, "smoothing T1 volume with sigma = %2.3f\n",
              current_sigma) ;
      parms.n_averages = n_averages ;
      parms.l_tsmooth = l_tsmooth ;
      /*
        replace bright stuff such as eye sockets with 255. Simply zeroing it out
        would make the border always go through the sockets, and ignore subtle
        local minima in intensity at the border of the sockets. Will set to 0
        after border values have been computed so that it doesn't mess up gradients.
      */
      MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_LABEL, 255) ;
      MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_BORDER_LABEL, MID_GRAY) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_T1, "pial_masked.mgz") ;
      MRIScomputeBorderValues(mris, mri_T1, mri_smooth, max_gray,
                              max_gray_at_csf_border, min_gray_at_csf_border,
                              min_csf,(max_csf+max_gray_at_csf_border)/2,
                              current_sigma, 2*max_thickness, parms.fp,
                              GRAY_CSF, NULL, 0) ;
      MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_LABEL, 0) ;
      if (vavgs) {
        fprintf(stderr, "averaging target values for %d iterations...\n",vavgs) ;
        MRISaverageMarkedVals(mris, vavgs) ;
        if (Gdiag_no > 0) {
          VERTEX *v ;
          v = &mris->vertices[Gdiag_no] ;
          fprintf(stderr,"v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f\n",
                  Gdiag_no, v->val, v->mean, v->d) ;
        }
      }

      if (write_vals) {
        sprintf(fname, "./%s-gray%2.2f.w", hemi, current_sigma) ;
        MRISwriteValues(mris, fname) ;
      }
      if (!mri_smooth)
        mri_smooth = MRIcopy(mri_T1, NULL) ;
      MRISpositionSurface(mris, mri_T1, mri_smooth,&parms);
      /*    parms.l_nspring = 0 ;*/
      if (!n_averages)
        break ;
    }
  }

  sprintf(fname, "%s/%s/surf/%s.%s%s%s", sdir, sname, hemi, pial_name,
          output_suffix, suffix) ;
  fprintf(stderr, "writing pial surface to %s...\n", fname) ;
  MRISwrite(mris, fname) ;
  if (create)   /* write out curvature and area files */
  {
    MRIScomputeMetricProperties(mris) ;
    MRIScomputeSecondFundamentalForm(mris) ;
    MRISuseMeanCurvature(mris) ;
    MRISaverageCurvatures(mris, curvature_avgs) ;
    sprintf(fname, "%s.curv.pial%s%s",
            mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", output_suffix,
            suffix);
    fprintf(stderr, "writing smoothed curvature to %s\n", fname) ;
    MRISwriteCurvature(mris, fname) ;
    sprintf(fname, "%s.area.pial%s%s",
            mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", output_suffix,
            suffix);
#if 1
    fprintf(stderr, "writing smoothed area to %s\n", fname) ;
    MRISwriteArea(mris, fname) ;
#endif
    MRISprintTessellationStats(mris, stderr) ;
  }

  if (in_out_in_flag) {
    sprintf(parms.base_name, "%s%s%s", white_matter_name, output_suffix, suffix) ;
    MRIScomputeMetricProperties(mris) ;    /* recompute surface normals */
    MRISstoreMetricProperties(mris) ;
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    fprintf(stderr, "repositioning cortical surface to gray/white boundary\n");

    l_intensity = parms.l_intensity ;
    MRISsetVals(mris, -1) ;  /* clear white matter intensities */

    current_sigma = white_sigma ;
    MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_LABEL, 0) ;
    for (n_averages = max_white_averages, i = 0 ;
         n_averages >= min_white_averages ;
         n_averages /= 2, current_sigma /= 2, i++) {
      if (nowhite)
        break ;

      parms.sigma = current_sigma ;
      mri_kernel = MRIgaussian1d(current_sigma, 100) ;
      fprintf(stderr, "smoothing T1 volume with sigma = %2.3f\n",
              current_sigma) ;
      if (!mri_smooth)
        mri_smooth = MRIclone(mri_T1, NULL) ;
#if 0
      MRIconvolveGaussian(mri_T1, mri_smooth, mri_kernel) ;
#endif
      MRIfree(&mri_kernel) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
        char fname[STRLEN] ;
        sprintf(fname, "sigma%.0f.mgz", current_sigma) ;
        fprintf(stderr, "writing smoothed volume to %s...\n", fname) ;
        MRIwrite(mri_smooth, fname) ;
      }

      parms.n_averages = n_averages ;
      MRISprintTessellationStats(mris, stderr) ;
      MRIScomputeBorderValues(mris, mri_T1, mri_smooth, MAX_WHITE,
                              max_border_white, min_border_white,
                              min_gray_at_white_border, max_border_white /*max_gray*/,
                              current_sigma, 2*max_thickness, parms.fp,
                              GRAY_WHITE, NULL, 0) ;
      MRISfindExpansionRegions(mris) ;
      if (vavgs) {
        fprintf(stderr,"averaging target values for %d iterations...\n",vavgs);
        MRISaverageMarkedVals(mris, vavgs) ;
        if (Gdiag_no > 0) {
          VERTEX *v ;
          v = &mris->vertices[Gdiag_no] ;
          fprintf(stderr,
                  "v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f\n",
                  Gdiag_no, v->val, v->mean, v->d) ;
        }
      }

      if (write_vals) {
        sprintf(fname, "./%s-white%2.2f.w", hemi, current_sigma) ;
        MRISwriteValues(mris, fname) ;
      }
      MRISpositionSurface(mris, mri_T1, mri_smooth,&parms);
      if (!n_averages)
        break ;
    }
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ; /* gray/white surface */
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;  /* pial surface */
    sprintf(fname, "%s/%s/surf/%s.%s2%s", sdir, sname,hemi,white_matter_name,
            suffix);
    fprintf(stderr, "writing gray/white surface to %s...\n", fname) ;
    MRISwrite(mris, fname) ;
  }
  MRIfree(&mri_T1);

  /*  if (!(parms.flags & IPFLAG_NO_SELF_INT_TEST))*/
  {
    fprintf(stderr, "measuring cortical thickness...\n") ;
    if (longitudinal)
      MRISmeasureCorticalThickness(mris, nbhd_size, 5.0) ;
    else
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
      sprintf(fname, "%s/%s/surf/%s.%s%s", sdir, sname, hemi, GRAYMID_NAME,
              suffix) ;
      fprintf(stderr, "writing layer IV surface to %s...\n", fname) ;
      MRISwrite(mris, fname) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    }
  }
  msec = TimerStop(&then) ;
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
  } else if (!stricmp(option, "mode")) {
    use_mode = atoi(argv[2]) ;
    printf("%susing class modes instead of means...\n", use_mode ? "" : "NOT ") ;
    nargs = 1 ;
  } else if (!stricmp(option, "T1") || !stricmp(option, "gvol")) {
    strcpy(T1_name, argv[2]) ;
    printf("using %s as T1 volume...\n", T1_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "wvol")) {
    white_fname = argv[2] ;
    printf("using %s as volume for white matter deformation...\n", white_fname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "hires") || !stricmp(option, "highres")) {
    highres_label = LabelRead(NULL, argv[2]) ;
    if (!highres_label)
      ErrorExit(ERROR_NOFILE, "%s: could not read highres label %s", Progname, argv[2]) ;
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
  } else if (!stricmp(option, "median")) {
    apply_median_filter = 1 ;
  } else if (!stricmp(option, "max_border_white")) {
    max_border_white_set = 1 ;
    max_border_white = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "min_border_white")) {
    min_border_white_set = 1 ;
    min_border_white = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "scale_std")) {
    std_scale = atof(argv[2]);
    printf("scale the estimated WM and GM std by %g \n", std_scale) ;
    nargs = 1 ;
  } else if (!stricmp(option, "min_gray_at_white_border")) {
    min_gray_at_white_border_set = 1 ;
    min_gray_at_white_border = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "max_gray")) {
    max_gray_set = 1 ;
    max_gray = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "max_gray_at_csf_border")) {
    max_gray_at_csf_border_set = 1 ;
    max_gray_at_csf_border = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "min_gray_at_csf_border")) {
    min_gray_at_csf_border_set = 1 ;
    min_gray_at_csf_border = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "min_csf")) {
    min_csf_set = 1 ;
    min_csf = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "max_csf")) {
    max_csf_set = 1 ;
    max_csf = atof(argv[2]) ;
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
      fprintf(stderr, "reading original vertex positions from %s\n", orig_name);
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
      fprintf(stderr, "%screating area and curvature files for wm surface...\n",
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
  printf("%s [options] <subject name> <hemisphere>\n",Progname) ;
  printf("\n");
  printf("options\n");
  printf("  -T1 T1vol : default is %s\n",T1_name);
  printf("  -wvol whitevol <hires>\n");
  printf("  -long : longitudinal\n");
  printf("  -SDIR SUBJECTS_DIR \n");
  printf("  -pial pialsurfname \n");
  printf("  -white whitesurfname \n");
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
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static int
mrisFindMiddleOfGray(MRI_SURFACE *mris) {
  int     vno ;
  VERTEX  *v ;
  float   nx, ny, nz, thickness ;

  MRISaverageCurvatures(mris, 3) ;
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris);
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;
    thickness = 0.5 * v->curv ;
    v->x = v->origx + thickness * nx ;
    v->y = v->origy + thickness * ny ;
    v->z = v->origz + thickness * nz ;
  }
  return(NO_ERROR) ;
}

MRI *
MRIfillVentricle(MRI *mri_inv_lv, MRI *mri_T1, float thresh,
                 int out_label, MRI *mri_dst) {
  BUFTYPE   *pdst, *pinv_lv, out_val, T1_val, inv_lv_val, *pT1 ;
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
        T1_val = *pT1++ ;
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

static MRI *
MRIsmoothMasking(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, int mask_val,
                 int wsize) {
  int      width, height, depth, x, y, z, xi, yi, zi, xk, yk, zk, whalf,
  nvox, mean, avg ;
  BUFTYPE  *psrc, *pdst ;

  whalf = (wsize-1) / 2 ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIcopy(mri_src, NULL) ;

  for ( z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++) {
        mean = *psrc++ ;
        nvox = 1 ;

        /* this is a hack to prevent smoothing of non-white values */
        if (MRIvox(mri_mask, x, y, z) > WM_MIN_VAL) {
          avg = 0 ;  /* only average if a masked value is close to this one */
          for (zk = -whalf ; zk <= whalf ; zk++) {
            zi = mri_mask->zi[z+zk] ;
            for (yk = -whalf ; yk <= whalf ; yk++) {
              yi = mri_mask->yi[y+yk] ;
              for (xk = -whalf ; xk <= whalf ; xk++) {
                xi = mri_mask->xi[x+xk] ;
                if (MRIvox(mri_mask, xi, yi, zi) == mask_val) {
                  avg = 1 ;
                  break ;
                }
              }
              if (avg)
                break ;
            }
            if (avg)
              break ;
          }
          if (avg) {
            for (zk = -whalf ; zk <= whalf ; zk++) {
              zi = mri_mask->zi[z+zk] ;
              for (yk = -whalf ; yk <= whalf ; yk++) {
                yi = mri_mask->yi[y+yk] ;
                for (xk = -whalf ; xk <= whalf ; xk++) {
                  xi = mri_mask->xi[x+xk] ;
                  if (MRIvox(mri_mask, xi, yi, zi) >= WM_MIN_VAL) {
                    mean += MRIvox(mri_src, xi, yi, zi) ;
                    nvox++ ;
                  }
                }
              }
            }
          }

        }
        *pdst++ = (BUFTYPE)nint((float)mean/(float)nvox) ;
      }
    }
  }
  return(mri_dst) ;
}

int
MRISfindExpansionRegions(MRI_SURFACE *mris) {
  int    vno, num, n, num_long, total ;
  float  d, dsq, mean, std, dist ;
  VERTEX *v, *vn ;

  d = dsq = 0.0f ;
  for (total = num = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
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
    v = &mris->vertices[vno] ;
    v->curv = 0 ;
    if (v->ripflag || v->val <= 0)
      continue ;
    if (fabs(v->d) < mean+2*std)
      continue ;
    for (num_long = num = 1, n = 0 ; n < v->vnum ; n++) {
      vn = &mris->vertices[v->v[n]] ;
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
MRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm) {
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
  MRIdilateInvThreshLabel(mri_labeled, mri_T1, mri_labeled, BRIGHT_LABEL, 3, 0) ;

  MRIfree(&mri_tmp) ;
  return(mri_labeled) ;
}
static float
check_contrast_direction(MRI_SURFACE *mris,MRI *mri_T1) {
  int     vno, n ;
  VERTEX  *v ;
  Real    x, y, z, xw, yw, zw, val, mean_inside, mean_outside ;

  mean_inside = mean_outside = 0.0 ;
  for (n = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag != 0)
      continue ;
    x = v->x+0.5*v->nx ;
    y = v->y+0.5*v->ny ;
    z = v->z+0.5*v->nz ;
    MRIsurfaceRASToVoxel(mri_T1, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_T1, xw, yw, zw, &val) ;
    mean_outside += val ;

    x = v->x-0.5*v->nx ;
    y = v->y-0.5*v->ny ;
    z = v->z-0.5*v->nz ;
    MRIsurfaceRASToVoxel(mri_T1, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_T1, xw, yw, zw, &val) ;
    mean_inside += val ;
    n++ ;
  }
  mean_inside /= (float)n ;
  mean_outside /= (float)n ;
  printf("mean inside = %2.1f, mean outside = %2.1f\n", mean_inside, mean_outside) ;
  return(mean_inside - mean_outside) ;
}

