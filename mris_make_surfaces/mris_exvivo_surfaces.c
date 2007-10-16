/**
 * @file  mris_exvivo_surfaces.c
 * @brief tailored for the ex-vivo MEF data
 *
 * This version is tailored for the ex-vivo MEF data. There are 6 channels
 * available, but using only 2 of them; 30 and 5 channel; again, 5 is better
 * for gray/white surface, and 30 is better for pial
 * Exvivo data only have one hemi, assume it to be left (255)
 * This version also gets rid of mri_em_seg, instead, it uses filled.mgz to
 * compute WM and GM stats.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/10/16 23:01:24 $
 *    $Revision: 1.5 $
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

static char vcid[] = 
"$Id: mris_exvivo_surfaces.c,v 1.5 2007/10/16 23:01:24 fischl Exp $";

int main(int argc, char *argv[]) ;

#define BRIGHT_LABEL         130
#define BRIGHT_BORDER_LABEL  100

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  mrisFindMiddleOfGray(MRI_SURFACE *mris) ;

int MRISfindExpansionRegions(MRI_SURFACE *mris) ;

int MRISaverageMarkedValbaks(MRI_SURFACE *mris, int navgs);
static int  MRIcomputeClassStatistics_mef(MRI *mri_T1_30, MRI *mri_T1_5, 
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
                                 float max_thickness, FILE *log_fp, 
                                 int formalin, int first);

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

static int formalin = 0 ;
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

static int auto_detect_stats = 1 ;

static int in_out_in_flag = 0 ;
static int apply_median_filter = 0 ;

static int nbhd_size = 20 ;

static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE    1.0
static float base_dt_scale = BASE_DT_SCALE ;

static char *label_name = NULL ;
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

static char *orig_name = ORIG_NAME ;
static char *suffix = "" ;
static char *output_suffix = "" ;
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

static int MGZ = 1; // for use with MGZ format

static int longitudinal = 0;


int
main(int argc, char *argv[]) {
  char          **av, *hemi, *sname, *cp, fname[STRLEN], mdir[STRLEN];
  int           ac, nargs, i, label_val, replace_val, msec, n_averages, j ;
  MRI_SURFACE   *mris ;
  MRI           *mri_filled, *mri_T1_30, *mri_T1_5; // *mri_labeled;
  MRI           *mri_em_seg = NULL;
  float         max_len ;
  //need to be estimated from EM segmentation; 
  //need to mask out cerebellum though
  float         white_mean[2], white_std[2], gray_mean[2], gray_std[2];
  double        l_intensity, current_sigma ;
  struct timeb  then ;

  char cmdline[CMD_LINE_LEN] ;

  make_cmd_version_string 
    (argc, argv, 
     "$Id: mris_exvivo_surfaces.c,v 1.5 2007/10/16 23:01:24 fischl Exp $", 
     "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: mris_exvivo_surfaces.c,v 1.5 2007/10/16 23:01:24 fischl Exp $", 
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.projection = NO_PROJECTION ;
  parms.tol = 5e-3 ;
  parms.check_tol = 1 ;
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

  sprintf(fname, "%s/%s/mri/filled", sdir, sname) ;
  if (MGZ) strcat(fname, ".mgz");
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_filled = MRIread(fname) ;
  if (!mri_filled)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;
  ////////////////////////////// we can handle only conformed volumes
  // setMRIforSurface(mri_filled);

  if (!stricmp(hemi, "lh")) //should always be left
  { label_val = lh_label ;
    replace_val = rh_label ;
  } else {
    label_val = rh_label ;
    replace_val = lh_label ;
  }

  sprintf(fname, "%s/%s/mri/%s", sdir, sname, T1_30_name) ;
  if (MGZ) strcat(fname, ".mgz");
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_T1_30 = MRIread(fname) ;

  if (!mri_T1_30)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

  sprintf(fname, "%s/%s/mri/%s", sdir, sname, T1_5_name) ;
  if (MGZ) strcat(fname, ".mgz");
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_T1_5 = MRIread(fname) ;

  if (!mri_T1_5)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

#if 0
  sprintf(fname, "%s/%s/mri/%s", sdir, sname, em_name) ;
  if (MGZ) strcat(fname, ".mgz");
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_em_seg = MRIread(fname) ;

  if (!mri_em_seg)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;
#endif

  /////////////////////////////////////////
  // setMRIforSurface(mri_T1_30);
  // setMRIforSurface(mri_T1_5);

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
#if 0
  MRInormalizeMEF(mri_T1_30, mri_em_seg);
  MRInormalizeMEF(mri_T1_5, mri_em_seg);
#else
  MRInormalizeMEF(mri_T1_30, mri_filled);
  MRInormalizeMEF(mri_T1_5, mri_filled);

  MRIwrite(mri_T1_30, "normalized_30.mgz");
  MRIwrite(mri_T1_5, "normalized_5.mgz");
#endif

  /* remove other hemi */
  MRIdilateLabel(mri_filled, mri_filled, replace_val, 1) ;

  MRImask(mri_T1_30, mri_filled, mri_T1_30, replace_val,0) ;
  MRImask(mri_T1_5, mri_filled, mri_T1_5, replace_val,0) ;


  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_T1_30, "r30.mgz") ;

  //the following two steps might as well be ignored
  //  MRIsmoothBrightWM(mri_T1_30, mri_wm) ;
  //  mri_labeled = MRIfindBrightNonWM(mri_T1_30, mri_wm) ;

#if 0
  MRIcomputeClassStatistics_mef(mri_T1_30, mri_T1_5, mri_em_seg, 
                                white_mean, white_std, gray_mean, gray_std) ;
#else
  MRIcomputeClassStatistics_mef(mri_T1_30, mri_T1_5, mri_filled, 
                                white_mean, white_std, gray_mean, gray_std) ;
#endif

  printf("WM: mean = (%g, %g) std = (%g, %g)\n", 
         white_mean[0], white_mean[1], white_std[0], white_std[1]);
  printf("GM: mean = (%g, %g) std = (%g, %g)\n", 
         gray_mean[0], gray_mean[1], gray_std[0], gray_std[1]);

  sprintf(fname, "%s/%s/surf/%s.%s%s", sdir, sname, hemi, orig_name, suffix) ;
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
    MRISsetNeighborhoodSize(mris, nbrs) ;

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
  l_intensity = parms.l_intensity ;

  //what if I need a vector-values intensity??
  MRISsetVals(mris, -1) ;  /* clear white matter intensities */
  MRIScopyValToValBak(mris);

  current_sigma = white_sigma ;

  if (label_name)
  {
    LABEL *area ;

    area = LabelRead(NULL, label_name) ;
    if (area == NULL)
      exit(Gerror) ;
    LabelRipRestOfSurface(area, mris) ;
    LabelFree(&area) ;
  }
  for (n_averages = max_white_averages, i = 0 ;
       n_averages >= min_white_averages ;
       n_averages /= 2, current_sigma /= 2, i++) {
    if (nowhite)
      break ;

    parms.sigma = current_sigma ;

    parms.n_averages = n_averages ;
    MRISprintTessellationStats(mris, stderr) ;

    //This following function needs major modification
    MRIScomputeBorderValues_MEF_WHITE(mris, mri_T1_30, mri_T1_30, mri_T1_5,
                                      white_mean, white_std, gray_mean,
                                      gray_std,
                                      current_sigma,
                                      max_thickness, parms.fp) ;

    //what does this do?
    MRISfindExpansionRegions(mris) ;

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
    MRISpositionSurface_mef(mris, mri_T1_30, mri_T1_5, &parms, 0.3, 0.7);

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
    msec = TimerStop(&then) ;
    fprintf(stderr,
            "refinement took %2.1f minutes\n", (float)msec/(60*1000.0f));
    MRIfree(&mri_T1_30);
    MRIfree(&mri_T1_5);
    if (mri_em_seg)
      MRIfree(&mri_em_seg);
    exit(0) ;
  }

  ///////////////////////////////////////////////////////////////////////////
  // pial surface
  ///////////////////////////////////////////////////////////////////////////
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
                                          was caused*/
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
        after border values have been computed so 
        that it doesn't mess up gradients.
      */
      //   MRImask(mri_T1_30, mri_labeled, mri_T1_30, BRIGHT_LABEL, 255) ;
      //  MRImask(mri_T1_30, mri_labeled, mri_T1_30, 
      //BRIGHT_BORDER_LABEL, MID_GRAY) ;
      //  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      //  MRIwrite(mri_T1_30, "pial_masked.mgz") ;

      //The following function need major modification
      MRIScomputeBorderValues_MEF_PIAL
        (mris, mri_T1_30, mri_T1_30, mri_T1_5,
         white_mean, white_std, gray_mean, gray_std,
         current_sigma, 2*max_thickness, parms.fp, formalin, j) ;

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
      if (formalin)
        MRISpositionSurface_mef(mris, mri_T1_30, mri_T1_5, &parms, 1.0, 0.0);
      else
        MRISpositionSurface_mef(mris, mri_T1_30, mri_T1_5, &parms, .95, 0.05);
      /*    parms.l_nspring = 0 ;*/
      if (!n_averages)
        break ;
    }
  }

  sprintf(fname, "%s/%s/surf/%s.%s%s%s", sdir, sname, hemi, pial_name,
          output_suffix, suffix) ;
  fprintf(stderr, "writing pial surface to %s...\n", fname) ;
  MRISwrite(mris, fname) ;

  MRIfree(&mri_filled) ;
  MRIfree(&mri_T1_30);
  MRIfree(&mri_T1_5);
  if (mri_em_seg != NULL)
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
  }
  else if (!stricmp(option, "formalin")) {
    formalin = atoi(argv[2]) ;
    fprintf(stderr,  "assuming hemisphere is %sembedded in formalin\n",
            formalin ? "" : "not ") ;
    nargs = 1 ;
  } else if (!stricmp(option, "wvol")) {
    white_fname = argv[2] ;
    printf("using %s as volume for white matter deformation...\n", 
           white_fname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "hires") || !stricmp(option, "highres")) {
    highres_label = LabelRead(NULL, argv[2]) ;
    if (!highres_label)
      ErrorExit(ERROR_NOFILE, "%s: could not read highres label %s", 
                Progname, argv[2]) ;
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
      fprintf(stderr, "reading original vertex positions from %s\n", 
              orig_name);
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
  case 'L':
    label_name = argv[2] ;
    nargs = 1 ;
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
  fprintf(stderr,
          "-formalin <0,1>  assume hemi is in formalin\n") ;
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

#define STEP_SIZE 0.1

int
MRIScomputeBorderValues_MEF_WHITE(MRI_SURFACE *mris, 
                                  MRI *mri_em_combined, MRI *mri_30,
                                  MRI *mri_5, 
                                  float wm_mean[2], float wm_std[2],
                                  float gm_mean[2], float gm_std[2],
                                  double sigma,
                                  float max_thickness, FILE *log_fp) {
  Real    val30, val5, x, y, z, max_mag_val30, 
    max_mag_val5, xw, yw, zw,mag,max_mag, max_mag_dist=0.0f,
    min_val30, min_val5, inward_dist,outward_dist,xw1,yw1,zw1,
    min_val_dist, orig_dist, dx, dy, dz;
  Real low30, high30, low5, high5;
  Real previous_mag, next_mag, previous_mag30, 
    previous_mag5, next_mag30, next_mag5;
  Real mag30, mag5, previous_val30, previous_val5, next_val30;
  int     total_vertices, vno, nmissing = 0, 
    nout = 0, nin = 0, nfound = 0,
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
  low30 = gm_mean[0] - 2*gm_std[0];
  low5 = wm_mean[1] - wm_std[1];
  high30 = wm_mean[0] + 2*wm_std[0];
  high5 = gm_mean[1] + gm_std[1];

  for (total_vertices = vno = 0 ; vno < mris->nvertices ; vno++) {

    if (Gdiag_no > 0) {
      if (vno == Gdiag_no) {
        printf("DEBUG\n");
      }
    }

    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x ;
    y = v->y ;
    z = v->z ;
    MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
    x = v->x + v->nx ;
    y = v->y + v->ny ;
    z = v->z + v->nz ;
    MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw1, &yw1, &zw1) ;
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
        MRIsurfaceRASToVoxelCached(mri_30, x, y, z, &xw, &yw, &zw) ;
        //       MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, nx, ny,nz,&mag30,  current_sigma);

        // if (mag30 >= 0.0) //this condition doesn't make much sense to me
        //  break ;


        MRIsampleVolume(mri_30, xw, yw, zw, &val30) ;
        MRIsampleVolume(mri_5, xw, yw, zw, &val5) ;

        if (Gdiag_no > 0) {
          if (vno == Gdiag_no) {
            printf("dist = %g, (x,y,z) = (%g, %g, %g), "
                   "val30 =%g, val5 = %g \n", dist, xw, yw, zw, val30, val5);
          }
        }


        //exvivo data, WM is darker than GM in flash30 as well
        if (val30 < low30 || val30 > high30 || val5 < low5 )
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
        MRIsurfaceRASToVoxelCached(mri_30, x, y, z, &xw, &yw, &zw) ;
        //       MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, nx, ny,nz, &mag30, current_sigma);
        // MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, nx, ny,nz, &mag5, current_sigma);

        //       if (mag30 >= 0.0)
        // break ;
        MRIsampleVolume(mri_30, xw, yw, zw, &val30) ;
        MRIsampleVolume(mri_5, xw, yw, zw, &val5) ;
        if (val30 < low30 || val5 > high5 )
          break ;
      }
      outward_dist = dist-step_size/2 ;
      if (!finite(outward_dist))
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
      MRIsampleVolume(mri_5, xw, yw, zw, &previous_val5) ;

      if (Gdiag_no > 0) {
        if (vno == Gdiag_no) {
          printf("dist = %g, previous30 = %g, prev5 = %g \n", 
                 dist, previous_val30, previous_val5);
        }
      }

      /* the previous point was inside the surface */
      if (previous_val30 > low30 && previous_val5 < gm_mean[1]) {
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

        if (Gdiag_no > 0) {
          if (vno == Gdiag_no) {
            printf("dist = %g, (x,y,z) = (%g, %g, %g), "
                   "mag5 = %g, pre_mag5 =%g, next_mag5 = %g \n", 
                   dist, xw, yw, zw, mag5, previous_mag5, next_mag5);
          }
        }


        //flash30 and flash5 have same contrast at 
        //gray/white boundary for ex-vivo; flash5 is better
        // the weights are arbitrary for now,
        mag = - mag5;
        previous_mag = - previous_mag5;
        next_mag = - next_mag5;
        if (val5 < min_val5) {
          min_val30 = val30;
          min_val5 = val5;
          min_val_dist = dist;
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
          (val5 <= (gm_mean[1]+ 1.2*gm_std[1]))
#if 0
          && 
          (val5 >= wm_mean[1]-0.5*wm_std[1]) //is this too restrictive??
#endif
          )
        {
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
          if ((((local_max_found == 0) &&
                (fabs(mag) > max_mag)) ||
               ((fabs(mag) > 2*max_mag)))
              &&
              (val5 <= (gm_mean[1] + gm_std[1])) &&
              (val5 >= (wm_mean[1] - 0.5*wm_mean[1]) )
             ) {
            local_max_found = 0;
            x = v->x + v->nx*(dist+1) ;
            y = v->y + v->ny*(dist+1) ;
            z = v->z + v->nz*(dist+1) ;
            MRISsurfaceRASToVoxel(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
            MRIsampleVolume(mri_30, xw, yw, zw, &next_val30) ;
            // MRIsampleVolume(mri_5, xw, yw, zw, &next_val5) ;
            if ((next_val30 >= (gm_mean[0] - 2*gm_std[0])) &&
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
         "v %d, target value = 30:%2.1f, 5:%2.1f, mag = %2.1f, dist = %2.2f, %s\n",
         Gdiag_no, v->val, v->valbak, v->mean, v->d,
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
#define F30IND 0
#define F5IND  1


int
MRIScomputeBorderValues_MEF_PIAL(MRI_SURFACE *mris, 
                                 MRI *mri_em_combined, MRI *mri_30,
                                 MRI *mri_5, float wm_mean[2], float wm_std[2],
                                 float gm_mean[2], float gm_std[2],
                                 double sigma,
                                 float max_thickness, FILE *log_fp,
                                 int formalin, int first) {
  //for pial surface, dura is a problem if I still use 
  //original images, really need to use the membership functions
  //or equivalently, the EM_combined
  //dura may be still low at both flash30 and flash5, but not that low
  Real  val,  val30, val5, x, y, z, 
    max_mag_val30, max_mag_val5, xw, yw, zw,mag,max_mag, max_mag_dist=0.0f,
    min_val30, min_val5, inward_dist,outward_dist,xw1,yw1,zw1,
    min_val_dist, orig_dist, dx, dy, dz;
  Real previous_mag, next_mag, 
    previous_mag30, previous_mag5, next_mag30, next_mag5;
  Real mag30, mag5, previous_val30, next_val30;
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
      if (first)
        dist = 0 ;
      else
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
          
          //The inside check can be restrictive, since pial surface have to 
          //move outside anyway. No not necessary for longitudinal method
          if (mag30 >= 0.0 || mag5 > 0.0)
            break ;
          MRIsampleVolume(mri_30, xw, yw, zw, &val30) ;
          MRIsampleVolume(mri_5, xw, yw, zw, &val5) ;
          
          if (val30 > wm_mean[0])
            break ;
        }
      inward_dist = dist+step_size/2 ;
      for (dist = first ? 1 : 0 ; dist < max_thickness ; dist += step_size) {
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
        //  MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, nx, ny,nz, &mag30, current_sigma);
        // MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, nx, ny,nz, &mag5, current_sigma);

        //       if (mag30 >= 0.0)
        // break ;
        MRIsampleVolume(mri_30, xw, yw, zw, &val30) ;
        // MRIsampleVolume(mri_em_combined, xw, yw, zw, &val5) ; //borrowed
        if (val30 < (gm_mean[0] - 3*gm_std[0]))
          break ;
      }
      outward_dist = dist;
      if (!finite(outward_dist))
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
      MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
      MRIsampleVolume(mri_30, xw, yw, zw, &previous_val30) ;
      //   MRIsampleVolume(mri_em_combined, xw, yw, zw, &previous_val5) ; //borrowed

      /* the previous point was inside the surface */
      //   if (previous_val30 < wm_mean[0] &&  previous_val5 > 50)
      if (previous_val30 < wm_mean[0] && 
          previous_val30 > (gm_mean[0] - 3*gm_std[0]) ) {
        /* see if we are at a local maximum in the gradient magnitude */
        x = v->x + v->nx*dist ;
        y = v->y + v->ny*dist ;
        z = v->z + v->nz*dist ;
        MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolume(mri_30, xw, yw, zw, &val30) ;
        MRIsampleVolume(mri_5, xw, yw, zw, &val5) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, 
                                       nx, ny, nz,&mag30,sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, 
                                       nx, ny, nz,&mag5,sigma);

        x = v->x + v->nx*(dist+STEP_SIZE) ;
        y = v->y + v->ny*(dist+STEP_SIZE) ;
        z = v->z + v->nz*(dist+STEP_SIZE) ;
        MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, nx, ny, nz,
                                       &next_mag30, sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, nx, ny, nz,
                                       &next_mag5, sigma);

        x = v->x + v->nx*(dist-STEP_SIZE) ;
        y = v->y + v->ny*(dist-STEP_SIZE) ;
        z = v->z + v->nz*(dist-STEP_SIZE) ;
        MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolumeDerivativeScale(mri_30, xw, yw, zw, nx, ny, nz,
                                       &previous_mag30, sigma);
        MRIsampleVolumeDerivativeScale(mri_5, xw, yw, zw, nx, ny, nz,
                                       &previous_mag5, sigma);


        //flash30 and flash5 have same contrast at gray/csf boundary
        // the weights are arbitrary for now,
        if (formalin)
        {
          mag = mag30 ;
          previous_mag = previous_mag30 ;
          next_mag = next_mag30 ;
        }
        else
        {
          mag = mag30*0.8 + mag5*0.2;
          previous_mag = previous_mag30*0.8 + previous_mag5*0.2;
          next_mag = next_mag30*0.8 + next_mag5*0.2;
        }

        if (val30 < min_val30) {
          min_val30 = val30 ;  /* used if no gradient max is found */
          min_val5 = val5;
          min_val_dist = dist ;
        }


        /*
          sample the next val we would process. If it is too low, then we
          have definitely reached the border, and the current gradient
          should be considered a local max.
        */
        x = v->x + v->nx*(dist+STEP_SIZE) ;
        y = v->y + v->ny*(dist+STEP_SIZE) ;
        z = v->z + v->nz*(dist+STEP_SIZE) ;
        MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        // MRIsampleVolume(mri_em_combined, xw, yw, zw, &next_val5) ; //borrowed
        MRIsampleVolume(mri_30, xw, yw, zw, &next_val30) ;
        if (next_val30 < 50)
          next_mag = 0 ;

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

          MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
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
            MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
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
        MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw, &yw, &zw) ;
        //       MRIsampleVolume(mri_em_combined, xw, yw, zw, &val) ; //borrowed
        MRIsampleVolume(mri_30, xw, yw, zw, &val) ; //borrowed
        if (val < gm_mean[0]-3*gm_std[0]) {
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
        //may better only use flash30 for intensity term
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

#define STD_TARGET 6

static int  MRInormalizeMEF(MRI *mri_src, MRI *mri_wm) {
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
        val1 = MRIgetVoxVal(mri_wm, x, y, z,0);
        val2 = MRIgetVoxVal(mri_src, x, y, z, 0);
        if (val2 <= 1e-10 || val1 < 10) continue;

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
        if (tmpval < 0) tmpval = 0;
        MRIsetVoxVal(mri_src,x,y,z,0,tmpval);
      }

  return 0;
}


static int  MRIcomputeClassStatistics_mef(MRI *mri_T1_30, MRI *mri_T1_5, 
                                          MRI *mri_fill, 
                                          float *white_mean, float *white_std, 
                                          float *gray_mean, float *gray_std) 
{
  // this version will follow MRIcomputeClassStatistics()
  int channel;
  MRI     *mri_border, *mri_tmp;
  float val30, val5;
  int x, y,z, width, height, depth;
  double sum_wm, sum_gm;
  BUFTYPE border_label;

  mri_tmp = MRIbinarize(mri_fill, NULL, 10, MRI_NOT_WHITE, MRI_WHITE);
  mri_border = MRImarkBorderVoxels(mri_tmp, NULL) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    MRIwrite(mri_border, "border.mgz") ;

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
  for (z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++) {
        val30 = MRIgetVoxVal(mri_T1_30, x, y, z,0);
        val5 = MRIgetVoxVal(mri_T1_5, x, y, z,0);
        border_label = MRIvox(mri_border, x, y, z) ;

        if (val30 <= 1e-10 || 
            val5 <= 1e-10 || 
            border_label == MRI_AMBIGUOUS) continue;
        if (border_label == MRI_WHITE) {
          if (val30 > 70) {
            white_mean[0] += val30;
            white_mean[1] += val5;

            white_std[0] += val30*val30;
            white_std[1] += val5*val5;

            sum_wm += 1;
          }
        } else if (border_label == MRI_NOT_WHITE) {
          if (val30 > 30) {
            gray_mean[0] += val30;
            gray_mean[1] += val5;

            gray_std[0] += val30*val30;
            gray_std[1] += val5*val5;

            sum_gm += 1;
          }
        }
      }

  if (sum_wm < 10 || sum_gm < 10) {
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

  MRIfree(&mri_border);
  MRIfree(&mri_tmp);
  return 0;
}



int MRISaverageMarkedValbaks(MRI_SURFACE *mris, int navgs) {
  int    i, vno, vnb, *pnb, vnum ;
  float  val, num ;
  VERTEX *v, *vn ;

  for (i = 0 ; i < navgs ; i++) {
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->marked == 0)
        continue ;
      val = v->valbak ;
      pnb = v->v ;
      vnum = v->vnum ;
      for (num = 0.0f, vnb = 0 ; vnb < vnum ; vnb++) {
        vn = &mris->vertices[*pnb++] ;    /* neighboring vertex pointer */
        if (vn->ripflag || v->marked == 0 )
          continue ;
        num++ ;
        val += vn->valbak ;
      }
      num++ ;  /*  account for central vertex */
      v->tdx = val / num ;
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->ripflag ||  v->marked == 0)
        continue ;
      v->valbak = v->tdx ;
    }
  }
  return(NO_ERROR) ;
}
