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


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "gca.h"
#include "cma.h"
#include "mrinorm.h"
#include "version.h"

const char         *Progname ;
static MORPH_PARMS  parms ;

static int mean_flag = 0 ;
static int fill_val = 0 ;
static int radius = 7 ;
static char *mask_fname = NULL ;
static char *norm_fname = NULL ;

static char *example_T1 = NULL ;
static char *example_segmentation = NULL ;

static double TR = -1 ;
static double alpha = -1 ;
static double TE = -1 ;

static char *sample_fname = NULL ;
static char *transformed_sample_fname = NULL ;
static char *normalized_transformed_sample_fname = NULL ;
static char *ctl_point_fname = NULL ;
static int novar = 0 ;

#define MAX_SPACING  8
static int max_spacing = MAX_SPACING ;
static int nscales = 1 ;

static char *xform_fname = NULL ;
static int use_contrast = 0 ;
static float min_prior = MIN_PRIOR ;
static double tol ;
static double tx = 0.0 ;
static double ty = 0.0 ;
static double tz = 0.0 ;
static double rzrot = 0.0 ;
static double rxrot = 0.0 ;
static double ryrot = 0.0 ;

static FILE *diag_fp = NULL ;

static LTA *Glta = NULL ;
static TRANSFORM *transform = NULL ;

double compute_mah_dist(float *wm_means,
                        MATRIX *m_wm_covar,
                        float *vals,
                        int ninputs) ;
static double find_optimal_linear_xform(GCA *gca,
                                        GCA_SAMPLE *gcas,
                                        MRI *mri,
                                        int nsamples,
                                        MATRIX *m_L, MATRIX *m_origin,
                                        float min_angle, float max_angle,
                                        float min_scale, float max_scale,
                                        float min_trans, float max_trans,
                                        float angle_steps, float scale_steps,
                                        float trans_steps,
                                        int nreductions);
static MRI * fill_brain_volume(MRI *mri,
                               GCA *gca,
                               TRANSFORM *transform,
                               int radius) ;
MRI *MRIremoveFace(MRI *mri_src,
                   MRI *mri_dst,
                   LTA *lta,
                   GCA *gca,
                   GCA *gca_face,
                   int radius,
                   int fill_val) ;
float compareLogSampleProbability(GCA *gca, GCA_SAMPLE *gcas,
                                  MRI *mri_inputs, TRANSFORM *transform,
                                  int nsamples) ;
#if 0
static int brain_in_region(GCA *gca,
                           MRI *mri,
                           TRANSFORM *transform,
                           int x, int y, int z, int whalf) ;
static int recompute_labels(MRI *mri, GCA *gca, GCA_SAMPLE *gcas,
                            int nsamples, MATRIX *m_L) ;
static float compare_transform(MRI *mri,
                               GCA *gca,
                               GCA_SAMPLE *gcas,
                               int nsamples,
                               MATRIX *m_L) ;
#endif
static int translation_only = 0 ;
static int get_option(int argc, char *argv[]) ;
static int register_mri(MRI *mri_in,
                        GCA *gca,
                        MP *parms,
                        int passno,
                        int spacing) ;

static char *renormalization_fname = NULL ;
static char *tissue_parms_fname = NULL ;
static int center = 1 ;
static int nreductions = 1 ;
static int noscale = 0 ;
static int noiscale = 0 ;
static int num_xforms = 1 ;
static int transform_loaded = 0 ;
static char *gca_mean_fname = NULL ;

#if 0
static double find_optimal_3x4(GCA *gca,
                               GCA_SAMPLE *gcas, MRI *mri, int nsamples,
                               MATRIX *m_L, float scale, int nsteps) ;
static MATRIX *
update_optimal_transform(MRI *mri, GCA *gca,
                         GCA_SAMPLE *gcas, int nsamples,
                         MATRIX *m_L,
                         double scale, int write_iterations, int nsteps,
                         int nreductions) ;
#endif
static MATRIX *find_optimal_transform(MRI *mri_in, GCA *gca, GCA_SAMPLE *gcas,
                                      int nsamples, MATRIX *m_L, int passno,
                                      int write_iterations, int spacing) ;
static double find_optimal_translation(GCA *gca, GCA_SAMPLE *gcas, MRI *mri,
                                       int nsamples, MATRIX *m_L,
                                       float min_trans, float max_trans,
                                       float trans_steps, int nreductions) ;
#if 0
static double find_optimal_scaling(GCA *gca, GCA_SAMPLE *gcas, MRI *mri,
                                   int nsamples, MATRIX *m_L, MATRIX *m_origin,
                                   float min_scale, float max_scale,
                                   float scale_steps, int nreductions) ;
static double find_optimal_rotation(GCA *gca, GCA_SAMPLE *gcas, MRI *mri,
                                    int nsamples,
                                    MATRIX *m_L, MATRIX *m_origin,
                                    float min_angle, float max_angle,
                                    float angle_steps,
                                    int nreductions) ;
static double find_optimal_scaling_and_rotation(GCA *gca, GCA_SAMPLE *gcas,
    MRI *mri,
    int nsamples,
    MATRIX *m_L, MATRIX *m_origin,
    float min_angle,
    float max_angle,
    float min_scale,
    float max_scale,
    float angle_steps,
    float scale_steps,
    int nreductions) ;
#endif

static double blur_sigma = 0.0f ;
double local_GCAcomputeLogSampleProbability(GCA *gca,
    GCA_SAMPLE *gcas,
    MRI *mri,
    MATRIX *m_L,
    int nsamples, double clamp);
/*
   command line consists of three inputs:

   argv[1]  - directory containing 'canonical' brain
   argv[2]  - directory containing brain to be registered
   argv[3]  - directory in which to write out registered brain.
*/

#define NPARMS           12
#define NSAMPLES        (NPARMS*20)
#define DEFAULT_CTL_POINT_PCT   .25
static double ctl_point_pct = DEFAULT_CTL_POINT_PCT ;
static int nsamples = NSAMPLES ;

static int exclude_list[MAX_CMA_LABEL+1] ;

int
main(int argc, char *argv[]) {
  char         *gca_fname, *in_fname, *out_fname,
  fname[STRLEN], **av, *gca_face_fname ;
  MRI          *mri_in, *mri_out, *mri_orig ;
  GCA          *gca /*, *gca_tmp, *gca_reduced*/, *gca_face ;
  int          ac, nargs, i ;
  int          msec, minutes, seconds, scale, spacing ;
  Timer start ;
  float        old_log_p, log_p ;

  // disable license checking
  sprintf(fname, "S%sER%sRONT%sOR", "URF", "_F", "DO") ;
  setenv(fname,"1",0);

  for (i = 0 ; i < MAX_CMA_LABEL ; i++) {
    switch (i) {
    case Unknown:
    case Left_Hippocampus:
    case Right_Hippocampus:
    case Left_Amygdala:
    case Right_Amygdala:
    case Left_Caudate:
    case Right_Caudate:
    case Left_Pallidum:
    case Right_Pallidum:
    case Left_Putamen:
    case Right_Putamen:
    case Left_Thalamus:
    case Right_Thalamus:
    case Left_Accumbens_area:
    case Right_Accumbens_area:
      exclude_list[i] = 1 ;
      break ;
    case Left_Cerebellum_Exterior:
    case Left_Cerebellum_Cortex:
    case Left_Cerebellum_White_Matter:
    case Right_Cerebellum_Exterior:
    case Right_Cerebellum_Cortex:
    case Right_Cerebellum_White_Matter:
      exclude_list[i] = 0 ;
      break ;
    default:
      exclude_list[i] = 0 ;
      break ;
    }
  }
  parms.l_intensity = 1.0f ;
  parms.niterations = 100 ;
  parms.levels = -1 ;   /* use default */
  parms.dt = 1e-6 ;  /* was 5e-6 */
  parms.tol = INTEGRATION_TOL*5 ;

  parms.max_levels = 0 ;
  parms.dt = 5e-6 ;  /* was 5e-6 */
  tol = parms.tol = 1e-3 ;
  parms.momentum = 0.8 ;
  parms.niterations = 25 ;
  Progname = argv[0] ;


  setRandomSeed(-1L) ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  nargs = handleVersionOption(argc, argv, "mri_deface");
  argc -= nargs ;
  if (1 == argc)
    ErrorExit
    (ERROR_BADPARM,
     "usage: %s <in volume> <brain template> "
     "<face template> <defaced output volume>\n", argv[0]) ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
    ErrorExit
    (ERROR_BADPARM,
     "usage: %s <in volume> <brain template> "
     "<face template> <defaced output volume>\n", argv[0]) ;

  in_fname = argv[1] ;
  gca_fname = argv[2] ;
  gca_face_fname = argv[3] ;
  out_fname = argv[4] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
  Gdiag |= DIAG_WRITE ;
  printf("logging results to %s.log\n", parms.base_name) ;

  start.reset() ;
  printf("reading '%s'...\n", gca_fname) ;
  fflush(stdout) ;
  gca = GCAread(gca_fname) ;
  if (gca == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open GCA %s.\n",
              Progname, gca_fname) ;
  parms.vgca = gca ;  /* for diagnostics in MRIemAlign */
  if (Ggca_label >= 0) {
    float means[MAX_GCA_INPUTS] ;
    int   i ;

    GCAlabelMean(gca, Ggca_label, means) ;
    printf("label %s: mean = ", cma_label_to_name(Ggca_label)) ;
    for (i = 0 ; i < gca->ninputs ; i++)
      printf("%2.1f ", means[i]) ;
    printf("\n") ;
  }

  gca_face = GCAread(gca_face_fname) ;
  if (gca_face == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open face GCA %s.\n",
              Progname, gca_face_fname) ;

  if (novar)
    GCAunifyVariance(gca) ;

  if (renormalization_fname) {
    FILE   *fp ;
    int    *labels, nlines, i ;
    float  *intensities, f1, f2 ;
    char   *cp, line[STRLEN] ;

    fp = fopen(renormalization_fname, "r") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not read %s",
                Progname, renormalization_fname) ;

    cp = fgetl(line, 199, fp) ;
    nlines = 0 ;
    while (cp) {
      nlines++ ;
      cp = fgetl(line, 199, fp) ;
    }
    rewind(fp) ;
    printf("reading %d labels from %s...\n", nlines,renormalization_fname) ;
    labels = (int *)calloc(nlines, sizeof(int)) ;
    intensities = (float *)calloc(nlines, sizeof(float)) ;
    cp = fgetl(line, 199, fp) ;
    for (i = 0 ; i < nlines ; i++) {
      sscanf(cp, "%e  %e", &f1, &f2) ;
      labels[i] = (int)f1 ;
      intensities[i] = f2 ;
      if (labels[i] == Left_Cerebral_White_Matter)
        DiagBreak() ;
      cp = fgetl(line, 199, fp) ;
    }
    GCArenormalizeIntensities(gca, labels, intensities, nlines) ;
    free(labels) ;
    free(intensities) ;
  }

  printf("reading '%s'...\n", in_fname) ;
  fflush(stdout) ;
  mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
              Progname, in_fname) ;
  mri_orig = MRIcopy(mri_in, NULL) ;
  if (mri_in->type!=MRI_UCHAR) {
    MRI *mri_tmp ;

    printf("changing type of input volume to 8 bits/voxel...\n") ;
    mri_tmp = MRIconform(mri_in) ;
    mri_in = mri_tmp ;
  }
  if (parms.write_iterations != 0) {
    char fname[STRLEN] ;
    MRI  *mri_gca ;
    mri_gca = MRIclone(mri_in, NULL) ;
    GCAbuildMostLikelyVolume(gca, mri_gca) ;
    sprintf(fname, "%s_target", parms.base_name) ;
    MRIwriteImageViews(mri_gca, fname, IMAGE_SIZE) ;
    sprintf(fname, "%s_target.mgh", parms.base_name) ;
    printf("writing target volume to %s...\n", fname) ;
    MRIwrite(mri_gca, fname) ;
  }
  if (mask_fname) {
    MRI *mri_mask ;

    mri_mask = MRIread(mask_fname) ;
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                Progname, mask_fname) ;

    MRImask(mri_in, mri_mask, mri_in, 0, 0) ;
    MRIfree(&mri_mask) ;
  }

  if (alpha > 0)
    mri_in->flip_angle = alpha ;
  if (TR > 0)
    mri_in->tr = TR ;
  if (TE > 0)
    mri_in->te = TE ;

  if (example_T1) {
    MRI *mri_T1, *mri_seg ;

    mri_seg = MRIread(example_segmentation) ;
    if (!mri_seg)
      ErrorExit
      (ERROR_NOFILE,"%s: could not read example segmentation from %s",
       Progname, example_segmentation) ;
    mri_T1 = MRIread(example_T1) ;
    if (!mri_T1)
      ErrorExit(ERROR_NOFILE,"%s: could not read example T1 from %s",
                Progname, example_T1) ;
    printf("scaling atlas intensities using specified examples...\n") ;
    MRIeraseBorderPlanes(mri_seg, 1) ;
    GCArenormalizeToExample(gca, mri_seg, mri_T1) ;
    MRIfree(&mri_seg) ;
    MRIfree(&mri_T1) ;
  }

  if (tissue_parms_fname)   /* use FLASH forward model */
    GCArenormalizeToFlash(gca, tissue_parms_fname, mri_in) ;

  if (!FZERO(tx) || !FZERO(ty) || !FZERO(tz)) {
    MRI *mri_tmp ;

    printf("translating second volume by (%2.1f, %2.1f, %2.1f)\n",
           tx, ty, tz) ;
    mri_tmp = MRItranslate(mri_in, NULL, tx, ty, tz) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }

  if (!transform_loaded)   /* wasn't preloaded */
  {
    parms.transform = transform = TransformAlloc(LINEAR_VOX_TO_VOX, mri_in) ;
    Glta = parms.lta = (LTA *)transform->xform ;
  }

  if (!FZERO(blur_sigma)) {
    MRI *mri_tmp, *mri_kernel ;

    mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
    mri_tmp = MRIconvolveGaussian(mri_in, NULL, mri_kernel) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }

  parms.lta->xforms[0].m_L = MatrixIdentity(4, NULL) ;
  if (parms.write_iterations != 0) {
    char fname[STRLEN] ;
    MRI *mri_aligned ;

    mri_aligned =
      MRIlinearTransform(mri_in, NULL, parms.lta->xforms[0].m_L) ;
    sprintf(fname, "%s_conformed.mgh", parms.base_name) ;
    MRIwrite(mri_aligned, fname) ;
    sprintf(fname, "%s_conformed", parms.base_name) ;
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
    MRIfree(&mri_aligned) ;
  }

  i = 0 ;
#if 0
  if (gca_reduced != gca)
    GCAfree(&gca_reduced) ;
#endif
  for (spacing = max_spacing, scale = 0 ;
       scale < nscales ;
       scale++, spacing /= 2) {
    if (use_contrast)
      parms.gcas = GCAfindContrastSamples(gca,&nsamples, spacing,min_prior);
    else
      parms.gcas = GCAfindStableSamples
	(gca, &nsamples,spacing, min_prior, exclude_list,1,-1);
    printf("spacing=%d, using %d sample points, tol=%2.2e...\n",
           spacing, nsamples, parms.tol) ;
    parms.nsamples = nsamples ;
    if (sample_fname) {
      GCAwriteSamples(gca, mri_in, parms.gcas, nsamples, sample_fname) ;
      printf("samples written\n") ;
    }
    old_log_p = GCAcomputeLogSampleProbability
                (gca, parms.gcas, mri_in, transform, nsamples, DEFAULT_CLAMP) ;
    register_mri(mri_in, gca, &parms,i, spacing) ;
    log_p = GCAcomputeLogSampleProbability
            (gca, parms.gcas, mri_in, transform, nsamples, DEFAULT_CLAMP) ;

    printf("pass %d, spacing %d: log(p) = %2.1f (old=%2.1f)\n",
           i+1, spacing, log_p, old_log_p) ;
    free(parms.gcas) ;
    parms.tol *= 10 ;
    i++ ;
  }

  parms.gcas = GCAfindAllSamples(gca, &nsamples, exclude_list,1) ;
  parms.nsamples = nsamples ;

  printf("computing final MAP estimate of "
         "linear transform using %d samples...\n", nsamples) ;

  parms.mri_in = mri_in ;  /* for diagnostics */
  if ((Gdiag & DIAG_WRITE) && (parms.write_iterations != 0)) {
    MRI *mri_aligned ;
    char fname[STRLEN] ;

    mri_aligned = MRIlinearTransform(mri_in, NULL, parms.lta->xforms[0].m_L);
    sprintf(fname, "%s_before_final_alignment.mgh", parms.base_name) ;
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
    MRIfree(&mri_aligned) ;
  }
  parms.tol = tol ;  /* otherwise NRC can get upset */
  MRIemAlign(mri_in, gca, &parms, parms.lta->xforms[0].m_L) ;

#if 0
  printf("final EM transform:\n") ;
  MatrixPrint(stdout, parms.lta->xforms[0].m_L) ;
  printf("\n") ;

#define DELTA_SCALE 0.01
  update_optimal_transform
  (mri_in, gca, parms.gcas, nsamples,
   parms.lta->xforms[0].m_L, DELTA_SCALE, parms.write_iterations, 3, 3) ;
  compare_transform
  (mri_in, gca, parms.gcas, nsamples, parms.lta->xforms[0].m_L) ;
  recompute_labels
  (mri_in, gca, parms.gcas, nsamples, parms.lta->xforms[0].m_L) ;
  compare_transform
  (mri_in, gca, parms.gcas, nsamples, parms.lta->xforms[0].m_L) ;
#endif
  printf("final transform:\n") ;
  MatrixPrint(stdout, parms.lta->xforms[0].m_L) ;
  printf("\n") ;

  if (xform_fname != NULL) {
    printf("writing transform to %s...\n", xform_fname) ;
    LTAwrite(parms.lta, xform_fname) ;
  }

  if ((Gdiag & DIAG_WRITE) && (parms.write_iterations != 0)) {
    MRI *mri_aligned ;
    char fname[STRLEN] ;

    mri_aligned = MRIlinearTransform(mri_in, NULL, parms.lta->xforms[0].m_L);
    sprintf(fname, "%s_after_final_alignment", parms.base_name) ;
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
    MRIfree(&mri_aligned) ;
  }

  if (transformed_sample_fname) {
    GCAtransformAndWriteSamples(gca, mri_in, parms.gcas, nsamples,
                                transformed_sample_fname, transform) ;
    printf("samples written\n") ;
  }
  printf("anonymizing volume...\n") ;
  MRIreplaceValues(mri_in, mri_in, 255, 254) ;
  mri_out =
    MRIremoveFace(mri_in, NULL, parms.lta, gca, gca_face, radius, 255) ;
  {
    MRI *mri_tmp ;

    printf("resampling to original coordinate system...\n");
    mri_tmp = MRIresample(mri_out, mri_orig, SAMPLE_NEAREST) ;
    MRIcopyHeader(mri_orig, mri_tmp) ;
    MRIfree(&mri_out) ;
    mri_out = mri_tmp ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_out, "after_resampling.mgh") ;
  }

  if (mean_flag > 0)
    MRImeanMask(mri_orig, mri_out, mri_orig, 255, mean_flag) ;
  else
    MRImask(mri_orig, mri_out, mri_orig, 255, fill_val) ;
  printf("writing anonymized volume to %s...\n", out_fname) ;
  MRIwrite(mri_orig, out_fname) ;

#if 0
  if (gca)
    GCAfree(&gca) ;
#endif
  if (mri_in)
    MRIfree(&mri_in) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("deidentification took %d minutes and %d seconds.\n",
         minutes, seconds) ;
  if (diag_fp)
    fclose(diag_fp) ;
  exit(0) ;
  return(0) ;
}


static int
register_mri(MRI *mri_in,
             GCA *gca,
             MORPH_PARMS *parms,
             int passno, int spacing) {
  MATRIX  *m_L ;

#if 0
  MRIscaleMeanIntensities(mri_in, mri_ref, mri_in);
  printf("initializing alignment using PCA...\n") ;
#endif


#if 0
  if (passno == 0)
    m_L = MatrixIdentity(4, NULL) ;
  else
#endif
    m_L = MatrixCopy(parms->lta->xforms[0].m_L, NULL) ;

  find_optimal_transform(mri_in, gca, parms->gcas, parms->nsamples,m_L,passno,
                         parms->write_iterations, spacing);

  /* make sure transform and lta are the same (sorry - retrofitting!) */
  if (!parms->lta) {
    parms->transform = transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
    Glta = parms->lta = (LTA *)transform->xform ;
  } else {
    parms->transform = transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
    transform->xform = (void *)parms->lta ;
  }

  MatrixCopy(m_L, parms->lta->xforms[0].m_L) ;
  if (Gdiag & DIAG_SHOW) {
    printf("global search transform:\n") ;
    MatrixPrint(stdout, m_L) ;
  }

  /*    parms->start_t++ ;*/
  printf("computing MAP estimate of linear transform...\n") ;

  parms->mri_in = mri_in ;  /* for diagnostics */
  MRIemAlign(mri_in, gca, parms, m_L) ;

#if 0
  printf("final transform:\n") ;
  MatrixPrint(stdout, parms->lta->xforms[0].m_L) ;
  printf("\n") ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRI *mri_aligned ;
    char fname[STRLEN] ;

    mri_aligned =
      MRIlinearTransform(mri_in, NULL, parms->lta->xforms[0].m_L) ;
    sprintf(fname, "%s_after_alignment", parms->base_name) ;
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
    MRIfree(&mri_aligned) ;
  }
#endif

  return(NO_ERROR) ;
}

#define DEFAULT_MAX_STEPS 5
static double MAX_ANGLES = DEFAULT_MAX_STEPS ;
#define MAX_ANGLE       RADIANS(15)
#define MIN_ANGLE       RADIANS(2)

#define MAX_SCALE       2.0
#define MIN_SCALE       0.5

static int max_angles = DEFAULT_MAX_STEPS ;


static int MAX_TRANS_STEPS = DEFAULT_MAX_STEPS ;
static double MAX_TRANS = 30 ;

static MATRIX *
find_optimal_transform(MRI *mri,
                       GCA *gca,
                       GCA_SAMPLE *gcas,
                       int nsamples,
                       MATRIX *m_L,
                       int passno,
                       int write_iterations,
                       int spacing) {
  MATRIX   *m_origin ;
  MRI      *mri_gca  ;
  double   gca_means[3], max_log_p, old_max,
  max_angle, angle_steps, min_scale, max_scale, scale_steps, scale,
  delta, mean ;
  int      niter, good_step, done, nscales, scale_samples ;
  float    min_search_scale ;
#if  0
  float      min_real_val, fmax, fmin ;
  int        min_real_bin, mri_peak ;
  double     in_means[3], dx, dy, dz ;
  MRI_REGION box, gca_box ;
  HISTOGRAM *h_mri, *h_smooth ;
#endif

#define MIN_SEARCH_SCALE 0.1
  min_search_scale = MIN_SEARCH_SCALE ;

  if (spacing >= 16)
    scale_samples = 5 ;
  else
    scale_samples = 3 ;

  if (spacing >= 8)
    min_search_scale /= 4;

  max_log_p =
    local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples, DEFAULT_CLAMP) ;

  mri_gca = MRIclone(mri, NULL) ;
  GCAmri(gca, mri_gca) ;
  MRIcenterOfMass(mri_gca, gca_means, 0) ;
  m_origin = MatrixIdentity(4, NULL) ;
  *MATRIX_RELT(m_origin, 1, 4) = gca_means[0]*(float)center ;
  *MATRIX_RELT(m_origin, 2, 4) = gca_means[1]*(float)center ;
  *MATRIX_RELT(m_origin, 3, 4) = gca_means[2]*(float)center ;
  *MATRIX_RELT(m_origin, 4, 4) = 1 ;

  if (passno == 0)   /* only first time*/
  {
    if (Gdiag & DIAG_WRITE && write_iterations > 0) {
      char fname[STRLEN] ;

      sprintf(fname, "%s_before_intensity.mgh", parms.base_name) ;
      printf("writing snapshot to %s...\n", fname) ;
      MRIwrite(mri, fname) ;
      /*    MRIwriteImageViews(mri, "before_intensity", IMAGE_SIZE) ;*/
    }

    if (!noiscale)
      GCAhistoScaleImageIntensities(gca, mri, 0) ;

    if (Gdiag & DIAG_WRITE && write_iterations > 0) {
      char fname[STRLEN] ;

      Glta->xforms[0].m_L = m_L ;
      sprintf(fname, "%s_after_intensity.mgh", parms.base_name) ;
      printf("writing snapshot to %s...\n", fname) ;
      MRIwrite(mri, fname) ;
      sprintf(fname, "%s_centering0", parms.base_name) ;
      MRIwriteImageViews(mri, fname, IMAGE_SIZE) ;
      sprintf(fname, "%s_fsamples_centering0.mgh", parms.base_name) ;
      GCAtransformAndWriteSamples(gca, mri, gcas, nsamples,
                                  fname, transform) ;
    }
    /* first align centroids */
    if (gca_mean_fname) {
      printf("writing gca volume to %s...\n", gca_mean_fname) ;
      MRIwrite(mri_gca, gca_mean_fname) ;
      printf("done\n") ;
    }

    printf("initial log_p = %2.1f\n", max_log_p) ;

#if 0
    MRIvalRange(mri, &fmin, &fmax) ;
    h_mri = MRIhistogram(mri, nint(fmax-fmin+1)) ;
    h_mri->counts[0] = 0 ; /* ignore background */
    h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
    mri_peak = HISTOfindHighestPeakInRegion(h_smooth, 0, h_smooth->nbins/3) ;
    min_real_bin = HISTOfindEndOfPeak(h_smooth, mri_peak, .25) ;
    min_real_val = h_smooth->bins[min_real_bin] ;
    printf("using real data threshold=%2.1f\n", min_real_val) ;
    MRIfindApproximateSkullBoundingBox(mri, min_real_val, &box) ;
    HISTOfree(&h_mri) ;
    HISTOfree(&h_smooth) ;
    printf("input bounding box (%d, %d, %d) --> (%d, %d, %d)\n",
           box.x, box.y, box.z, box.x+box.dx, box.y+box.dy, box.z+box.dz) ;

    MRIvalRange(mri_gca, &fmin, &fmax) ;
    h_mri = MRIhistogram(mri_gca, nint(fmax-fmin+1)) ;
    h_mri->counts[0] = 0 ; /* ignore background */
    h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
    mri_peak = HISTOfindHighestPeakInRegion(h_smooth, 0, h_smooth->nbins/3) ;
    min_real_bin = HISTOfindEndOfPeak(h_smooth, mri_peak, .25) ;
    min_real_val = h_smooth->bins[min_real_bin] ;
    printf("using GCA real data threshold=%2.1f\n", min_real_val) ;
    MRIfindApproximateSkullBoundingBox(mri_gca, min_real_val, &gca_box) ;
    HISTOfree(&h_mri) ;
    HISTOfree(&h_smooth) ;
    printf("gca bounding box (%d, %d, %d) --> (%d, %d, %d)\n",
           box.x, box.y, box.z, box.x+box.dx, box.y+box.dy, box.z+box.dz) ;

    /*                MRIcenterOfMass(mri, in_means, 0) ;*/
    in_means[0] = box.x+box.dx*0.5 ;
    in_means[2] = box.z+box.dz*0.5 ;
    printf("input centroid (%2.1f, %2.1f, %2.1f), "
           "gca centroid (%2.1f, %2.1f, %2.1f)\n",
           in_means[0], in_means[1], in_means[2],
           gca_means[0], gca_means[1], gca_means[2]) ;
    in_means[1] = box.y+box.dx*0.55 ;
    printf("resetting superior/inferior centroid to %2.1f\n", in_means[1]) ;

    /* now apply translation to take in centroid to ref centroid */
    dx = gca_means[0] - in_means[0] ;
    /* use top of skull as estimate of
       superior-inferior offset (about 1.5cm from brain) */
    dy = gca_box.y-(box.y+15*mri->ysize) ;
    dz = gca_means[2] - in_means[2] ;
    if (passno == 0) {
      *MATRIX_RELT(m_L, 1, 4) = dx ;
      *MATRIX_RELT(m_L, 2, 4) = dy ;
      *MATRIX_RELT(m_L, 3, 4) = dz ;
    }

    max_log_p =
      local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples, DEFAULT_CLAMP) ;
    printf("initial translation: (%2.1f, %2.1f, %2.1f): log p = %2.1f\n",
           dx,dy,dz, max_log_p) ;
#else
    max_log_p =
      find_optimal_translation
      (gca, gcas, mri, nsamples, m_L, -100, 100, 11, 3) ;
    max_log_p = local_GCAcomputeLogSampleProbability
                (gca, gcas, mri, m_L,nsamples, DEFAULT_CLAMP) ;
    printf("after initial translation: (%2.1f, %2.1f, %2.1f): "
           "log p = %2.1f\n",
           *MATRIX_RELT(m_L, 1, 4),
           *MATRIX_RELT(m_L, 2, 4),
           *MATRIX_RELT(m_L, 3, 4),
           max_log_p) ;
#endif

    if (write_iterations != 0) {
      char fname[STRLEN] ;
      MRI *mri_aligned ;

      mri_aligned = MRIlinearTransform(mri, NULL, m_L) ;
      sprintf(fname, "%s_centering1", parms.base_name) ;
      MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
      sprintf(fname, "%s_after_centering.mgh", parms.base_name) ;
      printf("writing image after centering to %s...\n", fname) ;
#if 0
      MRIwrite(mri_aligned, fname) ;
#else
      Glta->xforms[0].m_L = m_L ;
      sprintf(fname, "%s_fsamples_centering1.mgh", parms.base_name) ;
      printf("writing samples after centering to %s...\n", fname) ;
      GCAtransformAndWriteSamples(gca, mri, gcas, nsamples,
                                  fname, transform) ;
#endif
      MRIfree(&mri_aligned) ;
    }
    MRIfree(&mri_gca) ;
  }

  max_angle = MAX_ANGLE ;
  angle_steps = max_angles ;
  max_scale = MAX_SCALE ;
  min_scale = MIN_SCALE ;
  scale_steps = max_angles ;

#define MIN_SCALES 2
  niter = nscales = 0 ;
  scale = 1.0 ;
  good_step = 0 ;
  done = 0 ;
  do {
    old_max = max_log_p ;

    max_log_p = find_optimal_linear_xform
                (gca, gcas, mri, nsamples,
                 m_L, m_origin,
                 -RADIANS(2*spacing*scale),
                 RADIANS(2*spacing*scale),
                 1-.25*(spacing/16.0)*scale, 1+.25*(spacing/16.0)*scale,
                 -scale*(spacing/16.0)*MAX_TRANS, scale*(spacing/16.0)*MAX_TRANS,
                 scale_samples, 3, 3, 2);

#if 0
    /* it's more likely that there will be large-scale scaling
       than large-scale rotations, so look for scaling first */
    max_log_p = find_optimal_scaling_and_rotation
                (gca, gcas, mri, nsamples,
                 m_L, m_origin,
                 -RADIANS(5*scale),
                 RADIANS(5*scale),
                 1-.1*scale, 1+.1*scale,
                 max_angles/2,max_angles/2,3);
#if 1
    max_log_p = find_optimal_scaling(gca, gcas, mri, nsamples, m_L, m_origin,
                                     min_scale, max_scale, scale_steps, 3) ;
    max_log_p = find_optimal_rotation(gca, gcas, mri, nsamples, m_L,m_origin,
                                      -scale*max_angle, scale*max_angle,
                                      angle_steps, 3) ;
#endif
    max_log_p = find_optimal_translation(gca, gcas, mri, nsamples, m_L,
                                         -scale*MAX_TRANS, scale*MAX_TRANS,
                                         MAX_TRANS_STEPS, 2) ;
#endif

    if (write_iterations != 0) {
      char fname[STRLEN] ;
      MRI *mri_aligned ;

      mri_aligned = MRIlinearTransform(mri, NULL, m_L) ;
      sprintf(fname, "%s%03d", parms.base_name, parms.start_t+niter+1) ;
      MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
      sprintf(fname, "%s%03d.mgh",
              parms.base_name, parms.start_t+niter+1) ;
#if 0
      MRIwrite(mri_aligned, fname) ;
#else
      Glta->xforms[0].m_L = m_L ;
      GCAtransformAndWriteSamples(gca, mri, gcas, nsamples,
                                  fname, transform) ;
#endif
      MRIfree(&mri_aligned) ;
    }
    printf("scale %2.3f: max=%2.1f, old_max =%2.1f (thresh=%2.1f)\n",
           scale,max_log_p, old_max, old_max+fabs(tol*old_max)) ;

    /* search a finer nbhd (if do-while continues) */
    if ((spacing < 8) ||
        (max_log_p < old_max-tol*old_max)) /* couldn't take a step */
    {
      if ((spacing < 8) || good_step) {
        scale *= 0.25 ;
        if (scale < min_search_scale)
          break ;
        mean = (max_scale + min_scale)/2 ;
        delta = (max_scale - min_scale)/2 ;
        max_scale = 1.0 + delta*scale ;
        min_scale = 1.0 - delta*scale ;
        good_step = 0 ;
        printf("reducing scale to %2.4f\n", scale) ;
        nscales++ ;
      } else
        done = 1 ;
    } else
      good_step = 1 ; /* took at least one good step at this scale */

    niter++ ;
  } while (nscales++ < MIN_SCALES || (done == FALSE)) ;

  parms.start_t += niter ;
  MatrixFree(&m_origin) ;
  return(m_L) ;
}

#if 0
static double
find_optimal_rotation(GCA *gca, GCA_SAMPLE *gcas, MRI *mri, int nsamples,
                      MATRIX *m_L, MATRIX *m_origin,
                      float min_angle, float max_angle, float angle_steps,
                      int nreductions) {
  MATRIX   *m_rot, *m_x_rot, *m_y_rot, *m_z_rot, *m_tmp,*m_L_tmp,*m_origin_inv,
  *m_tmp2 ;
  double   x_angle, y_angle, z_angle, x_max, y_max, z_max, delta,
  log_p, max_log_p, mean_angle ;
  int      i ;

  m_origin_inv = MatrixCopy(m_origin, NULL) ;
  *MATRIX_RELT(m_origin_inv, 1, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 2, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 3, 4) *= -1 ;
  m_L_tmp = m_x_rot = m_y_rot = m_z_rot = m_rot = m_tmp = m_tmp2 = NULL ;
  x_max = y_max = z_max = 0.0 ;
  max_log_p =
    local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples, DEFAULT_CLAMP) ;
  for (i = 0 ; i <= nreductions ; i++) {
    delta = (max_angle-min_angle) / angle_steps ;
    if (FZERO(delta))
      return(max_log_p) ;

    if (Gdiag & DIAG_SHOW)
      printf(
        "scanning %2.2f degree nbhd (%2.1f) ",
        (float)DEGREES(max_angle), (float)DEGREES(delta)) ;

    for (x_angle = min_angle ; x_angle <= max_angle ; x_angle += delta) {
      m_x_rot = MatrixReallocRotation(4, x_angle, X_ROTATION, m_x_rot) ;
      for (y_angle = min_angle ; y_angle <= max_angle ; y_angle += delta) {
        m_y_rot =
          MatrixReallocRotation(4, y_angle, Y_ROTATION, m_y_rot) ;
        m_tmp =
          MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
        for (z_angle= min_angle ;
             z_angle <= max_angle ;
             z_angle += delta) {
          m_z_rot =
            MatrixReallocRotation(4, z_angle, Z_ROTATION, m_z_rot) ;
          m_rot = MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
          m_tmp2 = MatrixMultiply(m_rot, m_origin_inv, m_tmp2) ;
          MatrixMultiply(m_origin, m_tmp2, m_rot) ;

          m_L_tmp = MatrixMultiply(m_rot, m_L, m_L_tmp) ;
          log_p =
            local_GCAcomputeLogSampleProbability
            (gca, gcas, mri, m_L_tmp,nsamples, DEFAULT_CLAMP) ;
          if (log_p > max_log_p) {
            max_log_p = log_p ;
            x_max = x_angle ;
            y_max = y_angle ;
            z_max = z_angle ;
#if 0
            printf("new max p %2.1f found at "
                   "(%2.1f, %2.1f, %2.1f)\n",
                   max_log_p,
                   (float)DEGREES(x_angle),
                   (float)DEGREES(y_angle),
                   (float)DEGREES(z_angle)) ;
#endif
          }
        }
      }

    }

    if (Gdiag & DIAG_SHOW)
      printf(
        "max log p = %2.1f @ (%2.1f, %2.1f, %2.1f)\n",
        (float)max_log_p, (float)DEGREES(x_max), (float)DEGREES(y_max),
        (float)DEGREES(z_max)) ;

    /* update L to reflect new maximum and search around it */
    MatrixReallocRotation(4, x_max, X_ROTATION, m_x_rot) ;
    MatrixReallocRotation(4, y_max, Y_ROTATION, m_y_rot) ;
    MatrixReallocRotation(4, z_max, Z_ROTATION, m_z_rot) ;
    MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
    MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
    m_tmp2 = MatrixMultiply(m_rot, m_origin_inv, m_tmp2) ;
    MatrixMultiply(m_origin, m_tmp2, m_rot) ;
    MatrixMultiply(m_rot, m_L, m_L) ;

    x_max = y_max = z_max = 0.0 ;  /* we've rotated transform to old max */

    mean_angle = (max_angle + min_angle) / 2 ;
    delta = (max_angle-min_angle)/4 ;
    min_angle = mean_angle - delta ;
    max_angle = mean_angle + delta ;
#if 0
    delta = (max_angle-min_angle) / angle_steps ;
    min_angle -= delta/2 ;
    max_angle += delta/2 ;
#endif
  }

  printf("\n") ;

  MatrixFree(&m_x_rot) ;
  MatrixFree(&m_y_rot) ;
  MatrixFree(&m_z_rot) ;
  MatrixFree(&m_rot) ;
  MatrixFree(&m_tmp) ;
  MatrixFree(&m_origin_inv) ;
  MatrixFree(&m_tmp2) ;
  return(max_log_p) ;
}

static double
find_optimal_scaling(GCA *gca, GCA_SAMPLE *gcas, MRI *mri, int nsamples,
                     MATRIX *m_L, MATRIX *m_origin,
                     float min_scale, float max_scale, float scale_steps,
                     int nreductions) {
  MATRIX   *m_scale, *m_tmp,*m_L_tmp,*m_origin_inv ;
  double   x_scale, y_scale, z_scale, x_max, y_max, z_max, delta,
  log_p, max_log_p, mean_scale ;
  int      i ;

  m_origin_inv = MatrixCopy(m_origin, NULL) ;
  *MATRIX_RELT(m_origin_inv, 1, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 2, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 3, 4) *= -1 ;
  m_scale = MatrixIdentity(4, NULL) ;
  m_L_tmp = m_tmp = NULL ;
  x_max = y_max = z_max = 1.0 ;
  max_log_p =
    local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples, DEFAULT_CLAMP) ;
  for (i = 0 ; i <= nreductions ; i++) {
    delta = (max_scale-min_scale) / scale_steps ;
    if (FZERO(delta))
      return(max_log_p) ;
    if (Gdiag & DIAG_SHOW) {
      printf("scanning scales %2.3f->%2.3f (%2.3f) ",
             min_scale,max_scale, delta) ;
      fflush(stdout) ;
    }
    for (x_scale = min_scale ; x_scale <= max_scale ; x_scale += delta) {
      *MATRIX_RELT(m_scale, 1, 1) = x_scale ;
      for (y_scale = min_scale ; y_scale <= max_scale ; y_scale += delta) {
        *MATRIX_RELT(m_scale, 2, 2) = y_scale ;
        for (z_scale= min_scale ;
             z_scale <= max_scale ;
             z_scale += delta) {
          *MATRIX_RELT(m_scale, 3, 3) = z_scale ;

          /* reset translation values */
          *MATRIX_RELT(m_scale, 1, 4) =
            *MATRIX_RELT(m_scale, 2, 4) =
              *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
          m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
          MatrixMultiply(m_origin, m_tmp, m_scale) ;

          m_L_tmp = MatrixMultiply(m_scale, m_L, m_L_tmp) ;
          log_p =
            local_GCAcomputeLogSampleProbability
            (gca, gcas, mri, m_L_tmp,nsamples, DEFAULT_CLAMP) ;
          if (log_p > max_log_p) {
            max_log_p = log_p ;
            x_max = x_scale ;
            y_max = y_scale ;
            z_max = z_scale ;
#if 0
            printf("new max p %2.1f found at "
                   "(%2.3f, %2.3f, %2.3f)\n",
                   max_log_p, x_scale, y_scale, z_scale) ;
#endif
          }
        }
      }

    }

    if (Gdiag & DIAG_SHOW)
      printf("max log p = %2.1f @ (%2.3f, %2.3f, %2.3f)\n",
             max_log_p, x_max, y_max, z_max) ;

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_scale, 1, 4) =
      *MATRIX_RELT(m_scale, 2, 4) = *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
    *MATRIX_RELT(m_scale,1,1) = x_max ;
    *MATRIX_RELT(m_scale,2,2) = y_max ;
    *MATRIX_RELT(m_scale,3,3) = z_max ;
    m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
    MatrixMultiply(m_origin, m_tmp, m_scale) ;
    MatrixMultiply(m_scale, m_L, m_L_tmp) ;
    MatrixCopy(m_L_tmp, m_L) ;

    x_max = y_max = z_max = 1.0 ;  /* we've scaled transform by old maxs */

    mean_scale = (max_scale + min_scale) / 2 ;
    delta = (max_scale-min_scale)/4 ;
    min_scale = mean_scale - delta ;
    max_scale = mean_scale + delta ;
  }

  printf("\n") ;

  MatrixFree(&m_scale) ;
  MatrixFree(&m_tmp) ;
  MatrixFree(&m_origin_inv) ;
  return(max_log_p) ;
}

#endif
static double
find_optimal_translation(GCA *gca, GCA_SAMPLE *gcas, MRI *mri, int nsamples,
                         MATRIX *m_L, float min_trans, float max_trans,
                         float trans_steps, int nreductions) {
  MATRIX   *m_trans, *m_L_tmp ;
  double   x_trans, y_trans, z_trans, x_max, y_max, z_max, delta,
  log_p, max_log_p, mean_trans ;
  int      i ;

  delta = (max_trans-min_trans) / trans_steps ;
  m_L_tmp = NULL ;
  m_trans = MatrixIdentity(4, NULL) ;
  x_max = y_max = z_max = 0.0 ;
  max_log_p =
    local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples, DEFAULT_CLAMP) ;

  for (i = 0 ; i <= nreductions ; i++) {
    delta = (max_trans-min_trans) / trans_steps ;
    if (FZERO(delta))
      return(max_log_p) ;
    if (Gdiag & DIAG_SHOW) {
      printf(
        "scanning translations %2.2f->%2.2f (%2.1f) ",
        min_trans,max_trans, delta) ;
      fflush(stdout) ;
    }
    for (x_trans = min_trans ; x_trans <= max_trans ; x_trans += delta) {
      *MATRIX_RELT(m_trans, 1, 4) = x_trans ;
      for (y_trans = min_trans ; y_trans <= max_trans ; y_trans += delta) {
        *MATRIX_RELT(m_trans, 2, 4) = y_trans ;
        for (z_trans= min_trans ;
             z_trans <= max_trans ;
             z_trans += delta) {
          *MATRIX_RELT(m_trans, 3, 4) = z_trans ;
          m_L_tmp = MatrixMultiply(m_trans, m_L, m_L_tmp) ;
          log_p =
            local_GCAcomputeLogSampleProbability
            (gca, gcas, mri, m_L_tmp,nsamples, DEFAULT_CLAMP) ;
          if (log_p > max_log_p) {
            max_log_p = log_p ;
            x_max = x_trans ;
            y_max = y_trans ;
            z_max = z_trans ;
#if 0
            printf("new max p %2.1f found at "
                   "(%2.1f, %2.1f, %2.1f)\n",
                   max_log_p, x_trans, y_trans, z_trans) ;
#endif
          }
        }
      }
    }

    if (Gdiag & DIAG_SHOW)
      printf(
        "max log p = %2.1f @ (%2.1f, %2.1f, %2.1f)\n",
        max_log_p, x_max, y_max, z_max) ;

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_trans, 1, 4) = x_max ;
    *MATRIX_RELT(m_trans, 2, 4) = y_max ;
    *MATRIX_RELT(m_trans, 3, 4) = z_max ;
    MatrixMultiply(m_trans, m_L, m_L_tmp) ;
    MatrixCopy(m_L_tmp, m_L) ;
    x_max = y_max = z_max = 0.0 ;  /* we've translated
                                                        transform by old maxs */

    mean_trans = (max_trans + min_trans) / 2 ;
    delta = (max_trans-min_trans)/4 ;
    min_trans = mean_trans - delta ;
    max_trans = mean_trans + delta ;
  }

  printf("\n") ;

  MatrixFree(&m_trans) ;
  return(max_log_p) ;
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
  StrUpper(option) ;
  if (!stricmp(option, "DIST") || !stricmp(option, "DISTANCE")) {
    parms.l_dist = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_dist = %2.2f\n", parms.l_dist) ;
  } else if (!stricmp(option, "SAMPLES")) {
    sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing control points to %s...\n", sample_fname) ;
  } else if (!stricmp(option, "FILL")) {
    fill_val = atoi(argv[2]) ;
    nargs = 1 ;
    printf("filling defaced regions with %d\n", fill_val) ;
  } else if (!stricmp(option, "RADIUS")) {
    radius = atoi(argv[2]) ;
    nargs = 1 ;
    printf("erasing everything more than %d mm from possible brain\n",
           radius) ;
  } else if (!stricmp(option, "MASK")) {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  } else if (!stricmp(option, "MEAN")) {
    mean_flag = atoi(argv[2]) ;
    nargs = 1 ;
    printf("replacing values with local %dx%dx%d mean\n",
           mean_flag, mean_flag, mean_flag) ;
  } else if (!stricmp(option, "DEBUG_LABEL")) {
    Ggca_label = atoi(argv[2]) ;
    printf("debugging label %s (%d)\n",
           cma_label_to_name(Ggca_label), Ggca_label) ;
    nargs = 1 ;
  } else if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
  } else if (!stricmp(option, "DIAG")) {
    diag_fp = fopen(argv[2], "w") ;
    if (!diag_fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open diag file %s for writing",
                Progname, argv[2]) ;
    printf("opening diag file %s for writing\n", argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "TR")) {
    TR = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TR=%2.1f msec\n", TR) ;
  } else if (!stricmp(option, "EXAMPLE")) {
    example_T1 = argv[2] ;
    example_segmentation = argv[3] ;
    printf("using %s and %s as example T1 and segmentations respectively.\n",
           example_T1, example_segmentation) ;
    nargs = 2 ;
  } else if (!stricmp(option, "TE")) {
    TE = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TE=%2.1f msec\n", TE) ;
  } else if (!stricmp(option, "ALPHA")) {
    nargs = 1 ;
    alpha = RADIANS(atof(argv[2])) ;
    printf("using alpha=%2.0f degrees\n", DEGREES(alpha)) ;
  } else if (!stricmp(option, "FSAMPLES") || !stricmp(option, "ISAMPLES")) {
    transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing transformed control points to %s...\n",
           transformed_sample_fname) ;
  } else if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else if (!stricmp(option, "NSAMPLES")) {
    normalized_transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing  transformed normalization control points to %s...\n",
           normalized_transformed_sample_fname) ;
  } else if (!stricmp(option, "CONTRAST")) {
    use_contrast = 1 ;
    printf("using contrast to find labels...\n") ;
  } else if (!stricmp(option, "RENORM")) {
    renormalization_fname = argv[2] ;
    nargs = 1 ;
    printf("renormalizing using predicted intensity values in %s...\n",
           renormalization_fname) ;
  } else if (!stricmp(option, "FLASH")) {
    tissue_parms_fname = argv[2] ;
    nargs = 1 ;
    printf("using FLASH forward model and tissue parms in %s to predict"
           " intensity values...\n", tissue_parms_fname) ;
  } else if (!stricmp(option, "TRANSONLY")) {
    translation_only = 1 ;
    printf("only computing translation parameters...\n") ;
  } else if (!stricmp(option, "WRITE_MEAN")) {
    gca_mean_fname = argv[2] ;
    nargs = 1 ;
    printf("writing gca means to %s...\n", gca_mean_fname) ;
  } else if (!stricmp(option, "PRIOR")) {
    min_prior = atof(argv[2]) ;
    nargs = 1 ;
    printf("using prior threshold %2.2f\n", min_prior) ;
  } else if (!stricmp(option, "SPACING")) {
    max_spacing = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using max GCA spacing %d...\n", max_spacing) ;
  } else if (!stricmp(option, "NOVAR")) {
    novar = 1 ;
    printf("not using variance estimates\n") ;
  } else if (!stricmp(option, "DT")) {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    printf("dt = %2.2e\n", parms.dt) ;
  } else if (!stricmp(option, "TOL")) {
    tol = parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("tol = %2.2e\n", parms.tol) ;
  } else if (!stricmp(option, "CENTER")) {
    center = 1 ;
    printf("using GCA centroid as origin of transform\n") ;
  } else if (!stricmp(option, "NOSCALE")) {
    noscale = 1 ;
    printf("disabling scaling...\n") ;
  } else if (!stricmp(option, "NOISCALE")) {
    noiscale = 1 ;
    printf("disabling intensity scaling...\n") ;
  } else if (!stricmp(option, "NUM")) {
    num_xforms = atoi(argv[2]) ;
    nargs = 1 ;
    printf("finding a total of %d linear transforms\n", num_xforms) ;
  } else if (!stricmp(option, "AREA")) {
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_area = %2.2f\n", parms.l_area) ;
  } else if (!stricmp(option, "NLAREA")) {
    parms.l_nlarea = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_nlarea = %2.2f\n", parms.l_nlarea) ;
  } else if (!stricmp(option, "LEVELS")) {
    parms.levels = atoi(argv[2]) ;
    nargs = 1 ;
    printf("levels = %d\n", parms.levels) ;
  } else if (!stricmp(option, "INTENSITY") || !stricmp(option, "CORR")) {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_intensity = %2.2f\n", parms.l_intensity) ;
  } else if (!stricmp(option, "reduce")) {
    nreductions = atoi(argv[2]) ;
    nargs = 1 ;
    printf("reducing input images %d times before aligning...\n",
           nreductions) ;
  } else if (!stricmp(option, "nsamples")) {
    nsamples = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d samples of GCA...\n", nsamples) ;
  } else if (!stricmp(option, "XFORM")) {
    xform_fname = argv[2] ;
    nargs = 1 ;
    printf("writing transform to %s....\n", xform_fname) ;
  } else if (!stricmp(option, "norm")) {
    norm_fname = argv[2] ;
    nargs = 1 ;
    printf("intensity normalizing and writing to %s...\n",norm_fname);
  } else if (!stricmp(option, "scales") || !stricmp(option, "nscales")) {
    nscales = atoi(argv[2]) ;
    nargs = 1 ;
    printf("finding optimal linear transform over %d scales...\n", nscales);
  } else if (!stricmp(option, "steps")) {
    max_angles = atoi(argv[2]) ;
    nargs = 1 ;
    printf("taking %d angular steps...\n", max_angles) ;
  } else switch (*option) {
    case 'F':
      ctl_point_fname = argv[2] ;
      nargs = 1 ;
      printf("reading manually defined control points from %s\n",
             ctl_point_fname) ;
      break ;
    case 'D':
      tx = atof(argv[2]) ;
      ty = atof(argv[3]) ;
      tz = atof(argv[4]) ;
      nargs = 3 ;
      break ;
    case 'R':
      rxrot = RADIANS(atof(argv[2])) ;
      ryrot = RADIANS(atof(argv[3])) ;
      rzrot = RADIANS(atof(argv[4])) ;
      nargs = 3 ;
      break ;
    case 'T':
      parms.lta = LTAread(argv[2]) ;
      if (!parms.lta)
        ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s",
                  Progname, argv[2]) ;
      nargs = 1 ;
      printf("using previously computed transform %s\n", argv[2]) ;
      transform_loaded = 1 ;
      break ;
    case 'B':
      blur_sigma = atof(argv[2]) ;
      nargs = 1 ;
      printf("blurring input image with sigma=%2.3f\n", blur_sigma);
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'S':
#if 0
      parms.sigma = atof(argv[2]) ;
      printf("using sigma=%2.3f as upper bound on blurring.\n",
             parms.sigma) ;
      nargs = 1 ;
#else
      MAX_ANGLES = MAX_TRANS_STEPS = max_angles = (float)atoi(argv[2]) ;
      nargs = 1 ;
      printf("examining %2.0f different trans/rot/scale values...\n",
             MAX_ANGLES);
#endif
      break ;
    case '?':
    case 'H':
    case 'U':
      printf("usage: %s <in volume> <brain template> "
             "<face template> <defaced output volume>\n",
             argv[0]) ;
      exit(1) ;
      break ;
    case 'N':
      parms.niterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("niterations = %d\n", parms.niterations) ;
      break ;
    case 'W':
      parms.write_iterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("write iterations = %d\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'P':
      ctl_point_pct = atof(argv[2]) ;
      nargs = 1 ;
      printf("using top %2.1f%% wm points as control points....\n",
             100.0*ctl_point_pct) ;
      break ;
    case 'M':
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      printf("momentum = %2.2f\n", parms.momentum) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}
#if 0
static double
find_optimal_scaling_and_rotation(GCA *gca, GCA_SAMPLE *gcas,
                                  MRI *mri,
                                  int nsamples, MATRIX *m_L, MATRIX *m_origin,
                                  float min_angle, float max_angle,
                                  float min_scale, float max_scale,
                                  float angle_steps, float scale_steps,
                                  int nreductions) {
  MATRIX   *m_rot, *m_x_rot, *m_y_rot, *m_z_rot, *m_tmp,*m_L_tmp,*m_origin_inv,
  *m_tmp2, *m_scale;
  double   x_angle, y_angle, z_angle, x_max_rot, y_max_rot, z_max_rot, delta,
  x_max_scale, y_max_scale, z_max_scale, delta_scale,
  log_p, max_log_p, mean_angle, x_scale, y_scale, z_scale, mean_scale;
  int      i ;

  m_origin_inv = MatrixCopy(m_origin, NULL) ;
  *MATRIX_RELT(m_origin_inv, 1, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 2, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 3, 4) *= -1 ;
  m_L_tmp = m_x_rot = m_y_rot = m_z_rot = m_rot = m_tmp = m_tmp2 = NULL ;
  x_max_rot = y_max_rot = z_max_rot = 0.0 ;
  x_max_scale = y_max_scale = z_max_scale = 1.0f ;
  m_scale = MatrixIdentity(4, NULL) ;
  max_log_p =
    local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples, DEFAULT_CLAMP) ;
  for (i = 0 ; i <= nreductions ; i++) {
    delta_scale = (max_scale-min_scale) / scale_steps ;
    delta = (max_angle-min_angle) / angle_steps ;
    if (Gdiag & DIAG_SHOW) {
      printf("scanning %2.2f degree nbhd (%2.1f) and "
             "scale %2.3f->%2.3f (%2.3f)\n",
             (float)DEGREES(max_angle), (float)DEGREES(delta),
             min_scale,max_scale, delta_scale);
      fflush(stdout) ;
    }

    for (x_scale = min_scale ; x_scale <= max_scale ; x_scale += delta_scale) {
      /*      printf("x_scale = %2.3f\n", x_scale) ;*/
      *MATRIX_RELT(m_scale, 1, 1) = x_scale ;
      for (y_scale = min_scale ;
           y_scale <= max_scale ;
           y_scale += delta_scale) {
        *MATRIX_RELT(m_scale, 2, 2) = y_scale ;
        for (z_scale= min_scale ;
             z_scale <= max_scale;
             z_scale += delta_scale) {
          *MATRIX_RELT(m_scale, 3, 3) = z_scale ;

          /* reset translation values */
          *MATRIX_RELT(m_scale, 1, 4) =
            *MATRIX_RELT(m_scale, 2, 4) =
              *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
          m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
          MatrixMultiply(m_origin, m_tmp, m_scale) ;

          for (x_angle = min_angle ;
               x_angle <= max_angle ;
               x_angle += delta) {
            m_x_rot = MatrixReallocRotation
                      (4, x_angle, X_ROTATION, m_x_rot) ;
            for (y_angle = min_angle ;
                 y_angle <= max_angle ;
                 y_angle += delta) {
              m_y_rot = MatrixReallocRotation
                        (4, y_angle, Y_ROTATION, m_y_rot);
              m_tmp = MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
              for (z_angle= min_angle;
                   z_angle <= max_angle;
                   z_angle += delta) {
                m_z_rot = MatrixReallocRotation
                          (4, z_angle,Z_ROTATION,m_z_rot);
                m_rot = MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
                m_tmp2 = MatrixMultiply
                         (m_rot, m_origin_inv, m_tmp2) ;
                MatrixMultiply(m_origin, m_tmp2, m_rot) ;

                m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
                m_L_tmp = MatrixMultiply(m_tmp2, m_L, m_L_tmp) ;
                log_p =
                  local_GCAcomputeLogSampleProbability
                  (gca, gcas, mri, m_L_tmp, nsamples, DEFAULT_CLAMP) ;
                if (log_p > max_log_p) {
                  max_log_p = log_p ;
                  x_max_scale = x_scale ;
                  y_max_scale = y_scale ;
                  z_max_scale = z_scale ;
                  x_max_rot = x_angle ;
                  y_max_rot = y_angle ;
                  z_max_rot = z_angle ;
                }
              }
            }
          }
        }
      }
    }

    if (Gdiag & DIAG_SHOW) {
      printf("\tmax log p = %2.1f @ R=(%2.3f, %2.3f, %2.3f),"
             "S=(%2.3f,%2.3f,%2.3f)\n",
             max_log_p, DEGREES(x_max_rot), DEGREES(y_max_rot),
             DEGREES(z_max_rot),x_max_scale, y_max_scale, z_max_scale) ;
    }

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_scale, 1, 4) =
      *MATRIX_RELT(m_scale, 2, 4) = *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
    *MATRIX_RELT(m_scale,1,1) = x_max_scale ;
    *MATRIX_RELT(m_scale,2,2) = y_max_scale ;
    *MATRIX_RELT(m_scale,3,3) = z_max_scale ;
    m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
    MatrixMultiply(m_origin, m_tmp, m_scale) ;

    x_max_scale = y_max_scale = z_max_scale = 1.0 ;

    mean_scale = (max_scale + min_scale) / 2 ;
    delta_scale = (max_scale-min_scale)/4 ;
    min_scale = mean_scale - delta ;
    max_scale = mean_scale + delta ;

    /* update L to reflect new maximum and search around it */
    MatrixReallocRotation(4, x_max_rot, X_ROTATION, m_x_rot) ;
    MatrixReallocRotation(4, y_max_rot, Y_ROTATION, m_y_rot) ;
    MatrixReallocRotation(4, z_max_rot, Z_ROTATION, m_z_rot) ;
    MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
    MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
    m_tmp2 = MatrixMultiply(m_rot, m_origin_inv, m_tmp2) ;
    MatrixMultiply(m_origin, m_tmp2, m_rot) ;

    m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
    m_L_tmp = MatrixMultiply(m_tmp2, m_L, m_L_tmp) ;
    MatrixCopy(m_L_tmp, m_L) ;

    /* we've rotated transform to old max */
    x_max_rot = y_max_rot = z_max_rot = 0.0 ;

    mean_angle = (max_angle + min_angle) / 2 ;
    delta = (max_angle-min_angle)/4 ;
    min_angle = mean_angle - delta ;
    max_angle = mean_angle + delta ;
#if 0
    delta = (max_angle-min_angle) / angle_steps ;
    min_angle -= delta/2 ;
    max_angle += delta/2 ;
#endif
  }

  printf("\n") ;

  MatrixFree(&m_x_rot) ;
  MatrixFree(&m_y_rot) ;
  MatrixFree(&m_z_rot) ;
  MatrixFree(&m_rot) ;
  MatrixFree(&m_tmp) ;
  MatrixFree(&m_origin_inv) ;
  MatrixFree(&m_tmp2) ;
  return(max_log_p) ;
}
#endif


/*
  search 9-dimensional parameter space
*/
static double
find_optimal_linear_xform(GCA *gca, GCA_SAMPLE *gcas,
                          MRI *mri,
                          int nsamples, MATRIX *m_L, MATRIX *m_origin,
                          float min_angle, float max_angle,
                          float min_scale, float max_scale,
                          float min_trans, float max_trans,
                          float angle_steps, float scale_steps,
                          float trans_steps,
                          int nreductions) {
  MATRIX   *m_rot, *m_x_rot, *m_y_rot, *m_z_rot, *m_tmp,*m_L_tmp,*m_origin_inv,
  *m_tmp2, *m_scale, *m_trans, *m_tmp3 = NULL ;
  double   x_angle, y_angle, z_angle,
  x_max_rot, y_max_rot, z_max_rot, delta_rot,
  x_max_scale, y_max_scale, z_max_scale, delta_scale,
  x_trans, delta_trans, y_trans, z_trans,
  log_p, max_log_p, mean_angle, x_scale, y_scale, z_scale, mean_scale,
  x_max_trans, y_max_trans, z_max_trans, mean_trans ;
  int      i ;

  m_trans = MatrixIdentity(4, NULL) ;
  m_origin_inv = MatrixCopy(m_origin, NULL) ;
  *MATRIX_RELT(m_origin_inv, 1, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 2, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 3, 4) *= -1 ;
  m_L_tmp = m_x_rot = m_y_rot = m_z_rot = m_rot = m_tmp = m_tmp2 = NULL ;
  x_max_trans = y_max_trans = z_max_trans =
                                x_max_rot = y_max_rot = z_max_rot = 0.0 ;
  x_max_scale = y_max_scale = z_max_scale = 1.0f ;
  m_scale = MatrixIdentity(4, NULL) ;
  max_log_p = local_GCAcomputeLogSampleProbability
              (gca, gcas, mri, m_L,nsamples, DEFAULT_CLAMP) ;
  for (i = 0 ; i < nreductions ; i++) {
    delta_trans = (max_trans-min_trans) / (trans_steps-1) ;
    delta_scale = (max_scale-min_scale) / (scale_steps-1) ;
    delta_rot = (max_angle-min_angle) / (angle_steps-1) ;
    if (Gdiag & DIAG_SHOW) {
      printf("scanning %2.2f degree nbhd (%2.1f), "
             "scale %2.3f->%2.3f (%2.3f), "
             "and trans %2.2f->%2.2f (%2.2f)\n",
             (float)DEGREES(max_angle), (float)DEGREES(delta_rot),
             min_scale,max_scale, delta_scale,
             min_trans, max_trans, delta_trans);
      fflush(stdout) ;
    }

    for (x_scale = min_scale ;
         x_scale <= max_scale ;
         x_scale += delta_scale) {
      /*      printf("x_scale = %2.3f\n", x_scale) ;*/
      *MATRIX_RELT(m_scale, 1, 1) = x_scale ;
      for (y_scale = min_scale ;
           y_scale <= max_scale ;
           y_scale += delta_scale) {
        *MATRIX_RELT(m_scale, 2, 2) = y_scale ;
        for (z_scale= min_scale ;
             z_scale <= max_scale;
             z_scale += delta_scale) {
          *MATRIX_RELT(m_scale, 3, 3) = z_scale ;

          /* reset translation values */
          *MATRIX_RELT(m_scale, 1, 4) =
            *MATRIX_RELT(m_scale, 2, 4) =
              *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
          m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
          MatrixMultiply(m_origin, m_tmp, m_scale) ;

          for (x_angle = min_angle ;
               x_angle <= max_angle ;
               x_angle += delta_rot) {
            m_x_rot = MatrixReallocRotation
                      (4, x_angle, X_ROTATION, m_x_rot) ;
            for (y_angle = min_angle ;
                 y_angle <= max_angle ;
                 y_angle += delta_rot) {
              m_y_rot = MatrixReallocRotation
                        (4, y_angle, Y_ROTATION, m_y_rot);
              m_tmp = MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
              for (z_angle= min_angle;
                   z_angle <= max_angle;
                   z_angle += delta_rot) {
                m_z_rot = MatrixReallocRotation
                          (4, z_angle,Z_ROTATION,m_z_rot);
                m_rot = MatrixMultiply
                        (m_z_rot, m_tmp, m_rot) ;
                m_tmp2 = MatrixMultiply
                         (m_rot, m_origin_inv, m_tmp2) ;
                MatrixMultiply(m_origin, m_tmp2, m_rot) ;

                m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
                m_tmp3 = MatrixMultiply(m_tmp2, m_L, m_tmp3) ;

                for (x_trans = min_trans ;
                     x_trans <= max_trans ;
                     x_trans += delta_trans) {
                  *MATRIX_RELT(m_trans, 1, 4) = x_trans ;
                  for (y_trans = min_trans ;
                       y_trans <= max_trans ;
                       y_trans += delta_trans) {
                    *MATRIX_RELT(m_trans, 2, 4) = y_trans ;
                    for (z_trans= min_trans ;
                         z_trans <= max_trans ;
                         z_trans += delta_trans) {
                      *MATRIX_RELT(m_trans, 3, 4) =
                        z_trans ;

                      m_L_tmp = MatrixMultiply
                                (m_trans, m_tmp3, m_L_tmp) ;
                      log_p =
                        local_GCAcomputeLogSampleProbability
                        (gca,gcas,mri,m_L_tmp,nsamples, DEFAULT_CLAMP);
                      if (log_p > max_log_p) {
                        max_log_p = log_p ;
                        x_max_scale = x_scale ;
                        y_max_scale = y_scale ;
                        z_max_scale = z_scale ;
                        x_max_rot = x_angle ;
                        y_max_rot = y_angle ;
                        z_max_rot = z_angle ;
                        x_max_trans = x_trans ;
                        y_max_trans = y_trans ;
                        z_max_trans = z_trans ;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (Gdiag & DIAG_SHOW) {
      printf("\tmax log p = %2.1f @ R=(%2.3f, %2.3f, %2.3f),"
             "S=(%2.3f,%2.3f,%2.3f), T=(%2.1f, %2.1f, %2.1f)\n",
             max_log_p, DEGREES(x_max_rot), DEGREES(y_max_rot),
             DEGREES(z_max_rot),x_max_scale, y_max_scale, z_max_scale,
             x_max_trans, y_max_trans,z_max_trans) ;
    }

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_scale, 1, 4) =
      *MATRIX_RELT(m_scale, 2, 4) = *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
    *MATRIX_RELT(m_scale,1,1) = x_max_scale ;
    *MATRIX_RELT(m_scale,2,2) = y_max_scale ;
    *MATRIX_RELT(m_scale,3,3) = z_max_scale ;
    m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
    MatrixMultiply(m_origin, m_tmp, m_scale) ;

    x_max_scale = y_max_scale = z_max_scale = 1.0 ;

    mean_scale = (max_scale + min_scale) / 2 ;
    delta_scale = (max_scale-min_scale)/4 ;
    min_scale = mean_scale - delta_scale ;
    max_scale = mean_scale + delta_scale ;

    /* update L to reflect new maximum and search around it */
    MatrixReallocRotation(4, x_max_rot, X_ROTATION, m_x_rot) ;
    MatrixReallocRotation(4, y_max_rot, Y_ROTATION, m_y_rot) ;
    MatrixReallocRotation(4, z_max_rot, Z_ROTATION, m_z_rot) ;
    MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
    MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
    m_tmp2 = MatrixMultiply(m_rot, m_origin_inv, m_tmp2) ;
    MatrixMultiply(m_origin, m_tmp2, m_rot) ;

    m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
    m_tmp3 = MatrixMultiply(m_tmp2, m_L, m_tmp3) ;

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_trans, 1, 4) = x_max_trans ;
    *MATRIX_RELT(m_trans, 2, 4) = y_max_trans ;
    *MATRIX_RELT(m_trans, 3, 4) = z_max_trans ;
    MatrixMultiply(m_trans, m_tmp3, m_L_tmp) ;

    MatrixCopy(m_L_tmp, m_L) ;

    x_max_trans = y_max_trans = z_max_trans = 0.0 ;  /* we've translated
                                                                          transform by
                                                                          old maxs */
    mean_trans = (max_trans + min_trans) / 2 ;
    delta_trans = (max_trans-min_trans)/4 ;
    min_trans = mean_trans - delta_trans ;
    max_trans = mean_trans + delta_trans ;

    /* we've rotated transform to old max */
    x_max_rot = y_max_rot = z_max_rot = 0.0 ;

    mean_angle = (max_angle + min_angle) / 2 ;
    delta_rot = (max_angle-min_angle)/4 ;
    min_angle = mean_angle - delta_rot ;
    max_angle = mean_angle + delta_rot ;
  }

  printf("\n") ;

  MatrixFree(&m_x_rot) ;
  MatrixFree(&m_y_rot) ;
  MatrixFree(&m_z_rot) ;
  MatrixFree(&m_rot) ;
  MatrixFree(&m_tmp) ;
  MatrixFree(&m_origin_inv) ;
  MatrixFree(&m_tmp2) ;
  MatrixFree(&m_trans) ;
  MatrixFree(&m_tmp3) ;
  return(max_log_p) ;
}
double
local_GCAcomputeLogSampleProbability(GCA *gca,
                                     GCA_SAMPLE *gcas,
                                     MRI *mri,
                                     MATRIX *m_L,
                                     int nsamples,
                                     double clamp) {
  static TRANSFORM *transform = NULL ;

  if (!transform)
    transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
  ((LTA *)transform->xform)->xforms[0].m_L = m_L ;
  return(GCAcomputeLogSampleProbability(gca, gcas, mri, transform, nsamples, clamp)) ;
}

#if 0
static MATRIX *
update_optimal_transform(MRI *mri, GCA *gca, GCA_SAMPLE *gcas, int nsamples,
                         MATRIX *m_L, int write_iterations,
                         double min_trans, double max_trans,
                         double min_angle, double max_angle,
                         double min_scale, double max_scale,
                         int nsteps) {
  MATRIX   *m_origin ;
  MRI      *mri_gca  ;
  double   gca_means[3], max_log_p, old_max, scale, delta, mean ;
  int      niter ;

  mri_gca = MRIclone(mri, NULL) ;
  GCAmri(gca, mri_gca) ;

  max_log_p =
    local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples, DEFAULT_CLAMP) ;
  printf("initial log_p = %2.1f\n", max_log_p) ;

  MRIcenterOfMass(mri_gca, gca_means, 0) ;
  MRIfree(&mri_gca) ;

  m_origin = MatrixIdentity(4, NULL) ;
  *MATRIX_RELT(m_origin, 1, 4) = gca_means[0]*(float)center ;
  *MATRIX_RELT(m_origin, 2, 4) = gca_means[1]*(float)center ;
  *MATRIX_RELT(m_origin, 3, 4) = gca_means[2]*(float)center ;
  *MATRIX_RELT(m_origin, 4, 4) = 1 ;

  if (write_iterations != 0) {
    char fname[STRLEN] ;
    MRI *mri_aligned ;

    mri_aligned = MRIlinearTransform(mri, NULL, m_L) ;
    sprintf(fname, "%s_after_centering", parms.base_name) ;
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
    sprintf(fname, "%s_after_centering.mgh", parms.base_name) ;
    printf("writing image after centering to %s...\n", fname) ;
#if 0
    MRIwrite(mri_aligned, fname) ;
#else
    Glta->xforms[0].m_L = m_L ;
    GCAtransformAndWriteSamples(gca, mri, gcas, nsamples,
                                fname, transform) ;
#endif
    MRIfree(&mri_aligned) ;

  }
#define MIN_ITER 1
  niter = 0 ;
  scale = 1.0 ;
  do {
    old_max = max_log_p ;
    max_log_p = find_optimal_translation(gca, gcas, mri, nsamples, m_L,
                                         -scale*max_trans, scale*max_trans,
                                         nsteps, 2) ;

    if (!translation_only)
      max_log_p = find_optimal_rotation
                  (gca, gcas, mri, nsamples, m_L,m_origin,
                   -scale*max_angle, scale*max_angle,
                   nsteps, 3) ;

    if (!noscale && !translation_only) {
      max_log_p = find_optimal_scaling
                  (gca, gcas, mri, nsamples, m_L, m_origin,
                   min_scale, max_scale, nsteps, 3) ;
      max_log_p = find_optimal_scaling_and_rotation
                  (gca, gcas, mri, nsamples,
                   m_L, m_origin,
                   -RADIANS(2*scale),
                   RADIANS(2*scale),
                   1-.025*scale, 1+.025*scale,
                   nsteps/3,nsteps/3,3);
    }
    if (write_iterations != 0) {
      char fname[STRLEN] ;
      MRI *mri_aligned ;

      mri_aligned = MRIlinearTransform(mri, NULL, m_L) ;
      sprintf(fname, "%s%03d", parms.base_name, niter+1) ;
      MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
      sprintf(fname, "%s%03d.mgh", parms.base_name, niter+1) ;
      printf("writing image after centering to %s...\n", fname) ;
#if 0
      MRIwrite(mri_aligned, fname) ;
#else
      Glta->xforms[0].m_L = m_L ;
      GCAtransformAndWriteSamples(gca, mri, gcas, nsamples,
                                  fname, transform) ;
#endif
      MRIfree(&mri_aligned) ;
    }

    printf("scale %2.3f: max=%2.1f, old_max =%2.1f (thresh=%2.1f)\n",
           scale,max_log_p, old_max, old_max+fabs(tol*old_max)) ;

    /* search a finer nbhd (if do-while continues) */
    scale *= 0.5 ;
    mean = (max_scale + min_scale)/2 ;
    delta = (max_scale - min_scale)/2 ;
    max_scale = mean + delta*scale ;
    min_scale = mean - delta*scale ;
  } while (niter++ < MIN_ITER || (max_log_p > old_max+fabs(tol*old_max))) ;

  MatrixFree(&m_origin) ;
  return(m_L) ;
}
#else
#if 0
static MATRIX *
update_optimal_transform(MRI *mri, GCA *gca, GCA_SAMPLE *gcas, int nsamples,
                         MATRIX *m_L,
                         double scale, int write_iterations, int nsteps,
                         int nreductions) {
  double   max_log_p, old_max ;
  int      niter ;

  max_log_p =
    local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples, DEFAULT_CLAMP) ;
  printf("initial log_p = %2.1f\n", max_log_p) ;

  if (write_iterations != 0) {
    char fname[STRLEN] ;
    MRI *mri_aligned ;

    mri_aligned = MRIlinearTransform(mri, NULL, m_L) ;
    sprintf(fname, "%s_after_centering", parms.base_name) ;
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
    sprintf(fname, "%s_after_centering.mgh", parms.base_name) ;
    printf("writing image after centering to %s...\n", fname) ;
#if 0
    MRIwrite(mri_aligned, fname) ;
#else
Glta->xforms[0].m_L = m_L ;
GCAtransformAndWriteSamples(gca, mri, gcas, nsamples,
                            fname, transform) ;
#endif
    MRIfree(&mri_aligned) ;
  }

  niter = 0 ;
  do {
    old_max = max_log_p ;
    max_log_p =
      find_optimal_3x4(gca, gcas, mri, nsamples, m_L, scale, nsteps) ;

    if (write_iterations != 0) {
      char fname[STRLEN] ;
      MRI *mri_aligned ;

      mri_aligned = MRIlinearTransform(mri, NULL, m_L) ;
      sprintf(fname, "%s%03d", parms.base_name, niter+1) ;
      MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
      sprintf(fname, "%s%03d.mgh", parms.base_name, niter+1) ;
      printf("writing image after centering to %s...\n", fname) ;
#if 0
      MRIwrite(mri_aligned, fname) ;
#else
Glta->xforms[0].m_L = m_L ;
GCAtransformAndWriteSamples(gca, mri, gcas, nsamples,
                            fname, transform) ;
#endif
      MRIfree(&mri_aligned) ;
    }

    printf("scale %2.3f: max=%2.1f, old_max =%2.1f (thresh=%2.1f)\n",
           scale,max_log_p, old_max, old_max+fabs(tol*old_max)) ;

    /* search a finer nbhd (if do-while continues) */
    scale *= 0.5 ;
  } while (niter++ < nreductions ||
           (max_log_p > old_max+fabs(tol*old_max))) ;

  return(m_L) ;
}

#define MAX_ITER 10
#define MAX_DIM  5
static double
find_optimal_3x4(GCA *gca, GCA_SAMPLE *gcas, MRI *mri, int nsamples,
                 MATRIX *m_L, float scale, int nsteps) {
  MATRIX   *m_L_best, *m_L_tmp ;
  double   log_p, max_log_p, sf, old_log_p,
  min_vals[3*4], max_vals[3*4], best_vals[3*4], deltas[3*4],
  m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34 ;
  int      row, col, j, i, ndim, indices[4*3], k ;

  m_L_tmp = MatrixCopy(m_L, NULL) ;
  m_L_best = MatrixCopy(m_L, NULL) ;

  for (i = 0 ; i < MAX_ITER ; i++) {
    max_log_p =
      local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples, DEFAULT_CLAMP) ;
    old_log_p = max_log_p ;

    ndim = 0 ;
    for (j = 0 ; j < 4*3 ; j++)
      indices[j] = j ;

    /* now permute them */
    for (j = 0 ; j < 4*3 ; j++) {
      int tmp ;
      k = nint(randomNumber(0,4*3-.1)) ;
      tmp = indices[j] ;
      indices[j] = k ;
      indices[k] = tmp ;
    }

    for (k = 0 ; k < 4*3 ; k++) {
      j = indices[k] ;
      col = j%4+1 ;
      row = j/4+1 ;
      if (col == 4)   /* translations have much bigger values */
        sf = 100 ;
      else
        sf = 1 ;
      if (k < MAX_DIM)   /* include this dimension */
      {
        best_vals[j] = *MATRIX_RELT(m_L,row,col) ;
        min_vals[j] = best_vals[j] - sf*scale*fabs(best_vals[j]) ;
        max_vals[j] = best_vals[j] + sf*scale*fabs(best_vals[j]) ;
        deltas[j] = (max_vals[j] - min_vals[j]) / (double)nsteps ;
        if (FZERO(deltas[j]))
          return(max_log_p) ;
        ndim++ ;
      } else   /* leave this dimension out */
      {
        best_vals[j] = *MATRIX_RELT(m_L,row,col) ;
        min_vals[j] = max_vals[j] = best_vals[j] ;
        deltas[j] = 10000 ;
      }
    }

    printf("computing optimal xform for %d dimensions...\n", ndim) ;
    for (m11 = min_vals[0] ;
         m11 <= max_vals[0] ;
         m11 += deltas[0]) {
      *MATRIX_RELT(m_L_tmp, 1, 1) = m11 ;
      for (m12 = min_vals[1] ;
           m12 <= max_vals[1] ;
           m12 += deltas[1]) {
        *MATRIX_RELT(m_L_tmp, 1, 2) = m12 ;
        for (m13 = min_vals[2] ;
             m13 <= max_vals[2] ;
             m13 += deltas[2]) {
          *MATRIX_RELT(m_L_tmp, 1, 3) = m13 ;
          for (m14 = min_vals[3] ;
               m14 <= max_vals[3] ;
               m14 += deltas[3]) {
            *MATRIX_RELT(m_L_tmp, 1, 4) = m14 ;
            for (m21 = min_vals[4] ;
                 m21 <= max_vals[4] ;
                 m21 += deltas[4]) {
              *MATRIX_RELT(m_L_tmp, 2, 1) = m21 ;
              for (m22 = min_vals[5] ;
                   m22 <= max_vals[5] ;
                   m22 += deltas[5]) {
                *MATRIX_RELT(m_L_tmp, 2,2) = m22 ;
                for (m23 = min_vals[6] ;
                     m23 <= max_vals[6] ;
                     m23 += deltas[6]) {
                  *MATRIX_RELT(m_L_tmp, 2,3) = m23 ;
                  for (m24 = min_vals[7] ;
                       m24 <= max_vals[7] ;
                       m24 += deltas[7]) {
                    *MATRIX_RELT(m_L_tmp, 2,4) = m24 ;
                    for (m31 = min_vals[8] ;
                         m31 <= max_vals[8] ;
                         m31 += deltas[8]) {
                      *MATRIX_RELT(m_L_tmp, 3,1) = m31 ;
                      for (m32 = min_vals[9] ;
                           m32 <= max_vals[9] ;
                           m32 += deltas[9]) {
                        *MATRIX_RELT(m_L_tmp, 3,2) =
                          m32 ;
                        for (m33 = min_vals[10] ;
                             m33 <= max_vals[10] ;
                             m33 += deltas[10]) {
                          *MATRIX_RELT(m_L_tmp, 3,3) =
                            m33 ;
                          for (m34 = min_vals[11] ;
                               m34 <= max_vals[11] ;
                               m34 += deltas[11]) {
                            *MATRIX_RELT(m_L_tmp, 3,4) = m34 ;
                            {
                              log_p =
                                local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L_tmp,nsamples, DEFAULT_CLAMP) ;
                              if (log_p > max_log_p) {
                                printf
                                ("new optimal "
                                 "matrix found "
                                 "(%2.3f --> "
                                 "%2.3f)\n",
                                 max_log_p,
                                 log_p) ;
                                max_log_p = log_p ;
                                MatrixCopy
                                (m_L_tmp,
                                 m_L_best) ;
                                MatrixPrint
                                (stdout,
                                 m_L_best) ;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if ((Gdiag & DIAG_SHOW) && (old_log_p < max_log_p)) {
      printf("iter %d, scale %2.3f: "
             "max log p = %2.1f with matrix:\n", i, scale, max_log_p) ;
      MatrixPrint(stdout, m_L_best) ;
    }

    MatrixCopy(m_L_best, m_L) ;
  }

  MatrixFree(&m_L_tmp) ;
  MatrixCopy(m_L_best, m_L) ;
  MatrixFree(&m_L_best) ;
  return(max_log_p) ;
}

#endif
#endif

#if 0
static float
compare_transform(MRI *mri, GCA *gca, GCA_SAMPLE *gcas, int nsamples,
                  MATRIX *m_L) {
  static TRANSFORM *transform = NULL ;
  double  log_p1, log_p2 ;

  if (!transform)
    transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
  ((LTA *)transform->xform)->xforms[0].m_L = m_L ;

  log_p1 = GCAcomputeLogSampleProbability(gca, gcas, mri, transform,nsamples,DEFAULT_CLAMP);
  ((LTA *)transform->xform)->xforms[0].m_L = MatrixIdentity(4, NULL) ;
  log_p2 = compareLogSampleProbability(gca, gcas, mri, transform,nsamples);
  MatrixFree(&((LTA *)transform->xform)->xforms[0].m_L) ;
  return(log_p1 - log_p2) ;
}

float
compareLogSampleProbability(GCA *gca, GCA_SAMPLE *gcas,
                            MRI *mri_inputs, TRANSFORM *transform,
                            int nsamples) {
  int        x, y, z, width, height, depth, val,
  i, xp, yp, zp, max_i ;
  float      dist ;
  double     total_log_p, log_p, max_diff, diff ;

  /* go through each GC in the sample and compute the probability of
     the image at that point.
  */
  width = mri_inputs->width ;
  height = mri_inputs->height;
  depth = mri_inputs->depth ;
  TransformInvert(transform, mri_inputs) ;
  max_diff = -1 ;
  max_i = -1 ;
  for (total_log_p = 0.0, i = 0 ; i < nsamples ; i++) {
    if (i == Gdiag_no)
      DiagBreak() ;
    if (Gdiag_no == gcas[i].label)
      DiagBreak() ;

    xp = gcas[i].xp ;
    yp = gcas[i].yp ;
    zp = gcas[i].zp ;
    GCApriorToSourceVoxel
    (gca,mri_inputs, transform, xp, yp, zp, &x, &y, &z) ;
    gcas[i].x = x ;
    gcas[i].y = y ;
    gcas[i].z = z ;
    val = MRIvox(mri_inputs, x, y, z) ;

    dist = (val-gcas[i].mean) ;
    log_p =
      -log(sqrt(gcas[i].var)) -
      0.5 * (dist*dist/gcas[i].var) +
      log(gcas[i].prior) ;
    total_log_p += log_p ;
    /* current one should be better (bigger) than one in struct */
    diff = gcas[i].log_p - log_p  ;
    if (diff > max_diff) {
      max_diff = diff ;
      max_i = i ;
    }

    /*    gcas[i].log_p = log_p ;*/

    if (!check_finite("2", total_log_p)) {
      fprintf(stderr,
              "total log p not finite at (%d, %d, %d) var=%2.2f\n",
              x, y, z, gcas[i].var) ;
      DiagBreak() ;
    }
  }

  return((float)total_log_p) ;
}

static int
recompute_labels(MRI *mri, GCA *gca, GCA_SAMPLE *gcas,
                 int nsamples, MATRIX *m_L) {
  int            xp,yp,zp,x, y, z, label, n, i, nchanged = 0 ;
  double         val ;
  GCA_PRIOR      *gcap ;
  float          log_p ;
  GC1D           *gc ;

  for (i = 0 ; i < nsamples ; i++) {
    if (i == Gdiag_no)
      DiagBreak() ;
    xp = gcas[i].xp ;
    yp = gcas[i].yp ;
    zp = gcas[i].zp ;
    GCApriorToSourceVoxel(gca,mri, transform, xp, yp, zp, &x, &y, &z) ;
    val = MRIvox(mri, x, y, z) ;
    label = GCAcomputeMAPlabelAtLocation(gca, xp,yp,zp,val,&n,&log_p);
    gcap = &gca->priors[xp][yp][zp] ;
    if (n >= 0) {
      if (label != gcas[i].label)
        nchanged++ ;
      gcas[i].label = label ;
      gcas[i].prior = gcap->priors[n] ;
      gc = GCAfindPriorGC(gca, xp, yp, zp, label) ;
      if (gc) {
        gcas[i].var = sqrt(gc->var) ;
        gcas[i].mean = gc->mean ;
      } else   /* probably out of field of view */
      {
        gcas[i].var = 1 ;
        gcas[i].mean = 0 ;
      }

    }
  }

  return(NO_ERROR) ;
}

#endif

MRI *
MRIremoveFace(MRI *mri_src,
              MRI *mri_dst,
              LTA *lta,
              GCA *gca,
              GCA *gca_face,
              int radius, int fill_val) {
  int       x, y, z, frame, n, nerased = 0, nskipped = 0, i ;
  TRANSFORM *transform ;
  GCA_PRIOR *gcap_face ;
  MRI       *mri_brain, *mri_brain2 ;
  float     wm_means[MAX_GCA_INPUTS], gm_means[MAX_GCA_INPUTS],
  means[MAX_GCA_INPUTS], mah_dist ;
  float     vals[MAX_GCA_INPUTS] ;
  MATRIX    *m_wm_covar, *m_gm_covar, *m ;

  transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
  ((LTA *)transform->xform)->xforms[0].m_L = lta->xforms[0].m_L ;

  mri_dst = MRIcopy(mri_src, mri_dst) ;

  mri_brain = fill_brain_volume(mri_dst, gca, transform, radius) ;
  mri_brain2 = MRIcopy(mri_brain, NULL) ;
  for (n = 0 ; n < radius ; n++)
    MRIdilate(mri_brain2, mri_brain2) ;

  GCAlabelMeanFromImage
  (gca, transform, mri_src, Left_Cerebral_White_Matter, wm_means) ;
  GCAlabelMeanFromImage
  (gca, transform, mri_src, Right_Cerebral_White_Matter, means) ;
  for (i = 0 ; i < gca->ninputs ; i++)
    wm_means[i] = (wm_means[i] + means[i]) / 2 ;
  GCAlabelMeanFromImage
  (gca, transform, mri_src, Left_Cerebral_Cortex, gm_means) ;
  GCAlabelMeanFromImage
  (gca, transform, mri_src, Right_Cerebral_Cortex, means) ;
  for (i = 0 ; i < gca->ninputs ; i++)
    gm_means[i] = (gm_means[i] + means[i]) / 2 ;

  m_wm_covar = GCAlabelCovariance(gca, Left_Cerebral_White_Matter, NULL) ;
  m = GCAlabelCovariance(gca, Right_Cerebral_White_Matter, NULL) ;
  MatrixAdd(m_wm_covar, m, m_wm_covar) ;

  m_gm_covar = GCAlabelCovariance(gca, Left_Cerebral_Cortex, NULL) ;
  m = GCAlabelCovariance(gca, Right_Cerebral_Cortex, NULL) ;
  MatrixAdd(m_gm_covar, m, m_gm_covar) ;

  printf("using wm = %2.1f, gm = %2.1f\n",wm_means[0], gm_means[0]) ;
  printf("wm covar:\n") ;
  MatrixPrint(stdout, m_wm_covar) ;
  printf("gm covar:\n") ;
  MatrixPrint(stdout, m_gm_covar) ;
  for (x = 0 ; x < mri_src->width ; x++) {
    if (((x+1)%10) == 0) {
      printf(".") ;
      fflush(stdout) ;
    }
    for (y = 0 ; y < mri_src->height ; y++) {
      for (z = 0 ; z < mri_src->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        if (MRIvox(mri_brain, x, y, z)  > 0)
          continue ;
        if (MRIvox(mri_brain2, x, y, z) > 0)   /* within 2*radius -
                                                                                                test image value */
        {
          load_vals(mri_src, x, y, z, vals, gca->ninputs) ;
          mah_dist = compute_mah_dist
                     (wm_means, m_wm_covar, vals, gca->ninputs) ;
          if (mah_dist < 2) {
            nskipped++ ;
            continue ;
          }
          mah_dist = compute_mah_dist
                     (gm_means, m_gm_covar, vals, gca->ninputs) ;
          if (mah_dist < 2) {
            nskipped++ ;
            continue ;
          }
        }

        gcap_face = getGCAP(gca_face, mri_src, transform, x, y, z) ;
        if (gcap_face == NULL)
          continue ;
        for (n = 0 ; n < gcap_face->nlabels ; n++) {
          if (!IS_UNKNOWN(gcap_face->labels[n]) &&
              (gcap_face->priors[n] > 0.0001)) {
            for (frame = 0 ; frame < mri_dst->nframes ; frame++)
              MRIseq_vox(mri_dst, x, y, z, frame) = fill_val ;
            nerased++ ;
            break ;
          }
        }
      }
    }
  }

  printf("\n%d face voxels erased, %d ambiguous voxels retained\n",
         nerased, nskipped) ;
  MRIfree(&mri_brain) ;
  return(mri_dst)  ;
}

#if 0
static int
brain_in_region(GCA *gca, MRI *mri, TRANSFORM *transform,
                int x, int y, int z, int whalf) {
  int xk, yk, zk, xi, yi, zi, n ;
  GCA_PRIOR *gcap ;

  for (xk = -whalf ; xk <= whalf ; xk++) {
    xi = mri->xi[x+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++) {
      yi = mri->yi[y+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++) {
        zi = mri->zi[z+zk] ;
        gcap = getGCAP(gca, mri, transform, xi, yi, zi) ;
        if (gcap == NULL)
          continue ;
        for (n = 0 ; n < gcap->nlabels ; n++)
          if (IS_BRAIN(gcap->labels[n]))
            return(1) ;
      }
    }
  }
  return(0) ;
}
#endif

static MRI *
fill_brain_volume(MRI *mri, GCA *gca, TRANSFORM *transform, int radius) {
  int       i, x, y, z, n ;
  MRI       *mri_brain ;
  GCA_PRIOR *gcap ;

  mri_brain = MRIclone(mri, NULL) ;

  for (x = 0 ; x < mri->width ; x++) {
    for (y = 0 ; y < mri->height ; y++) {
      for (z = 0 ; z < mri->depth ; z++) {
        MRIvox(mri_brain, x, y, z) = 0 ;
        gcap = getGCAP(gca, mri, transform, x, y, z) ;
        if (gcap == NULL)
          continue ;
        for (n = 0 ; n < gcap->nlabels ; n++) {
          if (IS_BRAIN(gcap->labels[n])) {
            MRIvox(mri_brain, x, y, z) = 128 ;
            break ;
          }
        }
      }
    }
  }
  for (i = 0 ; i < radius ; i++)
    MRIdilate(mri_brain, mri_brain) ;
  return(mri_brain) ;
}

double
compute_mah_dist(float *means, MATRIX *m_cov, float *vals, int ninputs) {
  double  dsq ;
  VECTOR  *v_vals, *v_tmp ;
  int     i ;
  MATRIX  *m_cov_inv = NULL ;

  if (ninputs == 1) {
    double v ;
    v = vals[0] - means[0] ;
    dsq = v*v / *MATRIX_RELT(m_cov, 1, 1) ;
    return(sqrt(dsq)) ;
  }

  v_vals = VectorAlloc(ninputs, MATRIX_REAL) ;

  for (i = 0 ; i < ninputs ; i++)
    VECTOR_ELT(v_vals, i+1) = means[i] - vals[i] ;

  m_cov_inv = MatrixInverse(m_cov, m_cov_inv) ;
  if (!m_cov_inv)
    ErrorExit(ERROR_BADPARM, "singular covariance matrix!") ;
  MatrixSVDInverse(m_cov, m_cov_inv) ;

  v_tmp = MatrixMultiply(m_cov_inv, v_vals, NULL) ;  /* v_means is now
                                                                inverse(cov) *
                                                                v_vals */
  dsq = VectorDot(v_tmp, v_vals) ;

  MatrixFree(&m_cov_inv) ;
  VectorFree(&v_vals) ;
  return(sqrt(dsq)) ;
}
