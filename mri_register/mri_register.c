/**
 * @file  mri_ca_register.c
 * @brief high-dimensional alignment with canonical atlas
 *
 * Example usage:
 *  mri_ca_register -align -mask brainmask.mgz \
 *    -T transforms/talairach.lta norm.mgz \
 *    $FREESURFER_HOME/average/RB_all_2006-02-15.gca \
 *    transforms/talairach.m3z
 *
 * Inputs:
 *    brainmask.mgz
 *    transforms/talairach.lta
 *    norm.mgz
 *
 * Outputs:
 *    transforms/talairach.m3z
 *
 * Reference:
 *   "Automatically Parcellating the Human Cerebral Cortex", Fischl et al.
 *   (2004). Cerebral Cortex, 14:11-22.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:24 $
 *    $Revision: 1.25 $
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "utils.h"
#include "gca.h"
#include "cma.h"
#include "mrinorm.h"
#include "gcamorph.h"
#include "transform.h"
#include "mrisegment.h"
#include "version.h"

static int nogibbs=0 ;
static int remove_bright =0 ;
static int map_to_flash = 0 ;
static double TRs[MAX_GCA_INPUTS] ;
static double fas[MAX_GCA_INPUTS] ;
static double TEs[MAX_GCA_INPUTS] ;

char         *Progname ;
static GCA_MORPH_PARMS  parms ;

static int avgs = 0 ;  /* for smoothing conditional densities */
static int read_lta = 0 ;
static char *mask_fname = NULL ;
static char *norm_fname = NULL ;
static int renormalize = 0 ;
static int renormalize_new = 0 ;
static int renormalize_align = 0 ;
static int renormalize_align_after = 0 ;
static int demons = 1 ;
static int use_pvals = 0 ;

static char *long_reg_fname = NULL ;
//static int inverted_xform = 0 ;

static int max_rounds = 5 ;
static int max_reductions = 3 ;
static char *write_gca_fname = NULL ;
static float regularize = 0 ;
static float regularize_mean = 0 ;
static char *example_T1 = NULL ;
static char *example_segmentation = NULL ;
static int register_wm_flag = 0 ;

static double TR = -1 ;
static double alpha = -1 ;
static double TE = -1 ;
static char *tl_fname = NULL ;

static int ninsertions = 0 ;
static int insert_labels[MAX_INSERTIONS] ;
static int insert_intensities[MAX_INSERTIONS] ;
static int insert_coords[MAX_INSERTIONS][3] ;
static int insert_whalf[MAX_INSERTIONS] ;

static char *sample_fname = NULL ;
static char *transformed_sample_fname = NULL ;
static char *normalized_transformed_sample_fname = NULL ;
static char *ctl_point_fname = NULL ;
static int novar = 1 ;
static int reinit = 0 ;

static int use_contrast = 0 ;
static float min_prior = MIN_PRIOR ;
static int reset = 0 ;

static FILE *diag_fp = NULL ;

static int translation_only = 0 ;
static int get_option(int argc, char *argv[]) ;
static int write_vector_field(MRI *mri, GCA_MORPH *gcam, char *vf_fname) ;
static int remove_bright_stuff(MRI *mri, GCA *gca, TRANSFORM *transform) ;

static char *renormalization_fname = NULL ;
static char *tissue_parms_fname = NULL ;
static int center = 1 ;
static int nreductions = 1 ;
static char *xform_name = NULL ;
static int noscale = 0 ;
static int transform_loaded = 0 ;
static char *gca_mean_fname = NULL ;
static TRANSFORM  *transform = NULL ;
static char *vf_fname = NULL ;

static double blur_sigma = 0.0f ;

static int handle_expanded_ventricles = 0;

static int do_secondpass_renorm = 0;

/*
   command line consists of three inputs:

   argv[1]  - directory containing 'canonical' brain
   argv[2]  - directory containing brain to be registered
   argv[3]  - directory in which to write out registered brain.
*/

#define NPARMS           12
#define DEFAULT_CTL_POINT_PCT   .25
static double ctl_point_pct = DEFAULT_CTL_POINT_PCT ;

//static double sigmas[] = { 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 17.0, 20.0, 25.0, 30.0} ;
static double sigmas[] = { 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 25.0 } ;
#define NSIGMAS (sizeof(sigmas)/sizeof(sigmas[0]))
static int nsigmas = NSIGMAS ;

GCAM_MULTISCALE *GCAMMSalloc(GCA_MORPH *gcam, GCA *gca, double *sigmas, int nsigmas, char *write_ms_name) ;
GCAM_MULTISCALE *GCAMMSread(char *base_fname, GCA_MORPH *gcam, int nsigmas);
int GCAMMSwrite(GCAM_MULTISCALE *gcam_ms, char *base_fname);
 
GCAM_MULTISCALE *
GCAMMSalloc(GCA_MORPH *gcam, GCA *gca, double *sigmas, int nsigmas, char *write_ms_name)
{
  GCAM_MULTISCALE *gcam_ms ;
  int             s ;
  GCA             *gca_tmp ;
  float           base_sigma = 1.0 ;

  gcam_ms = (GCAM_MS *)calloc(1, sizeof(GCAM_MULTISCALE)) ;
  gcam_ms->nsigmas = nsigmas+1 ;  // including original atlas (unsmoothed)
  
  gcam_ms->gcam = gcam ;
  gcam_ms->sigmas = (double *)calloc(nsigmas+1, sizeof(sigmas[0])) ;
  gcam_ms->gcas = (GCA **)calloc(nsigmas+1, sizeof(GCA*)) ;
  gcam_ms->gcas[0] = gca ;
  gcam_ms->sigmas[0] = 0 ;
  for (s = 0 ; s < nsigmas ; s++)
  {
    if (s == Gdiag_no)
      DiagBreak() ;
    printf("generating atlas at smoothness level s=%2.3f\n", sigmas[s]) ;
    gcam_ms->sigmas[s+1] = sigmas[s] ;
    gca_tmp = GCAsmooth(gcam_ms->gcas[s], gcam_ms->sigmas[s+1]-gcam_ms->sigmas[s]) ;
    if (0 && sigmas[s]/base_sigma > 4)   // change node spacing
    {
      printf("setting node spacing to be %2.0f mm\n", gca_tmp->node_spacing*2) ;
      base_sigma = sigmas[s] ;
      gcam_ms->gcas[s+1] = GCAnodeDownsample2(gca_tmp) ;
      GCAfree(&gca_tmp) ;
    }
    else
      gcam_ms->gcas[s+1] = gca_tmp ;
    if (write_ms_name)
    {
      gcam_ms->nsigmas = s+2 ;
      GCAMMSwrite(gcam_ms, write_ms_name) ;
    }
  }
  printf("multi-scale atlas generation complete\n") ;
  gcam->vgcam_ms = gcam_ms ;
  return(gcam_ms) ;
}

int
GCAMMSwrite(GCAM_MULTISCALE *gcam_ms, char *base_fname)
{
  int  s ;
  char fname[STRLEN] ;
  FILE *fp ;

  sprintf(fname, "%s_ms.dat", base_fname) ;
  fp = fopen(fname, "w") ;
  if (fp == NULL)
    ErrorExit(ERROR_NOFILE, "GCAMMSwrite(%s): could not open %s",
              base_fname, fname) ;
  fprintf(fp, "nlevels=%d\n", gcam_ms->nsigmas) ;
  for (s = 0 ; s < gcam_ms->nsigmas ; s++)
  {
    sprintf(fname, "%s_sigma%.2f.gca", base_fname, gcam_ms->sigmas[s]) ;
    GCAwrite(gcam_ms->gcas[s], fname) ;
    fprintf(fp, "sigma=%f\n",gcam_ms->sigmas[s]) ;
  }
  fclose(fp) ;
  return(NO_ERROR) ;
}

GCAM_MULTISCALE *
GCAMMSread(char *base_fname, GCA_MORPH *gcam, int nsigmas)
{
  int             s, n ;
  char            fname[STRLEN] ;
  FILE            *fp ;
  double          sigmas[1000];
  GCAM_MULTISCALE *gcam_ms;

  sprintf(fname, "%s_ms.dat", base_fname) ;
  fp = fopen(fname, "r") ;
  if (fp == NULL)
    ErrorExit(ERROR_NOFILE, "GCAMMSread(%s): could not open %s",
              base_fname, fname) ;
  fscanf(fp, "nlevels=%d\n", &n) ;
  if (nsigmas <= 0)
    nsigmas = n ;
  for (s = 0 ; s < nsigmas ; s++)
    fscanf(fp, "sigma=%lf\n", &sigmas[s]) ;
  fclose(fp) ;

  gcam_ms = (GCAM_MS *)calloc(1, sizeof(GCAM_MULTISCALE)) ;
  gcam_ms->nsigmas = nsigmas ;  // including original atlas (unsmoothed)
  
  gcam_ms->gcam = gcam ;
  gcam_ms->sigmas = (double *)calloc(nsigmas, sizeof(sigmas[0])) ;
  gcam_ms->gcas = (GCA **)calloc(nsigmas, sizeof(GCA*)) ;
  for (s = 0 ; s < nsigmas ; s++)
  {
    sprintf(fname, "%s_sigma%.2f.gca", base_fname, sigmas[s]) ;
    printf("reading atlas at smoothness level s=%2.2f\n", sigmas[s]) ;
    gcam_ms->sigmas[s] = sigmas[s] ;
    gcam_ms->gcas[s] = GCAread(fname) ;
  }
  printf("multi-scale atlas reading complete\n") ;
  gcam_ms->nsigmas = nsigmas ;
  gcam->vgcam_ms = gcam_ms ;
  return(gcam_ms) ;
}


static char *read_ms_name = NULL ;
static char *write_ms_name = NULL ;
static char *read_dist_name = NULL ;
static char *write_dist_name = NULL ;

int
main(int argc, char *argv[]) {
  char         *gca_fname, *in_fname, *out_fname, fname[STRLEN], **av, *dot ;
  MRI          *mri_inputs, *mri_tmp, *mri_atlas_dist_map, *mri_atlas_labels ;
  GCA          *gca /*, *gca_tmp, *gca_reduced*/ ;
  int          ac, nargs, ninputs, input, extra = 0 ;
  int          msec, hours, minutes, seconds, round, done, nvox, reductions ;
  struct timeb start ;
  GCA_MORPH    *gcam ;
  GCAM_MS      *gcam_ms ;
  double       log_posterior, last_log_posterior, mean_scale, last_mean_scale,
               gca_log_posterior, last_gca_log_posterior, min_scale ;

	// for GCA Renormalization with Alignment (if called sequentially)
  float        label_scales[MAX_CMA_LABELS], label_offsets[MAX_CMA_LABELS];
	float        label_peaks[MAX_CMA_LABELS];
  int          label_computed[MAX_CMA_LABELS];
	int          got_scales =0;

  parms.l_log_likelihood = 0.0 ;
  parms.l_multiscale = 1.0 ;
  parms.niterations = 500 ;
  parms.levels = 6 ;
  parms.scale_smoothness = 1 ;
  parms.uncompress = 0 ;
  parms.npasses = 1 ;
  parms.diag_write_snapshots = 1 ;
  parms.diag_sample_type = SAMPLE_TRILINEAR ;
  parms.relabel_avgs = -1 ;  /* never relabel, was 1 */
  parms.reset_avgs = 0 ;  /* reset metric properties when navgs=0 */
  parms.dt = 0.05 ;  /* was 5e-6 */
  parms.momentum = 0.9 ;
  parms.tol = 1 ;  /* at least 1% decrease in sse */
  parms.l_jacobian = 1.0 ;
  parms.l_label = 1.0 ;
  parms.l_map = 0.0 ;
  parms.label_dist = 10.0 ;
  parms.l_smoothness = 1 ;
  parms.start_t = 0 ;
  parms.max_grad = 0.3 ;
  parms.sigma = 2.0f ;
  parms.exp_k = 20 ;
  parms.min_avgs = 0 ;
  parms.navgs = 256 ;
  parms.noneg = True ;
  parms.log_fp = NULL ;
  parms.ratio_thresh = 0.1 ;
  parms.nsmall = 1 ;
  parms.integration_type = GCAM_INTEGRATE_BOTH ;

  parms.l_multiscale = 0.0 ;
  parms.l_dtrans = 1.0 ;   parms.l_label = 0.0 ;
  Progname = argv[0] ;
  setRandomSeed(-1L) ;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  nargs = handle_version_option 
    (argc, argv, 
     "$Id: mri_register.c,v 1.25 2011/03/02 00:04:24 nicks Exp $", 
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4) {
    printf("mri_register [options] inbrain template outputfilename\n");
    printf("\n");
    printf("  -align\n");
    printf("  -T lta \n");
    printf("  -mask brainmask\n");
    printf("  -invert-and-save gcamfile\n");
    printf("\n");
    printf("  -dist distance\n");
    printf("  -regularize regularize\n");
    printf("  -regularize-mean regularizemean\n");
    printf("  -scale_smoothness smoothness\n");
    printf("  -nobright \n");
    printf("  -renormalize\n");
    printf("  -read_lta \n");
    printf("  -smoothness smoothness\n");
    printf("  -samples samples\n");
    printf("  -nsmall nsmall\n");
    printf("  -fixed\n");
    printf("  -optimal\n");
    printf("  -noneg noneg\n");
    printf("  -wm\n");
    printf("  -min_avgs <n>\n");
    printf("  -tl\n");
    printf("  -relabel relable\n");
    printf("  -relabel_avgs relable_avgs\n");
    printf("  -reset_avgs reset_avgs\n");
    printf("  -vf vf_fname\n");
    printf("  -diag diagfname\n");
    printf("  -tr tr\n");
    printf("  -te te\n");
    printf("  -example T1 seg\n");
    printf("  -<no>bigventricles\n");
    printf("  -uncompress\n");
    printf("  -secondpassrenorm\n");
    printf("\n");
    printf("\n");
    exit(1);
  }

  ninputs = argc-3 ;
  printf("reading %d input volumes...\n", ninputs) ;
  in_fname = argv[1] ;
  gca_fname = argv[ninputs+1] ;
  out_fname = argv[ninputs+2] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
  //  Gdiag |= DIAG_WRITE ;
  printf("logging results to %s.log\n", parms.base_name) ;

  TimerStart(&start) ;

  // build frames from ninputs ////////////////////////////////
  for (input = 0 ; input < ninputs ; input++) {
    in_fname = argv[1+input] ;
    printf("reading input volume '%s'...\n", in_fname) ;
    fflush(stdout) ;
    mri_tmp = MRIread(in_fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
                Progname, in_fname) ;

    TRs[input] = mri_tmp->tr ;
    fas[input] = mri_tmp->flip_angle ;
    TEs[input] = mri_tmp->te ;

    // -mask option
    if (mask_fname) {
      MRI *mri_mask ;

      mri_mask = MRIread(mask_fname) ;
      if (!mri_mask)
        ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                  Progname, mask_fname) ;
      // if mask == 0, then set dst as 0
      MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
      MRIfree(&mri_mask) ;
    }
    if (alpha > 0)
      mri_tmp->flip_angle = alpha ;
    if (TR > 0)
      mri_tmp->tr = TR ;
    if (TE > 0)
      mri_tmp->te = TE ;
    if (input == 0) {
      mri_inputs = MRIallocSequence(mri_tmp->width,
                                    mri_tmp->height,
                                    mri_tmp->depth,
                                    mri_tmp->type,
                                    ninputs+extra) ;
      // first one's header is copied
      MRIcopyHeader(mri_tmp, mri_inputs) ;
    }
    MRIcopyFrame(mri_tmp, mri_inputs, 0, input) ;
    MRIfree(&mri_tmp) ;
  }
  //
  nvox = mri_inputs->width*mri_inputs->height*mri_inputs->depth ;
  printf("reading GCA '%s'...\n", gca_fname) ;
  fflush(stdout) ;
  gca = GCAread(gca_fname) ;
  if (gca == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open GCA %s.\n",
              Progname, gca_fname) ;
  /////////////////////////////////////////////////////////////////
  // Remapping GCA
  /////////////////////////////////////////////////////////////////
  // GCA from (T1, PD) needs to map to the current input
  if (map_to_flash || gca->type == GCA_PARAM) {
    GCA *gca_tmp ;

    printf("mapping GCA into %d-dimensional FLASH space...\n",
           mri_inputs->nframes) ;
    gca_tmp = GCAcreateFlashGCAfromParameterGCA
              (gca, TRs, fas, TEs, 
               mri_inputs->nframes, GCA_DEFAULT_NOISE_PARAMETER) ;
    GCAfree(&gca) ;
    gca = gca_tmp ;
    if (ninputs != gca->ninputs)
      ErrorExit(ERROR_BADPARM,
                "%s: must specify %d inputs, not %d for this atlas\n",
                Progname, gca->ninputs, ninputs) ;
    GCAhistoScaleImageIntensities(gca, mri_inputs, 1) ;
    if (novar)
      GCAunifyVariance(gca) ;
  }
  // GCA from flash needs to map to the current input
  else if (gca->type == GCA_FLASH) {
    GCA *gca_tmp ;

    int need_map_flag = 0;
    int n;

    if (gca->ninputs != ninputs) need_map_flag = 1;
    else {
      for (n = 0 ; n < mri_inputs->nframes; n++) {
        if (!FZERO(gca->TRs[n] - TRs[n])) need_map_flag = 1;
        if (!FZERO(gca->FAs[n] - fas[n])) need_map_flag = 1;
        if (!FZERO(gca->TEs[n] - TEs[n])) need_map_flag = 1;
      }
    }

    if (need_map_flag) {
      printf("mapping %d-dimensional flash atlas into %d-dimensional "
             "input space\n", gca->ninputs, ninputs) ;

      gca_tmp = GCAcreateFlashGCAfromFlashGCA
                (gca, TRs, fas, TEs, mri_inputs->nframes) ;
      GCAfree(&gca) ;
      gca = gca_tmp ;
    }

    GCAhistoScaleImageIntensities(gca, mri_inputs, 1) ;// added by tosa
  }

  if (gca->flags & GCA_XGRAD)
    extra += ninputs ;
  if (gca->flags & GCA_YGRAD)
    extra += ninputs ;
  if (gca->flags & GCA_ZGRAD)
    extra += ninputs ;

  if ((ninputs+extra) != gca->ninputs)
    ErrorExit(ERROR_BADPARM,
              "%s: must specify %d inputs, not %d for this atlas\n",
              Progname, gca->ninputs, ninputs) ;

  /////////////////////////////////////////////////////////////////
  // clear six neighborhood information
  if (nogibbs)
  {
    printf("freeing gibbs priors...") ;
    GCAfreeGibbs(gca) ;
    printf("done.\n") ;
  }

  ////////////////////////////////////////////////
  // -example T1 T1seg option
  if (example_T1) {
    MRI *mri_T1, *mri_seg ;

    mri_seg = MRIread(example_segmentation) ;
    if (!mri_seg)
      ErrorExit(ERROR_NOFILE,"%s: could not read example segmentation from %s",
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

  /////////////////////////////////////////////////
  // -flash_parms fname option
  if (tissue_parms_fname)   /* use FLASH forward model */
    GCArenormalizeToFlash(gca, tissue_parms_fname, mri_inputs) ;

  /////////////////////////////////////////////////
  // -T transform option
  // transform is loaded at get_opt() with -T using TransformRead()
  // assumed to be vox-to-vox
  if (!transform_loaded)   /* wasn't preloaded */
    transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
  else
    // calculate inverse and cache it
    TransformInvert(transform, mri_inputs) ;

  /////////////////////////////////////////////////
  // -novar option  (default novar = 1)
  if (novar)
    GCAunifyVariance(gca) ;

  /////////////////////////////////////////////////
  // XGRAD or YGRAD or ZGRAD set
  // store (x,y,z)gradient info into mri_inputs
  if (gca->flags & GCA_GRAD) {
    int i, start = ninputs ;
    MRI *mri_kernel, *mri_smooth, *mri_grad, *mri_tmp ;

    mri_kernel = MRIgaussian1d(1.0, 30) ;
    mri_smooth = MRIconvolveGaussian(mri_inputs, NULL, mri_kernel) ;

    if (mri_inputs->type != MRI_FLOAT) {
      // change data to float
      mri_tmp = MRISeqchangeType(mri_inputs, MRI_FLOAT, 0, 0, 1) ;
      MRIfree(&mri_inputs) ;
      mri_inputs = mri_tmp ;
    }
    start = ninputs ;
    if (gca->flags & GCA_XGRAD) {
      for (i = 0 ; i < ninputs ; i++) {
        mri_grad = MRIxSobel(mri_smooth, NULL, i) ;
        MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
        MRIfree(&mri_grad) ;
      }
      start += ninputs ;
    }
    if (gca->flags & GCA_YGRAD) {
      for (i = 0 ; i < ninputs ; i++) {
        mri_grad = MRIySobel(mri_smooth, NULL, i) ;
        MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
        MRIfree(&mri_grad) ;
      }
      start += ninputs ;
    }
    if (gca->flags & GCA_ZGRAD) {
      for (i = 0 ; i < ninputs ; i++) {
        mri_grad = MRIzSobel(mri_smooth, NULL, i) ;
        MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
        MRIfree(&mri_grad) ;
      }
      start += ninputs ;
    }

    MRIfree(&mri_kernel) ;
    MRIfree(&mri_smooth) ;
  }

  ///////////////////////////////////////////////////////////
  // -nobright option
  if (remove_bright)
    remove_bright_stuff(mri_inputs, gca, transform) ;

  ///////////////////////////////////////////////////////////
  // -B blur option (default = 0)
  if (!FZERO(blur_sigma)) {
    MRI *mri_tmp, *mri_kernel ;

    mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
    mri_tmp = MRIconvolveGaussian(mri_inputs, NULL, mri_kernel) ;
    MRIfree(&mri_inputs) ;
    mri_inputs = mri_tmp ;
  }


  //////////////////////////////////////////////////////////
  // -debug_voxel Gvx Gvy Gvz option
  if (transform && Gvx > 0) {
    float xf, yf, zf ;

    TransformSample(transform, Gvx, Gvy, Gvz, &xf, &yf, &zf) ;
    Gsx = nint(xf) ;
    Gsy = nint(yf) ;
    Gsz = nint(zf) ;
    printf("mapping by transform (%d, %d, %d) --> "
           "(%d, %d, %d) for rgb writing\n",
           Gvx, Gvy, Gvz, Gsx, Gsy, Gsz) ;
  }

  //////////////////////////////////////////////////////////
  // -regularize val option (default = 0)
  if (regularize > 0)
    GCAregularizeCovariance(gca, regularize) ;

  //////////////////////////////////////////////////////////
  // -X prev.m3d option
  if (xform_name) {
    gcam = GCAMread(xform_name) ;
    if (!gcam)
      ErrorExit(ERROR_NOFILE,
                "%s: could not read transform from %s", Progname, xform_name) ;
    if (long_reg_fname) {
      TRANSFORM *transform_long ;

      transform_long = TransformRead(long_reg_fname) ;
      if (transform_long == NULL)
        ErrorExit(ERROR_NOFILE,
                  "%s: could not read longitudinal registration file %s",
                  Progname, long_reg_fname) ;

      //      if (inverted_xform)
      // {
      //   TransformInvert(transform_long, mri_inputs) ;
      //   TransformSwapInverse(transform_long) ;
      // }
      GCAMapplyTransform(gcam, transform_long) ;
      TransformFree(&transform_long) ;
      //      GCAMwrite(gcam, "combined_gcam.m3z");
    }
  } else // default is to create one
    gcam = GCAMalloc(gca->prior_width, gca->prior_height, gca->prior_depth) ;

  
  //////////////////////////////////////////////////////////////////
  // GCM initialization
  if (!xform_name)  /* only if the transform wasn't previously created */
    GCAMinit(gcam, mri_inputs, gca, transform,
             parms.relabel_avgs >= parms.navgs) ;
  else {
    // added by xhan
    int x, y, z, n, label, max_n, max_label;
    float max_p;
    GC1D *gc;
    GCA_MORPH_NODE  *gcamn ;
    GCA_PRIOR *gcap;

    gcam->ninputs = mri_inputs->nframes ;
    getVolGeom(mri_inputs, &gcam->image);
    GCAsetVolGeom(gca, &gcam->atlas);
    gcam->gca = gca ;
    gcam->spacing = gca->prior_spacing;

    // use gca information
    for (x = 0 ; x < gcam->width ; x++) {
      for (y = 0 ; y < gcam->height ; y++) {
        for (z = 0 ; z < gcam->depth ; z++) {
          gcamn = &gcam->nodes[x][y][z] ;
          gcap = &gca->priors[x][y][z] ;
          max_p = 0 ;
          max_n = -1 ;
          max_label = 0 ;

          // find the label which has the max p
          for (n = 0 ; n < gcap->nlabels ; n++) {
            label = gcap->labels[n] ;   // get prior label
            if (label == Gdiag_no)
              DiagBreak() ;
            if (label >= MAX_CMA_LABEL) {
              printf("invalid label %d at (%d, %d, %d) in prior volume\n",
                     label, x, y, z);
            }
            if (gcap->priors[n] >= max_p) // update the max_p and max_label
            {
              max_n = n ;
              max_p = gcap->priors[n] ;
              max_label = gcap->labels[n] ;
            }
          }

          gcamn->label = max_label ;
          gcamn->n = max_n ;
          gcamn->prior = max_p ;
          gc = GCAfindPriorGC(gca, x, y, z, max_label) ;
          // gc can be NULL
          gcamn->gc = gc ;
          gcamn->log_p = 0 ;

        }
      }
    }

    GCAMcomputeOriginalProperties(gcam) ;
    if (parms.relabel_avgs >= parms.navgs)
      GCAMcomputeLabels(mri_inputs, gcam) ;
    else
      GCAMcomputeMaxPriorLabels(gcam) ;
  }
  GCAMmarkNegativeNodesInvalid(gcam) ;
  
  //note that transform is meaningless when -L option is used! A bug! -xh
  //  if (renormalize)
  //  GCAmapRenormalize(gcam->gca, mri_inputs, transform) ;
  if (renormalize) {
    if (!xform_name)
      GCAmapRenormalize(gcam->gca, mri_inputs, transform) ;
    else {
      TRANSFORM *trans ;
      trans = (TRANSFORM *)calloc(1, sizeof(TRANSFORM)) ;
      trans->type = TransformFileNameType(xform_name);
      trans->xform = (void *)gcam;
      GCAmapRenormalize(gcam->gca, mri_inputs, trans) ;
      free(trans);
    }
  } else if (renormalize_new) {
    if (!xform_name)
      GCAmapRenormalizeByClass(gcam->gca, mri_inputs, transform) ;
    else {
      TRANSFORM *trans ;
      trans = (TRANSFORM *)calloc(1, sizeof(TRANSFORM)) ;
      trans->type = TransformFileNameType(xform_name);
      trans->xform = (void *)gcam;
      GCAmapRenormalizeByClass(gcam->gca, mri_inputs, trans) ;
      free(trans);
    }
  } else if (renormalize_align) {
    LTA _lta, *lta = &_lta ;

    lta->num_xforms = 0 ;
#if 0
    sprintf(fname, "%s.gca", parms.base_name) ;
    gca = GCAread(fname) ;
    sprintf(fname, "%s.lta", parms.base_name) ;
    lta = LTAread(fname) ;

    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN] ;
      sprintf(fname, "%s.log", parms.base_name) ;
      parms.log_fp = fopen(fname, "w") ;
    }
    {
      MRI *mri_seg, *mri_aligned ;
      int l ;
      mri_seg = MRIclone(mri_inputs, NULL) ;
      l = lta->xforms[0].label ;
      GCAbuildMostLikelyVolumeForStructure(gca, mri_seg, l, 0, transform,NULL);
      mri_aligned = MRIlinearTransform(mri_seg, NULL, lta->xforms[0].m_L) ;
      MRIwrite(mri_seg, "s.mgh")  ;
      MRIwrite(mri_aligned, "a.mgh") ;
      MRIfree(&mri_seg) ;
      MRIfree(&mri_aligned) ;
    }
#else
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN] ;
      sprintf(fname, "%s.log", parms.base_name) ;
      parms.log_fp = fopen(fname, "w") ;
    }
    if (read_lta) {
      sprintf(fname, "%s_array.lta", parms.base_name) ;
      lta = LTAread(fname) ;
    } else
      lta = NULL ;
    if (!xform_name) // normal (cross-sectional) processing
    {
      //      GCAmapRenormalize(gcam->gca, mri_inputs, transform) ;
      //   if (read_lta == 0)
      {
        int old_diag ;
        MRI *mri_morphed ;

        old_diag = Gdiag ;
        if (parms.write_iterations == 0)
          Gdiag &= ~DIAG_WRITE ;

#if 0
        {
          parms.tol *= 10 ;
          parms.l_smoothness *= 10 ;
          printf("---------------- doing initial registration ----------------------\n") ;
          GCAMregister(gcam, mri_inputs, &parms) ;
          parms.tol /= 10 ;
          parms.l_smoothness /= 10 ;
          mri_morphed = GCAMmorphToAtlas(mri_inputs, gcam, NULL, -1) ;
          printf("---------------- initial registration complete ----------------------\n") ;
        }
#else
        mri_morphed = mri_inputs ;
#endif

        // GCA Renormalization with Alignment: 
        if (!do_secondpass_renorm) // just run it once
				{
          // initial call (returning the label_* infos)
					// passing lta
          GCAcomputeRenormalizationWithAlignment
          (gcam->gca,
					 mri_morphed,
					 transform,
           parms.log_fp,
					 parms.base_name,
					 &lta,
					 handle_expanded_ventricles,
           label_scales,label_offsets,label_peaks,label_computed) ;
					 got_scales = 1;
				}
        else // run it twice
        {
          // initial call (returning the label_* infos)
					// not passing lta 
          GCAcomputeRenormalizationWithAlignment
          (gcam->gca, mri_morphed, transform,
           parms.log_fp, parms.base_name, NULL, handle_expanded_ventricles,
           label_scales,label_offsets,label_peaks,label_computed) ;

          // sequential call gets passed the results from first call
          // will overwrite the intensity.txt file with combinded results 
          printf("2nd pass renormalization with updated "
                 "intensity distributions\n");
          GCAseqRenormalizeWithAlignment
          (gcam->gca, mri_morphed, transform,
           parms.log_fp, parms.base_name, &lta, handle_expanded_ventricles,
           label_scales,label_offsets,label_peaks,label_computed) ;
					 got_scales = 1;
        }

        Gdiag = old_diag ;
        if (write_gca_fname) {
          printf("writing normalized gca to %s...\n", write_gca_fname) ;
          GCAwrite(gcam->gca, write_gca_fname) ;
        }
      }
    }
    else  // for longitudinal processing
    {
      TRANSFORM *trans ;
      trans = (TRANSFORM *)calloc(1, sizeof(TRANSFORM)) ;
      trans->type = TransformFileNameType(xform_name);
      trans->xform = (void *)gcam;

      /*The following inversion is necessary;
      but do I need to release the memory for the
      inverse transform after the mapRenormalize is done -xhan? */
      TransformInvert(trans, mri_inputs);

        // GCA Renormalization with Alignment: 
        if (!do_secondpass_renorm) // just run it once
        {
          // initial call (returning the label_* infos)
					// passing lta
          GCAcomputeRenormalizationWithAlignment
          (gcam->gca,
					 mri_inputs,
					 trans,
           parms.log_fp,
					 parms.base_name,
					 &lta,
					 handle_expanded_ventricles,
           label_scales,label_offsets,label_peaks,label_computed) ;
					 got_scales = 1;
				}
        else // run it twice (ensure correct output of label intensities in sequential run)
        {
          // initial call (returning the label_* infos)
					// not passing lta 
          GCAcomputeRenormalizationWithAlignment
          (gcam->gca, mri_inputs, trans,
           parms.log_fp, parms.base_name, NULL, handle_expanded_ventricles,
           label_scales,label_offsets,label_peaks,label_computed) ;

          // sequential call gets passed the results from first call
          // will overwrite the intensity.txt file with combinded results 
          printf("2nd pass renormalization with updated "
                 "intensity distributions\n");
          GCAseqRenormalizeWithAlignment
          (gcam->gca, mri_inputs, trans,
           parms.log_fp, parms.base_name, &lta, handle_expanded_ventricles,
           label_scales,label_offsets,label_peaks,label_computed) ;
					 got_scales = 1;
        }

      free(trans);
    }
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      sprintf(fname, "%s.gca", parms.base_name) ;
      printf("writing gca to %s...\n", fname) ;
      GCAwrite(gca, fname) ;
    }
    if (lta && !read_lta) {
      sprintf(fname, "%s_array.lta", parms.base_name) ;

      // should put volume geometry into file,
      // and probably change to RAS->RAS xform
      LTAwrite(lta, fname) ;
    }
    if (DIAG_VERBOSE_ON) {
      MRI *mri_seg, *mri_aligned ;
      int l ;
      lta = LTAread("gcam.lta") ;
      mri_seg = MRIclone(mri_inputs, NULL) ;
      l = lta->xforms[0].label ;
      GCAbuildMostLikelyVolumeForStructure
      (gca, mri_seg, l, 0, transform, NULL) ;
      LTAinvert(lta) ;
      mri_aligned = MRIlinearTransform(mri_seg, NULL, lta->xforms[0].m_L) ;
      MRIwrite(mri_seg, "s.mgh")  ;
      MRIwrite(mri_aligned, "a.mgh") ;
      MRIfree(&mri_seg) ;
      MRIfree(&mri_aligned) ;
    }
#endif
    if (reinit && (xform_name == NULL) && (lta != NULL))
      GCAMreinitWithLTA(gcam, lta, mri_inputs, &parms) ;
    if (DIAG_VERBOSE_ON) {
      MRI *mri_seg ;
      int l ;
      l = lta->xforms[0].label ;
      mri_seg = MRIclone(mri_inputs, NULL) ;
      GCAbuildMostLikelyVolumeForStructure
      (gca, mri_seg, l, 0, transform,NULL) ;
      MRIwrite(mri_seg, "sa.mgh") ;
      MRIfree(&mri_seg) ;
    }
    if (lta) LTAfree(&lta) ;
  }

  if (regularize_mean > 0)
    GCAregularizeConditionalDensities(gca, regularize_mean) ;

  if (parms.write_iterations != 0) {
    char fname[STRLEN] ;
    MRI  *mri_gca /*, *mri_tmp*/ ;
    mri_gca = MRIclone(mri_inputs, NULL) ;
    GCAMbuildMostLikelyVolume(gcam, mri_gca) ;
#if 0
    if (mri_gca->nframes > 1) {
      printf("careg: extracting %dth frame\n", mri_gca->nframes-1) ;
      mri_tmp = MRIcopyFrame(mri_gca, NULL, mri_gca->nframes-1, 0) ;
      MRIfree(&mri_gca) ;
      mri_gca = mri_tmp ;
    }
#endif
    sprintf(fname, "%s_target", parms.base_name) ;
    MRIwriteImageViews(mri_gca, fname, IMAGE_SIZE) ;
    sprintf(fname, "%s_target.mgz", parms.base_name) ;
    printf("writing target volume to %s...\n", fname) ;
    MRIwrite(mri_gca, fname) ;
    MRIfree(&mri_gca) ;
  }

  ///////////////////////////////////////////////////////////////////
  // -reset option
  if (reset) {
    GCAMcopyNodePositions(gcam, CURRENT_POSITIONS, ORIGINAL_POSITIONS) ;
    GCAMstoreMetricProperties(gcam) ;
  }
  mri_atlas_labels = MRIcloneDifferentType(mri_inputs, MRI_FLOAT) ;
  GCAMbuildLabelVolume(gcam, mri_atlas_labels) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRI *mri_morphed ;
    char fname[STRLEN] ;
    
    mri_morphed = GCAMmorphToAtlas(mri_inputs, gcam, NULL, 0, SAMPLE_TRILINEAR) ;
    sprintf(fname, "src_to_atlas%d.mgh", 0) ;
    MRIwrite(mri_morphed, fname) ; MRIfree(&mri_morphed) ;
    
    mri_morphed = GCAMmorphFromAtlas(mri_atlas_labels, gcam, NULL, SAMPLE_NEAREST);
    sprintf(fname, "atlas_to_src%d.mgh", 0) ;
    MRIwrite(mri_morphed, fname) ; MRIfree(&mri_morphed) ;
  }

  if (renormalize_align_after) // 1st morph should be smooth
  {
    int old_diag ;
    TRANSFORM  _transform, *transform = &_transform ;
    GCA_MORPH_PARMS  saved_parms ;
    int       start_t ;

    *(&saved_parms) = *(&parms) ;

    printf("performing initial intensity-based alignment for atlas renormalization\n");
    parms.l_dtrans = 0 ;
    parms.tol *= 5 ;
    parms.l_smoothness *= 5 ;
    parms.levels = 4 ;
    parms.l_log_likelihood = 0.2f ;
    parms.l_label = 1.0 ;
    GCAMregister(gcam, mri_inputs, &parms) ;
    transform->type = MORPH_3D_TYPE ;
    transform->xform = (void *)gcam ;
    old_diag = Gdiag ;
    if (parms.write_iterations == 0)
      Gdiag &= ~DIAG_WRITE ;
    TransformInvert(transform, mri_inputs) ;
    
    Gdiag = old_diag ;
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN] ;
      sprintf(fname, "%s.log", parms.base_name) ;
      parms.log_fp = fopen(fname, "a") ;
    }
		
    // GCA Renormalization with Alignment: 
		// check wether or not this is a sequential call
		if (!got_scales)
		// this is the first (and also last) call
      GCAcomputeRenormalizationWithAlignment
      (gcam->gca,
       mri_inputs,
       transform,
       parms.log_fp,
       parms.base_name,
       NULL,
       handle_expanded_ventricles,
			 label_scales,label_offsets,label_peaks,label_computed) ;
		else // this is a sequential call, pass scales..
      GCAseqRenormalizeWithAlignment
      (gcam->gca,
       mri_inputs,
       transform,
       parms.log_fp,
       parms.base_name,
       NULL,
       handle_expanded_ventricles,
			 label_scales,label_offsets,label_peaks,label_computed) ;
	 // (MR) Note that now label_scales and label_offsets are the sequential results
	 // not the results of this single call, this might mess up things furhter down
	 // when they are used? Needs to be checked!
		got_scales = 1;
		    
    Gdiag = old_diag ;
    if (write_gca_fname) {
      printf("writing normalized gca to %s...\n", write_gca_fname) ;
      GCAwrite(gcam->gca, write_gca_fname) ;
    }
    start_t = parms.start_t ;
    *(&parms) = *(&saved_parms) ;
    parms.start_t = start_t ;
#if 0
    parms.l_dtrans = l_dtrans ;
    parms.tol /= 5 ;  // reset parameters to previous level
    parms.l_log_likelihood = 0.0f ;
    parms.l_smoothness /= 5 ;
#endif
  }
  
  // build multi-scale atlases
  FileNameOnly(out_fname, fname) ;
  dot = strrchr(fname, '.') ;
  if (dot)
    *dot = 0 ;
  GCAinsertLabels(gca, mri_inputs, transform, ninsertions, 
                  insert_labels, insert_intensities, insert_coords, insert_whalf) ;
  if (read_ms_name)
  {
    gcam_ms = GCAMMSread(read_ms_name, gcam, nsigmas != NSIGMAS? nsigmas : -1) ;
    if (nogibbs)
    {
      int i ;
      for (i = 0 ; i < gcam_ms->nsigmas ; i++)
        GCAfreeGibbs(gcam_ms->gcas[i]) ;
    }

    if (nsigmas < NSIGMAS)
    {
      int i ;
      for (i = nsigmas ; i < gcam_ms->nsigmas ; i++)
        GCAfree(&gcam_ms->gcas[i]) ;
      gcam_ms->nsigmas = nsigmas ;
    }
  }
  else
    gcam_ms = GCAMMSalloc(gcam, gca, sigmas, nsigmas, write_ms_name) ;
  if (renormalize_align_after) // renormalize gcamms
  {
    int i ;
	 // (MR) Note that label_scales and label_offsets are possibly the sequential results
	 // not the results of the single call in align_after,
	 // this might mess up things here? Needs to be checked!
    for (i = 0 ; i < gcam_ms->nsigmas ; i++)
      GCAapplyRenormalization(gcam_ms->gcas[i], label_scales, label_offsets, 0) ;
  }

  if (write_ms_name)
  {
    GCAMMSwrite(gcam_ms, write_ms_name) ;
    exit(0) ;
  }

  //////////////////////////////////////////////////////////////
  // -renorm fname option
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
    for (i = 0 ; i < gcam_ms->nsigmas ; i++)
      GCArenormalizeIntensities(gcam_ms->gcas[i], labels, intensities, nlines) ;
      
    free(labels) ;
    free(intensities) ;
  }

  parms.ndtrans = NDTRANS_LABELS ;
  parms.dtrans_labels = dtrans_labels ;
  if (ninsertions > 0)
  {
    int i, found, l  ;

#if 0
    for (i = 0 ; i < gcam_ms->nsigmas ; i++)
      GCAinsertLabels(gcam_ms->gcas[i], mri_inputs, transform, ninsertions, 
                      insert_labels, insert_intensities, insert_coords, insert_whalf) ;
#endif

    // now add label to dtrans list if it's not already in there
    for (found = i = 0 ; i < parms.ndtrans ; i++)
    {
      for (l = 0 ; l < ninsertions ; l++)
        if (insert_labels[l] == parms.dtrans_labels[i])
          found++ ;
    }
    if (found < ninsertions)
    {
      int *dtrans_labels ;
      dtrans_labels = (int *)calloc(parms.ndtrans+(ninsertions-found), sizeof(int)) ;
      for (i = 0 ; i < parms.ndtrans ; i++)
        dtrans_labels[i] = parms.dtrans_labels[i] ;

      for (l = 0 ; l < ninsertions ; l++)
      {
        found = 0 ; // see if it's already in the list
        for (i = 0 ; i < parms.ndtrans ; i++)
          if (insert_labels[l] == parms.dtrans_labels[i])
          {
            found = 1 ;
            break ;
          }
        if (found == 0)
        {
          dtrans_labels[parms.ndtrans] = insert_labels[l] ;
          parms.ndtrans++ ;
        }
      }
      parms.dtrans_labels = dtrans_labels ; // new expanded list
    }
  }
  
  //////////////////////////////////////////////////////////////////
  // here is the main work force
  //  log_posterior = GCAMlogPosterior(gcam, mri_inputs)/(gcam->width*gcam->height*gcam->depth) ;
  parms.mri_diag = mri_inputs ;
  min_scale = mean_scale = gcam_ms->sigmas[gcam_ms->nsigmas-1] ;
  round = reductions = 0 ;
  replace_labels(mri_atlas_labels, mri_atlas_labels, combine_labels, NCOMBINE_LABELS, &parms) ;
  mri_atlas_dist_map = MRIcreateDistanceTransforms(mri_atlas_labels, NULL,
                                                   2*sigmas[nsigmas-1],
                                                   parms.dtrans_labels,
                                                   parms.ndtrans) ;
  done = 0 ;
  gca_log_posterior = 0.0 ;
  do
  {
    TRANSFORM _transform, *gcam_transform = &_transform ;
    char      fname[STRLEN] ;
    MRI       *mri_aseg, *mri_pvals, *mri_sigma, *mri_aseg_orig, *mri_s_index ;
    gcam_transform->type = MORPH_3D_TYPE;
    gcam_transform->xform = (void *)gcam ;
    if (round  > 0 || 1)
    {
      double maxp ;

      mri_aseg_orig = MRIallocSequence(mri_inputs->width, 
                                       mri_inputs->height, mri_inputs->depth, MRI_FLOAT, 4) ;
      MRIcopyHeader(mri_inputs, mri_aseg_orig);
      mri_pvals = MRIcloneDifferentType(mri_inputs, MRI_FLOAT) ;
      mri_s_index = MRIcloneDifferentType(mri_inputs, MRI_UCHAR) ;
      mri_sigma = GCAMMScomputeOptimalScale(gcam_ms, gcam_transform, mri_inputs,NULL, 
                                            mri_aseg_orig, mri_pvals, mri_s_index);
      
      GCAMSreclassifyUsingGibbsPriors(mri_inputs,
                                      gcam_ms,
                                      mri_aseg_orig,
                                      gcam_transform,
                                      5, 0,
                                      NULL, mri_s_index) ;
      log_posterior = GCAMMSgibbsImageLogPosterior(gcam_ms, mri_aseg_orig, mri_inputs,
                                                   gcam_transform, mri_s_index);
      log_posterior /= (mri_inputs->width*mri_inputs->height*mri_inputs->depth) ;
      mri_aseg = replace_labels(mri_aseg_orig, NULL, combine_labels, NCOMBINE_LABELS, NULL) ;
#if 0
      gca_log_posterior = 
        GCAimagePosteriorLogProbability(gcam_ms->gcas[0], mri_aseg_orig, 
                                        mri_inputs, gcam_transform)/nvox ;
#endif
      sprintf(fname, "%s_aseg%d.mgh", parms.base_name, round) ; 
      MRIwriteFrame(mri_aseg, fname, 0) ;
      if (mri_aseg->nframes > 1)
      {
        sprintf(fname, "%s_aseg%d_p1.mgh", parms.base_name, round) ; 
        MRIwriteFrame(mri_aseg, fname, 1) ;
        sprintf(fname, "%s_aseg%d_label2.mgh", parms.base_name, round) ; 
        MRIwriteFrame(mri_aseg, fname, 2) ;
        sprintf(fname, "%s_aseg%d_p2.mgh", parms.base_name, round) ; 
        MRIwriteFrame(mri_aseg, fname, 3) ;
      }

      sprintf(fname, "%s_sigma%d.mgh", parms.base_name, round) ; 
      MRIwrite(mri_sigma, fname) ;
      sprintf(fname, "%s_pvals%d.mgh", parms.base_name, round) ; 
      MRIwrite(mri_pvals, fname) ;
      maxp = MRImaxNorm(mri_pvals) ;
      MRIscalarMul(mri_pvals, mri_pvals, 1.0/maxp) ;
      sprintf(fname, "%s_norm_pvals%d.mgh", parms.base_name, round) ; 
      MRIwrite(mri_pvals, fname) ;
    }
    else
    {
      sprintf(fname, "%s_aseg%d.mgh", parms.base_name, round) ; 
      mri_aseg = MRIread(fname) ;
      sprintf(fname, "%s_sigma%d.mgh", parms.base_name, round) ; 
      mri_sigma = MRIread(fname) ;
      sprintf(fname, "%s_pvals%d.mgh", parms.base_name, round) ; 
      mri_pvals = MRIread(fname) ;
    }
    last_mean_scale = mean_scale;
    mean_scale = MRImeanFrame(mri_sigma, 0) ;
    mean_scale = MRImeanFrameNonzeroMask(mri_sigma, 0, mri_aseg) ;
    printf("round %d, log posterior = %f, mean_scale = %2.2f (min %2.2f)\n", 
           round, log_posterior, mean_scale, min_scale) ;
    if (min_scale > mean_scale)
      mean_scale = min_scale ;
    parms.sigma = mean_scale ;
    min_scale = .95*mean_scale ; // ensure it decreases

    last_log_posterior = log_posterior ;
    last_gca_log_posterior = gca_log_posterior ;
    if (Gdiag & DIAG_WRITE)
    {
      sprintf(fname, "%s.log", parms.base_name) ;
      if (parms.log_fp == NULL)
      {
        if (parms.start_t == 0)
          parms.log_fp = fopen(fname, "w") ;
        else
          parms.log_fp = fopen(fname, "a") ;
      }
      fprintf(parms.log_fp, "round %d, log posterior = %f, mean scale=%2.3f (%2.3f)\n", 
              round, log_posterior, mean_scale, min_scale) ;
      fflush(parms.log_fp) ;
    }
    else
      parms.log_fp = NULL ;
    if (demons)
    {
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        MRI *mri_tmp ;
        mri_tmp = TransformApplyType(transform, mri_aseg, NULL, SAMPLE_NEAREST) ;
        MRIwrite(mri_tmp, "sl.mgh") ;
        MRIfree(&mri_tmp) ;
      }
      if (use_pvals)
        GCAMdemonsRegister(gcam, mri_aseg, mri_atlas_labels, &parms,
                           2*sigmas[nsigmas-1], mri_pvals, mri_atlas_dist_map) ;
      else
        GCAMdemonsRegister(gcam, mri_aseg, mri_atlas_labels, &parms,
                           2*sigmas[nsigmas-1], NULL, mri_atlas_dist_map) ;
      if (Gdiag & DIAG_WRITE/* && DIAG_VERBOSE_ON*/)
      {
        MRI *mri_morphed ;
        char fname[STRLEN] ;

        sprintf(fname, "%s_round%d.m3z", parms.base_name, round+1) ;
        printf("writing %s\n", fname) ;
        GCAMwrite(gcam, fname) ;

        mri_morphed = GCAMmorphToAtlas(mri_inputs, gcam, NULL, 0, SAMPLE_TRILINEAR) ;
        sprintf(fname, "%s_src_to_atlas%d.mgh", parms.base_name,round+1) ;
        printf("writing %s\n", fname) ;
        MRIwrite(mri_morphed, fname) ; MRIfree(&mri_morphed) ;

        mri_morphed = GCAMmorphFromAtlas(mri_atlas_labels, gcam, NULL, SAMPLE_NEAREST);
        sprintf(fname, "%s_atlas_to_src%d.mgh", parms.base_name, round+1) ;
        printf("writing %s\n", fname) ;
        MRIwrite(mri_morphed, fname) ; MRIfree(&mri_morphed) ;
      }
    }
    else
      GCAMregister(gcam, mri_inputs, &parms) ;
#if 0
    log_posterior = GCAMlogPosterior(gcam, mri_inputs)/(gcam->width*gcam->height*gcam->depth) ;
#else
      log_posterior = GCAMMSgibbsImageLogPosterior(gcam_ms, mri_aseg, mri_inputs,
                                                   gcam_transform, mri_s_index);
      log_posterior /= (mri_inputs->width*mri_inputs->height*mri_inputs->depth) ;
#endif
    gca_log_posterior = 
      GCAimagePosteriorLogProbability(gcam_ms->gcas[0], mri_aseg_orig, 
                                      mri_inputs, gcam_transform)/nvox ;

    MRIfree(&mri_aseg) ; MRIfree(&mri_pvals) ; MRIfree(&mri_sigma); MRIfree(&mri_aseg_orig) ;
    MRIfree(&mri_s_index) ;
    if (Gdiag & DIAG_WRITE)
    {
      char fname[STRLEN] ;
      sprintf(fname, "%s.log", parms.base_name) ;

      fprintf(parms.log_fp, "round %d, old posterior %2.4f (%2.4f), new posterior %2.4f (%2.4f)\n",
              round, last_log_posterior, last_gca_log_posterior,
              log_posterior, gca_log_posterior) ;
      printf("round %d, old posterior %2.4f (%2.4f), new posterior %2.4f (%2.4f)\n",
              round, last_log_posterior, last_gca_log_posterior,
              log_posterior, gca_log_posterior) ;
      fflush(parms.log_fp) ; fflush(stdout) ;
      fclose(parms.log_fp) ;
      parms.log_fp = NULL ;
    }
    round++ ;
    if (round > max_rounds)
      done = 1 ;
    else if (log_posterior <= last_log_posterior+(last_log_posterior*2*parms.tol))
    {
      min_scale *= .9 ;
      reductions++ ;
      if (min_scale < 1 || reductions > max_reductions)
        done = 1;
      if (Gdiag & DIAG_WRITE)
      {
        char fname[STRLEN] ;
        sprintf(fname, "%s.log", parms.base_name) ;
        
        fprintf(parms.log_fp, "posterior increased, decreasing min_scale to %2.2f\n", min_scale) ;
        printf("posterior increased, decreasing min_scale to %2.2f\n", min_scale) ;
        fflush(parms.log_fp) ; fflush(stdout) ;
        fclose(parms.log_fp) ;
        parms.log_fp = NULL ;
      }
    }
    else
      done = 0 ;
    //    done = 1 ;  
  } while (!done) ;

  //record GCA filename to gcam
  strcpy(gcam->atlas.fname, gca_fname);
  printf("writing output transformation to %s...\n", out_fname) ;
  if (vf_fname)
    write_vector_field(mri_inputs, gcam, vf_fname) ;
  // GCAMwrite is used not MORPH3D
  if (GCAMwrite(gcam, out_fname) != NO_ERROR)
    ErrorExit(Gerror, "%s: GCAMwrite(%s) failed", Progname, out_fname) ;

#if 0
  if (gca)
    GCAfree(&gca) ;
#endif
  GCAMfree(&gcam) ;
  if (mri_inputs)
    MRIfree(&mri_inputs) ;
  if (diag_fp)
    fclose(diag_fp) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  hours = minutes / (60) ;
  minutes = minutes % 60 ;
  seconds = seconds % 60 ;
  printf("registration took %d hours, %d minutes and %d seconds.\n",
         hours, minutes, seconds) ;
  exit(0) ;
  return(0) ;
}



/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0, err ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  /*  StrUpper(option) ;*/
  if (!stricmp(option, "DIST") || !stricmp(option, "DISTANCE")) {
    parms.l_distance = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_dist = %2.2f\n", parms.l_distance) ;
  } else if (!stricmp(option, "REGULARIZE")) {
    regularize = atof(argv[2]) ;
    printf("regularizing variance to be sigma+%2.1fC(noise)\n", regularize) ;
    nargs = 1 ;
  } else if (!stricmp(option, "write_gca")) {
    write_gca_fname = argv[2] ;
    printf("writing gca to file name %s\n", write_gca_fname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "MIN_AVGS")) {
    parms.min_avgs = atoi(argv[2]) ;
    printf("setting min # of averages to %d\n", parms.min_avgs) ;
    nargs = 1 ;
  } else if (!stricmp(option, "READ_MS")) {
    read_ms_name = argv[2] ;
    printf("reading previously constructed multiscale GCA from %s\n", read_ms_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "READ_DIST")) {
    read_dist_name = argv[2] ;
    printf("reading previously constructed distance transforms %s\n", read_dist_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "WRITE_DIST")) {
    write_dist_name = argv[2] ;
    printf("writing distance transforms to %s\n", write_dist_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "READ_MS")) {
    read_ms_name = argv[2] ;
    printf("reading previously constructed multiscale GCA from %s\n", read_ms_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "MAX_ROUNDS")) {
    max_rounds = atoi(argv[2]) ;
    printf("setting max number of rounds to %d\n", max_rounds) ;
    nargs = 1 ;
  } else if (!stricmp(option, "MAX_REDUCTIONS")) {
    max_reductions = atoi(argv[2]) ;
    printf("setting max number of rounds to %d\n", max_rounds) ;
    nargs = 1 ;
  } else if (!stricmp(option, "DEMONS")) {
    demons = 1 ;
    printf("using demons algorithm\n") ;
  } else if (!stricmp(option, "USE_PVALS")) {
    use_pvals = 1 ;
    printf("using pvals in demons algorithm\n") ;
  } else if (!stricmp(option, "NSIGMAS")) {
    nsigmas = atoi(argv[2]) ;
    if (nsigmas > NSIGMAS)
      nsigmas = NSIGMAS ;
    printf("cropping sigmas to %d of %d\n", nsigmas, (int)NSIGMAS) ;
    nargs = 1 ;
  } else if (!stricmp(option, "WRITE_MS")) {
    write_ms_name = argv[2] ;
    printf("writing multiscale GCA from %s\n", write_ms_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "REGULARIZE_MEAN")) {
    regularize_mean = atof(argv[2]) ;
    printf("regularizing means to be %2.2f u(global) + %2.2f u(r)\n",
           regularize_mean, 1-regularize_mean) ;
    nargs = 1 ;
  } else if (!stricmp(option, "scale_smoothness")) {
    parms.scale_smoothness = atoi(argv[2]) ;
    parms.npasses = 2 ;
    printf("%sscaling smooothness coefficient (default=1), "
           "and setting npasses=%d\n",
           parms.scale_smoothness ? "" : "not ", parms.npasses) ;
    nargs = 1 ;
  } else if (!stricmp(option, "NOBRIGHT")) {
    remove_bright = 1 ;
    printf("removing bright non-brain structures...\n") ;
  } else if (!stricmp(option, "bigventricles")) {
    handle_expanded_ventricles = 1 ;
    printf("handling expanded ventricles...\n") ;
  } else if (!stricmp(option, "nogibbs")) {
    nogibbs = 1 ;
    printf("not using Gibbs priors in segmentation\n") ;
  } else if (!stricmp(option, "nobigventricles")) {
    handle_expanded_ventricles = 0 ;
    printf("not handling expanded ventricles...\n") ;
  } else if (!stricmp(option, "uncompress")) {
    parms.uncompress = 1 ;
    printf("parms.uncompress=1...\n") ;
  } else if (!stricmp(option, "secondpassrenorm")) {
    do_secondpass_renorm = 1 ;
    printf("performing 2nd-pass renormalization...\n") ;
  } else if (!stricmp(option, "renormalize")) {
    renormalize = 1 ;
    printf("renormalizing GCA to MAP estimate of means\n") ;
  } else if (!stricmp(option, "read_lta")) {
    read_lta = 1 ;
    printf("reading LTA from <base-name>.lta\n") ;
  } else if (!stricmp(option, "SMOOTH") || !stricmp(option, "SMOOTHNESS")) {
    parms.l_smoothness = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_smoothness = %2.2f\n", parms.l_smoothness) ;
  } else if (!stricmp(option, "SAMPLES")) {
    sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing control points to %s...\n", sample_fname) ;
  } else if (!stricmp(option, "SMALL") || !stricmp(option, "NSMALL")) {
    parms.nsmall = atoi(argv[2]) ;
    nargs = 1 ;
    printf("allowing %d small steps before terminating integration\n",
           parms.nsmall) ;
  } else if (!stricmp(option, "FIXED")) {
    parms.integration_type = GCAM_INTEGRATE_FIXED ;
    printf("using fixed time-step integration\n") ;
  } else if (!stricmp(option, "OPTIMAL")) {
    parms.integration_type = GCAM_INTEGRATE_OPTIMAL ;
    printf("using optimal time-step integration\n") ;
  } else if (!stricmp(option, "NONEG")) {
    parms.noneg = atoi(argv[2]) ;
    nargs = 1 ;
    printf("%s allowing temporary folds during numerical minimization\n",
           parms.noneg ? "not" : "") ;
  } else if (!stricmp(option, "ISIZE") || !stricmp(option, "IMAGE_SIZE")) {
    IMAGE_SIZE = atoi(argv[2]) ;
    nargs = 1 ;
    printf("setting diagnostic image size to %d\n", IMAGE_SIZE) ;
  } else if (!stricmp(option, "WM")) {
    register_wm_flag = 1 ;
    printf("registering white matter in initial pass...\n") ;
  } else if (!stricmp(option, "TL")) {
    tl_fname = argv[2] ;
    nargs = 1 ;
    printf("reading temporal lobe atlas from %s...\n", tl_fname) ;
  } else if (!stricmp(option, "RELABEL")) {
    parms.relabel = atoi(argv[2]) ;
    nargs = 1 ;
    printf("%srelabeling nodes with MAP estimates\n",
           parms.relabel ? "" : "not ") ;
  } else if (!stricmp(option, "RELABEL_AVGS")) {
    parms.relabel_avgs = atoi(argv[2]) ;
    nargs = 1 ;
    printf("relabeling nodes with MAP estimates at avgs=%d\n",
           parms.relabel_avgs) ;
  } else if (!stricmp(option, "RESET_AVGS")) {
    parms.reset_avgs = atoi(argv[2]) ;
    nargs = 1 ;
    printf("resetting metric properties at avgs=%d\n", parms.reset_avgs) ;
  } else if (!stricmp(option, "RESET")) {
    reset = 1 ;
    printf("resetting metric properties...\n") ;
  } else if (!stricmp(option, "VF")) {
    vf_fname = argv[2] ;
    nargs = 1 ;
    printf("writing vector field to %s...\n", vf_fname) ;
  } else if (!stricmp(option, "INSERT")) {
    if (ninsertions >= MAX_INSERTIONS)
      ErrorExit(ERROR_NOMEMORY, "%s: too many insertions (%d) specified\n", Progname, ninsertions) ;

    insert_labels[ninsertions] = atoi(argv[2]) ;
    insert_intensities[ninsertions] = atoi(argv[3]) ;
    insert_coords[ninsertions][0] = atoi(argv[4]) ;
    insert_coords[ninsertions][1] = atoi(argv[5]) ;
    insert_coords[ninsertions][2] = atoi(argv[6]) ;
    insert_whalf[ninsertions] = atoi(argv[7]) ;
    printf("inserting label %d (%s) at (%d, %d, %d) with intensity = %d, in %d voxel nbhd\n",
           insert_labels[ninsertions],
           cma_label_to_name(insert_labels[ninsertions]),
           insert_coords[ninsertions][0],
           insert_coords[ninsertions][1],
           insert_coords[ninsertions][2],
           insert_intensities[ninsertions],
           insert_whalf[ninsertions]) ;
           
    ninsertions++ ;
    nargs = 6 ;
  } else if (!stricmp(option, "MASK")) {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
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
    map_to_flash = 1 ;
    printf("using FLASH forward model to predict intensity values...\n") ;
  } else if (!stricmp(option, "FLASH_PARMS")) {
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
  } else if (!stricmp(option, "NOVAR")) {
    novar = 1 ;
    printf("not using variance estimates\n") ;
  } else if (!stricmp(option, "USEVAR")) {
    novar = 0 ;
    printf("using variance estimates\n") ;
  } else if (!stricmp(option, "DT")) {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    printf("dt = %2.2e\n", parms.dt) ;
  } else if (!stricmp(option, "DTRANS")) {
    parms.l_dtrans = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_dtrans = %2.2e\n", parms.l_dtrans) ;
  } else if (!stricmp(option, "TOL")) {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("tol = %2.2e\n", parms.tol) ;
  } else if (!stricmp(option, "CENTER")) {
    center = 1 ;
    printf("using GCA centroid as origin of transform\n") ;
  } else if (!stricmp(option, "NOSCALE")) {
    noscale = 1 ;
    printf("disabling scaling...\n") ;
  } else if (!stricmp(option, "LEVELS")) {
    parms.levels = atoi(argv[2]) ;
    nargs = 1 ;
    printf("levels = %d\n", parms.levels) ;
  } else if (!stricmp(option, "LIKELIHOOD")) {
    parms.l_likelihood = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_likelihood = %2.2f\n", parms.l_likelihood) ;
  } else if (!stricmp(option, "LOGLIKELIHOOD")) {
    parms.l_log_likelihood = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_log_likelihood = %2.2f\n", parms.l_log_likelihood) ;
  } else if (!stricmp(option, "LABEL")) {
    parms.l_label = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_label = %2.2f\n", parms.l_label) ;
  } else if (!stricmp(option, "MAP")) {
    parms.l_map = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_map = %2.2f\n", parms.l_map) ;
  } else if (!stricmp(option, "LDIST") || !stricmp(option, "LABEL_DIST")) {
    parms.label_dist = atof(argv[2]) ;
    nargs = 1 ;
    printf("label_dist = %2.2f\n", parms.label_dist) ;
  } else if (!stricmp(option, "REDUCE")) {
    nreductions = atoi(argv[2]) ;
    nargs = 1 ;
    printf("reducing input images %d times before aligning...\n",
           nreductions) ;
  } else if (!stricmp(option, "DEBUG_NODE")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging node (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx2 = atoi(argv[2]) ;
    Gy2 = atoi(argv[3]) ;
    Gz2 = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gvx, Gvy, Gvz) ;
  } else if (!stricmp(option, "norm")) {
    norm_fname = argv[2] ;
    nargs = 1 ;
    printf("intensity normalizing and writing to %s...\n",norm_fname);
  } else if (!stricmp(option, "avgs")) {
    avgs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,
            "applying mean filter %d times to conditional densities...\n",
            avgs) ;
  } else if (!stricmp(option, "cross-sequence") ||
             !stricmp(option, "cross_sequence")) {
    regularize = .5 ;
    avgs = 2 ;
    renormalize = 1 ;
    printf("registering sequences, equivalent to:\n") ;
    printf("\t-renormalize\n\t-avgs %d\n\t-regularize %2.3f\n",
           avgs, regularize) ;
  } else if (!stricmp(option, "align-cross-sequence") ||
             !stricmp(option, "align") ||
             !stricmp(option, "align-after")) {
    regularize = .5 ;
    reinit = 1 ;
    regularize_mean = .5 ;
    parms.ratio_thresh = 0.000001 ;
    //  avgs = 2 ;   not used anymore
    renormalize_align = 1 ;
    printf("renormalizing sequences with structure "
           "alignment, equivalent to:\n") ;
    printf("\t-renormalize\n\t-regularize_mean %2.3f\n\t-regularize %2.3f\n",
           regularize_mean, regularize) ;
    if (!stricmp(option, "align-after")) {
      renormalize_align = 0 ;
      renormalize_align_after = 1 ;
    }

  } else if (!stricmp(option, "no-re-init") ||
             !stricmp(option, "no-reinit") || 
             !stricmp(option, "no_re_init") ) {
    reinit = 0; //donot reinitialize GCAM with the multiple linear registration
  } else if (!stricmp(option, "cross-sequence-new") ||
             !stricmp(option, "cross_sequence_new")) {
    regularize = .5 ;
    avgs = 2 ;
    renormalize_new = 1 ;
    printf("registering sequences, equivalent to:\n") ;
    printf("\t-renormalize\n\t-avgs %d\n\t-regularize %2.3f\n",
           avgs, regularize) ;
  } else if (!stricmp(option, "area")) {
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    printf("using l_area=%2.3f\n", parms.l_area) ;
  } else if (!stricmp(option, "rthresh")) {
    parms.ratio_thresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("using compression ratio threshold = %2.3f...\n",
           parms.ratio_thresh) ;
  } else if (!stricmp(option, "invert-and-save")) {
    printf("Loading, Inverting, Saving, Exiting ...\n");
    err = GCAMwriteInverse(argv[2],NULL);
    exit(err);
  } else switch (toupper(*option)) {
    case 'L':   /* for longitudinal analysis */
      xform_name = argv[2] ;
      //invert is not needed if REG is from tp1 to current subject! -xh
      //   inverted_xform = 1 ;
      long_reg_fname = argv[3] ;
      nargs = 2 ;
      printf("reading previously computed atlas xform %s "
             "and applying inverse registration %s\n",
             xform_name, long_reg_fname) ;
      break ;
    case 'J':
      parms.l_jacobian = atof(argv[2]) ;
      nargs = 1 ;
      printf("using l_jacobian=%2.3f\n", parms.l_jacobian) ;
      break ;
    case 'A':
      parms.navgs = atoi(argv[2]) ;
      nargs = 1 ;
      printf("smoothing gradient with %d averages...\n", parms.navgs) ;
      break ;
    case 'F':
      ctl_point_fname = argv[2] ;
      nargs = 1 ;
      printf("reading manually defined control points from %s\n",
             ctl_point_fname) ;
      break ;
    case 'X':
      xform_name = argv[2] ;
      nargs = 1 ;
      printf("reading previous transform from %s...\n", xform_name) ;
      break ;
    case 'K':
      parms.exp_k = atof(argv[2]) ;
      printf("setting exp_k to %2.2f (default=%2.2f)\n",
             parms.exp_k, EXP_K) ;
      nargs = 1 ;
      break ;
    case 'T':
      transform = TransformRead(argv[2]) ;
      if (!transform)
        ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s",
                  Progname, argv[2]) ;
      if (transform->type == LINEAR_RAS_TO_RAS)
        ErrorExit(ERROR_BADPARM,
                  "%s: transform %s is RAS to RAS, cannot be used\n",
                  Progname, argv[2]);
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
      parms.sigma = atof(argv[2]) ;
      printf("using sigma=%2.3f as upper bound on blurring.\n",
             parms.sigma) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      printf("usage: %s <in volume> <template volume> <output transform>\n",
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

static int
write_vector_field(MRI *mri, GCA_MORPH *gcam, char *vf_fname) {
  FILE            *fp ;
  int             x, y, z ;
  GCA_MORPH_NODE  *gcamn ;

  fp = fopen(vf_fname, "w") ;

  for (x = 0 ; x < gcam->width ; x++) {
    for (y = 0 ; y < gcam->height ; y++) {
      for (z = 0 ; z < gcam->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
        fprintf(fp, "%f %f %f %f\n",
                gcamn->x-gcamn->origx,
                gcamn->y-gcamn->origy,
                gcamn->z-gcamn->origz,
                gcamn->gc ? gcamn->gc->means[0] : 0.0) ;
      }
    }
  }
  fclose(fp) ;
  return(NO_ERROR) ;
}

static int
remove_bright_stuff(MRI *mri, GCA *gca, TRANSFORM *transform) {
  HISTO            *h, *hs ;
  int              peak, num, end, x, y, z, xi, yi, zi,
  xk, yk, zk, i, n, erase, five_mm ;
  float            thresh ;
  Real             val, new_val ;
  MRI              *mri_tmp, *mri_nonbrain, *mri_tmp2 ;
  GCA_PRIOR        *gcap ;
  MRI_SEGMENTATION *mriseg ;
  MRI_SEGMENT      *mseg ;
  MSV              *msv ;

  if (gca->ninputs > 1)
    return(NO_ERROR) ;

  mri_tmp = MRIalloc(mri->width, mri->height, mri->depth, MRI_UCHAR) ;
  mri_nonbrain = MRIalloc(mri->width, mri->height, mri->depth, MRI_UCHAR) ;
  for (x = 0 ; x < mri->width ; x++) {
    for (y = 0 ; y < mri->height ; y++) {
      for (z = 0 ; z < mri->depth ; z++) {
        gcap = getGCAP(gca, mri, transform, x, y, z) ;
        if (gcap->nlabels == 0 ||
            (gcap->nlabels == 1 && IS_UNKNOWN(gcap->labels[0])))
          MRIvox(mri_nonbrain, x, y, z) = 1 ;
        else {
          MRIsampleVolume(mri, x, y, z, &val) ;
          if (FZERO(val))
            MRIvox(mri_nonbrain, x, y, z) = 128 ;
        }
      }
    }
  }
  /* dilate it by 0.5 cm */
  five_mm = nint(5.0*pow(mri->xsize*mri->ysize*mri->zsize, 1.0f/3.0f)) ;
  for (i = 0 ; i < five_mm ; i++) {
    MRIdilate(mri_nonbrain, mri_tmp) ;
    MRIcopy(mri_tmp, mri_nonbrain) ;
  }

  MRIclear(mri_tmp) ;
  h = MRIhistogram(mri, 0) ;
  h->counts[0] = 0 ;
  hs = HISTOsmooth(h, NULL, 2) ;

  peak = HISTOfindLastPeak(hs, 5, 0.1) ;
  end = HISTOfindEndOfPeak(hs, peak, 0.01) ;
  thresh = hs->bins[end] ;
  new_val = 0 ;

  printf("removing voxels brighter than %2.1f\n", thresh) ;

  for (num = x = 0 ; x < mri->width ; x++) {
    for (y = 0 ; y < mri->height ; y++) {
      for (z = 0 ; z < mri->depth ; z++) {
        if (x == Gvx && y == Gvy && z == Gvz)
          DiagBreak() ;
        MRIsampleVolume(mri, x, y, z, &val) ;
        if (val > thresh) {
          num++ ;
          MRIvox(mri_tmp, x, y, z) = 128 ;
          /* MRIsetVoxVal(mri, x, y, z, 0, (float)new_val) ;*/
        }
      }
    }
  }


  /* relax threshold somewhat, and reduce voxels that are above this thresh
     and nbrs of one above the more stringent one.
  */
  end = HISTOfindStartOfPeak(hs, peak, 0.1) ;
  thresh = hs->bins[end] ;
  mri_tmp2 = MRIcopy(mri_tmp, NULL) ;
  for (x = 0 ; x < mri->width ; x++) {
    for (y = 0 ; y < mri->height ; y++) {
      for (z = 0 ; z < mri->depth ; z++) {
        if (MRIvox(mri_tmp2, x, y, z) == 0)
          continue ;
        for (xk = -1 ; xk <= 1 ; xk++) {
          xi = mri_tmp->xi[x+xk] ;
          for (yk = -1 ; yk <= 1 ; yk++) {
            yi = mri_tmp->yi[y+yk] ;
            for (zk = -1 ; zk <= 1 ; zk++) {
              zi = mri_tmp->zi[z+zk] ;
              if (xi == Gvx && yi == Gvy && zi == Gvz)
                DiagBreak() ;
              MRIsampleVolume(mri, xi, yi, zi, &val) ;
              if (val > thresh) {
                num++ ;
                MRIvox(mri_tmp, xi, yi, zi) = 128 ;
                /* MRIsetVoxVal
                   (mri, xi, yi, zi, 0, (float)new_val) ;*/
              }
            }
          }
        }
      }
    }
  }

  MRIfree(&mri_tmp2) ;
  mriseg = MRIsegment(mri_tmp, 1, 255) ;
  printf("%d bright voxels found - %d segments\n", num, mriseg->nsegments) ;

  for (num = i = 0 ; i < mriseg->nsegments ; i++) {
    /* check to see that at least one voxel in
       segment is in nonbrain mask (i.e. it is within 1cm of
       nonbrain */
    mseg = &mriseg->segments[i] ;
    for (erase = 0, n = 0 ; n < mseg->nvoxels ; n++) {
      msv = &mseg->voxels[n] ;
      if (msv->x == Gvx && msv->y == Gvy && msv->z == Gvz)
        DiagBreak() ;
      if (MRIvox(mri_nonbrain, msv->x, msv->y, msv->z) > 0) {
        erase = 1 ;
        break ;
      }
    }
    if (erase) {
      if (DIAG_VERBOSE_ON)
        printf("erasing segment %d (%d voxels) with centroid "
               "at (%2.0f, %2.0f, %2.0f)\n",
               i, mseg->nvoxels, mseg->cx, mseg->cy, mseg->cz) ;
      for (n = 0 ; n < mseg->nvoxels ; n++) {
        msv = &mseg->voxels[n] ;
        if (msv->x == Gx && msv->y == Gy && msv->z == Gz)
          DiagBreak() ;
        MRIsetVoxVal(mri, msv->x, msv->y, msv->z, 0, 0.0f) ;
        num++ ;
      }
    }
  }

  printf("%d bright voxels erased\n", num) ;
  HISTOfree(&h) ;
  HISTOfree(&hs) ;
  MRIfree(&mri_tmp) ;
  MRIfree(&mri_nonbrain) ;
  MRIsegmentFree(&mriseg) ;
  return(NO_ERROR) ;
}


