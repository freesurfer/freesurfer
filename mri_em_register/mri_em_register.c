

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
#include "utils.h"
#include "gca.h"
#include "cma.h"
#include "mrinorm.h"

char         *Progname ;
static MORPH_PARMS  parms ;

static char *norm_fname = NULL ;

static double TR = -1 ;
static double alpha = -1 ;
static double TE = -1 ;

static char *sample_fname = NULL ;
static char *transformed_sample_fname = NULL ;
static int novar = 0 ;

#define MIN_SPACING   2.0
static float min_spacing = MIN_SPACING ;

static int use_contrast = 0 ;
static float min_prior = MIN_PRIOR ;
static double tol = 0.001 ;
static double tx = 0.0 ;
static double ty = 0.0 ;
static double tz = 0.0 ;
static double rzrot = 0.0 ;
static double rxrot = 0.0 ;
static double ryrot = 0.0 ;

static FILE *diag_fp = NULL ;

static LTA *Glta = NULL ;

static int translation_only = 0 ;
static int get_option(int argc, char *argv[]) ;
static int register_mri(MRI *mri_in, GCA *gca, MP *parms, int passno) ;

static char *renormalization_fname = NULL ;
static char *tissue_parms_fname = NULL ;
static int center = 1 ;
static int nreductions = 1 ;
static int noscale = 0 ;
static int num_xforms = 1 ;
static int transform_loaded = 0 ;
static char *gca_mean_fname = NULL ;

static MATRIX *find_optimal_transform(MRI *mri_in, GCA *gca, GCA_SAMPLE *gcas,
                                      int nsamples, MATRIX *m_L, int passno,
                                      int write_iterations) ;
static double find_optimal_translation(GCA *gca, GCA_SAMPLE *gcas, MRI *mri, 
                                       int nsamples, MATRIX *m_L, 
                                       float min_trans, float max_trans, 
                                       float trans_steps, int nreductions) ;
static double find_optimal_scaling(GCA *gca, GCA_SAMPLE *gcas, MRI *mri, 
                                   int nsamples, MATRIX *m_L, MATRIX *m_origin,
                                   float min_scale, float max_scale, 
                                   float scale_steps, int nreductions) ;
static double find_optimal_rotation(GCA *gca, GCA_SAMPLE *gcas, MRI *mri, 
                                 int nsamples, MATRIX *m_L, MATRIX *m_origin,
                                 float min_angle, float max_angle,
                                 float angle_steps,
                                 int nreductions) ;
static double find_optimal_scaling_and_rotation(GCA *gca, GCA_SAMPLE *gcas, 
                                                MRI *mri, 
                                 int nsamples, MATRIX *m_L, MATRIX *m_origin,
                                 float min_angle, float max_angle,
                                 float min_scale, float max_scale,
                                                float angle_steps, 
                                                float scale_steps,
                                                int nreductions) ;
static double blur_sigma = 0.0f ;

/* 
   command line consists of three inputs:

   argv[1]  - directory containing 'canonical' brain
   argv[2]  - directory containing brain to be registered
   argv[3]  - directory in which to write out registered brain.
*/

#define NPARMS     12
#define NSAMPLES   (NPARMS*20)

static int nsamples = NSAMPLES ;
int
main(int argc, char *argv[])
{
  char         *gca_fname, *in_fname, *out_fname, fname[STRLEN], **av ;
  MRI          *mri_in ;
  GCA          *gca /*, *gca_tmp, *gca_reduced*/ ;
  int          ac, nargs, i, done ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  float        old_log_p, log_p ;

  parms.l_intensity = 1.0f ;
  parms.niterations = 100 ;
  parms.levels = -1 ;   /* use default */
  parms.dt = 1e-6 ;  /* was 5e-6 */
  parms.tol = INTEGRATION_TOL*5 ;

  parms.max_levels = 0 ;
  parms.dt = 5e-6 ;  /* was 5e-6 */
  parms.tol = 1e-5 ;
  parms.momentum = 0.8 ;
  parms.niterations = 25 ;
  Progname = argv[0] ;


  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    ErrorExit(ERROR_BADPARM, 
              "usage: %s <in brain> <template> <output file name>\n",
              Progname) ;

  in_fname = argv[1] ;
  gca_fname = argv[2] ;
  out_fname = argv[3] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
  Gdiag |= DIAG_WRITE ;
  printf("logging results to %s.log\n", parms.base_name) ;

  TimerStart(&start) ;
  printf("reading '%s'...\n", gca_fname) ;
  fflush(stdout) ;
  gca = GCAread(gca_fname) ;
  if (gca == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open GCA %s.\n",
              Progname, gca_fname) ;


  if (novar)
    GCAunifyVariance(gca) ;

  if (renormalization_fname)
  {
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
    while (cp)
    {
      nlines++ ;
      cp = fgetl(line, 199, fp) ;
    }
    rewind(fp) ;
    printf("reading %d labels from %s...\n", nlines,renormalization_fname) ;
    labels = (int *)calloc(nlines, sizeof(int)) ;
    intensities = (float *)calloc(nlines, sizeof(float)) ;
    cp = fgetl(line, 199, fp) ;
    for (i = 0 ; i < nlines ; i++)
    {
      sscanf(cp, "%e  %e", &f1, &f2) ;
      labels[i] = (int)f1 ; intensities[i] = f2 ;
      if (labels[i] == Left_Cerebral_White_Matter)
        DiagBreak() ;
      cp = fgetl(line, 199, fp) ;
    }
    GCArenormalizeIntensities(gca, labels, intensities, nlines) ;
    free(labels) ; free(intensities) ;
  }


#if 0
  if (gca->spacing < min_spacing)
    printf(
            "reducing GCA to %d mm spacing before sampling "
            "interesting points...\n", (int)min_spacing) ;
  gca_reduced = gca ;
  while (gca_reduced->spacing < min_spacing)
  {
    gca_tmp = GCAreduce(gca_reduced) ;
    if (gca_reduced != gca)
      GCAfree(&gca_reduced) ;
    gca_reduced = gca_tmp ;
  }
  parms.gcas = GCAfindStableSamplesByLabel(gca_reduced, nsamples) ;
  if (gca_reduced != gca)
    GCAtransformSamples(gca_reduced, gca, parms.gcas, nsamples) ;
#else
  if (use_contrast)
    parms.gcas = GCAfindContrastSamples(gca,&nsamples,
                                        (int)min_spacing,min_prior);
  else
    parms.gcas = GCAfindStableSamples(gca, &nsamples,
                                      (int)min_spacing,min_prior);
  printf("using %d sample points...\n", nsamples) ;
#endif
  parms.nsamples = nsamples ;


  printf("reading '%s'...\n", in_fname) ;
  fflush(stdout) ;
  mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
              Progname, in_fname) ;

  if (alpha > 0)
    mri_in->flip_angle = alpha ;
  if (TR > 0)
    mri_in->tr = TR ;
  if (TE > 0)
    mri_in->te = TE ;

  if (tissue_parms_fname)   /* use FLASH forward model */
    GCArenormalizeToFlash(gca, tissue_parms_fname, mri_in) ;

  if (parms.write_iterations != 0)
  {
    char fname[STRLEN] ;
    MRI  *mri_gca ;
    mri_gca = MRIclone(mri_in, NULL) ;
    GCAmri(gca, mri_gca) ;
    sprintf(fname, "%s_target", parms.base_name) ;
    MRIwriteImageViews(mri_gca, fname, 400) ;
    sprintf(fname, "%s_target.mgh", parms.base_name) ;
    printf("writing target volume to %s...\n", fname) ;
    MRIwrite(mri_gca, fname) ;
  }
  if (sample_fname)
  {
    GCAwriteSamples(gca, mri_in, parms.gcas, nsamples, sample_fname) ;
    printf("samples written\n") ;
  }
#if 0
  if (gca_reduced != gca)
    GCAfree(&gca_reduced) ;
#endif

  if (!FZERO(tx) || !FZERO(ty) || !FZERO(tz))
  {
    MRI *mri_tmp ;
    
    printf("translating second volume by (%2.1f, %2.1f, %2.1f)\n",
            tx, ty, tz) ;
    mri_tmp = MRItranslate(mri_in, NULL, tx, ty, tz) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
#if 1
  if (!FZERO(rzrot))
  {
    MRI *mri_tmp ;
    
    printf(
            "rotating second volume by %2.1f degrees around Z axis\n",
            (float)DEGREES(rzrot)) ;
    mri_tmp = MRIrotateZ_I(mri_in, NULL, rzrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  if (!FZERO(rxrot))
  {
    MRI *mri_tmp ;
    
    printf(
            "rotating second volume by %2.1f degrees around X axis\n",
            (float)DEGREES(rxrot)) ;
    mri_tmp = MRIrotateX_I(mri_in, NULL, rxrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  if (!FZERO(ryrot))
  {
    MRI *mri_tmp ;
    
    printf(
            "rotating second volume by %2.1f degrees around Y axis\n",
            (float)DEGREES(ryrot)) ;
    mri_tmp = MRIrotateY_I(mri_in, NULL, ryrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
#else
  if (!FZERO(ryrot) || !FZERO(rxrot) || !FZERO(rzrot))
  {
    MRI *mri_tmp ;
    MATRIX *mX, *mY, *mZ, *mRot, *mTmp ;
    
    mX = MatrixAllocRotation(3, x_angle, X_ROTATION) ;
    mY = MatrixAllocRotation(3, y_angle, Y_ROTATION) ;
    mZ = MatrixAllocRotation(3, z_angle, Z_ROTATION) ;
    mTmp = MatrixMultiply(mX, mZ, NULL) ;
    mRot = MatrixMultiply(mY, mTmp, NULL)
      printf(
              "rotating second volume by (%2.1f, %2.1f, %2.1f) degrees\n",
              (float)DEGREES(rxrot), (float)DEGREES(ryrot)
              (float)DEGREES(rzrot)) ;
    
    mri_tmp = MRIrotate_I(mri_in, NULL, mRot, NULL) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
    
    MatrixFree(&mX) ; MatrixFree(&mY) ; MatrixFree(&mZ) ; 
    MatrixFree(&mTmp) ; MatrixFree(&mRot) ;
  }
#endif

  if (!transform_loaded)   /* wasn't preloaded */
    Glta = parms.lta = LTAalloc(1, mri_in) ;

  if (!FZERO(blur_sigma))
  {
    MRI *mri_tmp, *mri_kernel ;

    mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
    mri_tmp = MRIconvolveGaussian(mri_in, NULL, mri_kernel) ;
    MRIfree(&mri_in) ; mri_in = mri_tmp ;
  }

  parms.lta->xforms[0].m_L = MatrixIdentity(4, NULL) ;
#define MAX_ITER 3
  i = 0 ;
  do
  {
    old_log_p = 
      GCAcomputeLogSampleProbability(gca, parms.gcas, mri_in, 
                                     parms.lta->xforms[0].m_L, nsamples) ;
    register_mri(mri_in, gca, &parms,i) ;
    log_p = 
      GCAcomputeLogSampleProbability(gca, parms.gcas, mri_in, 
                                     parms.lta->xforms[0].m_L, nsamples) ;

    printf("pass %d: log(p) = %2.1f (old=%2.1f)\n", i+1, log_p, old_log_p) ;
    done = (((log_p - old_log_p) / fabs(old_log_p)) < tol) ;
    if (!done)
    {
      free(parms.gcas) ;
      parms.gcas = GCAfindStableSamples(gca, &nsamples,
                                        (int)min_spacing,min_prior);
      printf("using %d sample points...\n", nsamples) ;
      parms.nsamples = nsamples ;
    }
    if (i++ >= MAX_ITER)
      break ;
  } while (!done) ;
  
  if (transformed_sample_fname)
  {
    printf("writing transformed samples to %s...\n", 
            transformed_sample_fname) ;
    GCAtransformAndWriteSamples(gca, mri_in, parms.gcas, nsamples, 
                                transformed_sample_fname, parms.lta) ;
    printf("samples written\n") ;
  }

  if (norm_fname)
  {
    int   *ordered_indices, i, label, nused, nleft_cbm, nright_cbm, 
          nleft_used, nright_used ;
    MRI   *mri_norm ;

    /* make "unknowns" the bottom of the list */
    for (nleft_cbm = nright_cbm = nused = i = 0 ; i < nsamples ; i++)
    {
      label = parms.gcas[i].label ;
      if (label == Unknown || 
          ((label != Left_Cerebral_White_Matter ) &&
           (label != Right_Cerebral_White_Matter ) &&
           (label != Left_Cerebellum_White_Matter ) &&
           (label != Right_Cerebellum_White_Matter )))
        parms.gcas[i].log_p = -100000 ;
      else
      {
        nused++ ;
        if (label == Left_Cerebellum_White_Matter )
          nleft_cbm++ ;
        else if (label == Right_Cerebellum_White_Matter)
          nright_cbm++ ;
      }
    }

    /* rank samples by log probability */
    ordered_indices = (int *)calloc(nsamples, sizeof(int)) ;
    GCArankSamples(gca, parms.gcas, nsamples, ordered_indices) ;

#if 0
    for (i = 0 ; i < nsamples ; i++)
    {
      float pct ;

      pct = 100.0f*(float)(nsamples-i)/nsamples;
      if (pct < 50)
        parms.gcas[ordered_indices[i]].label = 0 ;
    }
#else

    /* remove the least likely samples */
    printf("sorting %d (%d/%d l/r cerebellum) white matter points by "
           "likelihood\n", nused, nleft_cbm, nright_cbm) ;
    for (nleft_used = nright_used = i = 0 ; i < nused/10 ; i++)
    {
      if (parms.gcas[ordered_indices[i]].label == Left_Cerebellum_White_Matter)
        nleft_used++ ;
      else if (parms.gcas[ordered_indices[i]].label == 
               Right_Cerebellum_White_Matter)
        nright_used++ ;
      if (diag_fp)
        fprintf(diag_fp, "%d %d %d %2.1f %d\n",
                parms.gcas[ordered_indices[i]].x,
                parms.gcas[ordered_indices[i]].y,
                parms.gcas[ordered_indices[i]].z,
                parms.gcas[ordered_indices[i]].mean,
                parms.gcas[ordered_indices[i]].label) ;

    }
#define MIN_CEREBELLUM  (30)
    printf("%d/%d (l/r) cerebellar points initially in top 10%%\n", 
           nleft_used,nright_used) ;
          
    for (i = nused/10 ; i < nsamples ; i++)
    {
      if ((parms.gcas[ordered_indices[i]].label == 
           Left_Cerebellum_White_Matter) && nleft_used < MIN_CEREBELLUM)
        nleft_used++ ;
      else  if ((parms.gcas[ordered_indices[i]].label == 
                 Right_Cerebellum_White_Matter) && nright_used< MIN_CEREBELLUM)
        nright_used++ ;
      else
        parms.gcas[ordered_indices[i]].label = 0 ;
      if (diag_fp && !IS_UNKNOWN(parms.gcas[ordered_indices[i]].label))
        fprintf(diag_fp, "%d %d %d %2.1f %d\n",
                parms.gcas[ordered_indices[i]].x,
                parms.gcas[ordered_indices[i]].y,
                parms.gcas[ordered_indices[i]].z,
                parms.gcas[ordered_indices[i]].mean,
                parms.gcas[ordered_indices[i]].label) ;

    }
#endif

#if 0
    GCAtransformAndWriteSamples(gca, mri_in, parms.gcas, nsamples, 
                                transformed_sample_fname, parms.lta) ;
#endif

#if 0
    /* replace sample label with pct rank so that we can write it out easily
       for diagnostic purposes (a bit of a hack) */
    for (i = 0 ; i < nsamples ; i++)
    {
      float pct ;

      pct = 100.0f*(float)(nsamples-i)/nsamples;
      parms.gcas[ordered_indices[i]].label = pct ;
    }
    GCAtransformAndWriteSamples(gca, mri_in, parms.gcas, nsamples, 
                                norm_fname, parms.lta) ;
#endif

    if (nused/10 > 0)
      mri_norm = GCAnormalizeSamples(mri_in, gca, parms.gcas, nsamples,
                                     parms.lta) ;
    printf("writing normalized volume to %s...\n", norm_fname) ;
    MRIwrite(mri_norm, norm_fname) ;
    MRIfree(&mri_norm) ;
  }


  printf("writing output transformation to %s...\n", out_fname) ;
#if 0
  MRIvoxelXformToRasXform(mri_in, mri_in, 
                          parms.lta->xforms[0].m_L, parms.lta->xforms[0].m_L) ;
#endif
  LTAwrite(parms.lta, out_fname) ;
  if (gca)
    GCAfree(&gca) ;
  if (mri_in)
    MRIfree(&mri_in) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("registration took %d minutes and %d seconds.\n", 
          minutes, seconds) ;
  if (diag_fp)
    fclose(diag_fp) ;
  exit(0) ;
  return(0) ;
}


static int
register_mri(MRI *mri_in, GCA *gca, MORPH_PARMS *parms, int passno)
{
  MATRIX  *m_L ;

#if 0
  MRIscaleMeanIntensities(mri_in, mri_ref, mri_in);
  printf("initializing alignment using PCA...\n") ;
#endif


  if (passno == 0)
    m_L = MatrixIdentity(4, NULL) ;
  else
    m_L = MatrixCopy(parms->lta->xforms[0].m_L, NULL) ;

  find_optimal_transform(mri_in, gca, parms->gcas, parms->nsamples,m_L,passno,
                         parms->write_iterations);

  if (!parms->lta)
    parms->lta = LTAalloc(1, NULL) ;
  MatrixCopy(m_L, parms->lta->xforms[0].m_L) ;
  if (Gdiag & DIAG_SHOW)
  {
    printf("global search transform:\n") ;
    MatrixPrint(stdout, m_L) ;
  }

  printf("computing MAP estimate of linear transform...\n") ;

  parms->mri_in = mri_in ;  /* for diagnostics */
  MRIemAlign(mri_in, gca, parms, m_L) ;

  printf("final transform:\n") ;
  MatrixPrint(stdout, parms->lta->xforms[0].m_L) ;
  printf("\n") ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRI *mri_aligned ;

    mri_aligned = 
      MRIapplyRASlinearTransform(mri_in, NULL, parms->lta->xforms[0].m_L) ;
    MRIwriteImageViews(mri_aligned, "after_alignment", 400) ;
    MRIfree(&mri_aligned) ;
  }

  return(NO_ERROR) ;
}

#define DEFAULT_MAX_STEPS 5
static double MAX_ANGLES = DEFAULT_MAX_STEPS ;
#define MAX_ANGLE       RADIANS(30)
#define MIN_ANGLE       RADIANS(2)

#define MAX_SCALE       1.3
#define MIN_SCALE       0.7

static int max_angles = DEFAULT_MAX_STEPS ;


static int MAX_TRANS_STEPS = DEFAULT_MAX_STEPS ;
static double MAX_TRANS = 50 ;

static MATRIX *
find_optimal_transform(MRI *mri, GCA *gca, GCA_SAMPLE *gcas, int nsamples,
                       MATRIX *m_L, int passno, int write_iterations)
{
  MATRIX   *m_origin ;
  MRI      *mri_gca  ;
  double   in_means[3], gca_means[3], dx, dy, dz, max_log_p, old_max,
           max_angle, angle_steps, min_scale, max_scale, scale_steps, scale,
           delta, mean, wm_mean ;
  int      niter, mri_peak ;
  HISTOGRAM *h_mri, *h_smooth ;
  MRI_REGION box ;


  if (!passno && Gdiag & DIAG_WRITE && write_iterations > 0)
  {
    char fname[STRLEN] ;

    sprintf(fname, "%s_before_intensity.mgh", parms.base_name) ;
    printf("writing snapshot to %s...\n", fname) ;
    MRIwrite(mri, fname) ;
    MRIwriteImageViews(mri, "before_intensity", 400) ;
  }

#define MIN_BIN  50
#define OFFSET_SIZE  25
#define BOX_SIZE   60   /* mm */
#define HALF_BOX   (BOX_SIZE/2)

  if (passno == 0)
  {
    float   x0, y0, z0 ;
    MRIfindApproximateSkullBoundingBox(mri, 50, &box) ;
    x0 = box.x+box.dx/2 ; y0 = box.y+box.dy/2 ; z0 = box.z+box.dz/2 ;
    printf("using (%.0f, %.0f, %.0f) as brain centroid...\n",x0, y0, z0) ;
    box.x = x0 - HALF_BOX*mri->xsize ; box.dx = BOX_SIZE*mri->xsize ;
    box.y = y0 - HALF_BOX*mri->ysize ; box.dy = BOX_SIZE*mri->ysize ;
    box.x = z0 - HALF_BOX*mri->zsize ; box.dz = BOX_SIZE*mri->zsize ;
    wm_mean = GCAlabelMean(gca, Left_Cerebral_White_Matter) ;
    wm_mean = (wm_mean + GCAlabelMean(gca, Right_Cerebral_White_Matter))/2.0 ;
    printf("mean wm in atlas = %2.0f, using box (%d,%d,%d) --> (%d, %d,%d) "
           "to find MRI wm\n", wm_mean, box.x, box.y, box.z, 
           box.x+box.dx-1,box.y+box.dy-1, box.z+box.dz-1) ;
    
    h_mri = MRIhistogramRegion(mri, 0, NULL, &box) ; 
    HISTOclearBins(h_mri, h_mri, 0, MIN_BIN) ;
    HISTOplot(h_mri, "mri.histo") ;
    mri_peak = HISTOfindLastPeak(h_mri, HISTO_WINDOW_SIZE,MIN_HISTO_PCT);
    printf("before smoothing, mri peak at %d\n", mri_peak) ;
    h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
    HISTOplot(h_smooth, "mri_smooth.histo") ;
    mri_peak = HISTOfindLastPeak(h_smooth, HISTO_WINDOW_SIZE,MIN_HISTO_PCT);
    printf("after smoothing, mri peak at %d, scaling input intensities "
           "by %2.1f\n", mri_peak, wm_mean/mri_peak) ;
    HISTOfree(&h_smooth) ; HISTOfree(&h_mri) ;
    MRIscalarMul(mri, mri, wm_mean/mri_peak) ;
  }

  
  /* first align centroids */
  mri_gca = MRIclone(mri, NULL) ;
  GCAmri(gca, mri_gca) ;

  if (!passno && Gdiag & DIAG_WRITE && write_iterations > 0)
  {
    char fname[STRLEN] ;

    Glta->xforms[0].m_L = m_L ;
    sprintf(fname, "%s_after_intensity.mgh", parms.base_name) ;
    printf("writing snapshot to %s...\n", fname) ;
    MRIwrite(mri, fname) ;
    MRIwriteImageViews(mri, "before_centering", 400) ;
    sprintf(fname, "%s_before_centering.mgh", parms.base_name) ;
    GCAtransformAndWriteSamples(gca, mri, gcas, nsamples, 
                                fname, Glta) ;
  }
  if (gca_mean_fname)
  {
    printf("writing gca volume to %s...\n", gca_mean_fname) ;
    MRIwrite(mri_gca, gca_mean_fname) ;
    printf("done\n") ;
  }

  max_log_p = GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples) ;
  printf("initial log_p = %2.1f\n", max_log_p) ;

  MRIcenterOfMass(mri, in_means, 0) ;
  MRIcenterOfMass(mri_gca, gca_means, 0) ;
  printf("input centroid (%2.1f, %2.1f, %2.1f), "
          "gca centroid (%2.1f, %2.1f, %2.1f)\n",
          in_means[0], in_means[1], in_means[2],
          gca_means[0], gca_means[1], gca_means[2]) ;



  MRIfree(&mri_gca) ;


  /* now apply translation to take in centroid to ref centroid */
  dx = gca_means[0] - in_means[0] ; dy = gca_means[1] - in_means[1] ;
  dz = gca_means[2] - in_means[2] ;
  if (passno == 0)
  {
    *MATRIX_RELT(m_L, 1, 4) = dx ; *MATRIX_RELT(m_L, 2, 4) = dy ;
    *MATRIX_RELT(m_L, 3, 4) = dz ;
  }

  m_origin = MatrixIdentity(4, NULL) ;
  *MATRIX_RELT(m_origin, 1, 4) = gca_means[0]*(float)center ; 
  *MATRIX_RELT(m_origin, 2, 4) = gca_means[1]*(float)center ;
  *MATRIX_RELT(m_origin, 3, 4) = gca_means[2]*(float)center ; 
  *MATRIX_RELT(m_origin, 4, 4) = 1 ;

  max_log_p = GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples) ;
  printf("initial translation: (%2.1f, %2.1f, %2.1f): log p = %2.1f\n",
          dx,dy,dz, max_log_p) ;

  max_angle = MAX_ANGLE ; angle_steps = max_angles ;
  max_scale = MAX_SCALE ; min_scale = MIN_SCALE ; scale_steps = max_angles ;

  if (write_iterations != 0)
  {
    char fname[STRLEN] ;
    MRI *mri_aligned ;

    mri_aligned = MRIlinearTransform(mri, NULL, m_L) ;
    sprintf(fname, "%s_after_centering", parms.base_name) ;
    MRIwriteImageViews(mri_aligned, fname, 400) ;
    sprintf(fname, "%s_after_centering.mgh", parms.base_name) ;
    printf("writing image after centering to %s...\n", fname) ;
#if 0
    MRIwrite(mri_aligned, fname) ;
#else
    Glta->xforms[0].m_L = m_L ;
    GCAtransformAndWriteSamples(gca, mri, gcas, nsamples, 
                                fname, Glta) ;
#endif
    MRIfree(&mri_aligned) ;
    
  }
#define MIN_ITER 1
  niter = 0 ; scale = 1.0 ;
  do
  {
    old_max = max_log_p ;
    max_log_p = find_optimal_translation(gca, gcas, mri, nsamples, m_L,
                                         -scale*MAX_TRANS, scale*MAX_TRANS, 
                                         MAX_TRANS_STEPS, 2) ;

    if (!translation_only)
      max_log_p = find_optimal_rotation(gca, gcas, mri, nsamples, m_L,m_origin,
                                        -scale*max_angle, scale*max_angle, 
                                        angle_steps, 3) ;

    if (!noscale && !translation_only)
    {
      max_log_p = find_optimal_scaling(gca, gcas, mri, nsamples, m_L, m_origin,
                                       min_scale, max_scale, scale_steps, 3) ;
      max_log_p = find_optimal_scaling_and_rotation(gca, gcas, mri, nsamples, 
                                                    m_L, m_origin,
                                                    -RADIANS(2*scale),
                                                    RADIANS(2*scale),
                                                  1-.025*scale, 1+.025*scale, 
                                                  max_angles/3,max_angles/3,3);
    }
    if (write_iterations != 0)
    {
      char fname[STRLEN] ;
      MRI *mri_aligned ;
      
      mri_aligned = MRIlinearTransform(mri, NULL, m_L) ;
      sprintf(fname, "%s%03d", parms.base_name, niter+1) ;
      MRIwriteImageViews(mri_aligned, fname, 400) ;
      sprintf(fname, "%s%03d.mgh", parms.base_name, niter+1) ;
      printf("writing image after centering to %s...\n", fname) ;
#if 0
      MRIwrite(mri_aligned, fname) ;
#else
      Glta->xforms[0].m_L = m_L ;
      GCAtransformAndWriteSamples(gca, mri, gcas, nsamples, 
                                  fname, Glta) ;
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

static double
find_optimal_rotation(GCA *gca, GCA_SAMPLE *gcas, MRI *mri, int nsamples, 
                   MATRIX *m_L, MATRIX *m_origin,
                   float min_angle, float max_angle, float angle_steps,
                   int nreductions)
{
  MATRIX   *m_rot, *m_x_rot, *m_y_rot, *m_z_rot, *m_tmp,*m_L_tmp,*m_origin_inv,
           *m_tmp2;
  double   x_angle, y_angle, z_angle, x_max, y_max, z_max, delta, 
           log_p, max_log_p, mean_angle ;
  int      i ;

  m_origin_inv = MatrixCopy(m_origin, NULL) ;
  *MATRIX_RELT(m_origin_inv, 1, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 2, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 3, 4) *= -1 ;
  m_L_tmp = m_x_rot = m_y_rot = m_z_rot = m_rot = m_tmp = m_tmp2 = NULL ;
  x_max = y_max = z_max = 0.0 ;
  max_log_p = GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples) ;
  for (i = 0 ; i <= nreductions ; i++)
  {
    delta = (max_angle-min_angle) / angle_steps ;
    if (FZERO(delta))
      return(max_log_p) ;

    if (Gdiag & DIAG_SHOW)
      printf(
              "scanning %2.2f degree nbhd (%2.1f) ",
              (float)DEGREES(max_angle), (float)DEGREES(delta)) ;

    for (x_angle = min_angle ; x_angle <= max_angle ; x_angle += delta)
    {
      m_x_rot = MatrixReallocRotation(4, x_angle, X_ROTATION, m_x_rot) ;
      for (y_angle = min_angle ; y_angle <= max_angle ; y_angle += delta)
      {
        m_y_rot = MatrixReallocRotation(4, y_angle, Y_ROTATION, m_y_rot) ;
        m_tmp = MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
        for (z_angle= min_angle ; z_angle <= max_angle ; z_angle += delta)
        {
          if (nint(DEGREES(x_angle)) == -9 && nint(DEGREES(y_angle)) == -5 &&
              nint(DEGREES(z_angle)) == -7)
            DiagBreak() ;
          m_z_rot = MatrixReallocRotation(4, z_angle, Z_ROTATION, m_z_rot) ;
          m_rot = MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
          m_tmp2 = MatrixMultiply(m_rot, m_origin_inv, m_tmp2) ;
          MatrixMultiply(m_origin, m_tmp2, m_rot) ;

          m_L_tmp = MatrixMultiply(m_rot, m_L, m_L_tmp) ;
          log_p = 
            GCAcomputeLogSampleProbability(gca, gcas, mri, m_L_tmp,nsamples) ;
          if (log_p > max_log_p)
          {
            max_log_p = log_p ;
            x_max = x_angle ; y_max = y_angle ; z_max = z_angle ;
#if 0
            printf("new max p %2.1f found at (%2.1f, %2.1f, %2.1f)\n",
                    max_log_p, (float)DEGREES(x_angle), 
                    (float)DEGREES(y_angle), (float)DEGREES(z_angle)) ;
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
    min_angle -= delta/2 ; max_angle += delta/2 ;
#endif
  }

  printf("\n") ;

  MatrixFree(&m_x_rot) ; MatrixFree(&m_y_rot) ; MatrixFree(&m_z_rot) ;
  MatrixFree(&m_rot) ;   MatrixFree(&m_tmp) ; MatrixFree(&m_origin_inv) ;
  MatrixFree(&m_tmp2) ;
  return(max_log_p) ;
}
static double
find_optimal_scaling(GCA *gca, GCA_SAMPLE *gcas, MRI *mri, int nsamples, 
                   MATRIX *m_L, MATRIX *m_origin,
                   float min_scale, float max_scale, float scale_steps,
                   int nreductions)
{
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
  max_log_p = GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples) ;
  for (i = 0 ; i <= nreductions ; i++)
  {
    delta = (max_scale-min_scale) / scale_steps ;
    if (FZERO(delta))
      return(max_log_p) ;
    if (Gdiag & DIAG_SHOW)
    {
      printf("scanning scales %2.3f->%2.3f (%2.3f) ",
              min_scale,max_scale, delta) ;
      fflush(stdout) ;
    }
    for (x_scale = min_scale ; x_scale <= max_scale ; x_scale += delta)
    {
      *MATRIX_RELT(m_scale, 1, 1) = x_scale ;
      for (y_scale = min_scale ; y_scale <= max_scale ; y_scale += delta)
      {
        *MATRIX_RELT(m_scale, 2, 2) = y_scale ;
        for (z_scale= min_scale ; z_scale <= max_scale ; z_scale += delta)
        {
          *MATRIX_RELT(m_scale, 3, 3) = z_scale ;

          /* reset translation values */
          *MATRIX_RELT(m_scale, 1, 4) = 
            *MATRIX_RELT(m_scale, 2, 4) = *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
          if (nint((x_scale)) == -9 && nint((y_scale)) == -5 &&
              nint((z_scale)) == -7)
            DiagBreak() ;

          m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
          MatrixMultiply(m_origin, m_tmp, m_scale) ;

          m_L_tmp = MatrixMultiply(m_scale, m_L, m_L_tmp) ;
          log_p = 
            GCAcomputeLogSampleProbability(gca, gcas, mri, m_L_tmp,nsamples) ;
          if (log_p > max_log_p)
          {
            max_log_p = log_p ;
            x_max = x_scale ; y_max = y_scale ; z_max = z_scale ;
#if 0
            printf("new max p %2.1f found at (%2.3f, %2.3f, %2.3f)\n",
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

  MatrixFree(&m_scale) ;   MatrixFree(&m_tmp) ; MatrixFree(&m_origin_inv) ;
  return(max_log_p) ;
}

static double
find_optimal_translation(GCA *gca, GCA_SAMPLE *gcas, MRI *mri, int nsamples, 
                         MATRIX *m_L, float min_trans, float max_trans, 
                         float trans_steps, int nreductions)
{
  MATRIX   *m_trans, *m_L_tmp ;
  double   x_trans, y_trans, z_trans, x_max, y_max, z_max, delta, 
           log_p, max_log_p, mean_trans ;
  int      i ;

  delta = (max_trans-min_trans) / trans_steps ;
  m_L_tmp = NULL ;
  m_trans = MatrixIdentity(4, NULL) ;
  x_max = y_max = z_max = 0.0 ;
  max_log_p = GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples) ;

  for (i = 0 ; i <= nreductions ; i++)
  {
    delta = (max_trans-min_trans) / trans_steps ;
    if (FZERO(delta))
      return(max_log_p) ;
    if (Gdiag & DIAG_SHOW)
    {
      printf(
              "scanning translations %2.2f->%2.2f (%2.1f) ",
              min_trans,max_trans, delta) ;
      fflush(stdout) ;
    }
    for (x_trans = min_trans ; x_trans <= max_trans ; x_trans += delta)
    {
      *MATRIX_RELT(m_trans, 1, 4) = x_trans ;
      for (y_trans = min_trans ; y_trans <= max_trans ; y_trans += delta)
      {
        *MATRIX_RELT(m_trans, 2, 4) = y_trans ;
        for (z_trans= min_trans ; z_trans <= max_trans ; z_trans += delta)
        {
          *MATRIX_RELT(m_trans, 3, 4) = z_trans ;
          if (nint((x_trans)) == -9 && nint((y_trans)) == -5 &&
              nint((z_trans)) == -7)
            DiagBreak() ;

          m_L_tmp = MatrixMultiply(m_trans, m_L, m_L_tmp) ;
          log_p = 
            GCAcomputeLogSampleProbability(gca, gcas, mri, m_L_tmp,nsamples) ;
          if (log_p > max_log_p)
          {
            max_log_p = log_p ;
            x_max = x_trans ; y_max = y_trans ; z_max = z_trans ;
#if 0
            printf("new max p %2.1f found at (%2.1f, %2.1f, %2.1f)\n",
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
    x_max = y_max = z_max = 0.0 ;  /* we've translated transform by old maxs */

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
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  StrUpper(option) ;
  if (!strcmp(option, "DIST") || !strcmp(option, "DISTANCE"))
  {
    parms.l_dist = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_dist = %2.2f\n", parms.l_dist) ;
  }
  else if (!strcmp(option, "SAMPLES"))
  {
    sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing control points to %s...\n", sample_fname) ;
  }
  else if (!strcmp(option, "DIAG"))
  {
    diag_fp = fopen(argv[2], "w") ;
    if (!diag_fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open diag file %s for writing",
                Progname, argv[2]) ;
    printf("opening diag file %s for writing\n", argv[2]) ;
    nargs = 1 ;
  }
  else if (!strcmp(option, "TR"))
  {
    TR = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TR=%2.1f msec\n", TR) ;
  }
  else if (!strcmp(option, "TE"))
  {
    TE = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TE=%2.1f msec\n", TE) ;
  }
  else if (!strcmp(option, "ALPHA"))
  {
    nargs = 1 ;
    alpha = RADIANS(atof(argv[2])) ;
    printf("using alpha=%2.0f degrees\n", DEGREES(alpha)) ;
  }
  else if (!strcmp(option, "FSAMPLES") || !strcmp(option, "ISAMPLES"))
  {
    transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing transformed control points to %s...\n", 
            transformed_sample_fname) ;
  }
  else if (!strcmp(option, "CONTRAST"))
  {
    use_contrast = 1 ;
    printf("using contrast to find labels...\n") ;
  }
  else if (!strcmp(option, "RENORM"))
  {
    renormalization_fname = argv[2] ;
    nargs = 1 ;
    printf("renormalizing using predicted intensity values in %s...\n",
           renormalization_fname) ;
  }
  else if (!strcmp(option, "FLASH"))
  {
    tissue_parms_fname = argv[2] ;
    nargs = 1 ;
    printf("using FLASH forward model and tissue parms in %s to predict"
           " intensity values...\n", tissue_parms_fname) ;
  }
  else if (!strcmp(option, "TRANSONLY"))
  {
    translation_only = 1 ;
    printf("only computing translation parameters...\n") ;
  }
  else if (!strcmp(option, "WRITE_MEAN"))
  {
    gca_mean_fname = argv[2] ;
    nargs = 1 ;
    printf("writing gca means to %s...\n", gca_mean_fname) ;
  }
  else if (!strcmp(option, "PRIOR"))
  {
    min_prior = atof(argv[2]) ;
    nargs = 1 ;
    printf("using prior threshold %2.2f\n", min_prior) ;
  }
  else if (!strcmp(option, "SPACING"))
  {
    min_spacing = atof(argv[2]) ;
    nargs = 1 ;
    printf("using min GCA spacing %2.0f...\n", min_spacing) ;
  }
  else if (!stricmp(option, "NOVAR"))
  {
    novar = 1 ;
    printf("not using variance estimates\n") ;
  }
  else if (!strcmp(option, "DT"))
  {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    printf("dt = %2.2e\n", parms.dt) ;
  }
  else if (!strcmp(option, "TOL"))
  {
    tol = parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("tol = %2.2e\n", parms.tol) ;
  }
  else if (!strcmp(option, "CENTER"))
  {
    center = 1 ;
    printf("using GCA centroid as origin of transform\n") ;
  }
  else if (!strcmp(option, "NOSCALE"))
  {
    noscale = 1 ;
    printf("disabling scaling...\n") ;
  }
  else if (!strcmp(option, "NUM"))
  {
    num_xforms = atoi(argv[2]) ;
    nargs = 1 ;
    printf("finding a total of %d linear transforms\n", num_xforms) ;
  }
  else if (!strcmp(option, "AREA"))
  {
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_area = %2.2f\n", parms.l_area) ;
  }
  else if (!strcmp(option, "NLAREA"))
  {
    parms.l_nlarea = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_nlarea = %2.2f\n", parms.l_nlarea) ;
  }
  else if (!strcmp(option, "LEVELS"))
  {
    parms.levels = atoi(argv[2]) ;
    nargs = 1 ;
    printf("levels = %d\n", parms.levels) ;
  }
  else if (!strcmp(option, "INTENSITY") || !strcmp(option, "CORR"))
  {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_intensity = %2.2f\n", parms.l_intensity) ;
  }
  else if (!stricmp(option, "reduce"))
  {
    nreductions = atoi(argv[2]) ;
    nargs = 1 ;
    printf("reducing input images %d times before aligning...\n",
            nreductions) ;
  }
  else if (!stricmp(option, "nsamples"))
  {
    nsamples = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d samples of GCA...\n", nsamples) ;
  }
  else if (!stricmp(option, "norm"))
  {
    norm_fname = argv[2] ;
    nargs = 1 ;
    printf("intensity normalizing and writing to %s...\n",norm_fname);
  }
  else if (!stricmp(option, "steps"))
  {
    max_angles = atoi(argv[2]) ;
    nargs = 1 ;
    printf("taking %d angular steps...\n", max_angles) ;
  }
  else switch (*option)
  {
  case 'D':
    tx = atof(argv[2]) ; ty = atof(argv[3]) ; tz = atof(argv[4]) ;
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
    printf("examining %2.0f different trans/rot/scale values...\n",MAX_ANGLES);
#endif
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
static double 
find_optimal_scaling_and_rotation(GCA *gca, GCA_SAMPLE *gcas, 
                                  MRI *mri, 
                                  int nsamples, MATRIX *m_L, MATRIX *m_origin,
                                  float min_angle, float max_angle,
                                  float min_scale, float max_scale,
                                  float angle_steps, float scale_steps,
                                  int nreductions)
{
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
  max_log_p = GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples) ;
  for (i = 0 ; i <= nreductions ; i++)
  {
    delta_scale = (max_scale-min_scale) / scale_steps ;
    delta = (max_angle-min_angle) / angle_steps ;
    if (Gdiag & DIAG_SHOW)
    {
      printf("scanning %2.2f degree nbhd (%2.1f) and "
             "scale %2.3f->%2.3f (%2.3f)\n",
              (float)DEGREES(max_angle), (float)DEGREES(delta),
             min_scale,max_scale, delta_scale);
      fflush(stdout) ;
    }

    for (x_scale = min_scale ; x_scale <= max_scale ; x_scale += delta_scale)
    {
      /*      printf("x_scale = %2.3f\n", x_scale) ;*/
      *MATRIX_RELT(m_scale, 1, 1) = x_scale ;
      for (y_scale = min_scale ; y_scale <= max_scale ; y_scale += delta_scale)
      {
        *MATRIX_RELT(m_scale, 2, 2) = y_scale ;
        for (z_scale= min_scale ; z_scale <= max_scale; z_scale += delta_scale)
        {
          *MATRIX_RELT(m_scale, 3, 3) = z_scale ;

          /* reset translation values */
          *MATRIX_RELT(m_scale, 1, 4) = 
            *MATRIX_RELT(m_scale, 2, 4) = *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
          if (nint((x_scale)) == -9 && nint((y_scale)) == -5 &&
              nint((z_scale)) == -7)
            DiagBreak() ;

          m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
          MatrixMultiply(m_origin, m_tmp, m_scale) ;

          for (x_angle = min_angle ; x_angle <= max_angle ; x_angle += delta)
          {
            m_x_rot = MatrixReallocRotation(4, x_angle, X_ROTATION, m_x_rot) ;
            for (y_angle = min_angle ; y_angle <= max_angle ; y_angle += delta)
            {
              m_y_rot = MatrixReallocRotation(4, y_angle, Y_ROTATION, m_y_rot);
              m_tmp = MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
              for (z_angle= min_angle; z_angle <= max_angle; z_angle += delta)
              {
                m_z_rot = MatrixReallocRotation(4, z_angle,Z_ROTATION,m_z_rot);
                m_rot = MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
                m_tmp2 = MatrixMultiply(m_rot, m_origin_inv, m_tmp2) ;
                MatrixMultiply(m_origin, m_tmp2, m_rot) ;

                
                m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
                m_L_tmp = MatrixMultiply(m_tmp2, m_L, m_L_tmp) ;
                log_p = 
                  GCAcomputeLogSampleProbability(gca, gcas, mri, m_L_tmp,
                                                 nsamples) ;
                if (log_p > max_log_p)
                {
                  max_log_p = log_p ;
                  x_max_scale = x_scale ; y_max_scale = y_scale ; 
                  z_max_scale = z_scale ;
                  x_max_rot = x_angle ; y_max_rot = y_angle ; 
                  z_max_rot = z_angle ;
                }
              }
            }
          }
        }
      }
      
    }

    if (Gdiag & DIAG_SHOW)
    {
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
    min_angle -= delta/2 ; max_angle += delta/2 ;
#endif
  }

  printf("\n") ;

  MatrixFree(&m_x_rot) ; MatrixFree(&m_y_rot) ; MatrixFree(&m_z_rot) ;
  MatrixFree(&m_rot) ;   MatrixFree(&m_tmp) ; MatrixFree(&m_origin_inv) ;
  MatrixFree(&m_tmp2) ;
  return(max_log_p) ;
}
