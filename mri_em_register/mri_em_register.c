

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

char         *Progname ;
static MORPH_PARMS  parms ;

static char *sample_fname = NULL ;
static char *transformed_sample_fname = NULL ;

#define MIN_SPACING   16.0
static float min_spacing = MIN_SPACING ;

static int use_contrast = 0 ;
static float min_prior = MIN_PRIOR ;
static double tol = 0.01 ;
static double tx = 0.0 ;
static double ty = 0.0 ;
static double tz = 0.0 ;
static double rzrot = 0.0 ;
static double rxrot = 0.0 ;
static double ryrot = 0.0 ;

static int get_option(int argc, char *argv[]) ;
static int register_mri(MRI *mri_in, GCA *gca, MP *parms) ;

static int center = 1 ;
static int nreductions = 1 ;
static int noscale = 0 ;
static int num_xforms = 1 ;
static int transform_loaded = 0 ;
static char *gca_mean_fname = NULL ;

static MATRIX *find_optimal_transform(MRI *mri_in, GCA *gca, GCA_SAMPLE *gcas,
                                      int nsamples) ;
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
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  parms.l_intensity = 1.0f ;
  parms.niterations = 100 ;
  parms.levels = -1 ;   /* use default */
  parms.dt = 1e-6 ;  /* was 5e-6 */
  parms.tol = INTEGRATION_TOL*5 ;

  parms.max_levels = 0 ;
  parms.dt = 5e-6 ;  /* was 5e-6 */
  parms.tol = 1e-3 ;
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
  fprintf(stderr, "logging results to %s.log\n", parms.base_name) ;

  TimerStart(&start) ;
  fprintf(stderr, "reading '%s'...\n", gca_fname) ;
  fflush(stderr) ;
  gca = GCAread(gca_fname) ;
  if (gca == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open GCA %s.\n",
              Progname, gca_fname) ;


#if 0
  if (gca->spacing < min_spacing)
    fprintf(stderr, 
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




  fprintf(stderr, "reading '%s'...\n", in_fname) ;
  fflush(stderr) ;
  mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
              Progname, in_fname) ;

  if (sample_fname)
  {
    fprintf(stderr, "writing samples to %s...\n", sample_fname) ;
    GCAwriteSamples(gca, mri_in, parms.gcas, nsamples, sample_fname) ;
    fprintf(stderr, "samples written\n") ;
  }
#if 0
  if (gca_reduced != gca)
    GCAfree(&gca_reduced) ;
#endif

  if (!FZERO(tx) || !FZERO(ty) || !FZERO(tz))
  {
    MRI *mri_tmp ;
    
    fprintf(stderr, "translating second volume by (%2.1f, %2.1f, %2.1f)\n",
            tx, ty, tz) ;
    mri_tmp = MRItranslate(mri_in, NULL, tx, ty, tz) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
#if 1
  if (!FZERO(rzrot))
  {
    MRI *mri_tmp ;
    
    fprintf(stderr, 
            "rotating second volume by %2.1f degrees around Z axis\n",
            (float)DEGREES(rzrot)) ;
    mri_tmp = MRIrotateZ_I(mri_in, NULL, rzrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  if (!FZERO(rxrot))
  {
    MRI *mri_tmp ;
    
    fprintf(stderr, 
            "rotating second volume by %2.1f degrees around X axis\n",
            (float)DEGREES(rxrot)) ;
    mri_tmp = MRIrotateX_I(mri_in, NULL, rxrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  if (!FZERO(ryrot))
  {
    MRI *mri_tmp ;
    
    fprintf(stderr, 
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
      fprintf(stderr, 
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
    parms.lta = LTAalloc(1, mri_in) ;

  if (!FZERO(blur_sigma))
  {
    MRI *mri_tmp, *mri_kernel ;

    mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
    mri_tmp = MRIconvolveGaussian(mri_in, NULL, mri_kernel) ;
    MRIfree(&mri_in) ; mri_in = mri_tmp ;
  }
  register_mri(mri_in, gca, &parms) ;
  
  if (transformed_sample_fname)
  {
    fprintf(stderr, "writing transformed samples to %s...\n", 
            transformed_sample_fname) ;
    GCAtransformAndWriteSamples(gca, mri_in, parms.gcas, nsamples, 
                                transformed_sample_fname, parms.lta) ;
    fprintf(stderr, "samples written\n") ;
  }

  fprintf(stderr, "writing output transformation to %s...\n", out_fname) ;
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
  fprintf(stderr, "registration took %d minutes and %d seconds.\n", 
          minutes, seconds) ;
  exit(0) ;
  return(0) ;
}


static int
register_mri(MRI *mri_in, GCA *gca, MORPH_PARMS *parms)
{
  MATRIX  *m_L ;

#if 0
  MRIscaleMeanIntensities(mri_in, mri_ref, mri_in);
  fprintf(stderr, "initializing alignment using PCA...\n") ;
#endif
  if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
  {
    MRIwriteImageViews(mri_in, "before_pca", 400) ;
  }

#if 0
  m_L = MatrixIdentity(4, NULL) ;
#else
  m_L = find_optimal_transform(mri_in, gca, parms->gcas, parms->nsamples) ;
#endif

  if (!parms->lta)
    parms->lta = LTAalloc(1, NULL) ;
  MatrixCopy(m_L, parms->lta->xforms[0].m_L) ;
  if (Gdiag & DIAG_SHOW)
  {
    printf("global search transform:\n") ;
    MatrixPrint(stdout, m_L) ;
  }

#if 0
  fprintf(stderr, "computing MAP estimate of linear transform...\n") ;

  parms->mri_in = mri_in ;  /* for diagnostics */
  MRIemAlign(mri_in, gca, parms, m_L) ;

  fprintf(stderr, "final transform:\n") ;
  MatrixPrint(stderr, parms->lta->xforms[0].m_L) ;
  fprintf(stderr, "\n") ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRI *mri_aligned ;

    mri_aligned = 
      MRIapplyRASlinearTransform(mri_in, NULL, parms->lta->xforms[0].m_L) ;
    MRIwriteImageViews(mri_aligned, "after_alignment", 400) ;
    MRIfree(&mri_aligned) ;
  }

#endif

  return(NO_ERROR) ;
}

#define MAX_ANGLES      15
#define MAX_ANGLE       RADIANS(30)
#define MIN_ANGLE       RADIANS(2)

#define MAX_SCALE       1.3
#define MIN_SCALE       0.7

static int max_angles = MAX_ANGLES ;


#define MAX_TRANS       15

static MATRIX *
find_optimal_transform(MRI *mri, GCA *gca, GCA_SAMPLE *gcas, int nsamples)
{
  MATRIX   *m_L, *m_origin ;
  MRI      *mri_gca ;
  double   in_means[3], gca_means[3], dx, dy, dz, max_log_p, old_max,
           max_angle, angle_steps, min_scale, max_scale, scale_steps, scale,
           delta, mean ;

  m_L = MatrixIdentity(4, NULL) ;

  /* first align centroids */
  mri_gca = MRIclone(mri, NULL) ;
  GCAmri(gca, mri_gca) ;

  if (gca_mean_fname)
  {
    fprintf(stderr, "writing gca volume to %s...\n", gca_mean_fname) ;
    MRIwrite(mri_gca, gca_mean_fname) ;
    fprintf(stderr, "done\n") ;
  }

  MRIcenterOfMass(mri, in_means, 0) ;
  MRIcenterOfMass(mri_gca, gca_means, 0) ;
  fprintf(stderr, "input centroid (%2.1f, %2.1f, %2.1f), "
          "gca centroid (%2.1f, %2.1f, %2.1f)\n",
          in_means[0], in_means[1], in_means[2],
          gca_means[0], gca_means[1], gca_means[2]) ;

  MRIfree(&mri_gca) ;

  /* now apply translation to take in centroid to ref centroid */
  dx = gca_means[0] - in_means[0] ; dy = gca_means[1] - in_means[1] ;
  dz = gca_means[2] - in_means[2] ;
  *MATRIX_RELT(m_L, 1, 4) = dx ; *MATRIX_RELT(m_L, 2, 4) = dy ;
  *MATRIX_RELT(m_L, 3, 4) = dz ;

  m_origin = MatrixIdentity(4, NULL) ;
  *MATRIX_RELT(m_origin, 1, 4) = gca_means[0]*(float)center ; 
  *MATRIX_RELT(m_origin, 2, 4) = gca_means[1]*(float)center ;
  *MATRIX_RELT(m_origin, 3, 4) = gca_means[2]*(float)center ; 
  *MATRIX_RELT(m_origin, 4, 4) = 1 ;

  max_log_p = GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples) ;
  fprintf(stderr,"initial translation: (%2.1f, %2.1f, %2.1f): log p = %2.1f\n",
          dx,dy,dz, max_log_p) ;

  max_angle = MAX_ANGLE ; angle_steps = max_angles ;
  max_scale = MAX_SCALE ; min_scale = MIN_SCALE ; scale_steps = max_angles ;

  scale = 1.0 ;
  do
  {
    old_max = max_log_p ;
    max_log_p = find_optimal_translation(gca, gcas, mri, nsamples, m_L,
                                         -scale*MAX_TRANS, scale*MAX_TRANS, 
                                         MAX_TRANS, 2) ;

    max_log_p = find_optimal_rotation(gca, gcas, mri, nsamples, m_L, m_origin,
                                      -scale*max_angle, 
                                      scale*max_angle, angle_steps, 3) ;

    if (!noscale)
      max_log_p = find_optimal_scaling(gca, gcas, mri, nsamples, m_L, m_origin,
                                       min_scale, max_scale, scale_steps, 3) ;
    
    fprintf(stderr, "scale %2.3f: max=%2.1f, old_max =%2.1f (thresh=%2.1f)\n",
            scale,max_log_p, old_max, old_max+fabs(tol*old_max)) ;

    /* search a finer nbhd (if do-while continues) */
    scale *= 0.75 ;
    mean = (max_scale + min_scale)/2 ;
    delta = (max_scale - min_scale)/2 ;
    max_scale = mean + delta*scale ;
    min_scale = mean - delta*scale ;

  } while (max_log_p > old_max+fabs(tol*old_max)) ;

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
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, 
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
            fprintf(stderr, "new max p %2.1f found at (%2.1f, %2.1f, %2.1f)\n",
                    max_log_p, (float)DEGREES(x_angle), 
                    (float)DEGREES(y_angle), (float)DEGREES(z_angle)) ;
#endif
          }
        }
      }
      
    }

    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, 
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

  fprintf(stderr, "\n") ;

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
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "scanning scales %2.3f->%2.3f (%2.3f) ",
              min_scale,max_scale, delta) ;
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
            fprintf(stderr, "new max p %2.1f found at (%2.3f, %2.3f, %2.3f)\n",
                    max_log_p, x_scale, y_scale, z_scale) ;
#endif
          }
        }
      }
      
    }

    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "max log p = %2.1f @ (%2.3f, %2.3f, %2.3f)\n", 
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

  fprintf(stderr, "\n") ;

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

  m_L_tmp = NULL ;
  m_trans = MatrixIdentity(4, NULL) ;
  x_max = y_max = z_max = 0.0 ;
  max_log_p = GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples) ;
  *MATRIX_RELT(m_trans, 4, 4) = 1.0 ;
  for (i = 0 ; i <= nreductions ; i++)
  {
    delta = (max_trans-min_trans) / trans_steps ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, 
              "scanning translations %2.2f->%2.2f (%2.1f) ",
              min_trans,max_trans, delta) ;
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
            fprintf(stderr, "new max p %2.1f found at (%2.1f, %2.1f, %2.1f)\n",
                    max_log_p, x_trans, y_trans, z_trans) ;
#endif
          }
        }
      }
      
    }

    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, 
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

  fprintf(stderr, "\n") ;

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
    fprintf(stderr, "l_dist = %2.2f\n", parms.l_dist) ;
  }
  else if (!strcmp(option, "SAMPLES"))
  {
    sample_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "writing control points to %s...\n", sample_fname) ;
  }
  else if (!strcmp(option, "ISAMPLES"))
  {
    transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "writing transformed control points to %s...\n", 
            transformed_sample_fname) ;
  }
  else if (!strcmp(option, "CONTRAST"))
  {
    use_contrast = 1 ;
    fprintf(stderr, "using contrast to find labels...\n") ;
  }
  else if (!strcmp(option, "WRITE_MEAN"))
  {
    gca_mean_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "writing gca means to %s...\n", gca_mean_fname) ;
  }
  else if (!strcmp(option, "PRIOR"))
  {
    min_prior = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using prior threshold %2.2f\n", min_prior) ;
  }
  else if (!strcmp(option, "SPACING"))
  {
    min_spacing = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using min GCA spacing %2.0f...\n", min_spacing) ;
  }
  else if (!strcmp(option, "DT"))
  {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt = %2.2e\n", parms.dt) ;
  }
  else if (!strcmp(option, "TOL"))
  {
    tol = parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "tol = %2.2e\n", parms.tol) ;
  }
  else if (!strcmp(option, "CENTER"))
  {
    center = 1 ;
    fprintf(stderr, "using GCA centroid as origin of transform\n") ;
  }
  else if (!strcmp(option, "NOSCALE"))
  {
    noscale = 1 ;
    fprintf(stderr, "disabling scaling...\n") ;
  }
  else if (!strcmp(option, "NUM"))
  {
    num_xforms = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "finding a total of %d linear transforms\n", num_xforms) ;
  }
  else if (!strcmp(option, "AREA"))
  {
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_area = %2.2f\n", parms.l_area) ;
  }
  else if (!strcmp(option, "NLAREA"))
  {
    parms.l_nlarea = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nlarea = %2.2f\n", parms.l_nlarea) ;
  }
  else if (!strcmp(option, "LEVELS"))
  {
    parms.levels = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "levels = %d\n", parms.levels) ;
  }
  else if (!strcmp(option, "INTENSITY") || !strcmp(option, "CORR"))
  {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_intensity = %2.2f\n", parms.l_intensity) ;
  }
  else if (!stricmp(option, "reduce"))
  {
    nreductions = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "reducing input images %d times before aligning...\n",
            nreductions) ;
  }
  else if (!stricmp(option, "nsamples"))
  {
    nsamples = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using %d samples of GCA...\n", nsamples) ;
  }
  else if (!stricmp(option, "steps"))
  {
    max_angles = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "taking %d angular steps...\n", max_angles) ;
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
    fprintf(stderr, "using previously computed transform %s\n", argv[2]) ;
    transform_loaded = 1 ;
    break ;
  case 'B':
    blur_sigma = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "blurring input image with sigma=%2.3f\n", blur_sigma);
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'S':
    parms.sigma = atof(argv[2]) ;
    fprintf(stderr, "using sigma=%2.3f as upper bound on blurring.\n", 
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
    fprintf(stderr, "niterations = %d\n", parms.niterations) ;
    break ;
  case 'W':
    parms.write_iterations = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "write iterations = %d\n", parms.write_iterations) ;
    Gdiag |= DIAG_WRITE ;
    break ;
  case 'M':
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "momentum = %2.2f\n", parms.momentum) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
