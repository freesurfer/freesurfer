

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

char         *Progname ;
static MORPH_PARMS  parms ;

static double tx = 0.0 ;
static double ty = 0.0 ;
static double tz = 0.0 ;
static double rzrot = 0.0 ;
static double rxrot = 0.0 ;
static double ryrot = 0.0 ;

static int get_option(int argc, char *argv[]) ;
static int register_mri(MRI *mri_in, MRI *mri_ref, MP *parms) ;
static int order_eigenvectors(MATRIX *m_src_evectors, MATRIX *m_dst_evectors) ;
static float window_size = 0 ;
static unsigned char thresh_low = 40 ;
static int binarize = 1 ;

static char *var_fname = NULL ;

#if 0
static unsigned char thresh_hi = 120 ;
#endif

static MATRIX *pca_matrix(MATRIX *m_in_evectors, double in_means[3],
                         MATRIX *m_ref_evectors, double ref_means[3]) ;
static MATRIX *compute_pca(MRI *mri_in, MRI *mri_ref) ;
#if 0
static int init_scaling(MRI *mri_in, MRI *mri_ref, MATRIX *m_L) ;
static int init_translation(MRI *mri_in, MRI *mri_ref, MATRIX *m_L);
#endif

static int nreductions = 1 ;
static int num_xforms = 1 ;
static int transform_loaded = 0 ;

static double blur_sigma = 2.0f ;

/* 
   command line consists of three inputs:

   argv[1]  - directory containing 'canonical' brain
   argv[2]  - directory containing brain to be registered
   argv[3]  - directory in which to write out registered brain.
*/

int
main(int argc, char *argv[])
{
  char         ref_fname[100], *in_fname, *out_fname, fname[100], **av ;
  MRI          *mri_ref, *mri_in ;
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  parms.l_intensity = 1.0f ;
  parms.niterations = 100 ;
  parms.levels = -1 ;   /* use default */
  parms.dt = 1e-6 ;  /* was 5e-6 */
  parms.tol = INTEGRATION_TOL*5 ;

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
  strcpy(ref_fname, argv[2]) ;
#if 0
  if (strchr(ref_fname, '#') == NULL)
    strcat(ref_fname, "#-2") ;
#endif
  out_fname = argv[3] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
  Gdiag |= DIAG_WRITE ;
  fprintf(stderr, "logging results to %s.log\n", parms.base_name) ;

  TimerStart(&start) ;
  fprintf(stderr, "reading '%s'...\n", ref_fname) ;
  fflush(stderr) ;
  mri_ref = MRIread(ref_fname) ;
  if (!mri_ref)
    ErrorExit(ERROR_NOFILE, "%s: could not open reference volume %s.\n",
              Progname, ref_fname) ;

  if (var_fname)  /* read in a volume of standard deviations */
  {
    MRI *mri_var, *mri_tmp ;

    fprintf(stderr, "reading '%s'...\n", var_fname) ;
    mri_var = MRIread(var_fname) ;
    if (!mri_var)
      ErrorExit(ERROR_NOFILE, "%s: could not open variance volume %s.\n",
                Progname, var_fname) ;
    mri_tmp = MRIconcatenateFrames(mri_ref, mri_var, NULL) ;
    MRIfree(&mri_var) ; MRIfree(&mri_ref) ;
    mri_ref = mri_tmp ;
  }
  fprintf(stderr, "reading '%s'...\n", in_fname) ;
  fflush(stderr) ;
  mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
              Progname, in_fname) ;


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
  register_mri(mri_in, mri_ref, &parms) ;
  
  fprintf(stderr, "writing output transformation to %s...\n", out_fname) ;
  LTAwrite(parms.lta, out_fname) ;
  if (mri_ref)
    MRIfree(&mri_ref) ;
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
  else if (!strcmp(option, "DT"))
  {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt = %2.2e\n", parms.dt) ;
  }
  else if (!strcmp(option, "TOL"))
  {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "tol = %2.2e\n", parms.tol) ;
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
  else if (!strcmp(option, "WINDOW"))
  {
    window_size = atof(argv[2]) ;
    fprintf(stderr, "applying Hanning window (R=%2.1f) to images...\n",
            window_size) ;
    nargs = 1 ;
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
  else if (!stricmp(option, "thresh"))
  {
    thresh_low = atoi(argv[2]) ;
#if 1
    fprintf(stderr, "setting threshold to %d\n", thresh_low) ;
    nargs = 1 ;
#else
    thresh_hi = atoi(argv[3]) ;
    fprintf(stderr, "thresholds set to %d --> %d\n", thresh_low, thresh_hi) ;
    nargs = 2 ;
#endif
  }
  else if (!stricmp(option, "reduce"))
  {
    nreductions = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "reducing input images %d times before aligning...\n",
            nreductions) ;
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
  case 'V':
    var_fname = argv[2] ;
    fprintf(stderr, "reading variance image from %s...\n", var_fname) ;
    nargs = 1 ;
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

#if 1
static int
register_mri(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms)
{
  MRI     *mri_in_red, *mri_ref_red ;
  MRI     *mri_in_windowed, *mri_ref_windowed, *mri_in_tmp, *mri_ref_tmp ;
  int     i ;
  MATRIX  *m_L ;

  MRIscaleMeanIntensities(mri_in, mri_ref, mri_in);
  fprintf(stderr, "initializing alignment using PCA...\n") ;
  if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
  {
    MRIwriteImageViews(mri_ref, "ref", 400) ;
    MRIwriteImageViews(mri_in, "before_pca", 400) ;
  }

  m_L = compute_pca(mri_in, mri_ref) ; 

#if 0
  init_scaling(mri_in, mri_ref, m_L) ;
  init_translation(mri_in, mri_ref, m_L) ; /* in case PCA failed */
#endif

  /* convert it to RAS mm coordinates */
  MRIvoxelXformToRasXform(mri_in, mri_ref, m_L, m_L) ;

  if (Gdiag & DIAG_SHOW)
  {
    printf("initial transform:\n") ;
    MatrixPrint(stdout, m_L) ;
  }
  if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
  {
    MRI *mri_aligned ;
    
    mri_aligned = MRIapplyRASlinearTransform(mri_in, NULL, m_L) ;
    MRIwriteImageViews(mri_aligned, "after_pca", 400) ;
    MRIfree(&mri_aligned) ;
  }

  fprintf(stderr, "aligning volume with average...\n") ;

  if (window_size > 0)
  {
    double in_means[3], ref_means[3] ;

    MRIcenterOfMass(mri_in, in_means, 0) ;
    MRIcenterOfMass(mri_ref, ref_means, 0) ;
    printf("windowing ref around (%d, %d, %d) and input around (%d, %d, %d)\n",
           nint(ref_means[0]), nint(ref_means[1]), nint(ref_means[2]),
           nint(in_means[0]), nint(in_means[1]), nint(in_means[2])) ;
    mri_in_windowed = 
      MRIwindow(mri_in, NULL, WINDOW_HANNING,nint(in_means[0]),
                nint(in_means[1]), nint(in_means[2]),window_size);
    mri_ref_windowed = 
      MRIwindow(mri_ref,NULL,WINDOW_HANNING,nint(ref_means[0]),
                nint(ref_means[1]), nint(ref_means[2]),window_size);
    mri_in = mri_in_windowed ; mri_ref = mri_ref_windowed ;
  }


  mri_in_red = mri_in_tmp = MRIcopy(mri_in, NULL) ;
  mri_ref_red = mri_ref_tmp = MRIcopy(mri_ref, NULL) ;
  for (i = 0 ; i < nreductions ; i++)
  {
    mri_in_red = MRIreduceByte(mri_in_tmp, NULL) ;
    mri_ref_red = MRIreduceByte(mri_ref_tmp,NULL);
    MRIfree(&mri_in_tmp); MRIfree(&mri_ref_tmp) ;
    mri_in_tmp = mri_in_red ; mri_ref_tmp = mri_ref_red ;
  }
  parms->mri_ref = mri_ref ; 
  parms->mri_in = mri_in ;  /* for diagnostics */
  MRIrigidAlign(mri_in_red, mri_ref_red, parms, m_L) ;

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

  MRIfree(&mri_in_red) ; MRIfree(&mri_ref_red) ;

  return(NO_ERROR) ;
}
static MATRIX *
compute_pca(MRI *mri_in, MRI *mri_ref)
{
  int    row, col, i ; 
  float  dot ;
  MATRIX *m_ref_evectors = NULL, *m_in_evectors = NULL ;
  float  in_evalues[3], ref_evalues[3] ;
  double  ref_means[3], in_means[3] ;
#if 0
  MRI     *mri_in_windowed, *mri_ref_windowed ;

  mri_in_windowed = MRIwindow(mri_in, NULL, WINDOW_HANNING,127,127,127,100.0f);
  mri_ref_windowed = MRIwindow(mri_ref,NULL,WINDOW_HANNING,127,127,127,100.0f);
  if (Gdiag & DIAG_WRITE)
  {
    MRIwriteImageViews(mri_in_windowed, "in_windowed", 400) ;
    MRIwriteImageViews(mri_ref_windowed, "ref_windowed", 400) ;
  }
#endif

  if (!m_ref_evectors)
    m_ref_evectors = MatrixAlloc(3,3,MATRIX_REAL) ;
  if (!m_in_evectors)
    m_in_evectors = MatrixAlloc(3,3,MATRIX_REAL) ;

  if (binarize)
  {
    MRIbinaryPrincipleComponents(mri_ref, m_ref_evectors, ref_evalues, 
                                 ref_means, thresh_low);
    MRIbinaryPrincipleComponents(mri_in,m_in_evectors,in_evalues,in_means,
                                 thresh_low);
  }
  else
  {
    MRIprincipleComponents(mri_ref, m_ref_evectors, ref_evalues, ref_means, 
                           thresh_low);
    MRIprincipleComponents(mri_in,m_in_evectors,in_evalues,in_means,
                           thresh_low);
  }
  fprintf(stderr, "before ordering....\n") ;
  fprintf(stderr, "ref_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    fprintf(stderr, "\t\t%2.2f    %2.2f    %2.2f\n",
            m_ref_evectors->rptr[i][1],
            m_ref_evectors->rptr[i][2],
            m_ref_evectors->rptr[i][3]) ;

  fprintf(stderr, "\nin_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    fprintf(stderr, "\t\t%2.2f    %2.2f    %2.2f\n",
            m_in_evectors->rptr[i][1],
            m_in_evectors->rptr[i][2],
            m_in_evectors->rptr[i][3]) ;

  order_eigenvectors(m_in_evectors, m_in_evectors) ;
  order_eigenvectors(m_ref_evectors, m_ref_evectors) ;

  /* check to make sure eigenvectors aren't reversed */
  for (col = 1 ; col <= 3 ; col++)
  {
#if 0
    float theta ;
#endif

    for (dot = 0.0f, row = 1 ; row <= 3 ; row++)
      dot += m_in_evectors->rptr[row][col] * m_ref_evectors->rptr[row][col] ;

    if (dot < 0.0f)
    {
      fprintf(stderr, "WARNING: mirror image detected in eigenvector #%d\n",
              col) ;
      dot *= -1.0f ;
      for (row = 1 ; row <= 3 ; row++)
        m_in_evectors->rptr[row][col] *= -1.0f ;
    }
#if 0
    theta = acos(dot) ;
    fprintf(stderr, "angle[%d] = %2.1f\n", col, DEGREES(theta)) ;
#endif
  }
  fprintf(stderr, "ref_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    fprintf(stderr, "\t\t%2.2f    %2.2f    %2.2f\n",
            m_ref_evectors->rptr[i][1],
            m_ref_evectors->rptr[i][2],
            m_ref_evectors->rptr[i][3]) ;

  fprintf(stderr, "\nin_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    fprintf(stderr, "\t\t%2.2f    %2.2f    %2.2f\n",
            m_in_evectors->rptr[i][1],
            m_in_evectors->rptr[i][2],
            m_in_evectors->rptr[i][3]) ;

  return(pca_matrix(m_in_evectors, in_means,m_ref_evectors, ref_means)) ;
}

static int
order_eigenvectors(MATRIX *m_src_evectors, MATRIX *m_dst_evectors)
{
  int    row, col, xcol, ycol, zcol ;
  double mx ;

  if (m_src_evectors == m_dst_evectors)
    m_src_evectors = MatrixCopy(m_src_evectors, NULL) ;

  /* find columx with smallest dot product with unit x vector */
  mx = fabs(*MATRIX_RELT(m_src_evectors, 1, 1)) ; xcol = 1 ;
  for (col = 2 ; col <= 3 ; col++)
    if (fabs(*MATRIX_RELT(m_src_evectors, 1, col)) > mx)
    {
      xcol = col ;
      mx = fabs(*MATRIX_RELT(m_src_evectors, 1, col)) ;
    }

  mx = fabs(*MATRIX_RELT(m_src_evectors, 2, 1)) ; ycol = 1 ;
  for (col = 2 ; col <= 3 ; col++)
    if (*MATRIX_RELT(m_src_evectors, 2, col) > mx)
    {
      ycol = col ;
      mx = fabs(*MATRIX_RELT(m_src_evectors, 2, col)) ;
    }

  mx = fabs(*MATRIX_RELT(m_src_evectors, 3, 1)) ; zcol = 1 ;
  for (col = 2 ; col <= 3 ; col++)
    if (fabs(*MATRIX_RELT(m_src_evectors, 3, col)) > mx)
    {
      zcol = col ;
      mx = fabs(*MATRIX_RELT(m_src_evectors, 3, col)) ;
    }

  for (row = 1 ; row <= 3 ; row++)
  {
    *MATRIX_RELT(m_dst_evectors,row,1) = *MATRIX_RELT(m_src_evectors,row,xcol);
    *MATRIX_RELT(m_dst_evectors,row,2) = *MATRIX_RELT(m_src_evectors,row,ycol);
    *MATRIX_RELT(m_dst_evectors,row,3) = *MATRIX_RELT(m_src_evectors,row,zcol);
  }
  return(NO_ERROR) ;
}


static MATRIX *
pca_matrix(MATRIX *m_in_evectors, double in_means[3],
           MATRIX *m_ref_evectors, double ref_means[3])
{
  float   dx, dy, dz ;
  MATRIX  *mRot, *m_in_T, *mOrigin, *m_L, *m_R, *m_T, *m_tmp ;
  double  x_angle, y_angle, z_angle, r11, r21, r31, r32, r33, cosy ;
  int     row, col ;

  m_in_T = MatrixTranspose(m_in_evectors, NULL) ;
  mRot = MatrixMultiply(m_ref_evectors, m_in_T, NULL) ;

  r11 = mRot->rptr[1][1] ;
  r21 = mRot->rptr[2][1] ;
  r31 = mRot->rptr[3][1] ;
  r32 = mRot->rptr[3][2] ;
  r33 = mRot->rptr[3][3] ;
  y_angle = atan2(-r31, sqrt(r11*r11+r21*r21)) ;
  cosy = cos(y_angle) ;
  z_angle = atan2(r21 / cosy, r11 / cosy) ;
  x_angle = atan2(r32 / cosy, r33 / cosy) ;

#define MAX_ANGLE  (RADIANS(25))
  if (fabs(x_angle) > MAX_ANGLE || fabs(y_angle) > MAX_ANGLE || 
      fabs(z_angle) > MAX_ANGLE)
  {
    MATRIX *m_I ;

    /*    MatrixFree(&m_in_T) ; MatrixFree(&mRot) ;*/
    fprintf(stderr, 
       "eigenvector swap detected (%2.0f, %2.0f, %2.0f): ignoring rotational PCA...\n",
            DEGREES(x_angle), DEGREES(y_angle), DEGREES(z_angle)) ;

    m_I = MatrixIdentity(3, NULL) ;
    MatrixCopy(m_I, mRot) ;
    MatrixFree(&m_I) ;
    x_angle = y_angle = z_angle = 0.0 ;
  }

  mOrigin = VectorAlloc(3, MATRIX_REAL) ;
  mOrigin->rptr[1][1] = ref_means[0] ;
  mOrigin->rptr[2][1] = ref_means[1] ;
  mOrigin->rptr[3][1] = ref_means[2] ;

  fprintf(stderr, "reference volume center of mass at (%2.1f,%2.1f,%2.1f)\n",
          ref_means[0], ref_means[1], ref_means[2]) ;
  fprintf(stderr, "input volume center of mass at     (%2.1f,%2.1f,%2.1f)\n",
          in_means[0], in_means[1], in_means[2]) ;
  dx = ref_means[0] - in_means[0] ;
  dy = ref_means[1] - in_means[1] ;
  dz = ref_means[2] - in_means[2] ;

  fprintf(stderr, "translating volume by %2.1f, %2.1f, %2.1f\n",
          dx, dy, dz) ;
  fprintf(stderr, "rotating volume by (%2.2f, %2.2f, %2.2f)\n",
          DEGREES(x_angle), DEGREES(y_angle), DEGREES(z_angle)) ;

  /* build full rigid transform */
  m_R = MatrixAlloc(4,4,MATRIX_REAL) ;
  m_T = MatrixAlloc(4,4,MATRIX_REAL) ;
  for (row = 1 ; row <= 3 ; row++)
  {
    for (col = 1 ; col <= 3 ; col++)
    {
      *MATRIX_RELT(m_R,row,col) = *MATRIX_RELT(mRot, row, col) ;
    }
    *MATRIX_RELT(m_T,row,row) = 1.0 ;
  }
  *MATRIX_RELT(m_R, 4, 4) = 1.0 ;

  /* translation so that origin is at ref eigenvector origin */
  dx = -ref_means[0] ; dy = -ref_means[1] ; dz = -ref_means[2] ;
  *MATRIX_RELT(m_T, 1, 4) = dx ; *MATRIX_RELT(m_T, 2, 4) = dy ;
  *MATRIX_RELT(m_T, 3, 4) = dz ; *MATRIX_RELT(m_T, 4, 4) = 1 ;
  m_tmp = MatrixMultiply(m_R, m_T, NULL) ;
  *MATRIX_RELT(m_T, 1, 4) = -dx ; *MATRIX_RELT(m_T, 2, 4) = -dy ;
  *MATRIX_RELT(m_T, 3, 4) = -dz ; 
  MatrixMultiply(m_T, m_tmp, m_R) ;

  /* now apply translation to take in centroid to ref centroid */
  dx = ref_means[0] - in_means[0] ; dy = ref_means[1] - in_means[1] ;
  dz = ref_means[2] - in_means[2] ;
  *MATRIX_RELT(m_T, 1, 4) = dx ; *MATRIX_RELT(m_T, 2, 4) = dy ;
  *MATRIX_RELT(m_T, 3, 4) = dz ; *MATRIX_RELT(m_T, 4, 4) = 1 ;

  m_L = MatrixMultiply(m_R, m_T, NULL) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("m_T:\n") ;
    MatrixPrint(stdout, m_T) ;
    printf("m_R:\n") ;
    MatrixPrint(stdout, m_R) ;
    printf("m_L:\n") ;
    MatrixPrint(stdout, m_L) ;
  }
  MatrixFree(&m_R) ; MatrixFree(&m_T) ;

  MatrixFree(&mRot) ;
  VectorFree(&mOrigin) ;
  return(m_L) ;
}
#else
static int
register_mri(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms)
{
  MRI  *mri_in_red, *mri_ref_red ;

  mri_in_red = MRIreduceByte(mri_in, NULL) ;
  mri_ref_red = MRIreduceMeanAndStdByte(mri_ref,NULL);

  /*  parms->write_iterations = 0 ; */
  if (!parms->niterations)
    parms->niterations = 1000 ;
  if (transform_loaded)  /* don't recompute rotation based on neck */
  {
    if (MRIfindNeck(mri_in_red, mri_in_red, thresh_low, thresh_hi, NULL, 
                    -1,NULL) == NULL)
      ErrorExit(Gerror, "%s: could not find subject neck.\n", Progname) ;
    if (MRIfindNeck(mri_ref_red, mri_ref_red, thresh_low, thresh_hi, NULL,1,
                    NULL) == NULL)
      ErrorExit(Gerror, "%s: could not find neck in reference volume.\n", 
                Progname) ;
  }
  else
  {
    if (MRIfindNeck(mri_in_red, mri_in_red, thresh_low, thresh_hi, parms, -1,
                &parms->in_np) == NULL)
      ErrorExit(Gerror, "%s: could not find subject neck.\n", Progname) ;
    if (MRIfindNeck(mri_ref_red, mri_ref_red, thresh_low, thresh_hi, parms, 1,
                &parms->ref_np) == NULL)
      ErrorExit(Gerror, "%s: could not find neck in reference volume.\n", 
                Progname) ;
  }

  parms->mri_ref = mri_ref_red ; 
  parms->mri_in = mri_in_red ;  /* for diagnostics */
  while (parms->lta->num_xforms < num_xforms)
    LTAdivide(parms->lta, mri_in_red) ;
  fprintf(stderr,"computing %d linear transformation%s...\n",
          parms->lta->num_xforms, parms->lta->num_xforms>1?"s":"");
  MRIlinearAlign(mri_in_red, mri_ref_red, parms) ;

  MRIfree(&mri_in_red) ; MRIfree(&mri_ref_red) ;
  return(NO_ERROR) ;
}
#endif

#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_DX   1.2
#define MAX_DY   1.2
#define MAX_DZ   1.2
#define MIN_DX   (1.0/MAX_DX)
#define MIN_DY   (1.0/MAX_DY)
#define MIN_DZ   (1.0/MAX_DZ)
#define MAX_RATIO 1.2

static int
init_scaling(MRI *mri_in, MRI *mri_ref, MATRIX *m_L)
{
  MATRIX      *m_scaling ;
  float       sx, sy, sz, dx, dy, dz ;
  MRI_REGION  in_bbox, ref_bbox ;

  m_scaling = MatrixIdentity(4, NULL) ;

  MRIboundingBox(mri_in, 60, &in_bbox) ;
  MRIboundingBox(mri_ref, 60, &ref_bbox) ;
  sx = (float)ref_bbox.dx / (float)in_bbox.dx ;
  sy = (float)ref_bbox.dy / (float)in_bbox.dy ;
  sz = (float)ref_bbox.dz / (float)in_bbox.dz ;
  dx = (ref_bbox.x+ref_bbox.dx-1)/2 - (in_bbox.x+in_bbox.dx-1)/2 ;
  dy = (ref_bbox.y+ref_bbox.dy-1)/2 - (in_bbox.y+in_bbox.dy-1)/2 ;
  dz = (ref_bbox.z+ref_bbox.dz-1)/2 - (in_bbox.z+in_bbox.dz-1)/2 ;
  
  if (sx > MAX_DX)
    sx = MAX_DX ;
  if (sx < MIN_DX)
    sx = MIN_DX ;
  if (sy > MAX_DY)
    sy = MAX_DY ;
  if (sy < MIN_DY)
    sy = MIN_DY ;
  if (sz > MAX_DZ)
    sz = MAX_DZ ;
  if (sz < MIN_DZ)
    sz = MIN_DZ ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "initial scaling: (%2.2f, %2.2f, %2.2f) <-- "
            "(%d/%d,%d/%d,%d/%d)\n", 
            sx,sy,sz, ref_bbox.dx, in_bbox.dx, ref_bbox.dy, in_bbox.dy, 
            ref_bbox.dz, in_bbox.dz) ;
  *MATRIX_RELT(m_scaling, 1, 1) = sx ;
  *MATRIX_RELT(m_scaling, 2, 2) = sy ;
  *MATRIX_RELT(m_scaling, 3, 3) = sz ;

#if 0
  *MATRIX_RELT(m_L, 1, 4) = dx ;
  *MATRIX_RELT(m_L, 2, 4) = dy ;
  *MATRIX_RELT(m_L, 3, 4) = dz ;
#endif
  MatrixMultiply(m_scaling, m_L, m_L) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
init_translation(MRI *mri_in, MRI *mri_ref, MATRIX *m_L)
{
  MATRIX *m_translation ;
  float  in_means[4], ref_means[4] ;
  double dx, dy, dz ;

  m_translation = MatrixIdentity(4, NULL) ;
  MRIfindCenterOfBrain(mri_in, in_means, in_means+1, in_means+2) ;
  MRIfindCenterOfBrain(mri_ref, ref_means, ref_means+1, ref_means+2) ;
  dx = (double)(ref_means[0] - in_means[0]) * mri_in->thick ;
  dy = (double)(ref_means[1] - in_means[1]) * mri_in->thick ;
  dz = (double)(ref_means[2] - in_means[2]) * mri_in->thick ;
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stderr, "centering template around (%d,%d,%d) and input around"
            " (%d,%d,%d)\n", 
            (int)ref_means[0], (int)ref_means[1], (int)ref_means[2],
            (int)in_means[0], (int)in_means[1], (int)in_means[2]) ;
    fprintf(stderr, "initial translation = (%2.0f, %2.0f, %2.0f).\n",dx,dy,dz);
  }
  *MATRIX_RELT(m_translation, 1, 4) = dx ;
  *MATRIX_RELT(m_translation, 2, 4) = dy ;
  *MATRIX_RELT(m_translation, 3, 4) = dz ;
  MatrixMultiply(m_translation, m_L, m_L) ;
  MatrixFree(&m_translation) ;
  return(NO_ERROR) ;
}
#endif
