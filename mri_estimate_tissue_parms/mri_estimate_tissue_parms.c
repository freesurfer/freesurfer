#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "timer.h"

int main(int argc, char *argv[]) ;
static MRI *compute_residuals(MRI *mri_flash, MRI *mri_T1, MRI *mri_PD) ;
static int get_option(int argc, char *argv[]) ;
#if 0
static MRI *align_with_average(MRI *mri_src, MRI *mri_avg) ;
static MATRIX *align_pca(MRI *mri_src, MRI *mri_avg) ;
static MATRIX *pca_matrix(MATRIX *m_in_evectors, double in_means[3],
                         MATRIX *m_ref_evectors, double ref_means[3]) ;
static int thresh_low = 0 ;
#endif

char *Progname ;
static int align = 1 ;
static int window_flag = 0 ;
static MORPH_PARMS  parms ;

static void usage_exit(int code) ;
static int no_valid_data(MRI **mri_flash, int nvolumes, 
                         int x, int y, int z, float thresh) ;

static float thresh = 25 ;
static int conform = 0 ;
static int sinc_flag = 1;
static int sinchalfwindow = 3;

static char *residual_name = NULL ;

#define MAX_IMAGES 100
#define MIN_T1   5  /* avoid singularity at T1=0 */
#define MIN_ITER 5


static double dM_dT1(double flip_angle, double TR, double PD, double T1) ;
static double dM_dPD(double flip_angle, double TR, double PD, double T1) ;
static double FLASHforwardModel(double flip_angle, double TR, double PD, 
                                double T1) ;
static double estimateVoxelParameters(MRI **mri_flash, int nvolumes, int x, 
                                      int y, int z, MRI *mri_T1, MRI *mri_PD,
                                      double T1, double PD) ;

static double
computeVoxelSSE(MRI **mri_flash, int nflash, int x, 
                int y, int z, double PD, double T1) ;

static int write_iterations = 0 ;

int
main(int argc, char *argv[])
{
  char   **av, fname[STRLEN] ;
  int    ac, nargs, i ;
  MRI    *mri_flash[MAX_IMAGES], *mri_T1, *mri_PD ;
  char   *in_fname, *out_PD_fname, *out_T1_fname ;
  int          msec, minutes, seconds, nvolumes ;
  struct timeb start ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;
  parms.dt = 1e-6 ;
  parms.tol = 1e-5 ;
  parms.momentum = 0.0 ;
  parms.niterations = 20 ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;

  out_T1_fname = argv[argc-2] ;
  out_PD_fname = argv[argc-1] ;
  FileNameOnly(out_T1_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;

  nvolumes = argc-3 ;
  printf("reading %d FLASH volumes.\n", nvolumes) ;
  for (i = 1 ; i < argc-2 ; i++)
  {
    in_fname = argv[i] ;
    printf("reading %s...", in_fname) ;

    mri_flash[i-1] = MRIread(in_fname) ;
    if (!mri_flash[i-1])
      ErrorExit(Gerror, "%s: MRIread(%s) failed", Progname, in_fname) ;
    printf("TR = %2.2f, alpha = %2.2f\n", mri_flash[i-1]->tr, 
           mri_flash[i-1]->flip_angle) ;
    mri_flash[i-1]->flip_angle = RADIANS(mri_flash[i-1]->flip_angle) ;
    if (conform)
    {
      MRI *mri_tmp ;

      printf("embedding and interpolating volume\n") ;
      mri_tmp = MRIconform(mri_flash[i-1]) ;
      /*      MRIfree(&mri_src) ;*/
      mri_flash[i-1] = mri_tmp ;
    }

  }
  mri_T1 = MRIclone(mri_flash[0], NULL) ;
  mri_PD = MRIclone(mri_flash[0], NULL) ;


  {
    double   sse, last_T1, last_PD ;
    int      x, y, z, width, height, depth, total_vox, ignored ;
    struct timeb first_slice ;

    Progname = argv[0] ;
    ErrorInit(NULL, NULL, NULL) ;
    DiagInit(NULL, NULL, NULL) ;
    
    TimerStart(&first_slice) ;

    last_T1 = last_PD = 1000 ; sse = 0.0 ;
    width = mri_T1->width ; height = mri_T1->height ; depth = mri_T1->depth ;
    total_vox = width*depth*height ;
#if 0
    estimateVoxelParameters(mri_flash, nvolumes, width/2, height/2, depth/2,
                            mri_T1, mri_PD, last_T1, last_PD) ;
#endif
    for (ignored = z = 0 ; z < depth ; z++)
    {
      printf("z = %d, ignored = %d, last sse=%2.3f, T1=%2.0f, PD=%2.0f...\n", 
             z, ignored, sse, last_T1, last_PD) ;
      if (z*width*height - ignored > 0);
      {
        int processed = z*width*height - ignored, hours ;

        msec = TimerStop(&first_slice) ;
        seconds = nint((float)msec/1000.0f) ;
        minutes = seconds / 60 ;
        seconds = seconds % 60 ;
        hours = minutes / 60 ;
        minutes = minutes % 60 ;
        printf("%02d:%02d:%02d total processing time ... ", 
               hours,minutes,seconds);
        msec = (int)((float)(total_vox-2*ignored)*msec/(float)processed) ;
        seconds = nint((float)msec/1000.0f) ;
        minutes = seconds / 60 ;
        seconds = seconds % 60 ;
        hours = minutes / 60 ;
        minutes = minutes % 60 ;
        printf("estimate %02d:%02d:%02d remaining.\n", hours,minutes, seconds);
      }
      if (write_iterations > 0 && z > 0 && !(z%write_iterations))
      {
        printf("writing T1 esimates to %s...\n", out_T1_fname) ;
        printf("writing PD estimates to %s...\n", out_PD_fname) ;
        MRIwrite(mri_T1, out_T1_fname) ;
        MRIwrite(mri_PD, out_PD_fname) ;
        if (residual_name) for (i = 0 ; i < nvolumes ; i++)
        {
          MRI *mri_res ;

          mri_res = compute_residuals(mri_flash[i], mri_T1, mri_PD) ;
          sprintf(fname, "%s%d.mnc", residual_name, i) ;
          MRIwrite(mri_res, fname) ;
          MRIfree(&mri_res) ;
        }
      }


      for (y = 0 ; y < height ; y++)
      {
        if (y%32 == 0)  
        printf("z = %d, y = %d, last sse=%2.3f, T1=%2.0f, PD=%2.0f...\n", 
                 z, y, sse, last_T1, last_PD) ;
        for (x = 0 ; x < width ; x++)
        {
#if 0
          for (i = 0 ; i < nvolumes ; i++)
            if (MRISvox(mri_flash[i],x,y,z) > thresh)
              break ;
          if (i >= nvolumes)
#else
            if (no_valid_data(mri_flash, nvolumes, x, y, z, thresh))
#endif
          {
            ignored++ ;
            MRISvox(mri_T1, x, y, z) = MRISvox(mri_PD, x, y, z) = 0 ;
            last_T1 = last_PD = 1000 ;
            continue ;
          }
          sse = estimateVoxelParameters(mri_flash, nvolumes, x, y, z,
                                        mri_T1, mri_PD, last_T1, last_PD) ;

          if (sqrt(sse/(float)nvolumes) < 15)  /* good fit */
          {
#if 0
            last_T1 = MRISvox(mri_T1, x, y, z) ;
            last_PD = MRISvox(mri_PD, x, y, z) ;
#endif
          }
          else
            last_T1 = last_PD = 1000 ;
        }
      }
    }
  }


  printf("writing T1 esimates to %s...\n", out_T1_fname) ;
  printf("writing PD estimates to %s...\n", out_PD_fname) ;
  MRIwrite(mri_T1, out_T1_fname) ;
  MRIwrite(mri_PD, out_PD_fname) ;
  if (residual_name) for (i = 0 ; i < nvolumes ; i++)
  {
    MRI *mri_res ;
    
    mri_res = compute_residuals(mri_flash[i], mri_T1, mri_PD) ;
    sprintf(fname, "%s%d.mnc", residual_name, i) ;
    MRIwrite(mri_res, fname) ;
    MRIfree(&mri_res) ;
  }
  MRIfree(&mri_T1) ; MRIfree(&mri_PD) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("parameter estimation took %d minutes and %d seconds.\n", 
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
  if (!stricmp(option, "dt"))
  {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    printf("using dt = %2.3e\n", parms.dt) ;
  }
  else if (!stricmp(option, "tol"))
  {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("using tol = %2.3e\n", parms.tol) ;
  }
  else if (!stricmp(option, "conform"))
  {
    conform = 1 ;
    printf("interpolating volume to be isotropic 1mm^3\n") ;
  }
  else if (!stricmp(option, "sinc"))
  {
    sinchalfwindow = atoi(argv[2]);
    sinc_flag = 1;
    nargs = 1;
    printf("using sinc interpolation with windowwidth of %d\n",
      2*sinchalfwindow);
  }
  else if (!stricmp(option, "trilinear"))
  {
    sinc_flag = 0;
    printf("using trilinear interpolation\n");
  }
  else if (!stricmp(option, "window"))
  {
    window_flag = 1 ;
    printf("applying hanning window to volumes...\n") ;
  }
  else if (!stricmp(option, "noconform"))
  {
    conform = 0 ;
    printf("inhibiting isotropic volume interpolation\n") ;
  }
  else switch (toupper(*option))
  {
#if 0
  case 'W':
    parms.write_iterations = atoi(argv[2]) ;
    Gdiag |= DIAG_WRITE ;
    nargs = 1 ;
    printf("writing snapshots every %d iterations\n", 
            parms.write_iterations) ;
    break ;
#endif
  case 'W':
    write_iterations = atoi(argv[2]) ;
    printf("writing out intermediate volumes every %d slices...\n",
           write_iterations) ;
    nargs = 1 ;
    break ;
  case 'M':
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    printf("using momentum = %2.3f\n", parms.momentum) ;
    break ;
  case 'R':
    residual_name = argv[2] ;
    printf("writing out residuals to %s...\n", residual_name) ;
    nargs = 1 ;
    break ;
  case 'T':
    thresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("ignoring locations in which all images are less than %2.2f\n", 
           thresh) ;
    break ;
  case 'A':
    align = 1 ;
    printf("aligning volumes before averaging...\n") ;
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    printf("unknown option %s\n", argv[1]) ;
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
usage_exit(int code)
{
  printf("usage: %s [options] <volume> ... <output T1 volume> <output PD volume>\n", Progname) ;
  printf("\t-a    rigid alignment of input volumes before averaging\n") ;
  exit(code) ;
}

#if 0
static MRI *
align_with_average(MRI *mri_src, MRI *mri_avg)
{
  MRI     *mri_aligned, *mri_in_red, *mri_ref_red ;
  MRI     *mri_in_windowed, *mri_ref_windowed, *mri_in_tmp, *mri_ref_tmp ;
  int     i ;
  MATRIX  *m_L ;

  printf("initializing alignment using PCA...\n") ;
  if (Gdiag & DIAG_WRITE)
  {
    MRIwriteImageViews(mri_avg, "ref", 400) ;
    MRIwriteImageViews(mri_src, "before_pca", 400) ;
  }

  m_L = align_pca(mri_src, mri_avg) ; 
  if (Gdiag & DIAG_SHOW)
  {
    printf("initial transform:\n") ;
    MatrixPrint(stdout, m_L) ;
  }
  if (Gdiag & DIAG_WRITE)
  {
    if(sinc_flag)
      mri_aligned = MRIsincTransform(mri_src, NULL, m_L,sinchalfwindow) ;
    else
      mri_aligned = MRIlinearTransform(mri_src, NULL, m_L) ;
    MRIwriteImageViews(mri_aligned, "after_pca", 400) ;
    MRIfree(&mri_aligned) ;
  }

  printf("aligning volume with average...\n") ;

  if (window_flag)
  {
    mri_in_windowed = 
      MRIwindow(mri_src, NULL, WINDOW_HANNING,127,127,127,100.0f);
    mri_ref_windowed = 
      MRIwindow(mri_avg,NULL,WINDOW_HANNING,127,127,127,100.0f);
    mri_src = mri_in_windowed ; mri_avg = mri_ref_windowed ;
  }

  MRIscaleMeanIntensities(mri_src, mri_avg, mri_src);

  mri_in_red = mri_in_tmp = MRIcopy(mri_src, NULL) ;
  mri_ref_red = mri_ref_tmp = MRIcopy(mri_avg, NULL) ;
  for (i = 0 ; i < nreductions ; i++)
  {
    mri_in_red = MRIreduceByte(mri_in_tmp, NULL) ;
    mri_ref_red = MRIreduceByte(mri_ref_tmp,NULL);
    MRIfree(&mri_in_tmp); MRIfree(&mri_ref_tmp) ;
    mri_in_tmp = mri_in_red ; mri_ref_tmp = mri_ref_red ;
  }
  parms.mri_ref = mri_avg ; 
  parms.mri_in = mri_src ;  /* for diagnostics */
  MRIrigidAlign(mri_in_red, mri_ref_red, &parms, m_L) ;

  printf("transforming input volume...\n") ;
  MatrixPrint(stderr, parms.lta->xforms[0].m_L) ;
  printf("\n") ;

  if(sinc_flag)
    mri_aligned = MRIsincTransform(mri_src, NULL, parms.lta->xforms[0].m_L,sinchalfwindow) ;
  else
    mri_aligned = MRIlinearTransform(mri_src, NULL, parms.lta->xforms[0].m_L) ;
  if (Gdiag & DIAG_WRITE)
    MRIwriteImageViews(mri_aligned, "after_alignment", 400) ;
  MRIfree(&mri_in_red) ; MRIfree(&mri_ref_red) ;

  return(mri_aligned) ;
}


static MATRIX *
align_pca(MRI *mri_in, MRI *mri_ref)
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

  MRIprincipleComponents(mri_ref, m_ref_evectors, ref_evalues, 
                         ref_means, thresh_low);
  MRIprincipleComponents(mri_in,m_in_evectors,in_evalues,in_means,thresh_low);

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
      printf("WARNING: mirror image detected in eigenvector #%d\n",
              col) ;
      dot *= -1.0f ;
      for (row = 1 ; row <= 3 ; row++)
        m_in_evectors->rptr[row][col] *= -1.0f ;
    }
#if 0
    theta = acos(dot) ;
    printf("angle[%d] = %2.1f\n", col, DEGREES(theta)) ;
#endif
  }
  printf("ref_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    printf("\t\t%2.2f    %2.2f    %2.2f\n",
            m_ref_evectors->rptr[i][1],
            m_ref_evectors->rptr[i][2],
            m_ref_evectors->rptr[i][3]) ;

  printf("\nin_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    printf("\t\t%2.2f    %2.2f    %2.2f\n",
            m_in_evectors->rptr[i][1],
            m_in_evectors->rptr[i][2],
            m_in_evectors->rptr[i][3]) ;

  return(pca_matrix(m_in_evectors, in_means,m_ref_evectors, ref_means)) ;
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

#define MAX_ANGLE  (RADIANS(30))
  if (fabs(x_angle) > MAX_ANGLE || fabs(y_angle) > MAX_ANGLE || 
      fabs(z_angle) > MAX_ANGLE)
  {
    MatrixFree(&m_in_T) ; MatrixFree(&mRot) ;
    printf("eigenvector swap detected: ignoring PCA...\n") ;
    return(MatrixIdentity(4, NULL)) ;
  }

  mOrigin = VectorAlloc(3, MATRIX_REAL) ;
  mOrigin->rptr[1][1] = ref_means[0] ;
  mOrigin->rptr[2][1] = ref_means[1] ;
  mOrigin->rptr[3][1] = ref_means[2] ;

  printf("reference volume center of mass at (%2.1f,%2.1f,%2.1f)\n",
          ref_means[0], ref_means[1], ref_means[2]) ;
  printf("input volume center of mass at     (%2.1f,%2.1f,%2.1f)\n",
          in_means[0], in_means[1], in_means[2]) ;
  dx = ref_means[0] - in_means[0] ;
  dy = ref_means[1] - in_means[1] ;
  dz = ref_means[2] - in_means[2] ;

  printf("translating volume by %2.1f, %2.1f, %2.1f\n",
          dx, dy, dz) ;
  printf("rotating volume by (%2.2f, %2.2f, %2.2f)\n",
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
#endif



static double
computeVoxelSSE(MRI **mri_flash, int nflash, int x, 
                int y, int z, double PD, double T1)
{
  double    sse, estimate, err ;
  int       i ;
  MRI       *mri ;

  for (sse = 0.0, i = 0 ; i < nflash ; i++)
  {
    mri = mri_flash[i] ;
    estimate = FLASHforwardModel(mri->flip_angle, mri->tr, PD, T1) ;
    err = (MRISvox(mri,x,y,z)- estimate) ;
    sse += err*err ;
  }
  if (!finite(sse))
    DiagBreak() ;
  return(sse) ;
}

static double
FLASHforwardModel(double flip_angle, double TR, double PD, double T1)
{
  double FLASH, E1 ;
  double  CFA, SFA ;


  CFA = cos(flip_angle) ; SFA = sin(flip_angle) ;
  E1 = exp(-TR/T1) ;
      
  FLASH = PD * SFA ;
  if (!DZERO(T1))
    FLASH *= (1-E1)/(1-CFA*E1);
  return(FLASH) ;
}


static double
dM_dT1(double flip_angle, double TR, double PD, double T1)
{
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
dM_dPD(double flip_angle, double TR, double PD, double T1)
{
  double  dPD, E1 ;
  static double  CFA, SFA ;

  CFA = cos(flip_angle) ;
  SFA = sin(flip_angle) ;
  E1 = exp(TR/T1) ;
  dPD = (-1 + E1)*SFA/(E1 - CFA) ;

  return(dPD) ;
}
static double STEP_SIZE = 0.1 ;
static double momentum = 0.8 ;
static double tol = 0.0001 ;

static double
estimateVoxelParameters(MRI **mri_flash, int nvolumes, int x, 
                                      int y, int z, MRI *mri_T1, MRI *mri_PD,
                        double T1, double PD)
{
  double   sse, last_sse, dT1, dPD, ival, estimate, err,
           last_dT1, last_dPD ;
  int      i, niter ;
  MRI      *mri ;

  last_sse = sse = computeVoxelSSE(mri_flash, nvolumes, x, y, z,T1,PD);

  niter = 0 ; last_dT1 = last_dPD = 0.0 ;
  do
  {
    for (dT1 = dPD = 0.0, i = 0 ; i < nvolumes ; i++)
    {
      mri = mri_flash[i] ;
      ival = (double)MRISvox(mri, x, y, z) ;
      estimate = FLASHforwardModel(mri->flip_angle, mri->tr, PD, T1) ;
      err = (ival - estimate) ;
      dT1 += err*dM_dT1(mri->flip_angle, mri->tr, PD, T1) ;
      dPD += err*dM_dPD(mri->flip_angle, mri->tr, PD, T1) ;
    }

    dT1 = STEP_SIZE * dT1 + momentum * last_dT1 ;
    T1 += dT1 ; last_dT1 = dT1 ;
    dPD = STEP_SIZE * dPD + momentum * last_dPD ;
    PD += dPD ; last_dPD = dPD ;
    sse = computeVoxelSSE(mri_flash, nvolumes, x, y, z, PD, T1);
#if 0
    if (Gdiag_no == 1)
      printf("%03d: sse=%2.2f (d=%2.2f), dT1=%2.4f, dPD=%2.4f, T1=%2.1f, PD=%2.1f\n",
             niter, sse, sse-last_sse, dT1, dPD, T1, PD) ;
#endif
    if (T1 < MIN_T1)
    {
      T1 = MIN_T1 ;
      break ;
    }
    if ((last_sse < sse) || FZERO(sse))
      break ;
    if (((last_sse - sse)/last_sse < tol) && niter > MIN_ITER)
      break ;

    last_sse = sse ;
  } while (niter++ < 1000) ;
  MRISvox(mri_T1, x, y, z) = (short)(T1+0.5) ; ;
  MRISvox(mri_PD, x, y, z) = (short)(PD+0.5) ;
  if (MRISvox(mri_T1, x, y, z) == 0)
    DiagBreak() ;
  return(sse) ;
}

static int compare_sort_array(const void *pc1, const void *pc2) ;

static int
no_valid_data(MRI **mri_flash, int nvolumes, int x, int y, int z, float thresh)
{
  int sort_array[MAX_IMAGES], i, nsmall ;
  float median ;

  for (nsmall = 0, i = 0 ; i < nvolumes ; i++)
  {
    sort_array[i] = MRISvox(mri_flash[i], x, y, z) ;
    if (sort_array[i] < thresh/2)
      nsmall++ ;
  }
  if (nsmall >= nvolumes)
    return(1) ;

  qsort(sort_array, nvolumes, sizeof(int), compare_sort_array) ;
  if (ODD(nvolumes))
    median = sort_array[(nvolumes-1)/2] ;
  else
    median = (sort_array[(nvolumes-1)/2] + sort_array[nvolumes/2]) / 2  ;

  return(median < thresh/2) ;
}

static int
compare_sort_array(const void *pc1, const void *pc2)
{
  register int c1, c2 ;

  c1 = *(int *)pc1 ;
  c2 = *(int *)pc2 ;

/*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (c1 > c2)
    return(1) ;
  else if (c1 < c2)
    return(-1) ;

  return(0) ;
}

static MRI *
compute_residuals(MRI *mri_flash, MRI *mri_T1, MRI *mri_PD)
{
  MRI    *mri_res ;
  int    x, y, z, width, depth, height, T1, PD, val ;
  double prediction, rms ;

  mri_res = MRIclone(mri_flash, NULL) ;

  width = mri_PD->width ; height = mri_PD->height ; depth = mri_PD->depth ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        T1 = MRISvox(mri_T1, x, y, z) ; PD = MRISvox(mri_PD, x, y, z) ;
        prediction = FLASHforwardModel(mri_flash->flip_angle, mri_flash->tr,
                                       PD, T1) ;
        
        val = MRISvox(mri_flash, x, y, z) ;
        rms = val - prediction ; rms = sqrt(rms*rms) ;
        MRISvox(mri_res, x, y, z) = (short)nint(rms) ;
      }
    }
  }
  return(mri_res) ;
}
