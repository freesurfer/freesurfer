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

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static MRI *align_with_average(MRI *mri_src, MRI *mri_avg) ;
char *FileNameRemoveExtension(char *in_fname, char *out_fname) ;
static MRI *align_pca(MRI *mri_src, MRI *mri_avg) ;
static MRI *apply_pca(MRI *mri_in, MRI *mri_ref, MRI *mri_reg,
                         MATRIX *m_in_evectors, double in_means[3],
                         MATRIX *m_ref_evectors, double ref_means[3]) ;
MRI *MRIscaleMeanIntensities(MRI *mri_src, MRI *mri_ref, MRI *mri_dst) ;

char *Progname ;
static int align = 1 ;
static MORPH_PARMS  parms ;

static void usage_exit(int code) ;

static double tx = 0.0 ;
static double ty = 0.0 ;
static double tz = 0.0 ;
static double rzrot = 0.0 ;
static double rxrot = 0.0 ;
static double ryrot = 0.0 ;
static int thresh_low = 0 ;

static int conform = 1 ;


int
main(int argc, char *argv[])
{
  char   **av, fname[100] ;
  int    ac, nargs, i ;
  MRI    *mri_src, *mri_avg = NULL, *mri_tmp ;
  char   *in_fname, *out_fname ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.dt = 1e-6 ;
  parms.tol = 1e-6 ;
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

  if (argc < 3)
    usage_exit(1) ;

  out_fname = argv[argc-1] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;

  for (i = 1 ; i < argc-1 ; i++)
  {
    in_fname = argv[i] ;
    fprintf(stderr, "reading %s...\n", in_fname) ;

    mri_src = MRIread(in_fname) ;
    if (!mri_src)
      ErrorExit(Gerror, "%s: MRIread(%s) failed", Progname, in_fname) ;
    if (conform)
    {
      MRI *mri_tmp ;

      fprintf(stderr, "embedding and interpolating volume\n") ;
      mri_tmp = MRIconform(mri_src) ;
      /*      MRIfree(&mri_src) ;*/
      mri_src = mri_tmp ;
    }

    if (i == 2)
    {
      if (!FZERO(tx) || !FZERO(ty) || !FZERO(tz))
      {
        MRI *mri_tmp ;

        fprintf(stderr, "translating second volume by (%2.1f, %2.1f, %2.1f)\n",
                tx, ty, tz) ;
        mri_tmp = MRItranslate(mri_src, NULL, tx, ty, tz) ;
        MRIfree(&mri_src) ;
        mri_src = mri_tmp ;
      }
      if (!FZERO(rzrot))
      {
        MRI *mri_tmp ;

        fprintf(stderr,
                "rotating second volume by %2.1f degrees around Z axis\n",
                (float)DEGREES(rzrot)) ;
        mri_tmp = MRIrotateZ(mri_src, NULL, rzrot) ;
        MRIfree(&mri_src) ;
        mri_src = mri_tmp ;
      }
      if (!FZERO(rxrot))
      {
        MRI *mri_tmp ;

        fprintf(stderr,
                "rotating second volume by %2.1f degrees around X axis\n",
                (float)DEGREES(rxrot)) ;
        mri_tmp = MRIrotateX(mri_src, NULL, rxrot) ;
        MRIfree(&mri_src) ;
        mri_src = mri_tmp ;
      }
      if (!FZERO(ryrot))
      {
        MRI *mri_tmp ;

        fprintf(stderr,
                "rotating second volume by %2.1f degrees around Y axis\n",
                (float)DEGREES(ryrot)) ;
        mri_tmp = MRIrotateY(mri_src, NULL, ryrot) ;
        MRIfree(&mri_src) ;
        mri_src = mri_tmp ;
      }
#if 1
      if (!FZERO(rxrot) || !FZERO(ryrot) || !FZERO(rzrot))
        MRIwrite(mri_src, "/disk2/mri/tamily/mri/tmp") ;
#endif
    }
    mri_src->xsize = mri_src->ysize = mri_src->zsize = mri_src->thick = 1.0f ;
    mri_src->imnr0 = 1 ; mri_src->imnr1 = mri_src->depth ;
    if (align && mri_avg)  /* don't align the first time */
    {
      mri_tmp = align_with_average(mri_src, mri_avg) ;
      MRIfree(&mri_src) ;
      mri_src = mri_tmp ;
    }

    mri_avg = MRIaverage(mri_src, i-1, mri_avg) ;
    MRIfree(&mri_src) ;
  }
  fprintf(stderr, "writing to %s...\n", out_fname) ;
  MRIwrite(mri_avg, out_fname) ;
  MRIfree(&mri_avg) ;
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
    fprintf(stderr, "using dt = %2.3e\n", parms.dt) ;
  }
  else if (!stricmp(option, "tol"))
  {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using tol = %2.3e\n", parms.tol) ;
  }
  else if (!stricmp(option, "conform"))
  {
    conform = 1 ;
    fprintf(stderr, "interpolating volume to be isotropic 1mm^3\n") ;
  }
  else if (!stricmp(option, "noconform"))
  {
    conform = 0 ;
    fprintf(stderr, "inhibiting isotropic volume interpolation\n") ;
  }
  else switch (toupper(*option))
  {
  case 'W':
    parms.write_iterations = atoi(argv[2]) ;
    Gdiag |= DIAG_WRITE ;
    nargs = 1 ;
    fprintf(stderr, "writing snapshots every %d iterations\n",
            parms.write_iterations) ;
    break ;
  case 'T':
    tx = atof(argv[2]) ; ty = atof(argv[3]) ; tz = atof(argv[4]) ;
    nargs = 3 ;
    break ;
  case 'R':
    rxrot = RADIANS(atof(argv[2])) ;
    ryrot = RADIANS(atof(argv[3])) ;
    rzrot = RADIANS(atof(argv[4])) ;
    nargs = 3 ;
    break ;
  case 'M':
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using momentum = %2.3f\n", parms.momentum) ;
    break ;
  case 'A':
    align = 1 ;
    fprintf(stderr, "aligning volumes before averaging...\n") ;
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
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
  printf("usage: %s [options] <volume> ... <output volume>\n", Progname) ;
  printf("\t-a    rigid alignment of input volumes before averaging\n") ;
  exit(code) ;
}


static MRI *
align_with_average(MRI *mri_src, MRI *mri_avg)
{
  MRI     *mri_aligned, *mri_in_red, *mri_ref_red, *mri_tmp ;
  MRI     *mri_in_windowed, *mri_ref_windowed ;

  fprintf(stderr, "initializing alignment using PCA...\n") ;
  if (Gdiag & DIAG_WRITE)
  {
    MRIwriteImageViews(mri_avg, "ref", 400) ;
    MRIwriteImageViews(mri_src, "before_pca", 400) ;
  }

  mri_aligned = align_pca(mri_src, mri_avg) ;
  if (Gdiag & DIAG_WRITE)
    MRIwriteImageViews(mri_aligned, "after_pca", 400) ;

  fprintf(stderr, "aligning volume with average...\n") ;
#define WINDOW_IMAGES 0
#if WINDOW_IMAGES
  mri_in_windowed = MRIwindow(mri_aligned, NULL, WINDOW_HANNING,127,127,127,100.0f);
  mri_ref_windowed = MRIwindow(mri_avg,NULL,WINDOW_HANNING,127,127,127,100.0f);
#else
  mri_in_windowed = MRIcopy(mri_aligned, NULL) ;
  mri_ref_windowed = MRIcopy(mri_avg, NULL) ;
#endif

  MRIscaleMeanIntensities(mri_in_windowed, mri_ref_windowed, mri_in_windowed);
#define USE_REDUCED_IMAGES 0
#if USE_REDUCED_IMAGES
  mri_in_red = MRIreduceByte(mri_in_windowed, NULL) ;
  mri_ref_red = MRIreduceByte(mri_ref_windowed,NULL);
#else
  mri_in_red = MRIcopy(mri_in_windowed, NULL) ;
  mri_ref_red = MRIcopy(mri_ref_windowed,NULL);
#endif
  parms.mri_ref = mri_avg ;
  parms.mri_in = mri_aligned ;  /* for diagnostics */
  MRIrigidAlign(mri_in_red, mri_ref_red, &parms) ;

  fprintf(stderr, "transforming input volume...\n") ;
  MatrixPrint(stderr, parms.lta->xforms[0].m_L) ;
  fprintf(stderr, "\n") ;

  mri_tmp = MRIcopy(mri_aligned, NULL) ;
  mri_aligned = MRIlinearTransform(mri_tmp, NULL, parms.lta->xforms[0].m_L) ;
  MRIfree(&mri_tmp) ;
  MRIfree(&mri_in_red) ; MRIfree(&mri_ref_red) ;

  return(mri_aligned) ;
}


static MRI *
align_pca(MRI *mri_in, MRI *mri_ref)
{
  int    row, col ;
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

{
  int i ;

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
}
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
      fprintf(stderr, "WARNING: mirror image detected in %dth eigenvector!\n",
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
  return(apply_pca(mri_in, mri_ref, NULL,
                      m_in_evectors, in_means,
                      m_ref_evectors, ref_means)) ;
}

static MRI *
apply_pca(MRI *mri_in, MRI *mri_ref, MRI *mri_reg,
             MATRIX *m_in_evectors, double in_means[3],
             MATRIX *m_ref_evectors, double ref_means[3])
{
  float   dx, dy, dz ;
  MATRIX  *mRot, *m_in_T, *mOrigin ;
  double  x_angle, y_angle, z_angle, r11, r21, r31, r32, r33, cosy ;
  MRI     *mri_tmp ;

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

  mri_tmp = MRItranslate(mri_in, NULL, dx, dy, dz) ;

  fprintf(stderr, "rotating volume by (%2.2f, %2.2f, %2.2f)\n",
          DEGREES(x_angle), DEGREES(y_angle), DEGREES(z_angle)) ;
#if 0
  if (Gdiag & DIAG_WRITE)
    MRIwriteImageViews(mri_tmp, "after_trans", 400) ;
  fprintf(stderr, "rotation matrix:\n") ;
  MatrixPrint(stderr, mRot) ;
  mri_reg = MRIrotate_I(mri_tmp, mri_reg, mRot, mOrigin) ;
#else
  mri_reg = MRIcopy(mri_tmp, NULL) ;
#endif
  MRIfree(&mri_tmp) ;

  if (Gdiag & DIAG_WRITE)
    MRIwrite(mri_reg, "reg.mgh") ;

  MatrixFree(&mRot) ;
  VectorFree(&mOrigin) ;
  return(mri_reg) ;
}
char *
FileNameRemoveExtension(char *in_fname, char *out_fname)
{
  char *dot ;

  if (out_fname != in_fname)
    strcpy(out_fname, in_fname) ;
  dot = strrchr(out_fname, '.') ;
  if (dot)
    *dot = 0 ;
  return(out_fname) ;
}

MRI *
MRIscaleMeanIntensities(MRI *mri_src, MRI *mri_ref, MRI *mri_dst)
{
  int    width, height, depth, x, y, z, val ;
  double ref_mean, src_mean, nref_vox, nsrc_vox, scale ;

  mri_dst = MRIcopy(mri_src, mri_dst) ;

  width = mri_dst->width ; height = mri_dst->height ; depth = mri_dst->depth;

  nref_vox = nsrc_vox = src_mean = ref_mean = 0.0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (MRIvox(mri_ref,x,y,z) > 10)
        {
          nref_vox++ ;
          ref_mean += (double)MRIvox(mri_ref, x, y, z) ;
        }
        if (MRIvox(mri_src,x,y,z) > 10)
        {
          src_mean += (double)MRIvox(mri_src, x, y, z) ;
          nsrc_vox++ ;
        }
      }
    }
  }

  ref_mean /= nref_vox ; src_mean /= nsrc_vox ;
  fprintf(stderr, "mean brightnesses: ref = %2.1f, in = %2.1f\n",
          ref_mean, src_mean) ;
  scale = ref_mean / src_mean ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        val = MRIvox(mri_src, x, y, z) ;
        val = nint(val*scale) ;
        if (val > 255)
          val = 255 ;
        MRIvox(mri_src, x, y, z) = val ;
      }
    }
  }

  return(mri_dst) ;
}
