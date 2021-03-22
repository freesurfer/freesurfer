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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "proto.h"
#include "histo.h"
#include "transform.h"
#include "mrinorm.h"
#include "version.h"
#include "tags.h"
#include "flash.h"


int main(int argc, char *argv[]) ;

static int saturate_PD(MRI *mri_PD, float PDsat) ;
static int normalize_PD(MRI *mri_PD, float target) ;
static int discard_PD(MRI *mri_PD, short thresh, short target) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static int  transform_T1_values_using_joint_pdf(MRI *mri_T1, char *jpdf_name, int invert) ;
static void print_version(void) ;
static int  apply_bias_field(MRI *mri, int nbias, float *bias_coefs[3][2]) ;
#if 0
static double FLASHforwardModel(double flip_angle, double TR, double PD,
                                double T1) ;
#endif

static MRI *MRIsynthesizeWithFAF(MRI *mri_T1, MRI *mri_PD, MRI *mri_dst, double TR, double alpha, double TE,
                                 int nfaf, float *faf_coefs[3][2]) ;
MRI *MRIsynthesizeWeightedVolume(MRI *mri_T1, MRI *mri_PD, float w5, float TR5,
                                 float w30, float TR30, float target_wm, float TE) ;
static MRI *MRIsynthesize(MRI *mri_T1, MRI *mri_PD, MRI *mri_T2star, MRI *mri_dst, double TR, double alpha, double TE) ;
static int remap_T1(MRI *mri_T1, float mean, float scale) ;

const char *Progname ;
static int normalize = 0 ;
static int discard = 0 ;
static int nl_remap_T1 = 0 ;
static float nl_scale = 0.01 ;
static float nl_mean = 950 ;
static char *jpdf_name = NULL ;
static int invert = 0 ;
static float PDsat = 0.0 ;
static char *T2star_fname = NULL ;

static int extract = 0 ;

static int nbias = 0 ;
static float *bias_coefs[3][2] ;

static int nfaf = 0 ;
static float *faf_coefs[3][2] ;

/* optimal class separation (sort of) */
#if 0
static double w30 = 0.7718 ;
static double w5 = -0.6359 ;
#else
static double w30 = 2* 0.9527 ;
static double w5 =  2*-0.3039 ;
#endif

static int use_weighting = 0 ;

int
main(int argc, char *argv[]) {
  char        **av, *out_fname, *T1_fname, *PD_fname ;
  int         ac, nargs ;
  MRI         *mri_T1, *mri_PD, *mri_out, *mri_T2star = NULL ;
  float       TR, TE, alpha ;


  std::string cmdline = getAllInfo(argc, argv, "mri_synthesize");

  nargs = handleVersionOption(argc, argv, "mri_synthesize");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 7)
    usage_exit() ;

  TR = atof(argv[1]) ;
  alpha = atof(argv[2]) ;
  TE = atof(argv[3]) ;
  T1_fname = argv[4] ;
  PD_fname = argv[5] ;
  out_fname = argv[6] ;

  printf("reading T1 volume from %s...\n", T1_fname) ;
  mri_T1 = MRIread(T1_fname) ;
  if (!mri_T1)
    ErrorExit(ERROR_NOFILE, "%s: could not read T1 volume %s", Progname,
              T1_fname) ;
  if (extract) {
    MRI *mri_tmp ;
    int dx, dy, dz ;

    dx = mri_T1->width/2 ;
    dy = mri_T1->height/2 ;
    dz = mri_T1->depth/2 ;
    printf("extracting interior %dx%dx%d\n", dx, dy, dz) ;
    mri_tmp = MRIextract(mri_T1, NULL, dx/2, dy/2, dz/2, dx, dy, dz) ;
    MRIfree(&mri_T1) ;
    mri_T1 = mri_tmp ;
  }
  if (jpdf_name) {
    transform_T1_values_using_joint_pdf(mri_T1, jpdf_name, invert) ;
  }

  printf("reading PD volume from %s...\n", PD_fname) ;
  mri_PD = MRIread(PD_fname) ;
  if (!mri_PD)
    ErrorExit(ERROR_NOFILE, "%s: could not read PD volume %s", Progname,
              PD_fname) ;

  if (T2star_fname != NULL) {
    printf("reading T2* volume from %s...\n", T2star_fname) ;
    mri_T2star = MRIread(T2star_fname) ;
    if (!mri_T2star)
      ErrorExit(ERROR_NOFILE, "%s: could not read T2* volume %s", Progname,
                T2star_fname) ;
  }

  if (PDsat > 0)
    saturate_PD(mri_PD, PDsat) ;
  if (extract) {
    MRI *mri_tmp ;
    int dx, dy, dz ;

    dx = mri_PD->width/2 ;
    dy = mri_PD->height/2 ;
    dz = mri_PD->depth/2 ;
    mri_tmp = MRIextract(mri_PD, NULL, dx/2, dy/2, dz/2, dx, dy, dz) ;
    MRIfree(&mri_PD) ;
    mri_PD = mri_tmp ;
  }
  if (use_weighting) {
    mri_out = MRIsynthesizeWeightedVolume(mri_T1, mri_PD, w5, TR, w30, TR, 110,TE);
  } else {
    printf("synthesizing volume with TR=%2.1f msec, TE=%2.1f msec, and alpha=%2.2f degrees...\n",
           TR, TE, alpha) ;
    if (normalize)
      normalize_PD(mri_PD, 1000) ;
    if (discard)
      discard_PD(mri_PD, 250, 1500) ;
    if (nl_remap_T1)
      remap_T1(mri_T1, nl_mean, nl_scale) ;
    if (nfaf > 0)
      mri_out = MRIsynthesizeWithFAF(mri_T1, mri_PD, NULL, TR, RADIANS(alpha), TE, nfaf, faf_coefs) ;
    else
      mri_out = MRIsynthesize(mri_T1, mri_PD, mri_T2star, NULL, TR, RADIANS(alpha), TE) ;
  }

  if (nbias > 0)
    apply_bias_field(mri_out, nbias, bias_coefs) ;
  printf("writing output to %s.\n", out_fname) ;
  MRIaddCommandLine(mri_out, cmdline) ;
  MRIwrite(mri_out, out_fname) ;

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
  else if (!stricmp(option, "remap")) {
    nl_remap_T1 = 1 ;
    nl_mean = atof(argv[2]) ;
    nl_scale = atof(argv[3]) ;
    printf("remapping T1 with %2.0f * (tanh(%2.4f * (T1-%2.0f))+1.5)\n",
           nl_mean, nl_scale, nl_mean) ;
    nargs = 2 ;
  } else if (!stricmp(option, "PDsat")) {
    PDsat = atof(argv[2]) ;
    nargs = 1 ;
    printf("saturating PD with tanh(PD/(%2.1f*(max-min)))\n", PDsat) ;
  } else if (!stricmp(option, "T2") || !stricmp(option, "T2star")) {
    T2star_fname = argv[2] ;
    nargs = 1 ;
    printf("using T2* volume %s for synthesis\n", T2star_fname) ;
  } else if (!stricmp(option, "bias")) {
    int i, j ;

    nbias = atoi(argv[2]) ;
    nargs = 2*3*nbias+1 ;
    for (i = 0 ; i < 3 ; i++)
      for (j = 0 ; j < 2 ; j++) {
        bias_coefs[i][j] = (float *)calloc(nbias, sizeof(float)) ;
        if (!bias_coefs[i][j])
          ErrorExit(ERROR_NOMEMORY, "%s: could not allocate bias coefficient array (%d, %d)", Progname,i, j) ;
      }

    for (j = 3, i = 0 ; i < nbias ; i++) {
      bias_coefs[0][0][i] = atof(argv[j++]) ;
      bias_coefs[0][1][i] = atof(argv[j++]) ;
      printf("%2.3f cos(%d w0 x) + %2.3f sin(%d w0 x)\n",
             bias_coefs[0][0][i], i+1, bias_coefs[0][1][i], i+1) ;
    }
    for (i = 0 ; i < nbias ; i++) {
      bias_coefs[1][0][i] = atof(argv[j++]) ;
      bias_coefs[1][1][i] = atof(argv[j++]) ;
      printf("%2.3f cos(%d w0 y) + %2.3f sin(%d w0 y)\n",
             bias_coefs[1][0][i], i+1, bias_coefs[1][1][i], i+1) ;
    }
    for (i = 0 ; i < nbias ; i++) {
      bias_coefs[2][0][i] = atof(argv[j++]) ;
      bias_coefs[2][1][i] = atof(argv[j++]) ;
      printf("%2.3f cos(%d w0 z) + %2.3f sin(%d w0 z)\n",
             bias_coefs[2][0][i], i+1, bias_coefs[2][1][i], i+1) ;
    }
  } else if (!stricmp(option, "faf")) {
    int i, j ;

    nfaf = atoi(argv[2]) ;
    nargs = 2*3*nfaf+1 ;
    for (i = 0 ; i < 3 ; i++)
      for (j = 0 ; j < 2 ; j++) {
        faf_coefs[i][j] = (float *)calloc(nfaf, sizeof(float)) ;
        if (!faf_coefs[i][j])
          ErrorExit(ERROR_NOMEMORY, "%s: could not allocate faf coefficient array (%d, %d)", Progname,i, j) ;
      }

    for (j = 3, i = 0 ; i < nfaf ; i++) {
      faf_coefs[0][0][i] = atof(argv[j++]) ;
      faf_coefs[0][1][i] = atof(argv[j++]) ;
      printf("%2.3f cos(%d w0 x) + %2.3f sin(%d w0 x)\n",
             faf_coefs[0][0][i], i+1, faf_coefs[0][1][i], i+1) ;
    }
    for (i = 0 ; i < nfaf ; i++) {
      faf_coefs[1][0][i] = atof(argv[j++]) ;
      faf_coefs[1][1][i] = atof(argv[j++]) ;
      printf("%2.3f cos(%d w0 y) + %2.3f sin(%d w0 y)\n",
             faf_coefs[1][0][i], i+1, faf_coefs[1][1][i], i+1) ;
    }
    for (i = 0 ; i < nfaf ; i++) {
      faf_coefs[2][0][i] = atof(argv[j++]) ;
      faf_coefs[2][1][i] = atof(argv[j++]) ;
      printf("%2.3f cos(%d w0 z) + %2.3f sin(%d w0 z)\n",
             faf_coefs[2][0][i], i+1, faf_coefs[2][1][i], i+1) ;
    }
  }
#if 0
  else if (!stricmp(option, "faf")) {
    int i ;

    nfaf = atoi(argv[2]) ;
    nargs = 3*nfaf+1 ;
    faf_coefs[0][0] = (float *)calloc(nfaf, sizeof(float)) ;
    faf_coefs[0][1] = (float *)calloc(nfaf, sizeof(float)) ;
    faf_coefs[1][0] = (float *)calloc(nfaf, sizeof(float)) ;
    faf_coefs[1][1] = (float *)calloc(nfaf, sizeof(float)) ;
    faf_coefs[1][2] = (float *)calloc(nfaf, sizeof(float)) ;
    faf_coefs[2][2] = (float *)calloc(nfaf, sizeof(float)) ;
    if (!faf_coefs[0][0] || !faf_coefs[1][0] || !faf_coefs[2][0] ||
        !faf_coefs[0][1] || !faf_coefs[1][1] || !faf_coefs[2][1])
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate faf coefficient array", Progname) ;
    for (i = 0 ; i < nfaf ; i++) {
      faf_coefs[0][i] = atof(argv[3+3*i]) ;
      faf_coefs[1][i] = atof(argv[3+3*i+1]) ;
      faf_coefs[2][i] = atof(argv[3+3*i+2]) ;
      printf("%2.3f cos(%d w0 x), %2.3f cos(%d w0 y), %2.3f cos(%d w0 z)\n",
             faf_coefs[0][i], i+1, faf_coefs[1][i], i+1, faf_coefs[2][i], i+1) ;
    }
  }
#endif
  else if (!stricmp(option, "jpdf")) {
    jpdf_name = argv[2] ;
    nargs = 1 ;
    printf("remapping T1 values using joint pdf file %s...\n", jpdf_name) ;
  } else if (!stricmp(option, "ijpdf")) {
    jpdf_name = argv[2] ;
    nargs = 1 ;
    invert = 1 ;
    printf("remapping T1 values using inverse of joint pdf file %s...\n", jpdf_name) ;
  } else if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else if (!stricmp(option, "w5")) {
    w5 = atof(argv[2]) ;
    printf("setting 5 degree weight to %f\n", w5) ;
    nargs = 1 ;
  } else if (!stricmp(option, "w30")) {
    w30 = atof(argv[2]) ;
    printf("setting 30 degree weight to %f\n", w30) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'X':
      printf("extracting middle half of images\n") ;
      extract = 1 ;
      break ;
    case 'W':
      use_weighting = 1 ;
      break ;
    case 'D':
      discard = 1 ;
      printf("setting all PD values to constant...\n") ;
      break ;
    case 'N':
      normalize = 1 ;
      printf("normalizing PD before synthesizing image...\n") ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      print_usage();
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void) {
  print_usage() ;
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  printf("usage: %s [options] <TR> <alpha (deg)> <TE> <T1 volume> <PD volume> <output volume>\n",
         Progname) ;
  printf("The -w switch will use a fixed weighting in order to generate an output volume with\n"
         "optimal gray/white contrast\n") ;

}

static void
print_help(void) {
  fprintf(stderr,
          "\nThis program will synthesize a flash acquisition based on previously computed T1/PD maps\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

#if 0
static double
FLASHforwardModel(double flip_angle, double TR, double PD, double T1) {
  double  CFA = 1, SFA = 0 ;
  double  E1, FLASH ;

  CFA = cos(flip_angle) ;
  SFA = sin(flip_angle) ;

  E1 = exp(-TR/T1) ;

  FLASH = PD * SFA ;
  if (!DZERO(T1))
    FLASH *= (1-E1)/(1-CFA*E1);
  return(FLASH) ;
}
#endif
static MRI *
MRIsynthesize(MRI *mri_T1, MRI *mri_PD, MRI *mri_T2star, MRI *mri_dst, double TR, double alpha, double TE) {
  int   x, y, z, width, height, depth ;
  double flash, T1, PD ;


  if (!mri_dst)
    mri_dst = MRIclone(mri_T1, NULL) ;

  mri_dst->tr = TR ;
  mri_dst->flip_angle = alpha ;
  mri_dst->te = TE ;
  mri_dst->ti = 0 ;
  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;
  for (x = 0 ; x < width ; x++) {
    for (y = 0 ; y < height ; y++) {
      for (z = 0 ; z < depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        MRIsampleVolume(mri_T1, x, y, z, &T1) ;
        if (T1 <= 0)
          T1 = 1 ;
        if (T1 < 900 && T1 > 600)
          DiagBreak() ;
        MRIsampleVolume(mri_PD, x, y, z, &PD) ;
        if (mri_T2star) {
          double T2star ;
          MRIsampleVolume(mri_T2star, x, y, z, &T2star) ;
          flash = FLASHforwardModelT2star(T1, PD, T2star, TR, alpha, TE) ;
        } else
          flash = FLASHforwardModel(T1, PD, TR, alpha, TE) ;
        MRIsetVoxVal(mri_dst, x, y, z, 0, flash) ;
        if (!std::isfinite(flash))
          DiagBreak() ;
      }
    }
  }

  return(mri_dst) ;
}

static int
normalize_PD(MRI *mri_PD, float target) {
  double mean_PD, scale, val ;
  int    x, y, z ;

  for (mean_PD = 0.0, x = 0 ; x < mri_PD->width ; x++) {
    for (y = 0 ; y < mri_PD->height ; y++) {
      for (z = 0 ; z < mri_PD->depth ; z++) {
        mean_PD += (double)MRIgetVoxVal(mri_PD, x, y, z,0) ;
      }
    }
  }
  mean_PD /= (mri_PD->width * mri_PD->height * mri_PD->depth) ;
  scale = target / mean_PD ;
  printf("mean PD %2.0f, scaling by %2.2f to set mean to %2.0f\n",
         mean_PD, scale, target) ;
  for (mean_PD = 0.0, x = 0 ; x < mri_PD->width ; x++) {
    for (y = 0 ; y < mri_PD->height ; y++) {
      for (z = 0 ; z < mri_PD->depth ; z++) {
        val = (double)MRIgetVoxVal(mri_PD, x, y, z,0) ;
        val *= scale ;
        MRIsetVoxVal(mri_PD, x, y, z,0, val);
      }
    }
  }
  return(NO_ERROR) ;
}

static int
discard_PD(MRI *mri_PD, short thresh, short target) {
  int    x, y, z ;
  double  val ;

  for (x = 0 ; x < mri_PD->width ; x++) {
    for (y = 0 ; y < mri_PD->height ; y++) {
      for (z = 0 ; z < mri_PD->depth ; z++) {
        val = MRIgetVoxVal(mri_PD, x, y, z,0) ;
        if (val > thresh)
          MRIsetVoxVal(mri_PD, x, y, z,0, target);
        else
          MRIsetVoxVal(mri_PD, x, y, z,0,0);
      }
    }
  }
  return(NO_ERROR) ;
}
static int
remap_T1(MRI *mri_T1, float mean, float scale) {
  int    x, y, z ;
  float val ;

  for (x = 0 ; x < mri_T1->width ; x++) {
    for (y = 0 ; y < mri_T1->height ; y++) {
      for (z = 0 ; z < mri_T1->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val = (float)MRIgetVoxVal(mri_T1, x, y, z,0) ;
        val = mean * (tanh(scale * (val-mean))+1.5) ;
        MRIsetVoxVal(mri_T1, x, y, z,0, val) ;
      }
    }
  }
  return(NO_ERROR) ;
}
#define MIN_REAL_VAL 100  /* with PD scaled to be 1000 */
#define MIN_BIN  50
MRI *
MRIsynthesizeWeightedVolume(MRI *mri_T1, MRI *mri_PD, float w5, float TR5,
                            float w30, float TR30, float target_wm, float TE) {
  MRI        *mri_dst ;
  int        x, y, z, width, height, depth ;
  MRI        *mri30, *mri5 ;
  double     val30, val5, val, min_val ;
#if 0
  int        mri_peak, n, min_real_bin ;
  double    mean_PD ;
  MRI_REGION box ;
  HISTOGRAM *h_mri, *h_smooth ;
  float      x0, y0, z0, min_real_val ;
#endif

  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;
  mri_dst = MRIalloc(width, height, depth, MRI_FLOAT) ;
  MRIcopyHeader(mri_T1, mri_dst) ;
  mri30 = MRIsynthesize(mri_T1, mri_PD, NULL, NULL, TR30, RADIANS(30), TE) ;
  mri5 = MRIsynthesize(mri_T1, mri_PD, NULL, NULL, TR5, RADIANS(5), TE) ;
#if 0
  mean_PD = MRImeanFrame(mri_PD, 0) ;
  /*  MRIscalarMul(mri_PD, mri_PD, 1000.0f/mean_PD) ;*/


  h_mri = MRIhistogram(mri30, 100) ;
  h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
  mri_peak = HISTOfindHighestPeakInRegion(h_smooth, 0, h_smooth->nbins) ;
  min_real_bin = HISTOfindNextValley(h_smooth, mri_peak) ;
  min_real_val = h_smooth->bins[min_real_bin] ;

  MRIfindApproximateSkullBoundingBox(mri30, min_real_val, &box) ;
  x0 = box.x+box.dx/3 ;
  y0 = box.y+box.dy/3 ;
  z0 = box.z+box.dz/2 ;
  printf("using (%.0f, %.0f, %.0f) as brain centroid...\n",x0, y0, z0) ;
  box.dx /= 4 ;
  box.x = x0 - box.dx/2;
  box.dy /= 4 ;
  box.y = y0 - box.dy/2;
  box.dz /= 4 ;
  box.z = z0 - box.dz/2;


  printf("using box (%d,%d,%d) --> (%d, %d,%d) "
         "to find MRI wm\n", box.x, box.y, box.z,
         box.x+box.dx-1,box.y+box.dy-1, box.z+box.dz-1) ;

  h_mri = MRIhistogramRegion(mri30, 0, NULL, &box) ;
  for (n = 0 ; n < h_mri->nbins-1 ; n++)
    if (h_mri->bins[n+1] > min_real_val)
      break ;
  HISTOclearBins(h_mri, h_mri, 0, n) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    HISTOplot(h_mri, "mri.histo") ;
  mri_peak = HISTOfindLastPeak(h_mri, HISTO_WINDOW_SIZE,MIN_HISTO_PCT);
  mri_peak = h_mri->bins[mri_peak] ;
  printf("before smoothing, mri peak at %d\n", mri_peak) ;
  h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    HISTOplot(h_smooth, "mri_smooth.histo") ;
  mri_peak = HISTOfindLastPeak(h_smooth, HISTO_WINDOW_SIZE,MIN_HISTO_PCT);
  mri_peak = h_mri->bins[mri_peak] ;
  printf("after smoothing, mri peak at %d\n", mri_peak) ;
  HISTOfree(&h_smooth) ;
  HISTOfree(&h_mri) ;
#endif

  min_val = 0 ;
  for (x = 0 ; x < width ; x++) {
    for (y = 0 ; y < height ; y++) {
      for (z = 0 ; z < depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        MRIsampleVolumeType(mri30, x, y, z, &val30, SAMPLE_NEAREST) ;
        MRIsampleVolumeType(mri5, x, y, z, &val5, SAMPLE_NEAREST) ;
        val = w30*val30 + w5*val5 ;
        MRIFvox(mri_dst, x, y, z) = val ;
        if (val < min_val)
          min_val = val ;
      }
    }
  }

  for (x = 0 ; x < width ; x++) {
    for (y = 0 ; y < height ; y++) {
      for (z = 0 ; z < depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        MRIFvox(mri_dst, x, y, z) += min_val ;
      }
    }
  }

  MRIfree(&mri30) ;
  MRIfree(&mri5) ;
  return(mri_dst) ;
}
static int
transform_T1_values_using_joint_pdf(MRI *mri_T1, char *jpdf_name, int invert) {
  int   x, y, z, nbins, i, j, **jpdf, max_j, max_count ;
  double  T1 ;
  FILE  *fp ;
  float fstep, fmin, fmax, val ;
  char  line[STRLEN], *cp, var_name[STRLEN] ;

  fstep = 10 ;
  fmin = 10;
  fmax = 5000 ;
  nbins = 500 ;
  fp = fopen(jpdf_name, "r") ;
  if (fp == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read joint pdf file %s\n", Progname, jpdf_name) ;


  while ((cp = fgetl(line, STRLEN-1, fp)) != NULL) {
    sscanf(cp, "%s = %f", var_name, &val) ;
    if (stricmp(var_name, "nbins") == 0)
      nbins = nint(val) ;
    else if (stricmp(var_name, "fmin") == 0)
      fmin = val ;
    else if (stricmp(var_name, "fmax") == 0)
      fmax = val ;
    else if (stricmp(var_name, "fstep") == 0)
      fstep = val ;
    else if (stricmp(var_name, "joint_density") == 0)
      break ;
    else
      break ;
  }

  printf("using nbins = %d, fmin = %2.1f, fmax = %2.1f, fstep = %2.1f\n", nbins, fmin, fmax, fstep) ;
  jpdf = (int **)calloc(nbins, sizeof(int *)) ;
  if (!jpdf)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d x %d bin jpdf lookup table\n", Progname, nbins, nbins);
  for (i = 0 ; i < nbins ; i++) {
    jpdf[i] = (int *)calloc(nbins, sizeof(int)) ;
    if (!jpdf[i])
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d x %d bin jpdf lookup table\n", Progname, nbins, nbins);
    for (j = 0; j < nbins ; j++) {
      if (fscanf(fp, "%d", &jpdf[i][j]) != 1)
        ErrorExit(ERROR_BADFILE, "%s: could not scan element %d, %d from %s",
                  Progname, i, j, jpdf_name) ;
    }
    fscanf(fp, " ;\n") ;
  }
  fclose(fp) ;

  for (x = 0 ; x < mri_T1->width ; x++) {
    for (y = 0 ; y < mri_T1->height ; y++) {
      for (z = 0 ;  z < mri_T1->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak()  ;
        MRIsampleVolume(mri_T1, x, y, z, &T1) ;
        if (T1 <= 0)
          T1 = 1 ;
        if (invert) {
          i = nint((T1-fmin)/fstep) ;
          if (i >= nbins)
            i = nbins-1 ;
          if (i < 0)
            i = 0 ;
          max_count = jpdf[i][max_j=0] ;
          for (j = 1 ; j < nbins ; j++) {
            if (jpdf[i][j] > max_count) {
              max_count = jpdf[i][j] ;
              max_j = j ;
            }
          }
        } else {
          i = nint((T1-fmin)/fstep) ;
          if (i >= nbins)
            i = nbins-1 ;
          if (i < 0)
            i = 0 ;
          max_count = jpdf[max_j=0][i] ;
          for (j = 1 ; j < nbins ; j++) {
            if (jpdf[i][j] > max_count) {
              max_count = jpdf[j][i] ;
              max_j = j ;
            }
          }
        }

        T1 = max_j * fstep + fmin ;
        MRIsetVoxVal(mri_T1, x, y, z, 0, T1) ;
      }
    }
  }

  return(NO_ERROR) ;
}

static int
apply_bias_field(MRI *mri, int nbias, float *bias_coefs[3][2]) {
  int    x, y, z, n ;
  double xb, yb, zb, x0, y0, z0, w0x, w0y, w0z ;
  float  val ;

  x0 = mri->width/2 ;
  y0 = mri->height/2 ;
  z0 = mri->depth/2 ;
  w0x = 2/x0 ;
  w0y = 2/y0 ;
  w0z = 2/z0 ;
  for (x = 0 ; x < mri->width ; x++) {
    for (xb = 1.0, n=1 ; n <= nbias ; n++)
      xb += bias_coefs[0][0][n-1] * cos(w0x*n*(x-x0)) + bias_coefs[0][1][n-1] * sin(w0x*n*(x-x0)) ;
    for (y = 0 ; y < mri->height ; y++) {
      for (yb = 1.0, n=1 ; n <= nbias ; n++)
        yb += bias_coefs[1][0][n-1] * cos(w0y*n*(y-y0)) + bias_coefs[1][1][n-1] * sin(w0y*n*(y-y0)) ;
      for (z = 0 ; z < mri->depth ; z++) {
        for (zb = 1.0, n=1 ; n <= nbias ; n++)
          zb += bias_coefs[2][0][n-1] * cos(w0z*n*(z-z0)) + bias_coefs[2][1][n-1] * sin(w0z*n*(z-z0)) ;
        val = MRIgetVoxVal(mri, x, y, z, 0) ;
        val = val * xb * yb * zb ;
        MRIsetVoxVal(mri, x, y, z, 0, val) ;
      }
    }
  }

  return(NO_ERROR) ;
}

static MRI *
MRIsynthesizeWithFAF(MRI *mri_T1, MRI *mri_PD, MRI *mri_dst, double TR, double alpha, double TE, int nfaf,
                     float *faf_coefs[3][2]) {
  int   x, y, z, width, height, depth, n ;
  double flash, T1, PD ;
  double xb, yb, zb, x0, y0, z0, w0x, w0y, w0z ;

  x0 = mri_T1->width/2 ;
  y0 = mri_T1->height/2 ;
  z0 = mri_PD->depth/2 ;
  w0x = 2/x0 ;
  w0y = 2/y0 ;
  w0z = 2/z0 ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_T1, NULL) ;

  mri_dst->tr = TR ;
  mri_dst->flip_angle = alpha ;
  mri_dst->te = TE ;
  mri_dst->ti = 0 ;
  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;
  for (x = 0 ; x < width ; x++) {
    for (xb = 1.0, n=1 ; n <= nfaf ; n++)
      xb += faf_coefs[0][0][n-1] * cos(w0x*n*(x-x0)) + faf_coefs[0][1][n-1] * sin(w0x*n*(x-x0)) ;
    for (y = 0 ; y < height ; y++) {
      for (yb = 1.0, n=1 ; n <= nfaf ; n++)
        yb += faf_coefs[1][0][n-1] * cos(w0y*n*(y-y0)) + faf_coefs[1][1][n-1] * sin(w0y*n*(y-y0)) ;
      for (z = 0 ; z < depth ; z++) {
        for (zb = 1.0, n=1 ; n <= nfaf ; n++)
          zb += faf_coefs[2][0][n-1] * cos(w0z*n*(z-z0)) + faf_coefs[2][1][n-1] * sin(w0z*n*(z-z0)) ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        MRIsampleVolume(mri_T1, x, y, z, &T1) ;
        if (T1 <= 0)
          T1 = 1 ;
        if (T1 < 900 && T1 > 600)
          DiagBreak() ;
        MRIsampleVolume(mri_PD, x, y, z, &PD) ;
        flash = FLASHforwardModel(T1, PD, TR, xb*yb*zb*alpha, TE) ;
        MRIsetVoxVal(mri_dst, x, y, z, 0, flash) ;
      }
    }
  }

  return(mri_dst) ;
}
static int
saturate_PD(MRI *mri, float PDsat) {
  int  x,y , z ;
  float mx, mn, val ;

  printf("saturating PD (%2.3f)\n", PDsat) ;
  MRIvalRange(mri, &mn, &mx) ;

  for (x = 0 ; x < mri->width ; x++) {
    for (y = 0 ; y < mri->height ; y++) {
      for (z = 0 ; z < mri->depth ; z++) {
        val = MRIgetVoxVal(mri, x, y, z, 0) ;
        val = 1000 * tanh(val / (PDsat * (mx-mn))) ;
        MRIsetVoxVal(mri, x, y, z, 0, val) ;
      }
    }
  }

  return(NO_ERROR) ;
}

