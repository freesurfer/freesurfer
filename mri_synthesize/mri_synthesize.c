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

static char vcid[] = "$Id: mri_synthesize.c,v 1.6 2003/04/16 18:00:03 kteich Exp $";

int main(int argc, char *argv[]) ;

static int normalize_PD(MRI *mri_PD, float target) ;
static int discard_PD(MRI *mri_PD, short thresh, short target) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static double FLASHforwardModel(double flip_angle, double TR, double PD, 
                                double T1) ;

MRI *MRIsynthesizeWeightedVolume(MRI *mri_T1, MRI *mri_PD, float w5, float TR5,
                                 float w30, float TR30, float target_wm, float TE) ;
static MRI *MRIsynthesize(MRI *mri_T1, MRI *mri_PD, MRI *mri_dst, double TR, double alpha, double TE) ;
static int remap_T1(MRI *mri_T1, float mean, float scale) ;

char *Progname ;
static int normalize = 0 ;
static int discard = 0 ;
static int nl_remap_T1 = 0 ;
static float nl_scale = 0.01 ;
static float nl_mean = 950 ;

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
main(int argc, char *argv[])
{
  char        **av, *out_fname, *T1_fname, *PD_fname ;
  int         ac, nargs ;
  MRI         *mri_T1, *mri_PD, *mri_out ;
  float       TR, TE, alpha ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_synthesize.c,v 1.6 2003/04/16 18:00:03 kteich Exp $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)
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

  printf("reading PD volume from %s...\n", PD_fname) ;
  mri_PD = MRIread(PD_fname) ;
  if (!mri_PD)
    ErrorExit(ERROR_NOFILE, "%s: could not read PD volume %s", Progname, 
              PD_fname) ;

  if (use_weighting)
  {
    mri_out = MRIsynthesizeWeightedVolume(mri_T1, mri_PD, w5, TR, w30, TR, 110,TE);
  }
  else
  {
    printf("synthesizing volume with TR=%2.1f msec, TE=%2.1f msec, and alpha=%2.2f degrees...\n",
           TR, TE, alpha) ;
    if (normalize)
      normalize_PD(mri_PD, 1000) ;
    if (discard)
      discard_PD(mri_PD, 250, 1500) ;
    if (nl_remap_T1)
      remap_T1(mri_T1, nl_mean, nl_scale) ;
    mri_out = MRIsynthesize(mri_T1, mri_PD, NULL, TR, RADIANS(alpha), TE) ;
  }

  printf("writing output to %s.\n", out_fname) ;
  MRIwrite(mri_out, out_fname) ;

  exit(0) ;
  return(0) ;  /* for ansi */
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "remap"))
  {
    nl_remap_T1 = 1 ;
    nl_mean = atof(argv[2]) ;
    nl_scale = atof(argv[3]) ;
    printf("remapping T1 with %2.0f * (tanh(%2.4f * (T1-%2.0f))+1.5)\n",
           nl_mean, nl_scale, nl_mean) ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "w5"))
  {
    w5 = atof(argv[2]) ;
    printf("setting 5 degree weight to %f\n", w5) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "w30"))
  {
    w30 = atof(argv[2]) ;
    printf("setting 30 degree weight to %f\n", w30) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
  {
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
usage_exit(void)
{
  print_usage() ;
  print_help() ;
  exit(1) ;
}

static void
print_usage(void)
{
  printf("usage: %s [options] <TR> <alpha (deg)> <TE> <T1 volume> <PD volume> <output volume>\n",
				 Progname) ;
	printf("The -w switch will use a fixed weighting in order to generate an output volume with\n"
				 "optimal gray/white contrast\n") ;
	
}

static void
print_help(void)
{
  fprintf(stderr, 
          "\nThis program will synthesize a flash acquisition based on previously computed T1/PD maps\n");
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static double
FLASHforwardModel(double flip_angle, double TR, double PD, double T1)
{
  double  CFA = 1, SFA = 0 ;
  double  E1, FLASH ;

  CFA = cos(flip_angle) ; SFA = sin(flip_angle) ;

  E1 = exp(-TR/T1) ;
      
  FLASH = PD * SFA ;
  if (!DZERO(T1))
    FLASH *= (1-E1)/(1-CFA*E1);
  return(FLASH) ;
}

static MRI *
MRIsynthesize(MRI *mri_T1, MRI *mri_PD, MRI *mri_dst, double TR, double alpha, double TE)
{
  int   x, y, z, width, height, depth ;
  double flash, T1, PD ;
  

  if (!mri_dst)
    mri_dst = MRIclone(mri_T1, NULL) ;

	mri_dst->tr = TR ; mri_dst->flip_angle = alpha ; mri_dst->te = TE ; mri_dst->ti = 0 ;
  width = mri_T1->width ; height = mri_T1->height ; depth = mri_T1->depth ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        T1 = MRISvox(mri_T1, x, y, z) ;
        if (T1 < 900 && T1 > 600)
          DiagBreak() ;
        PD = MRISvox(mri_PD, x, y, z) ;
        flash = FLASHforwardModel(alpha, TR, PD, T1) ;
        MRISvox(mri_dst, x, y, z) = (short)nint(flash) ;
      }
    }
  }

  return(mri_dst) ;
}

static int
normalize_PD(MRI *mri_PD, float target)
{
  double mean_PD, scale, val ;
  int    x, y, z ;

  for (mean_PD = 0.0, x = 0 ; x < mri_PD->width ; x++)
  {
    for (y = 0 ; y < mri_PD->height ; y++)
    {
      for (z = 0 ; z < mri_PD->depth ; z++)
      {
        mean_PD += (double)MRISvox(mri_PD, x, y, z) ;
      }
    }
  }
  mean_PD /= (mri_PD->width * mri_PD->height * mri_PD->depth) ;
  scale = target / mean_PD ;
  printf("mean PD %2.0f, scaling by %2.2f to set mean to %2.0f\n",
         mean_PD, scale, target) ;
  for (mean_PD = 0.0, x = 0 ; x < mri_PD->width ; x++)
  {
    for (y = 0 ; y < mri_PD->height ; y++)
    {
      for (z = 0 ; z < mri_PD->depth ; z++)
      {
        val = (double)MRISvox(mri_PD, x, y, z) ;
        val *= scale ;
        MRISvox(mri_PD, x, y, z) = (short)val ;
      }
    }
  }
  return(NO_ERROR) ;
}

static int
discard_PD(MRI *mri_PD, short thresh, short target)
{
  int    x, y, z ;
  short  val ;

  for (x = 0 ; x < mri_PD->width ; x++)
  {
    for (y = 0 ; y < mri_PD->height ; y++)
    {
      for (z = 0 ; z < mri_PD->depth ; z++)
      {
        val = MRISvox(mri_PD, x, y, z) ;
        if (val > thresh)
          MRISvox(mri_PD, x, y, z) = target ;
        else
          MRISvox(mri_PD, x, y, z) = 0 ;
      }
    }
  }
  return(NO_ERROR) ;
}
static int
remap_T1(MRI *mri_T1, float mean, float scale)
{
  int    x, y, z ;
  float val ;

  for (x = 0 ; x < mri_T1->width ; x++)
  {
    for (y = 0 ; y < mri_T1->height ; y++)
    {
      for (z = 0 ; z < mri_T1->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val = (float)MRISvox(mri_T1, x, y, z) ;
        val = mean * (tanh(scale * (val-mean))+1.5) ;
        MRISvox(mri_T1, x, y, z) = val ;
      }
    }
  }
  return(NO_ERROR) ;
}
#define MIN_REAL_VAL 100  /* with PD scaled to be 1000 */
#define MIN_BIN  50
MRI *
MRIsynthesizeWeightedVolume(MRI *mri_T1, MRI *mri_PD, float w5, float TR5,
                            float w30, float TR30, float target_wm, float TE)
{
  MRI *mri_dst ;
  MRI_REGION box ;
  float      x0, y0, z0, min_real_val ;
  HISTOGRAM *h_mri, *h_smooth ;
  int        mri_peak, n, min_real_bin, x, y, z, width, height, depth ;
  MRI       *mri30, *mri5 ;
  double    mean_PD ;
  Real      val30, val5, val, min_val ;

  mean_PD = MRImeanFrame(mri_PD, 0) ;
  /*  MRIscalarMul(mri_PD, mri_PD, 1000.0f/mean_PD) ;*/
  mri30 = MRIsynthesize(mri_T1, mri_PD, NULL, TR30, RADIANS(30), TE) ;
  mri5 = MRIsynthesize(mri_T1, mri_PD, NULL, TR5, RADIANS(5), TE) ;
  width = mri30->width ; height = mri30->height ; depth = mri30->depth ;

  mri_dst = MRIalloc(width, height, depth, MRI_FLOAT) ;
  MRIcopyHeader(mri_T1, mri_dst) ;

  h_mri = MRIhistogram(mri30, 100) ;
  h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
  mri_peak = HISTOfindHighestPeakInRegion(h_smooth, 0, h_smooth->nbins) ;
  min_real_bin = HISTOfindNextValley(h_smooth, mri_peak) ;
  min_real_val = h_smooth->bins[min_real_bin] ;

  MRIfindApproximateSkullBoundingBox(mri30, min_real_val, &box) ;
  x0 = box.x+box.dx/3 ; y0 = box.y+box.dy/3 ; z0 = box.z+box.dz/2 ;
  printf("using (%.0f, %.0f, %.0f) as brain centroid...\n",x0, y0, z0) ;
  box.dx /= 4 ; box.x = x0 - box.dx/2;
  box.dy /= 4 ; box.y = y0 - box.dy/2;
  box.dz /= 4 ; box.z = z0 - box.dz/2;


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
  HISTOfree(&h_smooth) ; HISTOfree(&h_mri) ;

  min_val = 0 ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
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

  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        MRIFvox(mri_dst, x, y, z) += min_val ;
      }
    }
  }

  MRIfree(&mri30) ; MRIfree(&mri5) ;
  return(mri_dst) ;
}
