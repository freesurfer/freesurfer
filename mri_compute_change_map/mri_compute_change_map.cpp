/**
 * @brief computes longitudinal change map for a pair of (registered) images
 *
 * Program for computing a change map for a pair of registered images by
 * using a robust estimate of the noise.
 */
/*
 * Original Author: Bruce Fischl
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
#include "version.h"
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;


const char *Progname ;


#define MAX_LOGP 1000
#define MAX_SIMS 1000

static void usage_exit(int code) ;

static char *log_fname = NULL ;
static float smooth_sigma = 0.0 ;
static int mean_filter = 0 ;
static int erode = 0 ;
static int nbhd_size = 0 ;
static int bonferroni = 0 ;

MRI *MRIcomputeChangeMap(MRI *mri1, MRI *mri2, TRANSFORM *transform, MRI *mri_change, float *pstd) ;
MRI *MRIcomputeNbhdPvalues(MRI *mri_src, int nbhd_size, MRI *mri_dst) ;
MRI *MRIbonferroniCorrect(MRI *mri_src, MRI *mri_dst, int is_log) ; 

int MRImaskRegisteredPair(MRI *mri1, MRI *mri2, TRANSFORM *transform, MRI *mri_mask) ;

static MRI *mri_mask = NULL ;
int
main(int argc, char *argv[]) {
  char   *out_fname, **av ;
  int    ac, nargs ;
  MRI    *mri1, *mri2, *mri_change ;
  int     msec, minutes, seconds, i ;
  Timer start ;
  TRANSFORM     *transform ;
  float         std ;

  nargs = handleVersionOption(argc, argv, "mri_compute_change_map");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;

  mri1 = MRIread(argv[1]) ;
  if (mri1 == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read input volume %s\n",
              Progname,argv[1]);
  for (i = 0 ; i < erode ; i++)
    MRIerodeZero(mri1, mri1) ;
  mri2 = MRIread(argv[2]) ;
  if (mri2 == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read input volume %s\n",
              Progname,argv[2]);
  for (i = 0 ; i < erode ; i++)
    MRIerodeZero(mri2, mri2) ;
  transform = TransformRead(argv[3]) ;
  if (transform == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read transform %s\n",
              Progname,argv[3]);
	// if identity, set geometries (to allow conversion later
	// from vox2vox to ras2ras etc) also set vox2vox as that
	// is needed below (avoid conversion)
  if (0 == strcmp(argv[3], "identity.nofile"))
  {
	  transform->type =  LINEAR_VOX_TO_VOX;
		LTA* lta = (LTA *)transform->xform;
		lta->type =  LINEAR_VOX_TO_VOX;
	  getVolGeom(mri2,&lta->xforms[0].src);
	  getVolGeom(mri1,&lta->xforms[0].dst);
		//LTAwrite(lta,"identity.lta");
	}

  if (smooth_sigma > 0)
  {
    MRI *mri_smooth, *mri_kernel = MRIgaussian1d(smooth_sigma, 100), *mri_mask ;
    int i, erodes ;

    erodes = MAX((int)ceil(smooth_sigma*2),1) ;
    printf("eroding input images %d times to remove border effects\n", erodes) ;

    mri_smooth = MRIconvolveGaussian(mri1, NULL, mri_kernel) ;
    MRIbinarize(mri1, mri1, 1, 0, 1) ;
    for (i = 0 ; i < erodes ; i++)
    {
      mri_mask = MRIerodeZero(mri1, NULL) ;
      MRIcopy(mri_mask, mri1) ;
      MRIfree(&mri_mask) ;
    }
    MRImask(mri_smooth, mri1, mri_smooth, 0, 0) ;
    MRIfree(&mri1) ; mri1 = mri_smooth ;


    mri_smooth = MRIconvolveGaussian(mri2, NULL, mri_kernel) ;
    MRIbinarize(mri2, mri2, 1, 0, 1) ;
    for (i = 0 ; i < erodes ; i++)
    {
      mri_mask = MRIerodeZero(mri2, NULL) ;
      MRIcopy(mri_mask, mri2) ;
      MRIfree(&mri_mask) ;
    }
    MRImask(mri_smooth, mri2, mri_smooth, 0, 0) ;
    MRIfree(&mri2) ; mri2 = mri_smooth ;

    MRIfree(&mri_kernel) ;
  }

  if (mean_filter)
  {
    MRI *mri_smooth, *mri_mask ;
    int i ;
    mri_smooth = MRImean(mri1, NULL, mean_filter) ;
    MRIbinarize(mri1, mri1, 1, 0, 1) ;
    for (i = 0 ; i < (mean_filter-1)/2 ; i++)
    {
      mri_mask = MRIerodeZero(mri1, NULL) ;
      MRIcopy(mri_mask, mri1) ;
      MRIfree(&mri_mask) ;
    }
    MRImask(mri_smooth, mri1, mri_smooth, 0, 0) ;
    MRIfree(&mri1) ; mri1 = mri_smooth ; 
    mri_smooth = MRImean(mri2, NULL, mean_filter) ;
    MRIbinarize(mri2, mri2, 1, 0, 1) ;
    for (i = 0 ; i < (mean_filter-1)/2 ; i++)
    {
      mri_mask = MRIerodeZero(mri2, NULL) ;
      MRIcopy(mri_mask, mri2) ;
      MRIfree(&mri_mask) ;
    }
    MRImask(mri_smooth, mri2, mri_smooth, 0, 0) ;
    MRIfree(&mri2) ; mri2 = mri_smooth ; 
  }

  if (mri_mask)
    MRImaskRegisteredPair(mri1, mri2, transform, mri_mask) ;

  mri_change = MRIcomputeChangeMap(mri1, mri2, transform, NULL, &std) ;

#if 0
  if (smooth_sigma > 0)  // correct for multiple comparisons using GRF and monte carlo sims
  {
    MRI          *mri_kernel = MRIgaussian1d(smooth_sigma, 100), *mri_n1, *mri_n2, 
                 *mri_n1_smoothed, *mri_n2_smoothed, *mri_noise_change;
    int          nsim, b ;
    HISTOGRAM    *h = NULL ;
    float        min_val, max_val ;
    
    mri_n1 = MRIcloneDifferentType(mri1, MRI_FLOAT) ;
    mri_n2 = MRIcloneDifferentType(mri2, MRI_FLOAT) ;
    mri_n1_smoothed = MRIclone(mri_n1, NULL) ;
    mri_n2_smoothed = MRIclone(mri_n2, NULL) ;
    h = HISTOalloc(MAX_LOGP) ; HISTOinit(h, MAX_LOGP, 0, MAX_LOGP) ;
    mri_noise_change = MRIclone(mri_change, NULL) ;
    for (nsim = 0 ; nsim < MAX_SIMS ; nsim++)
    {
      if (!(nsim % 10))
        printf("\rsimulation %d", nsim) ;
      MRIrandn(mri_n1->width, mri_n1->height, mri_n1->depth, mri_n1->nframes, 0.0, std, mri_n1) ;
      MRIrandn(mri_n2->width, mri_n2->height, mri_n2->depth, mri_n2->nframes, 0.0, std, mri_n2) ;
      MRIconvolveGaussian(mri_n1, mri_n1_smoothed, mri_kernel) ;
      MRIconvolveGaussian(mri_n2, mri_n2_smoothed, mri_kernel) ;
      MRIcomputeChangeMap(mri_n1_smoothed, mri_n2_smoothed, transform, mri_noise_change, NULL) ;
      MRIvalRange(mri_noise_change, &min_val, &max_val) ;
      for (b = 0 ; b <= max_val ; b++)
        h->counts[b]++ ;
    }
    printf("\n") ;
    HISTOplot(h, "h.plt") ;

    MRIfree(&mri_kernel) ; MRIfree(&mri_n1) ; MRIfree(&mri_n2) ; MRIfree(&mri_noise_change) ;
    HISTOfree(&h) ;
  }
#endif

#if 0
  if (smooth_sigma > 0)
  {
    MRI *mri_smooth, *mri_kernel = MRIgaussian1d(smooth_sigma, 100) ;
    
    mri_smooth = MRIconvolveGaussian(mri_change, NULL, mri_kernel) ;
    MRIfree(&mri_change) ; mri_change = mri_smooth ;

    MRIfree(&mri_kernel) ;
  }
#endif
  if (bonferroni)
    MRIbonferroniCorrect(mri_change, mri_change, 1) ;

  if (nbhd_size > 0)
  {
    MRI *mri_tmp ;
    mri_tmp = MRIcomputeNbhdPvalues(mri_change, nbhd_size, NULL) ;
    MRIfree(&mri_change) ; mri_change = mri_tmp ;
  }

  out_fname = argv[4] ;
  fprintf(stderr, "writing to %s...\n", out_fname) ;
  if (log_fname)
  {
    FILE *fp ;
    fp = fopen(log_fname, "w") ;
    if (fp == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not optn %s", Progname, log_fname) ;
    fprintf(fp, "%f\n", std) ;
    fclose(fp) ;
  }
  MRIwrite(mri_change, out_fname) ;
  MRIfree(&mri_change) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "morphological processing took %d minutes and %d seconds.\n",
          minutes, seconds) ;
  exit(0) ;
  return(0) ;
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
  if (!stricmp(option, "mask")) 
  {
    mri_mask = MRIread(argv[2]) ;
    if (mri_mask == NULL)
      exit(Gerror) ;
    MRIbinarize(mri_mask, mri_mask, 1, 0, 1) ;
    nargs = 1 ;
  } else if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else switch (toupper(*option)) {
  case 'E':
    erode = atoi(argv[2]) ;
    nargs = 1 ;
    printf("eroding input volumes %d times to reduce effects of skull-strip differences\n", erode) ;
    break ;
  case 'M':
    mean_filter = atoi(argv[2]) ;
    nargs = 1 ;
    printf("smoothing difference image with mean filter w/wsize=%d\n", mean_filter) ;
    break ;
  case 'L':
    log_fname = argv[2] ;
    nargs = 1 ;
    printf("logging results to %s\n", log_fname) ;
    break ;
  case 'S':
    smooth_sigma = atof(argv[2]) ;
    printf("smoothing difference image with Gaussian w/sigma=%2.2f\n", smooth_sigma) ;
    nargs = 1 ;
    break ;
  case 'B':
    bonferroni = 1 ;
    printf("performing Bonferroni correction\n") ;
    break ;
  case 'N':
    nbhd_size = atoi(argv[2]) ;
    printf("computing joint p-value in %d voxel neighborhood\n", nbhd_size) ;
    nargs = 1 ;
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
usage_exit(int code) {
  printf("usage: %s [options] <volume1> <volume2> <transform> <out volume>\n", Progname) ;
  printf("\twhere <transform> should take volume2 coords into volume1 space\n") ;
  printf("\tvalid options are:\n") ;
  printf("\t-m                 : mean filter output before writing\n") ;
  printf("\t-s <sigma>           smooth with Gaussian filter before writing\n") ;
  printf("\tthe output map will be in register with <volume1>\n") ;
  exit(code) ;
}

int
MRImaskRegisteredPair(MRI *mri1, MRI *mri2, TRANSFORM *transform, MRI *mri_mask)
{
  int       x1, y1, z1, x2, y2, z2 ;
  float     xf, yf, zf;
  double    val ;

  for (x1 = 0 ; x1 < mri1->width ; x1++)
    for (y1 = 0 ; y1 < mri1->height ; y1++)
      for (z1 = 0 ; z1 < mri1->depth ; z1++)
      {
        if (x1 == Gx && y1 == Gy && z1 == Gz)
          DiagBreak() ;
        val = MRIgetVoxVal(mri_mask, x1, y1, z1, 0) ;
        if (FZERO(val) == 0)
          continue ;
        MRIsetVoxVal(mri1, x1, y1, z1, 0, 0) ;
      }

  for (x2 = 0 ; x2 < mri2->width ; x2++)
    for (y2 = 0 ; y2 < mri2->height ; y2++)
      for (z2 = 0 ; z2 < mri2->depth ; z2++)
      {
        TransformSampleInverse(transform, (float)x2, (float)y2, (float)z2, &xf, &yf, &zf) ;
        MRIsampleVolume(mri_mask, xf, yf, zf, &val) ;
        if (FZERO(val) == 0)
          continue ;
        MRIsetVoxVal(mri2, x2, y2, z2, 0, 0) ;
      }
  return(NO_ERROR) ;
}
MRI *
MRIcomputeChangeMap(MRI *mri1, MRI *mri2, TRANSFORM *transform, MRI *mri_change, float *pstd)
{
  int       x1, y1, z1, nvox ;
  float     x2, y2, z2, val1, min1, max1, min2, max2, min_change, 
    max_change, dif, mask1 ;
  double    val2, p, logp, mask2, sigma, mean ;
  MRI       *mri_dif, *mri_mean, *mri_used, *mri_mask1, *mri_mask2, *mri_big ;
  HISTOGRAM *h, *hs, *hg ;
  MRI       *mri_xformed ;

  MRIvalRange(mri1, &min1, &max1) ;
  MRIvalRange(mri2, &min2, &max2) ;
  TransformInvert(transform,mri1) ;
  mri_change = MRIcloneDifferentType(mri1, MRI_FLOAT) ;
  mri_dif = MRIcloneDifferentType(mri1, MRI_FLOAT) ;
  mri_mean = MRIcloneDifferentType(mri1, MRI_FLOAT) ;
  mri_used = MRIcloneDifferentType(mri1, MRI_UCHAR) ;
  mri_big = MRIcloneDifferentType(mri1, MRI_UCHAR) ; // for keeping track of underflow values

  mri_mask1 = MRIerodeZero(mri1, NULL) ;
  mri_mask2 = MRIerodeZero(mri2, NULL) ;
  nvox = 0 ;

  mri_xformed = MRIclone(mri1, NULL) ;
  // compute mode of distribution and use it for centering
  for (x1 = 0 ; x1 < mri1->width ; x1++)
    for (y1 = 0 ; y1 < mri1->height ; y1++)
      for (z1 = 0 ; z1 < mri1->depth ; z1++)
      {
        if (x1 == Gx && y1 == Gy && z1 == Gz)
          DiagBreak() ;
        TransformSampleInverse(transform, (float)x1, (float)y1, (float)z1, &x2, &y2, &z2) ;
        val1 = MRIgetVoxVal(mri1, x1, y1, z1, 0) ;
        MRIsampleVolume(mri2, x2, y2, z2, &val2) ;
        mask1 = MRIgetVoxVal(mri_mask1, x1, y1, z1, 0) ;
        MRIsampleVolume(mri_mask2, x2, y2, z2, &mask2) ;
        MRIsetVoxVal(mri_xformed, x1, y1, z1, 0, val2) ;
        if (!FZERO(val1) && !FZERO(val2))
        {
          dif = val2-val1 ;
          MRIsetVoxVal(mri_dif, x1, y1, z1, 0, dif) ;
          MRIsetVoxVal(mri_used, x1, y1, z1, 0, 1) ;
        }
      }

  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] = "xformed.mgz" ;
    printf("writing xformed volume to %s\n",  fname) ;
    MRIwrite(mri_xformed, fname) ;
  }
  MRIfree(&mri_xformed) ;
  MRIvalRange(mri_dif, &min_change, &max_change) ;
  if (max_change-min_change <= 0)
    ErrorExit(ERROR_BADPARM, "%s: input volumes are identical", Progname) ;
  h = MRIhistogramLabel(mri_dif, mri_used, 1, 5*(max_change-min_change)) ;
  HISTOfillHoles(h) ;
  hs = HISTOsmooth(h, NULL, 20) ;
  HISTOrobustGaussianFit(hs, 0.5, &mean, &sigma) ;
  hs = HISTOsmooth(h, NULL, 10) ;
  HISTOrobustGaussianFit(hs, 0.5, &mean, &sigma) ;
  hs = HISTOsmooth(h, NULL, 5) ;
  HISTOrobustGaussianFit(hs, 0.5, &mean, &sigma) ;
  HISTOrobustGaussianFit(h, 0.5, &mean, &sigma) ;
  if (Gdiag & DIAG_WRITE)
  {
    HISTOplot(h, "h.plt") ; 
    HISTOplot(hs, "hs.plt") ; 
  }
  HISTOfree(&h) ; HISTOfree(&hs) ;
  printf("difference image distributed as %2.2f +- %2.2f\n", mean, sigma) ;

  //  sigma = 2.0 ;
  if (pstd)
  {
    *pstd = sigma ;
    printf("setting std = %2.2f\n", sigma) ;
  }

  hg = HISTOgaussianCDF(NULL, mean,sigma, 10000);
  for (x1 = 0 ; x1 < mri1->width ; x1++)
    for (y1 = 0 ; y1 < mri1->height ; y1++)
      for (z1 = 0 ; z1 < mri1->depth ; z1++)
      {
        if (x1 == Gx && y1 == Gy && z1 == Gz)
          DiagBreak() ;
        TransformSampleInverse(transform, (float)x1, (float)y1, (float)z1, &x2, &y2, &z2) ;
        val1 = MRIgetVoxVal(mri1, x1, y1, z1, 0) ;
        MRIsampleVolume(mri2, x2, y2, z2, &val2) ;
        mask1 = MRIgetVoxVal(mri_mask1, x1, y1, z1, 0) ;
        MRIsampleVolume(mri_mask2, x2, y2, z2, &mask2) ;
        if (!FZERO(val1) && !FZERO(val2))
        {
          dif = ((val2-mean)-val1)  ;
          p = (1.0 / (sqrt(2*M_PI)*sigma)) * exp(-0.5 * (dif*dif) / (2*sigma*sigma)) ;
          p = HISTOgetCount(hg, fabs(dif));
          p = 1.0 - 1.0*p ;
          logp = -log10(p) ;
          if (std::isfinite(logp) == 0 || (DZERO(logp) && p < .1))
            MRIsetVoxVal(mri_big, x1, y1, z1, 0, 1) ;
          else
            MRIsetVoxVal(mri_change, x1, y1, z1, 0, logp) ;
        }
      }
  MRIvalRange(mri_dif, &min_change, &max_change) ;
  
  if (pstd)
    printf("setting underflow voxels to %2.1f\n", 2*max_change) ;
  for (x1 = 0 ; x1 < mri1->width ; x1++)
    for (y1 = 0 ; y1 < mri1->height ; y1++)
      for (z1 = 0 ; z1 < mri1->depth ; z1++)
      {
        if (MRIgetVoxVal(mri_big, x1, y1, z1, 0) == 0)
          continue ;
        MRIsetVoxVal(mri_change, x1, y1, z1, 0, 2*max_change) ;
      }

  HISTOfree(&hs) ; HISTOfree(&h) ; MRIfree(&mri_dif) ; MRIfree(&mri_mean) ; 
  MRIfree(&mri_mask1) ; MRIfree(&mri_mask2) ;
  return(mri_change) ;
}

MRI *
MRIcomputeNbhdPvalues(MRI *mri_src, int nbhd_size, MRI *mri_dst)
{
  int    x, y, z, xk, yk, zk, xi, yi, zi ;
  double pval, log_pnbhd, logp ;

  mri_dst = MRIclone(mri_src, mri_dst) ;

  for (x = 0 ; x < mri_src->width ; x++)
    for (y = 0 ; y < mri_src->height ; y++)
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        log_pnbhd = 0 ;
        for (xk = -nbhd_size ; xk <= nbhd_size ; xk++)
        {
          xi = mri_src->xi[x+xk] ;
          for (yk = -nbhd_size ; yk <= nbhd_size ; yk++)
          {
            yi = mri_src->yi[y+yk] ;
            for (zk = -nbhd_size ; zk <= nbhd_size ; zk++)
            {
              zi = mri_src->zi[z+zk] ;
              logp = MRIgetVoxVal(mri_src, xi, yi, zi, 0) ;
              pval = pow(10.0, -logp) ;
              logp = -log10(1.0-pval) ;
              
              log_pnbhd += logp ;
            }
          }
        }
        pval = pow(10.0, -log_pnbhd) ;
        logp = -log10(1.0-pval) ;
        MRIsetVoxVal(mri_dst, x, y, z, 0, logp) ;
      }

  return(mri_dst) ;
}
MRI *
MRIbonferroniCorrect(MRI *mri_src, MRI *mri_dst, int is_log) 
{
  int    x, y, z, nvox ;
  double logp, pval ;

  if (mri_dst == NULL)
    mri_dst = MRIclone(mri_src, NULL) ;


  for (nvox = x = 0 ; x < mri_src->width ; x++)
    for (y = 0 ; y < mri_src->height ; y++)
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        logp = MRIgetVoxVal(mri_src, x, y, z, 0) ;
        if (FZERO(logp) != 0)
          nvox++ ;
      }

  printf("using %d voxels for Bonferroni correction\n", nvox) ;
  for (x = 0 ; x < mri_src->width ; x++)
    for (y = 0 ; y < mri_src->height ; y++)
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        if (is_log)
        {
          logp = MRIgetVoxVal(mri_src, x, y, z, 0) ;
          pval = 1-pow(10.0, -logp) ;
          logp = -log10(pval) ;
          logp *= (float)nvox ;
          pval = 1-pow(10.0, -logp) ;
          logp = -log10(pval) ;
          MRIsetVoxVal(mri_src, x, y, z, 0, logp) ;
        }
      }

  return(mri_dst) ;
}

