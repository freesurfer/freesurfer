/**
 * @file  mri_apply_autoencoder.c
 * @brief main program for applying a stacked autoencoder for feature extraction.
 *
H.-C. Shin, M. R. Orton, D. J. Collins, S. J. Doran, and M. O. Leach,
"Stacked Autoencoders for
Unsupervised Feature Learning and Multiple Organ Detectionin a Pilot Study
Using 4D Patient Data,"
IEEE Transaction on Pattern Analysis and Machine Intelligence, 2012.

 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2013/11/22 19:41:39 $
 *    $Revision: 1.2 $
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



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "error.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "autoencoder.h"

//static VECTOR *extract_neighborhood(MRI *mri, int  whalf, int x0, int y0, int z0, VECTOR *v)  ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

static int synthesize = 0 ;

static int x0 = -1 ;
static int x1 ;
// for some reason y0 and y1 are defined in /usr/include/bits/mathcalls.h so have to include annoying _
static int y0_ ;
static int y1_ ;

static int z0 ;
static int z1 ;

#define MAX_PYR_LEVELS 100 
static int nlevels = 4 ;

static MRI *mri_train ;

static int ras_point_set = 0 ;
static double Gr = -1 ;
static double Ga = -1 ;
static double Gs = -1 ;
static int read_flag = 0 ;

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs ;
  char   *in_fname, *out_fname, *sae_fname ;
  int          msec, minutes, seconds, n ;
  Timer start ;
  MRI          *mri, *mri_orig, *mri_scaled, *mri_ae_out, *mri_dot_out, *mri_orig_cropped, *mri_pyramid[MAX_PYR_LEVELS],
               *mri_train_pyramid[MAX_PYR_LEVELS]   ;
  SAE          *sae ;
  double       mean, train_mean ;
  AE           *ae_last ; // deepest layer

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_apply_autoencoder.c,v 1.2 2013/11/22 19:41:39 fischl Exp $", "$Name:  $");
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

  in_fname = argv[1] ;
  sae_fname = argv[2] ;
  out_fname = argv[3] ;
  if (argc < 4)
    usage_exit(1) ;
//  setRandomSeed(-1L) ;

  mri = MRIread(in_fname) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume from %s", Progname, in_fname) ;

  if (ras_point_set >= 0)
  {
    double x, y, z ;
    MRIsurfaceRASToVoxel(mri_train, Gr, Ga, Gs,&x, &y, &z) ;
    Gx = nint(x) ; Gy = nint(y) ; Gz = nint(z) ;
    printf("RAS point (%2.0f, %2.0f, %2.0f) maps to voxel (%d, %d, %d)\n", Gr, Ga, Gs, Gx, Gy, Gz) ;
  }

  if (read_flag)
  {
    LABEL *area ;
    char   fname[STRLEN] ;
//    HISTOGRAM *h, *hcdf ;
    int        bin ;
    double     thresh, xv, yv, zv ;
    float      fmin, fmax ;
    MRI        *mri_ae_p ;

    sprintf(fname, "%s.AE.p.mgz", out_fname) ;
    printf("reading AE p-val product from  %s\n", fname) ;
    mri_ae_p = MRIread(fname) ;
    if (mri_ae_p == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read precomputed AE p-value map from %s\n", Progname, fname) ;

    in_fname = "cube.inputs.label" ;
    out_fname = "cube.outputs.label" ;
    if (FileExists(in_fname) == 0)
    {
      area = LabelAlloc(1, NULL, in_fname) ;
      if (area == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not create label file %s\n", Progname,in_fname) ;
      LabelToVoxel(area, mri_train, area) ;
    }
    else
    {
      area = LabelRead(NULL, in_fname) ;
      if (area == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read label file %s\n", Progname,in_fname) ;
      area = LabelRealloc(area, area->n_points+1) ;
    }

    bin = area->n_points++ ;
    area->lv[bin].vno = -1 ;
    area->lv[bin].deleted = 0 ;
    area->lv[bin].stat = 0 ;
    area->lv[bin].x = Gx ;
    area->lv[bin].y = Gy ;
    area->lv[bin].z = Gz ;
    printf("writing label file with %d points\n", area->n_points) ;
    LabelWrite(area, in_fname) ;

    if (FileExists(out_fname) == 0)
    {
      area = LabelAlloc(1, NULL, out_fname) ;
      if (area == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not create label file %s\n", Progname, out_fname) ;
      LabelToVoxel(area, mri_train, area) ;
    }
    else
    {
      area = LabelRead(NULL, out_fname) ;
      if (area == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read label file %s\n", Progname, out_fname) ;
      area = LabelRealloc(area, area->n_points+1) ;
    }

#if 0
    h = MRIhistogram(mri_ae_p, 100000) ;
    hcdf = HISTOmakeCDF(h, NULL) ;
    bin = HISTOfindBinWithCount(hcdf, .99999) ;
    thresh = h->bins[bin] ;
#endif
    MRIvalRange(mri_ae_p, &fmin, &fmax) ;
    thresh = fmax ;
    MRIthreshold(mri_ae_p, mri_ae_p, thresh) ;
    MRIcomputeCentroid(mri_ae_p, &xv, &yv, &zv) ;
    printf("thresholding SAE output at %f, centroid at (%2.1f, %2.1f, %2.1f)\n", thresh, xv, yv, zv) ;
    bin = area->n_points++ ;
    area->lv[bin].deleted = 0 ;
    area->lv[bin].stat = thresh ;
    area->lv[bin].vno = -1 ;
    area->lv[bin].x = xv ;
    area->lv[bin].y = yv ;
    area->lv[bin].z = zv ;
    LabelWrite(area, out_fname) ;
    exit(0) ;
  }

  mri_orig = MRIcopy(mri, NULL) ;
  if (mri->type == MRI_UCHAR)
  {
    MRI *mri_tmp = MRIcloneDifferentType(mri, MRI_FLOAT) ;
    MRIscalarMul(mri, mri_tmp, 1.0/255.0) ;
    MRIfree(&mri) ; mri = mri_tmp ;
  }
  mri_scaled = MRIcopy(mri, NULL) ;
  mean = MRImeanFrame(mri, 0) ;
  MRIaddScalar(mri, mri, -mean) ;

  if (x0 >= 0)
  {
    MRI *mri_tmp ;
    mri_tmp = MRIextract(mri, NULL, x0, y0_, z0, x1-x0+1, y1_-y0_+1, z1-z0+1) ;
    MRIfree(&mri) ;
    mri = mri_tmp ;
    mri_orig_cropped = MRIextract(mri, NULL, x0, y0_, z0, x1-x0+1, y1_-y0_+1, z1-z0+1) ;
  }
  else
    mri_orig_cropped = mri_orig ;

  mri_pyramid[0] = mri ;
  for (n = 1 ; n < nlevels ; n++)
    mri_pyramid[n] = MRIreduce(mri_pyramid[n-1], NULL) ;

  if (mri_train->type == MRI_UCHAR)
  {
    MRI *mri_tmp = MRIcloneDifferentType(mri_train, MRI_FLOAT) ;
    MRIscalarMul(mri_train, mri_tmp, 1.0/255.0) ;
    MRIfree(&mri_train) ; mri_train = mri_tmp ;
  }
  train_mean = MRImeanFrame(mri_train, 0) ;
  MRIaddScalar(mri_train, mri_train, -train_mean) ;
  mri_train_pyramid[0] = mri_train ;
  for (n = 1 ; n < nlevels ; n++)
    mri_train_pyramid[n] = MRIreduce(mri_train_pyramid[n-1], NULL) ;
  sae = SAEread(sae_fname) ;
  if (sae == NULL)
    ErrorExit(Gerror, "") ;
  ae_last = SAEfindLastLayer(sae, NULL) ;
  {
    int     x, y, z ;
    VECTOR  *v_src, *v_dst, *v_ae_src, *v_ae_dst, *v_src_orig, *v_ae_src_orig ;
    char    fname[STRLEN] ;
    double  rms ;
    MRI     *mri_ae_p, *mri_p ;

    if (x0 >= 0)
    {
      Gx -= x0 ; Gy -= y0_ ; Gz -= z0 ;
    }

    mri_ae_out = MRIclone(mri, NULL) ;
    mri_dot_out = MRIclone(mri, NULL) ;
    mri_ae_p = MRIclone(mri, NULL) ;
    mri_p = MRIclone(mri, NULL) ;
    printf("computing feature similarity using auto-encoder...\n") ;

    SAEfillInputVector(mri_train_pyramid, sae->nlevels, Gx, Gy, Gz, sae->whalf, sae->first->v_input) ;
    SAEactivateNetwork(sae) ;
    v_ae_src_orig = VectorCopy(ae_last->v_hidden, NULL) ; 
    v_ae_src = VectorZeroMean(v_ae_src_orig, NULL) ; 

    v_src_orig = SAEfillInputVector(mri_train_pyramid, sae->nlevels, Gx, Gy, Gz, sae->whalf, NULL) ;
    v_src = VectorZeroMean(v_src_orig, NULL) ; 
    v_dst = VectorClone(v_src) ; 
    v_ae_dst = VectorClone(v_ae_src) ; 

    for (x = sae->whalf ; x < mri->width-sae->whalf; x++)
      for (y = sae->whalf ; y < mri->height-sae->whalf ; y++)
	for (z = sae->whalf ; z < mri->depth-sae->whalf ; z++)
	{
	  if (x == Gx && y == Gy && z == Gz)
	    DiagBreak() ;
	  if (MRIgetVoxVal(mri_orig_cropped, x, y, z, 0) == 0)
	    continue ;
	  SAEfillInputVector(mri_pyramid, sae->nlevels, x, y, z, sae->whalf, sae->first->v_input) ;
	  SAEactivateNetwork(sae) ;
	  rms = VectorRMS(ae_last->v_hidden, v_ae_src_orig) ;
#define SIGMA_AE .1
	  MRIsetVoxVal(mri_ae_p, x, y, z, 0, exp(-rms/SIGMA_AE)) ;
	  MRIsetVoxVal(mri_ae_p, x, y, z, 0, (1-rms)) ;

	  VectorZeroMean(ae_last->v_hidden, v_ae_dst) ;
	  if (VectorNormalizedDot(v_ae_src, v_ae_dst) < 0)
	    DiagBreak() ;
	  MRIsetVoxVal(mri_ae_out, x, y, z, 0, VectorNormalizedDot(v_ae_src, v_ae_dst)) ;

	  SAEfillInputVector(mri_pyramid, sae->nlevels, x, y, z, sae->whalf, v_dst) ;
	  rms = 1+VectorRMS(v_src_orig, v_dst) ;
	  MRIsetVoxVal(mri_p, x, y, z, 0, exp(-rms/SIGMA_AE)) ;
	  MRIsetVoxVal(mri_p, x, y, z, 0, (1-rms)) ;
	  VectorZeroMean(v_dst, v_dst) ;
	  MRIsetVoxVal(mri_dot_out, x, y, z, 0, VectorNormalizedDot(v_src, v_dst)) ;
	} 
    VectorFree(&v_src) ;
    VectorFree(&v_ae_src) ;
    VectorFree(&v_dst) ;
    VectorFree(&v_ae_dst) ;
    if (x0 < 0)
    {
      x0 = y0_ = z0 = 0 ;
      x1 = mri->width-1 ;
      y1_ = mri->height-1 ;
      z1 = mri->depth-1 ;
    }
    if (x0 >= 0)
    {
      MRI *mri_tmp = MRIclone(mri_scaled, NULL) ;
      MRIextractInto(mri_ae_out, mri_tmp, 0, 0, 0, mri_ae_out->width, mri_ae_out->height, mri_ae_out->depth, x0, y0_, z0) ;
      MRIfree(&mri_ae_out) ; mri_ae_out = mri_tmp ;

      mri_tmp = MRIclone(mri_scaled, NULL) ;
      MRIextractInto(mri_dot_out, mri_tmp, 0, 0, 0, mri_dot_out->width, mri_dot_out->height, mri_dot_out->depth, x0, y0_, z0) ;
      MRIfree(&mri_dot_out) ; mri_dot_out = mri_tmp ;

    }
    sprintf(fname, "%s.dot.mgz", out_fname) ;
    printf("writing dot product to %s\n", fname) ;
//    MRIwrite(mri_dot_out, fname) ;
    sprintf(fname, "%s.p.mgz", out_fname) ;
    printf("writing p-val to %s\n", fname) ;
//    MRIwrite(mri_p, fname) ;

    sprintf(fname, "%s.AE.mgz", out_fname) ;
    printf("writing AE dot product to %s\n", fname) ;
//    MRIwrite(mri_ae_out, fname) ;
    sprintf(fname, "%s.AE.p.mgz", out_fname) ;
    printf("writing AE p-val product to %s\n", fname) ;
    MRIwrite(mri_ae_p, fname) ;

  }
  if (synthesize)
  {
    int x, y, z, wsize = 2*sae->whalf+1, ind  ;
    MRI *mri_tmp, *mri_total, *mri_tmp2 ;
    VECTOR *v ;
    char   path[STRLEN], fname[STRLEN] ;

    printf("synthesizing output volume using auto-encoder...\n") ;

    ind = wsize*wsize*wsize /2 + 1 ;
    FileNameRemoveExtension(out_fname, path) ;
    mri_tmp = MRIalloc(wsize, wsize, wsize, mri->type) ;
    mri_tmp2 = MRIclone(mri, NULL) ;
    mri_total = MRIcopy(mri, NULL) ;
//    MRIscalarMul(mri_total, mri_total, (float)(wsize*wsize*wsize)) ;
    for (x = sae->whalf ; x < mri->width-sae->whalf; x++)
      for (y = sae->whalf ; y < mri->height-sae->whalf ; y++)
	for (z = sae->whalf ; z < mri->depth-sae->whalf ; z++)
	{
	  if (x == Gx && y == Gy && z == Gz)
	    DiagBreak() ;
	  if (FZERO(MRIgetVoxVal(mri, x, y, z, 0)))
	    continue ;
	  SAEfillInputVector(mri_pyramid, sae->nlevels, x, y, z, sae->whalf, ae_last->v_input) ;
	  v = SAEactivateNetwork(sae) ;
	  SAEvectorToMRI(v, sae->nlevels, sae->whalf, mri_tmp) ;
	  MRIclear(mri_tmp2) ;
	  MRIextractInto(mri_tmp, mri_tmp2,0,0,0, mri_tmp2->width, mri_tmp2->height, mri_tmp2->depth, x, y, z) ;
	  MRIcopyHeader(mri_total, mri_tmp2) ;
//	  MRIadd(mri_tmp2, mri_total, mri_total) ;
	  MRIsetVoxVal(mri_total, x, y, z, 0, ae_last->v_output->rptr[ind][1]) ;
	} 
//   MRIscalarMul(mri_total, mri_total, 1.0/(float)(wsize*wsize*wsize)) ;
    if (x0 < 0)
    {
      x0 = y0_ = z0 = 0 ;
      x1 = mri->width-1 ;
      y1_ = mri->height-1 ;
      z1 = mri->depth-1 ;
    }
    if (mri_orig->type == MRI_UCHAR)
    {
      MRIwrite(mri_total, "t1.mgz") ;
      MRIextractInto(mri_total, mri_scaled, 0, 0, 0, mri_total->width, mri_total->height, mri_total->depth, x0, y0_, z0) ;
      MRIwrite(mri_scaled, "synth.mgz") ;
      MRIscalarMul(mri_total, mri_total, 255.0f) ;
      mri_total = MRIchangeType(mri_total, MRI_UCHAR, 0, 255, 1) ;
    }
    MRIextractInto(mri_total, mri_orig, 0, 0, 0, mri_total->width, mri_total->height, mri_total->depth, x0, y0_, z0) ;
    sprintf(fname, "%s.out.mgz", path) ;
    printf("writing synthesized output to %s\n", fname) ;
    MRIwrite(mri_orig, fname) ;
    MRIfree(&mri_total) ; MRIfree(&mri_tmp) ;
  }

  SAEfree(&sae) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("autoeencoder application took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "RAS"))
  {
    ras_point_set = 1 ;
    Gr = atoi(argv[2]) ;
    Ga = atoi(argv[3]) ;
    Gs = atoi(argv[4]) ;
    nargs = 3 ;
    printf("applying SAE at TK RAS point (%2.1f, %2.1f, %2.1f)\n", Gr, Ga, Gs) ;
  }
  else switch (toupper(*option)) {
    case 'T':
      mri_train = MRIread(argv[2]) ;
      if (mri_train == NULL)
	ErrorExit(ERROR_NOFILE, "") ;
      printf("using training volume %s to compute similarity measures\n", argv[2]) ;
      nargs = 1 ;
      break ;
    case 'R':
      read_flag = 1 ;
      printf("reading in pre-computed results and creating labels\n") ;
      break ;
    case 'S':
      synthesize = 1 ;
      printf("synthesizing volume using autoencoder\n") ;
      break ;
    case 'X':
      x0 = atoi(argv[2]) ;
      x1 = atoi(argv[3]) ;
      y0_ = atoi(argv[4]) ;
      y1_ = atoi(argv[5]) ;
      z0 = atoi(argv[6]) ;
      z1 = atoi(argv[7]) ;
      nargs = 6 ;
      printf("extracting subregion X %d -->%d, y %d --> %d, z %d --> %d\n", x0, x1, y0_, y1_, z0, z1) ;
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
  printf("usage: %s [options] <input image> <autoencoder> <output image>\n", Progname) ;
  exit(code) ;
}



#if 0
static VECTOR *
extract_neighborhood(MRI **mri, int nlevels, int  whalf, int x0, int y0, int z0, VECTOR *v) 
{
  int wsize = 2*whalf+1, i, x, y, z, xk, yk, zk ;

  if (v == NULL)
    v = VectorAlloc(wsize*wsize*wsize, MATRIX_REAL) ;


  for (i = 1, xk = -whalf ; xk <= whalf ; xk++)
  {
    x = mri->xi[x0+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      y = mri->yi[y0+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++, i++)
      { 
	z = mri->zi[z0+zk] ;
	VECTOR_ELT(v, i) = MRIgetVoxVal(mri, x, y, z, 0) ;
      }
    }
  }

  return(v) ;
}

#endif
