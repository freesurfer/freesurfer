/**
 * @brief main program for registring a set of OCT tiles into a mosaic
 *
 * program for registering a set of OCT subimages (tiles) into a single mosaic.
 * Will peform registration and potentially optical distortion correction
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

#include "romp_support.h"


#include "mri.h"
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
#include "numerics.h"
#include "histo.h"
#include "fsinit.h"

#define DX_I  1
#define DY_I  2
#define D_I   3
#define B_I   4
#define C_I   5
#define A_I   6
#define AX0_I 7
#define AY0_I 8
#define AX1_I 9
#define AY1_I 10
//#define NPARMS 10
#define NPARMS 2


static int downsample = 0 ;
static int undistort = 0 ;
static char Gout_name[STRLEN] ;
static double blur_sigma = 0.0 ;

static int pad = 0 ;

static float intensity_tolerance = 0.1 ;
static float translation_tolerance = 0.025 ;

static MRI *mosaic_images_with_xforms(MRI **mri, int nimages, MRI *mri_weights, int pad)  ;
static int insert_image_into_mosaic(MRI *mri_mosaic, MRI *mri, double x0, double y0) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static MRI *mosaic_images(MRI **mri, int *x0, int *y0, int nimages, int pad) ;
//static MRI *add_image_to_mosaic(MRI *mri_orig_mosaic, MRI *mri, double x0, double y0, double ax, double ay, double a, double b, double c, double d) ;
static MRI *undistort_and_mosaic_images(MRI **mri, double *x0, double *y0, int nimages, double *ax, double *ay, double a, double b, double c, double d, double *penergy) ;
static double compute_mosaic_energy(MRI *mri_mosaic, int nimages)  ;

#if 0
static int    jacobian_correct(MRI *mri_src, MRI *mri_dst, double ax, double ay, double a, double b, double c, double d) ;
#endif

const char *Progname ;
static void usage_exit(int code) ;

static double compute_pairwise_deformation_energy(MRI *mri1, MRI *mri2, double dx, double dy, double *ax, double *ay, double a, double b, double c, double d) ;
static int undistorted_coords(MRI *mri, double x, double ax, double ay, double y, double a, double b, double c, double d, double *px, double *py) ;
static float compute_powell_sse(float *p) ;
static float compute_powell_global_sse(float *p) ;
static MRI *build_best_mosaic(MRI **mri, MRI **mri_smooth, double *x0d, double *y0d, int nimages) ;
static int powell_minimize(MRI *mri1, MRI *mri2, double dx, double dy, double a, double b, double c, double d, double *ax, double *ay, 
			   double *pa, double *pb, double *pc, double *pd, double *pdx, double *pdy) ;

static double a = 0.0, b = 0.0, c = 0.0, d = 1.0 ;

static double xsize = 0, ysize = 0, zsize = 0 ;

static int niters = 100 ;

#define MAX_IMAGES 10000
#define MAX_TRANS 3

static char *weight_fname = NULL ;

int
main(int argc, char *argv[]) {
  char       **av, *mosaic_input_fname, *out_fname ;
  int         ac, nargs, msec, minutes, seconds, nimages, i ;
  int          dx, dy, n, *image_numbers, frame, dx_best, dy_best, skip ;
  Timer start ;
  MRI          *mri[MAX_IMAGES], *mri_orig[MAX_IMAGES], *mri_mosaic, *mri_smooth[MAX_IMAGES],
               *mri_weights  ;
  FILE         *fp = NULL ;
  double       energy, best_energy, prev_best_energy, dxd, dyd ;
  int          x0[MAX_IMAGES], y0[MAX_IMAGES], xbest[MAX_IMAGES], ybest[MAX_IMAGES] ; 
  double       x0d[MAX_IMAGES], y0d[MAX_IMAGES], ax[MAX_IMAGES], ay[MAX_IMAGES] ;
//  HISTOGRAM    *h0 = NULL, *h1 ;

  nargs = handleVersionOption(argc, argv, "oct_register_mosaic");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

#ifdef HAVE_OPENMP
  if (getenv("OMP_NUM_THREADS") == NULL)
    omp_set_num_threads(1);
#endif
//  FSinit() ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  setRandomSeed(-1L) ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  mosaic_input_fname = argv[1] ;
  out_fname = argv[argc-1] ;
  if (argc < 2)
    usage_exit(1) ;

  if (argc == 3 && 0)  // disable for now
  {
    FileNameRemoveExtension(out_fname, Gout_name) ;
    printf("reading mosaic file names and initial positions from %s\n", mosaic_input_fname) ;
    
    nimages = FileNumberOfEntries(mosaic_input_fname) ;

    fp = fopen(mosaic_input_fname, "r") ;
    if (fp == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not open input list file %s\n", Progname, mosaic_input_fname) ;
  }
  else
  {
    nimages = argc-2 ;
  }
  printf("mosaicing %d images...\n", nimages) ;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) firstprivate(downsample, blur_sigma)  shared(x0, y0, xbest, ybest, mri, mri_orig, ax, ay, mri_smooth)
#endif
  for (i = 0 ; i < nimages ; i++)
  {
    ROMP_PFLB_begin
    
    VOL_GEOM   vg;
    char      *cp, fname[STRLEN], line[STRLEN] ;
    if (fp)
    {
      cp = fgetl(line, 199, fp) ;
      sscanf(cp, "%d %d %s", &x0[i], &y0[i], fname) ;
      printf("reading input volume %d of %d: %s at (%d, %d)\n", i+1, nimages, fname, x0[i], y0[i]) ;
    }
    else
      strcpy(fname, argv[i+1]) ;
  

    xbest[i] = x0[i] ; ybest[i] = y0[i] ;
    mri[i] = MRIread(fname) ;
    if (mri[i] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s\n", Progname, fname) ;
    getVolGeom(mri[i], &vg);
    if (fp == NULL)
      printf("reading input volume %d of %d %s at (%2.1f, %2.1f, %2.1f)\n", i+1, nimages, fname, vg.c_r, vg.c_a, vg.c_s) ;
    if (xsize > 0)
      vg.xsize = xsize ;
    if (ysize > 0)
      vg.ysize = ysize ;
    if (zsize > 0)
      vg.zsize = zsize ;
    useVolGeomToMRI(&vg,mri[i]);

    if (downsample)
    {
      MRI *mri_tmp ;
      int d ;

      for (d = 0 ; d < downsample ; d++)
      {
#if 0
	mri_tmp = MRIdownsample2(mri[i], NULL) ;
#else
	mri_tmp = MRIreduce(mri[i], NULL) ;
#endif
	MRIfree(&mri[i]) ;
	mri[i] = mri_tmp ;
      }
    }
#if 0
    vg.c_r += (x0[i]*vg.xsize) ;
    vg.c_s -= (y0[i]*vg.ysize) ;
#endif

    mri_orig[i] = MRIcopy(mri[i],NULL) ;
    sprintf(fname, "v%d.mgz", i);
#if 0
    if (i == 0)
      h0 = MRIhistogram(mri[i], 256) ;
    else
    {
      float      a, b ;
      if (h0 == NULL)
	h0 = MRIhistogram(mri[i], 256) ;
      h1 = MRIhistogram(mri[i], 256) ;
      HISTOfindLinearFit(h1, h0, .8, 1.2, -10, 10, &a, &b) ;
      printf("scaling image %d intensities by %2.3f x + %2.0f to match baseline image\n", i, a, b) ;
      MRIscaleIntensities(mri[i], mri[i], a, b) ;
    }
#endif
    ax[i] = (mri[i]->width-1.0)/2.0 ;
    ay[i] = (mri[i]->height-1.0)/2.0 ;
    if (blur_sigma > 0)
    {
      MRI *mri_blur = MRIgaussian1d(blur_sigma, 0) ;
      mri_smooth[i] = MRIconvolveGaussian(mri[i], NULL, mri_blur) ;
      MRIfree(&mri_blur) ;
    }
    else
      mri_smooth[i] = MRIcopy(mri[i], NULL) ;

//    MRIwrite(mri[i], fname) ;
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  if (fp)
    fclose(fp) ;

  if (downsample > 0)
      printf("final downsampling image size = %d x %d x %d, size %2.3f x %2.3f x %2.3f um^3\n", 
	     mri[0]->width, mri[0]->height, mri[0]->depth,
	     1000*mri[0]->xsize, 1000*mri[0]->ysize, 1000*mri[0]->zsize) ;

  {
    char fname[STRLEN], ext[STRLEN], orig_fname[STRLEN] ;
    
    if (fp)
    {
      FileNameExtension(out_fname, ext) ;
      FileNameRemoveExtension(out_fname, fname) ;
    } else {
      sprintf(fname, "mosaic") ;
      int req = snprintf(orig_fname, STRLEN, "%s.orig.%s", fname, ext) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    for (i = 0 ; i < nimages ; i++)
    {
      x0d[i] = x0[i] ;
      y0d[i] = y0[i] ;
    }
    if (weight_fname)
    {
      mri_weights = MRIread(weight_fname) ;
      if (mri_weights == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read weights from %s", weight_fname) ;
    }
    else
      mri_weights = NULL ;
    mri_mosaic = mosaic_images_with_xforms(mri, nimages, mri_weights, 0) ;
    if (undistort)  // since this isn't the final image, save a copy
    {
      printf("writing mosaic with original translations to %s\n", orig_fname) ;
      MRIwrite(mri_mosaic, orig_fname) ;
    }
  }

  if (undistort)
  {
    if (mri_mosaic)
      MRIfree(&mri_mosaic) ;
    mri_mosaic = mosaic_images(mri, x0, y0, nimages, pad) ;
    MRIcopyHeader(mri[0], mri_mosaic) ;
    MRIclear(mri_mosaic) ;
    for (i = 0 ; i < nimages ; i++)
    {
      x0[i] += pad ;
      y0[i] += pad ;
      x0d[i] += pad ;
      y0d[i] += pad ;
    }
    
    mri_mosaic = build_best_mosaic(mri, mri_smooth, x0d,y0d, nimages) ;
  }
  printf("writing final mosaic to %s\n", out_fname) ;
  MRIwrite(mri_mosaic, out_fname) ;
  exit(0) ;
  
  insert_image_into_mosaic(mri_mosaic, mri[0], x0d[0], y0d[0]) ;
  for (i = 1 ; i < nimages ; i++)
  {
    powell_minimize(mri_mosaic, mri[i], x0d[i], y0d[i],  a,  b,  c,  d, ax, ay, &a, &b, &c, &d, &dxd, &dyd) ;

    if (NPARMS == 2)
      printf("%d: best parameters found at dr = (%2.1f, %2.1f)\n", i, dxd, dyd) ;
    else
    {
      printf("%d: best parameters found at r = %2.4f r^4 + %2.4f r^3 + %2.4f r^2 + %2.4f r, dr = (%2.1f, %2.1f)\n", i, a, b, c, d, dxd, dyd) ;
      printf("optic axes:\n") ;
      for (i = 0 ; i < nimages ; i++)
	printf("(%2.1f, %2.1f)\n", ax[i], ay[i]) ;
    }
    x0d[i] = dxd ;
    y0d[i] = dyd;
    insert_image_into_mosaic(mri_mosaic, mri[i], x0d[i], y0d[i]) ;

    {
      char fname[STRLEN], ext[STRLEN], mosaic_fname[STRLEN] ;
      
      FileNameExtension(out_fname, ext) ;
      FileNameRemoveExtension(out_fname, fname) ;
      sprintf(mosaic_fname, "%s.m%d.%s", fname, i, ext) ;
      MRIwrite(mri_mosaic, mosaic_fname) ;
    }
  }
  printf("final positions\n") ;
  for (i = 0 ; i < nimages ; i++)
  {
    printf("(%2.0f, %2.0f)\n", x0d[i], y0d[i]) ;
  }
  printf("writing final mosaic to %s\n", out_fname) ;

  mri_mosaic = undistort_and_mosaic_images(mri, x0d, y0d, nimages, ax, ay, a, b, c, d, NULL) ;
  MRIwrite(mri_mosaic, out_fname) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ; minutes = seconds / 60 ; seconds = seconds % 60 ;
  fprintf(stderr, "OCT mosaicing  took %d minutes and %d seconds.\n",  minutes, seconds) ;
  exit(0) ;

  mri_mosaic = mosaic_images(mri,x0, y0, nimages, 0) ;
  best_energy = energy = compute_mosaic_energy(mri_mosaic, nimages) ;
  printf("original mosaic energy = %2.3f\n", energy) ;
  for (skip = 1 ; skip > 0 ; skip--)
  {
    int max_trans = skip*MAX_TRANS ;
    for (n = 0 ; n < 10 ; n++)
    {
      image_numbers = compute_permutation(nimages, NULL)  ;
      prev_best_energy = best_energy ;
      for (i = 0 ; i < nimages ; i++)
      {
	frame = image_numbers[i] ; frame = i ;
	dx_best = dy_best = 0 ;
	for (dx = -max_trans ; dx <= max_trans ; dx += skip)
	  for (dy = -max_trans; dy <= max_trans ; dy += skip)
	  {
	    x0[frame] = xbest[frame]+dx ; y0[frame] = ybest[frame]+dy ;
	    mri_mosaic = mosaic_images(mri,x0, y0, nimages, 0) ;
	    energy = compute_mosaic_energy(mri_mosaic, nimages) ;
	    if (energy < best_energy)
	    {
	      printf("new best energy %2.5f found for image %d at (%d, %d)\n", energy, frame, dx, dy) ;
	      best_energy = energy ;
	      dx_best = dx ;
	      dy_best = dy ;
	    }
	    MRIfree(&mri_mosaic) ;
	  }
	xbest[frame] += dx_best ;
	ybest[frame] += dy_best ;
      }
      printf("skip %d, completed iter %d, best energy = %2.3f\n", skip, n, best_energy) ;
      mri_mosaic = mosaic_images(mri,xbest, ybest, nimages, 0) ;
      MRIwrite(mri_mosaic, out_fname) ;
      free(image_numbers) ;
      if (FEQUAL(best_energy, prev_best_energy))
	break ;
    }
    printf("skip %d completed best energy = %2.3f\n", skip, best_energy) ;
  }
  printf("final positions\n") ;
  for (i = 0 ; i < nimages ; i++)
  {
    printf("(%d, %d)\n", xbest[i], ybest[i]) ;
  }

  mri_mosaic = mosaic_images(mri,xbest, ybest, nimages, 0) ;
  MRIwrite(mri_mosaic, out_fname) ;
  
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ; minutes = seconds / 60 ; seconds = seconds % 60 ;
  fprintf(stderr, "OCT mosaicing  took %d minutes and %d seconds.\n",  minutes, seconds) ;
  exit(0) ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#include "mrisegment.h"
#include "label.h"
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
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
  }
#if 0
  else if (!stricmp(option, "hough"))
  {
    MRI *mri_hough ; 
    MRI_SEGMENTATION *mseg ;

    mri_hough = MRIread(argv[2]) ;
    printf("reading in Hough transform of 25um grid from %s...\n", argv[2]) ;
    if (mri_hough == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read in Hough transformed grid from %s...\n", argv[2]) ;
    nargs = 1 ;
    mseg = MRIsegment(mri_hough, 700, 2000) ;
    MRIremoveSmallSegments(mseg, 20) ;
    printf("%d Hough segments found\n", mseg->nsegments) ;
    exit(1) ;
  }
  else if (!stricmp(option, "unwarp"))
  {
    LABEL *grid_points ;

    grid_points = MRIread(argv[2]) ;
    printf("reading grid points from %s...\n", argv[2]) ;
    if (mri_hough == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read grid from %s...\n", argv[2]) ;
    nargs = 1 ;
    compute_grid_deformation(grid_points) ;
    exit(0) ;
  }
#endif
  else if (!stricmp(option, "xsize"))
  {
    xsize = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting x voxel size to %2.6f\n", xsize) ;
  }
  else if (!stricmp(option, "ysize"))
  {
    ysize = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting y voxel size to %2.6f\n", ysize) ;
  }
  else if (!stricmp(option, "zsize"))
  {
    zsize = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting z voxel size to %2.6f\n", zsize) ;
  }
  else switch (toupper(*option)) {
    case 'D':
      downsample = atoi(argv[2]) ;
      nargs = 1 ;
      printf("downsampling input volumes %d times (a factor of %d) before mosaicing\n", downsample, (int)exp2(downsample)) ;
      break ;
    case 'S':
      blur_sigma = atof(argv[2]) ;
      printf("blurring input with sigma=%2.2f\n", blur_sigma) ;
      nargs = 1 ;
      break ;
    case 'U':
      undistort = 1 ;
      printf("undistorting images before mosaicing\n") ;
      break ;
    case 'W':
      weight_fname = argv[2] ;
      printf("using weighting file %s\n", weight_fname) ;
      nargs = 1 ;
      break ;
    case 'I':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      printf("debugging image # %d\n", Gdiag_no) ;
      break ;
    case 'P':
      pad = atoi(argv[2]) ;
      nargs = 1 ;
      printf("padding final mosaic with %d voxels\n", pad) ;
      break ;
    case 'N':
      niters = atoi(argv[2]) ;
      printf("setting max iters = %d\n", niters) ;
      nargs = 1 ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      usage_exit(1) ;
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
  printf("usage 1: %s [options] <tile 1> <tile 2> ...  <output volume>\n", Progname) ;
  printf("usage 2: %s [options] <mosaic list file>  <output volume>\n", Progname) ;
  printf("\toptions:\n") ;
  printf("\t-D <downsample>\t-\tuse Gaussian downsampling <downsample> times\n") ;
  printf("\t-W <weight fname>\t-\tread tile of weights from <weight fname> and use it in tile averaging\n") ;
  exit(code) ;
}
static MRI *
mosaic_images(MRI **mri, int *x0, int *y0, int nimages, int pad) 
{
  int      xmin, ymin, xmax, ymax, x1, y1, x2, y2 ;
  MRI      *mri_mosaic, *mri_counts ;
  int      i, width, height, x ;

  ymin = xmin = 10000000 ;
  ymax = xmax = -ymin ;

  for (i =  0 ; i < nimages ; i++)
  {
    x1 = x0[i]-pad ; y1 = y0[i]-pad ; x2 = 2*pad+x1 + mri[i]->width-1 ; y2 = 2*pad+y1 + mri[i]->height-1 ;
    if (x1 < xmin)
      xmin = x1 ;
    if (y1 < ymin)
      ymin = y1 ;
    if (x2 > xmax)
      xmax = x2 ;
    if (y2 > ymax)
      ymax = y2 ;
  }
  
  width = xmax - xmin + 1 ;
  height = ymax - ymin + 1 ;
  mri_mosaic = MRIallocSequence(width, height, 1, MRI_FLOAT, nimages+1) ;
  MRIcopyHeader(mri[0], mri_mosaic) ;

  mri_counts = MRIclone(mri_mosaic, NULL) ;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) 
#endif
  for (i =  0 ; i < nimages ; i++)
  {
    ROMP_PFLB_begin
    int x, y, count, xm, ym ;
    float  val, sum ;
    for (x = 0 ; x < mri[i]->width ; x++)
    {
      xm = x + x0[i]-xmin ;
      for (y = 0 ; y < mri[i]->height ; y++)
      {
	if ( x== Gx && y == Gy)
	  DiagBreak() ;
	ym = y + y0[i]-ymin ;
	count = (int)MRIgetVoxVal(mri_counts, xm, ym, 0, 0) ;
	count++ ;
	MRIsetVoxVal(mri_counts, xm, ym, 0, 0, count) ;
	sum = MRIgetVoxVal(mri_mosaic, xm, ym, 0, 0) ;
	val = MRIgetVoxVal(mri[i], x, y, 0, 0) ;
	MRIsetVoxVal(mri_mosaic, xm, ym, 0, 0, val+sum) ;
	MRIsetVoxVal(mri_mosaic, xm, ym, 0, i+1, val) ;
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) 
#endif
  for (x = 0 ; x < mri_mosaic->width ; x++)
  {
    ROMP_PFLB_begin
    
    float sum, count ;
    int   y ;
    for (y = 0 ; y < mri_mosaic->height ; y++)
    {
      count = MRIgetVoxVal(mri_counts, x, y, 0, 0) ;
      if (count <= 1)
	continue ;
      sum = MRIgetVoxVal(mri_mosaic, x, y, 0, 0) ;
      MRIsetVoxVal(mri_mosaic, x, y, 0, 0, sum/count) ;
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

  MRIfree(&mri_counts) ;
  return(mri_mosaic) ;
}


static double
compute_mosaic_energy(MRI *mri_mosaic, int nimages) 
{
  int    x, overlap_voxels, min_overlap_voxels, nonzero ;
  double energy ;


  energy =  0.0 ;
  nonzero = overlap_voxels = 0 ;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) reduction (+:overlap_voxels,nonzero,energy)  
#endif
  for (x = 0 ; x < mri_mosaic->width ; x++)
  {
    ROMP_PFLB_begin
    int    y, count, frame ;
    double var, dif, val, mean, vals[MAX_IMAGES] ; ;
    for (y = 0 ; y < mri_mosaic->height ; y++)
    {
      if (x == Gx && y == Gy)
	DiagBreak() ;
      for (count = 0, mean = 0.0, frame = 1 ; frame <= nimages ; frame++)
      {
	val = MRIgetVoxVal(mri_mosaic, x, y, 0, frame) ;
	vals[frame-1] = val ;
	if (val > 0)
	{
	  mean += val ;
	  count++ ;
	}
      }
      if (count == 0)
	continue ;  // no images at this location
      mean /= (double)count ;
      for (var = 0.0, frame = 0 ; frame < nimages ; frame++)
      {
	if (!FZERO(vals[frame]) && (vals[frame] > 0))
	{
	  dif = vals[frame]-mean ;
	  var += dif*dif ;
	}
      }

      var /= (double)count ;
      nonzero++ ;
      if (count > 1)
      {
	overlap_voxels++ ;
	energy += ((var/(double)count)) ;
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  min_overlap_voxels = nint(nonzero*.1) ;
  energy = (sqrt(energy/(double)overlap_voxels)) ;
//  printf("%d overlapping voxels (%2.3f%%) - energy = %2.4f\n", overlap_voxels, 100.0f*overlap_voxels/(nonzero), energy) ;
#if 0
  if (overlap_voxels < min_overlap_voxels)
    energy = 1e10 ;  // disallow certain configurations
#endif

  return(energy) ;
}

/*
  to compute voxel coords in mri2, add (dx,dy) to the mri1 coords

  integration will take place in space of MRI1 image
*/
static MRI *Gmri1, *Gmri2 ;
static double Gdx, Gdy ;

#define VALID(a,b,c,d)   (fabs(a)<.1 && fabs(b)<.1 && fabs(c)<.1 && fabs(1-d)<.1)
static double
compute_pairwise_deformation_energy(MRI *mri1, MRI *mri2, double dx, double dy, double *ax, double *ay, double a, double b, double c, double d)
{
  double  energy ;
  int     xmin, ymin, xmax, ymax, x1, nvox ;


  xmin = dx ; ymin = dy ; xmax = mri1->width-1 ; ymax = mri1->height-1;
  if (ymin < 0 || xmin < 0 || xmax >= mri1->width || ymax >= mri1->height || !devFinite(a) || !devFinite(b) || !devFinite(c) || !devFinite(d))
    return(1e10) ;

  if (fabs(Gdx-dx) > 20 || fabs(Gdy-dy)>20)
    return(1e10) ;

  if (!VALID(a,b,c,d))
    return(1e10) ;

  nvox = 0 ; energy = 0.0 ;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) reduction(+:energy,nvox)
#endif
  for (x1 = xmin ; x1 <= xmax ; x1++)
  {
    ROMP_PFLB_begin
    double  val1, val2, x1u, y1u, x2u, y2u ;
    int x2, y1, y2 ;

    x2 = x1 - dx ;
    if (x2 < 0 || x2 >= mri2->width)
      continue ;
    for (y1 = ymin ; y1 <= ymax ; y1++)
    {
      y2 = y1-dy ;
      if (y2 < 0 || y2 >= mri2->height)
	continue ;
      undistorted_coords(mri1, x1, y1, ax[0], ay[0], a, b, c, d, &x1u, &y1u) ;
      undistorted_coords(mri2, x2, y2, ax[1], ay[1], a, b, c, d, &x2u, &y2u) ;
      if ((MRIindexNotInVolume(mri1, x1u, y1u, 0) == 0) &&  (MRIindexNotInVolume(mri1, x1u, y1u, 0) == 0))
      {
	nvox++ ;
	MRIsampleVolumeType(mri1, x1u, y1u, 0, &val1, SAMPLE_TRILINEAR) ;
	MRIsampleVolumeType(mri2, x2u, y2u, 0, &val2, SAMPLE_TRILINEAR) ;
	energy += SQR(val1-val2) ;
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  if (nvox == 0 || (xmax-xmin ) < 100 || (ymax - ymin) < 100)
    return(1e10) ;
  energy = sqrt(energy/(double)nvox) ;
  if (energy < 7)
    DiagBreak() ;
  return(energy) ;
}


static int
undistorted_coords(MRI *mri, double x, double y, double ax, double ay, double a, double b, double c, double d, double *px, double *py)
{
  double rdst, x0, y0, r0, r2, r3, r4, theta, r, dx, dy ;

  x0 = (mri->width-1)/2.0 ; y0 = (mri->height-1)/2.0 ;
  r0 = MIN(x0,y0) ;
  dx = x-ax ; dy = y-ay ;
  r = sqrt(SQR(dx) + SQR(dy)) ;
  r /= r0 ;
  r2 = r*r ; r3 = r * r2 ; r4 = r*r3 ;
  rdst = a*r4 + b * r3 + c*r2 + d*r ;
  theta = atan2(dy, dx) ;
  rdst *= r0 ;
  *px = rdst*cos(theta) + ax ;
  *py = rdst*sin(theta) + ay ;
  return(NO_ERROR);
}
  
static int
powell_minimize(MRI *mri1, MRI *mri2, double dx, double dy, double a, double b, double c, double d, double *ax, double *ay,
		double *pa, double *pb, double *pc, double *pd, double *pdx, double *pdy) 

{
  float *p, **xi, fret, fstart, min_sse, orig_sse;
  int   iter, row, col ;

  Gmri1 = mri1 ;
  Gmri2 = mri2 ;
  Gdx = dx ; Gdy = dy ;

  p = vector(1, NPARMS) ;
  p[DX_I] = dx ; p[DY_I]  = dy; 
  if (NPARMS > 2)
  {
    p[A_I] = a ; p[B_I] = b ; p[C_I] = c ; p[D_I] = d ; 
    p[AX0_I] = ax[0] ; p[AY0_I] = ay[0] ;
    p[AX1_I] = ax[1] ; p[AY1_I] = ay[1] ; // optic axis
  }
  xi = matrix(1, NPARMS, 1, NPARMS) ;
  
  orig_sse = min_sse = compute_powell_sse(p) ;

  for (row = 1 ; row <= NPARMS ; row++)
  {
    for (col = 1 ; col <= NPARMS ; col++)
    {
      xi[row][col] = row == col ? 1 : 0 ;
    }
  }
  OpenPowell2(p, xi, NPARMS, intensity_tolerance, translation_tolerance, niters, &iter, &fret, compute_powell_sse);
  printf("best energy after powell: %2.3f (%d steps)\n", fret, iter) ;
  do
  {
    // reinitialize powell directions
    for (row = 1 ; row <= NPARMS ; row++)
    {
      for (col = 1 ; col <= NPARMS ; col++)
      {
        xi[row][col] = row == col ? 1 : 0 ;
      }
    }

    fstart = fret ;
    OpenPowell2(p, xi, NPARMS, intensity_tolerance, translation_tolerance, niters, &iter, &fret, compute_powell_sse);
//    OpenPowell(p, xi, NPARMS, translation_tolerance, &iter, &fret, compute_powell_sse);
    dx = p[DX_I] ; dy = p[DY_I] ; 
    if (NPARMS > 2)
    {
      a = p[A_I] ; b = p[B_I] ; c = p[C_I]  ; d = p[D_I] ; 
      ax[0] = p[AX0_I] ;  ay[0] = p[AY0_I] ; 
      ax[1] = p[AX1_I] ;  ay[1] = p[AY1_I] ; 
    }
    else
    {
      a = b = c = 0 ; d = 1 ;
    }
    printf("best energy after powell: %2.3f (%d steps)\n", fret, iter) ;
    if ((fstart-fret)/fstart < TOL)
      break ;
  }
  while (fret < fstart) ;

  min_sse = compute_powell_sse(p) ;
  *pa = a ; *pb = b ; *pc = c ; *pd = d ; *pdx = dx ; *pdy = dy ;
  free_matrix(xi, 1, NPARMS, 1, NPARMS) ;
  free_vector(p, 1, NPARMS) ;
  return((orig_sse-min_sse)/orig_sse < TOL) ;
}

static float
compute_powell_sse(float *p)
{
  double a, b, c, d, dx, dy, energy, ax[2], ay[2] ;

  dx = p[DX_I] ; dy = p[DY_I] ; 
  if (NPARMS > 2)
  {
    a = p[A_I] ; b = p[B_I]  ; c = p[C_I]  ; d = p[D_I] ; 
    ax[0] = p[AX0_I] ; ay[0] = p[AY0_I] ;
    ax[1] = p[AX1_I] ; ay[1] = p[AY1_I] ;
  }
  else
  {
    int i ;
    a = b = c = 0 ; d = 1 ;
    for (i = 0 ; i < 2 ; i++)
    {
      ax[i] = (Gmri1->width-1.0)/2.0 ;
      ay[i] = (Gmri1->height-1.0)/2.0 ;
    }
  }
  energy = compute_pairwise_deformation_energy(Gmri1, Gmri2, dx, dy, ax, ay, a, b, c, d) ;
  return(energy) ;
}
#if 0
static MRI *
add_image_to_mosaic(MRI *mri_orig_mosaic, MRI *mri, double x0, double y0, double ax, double ay, double a, double b, double c, double d)
{
  int      xmin, ymin, xmax, ymax ;
  MRI      *mri_mosaic ;
  int      width, height, x ;
  double   x1, y1, x2, y2 ;

  x1 = 0 ; y1 = 0 ; 
  x2 = MAX(mri_orig_mosaic->width-1, x0 + mri->width-1) ; 
  y2 = MAX(mri_orig_mosaic->height-1, y0 + mri->height-1) ;
  xmin = floor(x1) ; ymin = floor(y1) ;
  xmax = ceil(x2) ;  ymax = ceil(y2) ;
  
  width = xmax - xmin + 1 ;
  height = ymax - ymin + 1 ;
  mri_mosaic = MRIallocSequence(width, height, 1, MRI_FLOAT, 2) ;
  MRIcopyHeader(mri_orig_mosaic, mri_mosaic) ;
  MRIextractInto(mri_orig_mosaic, mri_mosaic, 0, 0, 0, 
		 mri_orig_mosaic->width,mri_orig_mosaic->height,mri_orig_mosaic->depth,0,0,0) ;

//#ifdef HAVE_OPENMP
//#pragma omp parallel for if_ROMP(experimental) 
//#endif
  {
    int x, y, count, xstart, xend, ystart, yend ;
    double   xu, yu, xf, yf ;
    double  val, sum ;

    // x and y are in the (larger) mosaic coords, while xf and yf are in the indidual tile coords
    xstart = MAX(0,x0-xmin) ;
    xend = MIN(mri->width-1+x0-xmin, mri_mosaic->width-1) ;
    ystart = MAX(0,y0-ymin) ;
    yend = MIN(mri->height-1+y0-ymin, mri_mosaic->height-1) ;
    for (x = xstart ; x <= xend ; x++)
    {
      xf = x - x0 + xmin ;
      for (y = ystart ; y <= yend ; y++)
      {
	if ( x== Gx && y == Gy)
	  DiagBreak() ;
	yf = y - y0 + ymin ;
	count = (int)MRIgetVoxVal(mri_mosaic, x, y, 0, 1) ;
	sum = MRIgetVoxVal(mri_mosaic, x, y, 0, 0) ;
	count++ ;
	MRIsetVoxVal(mri_mosaic, x, y, 0, 1, count) ;

	undistorted_coords(mri, xf, yf, ax, ay, a, b, c, d, &xu, &yu) ;
	MRIsampleVolumeType(mri, xu, yu, 0, &val, SAMPLE_TRILINEAR) ;

	MRIsetVoxVal(mri_mosaic, x, y, 0, 0, val+sum) ;
//	MRIsetVoxVal(mri_mosaic, x, y, 0, 1, val) ;
      }
    }
  }

  MRIwrite(mri_mosaic, "m.mgz");

#ifdef HAVE_OPENMP
#pragma omp parallel for if_ROMP(experimental) 
#endif
  for (x = 0 ; x < mri_mosaic->width ; x++)
  {
    float sum, count ;
    int   y ;
    for (y = 0 ; y < mri_mosaic->height ; y++)
    {
      if (x == Gx && y == Gy)
	DiagBreak() ;
      count = MRIgetVoxVal(mri_mosaic, x, y, 0, 1) ;
      if (count <= 1)
	continue ;
      sum = MRIgetVoxVal(mri_mosaic, x, y, 0, 0) ;
      MRIsetVoxVal(mri_mosaic, x, y, 0, 0, sum/count) ;
    }
  }

  return(mri_mosaic) ;
}
#endif

static MRI *
undistort_and_mosaic_images(MRI **mri, double *x0, double *y0, int nimages, double *ax, double *ay, double a, double b, double c, double d, double *penergy) 
{
  int      xmin, ymin, xmax, ymax, nvoxels ;
  MRI      *mri_mosaic, *mri_counts, *mri_vars ;
  int      i, width, height, x ;
  MRI      *mri_jac ;
  double   x1, y1, x2, y2, energy ;

  ymin = xmin = 10000000 ;
  ymax = xmax = -ymin ;

  for (i =  0 ; i < nimages ; i++)
  {
    x1 = x0[i] ; y1 = y0[i] ; x2 = x1 + mri[i]->width-1 ; y2 = y1 + mri[i]->height-1 ;
    if (x1 < xmin)
      xmin = floor(x1) ;
    if (y1 < ymin)
      ymin = floor(y1) ;
    if (x2 > xmax)
      xmax = ceil(x2) ;
    if (y2 > ymax)
      ymax = ceil(y2) ;
  }
  
  width = xmax - xmin + 1 ;
  height = ymax - ymin + 1 ;
  mri_mosaic = MRIallocSequence(width, height, 1, MRI_FLOAT, nimages+1) ;
  mri_jac = MRIallocSequence(mri[0]->width, mri[0]->height, 1, MRI_FLOAT, nimages) ;
  mri_vars = MRIallocSequence(width, height, 1, MRI_FLOAT, 1) ;
  MRIcopyHeader(mri[0], mri_mosaic) ;
  MRIcopyHeader(mri[0], mri_jac) ;
  MRIcopyHeader(mri_mosaic, mri_vars) ;

  mri_counts = MRIclone(mri_mosaic, NULL) ;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (i =  0 ; i < nimages ; i++)
  {
    ROMP_PFLB_begin
    int x, y, count, xstart, xend, ystart, yend ;
    double   xu, yu, xf, yf ;
    double  val, sum ;

    // x and y are in the (larger) mosaic coords, while xf and yf are in the indidual tile coords
    xstart = MAX(0,x0[i]-xmin) ;
    xend = MIN(mri[i]->width-1+x0[i]-xmin, mri_mosaic->width-1) ;
    ystart = MAX(0,y0[i]-ymin) ;
    yend = MIN(mri[i]->height-1+y0[i]-ymin, mri_mosaic->height-1) ;
    for (x = xstart ; x <= xend ; x++)
    {
      xf = x - x0[i] + xmin ;
      for (y = ystart ; y <= yend ; y++)
      {
	if ( x== Gx && y == Gy)
	  DiagBreak() ;
	yf = y - y0[i] + ymin ;
//	undistorted_coords(mri[i], xf, yf, ax[i], ay[i], a, b, c, d, &xu, &yu) ;
	xu = xf ; yu = yf ;  //disable unwarping for now
	if (MRIindexNotInVolume(mri[i], xu, yu, 0))
	  continue ;
	MRIsampleVolumeType(mri[i], xu, yu, 0, &val, SAMPLE_TRILINEAR) ;
	count = (int)MRIgetVoxVal(mri_counts, x, y, 0, 0) ;
	count++ ;
	MRIsetVoxVal(mri_counts, x, y, 0, 0, count) ;
	sum = MRIgetVoxVal(mri_mosaic, x, y, 0, 0) ;
	
#if 0
	{
	  double xp1, yp1, jac;
	  undistorted_coords(mri[i], x+1, y+1, ax[i], ay[i], a, b, c, d, &xp1, &yp1) ;
	  jac = fabs((xp1-xu)*(yp1-yu)) ;
	  if (jac > 10)
	    DiagBreak() ;
	  MRIsetVoxVal(mri_jac, x, y, 0, i, jac) ;
	  val /= jac ;
	}
#endif

	MRIsetVoxVal(mri_mosaic, x, y, 0, 0, val+sum) ;
	MRIsetVoxVal(mri_mosaic, x, y, 0, i+1, val) ;
	sum = MRIgetVoxVal(mri_vars, x, y, 0, 0) ;
	MRIsetVoxVal(mri_vars, x, y, 0, 0, sum+val*val) ;
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

//  MRIwrite(mri_jac, "jac.mgz") ;
  MRIfree(&mri_jac) ;
  energy = 0.0 ; nvoxels = 0 ;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) reduction (+:energy,nvoxels)
#endif
  for (x = 0 ; x < mri_mosaic->width ; x++)
  {
    ROMP_PFLB_begin
    
    float sum, count, var ;
    int   y ;
    for (y = 0 ; y < mri_mosaic->height ; y++)
    {
      count = MRIgetVoxVal(mri_counts, x, y, 0, 0) ;
      if (count <= 1)
	continue ;
      sum = MRIgetVoxVal(mri_mosaic, x, y, 0, 0) ;
      sum /= count ;
      MRIsetVoxVal(mri_mosaic, x, y, 0, 0, sum) ;
      var = MRIgetVoxVal(mri_vars, x, y, 0, 0) ;
      var = var / count - sum*sum ;
      energy += var ;
      nvoxels++ ;
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  energy = sqrt(energy / nvoxels) ;
  if (penergy)
    *penergy = energy ;
  MRIfree(&mri_counts) ;
  MRIfree(&mri_vars) ;
  return(mri_mosaic) ;
}

#if 0
static int
jacobian_correct(MRI *mri_src, MRI *mri_dst, double ax, double ay, double a, double b, double c, double d) 
{
  int x, y ;
  
  for (x = 0 ; x < mri_src->width ; x++)
  {
    double xp1, yp1, jac, xu, yu, val;

    for (y = 0 ; y < mri_src->height ; y++)
    {
      if ( x== Gx && y == Gy)
	DiagBreak() ;
      
      undistorted_coords(mri_src, x, y, ax, ay, a, b, c, d, &xu, &yu) ;
      val = MRIgetVoxVal(mri_src, x, y, 0, 0) ;
      undistorted_coords(mri_src, x+1, y+1, ax, ay, a, b, c, d, &xp1, &yp1) ;
      jac = fabs((xp1-xu)*(yp1-yu)) ;
      if (jac > 10)
	DiagBreak() ;
      val /= jac ;
      MRIsetVoxVal(mri_dst, x, y, 0, 0, val) ;
    }
  }
  return(NO_ERROR) ;
}
#endif

static int
insert_image_into_mosaic(MRI *mri_mosaic, MRI *mri, double x0, double y0)
{
  int   x, y ;
  double xi, yi, dof, vali, valm ;

  for (x = nint(x0) ; x < nint(x0)+mri->width; x++)
    for (y = nint(y0) ; y < nint(y0)+mri->height; y++)
      {
	xi = x - x0 ; yi = y - y0 ; 
	dof = MRIgetVoxVal(mri_mosaic, x, y, 0, 1) ;
	valm = MRIgetVoxVal(mri_mosaic, x, y, 0, 0) ;
	vali = MRIgetVoxVal(mri, xi, yi, 0, 0) ;
	valm = (dof * valm + vali) / (dof+1) ;
	dof += 1 ;
	MRIsetVoxVal(mri_mosaic, x, y, 0, 0, valm) ;
	MRIsetVoxVal(mri_mosaic, x, y, 0, 1, dof) ;
      }
  return(NO_ERROR) ;
}
static int Gnimages ;
static MRI **Gmri ;
static MRI **Gmri_smooth ;
static double Gx0 , Gy0 ;  // offset of first image that won't be modified by Powell

extern void (*powell_iteration_func)(float *p, int nparms) ;

static void
dump_parms(float *p, int nimages, const char *name) 
{
  int i ;

  printf("%s:\n",name) ;
  printf("image %2.2d: (%2.1f, %2.1f)\n", 0, Gx0, Gy0) ;
  for (i = 1 ; i < nimages ; i++)
    printf("image %2.2d: (%2.1f, %2.1f)\n", i, p[i], p[nimages-1+i]) ;
    
}
static void
powell_step_func(float *p, int nparms) 
{
  static int step = 0 ;
  int  i, nimages ;
  char fname[STRLEN] ;
  MRI *mri_mosaic ;
  double x0d[MAX_IMAGES], y0d[MAX_IMAGES] ;

  nimages = (nparms+2)/2 ;
  int req = snprintf(fname, STRLEN, "%s powell iter %d", Gout_name, step) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  dump_parms(p, nimages, fname) ;
  x0d[0] = Gx0 ; y0d[0] = Gy0 ;
  for (i = 1 ; i < Gnimages ; i++)
  {
    x0d[i] = p[i] ;
    y0d[i] = p[Gnimages-1+i] ;
  }
  mri_mosaic = undistort_and_mosaic_images(Gmri, x0d, y0d, nimages, 0, 0, 0, 0, 0, 1, NULL) ;
  req = snprintf(fname, STRLEN, "%s.powell.%2.2d.mgz", Gout_name, step) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("writing image after powell step %d to %s\n", step, fname) ;
  MRIwrite(mri_mosaic,fname) ;

  step++;
  MRIfree(&mri_mosaic) ;
}

static MRI *
build_best_mosaic(MRI **mri, MRI **mri_smooth, double *x0d, double *y0d, int nimages) 
{
  float *p, **xi, fret, fstart, min_sse, orig_sse;
  int   iter, row, col, nparms, i ;
  MRI   *mri_mosaic ;

  powell_iteration_func = powell_step_func ;
  Gnimages = nimages ;
  Gmri = mri ;
  Gmri_smooth = mri_smooth ;

  nparms = nimages*2-2 ;     // dx and dy for each image, but first image is fixed
  p = vector(1, nparms) ;   
  xi = matrix(1, nparms, 1, nparms) ;

  Gx0 = x0d[0] ; Gy0 = y0d[0] ;
  for (i = 1 ; i < Gnimages ; i++)
  {
    p[i] = x0d[i] ;
    p[Gnimages-1+i] = y0d[i] ;   // [0] slots are used to store first image offset but won't be modified
  }
  dump_parms(p, nimages, "original") ;
  orig_sse = min_sse = compute_powell_global_sse(p) ;

  for (row = 1 ; row <= nparms ; row++)
  {
    for (col = 1 ; col <= nparms ; col++)
    {
      xi[row][col] = row == col ? 1 : 0 ;
    }
  }

  OpenPowell2(p, xi, nparms, intensity_tolerance, translation_tolerance, niters, &iter, &fret, compute_powell_global_sse);
  dump_parms(p, nimages, "after powell") ;

  do
  {
    // reinitialize powell directions
    for (row = 1 ; row <= nparms ; row++)
    {
      for (col = 1 ; col <= nparms ; col++)
      {
        xi[row][col] = row == col ? 1 : 0 ;
      }
    }

    fstart = fret ;
    OpenPowell2(p, xi, nparms, intensity_tolerance, translation_tolerance, niters, &iter, &fret, compute_powell_global_sse);
    printf("best energy after powell: %2.3f (%d steps)\n", fret, iter) ;
    dump_parms(p, nimages, "after second powell") ;
  }
  while ((100*(fstart-fret)/fstart) < intensity_tolerance) ;

  for (i = 1 ; i < Gnimages ; i++)
  {
    x0d[i] = p[i] ;
    y0d[i] = p[Gnimages-1+i] ;
  }
  x0d[0] = Gx0 ; y0d[0] = Gy0 ;
  mri_mosaic = undistort_and_mosaic_images(mri, x0d, y0d, nimages, 0, 0, 0, 0, 0, 1, NULL) ;

  min_sse = compute_powell_global_sse(p) ;
  free_matrix(xi, 1, nparms, 1, nparms) ;
  free_vector(p, 1, nparms) ;
  return(mri_mosaic) ;
}

static float
compute_powell_global_sse(float *p)
{
  double    x0[MAX_IMAGES], y0[MAX_IMAGES], energy ;
  int       i ;
  MRI       *mri_mosaic ;
  
  for (i = 1 ; i < Gnimages ; i++)
  {
    x0[i] = p[i] ;
    y0[i] = p[Gnimages-1+i] ;
  }
  x0[0] = Gx0 ; y0[0] = Gy0 ;
  mri_mosaic = undistort_and_mosaic_images(Gmri_smooth, x0, y0, Gnimages, 0, 0, 0, 0, 0, 1, &energy) ;
//  energy = compute_mosaic_energy(mri_mosaic, Gnimages) ;
  MRIfree(&mri_mosaic) ;
  return(energy) ;
}

static MRI *
mosaic_images_with_xforms(MRI **mri, int nimages, MRI *mri_weights, int pad) 
{
  double   xmin, ymin, zmin, xmax, ymax, zmax, x1, y1, z1 ;
  MRI      *mri_mosaic, *mri_counts ;
  int      i, width, height, depth, x ;
  MATRIX   *m_vox2ras[MAX_IMAGES], *m_ras2vox ;
  VECTOR   *v_ras, *v_vox ;

  zmin = ymin = xmin = 10000000 ;
  zmax = ymax = xmax = -ymin ;

  v_ras = VectorAlloc(4, MATRIX_REAL) ;
  v_vox = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v_ras, 4) = VECTOR_ELT(v_vox, 4) = 1.0 ;
  for (i =  0 ; i < nimages ; i++)
  {
    int x, y, z ;
    m_vox2ras[i] = MRIgetVoxelToRasXform(mri[i]) ;

    for (x = 0 ; x < mri[i]->width ; x++)
    {
      for (y = 0 ; y < mri[i]->height ; y++)
      {
	for (z = 0 ; z < mri[i]->depth ; z++)
	{
	  if (x == Gx && y == Gy && z == Gz)
	    DiagBreak() ;
	  if (i == Gdiag_no && x == Gx && y == Gy && z == Gz)
	    DiagBreak() ;

	  V3_X(v_vox) = x ; V3_Y(v_vox) = y ; V3_Z(v_vox) = z ;
	  MatrixMultiply(m_vox2ras[i], v_vox, v_ras) ;
	  x1 = V3_X(v_ras) ; y1 = V3_Y(v_ras) ; z1 = V3_Z(v_ras) ;

	  if (x1-pad < xmin)
	    xmin = x1 - pad;
	  if (y1-pad < ymin)
	    ymin = y1-pad ;
	  if (z1-pad < zmin)
	    zmin = z1-pad ;

	  if (x1+pad > xmax)
	    xmax = x1+pad ;
	  if (y1+pad > ymax)
	    ymax = y1 + pad ;
	  if (z1+pad > zmax)
	    zmax = z1+pad ;
	}
      }
    }
  }
  
  width =  ceil((xmax - xmin ) / mri[0]->xsize) + 1 ;
  height = ceil((ymax - ymin ) / mri[0]->ysize) + 1 ;
  depth =  ceil((zmax - zmin ) / mri[0]->zsize) + 1 ;
  mri_mosaic = MRIallocSequence(width, height, depth, MRI_FLOAT, 1) ;
  MRIcopyHeader(mri[0], mri_mosaic) ;

  //  compute vox2ras xform by mapping corners of RAS rectangle to corners of image
  {
    MATRIX  *m_vox, *m_ras, *m_vox2ras, *m_pinv ;

    m_vox = MatrixAlloc(4, 8, MATRIX_REAL) ;
    m_ras = MatrixAlloc(4, 8, MATRIX_REAL) ;

    // (0,0,0)
    *MATRIX_RELT(m_vox, 1, 1) = 0 ; 
    *MATRIX_RELT(m_vox, 2, 1) = 0 ;
    *MATRIX_RELT(m_vox, 3, 1) = 0 ; 
    *MATRIX_RELT(m_vox, 4, 1) = 1 ;

    *MATRIX_RELT(m_ras, 1, 1) = xmin ; 
    *MATRIX_RELT(m_ras, 2, 1) = ymin ;
    *MATRIX_RELT(m_ras, 3, 1) = zmin ; 
    *MATRIX_RELT(m_ras, 4, 1) = 1 ;

    // (1, 1, 1)
    *MATRIX_RELT(m_vox, 1, 2) = width-1 ; 
    *MATRIX_RELT(m_vox, 2, 2) = height-1 ;
    *MATRIX_RELT(m_vox, 3, 2) = depth-1 ; 
    *MATRIX_RELT(m_vox, 4, 2) = 1 ;

    *MATRIX_RELT(m_ras, 1, 2) = xmax ; 
    *MATRIX_RELT(m_ras, 2, 2) = ymax ;
    *MATRIX_RELT(m_ras, 3, 2) = zmax ; 
    *MATRIX_RELT(m_ras, 4, 2) = 1 ;

    // (0, 1, 1)
    *MATRIX_RELT(m_vox, 1, 3) = 0 ; 
    *MATRIX_RELT(m_vox, 2, 3) = height-1 ;
    *MATRIX_RELT(m_vox, 3, 3) = depth-1 ; 
    *MATRIX_RELT(m_vox, 4, 3) = 1 ;

    *MATRIX_RELT(m_ras, 1, 3) = xmin ; 
    *MATRIX_RELT(m_ras, 2, 3) = ymax ;
    *MATRIX_RELT(m_ras, 3, 3) = zmax ; 
    *MATRIX_RELT(m_ras, 4, 3) = 1 ;

    // (1, 0, 1)
    *MATRIX_RELT(m_vox, 1, 4) = width-1 ; 
    *MATRIX_RELT(m_vox, 2, 4) = 0 ;
    *MATRIX_RELT(m_vox, 3, 4) = depth-1 ; 
    *MATRIX_RELT(m_vox, 4, 4) = 1 ;

    *MATRIX_RELT(m_ras, 1, 4) = xmax ; 
    *MATRIX_RELT(m_ras, 2, 4) = ymin ;
    *MATRIX_RELT(m_ras, 3, 4) = zmax ; 
    *MATRIX_RELT(m_ras, 4, 4) = 1 ;

    // (1, 1, 0)
    *MATRIX_RELT(m_vox, 1, 5) = width-1 ; 
    *MATRIX_RELT(m_vox, 2, 5) = height-1 ;
    *MATRIX_RELT(m_vox, 3, 5) = 0 ; 
    *MATRIX_RELT(m_vox, 4, 5) = 1 ;

    *MATRIX_RELT(m_ras, 1, 5) = xmax ; 
    *MATRIX_RELT(m_ras, 2, 5) = ymax ;
    *MATRIX_RELT(m_ras, 3, 5) = zmin ; 
    *MATRIX_RELT(m_ras, 4, 5) = 1 ;

    // (1, 0, 0)
    *MATRIX_RELT(m_vox, 1, 6) = width-1 ; 
    *MATRIX_RELT(m_vox, 2, 6) = 0 ;
    *MATRIX_RELT(m_vox, 3, 6) = 0 ; 
    *MATRIX_RELT(m_vox, 4, 6) = 1 ;

    *MATRIX_RELT(m_ras, 1, 6) = xmax ; 
    *MATRIX_RELT(m_ras, 2, 6) = ymin ;
    *MATRIX_RELT(m_ras, 3, 6) = zmin ; 
    *MATRIX_RELT(m_ras, 4, 6) = 1 ;

    // (0, 0, 1)
    *MATRIX_RELT(m_vox, 1, 7) = 0 ; 
    *MATRIX_RELT(m_vox, 2, 7) = 0 ;
    *MATRIX_RELT(m_vox, 3, 7) = depth-1 ; 
    *MATRIX_RELT(m_vox, 4, 7) = 1 ;

    *MATRIX_RELT(m_ras, 1, 7) = xmin ; 
    *MATRIX_RELT(m_ras, 2, 7) = ymin ;
    *MATRIX_RELT(m_ras, 3, 7) = zmax ; 
    *MATRIX_RELT(m_ras, 4, 7) = 1 ;

    // (0, 1, 0)
    *MATRIX_RELT(m_vox, 1, 8) = width-1 ; 
    *MATRIX_RELT(m_vox, 2, 8) = 0 ;
    *MATRIX_RELT(m_vox, 3, 8) = depth-1 ; 
    *MATRIX_RELT(m_vox, 4, 8) = 1 ;

    *MATRIX_RELT(m_ras, 1, 8) = xmax ; 
    *MATRIX_RELT(m_ras, 2, 8) = ymin ;
    *MATRIX_RELT(m_ras, 3, 8) = zmax ; 
    *MATRIX_RELT(m_ras, 4, 8) = 1 ;

    m_pinv = MatrixPseudoInverse(m_vox, NULL) ;
    if (m_pinv == NULL)
    {
      MatrixPrint(stderr, m_vox) ;
      MatrixPrint(stderr, m_ras) ;
      ErrorExit(ERROR_BADPARM, "%s: could not invert vox coord matrix\n", Progname) ;
    }
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    {
      MatrixPrint(stderr, m_vox) ;
      MatrixPrint(stderr, m_ras) ;
    }
    m_vox2ras = MatrixMultiply(m_ras, m_pinv, NULL) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    {
      printf("setting vox2ras of mosaic to:\n") ;
      MatrixPrint(stdout, m_vox2ras) ;
    }
    MRIsetVox2RASFromMatrix(mri_mosaic, m_vox2ras) ;

    MatrixFree(&m_ras) ; MatrixFree(&m_vox) ; MatrixFree(&m_pinv) ; MatrixFree(&m_vox2ras) ;
  }

  printf("mosaic size = (%d x %d x %d) : (%2.2f-->%2.2f, %2.2f-->%2.2f, %2.2f-->%2.2f \n", 
	 width, height, depth, xmin, xmax, ymin, ymax, zmin, zmax) ;

  m_ras2vox = MRIgetRasToVoxelXform(mri_mosaic) ;
  mri_counts = MRIclone(mri_mosaic, NULL) ;
//#ifdef HAVE_OPENMP
//#pragma omp parallel for if_ROMP(experimental) 
//#endif
  for (i =  0 ; i < nimages ; i++)
  {
    int x, y, z, xm, ym, zm ;
    float  val, sum, count, w ;
    MATRIX *m_vox2vox ;
    VECTOR *v_vox2 ;

    v_vox2 = VectorAlloc(4, MATRIX_REAL) ; VECTOR_ELT(v_vox2, 4) = 1.0 ;
    m_vox2vox = MRIgetVoxelToVoxelXform(mri[i], mri_mosaic) ;
    for (x = 0 ; x < mri[i]->width ; x++)
    {
      for (y = 0 ; y < mri[i]->height ; y++)
      {
	for (z = 0 ; z < mri[i]->depth ; z++)
	{
	  V3_X(v_vox) = x ; V3_Y(v_vox) = y ; V3_Z(v_vox) = z ;
	  MatrixMultiply(m_vox2vox, v_vox, v_vox2) ;
	  if ( x== Gx && y == Gy && z == Gz)
	    DiagBreak() ;
	  xm = nint(V3_X(v_vox2)) ; ym = nint(V3_Y(v_vox2)) ; zm = nint(V3_Z(v_vox2)) ;

	  if (xm == Gx && ym == Gy && zm == Gz)
	    DiagBreak() ;
	  if (MRIindexNotInVolume(mri_counts, xm, ym, zm) != 0)
	  {
	    static int first = 1 ;

	    if (first)
	    {
	      MatrixMultiply(m_vox2ras[i], v_vox, v_ras) ;
	      x1 = V3_X(v_ras) ; y1 = V3_Y(v_ras) ; z1 = V3_Z(v_ras) ;
	      printf("image %d: voxel coord (%d %d %d) -> (%d %d %d), RAS (%2.2f, %2.2f %2.2f) OOB!\n",
		     i, x, y, z, xm, ym, zm, x1, y1, z1) ;
	    }
	    first = 0 ;
	    continue ;
	  }
	  count = MRIgetVoxVal(mri_counts, xm, ym, zm, 0) ;
	  if (mri_weights)
	    w = MRIgetVoxVal(mri_weights, x, y, z, 0) ; /* use confocal weighting */
	  else
	    w = 1 ;

	  count = count + w ;
	  if (count > 1)
	    DiagBreak() ;
	  MRIsetVoxVal(mri_counts, xm, ym, zm, 0, count) ;
	  sum = MRIgetVoxVal(mri_mosaic, xm, ym, zm, 0) ;
	  val = MRIgetVoxVal(mri[i], x, y, z, 0) ;
	  MRIsetVoxVal(mri_mosaic, xm, ym, zm, 0, w*val+sum) ;
//	  MRIsetVoxVal(mri_mosaic, xm, ym, zm, i+1, val) ;
	}
      }
    }
    VectorFree(&v_vox2) ; MatrixFree(&m_vox2vox) ;
  }
  if (Gdiag & DIAG_WRITE)
  {  MRIwrite(mri_mosaic, "m.mgz") ; MRIwrite(mri_counts, "c.mgz") ;}

//#ifdef HAVE_OPENMP
//#pragma omp parallel for if_ROMP(experimental) 
//#endif
  for (x = 0 ; x < mri_mosaic->width ; x++)
  {
    float sum, count ;
    int   y, z ;
    for (y = 0 ; y < mri_mosaic->height ; y++)
    {
      for (z = 0 ; z < mri_mosaic->depth ; z++)
      {
	if ( x== Gx && y == Gy && z == Gz)
	  DiagBreak() ;
	count = MRIgetVoxVal(mri_counts, x, y, z, 0) ;
	if (count <= 1)
	  continue ;
	sum = MRIgetVoxVal(mri_mosaic, x, y, z, 0) ;
	MRIsetVoxVal(mri_mosaic, x, y, z, 0, sum/count) ;
      }
    }
  }

  if (Gdiag & DIAG_WRITE)
    MRIwrite(mri_mosaic, "m2.mgz") ; 

  for (i = 0 ; i < nimages ; i++)
    MatrixFree(&m_vox2ras[i]) ;
  VectorFree(&v_ras) ; VectorFree(&v_vox) ;
  MatrixFree(&m_ras2vox) ;
  MRIfree(&mri_counts) ;
  return(mri_mosaic) ;
}


