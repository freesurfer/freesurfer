#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "diag.h"
#include "error.h"
#include "macros.h"
#include "utils.h"
#include "proto.h"
#include "classify.h"
#include "mri.h"

static int extract = 0 ;
static int  verbose = 0 ;
static int wsize = 11 ;
static float pct = 0.8 ;
static float pslope = 1.0f ;
static float nslope = 1.0f ;
static float wm_low = 90 ;
static float wm_hi = 125 ;
static float gray_hi = 100 ;
static int niter = 1 ;

char *Progname ;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
MRI *MRIremoveWrongDirection(MRI *mri_src, MRI *mri_dst, int wsize,
                             float low_thresh, float hi_thresh) ;

int
main(int argc, char *argv[])
{
  MRI     *mri_src, *mri_dst, *mri_tmp, *mri_labeled ;
  char    *input_file_name, *output_file_name ;
  int     nargs, i ;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    ErrorExit(ERROR_BADPARM,
              "usage: %s <input volume> <output volume>", Progname);

  input_file_name = argv[1] ;
  output_file_name = argv[2] ;

  mri_src = MRIread(input_file_name) ;
  if (!mri_src)
    ErrorExit(ERROR_NOFILE, "%s: could not read source volume from %s",
              Progname, input_file_name) ;

  fprintf(stderr, "doing initial intensity segmentation...\n") ;
  mri_tmp = MRIintensitySegmentation(mri_src, NULL, wm_low, wm_hi, gray_hi);

  fprintf(stderr, "using local statistics to label ambiguous voxels...\n") ;
  MRIhistoSegment(mri_src, mri_tmp, wm_low, wm_hi, gray_hi, wsize, 3.0f) ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_tmp, "/tmp/int.mnc") ;
  fprintf(stderr, 
          "using local geometry to label remaining ambiguous voxels...\n") ;
  mri_labeled = MRIcpolvMedianCurveSegment(mri_src, mri_tmp, NULL, 5, 3);
  fprintf(stderr, 
          "\nreclassifying voxels using Gaussian border classifier...\n") ;

  /*
    now use the gray and white mattter border voxels to build a Gaussian
    classifier at each point in space and reclassify all voxels in the
    range [wm_low-5,gray_hi+5].
    */
  for (i = 0 ; i < niter ; i++)
  {
    MRIreclassify(mri_src, mri_labeled, mri_labeled, wm_low-5,gray_hi+5,wsize);
  }
  MRIfree(&mri_tmp) ;

  mri_dst = MRImaskLabels(mri_src, mri_labeled, NULL) ;
  fprintf(stderr, 
          "\nremoving voxels with positive offset direction...\n") ;
  MRIremoveWrongDirection(mri_dst, mri_dst, 3, wm_low-5, gray_hi+5) ;

  MRIfree(&mri_src) ;
  MRIfree(&mri_labeled) ; 
  fprintf(stderr, "writing output to %s...", output_file_name) ;
  MRIwrite(mri_dst, output_file_name) ;
  fprintf(stderr, "done.\n") ;

  MRIfree(&mri_dst) ;

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
  if (!stricmp(option, "slope"))
  {
    nslope = pslope = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using curvature slope = %2.2f\n", pslope) ;
  }
  else if (!stricmp(option, "pslope"))
  {
    pslope = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using curvature pslope = %2.2f\n", pslope) ;
  }
  else if (!stricmp(option, "nslope"))
  {
    nslope = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using curvature nslope = %2.2f\n", nslope) ;
  }
  else switch (toupper(*option))
  {
  case 'N':
    niter = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "running border classification %d times\n", niter) ;
    break ;
  case 'V':
    verbose = !verbose ;
    break ;
  case 'P':
    pct = atof(argv[2]) ;
    fprintf(stderr, "using %2.0f%% threshold\n", pct*100.0f) ;
    break ;
  case 'X':
    if (sscanf(argv[2], "%d", &extract) != 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                Progname, argv[2]) ;
    nargs = 1 ;
    break ;
  case 'W':
    wsize = atoi(argv[2]) ;
    fprintf(stderr, "using wsize = %d\n", wsize) ;
    break ;
  case '?':
  case 'U':
    fprintf(stderr,
            "usage: %s <classifier file> <input volume> <output volume>\n",
            Progname) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
#define BLUR_SIGMA 0.25f
MRI *
MRIremoveWrongDirection(MRI *mri_src, MRI *mri_dst, int wsize, 
                        float low_thresh, float hi_thresh)
{
  MRI  *mri_kernel, *mri_smooth ;
  int  x, y, z, width, height, depth, val, nchanged, ntested ;

  mri_kernel = MRIgaussian1d(BLUR_SIGMA, 100) ;
  fprintf(stderr, "smoothing T1 volume with sigma = %2.3f\n", BLUR_SIGMA) ;
  mri_smooth = MRIclone(mri_src, NULL) ;
  MRIconvolveGaussian(mri_src, mri_smooth, mri_kernel) ;
  MRIfree(&mri_kernel) ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  nchanged = ntested = 0 ;
  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ; 
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        val = MRIvox(mri_src, x, y, z) ;
        if (val >= low_thresh && val <= hi_thresh)
        {
          ntested++ ;
          if (MRIvoxelDirection(mri_smooth, x, y, z, wsize) > 0)
          {
            nchanged++ ;
            val = 0 ;
          }
        }
        MRIvox(mri_dst, x, y, z) = val ;
      }
    }
  }
    
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stderr, "               %8d voxels tested (%2.2f%%)\n",
            ntested, 100.0f*(float)ntested/ (float)(width*height*depth));
    fprintf(stderr, "               %8d voxels changed (%2.2f%%)\n",
            nchanged, 100.0f*(float)nchanged/ (float)(width*height*depth));
  }
  MRIfree(&mri_smooth) ;
  return(mri_dst) ;
}




