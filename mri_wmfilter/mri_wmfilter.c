 #include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <memory.h>

#include "filter.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "mri.h"
#include "region.h"

static char vcid[] = "$Id: mri_wmfilter.c,v 1.1 1997/07/23 03:41:01 fischl Exp $";

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

#define FILTER_WSIZE         5
#define DEFAULT_FRACTION_ON  0.6f

char *Progname ;

static int   filter_window_size = FILTER_WSIZE ;

static float frac_on = DEFAULT_FRACTION_ON ;
static float white_hilim = 100;
static float white_lolim = 60; 
static float gray_hilim = 70; 
static float threshold = 30;

static void read_parameter_file(char *fname) ;

int
main(int argc, char *argv[])
{
  char        **av ;
  int         ac, nargs, num_above_threshold ;
  char        *in_fname, *out_fname ;
  MRI         *mri_full, *mri_src, *mri_dst, *mri_cpolv, *mri_cpolv_count ;
  MRI_REGION  clip_region ;

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

  if (argc < 3)
    usage_exit() ;

  num_above_threshold = nint((float)filter_window_size * frac_on) ;
  read_parameter_file("wmfilter.dat") ;
  in_fname = argv[1] ;
  out_fname = argv[2] ;

  mri_full = MRIread(in_fname) ;
  if (!mri_full)
    ErrorExit(ERROR_NOFILE, "%s: could not read '%s'", Progname, in_fname) ;

  /* only work on the region with actual data */
  MRIboundingBox(mri_full, threshold, &clip_region) ;
  REGIONexpand(&clip_region, &clip_region, (filter_window_size+1)/2) ;
  mri_src = MRIextractRegion(mri_full, NULL, &clip_region) ;

  /* compute the plane of least variance */
  mri_cpolv = MRIplaneOfLeastVarianceNormal(mri_src, NULL, filter_window_size);
MRIwrite(mri_cpolv, "cpolv.mnc") ;

  /* find only those pixels that lie in with white-like intensity values */
  mri_cpolv_count = MRIpolvCount(mri_src, NULL, mri_cpolv, 
                                 filter_window_size, white_lolim, white_hilim);
MRIwrite(mri_cpolv_count, "order.mnc") ;
  mri_dst = MRIorderThreshold(mri_src,NULL,mri_cpolv_count,
                              num_above_threshold);

  /* also include other values that are in the correct range */
  MRIthresholdRangeInto(mri_src, mri_dst, gray_hilim+1, white_hilim) ;
  if (!mri_dst)
    ErrorExit(ERROR_NOFILE, "%s: could allocate space for destination image", 
              Progname) ;

  /* now put it back in the full size image */
  MRIextractIntoRegion(mri_dst, mri_full, 0, 0, 0, &clip_region);
  MRIfree(&mri_dst) ;
  MRIwrite(mri_full, out_fname) ;
  MRIfree(&mri_full) ;
  if (DIAG_VERBOSE_ON && (Gdiag & DIAG_SHOW))
    fprintf(stderr, "done.\n") ;
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
  else switch (toupper(*option))
  {
  case 'W':
    if (sscanf(argv[2], "%d", &filter_window_size) < 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan window size from '%s'",
                Progname, argv[2]) ;
    if (filter_window_size < 3)
      ErrorExit(ERROR_BADPARM, "%s: filter window size must be >= 3",Progname);
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <input image file> <output image file>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
#if 0
  fprintf(stderr, 
          "\nThis program will process the image contained in "
          "<input image file>\n"
          "using a nonlocal filter, and write the results to "
          "<output image file>.\n"
          "\nThe default filter is the median, but either Gaussian or mean "
          "filtering\n"
          "can be specified through command line options (see below).\n"
          "\nBy default the image is smoothed using a Gaussian kernel "
          "(sigma=%2.3f)\n"
          "before the offset field is calculated. This can be modified using\n"
          "the -blur command line option.\n"
          "\nThe input and output image formats are specified by the file name"
          " extensions.\n"
          "Supported formats are:\n\n"
          "   HIPS   (.hipl or .hips)\n"
          "   MATLAB (.mat)\n"
          "   TIFF   (.tif or .tiff).\n"
          "\nNote that 8-bit output images, which are generated if the input\n"
          "image is 8-bit, are scaled to be in the range 0-255.\n", 
          BLUR_SIGMA) ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, 
          "  -blur <sigma>     - specify sigma of blurring "
          "kernel (default=%2.3f).\n", BLUR_SIGMA) ;
  fprintf(stderr, 
          "  -gaussian <sigma> - filter with Gaussian instead of median.\n") ;
  fprintf(stderr, 
          "  -mean             - filter with mean instead of median.\n") ;
  fprintf(stderr, 
          "  -w <window size>  - specify window used for offset calculation "
          "(default=%d).\n", OFFSET_WSIZE) ;
  fprintf(stderr, 
          "  --version         - display version #.\n") ;
  fprintf(stderr, 
          "  --help            - display this help message.\n\n") ;
#endif
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static void
read_parameter_file(char *fname)
{
  FILE *fp ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorExit(ERROR_NOFILE, "File %s not found.\n",fname);

  if (fscanf(fp,"%*s %f",&white_hilim) != 1)
    ErrorExit(ERROR_BADFILE, "could not scan white hilim from %s", fname) ;

  if (fscanf(fp,"%*s %f",&white_lolim) != 1)
    ErrorExit(ERROR_BADFILE, "could not scan white lolim from %s", fname) ;
  if (fscanf(fp,"%*s %f",&gray_hilim) != 1)
    ErrorExit(ERROR_BADFILE, "could not scan gray hilim from %s", fname) ;
  fclose(fp);
  printf("white_hilim = %f, white_lolim = %f, gray_hilim = %f\n",
         white_hilim,white_lolim,gray_hilim);
}

