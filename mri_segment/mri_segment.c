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
static int wsize = 5 ;
static float pct = 0.8 ;
static float pslope = 1.0f ;
static float nslope = 1.0f ;
static float wm_low = 90 ;
static float wm_hi = 130 ;
static float gray_hi = 100 ;
char *Progname ;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

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
  MRIhistoSegment(mri_src, mri_tmp, wm_low, wm_hi, gray_hi, 11, 3.0f) ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_tmp, "/tmp/int.mnc") ;
  fprintf(stderr, 
          "using local geometry to label remaining ambiguous voxels...\n") ;
  mri_labeled = MRIcpolvMedianCurveSegment(mri_src, mri_tmp, NULL, 5, 2.5);
  fprintf(stderr, 
          "\nreclassifying voxels using Gaussian border classifier...\n") ;

  /*
    now use the gray and white mattter border voxels to build a Gaussian
    classifier at each point in space and reclassify all voxels in the
    range [wm_low-5,gray_hi+5].
    */
  for (i = 0 ; i < 1 ; i++)
  {
    MRIreclassify(mri_src, mri_labeled, mri_labeled, wm_low-5, gray_hi+5, 7);
  }
  MRIfree(&mri_tmp) ;

  mri_dst = MRImaskLabels(mri_src, mri_labeled, NULL) ;

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
