#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "timer.h"
#include "proto.h"
#include "mrinorm.h"


int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void  usage_exit(void) ;

char *Progname ;

static MRI_NORM_INFO  mni ;
static int verbose = 1 ;
static int num_3d_iter = 3 ;

static float intensity_above = 25 ;
static float intensity_below = 10 ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs, n ;
  MRI    *mri_src, *mri_dst = NULL ;
  char   *in_fname, *out_fname ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  mni.max_gradient = MAX_GRADIENT ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    usage_exit() ;

  if (argc < 1)
    ErrorExit(ERROR_BADPARM, "%s: no input name specified", Progname) ;
  in_fname = argv[1] ;

  if (argc < 2)
    ErrorExit(ERROR_BADPARM, "%s: no output name specified", Progname) ;
  out_fname = argv[2] ;

  if (verbose)
    fprintf(stderr, "reading from %s...\n", in_fname) ;
  mri_src = MRIread(in_fname) ;
  if (!mri_src)
    ErrorExit(ERROR_NO_FILE, "%s: could not open source file %s", 
              Progname, in_fname) ;

  if (verbose)
    fprintf(stderr, "normalizing image...\n") ;
  TimerStart(&start) ;
  MRInormInit(mri_src, &mni, 0, 0, 0, 0, 0.0f) ;
  mri_dst = MRInormalize(mri_src, NULL, &mni) ;
  if (!mri_dst)
    ErrorExit(ERROR_BADPARM, "%s: normalization failed", Progname) ;

  for (n = 0 ; n < num_3d_iter ; n++)
  {
    fprintf(stderr, "3d normalization pass %d of %d\n", n+1, num_3d_iter) ;
    MRI3dNormalize(mri_dst, NULL, DEFAULT_DESIRED_WHITE_MATTER_VALUE, mri_dst,
                   intensity_above, intensity_below);
  }

  if (verbose)
    fprintf(stderr, "writing output to %s\n", out_fname) ;
  MRIwrite(mri_dst, out_fname) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "3D bias adjustment took %d minutes and %d seconds.\n", 
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
  switch (toupper(*option))
  {
  case 'A':
    intensity_above = atof(argv[2]) ;
    fprintf(stderr, "using control point with intensity %2.1f above target.\n",
            intensity_above) ;
    nargs = 1 ;
    break ;
  case 'B':
    intensity_below = atof(argv[2]) ;
    fprintf(stderr, "using control point with intensity %2.1f below target.\n",
            intensity_below) ;
    nargs = 1 ;
    break ;
  case 'G':
    mni.max_gradient = atof(argv[2]) ;
    fprintf(stderr, "using max gradient = %2.3f\n", mni.max_gradient) ;
    nargs = 1 ;
    break ;
  case 'V':
    verbose = !verbose ;
    break ;
  case 'N':
    num_3d_iter = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "performing 3d normalization %d times\n", num_3d_iter) ;
    break ;
  case '?':
  case 'U':
    printf("usage: %s [input directory] [output directory]\n", argv[0]) ;
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
  fprintf(stderr, "usage: %s [options] <input volume> <output volume>\n", 
          Progname) ;
  fprintf(stderr, "\t-n <# of 3d normalization iterations>, default=5\n") ;
  fprintf(stderr, "\t-g <max intensity/mm gradient>, default=0.6\n") ;
  fprintf(stderr, "\t-s <sigma>, "
          "intensity interval pct for 3d, default=0.15\n") ;
  exit(0) ;
}

