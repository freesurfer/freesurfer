
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "tags.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"

static char vcid[] =
"$Id: mris_smooth.c,v 1.15 2006/11/15 19:57:16 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int normalize_flag = 0 ;
static char curvature_fname[STRLEN] = "curv" ;
static char area_fname[STRLEN] = "area" ;
static int nbrs = 2 ;
static int normalize_area = 0 ;
static int navgs = 10 ;
static int niterations = 10 ;
static int rescale = 0 ;
static int write_iterations = 0 ;
static double l_spring = 1.0 ;
static float momentum = 0.5 ;
static int no_write = 0 ;

// -g 20 8 works well for hippo
static double gaussian_norm = 0 ;
static int gaussian_avgs = 0 ;

int
main(int argc, char *argv[])
{
  char               **av, *in_fname, *out_fname, fname[STRLEN], path[STRLEN] ;
  int                ac, nargs ;
  MRI_SURFACE        *mris ;

  char cmdline[CMD_LINE_LEN] ;

  make_cmd_version_string
    (argc, argv,
     "$Id: mris_smooth.c,v 1.15 2006/11/15 19:57:16 fischl Exp $",
     "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
    (argc, argv,
     "$Id: mris_smooth.c,v 1.15 2006/11/15 19:57:16 fischl Exp $",
     "$Name:  $");
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

  if (argc < 3)
    print_help() ;

  in_fname = argv[1] ;
  out_fname = argv[2] ;
  FileNamePath(out_fname, path) ;

  mris = MRISfastRead(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;

  MRISaddCommandLine(mris, cmdline) ;
  MRISremoveTriangleLinks(mris) ;
  fprintf(stderr, "smoothing surface tessellation for %d iterations...\n",
          niterations);

  MRIScomputeMetricProperties(mris) ;
  MRISstoreMetricProperties(mris) ;
  MRISsetNeighborhoodSize(mris, nbrs) ;
  if (gaussian_norm > 0)
  {
    int i ;

    for (i = 0 ; i < niterations+5 ; i++)
    {
      if (i < niterations)
      {
        MRIScomputeSecondFundamentalForm(mris) ;
        MRISspringTermWithGaussianCurvature
          (mris, gaussian_norm, l_spring) ;
        MRISaverageGradients(mris, gaussian_avgs) ;
        MRISmomentumTimeStep(mris, momentum, 1, 1, 0) ;
        MRISclearGradient(mris) ;
        if ((write_iterations > 0) && (((i+1) % write_iterations) == 0))
        {
          char fname[STRLEN] ;
          sprintf(fname, "%s%04d", out_fname, i+1) ;
          printf("writing snapshot to %s...\n", fname) ;
          MRISwrite(mris, fname) ;
        }
      }
      else
      {
        MRISaverageVertexPositions(mris, 1) ;
        if (write_iterations > 0)
        {
          char fname[STRLEN] ;
          sprintf(fname, "%s%04d", out_fname, i+1) ;
          printf("writing snapshot to %s...\n", fname) ;
          MRISwrite(mris, fname) ;
        }
      }

    }
  }
  else
    MRISaverageVertexPositions(mris, niterations) ;

  fprintf(stderr, "smoothing complete - recomputing first and second "
          "fundamental forms...\n") ;
  MRIScomputeMetricProperties(mris) ;

  if (rescale)
    MRISscaleBrainArea(mris) ;
  MRIScomputeSecondFundamentalForm(mris) ;
  MRISuseMeanCurvature(mris) ;
  MRISaverageCurvatures(mris, navgs) ;
  if (normalize_flag)
    MRISnormalizeCurvature(mris) ;
  sprintf(fname, "%s.%s", mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh",
          curvature_fname);
  if (no_write == 0)
  {
    fprintf(stderr, "writing smoothed curvature to %s/%s\n", path,fname) ;
    MRISwriteCurvature(mris, fname) ;
    sprintf(fname, "%s.%s", mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh",
            area_fname);
    fprintf(stderr, "writing smoothed area to %s/%s\n", path, fname) ;
    MRISwriteArea(mris, fname) ;
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing smoothed surface to %s\n", out_fname) ;
  MRISwrite(mris, out_fname) ;
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
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size = %d\n", nbrs) ;
  }
  else if (!stricmp(option, "normalize"))
    normalize_flag = 1 ;
  else if (!stricmp(option, "nw"))
    no_write = 1 ;
  else if (!stricmp(option, "seed"))
  {
    setRandomSeed(atol(argv[2])) ;
    fprintf(stderr,"setting seed for random number generator to %d\n",
            atoi(argv[2])) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "area"))
  {
    normalize_area = 1 ;
    printf("normalizing area after smoothing\n") ;
  }
  else switch (toupper(*option))
  {
  case 'M':
    momentum = atof(argv[2]) ;
    printf("using momentum = %2.2f\n", momentum) ;
    nargs = 1 ;
    break ;
  case 'G':
    gaussian_norm = atof(argv[2]) ;
    gaussian_avgs = atoi(argv[3]) ;
    printf("using Gaussian curvature smoothing with norm %2.2f with %d smooth steps\n",
           gaussian_norm, gaussian_avgs) ;
    nargs = 2 ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    printf("debugging vertex %d\n", Gdiag_no) ;
    nargs = 1 ;
    break ;
  case 'W':
    write_iterations = atoi(argv[2]) ;
    printf("writing out snapshots every %d iterations\n", write_iterations) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'R':
    rescale = 1 ;
    fprintf(stderr, "rescaling brain area after smoothing...\n") ;
    break ;
  case 'C':
    strcpy(curvature_fname, argv[2]) ;
    nargs = 1 ;
    break ;
  case 'N':
    niterations = atoi(argv[2]) ;
    fprintf(stderr, "smoothing for %d iterations\n", niterations) ;
    nargs = 1 ;
    break ;
  case 'A':
    navgs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "averaging curvature for %d iterations\n", navgs) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
print_usage(void)
{
  fprintf(stderr,
          "usage: %s [options] <input surface> <output surface>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr,
          "\nThis program smooths the tessellation of a cortical surface and\n"
          "write out the first and second order properties after smoothing\n"
          "to the files $hemi.curv (mean curvature) and $hemi.area (area).\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-a <avgs>  "
          "specify # of curvature averaging iterations (def=10).\n") ;
  fprintf(stderr, "-n <niter> specify # of smoothing iterations (def=10).\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

