
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"

static char vcid[] = "$Id: mris_smooth.c,v 1.1 1998/03/26 05:19:26 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int normalize_flag = 0 ;
static char curvature_fname[100] = "curv" ;
static char area_fname[100] = "area" ;
static int nbrs = 2 ;
static int navgs = 10 ;
static int niterations = 10 ;

int
main(int argc, char *argv[])
{
  char               **av, *in_fname, *out_fname, fname[100], path[100] ;
  int                ac, nargs ;
  MRI_SURFACE        *mris ;

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

  in_fname = argv[1] ;
  out_fname = argv[2] ;
  FileNamePath(out_fname, path) ;

  mris = MRISfastRead(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;

  fprintf(stderr, "smoothing surface tessellation for %d iterations...\n",
          niterations);

  MRIScomputeMetricProperties(mris) ;
  MRISstoreMetricProperties(mris) ;
  MRISsetNeighborhoodSize(mris, nbrs) ;
  MRISaverageVertexPositions(mris, niterations) ;

  fprintf(stderr, "smoothing complete - recomputing first and second "
          "fundamental forms...\n") ;
  MRIScomputeMetricProperties(mris) ;
  MRISscaleBrainArea(mris) ;
  MRIScomputeSecondFundamentalForm(mris) ;
  MRISuseMeanCurvature(mris) ;
  MRISaverageCurvatures(mris, navgs) ;
  if (normalize_flag)
    MRISnormalizeCurvature(mris) ;
  sprintf(fname, "%s.%s", mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh",
          curvature_fname);
  fprintf(stderr, "writing smoothed curvature to %s/%s\n", path,fname) ;
  MRISwriteCurvature(mris, fname) ;
  sprintf(fname, "%s.%s", mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh",
          area_fname);
  fprintf(stderr, "writing smoothed area to %s/%s\n", path, fname) ;
  MRISwriteArea(mris, fname) ;
  
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
  else switch (toupper(*option))
  {
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
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
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <input surface> <sigma> <output curvature file name>\n", 
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program smooth the tessellation of a cortical surface and\n"
          "write out the first and second order properties after smoothing\n"
          "to the files $hemi.curv (mean curvature) and $hemi.area (area).\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-n    normalize output curvatures.\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

