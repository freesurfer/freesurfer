
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
#include "version.h"

static char vcid[] = "$Id: mris_density.c,v 1.1 2006/12/18 18:43:36 fischl Exp $";


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static double resolution = 1.0/8.0 ;
static double radius = 20 ;

int
main(int argc, char *argv[])
{
  char          **av, *out_fname, *in_fname ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_density.c,v 1.1 2006/12/18 18:43:36 fischl Exp $", "$Name:  $");
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
    usage_exit() ;

  in_fname = argv[1] ;
  out_fname = argv[2] ;

  fprintf(stderr, "reading surface from %s...\n", in_fname) ;
  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",Progname, in_fname) ;
  MRISmakeDensityMap(mris, resolution, radius) ;

  fprintf(stderr, "writing density map to curvature file %s...\n", out_fname) ;
  MRISwriteCurvature(mris, out_fname) ;
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
  else if (!stricmp(option, "radius"))
  {
    radius = atof(argv[2]) ;
    nargs =  1 ;
    printf("using radius = %2.3f\n", radius) ;
  }
  else switch (toupper(*option))
  {
  case 'R':
    resolution = (double)atof(argv[2]) ;
    nargs = 1 ;
    printf("setting resolution for intermediate calculations to %2.4f\n", resolution) ;
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
          "usage: %s [options] <input surface> <output map>\n", 
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program computes a density map and attaches it to each vertex on the surface.\n"
          "The resulting measurement are written into a 'curvature' file\n"
          "<output file>.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  /*  fprintf(stderr, "-n    normalize output curvatures.\n") ;*/
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


