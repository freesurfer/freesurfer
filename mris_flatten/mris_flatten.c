
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

static char vcid[] = "$Id: mris_flatten.c,v 1.1 1997/10/23 19:45:22 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static INTEGRATION_PARMS  parms ;

int
main(int argc, char *argv[])
{
  char         **av, *in_surf_fname, *in_patch_fname, *out_patch_fname, 
               path[100] ;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.projection = PROJECT_PLANE ;
  parms.tol = 100*TOL ;
  parms.n_averages = N_AVERAGES ;
  parms.l_angle = L_ANGLE ;
  parms.l_area = L_AREA ;
  parms.niterations = NITERATIONS ;
  parms.write_iterations = 1 /*WRITE_ITERATIONS */;
  parms.a = parms.b = parms.c = 0.0f ;  /* ellipsoid parameters */

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  in_surf_fname = argv[1] ;
  in_patch_fname = argv[2] ;
  out_patch_fname = argv[3] ;
  FileNamePath(in_surf_fname, path) ;

  mris = MRISread(in_surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_surf_fname) ;


  if (MRISreadPatch(mris, in_patch_fname) != NO_ERROR)
    ErrorExit(ERROR_BADPARM, "%s: could not read patch file %s",
              Progname, in_patch_fname) ;
  MRISreadTriangleProperties(mris, in_surf_fname) ;
  MRISflattenPatch(mris, path) ;
  /*  MRISwritePatch(mris, out_patch_fname) ; exit(0) ;*/
  if (parms.niterations > 0)
  {
    MRISunfold(mris, &parms) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "writing flattened patch to %s\n", out_patch_fname) ;
    MRISwritePatch(mris, out_patch_fname) ;
  }

#if 0
  sprintf(fname, "%s.area_error", out_fname) ;
  printf("writing area errors to %s\n", fname) ;
  MRISwriteAreaError(mris, fname) ;
  sprintf(fname, "%s.angle_error", out_fname) ;
  printf("writing angle errors to %s\n", fname) ;
  MRISwriteAngleError(mris, fname) ;
  MRISfree(&mris) ;
#endif

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
  else if (!stricmp(option, "area"))
  {
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
  }
  else if (!stricmp(option, "name"))
  {
    strcpy(parms.base_name, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using base name = %s\n", parms.base_name) ;
  }

  else if (!stricmp(option, "angle"))
  {
    sscanf(argv[2], "%f", &parms.l_angle) ;
    nargs = 1 ;
    fprintf(stderr, "using l_angle = %2.3f\n", parms.l_angle) ;
  }
  else if (!stricmp(option, "area"))
  {
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
  }
  else if (!stricmp(option, "no_proj"))
  {
    parms.projection = NO_PROJECTION ;
    fprintf(stderr, "disabling ellipsoid projection\n") ;
  }
  else if (!stricmp(option, "tol"))
  {
    sscanf(argv[2], "%e", &parms.tol) ;
    nargs = 1 ;
    fprintf(stderr, "using tol = %2.2e\n", parms.tol) ;
  }
  else switch (toupper(*option))
  {
  case 'W':
    sscanf(argv[2], "%d", &parms.write_iterations) ;
    nargs = 1 ;
    fprintf(stderr, "using write iterations = %d\n", parms.write_iterations) ;
    break ;
  case 'A':
    sscanf(argv[2], "%d", &parms.n_averages) ;
    nargs = 1 ;
    fprintf(stderr, "using n_averages = %d\n", parms.n_averages) ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'N':
    sscanf(argv[2], "%d", &parms.niterations) ;
    nargs = 1 ;
    fprintf(stderr, "using niterations = %d\n", parms.niterations) ;
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
          "usage: %s [options] <surface file> <patch file name> <output patch>"
          "\n", Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program will add a template into an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

