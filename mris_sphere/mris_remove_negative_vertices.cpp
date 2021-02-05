/*
 * Original Author: Bruce Fischl
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
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"

#include "mri.h"
#include "mrisurf_project.h"

#include "error.h"
#include "tags.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static INTEGRATION_PARMS  parms ;

int
main(int argc, char *argv[])
{
  char         **av, *in_surf_fname, *out_fname,
  fname[STRLEN], *cp ;
  int          ac, nargs, msec ;
  MRI_SURFACE  *mris ;
  Timer then ;

  std::string cmdline = getAllInfo(argc, argv, "mris_remove_negative_vertices");

  nargs = handleVersionOption(argc, argv, "mris_remove_negative_vertices");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag = DIAG_SHOW ;

  then.reset() ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.dt = 1 ;
  parms.tol = .5 /*1e-1*/ ;
  parms.min_averages = 0 ;
  parms.l_angle = 0.0 /* L_ANGLE */ ;
  parms.l_area = 0.0 /* L_AREA */ ;
  parms.l_neg = 0.0 ;
  parms.l_spring = 0.0 ;
  parms.l_boundary = 0.0 ;
  parms.l_curv = 0.0 ;
  parms.niterations = 1000 ;
  parms.a = parms.b = parms.c = 0.0f ;  /* ellipsoid parameters */
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = 0.0 ;
  parms.base_name[0] = 0 ;

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

  in_surf_fname = argv[1] ;
  out_fname = argv[2] ;

  if (parms.base_name[0] == 0)
  {
    FileNameOnly(out_fname, fname) ;
    cp = strchr(fname, '.') ;
    if (cp)
      strcpy(parms.base_name, cp+1) ;
    else
      strcpy(parms.base_name, "sphere") ;
  }

  mris = MRISread(in_surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_surf_fname) ;

  MRISaddCommandLine(mris, cmdline) ;

  MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
  MRISremoveOverlapWithSmoothing(mris,&parms) ;


  printf("writing output to %s\n", out_fname) ;
  MRISwrite(mris, out_fname) ;
  msec = then.milliseconds() ;
  fprintf(stderr, 
          "regularization of spherical transformation took %2.2f hours\n",
          (float)msec/(1000.0f*60.0f*60.0f));
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
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'M':
      parms.integration_type = INTEGRATE_MOMENTUM ;
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "momentum = %2.2f\n", (float)parms.momentum) ;
      break ;
    case 'W':
      Gdiag |= DIAG_WRITE ;
      sscanf(argv[2], "%d", &parms.write_iterations) ;
      nargs = 1 ;
      fprintf(stderr, "using write iterations = %d\n", 
              parms.write_iterations) ;
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
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

