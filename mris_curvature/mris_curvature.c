
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

static char vcid[] = "$Id: mris_curvature.c,v 1.3 1997/11/03 17:27:09 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int write_flag = 0 ;

int
main(int argc, char *argv[])
{
  char         **av, *in_fname, fname[100], hemi[10], path[100], name[100],*cp;
  int          ac, nargs, nhandles ;
  MRI_SURFACE  *mris ;
  double       ici, fi ;

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

  if (argc < 2)
    usage_exit() ;

  in_fname = argv[1] ;

  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;

#if 0
MRISprojectOntoCylinder(mris, 50.0f) ;
MRISwrite(mris,in_fname) ;
#endif

#if 0
MRISprojectOntoEllipsoid(mris, mris, 50.0f, 50.0f, 50.0f) ;
MRISwrite(mris, in_fname) ;
mrisComputeNormals(mris) ;
#endif

  MRIScomputeSecondFundamentalForm(mris) ;
  nhandles = nint(1.0 - mris->Ktotal / (4.0*M_PI)) ;
  fprintf(stderr, "total integrated curvature = %2.3f*4pi (%2.3f) --> "
          "%d handles\n", (float)mris->Ktotal/(4.0f*M_PI), 
          (float)mris->Ktotal, nhandles) ;
  
  MRIScomputeCurvatureIndices(mris, &ici, &fi);
  fprintf(stderr, "ICI = %2.1f, FI = %2.1f\n", ici, fi) ;

  if (write_flag)
  {
    FileNamePath(in_fname, path) ;
    FileNameOnly(in_fname, name) ;
    cp = strchr(name, '.') ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: could not scan hemisphere from '%s'",
                Progname, fname) ;
    strncpy(hemi, cp-2, 2) ;
    hemi[2] = 0 ;
    
    MRISuseGaussianCurvature(mris) ;
    sprintf(fname, "%s/%s.K", path,name) ; MRISwriteCurvature(mris, fname) ;
    MRISuseMeanCurvature(mris) ;
    sprintf(fname, "%s/%s.H", path,name) ; MRISwriteCurvature(mris, fname) ;
  }
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
    write_flag = 1 ;
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
  fprintf(stderr, "usage: %s [options] <input surface file>\n", Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
      "\nThis program will the second fundamental form of a cortical surface."
          "\nIt will create two new files <hemi>.H and <hemi>.K with the\n"
          "mean and Gaussian curvature respectively.");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

