
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

static char vcid[] = "$Id: mris_curvature.c,v 1.8 1998/01/27 00:46:23 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int write_flag = 0 ;
static int nbrs = 2 ;
static int navgs = 0 ;
static char *param_file = NULL ;
static int normalize = 0 ;

int
main(int argc, char *argv[])
{
  char         **av, *in_fname, fname[100], hemi[10], path[100], name[100],*cp;
  int          ac, nargs, nhandles ;
  MRI_SURFACE  *mris ;
  double       ici, fi, var ;

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

  MRISsetNeighborhoodSize(mris, nbrs) ;

  if (param_file)
  {
    MRI_SP *mrisp ;
    mrisp = MRISPread(param_file) ;
    MRISfromParameterization(mrisp, mris, 0) ;
    MRISPfree(&mrisp) ;
    if (normalize)
      MRISnormalizeCurvature(mris) ;
    FileNamePath(in_fname, path) ;
    FileNameOnly(in_fname, name) ;
    cp = strchr(name, '.') ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: could not scan hemisphere from '%s'",
                Progname, fname) ;
    strncpy(hemi, cp-2, 2) ;
    hemi[2] = 0 ;
    sprintf(fname, "%s/%s.param", path,name) ; 
    fprintf(stderr, "writing parameterized curvature to %s...", fname) ;
    MRISwriteCurvature(mris, fname) ;
    fprintf(stderr, "done.\n") ;
  }
  else
  {
    MRIScomputeSecondFundamentalForm(mris) ;
    nhandles = nint(1.0 - mris->Ktotal / (4.0*M_PI)) ;
    fprintf(stderr, "total integrated curvature = %2.3f*4pi (%2.3f) --> "
            "%d handles\n", (float)mris->Ktotal/(4.0f*M_PI), 
            (float)mris->Ktotal, nhandles) ;
    
    fprintf(stderr, "0: k1 = %2.3f, k2 = %2.3f, H = %2.3f, K = %2.3f\n",
            mris->vertices[0].k1, mris->vertices[0].k2, 
            mris->vertices[0].H, mris->vertices[0].K) ;
    fprintf(stderr, "0: vnum = %d, v2num = %d, total=%d, area=%2.3f\n",
            mris->vertices[0].vnum, mris->vertices[0].v2num,
            mris->vertices[0].vtotal,mris->vertices[0].area) ;
    MRIScomputeCurvatureIndices(mris, &ici, &fi);
    var = MRIStotalVariation(mris) ;
    fprintf(stderr,"ICI = %2.1f, FI = %2.1f, variation=%2.3f\n", ici, fi, var);

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
      MRISaverageCurvatures(mris, navgs) ;
      sprintf(fname, "%s/%s.K", path,name) ; 
      fprintf(stderr, "writing Gaussian curvature to %s...", fname) ;
      MRISwriteCurvature(mris, fname) ;
      MRISuseMeanCurvature(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      sprintf(fname, "%s/%s.H", path,name) ; 
      fprintf(stderr, "done.\nwriting mean curvature to %s...", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }
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
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
  }
  else switch (toupper(*option))
  {
  case 'N':
    normalize = 1 ;
    break ;
  case 'P':
    param_file = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "using parameterization file %s\n", param_file) ;
    break ;
  case 'A':
    navgs = atoi(argv[2]) ;
    fprintf(stderr, "averaging curvature patterns %d times.\n", navgs) ;
    nargs = 1 ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
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

