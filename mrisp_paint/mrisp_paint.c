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

static char vcid[] = "$Id: mrisp_paint.c,v 1.1 1998/10/27 01:00:26 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int normalize = 0 ;


int
main(int argc, char *argv[])
{
  char         **av, *surf_fname, *template_fname, *out_fname, *cp;
  int          ac, nargs, param_no = 0 ;
  MRI_SURFACE  *mris ;
  MRI_SP       *mrisp ;

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

  if (argc < 4)
    usage_exit() ;

  template_fname = argv[1] ;
  surf_fname = argv[2] ;
  out_fname = argv[3] ;

  fprintf(stderr, "reading surface from %s...\n", surf_fname) ;
  mris = MRISread(surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, surf_fname) ;

  if (normalize)
  {
    cp = strchr(template_fname, '#') ;
    if (cp)   /* # explicitly given */
    {
      param_no = atoi(cp+1) ;
      *cp = 0 ;
    }
    else
      param_no = 0 ;
  }
  else
  {
    cp = strchr(template_fname, '#') ;
    if (cp)   /* # explicitly given */
    {
      param_no = atoi(cp+1) ;
      *cp = 0 ;
    }
    else
      param_no = 0 ;
  }
  fprintf(stderr, "reading template parameterization from %s...\n",
          template_fname) ;
  mrisp = MRISPread(template_fname) ;
  if (!mrisp)
    ErrorExit(ERROR_NOFILE, "%s: could not open template file %s",
                Progname, template_fname) ;

  if (normalize)
    MRISnormalizeFromParameterization(mrisp, mris, param_no) ;
  else
    MRISfromParameterization(mrisp, mris, param_no) ;
  fprintf(stderr, "writing curvature file to %s...\n", out_fname) ;
  MRISwriteCurvature(mris, out_fname) ;

  MRISPfree(&mrisp) ;
  MRISfree(&mris) ;
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
  int    nargs = 0 ;
  char   *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option))
  {
  case 'N':
    normalize = 1 ;
    fprintf(stderr, "normalizing curvature by variance.\n") ;
    break ;
  case 'W':
    Gdiag |= DIAG_WRITE ;
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
       "usage: %s [options] <parameterization file> <input surface> <output name>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program paint a parameterization onto a surface and output.\n"
          "the results as a curvature file.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

