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

static char vcid[] = "$Id: mris_add_template.c,v 1.1 1997/09/12 21:49:33 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int nsurfaces=0 ;   /* total # of surfaces currently in average brain */

int
main(int argc, char *argv[])
{
  char         **av, *surf_fname, *template_fname ;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;
  MRI          *mri ;
  MRI_SP       *mrisp, *mrisp_blur ;
  float        sigma ;

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

  surf_fname = argv[1] ;
  template_fname = argv[2] ;
  if (argc > 3)
    sigma = atof(argv[3]) ;
  else
    sigma = 5.0f ;

  mris = MRISread(surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, surf_fname) ;

  mrisp = MRIStoParameterization(mris, NULL, 1) ;
  MRISfromParameterization(mrisp, mris) ;
  MRISwriteCurvature(mris, "surf/lh.curv") ;
#if 1
  mrisp_blur = MRISPblur(mrisp, NULL, sigma) ;
  MRISfromParameterization(mrisp_blur, mris) ;
  MRISwriteCurvature(mris, "surf/lh.curv_blur") ;
#endif
  exit(0) ;
  fprintf(stderr, "current average consists of %d surfaces\n", nsurfaces) ;

  if (!nsurfaces)  /* first time - create average surface */
  {
    mri = MRIalloc(256, 256, 256, MRI_FLOAT) ;
    if (!mri)
      ErrorExit(ERROR_NOFILE, "%s: could not allocate template brain") ;
  }
  else
  {
    mri = MRIread(template_fname) ;
    if (!mri)
      ErrorExit(ERROR_NOFILE, "%s: could not read template brain %s",
                template_fname) ;
  }

  fprintf(stderr, "averaging surface %s into %d averages %s\n",
          surf_fname, nsurfaces, template_fname) ;
  MRIwrite(mri, template_fname) ;

  MRISPfree(&mrisp) ;
  MRIfree(&mri) ;
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
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option))
  {
  case 'N':
    sscanf(argv[2], "%d", &nsurfaces) ;
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
          "usage: %s [options] <surface file> <average surface>\n",
          Progname) ;
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

