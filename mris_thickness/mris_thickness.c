
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

static char vcid[] = "$Id: mris_thickness.c,v 1.2 1999/01/21 18:29:58 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

int
main(int argc, char *argv[])
{
  char          **av, *out_fname, *sname, sdir[400], *cp, fname[500], *hemi ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;

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

  sname = argv[1] ;
  hemi = argv[2] ;
  out_fname = argv[3] ;
  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, 
              "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
  strcpy(sdir, cp) ;
  

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, GRAY_MATTER_NAME) ;
  fprintf(stderr, "reading gray matter surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  if (MRISreadOriginalProperties(mris, WHITE_MATTER_NAME) != NO_ERROR)
    ErrorExit(Gerror, "%s: could not read white matter surface", Progname) ;
  fprintf(stderr, "measuring gray matter thickness...\n") ;

  MRISmeasureCorticalThickness(mris) ;

#if 0
  sprintf(fname, "%s/%s/surf/%s", sdir, sname, out_fname) ;
  fprintf(stderr, "writing output surface to %s...\n", fname) ;
#endif
  fprintf(stderr, "writing thickness to curvature file %s...\n", out_fname) ;
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
  else switch (toupper(*option))
  {
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
          "usage: %s [options] <subject name> <hemi> <thickness file>\n", 
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program measures the thickness of the cortical surface\n"
          "and writes the resulting measurement into a 'curvature' file\n"
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

