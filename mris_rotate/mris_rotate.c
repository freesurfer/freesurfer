
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

static char vcid[] = "$Id: mris_rotate.c,v 1.1 1997/09/15 16:45:16 fischl Exp $";

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
  char         **av, *in_fname, *out_fname ;
  int          ac, nargs ;
  MRI_SURFACE  *mris_src, *mris_dst ;
  float        delta_theta, delta_phi ;

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

  if (argc < 5)
    usage_exit() ;

  in_fname = argv[1] ;
  if (sscanf(argv[2], "%f", &delta_phi) != 1)
    ErrorExit(ERROR_BADPARM, "%s: could not scan delta phi from %s",
              Progname, argv[2]) ;
  if (sscanf(argv[3], "%f", &delta_theta) != 1)
    ErrorExit(ERROR_BADPARM, "%s: could not scan delta phi from %s",
              Progname, argv[3]) ;
  out_fname = argv[4] ;

  mris_src = MRISread(in_fname) ;
  if (!mris_src)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;
  
  mris_dst = MRISrotate(mris_src,NULL,RADIANS(delta_phi),RADIANS(delta_theta));
  if (!mris_dst)
    ErrorExit(ERROR_NOFILE, "%s: could not rotate surface", Progname) ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing rotated surface to %s\n", out_fname) ;
  MRISwrite(mris_dst, out_fname) ;

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
          "usage: %s [options] <surface file> <delta phi> <delta theta>"
          " <output surface>\n", Progname) ;
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

