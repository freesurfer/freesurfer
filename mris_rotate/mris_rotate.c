/**
 * @file  mris_rotate.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:33 $
 *    $Revision: 1.6 $
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
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"

static char vcid[] = "$Id: mris_rotate.c,v 1.6 2011/03/02 00:04:33 nicks Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

int
main(int argc, char *argv[]) {
  char         **av, *in_fname, *out_fname ;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;
  float        alpha, beta, gamma ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_rotate.c,v 1.6 2011/03/02 00:04:33 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)
    usage_exit() ;

  in_fname = argv[1] ;
  if (sscanf(argv[2], "%f", &alpha) != 1)
    ErrorExit(ERROR_BADPARM, "%s: could not scan alpha from %s",
              Progname, argv[2]) ;
  if (sscanf(argv[3], "%f", &beta) != 1)
    ErrorExit(ERROR_BADPARM, "%s: could not scan beta from %s",
              Progname, argv[3]) ;
  if (sscanf(argv[4], "%f", &gamma) != 1)
    ErrorExit(ERROR_BADPARM, "%s: could not scan gamma from %s",
              Progname, argv[4]) ;
  out_fname = argv[5] ;

  mris = MRISfastRead(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;

  alpha = RADIANS(alpha) ;
  beta = RADIANS(beta) ;
  gamma = RADIANS(gamma) ;
  MRIScenter(mris, mris) ;
  MRISrotate(mris, mris, alpha, beta, gamma) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not rotate surface", Progname) ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing rotated surface to %s\n", out_fname) ;
  MRISwrite(mris, out_fname) ;

  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option)) {
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
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <surface file> <delta phi> <delta theta>"
          " <output surface>\n", Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will add a template into an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

