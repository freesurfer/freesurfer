/**
 * @file  mri_mosaic.c
 * @brief program for taking a set of images and constructing one large image 
 * that contains them all
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2015/09/25 16:05:20 $
 *    $Revision: 1.1 $
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

static int rectify = 0 ;

#define MAX_INPUT_IMAGES 10000

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs, i ;
  char   *in_fname, *out_fname ;
  int          msec, minutes, seconds, nimages ;
  Timer start ;
  MRI          *mri[MAX_INPUT_IMAGES], *mri_mosaic ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_mosaic.c,v 1.1 2015/09/25 16:05:20 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;

  out_fname = argv[argc-1] ;
  nimages = argc-2 ;
  printf("processing %d images and writing output to %s\n", nimages,out_fname) ;

  for (i = 0 ; i < nimages ; i++)
  {
    in_fname = argv[i+1] ;
    printf("reading image %d: %s\n", i, in_fname) ;
    mri[i] = MRIread(in_fname) ;
    if (mri[i] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read %dth input volume %s",
		Progname, i, in_fname) ;
  }
  mri_mosaic = MRImakeMosaic(mri, nimages, rectify) ;

  MRIwrite(mri_mosaic, out_fname) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "mosaicing took %d minutes and %d seconds.\n", 
	  minutes, seconds) ;
  exit(0) ;
  return(0) ;
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
  switch (toupper(*option)) {
  case 'R':
    rectify = 1 ;
    printf("rectifying input\n") ;
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s [options] <inverse operator> <EEG/MEG data file>",
         Progname) ;
  printf(
    "\tf <f low> <f hi> - apply specified filter (not implemented yet)\n"
  );
  printf("\tn - noise-sensitivity normalize inverse (default=1)") ;
  exit(code) ;
}

