/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static char *log_fname = NULL ;
static void usage_exit(int code) ;

static int quiet = 0 ;
static int all_flag = 0 ;

int
main(int argc, char *argv[]) {
  char   **av ;
  int    label, ac, nargs ;
  int          msec, minutes, seconds/*, wrong, total, correct*/ ;
  Timer start ;
  MRI    *mri_T1, *mri_labeled ;
  FILE   *log_fp ;
  HISTOGRAM *histo ;

  nargs = handleVersionOption(argc, argv, "mri_label_histo");
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

  if (argc < 4)
    usage_exit(1) ;


  mri_T1 = MRIread(argv[1]) ;
  if (!mri_T1)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s",Progname,
              argv[1]) ;
  mri_labeled = MRIread(argv[2]) ;
  if (!mri_labeled)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s",Progname,
              argv[2]) ;

  if (log_fname) {
    log_fp = fopen(log_fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, log_fname) ;
  } else
    log_fp = NULL ;


  label = atoi(argv[3]) ;
  printf("generating histogram for label %d...\n", label) ;
  histo = MRIhistogramLabel(mri_T1, mri_labeled, label, 0) ;

  HISTOplot(histo, argv[4]) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;

  if (DIAG_VERBOSE_ON)
    fprintf(stderr, "overlap calculation took %d minutes and %d seconds.\n",
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
  case 'Q':
    quiet = 1 ;
    break ;
  case 'A':
    all_flag = 1 ;
    break ;
  case 'L':
    log_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "logging results to %s\n", log_fname) ;
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
  printf("usage: %s [options] <T1 volume> <labeled volume> <label> <output>\n",
         Progname) ;
  exit(code) ;
}
