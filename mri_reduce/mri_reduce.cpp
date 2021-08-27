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
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
int reductions = 1 ;

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs, i ;
  MRI    *mri_src, *mri_dst = NULL ;
  char   *in_fname, *out_fname ;

  nargs = handleVersionOption(argc, argv, "mri_reduce");
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

  if (argc < 1)
    argc = 1 ;

  if (argc < 1)
    ErrorExit(ERROR_BADPARM, "%s: no input name specified", Progname) ;
  in_fname = argv[1] ;

  if (argc < 2)
    ErrorExit(ERROR_BADPARM, "%s: no output name specified", Progname) ;
  out_fname = argv[2] ;

  fprintf(stderr, "reading from %s...", in_fname) ;
  mri_src = MRIread(in_fname) ;
  if (mri_src == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read src image from %s\n", Progname,in_fname) ;

  i = 0 ;
  do {
    if (i)
      mri_src = MRIcopy(mri_dst, NULL) ;
    fprintf(stderr, "\nreducing by 2\n");

    mri_dst = MRIallocSequence(MAX(1,mri_src->width/2), MAX(1,mri_src->height/2), MAX(1,mri_src->depth/2), MRI_FLOAT, mri_src->nframes);
    MRIreduce(mri_src, mri_dst) ;
    MRIfree(&mri_src) ;
  } while (++i < reductions) ;

  fprintf(stderr, "\nwriting to %s", out_fname) ;
  MRIwrite(mri_dst, out_fname) ;
  fprintf(stderr, "\n") ;
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
  case 'N':
    sscanf(argv[2], "%d", &reductions) ;
    fprintf(stderr, "reducing %d times\n", reductions) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    printf("usage: %s [input directory] [output directory]\n", argv[0]) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
