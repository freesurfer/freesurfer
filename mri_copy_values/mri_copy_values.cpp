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

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs ;
  MRI    *mri_src, *mri_dst ;
  char   *in_fname, *out_fname ;
  int    label, nvox ;

  nargs = handleVersionOption(argc, argv, "mri_copy_values");
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

  if (argc < 2)
    ErrorExit(ERROR_BADPARM, "%s: no input name specified", Progname) ;
  in_fname = argv[1] ;

  if (argc < 3)
    ErrorExit(ERROR_BADPARM, "%s: no value specified", Progname) ;
  label = atoi(argv[2]) ;

  if (argc < 4)
    ErrorExit(ERROR_BADPARM, "%s: no output name specified", Progname) ;
  out_fname = argv[3] ;

  fprintf(stderr, "reading from %s...\n", in_fname) ;
  mri_src = MRIread(in_fname) ;
  if (!mri_src)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, in_fname) ;
  mri_dst = MRIread(out_fname) ;
  if (!mri_dst)
    ErrorExit(ERROR_NOFILE, "%s: could not read destination volume %s",
              Progname, out_fname) ;
  nvox = MRIcopyLabel(mri_src, mri_dst, label) ;
  fprintf(stderr, "%d voxels copied from input to output volume...\n", nvox);

  fprintf(stderr, "writing to %s...\n", out_fname) ;
  MRIwrite(mri_dst, out_fname) ;
  MRIfree(&mri_dst) ;
  MRIfree(&mri_src) ;
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
  case '?':
  case 'U':
    printf("usage: %s <input volume> <label> <output volume>\n", argv[0]) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
