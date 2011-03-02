/**
 * @file  mri_build_priors.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:13 $
 *    $Revision: 1.8 $
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
#include "mrinorm.h"
#include "mriclass.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;

static int verbose = 1 ;
#define PRIOR_SCALE 3

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs ;
  char   *training_file_name, source_fname[100], target_fname[100], *cp,
  line[250], *output_file_name ;
  FILE   *fp ;
  int    fno, nfiles ;
  MRI    *mri_src, *mri_target, *mri_wm, *mri_priors = NULL ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_build_priors.c,v 1.8 2011/03/02 00:04:13 nicks Exp $", "$Name:  $");
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
    training_file_name = "train.dat" ;
  else
    training_file_name = argv[1] ;

  if (argc < 3)
    output_file_name = "priors.mnc" ;
  else
    output_file_name = argv[2] ;

  fp = fopen(training_file_name, "r") ;
  if (!fp)
    ErrorExit(ERROR_NO_FILE, "%s: could not open file %s",
              Progname, training_file_name) ;

  nfiles = 0 ;
  while ((cp = fgetl(line, 299, fp)) != NULL)
    nfiles++ ;
  fprintf(stderr, "processing %d files\n", nfiles) ;
  rewind(fp) ;

  fno = 0 ;
  while ((cp = fgetl(line, 299, fp)) != NULL) {
    sscanf(cp, "%s %s", source_fname, target_fname) ;
    fprintf(stderr, "file[%d]: %s --> %s\n", fno, source_fname, target_fname);
    mri_src = MRIread(source_fname) ;
    if (!mri_src) {
      fprintf(stderr, "could not read MR image %s\n", source_fname) ;
      continue ;
    }

    mri_wm = MRIread(target_fname) ;
    if (!mri_wm) {
      fprintf(stderr, "could not read MR image %s\n", target_fname) ;
      MRIfree(&mri_src) ;
      continue ;
    }

    mri_target = MRICbuildTargetImage(mri_src, NULL, mri_wm, 0, 0) ;
    MRIfree(&mri_src) ;
    MRIfree(&mri_wm) ;
    fno++ ;

    mri_priors = MRICupdatePriors(mri_target, mri_priors, PRIOR_SCALE) ;
    MRIfree(&mri_target) ;
  }
#if 0
  if (fno > 0)
    MRIscalarMul(mri_priors, mri_priors, 1.0f / (float)fno) ;
#else
  MRInormalizePriors(mri_priors) ;
#endif
  MRIwrite(mri_priors, output_file_name) ;
  fclose(fp) ;
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
  case 'V':
    verbose = !verbose ;
    break ;
  case '?':
  case 'U':
    printf("usage: %s <training file> <output file>", argv[0]) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
