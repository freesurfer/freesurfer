/**
 * @file  mri_label_accuracy.c
 * @brief computes the accuracy of the labeling of two images -  one input and one reference.
 *
 * Program computes the accuracy of the labeling of two images -  one input and one reference.
 * A set of algorithms can be selected from, the default is the mean of the symmetric min distances
 * between the boundaries
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:22 $
 *    $Revision: 1.2 $
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
#include "timer.h"
#include "version.h"
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;


char *Progname ;

#define ERODE          1
#define DILATE         2
#define CLOSE          3
#define OPEN           4
#define DILATE_LABEL   5
#define MODE_FILTER    6


static void usage_exit(int code) ;

static int target_label = -1 ;

#define PAD 20

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs ;
  MRI    *mri_src, *mri_ref, *mri_tmp ;
  double accuracy ;
  MRI_REGION box ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_label_accuracy.c,v 1.2 2011/03/02 00:04:22 nicks Exp $", "$Name:  $");
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

  if (argc < 3)
    usage_exit(1) ;

  mri_src = MRIread(argv[1]) ;
  if (mri_src == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read input volume %s\n", Progname,argv[1]);
  MRIboundingBox(mri_src, 0, &box) ;
  mri_tmp = MRIextractRegionAndPad(mri_src, NULL, &box, PAD) ;
  MRIfree(&mri_src) ; mri_src = mri_tmp ;
  if (mri_src->type == MRI_SHORT)
  {
    mri_tmp = MRIchangeType(mri_src, MRI_FLOAT, 0, 0, 0) ;
    MRIfree(&mri_src) ; mri_src = mri_tmp ;
  }

  mri_ref = MRIread(argv[2]) ;
  if (mri_ref == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read reference volume %s\n", Progname,argv[1]);
  MRIboundingBox(mri_ref, 0, &box) ;
  mri_tmp = MRIextractRegionAndPad(mri_ref, NULL, &box, PAD) ;
  MRIfree(&mri_ref) ; mri_ref = mri_tmp ;

  accuracy = MRIcomputeLabelAccuracy(mri_src, mri_ref, MRI_MEAN_MIN_DISTANCE, stdout) ;

  if (Gdiag_fp)
    fclose(Gdiag_fp) ;
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
  if (!stricmp(option, "log")) {
    Gdiag_fp = fopen(argv[2], "w") ;
    if (Gdiag_fp == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not open log file %s", argv[2]) ;
    nargs = 1 ;
  }
  else switch (toupper(*option)) {
  case 'L':
    Gdiag_no = target_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("only applying operations to label %d\n", target_label) ;
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
  printf("usage: %s [options] <src volume> <ref volume>\n", Progname) ;
  printf("\tvalid options are:\n") ;
  printf("\t-l <label>  only apply operations to <label> instead of all nonzero voxels\n") ;
  exit(code) ;
}

