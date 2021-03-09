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
#include "utils.h"
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;


static int erode = 0 ;
static double cfmin = 10 ;
static double cfmax = 5000 ;
static double step = 10 ;


int
main(int argc, char *argv[]) {
  char         **av ;
  int          ac, nargs, i, j, nbins ;
  int          msec, minutes, seconds, x, y, z, **joint_density ;
  Timer start ;
  MRI          *mri1, *mri2, *mri_nonbrain, *mri_tmp ;
  FILE         *fp ;
  float        fmin1, fmax1, fmin2, fmax2 ;
  double       val1, val2 ;

  nargs = handleVersionOption(argc, argv, "mri_joint_density");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  if (argc < 3)
    usage_exit(1) ;

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

  printf("reading volume 1 from %s...\n", argv[1]) ;
  mri1 = MRIread(argv[1]) ;
  if (!mri1)
    ErrorExit(ERROR_NOFILE, "%s: could not read MR volume 1 from %s", Progname, argv[1]) ;
  printf("reading volume 2 from %s...\n", argv[2]) ;
  mri2 = MRIread(argv[2]) ;
  if (!mri2)
    ErrorExit(ERROR_NOFILE, "%s: could not read MR volume 2 from %s", Progname, argv[2]) ;

  MRIvalRange(mri1, &fmin1, &fmax1) ;
  MRIvalRange(mri2, &fmin2, &fmax2) ;
  printf("MR ranges [%2.1f %2.1f], [%2.1f %2.1f]\n", fmin1, fmax1, fmin2, fmax2) ;

  mri_nonbrain = MRIand(mri1, mri2, NULL, 1) ;
  MRIbinarize(mri_nonbrain, mri_nonbrain, 1, 128, 0) ;
  mri_tmp = MRIallocSequence(mri1->width, mri1->height, mri1->depth, MRI_UCHAR, mri1->nframes) ;
  MRIcopy(mri_nonbrain, mri_tmp) ;
  MRIfree(&mri_nonbrain) ;
  mri_nonbrain = mri_tmp ;
  MRIdilateLabel(mri_nonbrain, mri_nonbrain, 128, erode) ;

  /* cfmin = MIN(fmin1, fmin2) ; cfmax = MAX(fmax1, fmax2) ;*/

  nbins = nint(((cfmax - cfmin + 1) / step)+.99) ;
  printf("computing joint density with %d bins...\n", nbins) ;

  joint_density = (int **)calloc(nbins, sizeof(int *)) ;
  if (!joint_density)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d bin joint density\n", nbins) ;
  for (i = 0 ; i < nbins ; i++) {
    joint_density[i] = (int *)calloc(nbins, sizeof(int)) ;
    if (!joint_density[i])
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d bin joint density\n", nbins) ;
  }

  for (x = 0 ; x < mri1->width ; x++) {
    for (y = 0 ; y < mri1->height ; y++) {
      for (z = 0 ; z < mri1->depth ; z++) {
        MRIsampleVolume(mri_nonbrain, x, y, z, &val1) ;
        if (FZERO(val1) != 0)
          continue ;  /* close to border of brain */

        MRIsampleVolume(mri1, x, y, z, &val1) ;
        if (val1 > cfmax || val1 < cfmin)
          continue ;
        MRIsampleVolume(mri2, x, y, z, &val2) ;
        if (val2 > cfmax || val2 < cfmin)
          continue ;
        i = nint((val1-cfmin)/step) ;
        j = nint((val2-cfmin)/step) ;
        joint_density[i][j]++ ;
      }
    }
  }

  printf("writing joint density to %s...\n", argv[3]) ;
  fp = fopen(argv[3], "w") ;
  if (fp == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open output file %s", Progname, argv[3]) ;
  fprintf(fp, "nbins = %d;\n cfmin = %f ;\n cfmax = %f  ; \nfstep = %f;\n", nbins, cfmin, cfmax, step) ;
  fprintf(fp, "joint_density = [...\n") ;
  for (i = 0  ; i < nbins ; i++) {
    for (j = 0 ; j < nbins ; j++) {
      fprintf(fp, "%d ", joint_density[i][j]) ;
    }
    fprintf(fp, ";\n") ;
  }
  fprintf(fp, "];\n") ;
  fclose(fp) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("joint density estimation took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  if (!stricmp(option, "step")) {
    step = atof(argv[2]) ;
    nargs = 1 ;
    printf("using step size %2.2f in histogram\n", step) ;
  } else if (!stricmp(option, "max")) {
    cfmax = atof(argv[2]) ;
    nargs = 1 ;
    printf("using max val %2.2f in histogram\n", cfmax) ;
  } else if (!stricmp(option, "erode")) {
    erode = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d erosions to prevent alignment of brain with non-brain\n", erode) ;
  } else if (!stricmp(option, "min")) {
    cfmin = atof(argv[2]) ;
    nargs = 1 ;
    printf("using min val %2.2f in histogram\n", cfmin) ;
  } else switch (toupper(*option)) {
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
  printf("usage: %s [options] <vol1> <vol2> <output density file>", Progname) ;
  exit(code) ;
}
