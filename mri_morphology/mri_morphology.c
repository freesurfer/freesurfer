/**
 * @file  mri_morphology.c
 * @brief applies morphological operations to a volume.
 *
 * Program for applying morphological operations open, close, dilate, and erode either to an
 * entire volume or to only a label (specified with -l <label>) within it.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/01/05 17:24:11 $
 *    $Revision: 1.7 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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

static int label = -1 ;


int
main(int argc, char *argv[]) {
  char   *out_fname, **av ;
  int    ac, nargs, niter, operation, i ;
  MRI    *mri_src, *mri_dst, *mri_saved_src = NULL ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_morphology.c,v 1.7 2007/01/05 17:24:11 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
    usage_exit(1) ;

  out_fname = argv[argc-1] ;
  mri_src = MRIread(argv[1]) ;
  if (mri_src == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read input volume %s\n",
              Progname,argv[1]);

  if (!stricmp(argv[2], "dilate"))
    operation = DILATE ;
  else  if (!stricmp(argv[2], "open"))
    operation = OPEN ;
  else  if (!stricmp(argv[2], "close"))
    operation = CLOSE ;
  else  if (!stricmp(argv[2], "erode"))
    operation = ERODE ;
  else  if (!stricmp(argv[2], "open"))
    operation = OPEN ;
  else  if (!stricmp(argv[2], "mode"))
    operation = OPEN ;
  else {
    operation = 0 ;
    ErrorExit(ERROR_UNSUPPORTED, "morphological operation '%s'  is not supported", argv[2]) ;
  }

  if (label > 0)  // erase everything but label in src, and save orig vol
  {
    MRI *mri_tmp ;

    mri_saved_src = MRIcopy(mri_src, NULL) ;
    mri_tmp = MRIclone(mri_src, NULL);
    MRIcopyLabel(mri_src, mri_tmp,label) ;
    MRIfree(&mri_src) ;
    mri_src = mri_tmp ;
  }
    
  niter = atoi(argv[3]) ;
  switch (operation) {
  case MODE_FILTER: {
    MRI *mri_tmp ;
    if (mri_src->type != MRI_UCHAR) {
      mri_tmp = MRIchangeType(mri_src, MRI_UCHAR, 0, 1, 1) ;
      MRIfree(&mri_src) ;
      mri_src = mri_tmp ;
    }
    printf("applying mode filter %d times\n", niter) ;
    mri_dst = MRImodeFilter(mri_src, NULL, niter) ;
    break ;
  }
  case DILATE:
    mri_dst = NULL ;
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIdilate(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    break ;
  case CLOSE:
    mri_dst = NULL ;
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIdilate(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIerode(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    break ;
  case OPEN:
    mri_dst = NULL ;
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIerode(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIdilate(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    break ;
  case ERODE:
    mri_dst = NULL ;
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIerode(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    break ;
  default:
    break ;
  }
  if (label > 0)
  {
    MRIreplaceValues(mri_saved_src, mri_saved_src, label, 0) ;
    MRIcopyLabel(mri_dst, mri_saved_src, label) ;
    MRIcopy(mri_saved_src, mri_dst) ;
    MRIfree(&mri_saved_src) ;
  }

  fprintf(stderr, "writing to %s...\n", out_fname) ;
  MRIwrite(mri_dst, out_fname) ;
  MRIfree(&mri_dst) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "morphological processing took %d minutes and %d seconds.\n",
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
  if (!stricmp(option, "dt")) {}
  else switch (toupper(*option)) {
  case 'L':
    label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("only applying operations to label %d\n", label) ;
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
  printf("usage: %s [options] <volume> <operation> <# iter> <out volume>\n", Progname) ;
  printf("\twhere <operation> can be [open,close,dilate,erode,mode]\n");
  printf("\tvalid options are:\n") ;
  printf("\t-l <label>  only apply operations to <label> instead of all nonzero voxels\n") ;
  exit(code) ;
}

