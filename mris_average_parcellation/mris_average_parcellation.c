/**
 * @file  mris_average_parcellation.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:27 $
 *    $Revision: 1.4 $
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
#include "mrisurf.h"
#include "utils.h"
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;


char *Progname ;

static char sdir[STRLEN] = "" ;
static void usage_exit(int code) ;

int
main(int argc, char *argv[]) {
  char        **av, *cp, *annot_name, fname[STRLEN], *hemi ;
  int         **counts, ac, nargs, i, j, num, nvertices, vno ;
  MRI_SURFACE *mris ;
  char        *subject_name, *out_fname ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  VERTEX       *v ;

  counts=NULL;
  nvertices=0;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mris_average_parcellation.c,v 1.4 2011/03/02 00:04:27 nicks Exp $",
           "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  if (!strlen(sdir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
    usage_exit(1) ;

  annot_name = argv[1] ;
  hemi = argv[2] ;
  out_fname = argv[argc-1] ;

  printf("annot %s, hemi %s, out %s\n", annot_name, hemi, out_fname);
  for (num = 0, i = 3 ; i < argc-1 ; i++, num++) {
    subject_name = argv[i] ;
    sprintf(fname, "%s/%s/surf/%s.orig", sdir, subject_name, hemi) ;
    fprintf(stderr, "%d of %d: reading %s...\n", num+1, argc-4, fname) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(Gerror, "%s: MRISread(%s) failed", Progname, fname) ;
    if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
      ErrorExit(Gerror, "%s: MRISreadAnnotation(%s) failed", Progname, annot_name) ;
    if (mris->ct == NULL)
      ErrorExit(Gerror, "%s: MRISreadAnnotation(%s): no color table", Progname, annot_name) ;
    if (num == 0) // first one
    {
      nvertices = mris->nvertices ;
      counts = (int **)calloc(sizeof(int *), nvertices) ;
      if (counts == NULL)
        ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d vertex counts", Progname, nvertices) ;
      for (j = 0 ; j < nvertices ; j++) {
        counts[j] = (int *)calloc(sizeof(int), mris->ct->nentries) ;
        if (counts[j] == NULL)
          ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d counts at vertex %d",
                    Progname, mris->ct->nentries, nvertices) ;
      }
    }
    else if (mris->nvertices != nvertices)
      ErrorExit(ERROR_BADPARM, "%s: surface %d nvertices doesn't match previous ones (%d)",
                Progname, mris->nvertices, nvertices) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      int index ;
      v = &mris->vertices[vno] ;
      CTABfindAnnotation(mris->ct, v->annotation, &index);
      if (index < 0)
        continue ;
      counts[vno][index]++ ;
    }

    if (i < argc-2)
      MRISfree(&mris) ;
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    int max_j, max_count, annot ;

    v = &mris->vertices[vno] ;
    max_j = 0 ;
    max_count = counts[vno][0] ;
    for (j = 0 ; j < mris->ct->nentries ; j++) {
      if (counts[vno][j] > max_count) {
        max_count = counts[vno][j] ;
        max_j = j ;
      }
    }
    CTABannotationAtIndex(mris->ct, max_j, &annot) ;
    v->annotation = annot ;
  }
  printf("writing to %s...\n", out_fname) ;
  MRISwriteAnnotation(mris, out_fname) ;
  MRISfree(&mris) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "parcellation averaging took %d minutes and %d seconds.\n",
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
  if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
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
  printf("usage: %s [options] <annot name> <hemi> <subject> <subject>  ... <output annot>\n", Progname) ;
  exit(code) ;
}

