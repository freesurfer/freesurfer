/**
 * @file  mris_w_to_curv.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:34 $
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

static char vcid[] = "$Id: mris_w_to_curv.c,v 1.6 2011/03/02 00:04:34 nicks Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int navgs = 0 ;
static char sdir[STRLEN] ;
static int normalize_flag = 0 ;

int
main(int argc, char *argv[]) {
  char               **av, fname[STRLEN], *subject_name, *wfile_name,
  *cp, *curv_name, *hemi ;
  int                ac, nargs, vno ;
  MRI_SURFACE        *mris ;
  VERTEX             *v ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_w_to_curv.c,v 1.6 2011/03/02 00:04:34 nicks Exp $", "$Name:  $");
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

  if (argc < 4)
    print_help() ;

  subject_name = argv[1] ;
  hemi = argv[2] ;
  wfile_name = argv[3] ;
  curv_name = argv[4] ;

  if (strlen(sdir) == 0) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_UNSUPPORTED, "%s: must specifiy SUBJECTS_DIR in env",
                Progname) ;
    strcpy(sdir, cp) ;
  }
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject_name, hemi, ORIG_NAME) ;
  mris = MRISfastRead(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  MRISremoveTriangleLinks(mris) ;
  if (MRISreadValues(mris, wfile_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read w file %s",
              Progname, wfile_name) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = v->val ;
  }

  if (navgs > 0) {
    printf("smoothing surface tessellation for %d iterations...\n",navgs) ;
    MRISaverageCurvatures(mris, navgs) ;
  }

  if (normalize_flag) {
    printf("normalizing curvature before writing...\n") ;
    MRISnormalizeCurvature(mris, NORM_MEAN) ;
  }
  MRISwriteCurvature(mris, curv_name) ;

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
  else if (!stricmp(option, "sdir")) {
    strcpy(sdir, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using SUBJECTS_DIR=%s\n", sdir) ;
  } else switch (toupper(*option)) {
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "averaging curvature for %d iterations\n", navgs) ;
      break ;
    case 'N':
      normalize_flag = atoi(argv[2]) ;
      printf("normalizing curvature before writing\n") ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <subject> <hemi> <wfile> <curv file>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program converts w files to curv files, potentially "
          "smoothing (-a <avgs>) or normalizing (-n)\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-a <avgs>  "
          "specify # of curvature averaging iterations (def=0).\n") ;
  fprintf(stderr, "-n  normalize curvature before writing\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

