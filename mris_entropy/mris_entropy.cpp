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


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static char *log_fname = NULL ;
static int navgs = 0 ;
static char sdir[STRLEN] ;
static int normalize_flag = 0 ;

int
main(int argc, char *argv[]) {
  char               **av, fname[STRLEN], *subject_name, *wfile_name,
  *cp, *hemi ;
  int                ac, nargs, vno ;
  MRI_SURFACE        *mris ;
  VERTEX             *v ;
  double             entropy, total_len, min_w ;

  nargs = handleVersionOption(argc, argv, "mris_entropy");
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
    print_help() ;

  subject_name = argv[1] ;
  hemi = argv[2] ;
  wfile_name = argv[3] ;

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

  if (MRISreadValues(mris, wfile_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read w file %s",
              Progname, wfile_name) ;

  if (navgs > 0) {
    printf("smoothing surface values for %d iterations...\n",navgs) ;
    MRISaverageVals(mris, navgs) ;
  }

  min_w = fabs(mris->vertices[0].val) ;
  for (total_len = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val = fabs(v->val) ;
    total_len += (v->val*v->val) ;
    if (v->val < min_w)
      min_w = v->val ;
  }

  total_len = sqrt(total_len) ;
  if (FZERO(total_len))
    ErrorExit(ERROR_BADPARM, "total vector len = 0, entropy = 0 (trivial case)") ;

  for (entropy = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (DZERO(v->val))
      continue ;
    v->val = v->val / total_len ;
    entropy += v->val * log2(v->val) ;
  }

  entropy = -entropy ;
  printf("total entropy = %f\n", entropy) ;
  if (log_fname) {
    FILE *fp ;

    fp = fopen(log_fname, "a") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open log file %s\n", Progname,log_fname) ;
    fprintf(fp, "%f\n", entropy) ;
    fclose(fp) ;
  }
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
    case 'L':
      log_fname = argv[2] ;
      printf("logging results to %s...\n", log_fname) ;
      nargs = 1 ;
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
          "\nThis program computes the entropy of a surface activation pattern "
          "smoothing (-a <avgs>) or normalizing (-n)\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-a <avgs>  "
          "specify # of curvature averaging iterations (def=0).\n") ;
  fprintf(stderr, "-n  normalize curvature before writing\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

