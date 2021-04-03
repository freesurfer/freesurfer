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
#include "label.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  compute_thickness_stats(MRI_SURFACE *mris, LABEL *area,
                                    double *pmean_w, double *pvar_w,
                                    double *pmean_thick, double *pvar_thick,
                                    int *pn) ;

const char *Progname ;

static char subjects_dir[STRLEN] ;
static int  avgs = 0 ;
static FILE *log_fp = NULL ;

int
main(int argc, char *argv[]) {
  char               **av, fname[STRLEN], env_string[STRLEN],
  *cp, *subject_name, *hemi, *thickness_name, *wfile_name, *label_name ;
  int                ac, nargs, i, npoints ;
  MRI_SURFACE        *mris ;
  double             mean_w, var_w, mean_thick, var_thick ;
  LABEL              *area ;

  nargs = handleVersionOption(argc, argv, "mris_thickness_comparison");
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

  if (argc < 5)
    print_help() ;

  if (strlen(subjects_dir) == 0) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_UNSUPPORTED, "%s: must specify -sdir or SUBJECTS_DIR in env",Progname) ;
    strcpy(subjects_dir, cp) ;
  }

  int req = snprintf(env_string, STRLEN, "SUBJECTS_DIR=%s", subjects_dir) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  putenv(env_string) ;
  cp = getenv("SUBJECTS_DIR") ;
  subject_name = argv[1] ;
  hemi = argv[2] ;
  thickness_name = argv[3] ;
  wfile_name = argv[4] ;

  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", subjects_dir, subject_name, hemi, ORIG_NAME) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  if (MRISreadCurvature(mris, thickness_name) != NO_ERROR)
    ErrorExit(ERROR_BADFILE, "%s: could not read thickness data from %s",
              Progname, thickness_name) ;

  MRISaverageCurvatures(mris, avgs) ;
  if (MRISreadValues(mris, wfile_name) != NO_ERROR)
    ErrorExit(ERROR_BADFILE, "%s: could not read wfile from %s",
              Progname, wfile_name) ;

  for (i = 5 ; i < argc ; i++) {
    label_name = argv[i] ;
    printf("reading label file %s\n", label_name) ;
    area = LabelRead(subject_name, label_name) ;
    if (area == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not read label file %s\n", Progname,label_name) ;
    compute_thickness_stats(mris, area, &mean_w, &var_w,
                            &mean_thick, &var_thick, &npoints) ;
    printf("%s: manual thickness = %2.2f+-%2.2f (%d points)\n",
           label_name,mean_w, sqrt(var_w), npoints) ;
    printf("%s: auto   thickness = %2.2f+-%2.2f\n",  label_name,mean_thick, sqrt(var_thick)) ;
    printf("\n") ;
    if (log_fp)
      fprintf(log_fp, "%s %f %f %f %f %d\n",
              label_name, mean_w, sqrt(var_w), mean_thick, sqrt(var_thick), npoints) ;
    LabelFree(&area) ;
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
  else switch (toupper(*option)) {
    case 'L':
      log_fp = fopen(argv[2], "w") ;
      if (log_fp == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not open log file %s", Progname, argv[2]) ;
      nargs = 1 ;
      break ;
    case 'A':
      avgs = atoi(argv[2]) ;
      printf("smoothing auto thicknesses %d times\n", avgs) ;
      nargs = 1 ;
      break ;
    case 'S':
      strcpy(subjects_dir, argv[2]) ;
      nargs = 1 ;
      printf("using SUBJECTS_DIR=%s\n", subjects_dir) ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
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
          "usage: %s [options] <subject> <hemi> <thickness file> <w file> "
          "<label1> <label2>...\n", Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

static int  compute_thickness_stats(MRI_SURFACE *mris, LABEL *area,
                                    double *pmean_w, double *pvar_w,
                                    double *pmean_thick, double *pvar_thick,
                                    int *pn) {
  int    npoints ;
  double mean_w, var_w, mean_thick, var_thick ;
  VERTEX       *v ;
  int     vno, n ;

  mean_w = var_w = mean_thick = var_thick = 0 ;
  for (npoints = n = 0 ; n < area->n_points ; n++) {
    vno = area->lv[n].vno ;
    v = &mris->vertices[vno] ;
    if (v->ripflag || FZERO(v->val))
      continue ;
    npoints++ ;
    mean_w += v->val ;
    var_w += (v->val*v->val) ;
    mean_thick += v->curv ;
    var_thick += (v->curv*v->curv) ;
  }

  mean_w /= (float)npoints ;
  mean_thick /= (float)npoints ;
  var_w = var_w / npoints - mean_w*mean_w ;
  var_thick = var_thick / npoints - mean_thick*mean_thick ;

  *pmean_w = mean_w ;
  *pmean_thick = mean_thick ;
  *pvar_w = var_w ;
  *pvar_thick = var_thick ;
  *pn = npoints ;
  return(NO_ERROR) ;
}

