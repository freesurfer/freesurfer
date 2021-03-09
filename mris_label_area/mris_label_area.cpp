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

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "utils.h"
#include "timer.h"
#include "annotation.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static char *log_fname = NULL ;
static void usage_exit(int code) ;
static double MRISannotArea(MRI_SURFACE *mris, int label) ;
static int in_label = -1 ;
static int out_label = -1 ;

static int compute_pct = 0 ;

static char sdir[STRLEN] ;

int
main(int argc, char *argv[]) {
  char   **av, *subject_name, *cp, *hemi, *surf_name, *annot_name, fname[STRLEN];
  const char *name;
  int    ac, nargs, msec, minutes, label, seconds, i ;
  double area, total_area ;
  Timer start ;
  MRIS   *mris ;
  FILE   *log_fp ;

  nargs = handleVersionOption(argc, argv, "mris_label_area");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

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

  if (argc < 6)
    usage_exit(1) ;


  subject_name = argv[1] ;
  hemi = argv[2] ;
  surf_name = argv[3] ;
  annot_name = argv[4] ;

  if (strlen(sdir) == 0) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in env or cmd line",Progname) ;
    strcpy(sdir, cp) ;
  }
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject_name, hemi, surf_name) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s",Progname,
              fname) ;

  MRIScomputeMetricProperties(mris) ;
#if 0
  if (in_label >= 0)
    MRIreplaceValues(mri, mri, in_label, out_label) ;
#endif

  if (compute_pct)
    total_area = mris->total_area ;
  else
    total_area = 1 ;

  if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read annot file %s", Progname, annot_name) ;

  for (i = 5 ; i < argc ; i++) {
    label = atoi(argv[i]) ;
    name = annotation_to_name(index_to_annotation(label), NULL) ;
    printf("processing label %s (%d)...\n", name, label) ;


    area = MRISannotArea(mris, label) ;
    if (log_fname) {
      char fname[STRLEN] ;

      sprintf(fname, log_fname, label) ;
      printf("logging to %s...\n", fname) ;
      log_fp = fopen(fname, "a+") ;
      if (!log_fp)
        ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                  Progname, fname) ;
    } else
      log_fp = NULL ;

    if (compute_pct) {
      printf("%2.3f mm^2 in label %d (%s), %%%2.6f of total cortical area (%2.2f)\n",
             area, label, name,
             100.0*(float)area/(float)total_area,
             total_area) ;
      if (log_fp) {
        fprintf(log_fp,"%2.6f\n", 100.0*area/(float)total_area) ;
        fclose(log_fp) ;
      }
    } else {
      printf("%2.0f mm^2 in label %s (%d)\n", area, name, label) ;
      if (log_fp) {
        fprintf(log_fp,"%f\n", area) ;
        fclose(log_fp) ;
      }
    }
  }

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;

  if (DIAG_VERBOSE_ON)
    printf("area calculation took %d minutes and %d seconds.\n",
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
  if (!stricmp(option, "sdir")) {
    strcpy(sdir, argv[2]) ;
    printf("using SUBJECTS_DIR=%s\n", sdir) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'C':
      if (read_named_annotation_table(argv[2]) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read lookup table %s",
                  Progname, argv[2]) ;
      nargs = 1 ;
      printf("using lookup table %s...\n", argv[2]) ;
      break ;
    case 'T':
      in_label = atoi(argv[2]) ;
      out_label = atoi(argv[3]) ;
      nargs = 2 ;
      printf("translating label %d to label %d\n", in_label, out_label) ;
      break ;
    case 'L':
      log_fname = argv[2] ;
      nargs = 1 ;
      /*    fprintf(stderr, "logging results to %s\n", log_fname) ;*/
      break ;
    case 'P':
      compute_pct = 1 ;
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
  printf("usage: %s [options] <subject name> <hemi> <surf name> <annot name> <label 1> <label 2> ...\n", Progname) ;
  printf(
    "\tp             - compute brain area as a pct of all brain labels\n"
    "\tl <log fname> - log results to file (note %%d will include label #)\n"
    "\tb <brain vol> - load brain vol and use it to normalize areas\n"
  );
  exit(code) ;
}
static double
MRISannotArea(MRI_SURFACE *mris, int label) {
  int    vno, annotation ;
  VERTEX *v ;
  double area ;

  annotation = index_to_annotation(label) ;
  for (area = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->annotation == annotation)
      area += v->area ;
  }
  return(area) ;
}

