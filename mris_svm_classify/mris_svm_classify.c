/**
 * @file  mris_svm_classify.c
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "mrisurf.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "svm.h"
#include "label.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;

static char subjects_dir[STRLEN] ;
static float true_class = 0 ;
static char *log_fname = NULL ;
static int navgs = 0 ;
static LABEL *area = NULL ;
static char *label_name = NULL ;

#define MAX_ANNOTATIONS 100
static int annotations[MAX_ANNOTATIONS] ;
static char *anames[MAX_ANNOTATIONS] ;
static int nannotations = 0 ;
static char *annot_name = "aparc" ;

int
main(int argc, char *argv[]) {
  char         **av, fname[STRLEN], *input_name, *subject_name, *cp,*hemi,
  *svm_name, *surf_name, *output_subject_name ;
  int          ac, nargs, vno ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  MRI_SURFACE  *mris ;
  SVM          *svm ;
  double       classification ;
  float        *inputs ;
  MRI_SP       *mrisp ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_svm_classify.c,v 1.6 2011/03/02 00:04:34 nicks Exp $", "$Name:  $");
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

  if (!strlen(subjects_dir)) /* hasn't been set on command line */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
                Progname);
    strcpy(subjects_dir, cp) ;
  }
  if (argc < 7)
    usage_exit(1) ;

  subject_name = argv[1] ;
  hemi = argv[2] ;
  surf_name = argv[3] ;
  input_name = argv[4] ;
  output_subject_name = argv[5] ;
  svm_name = argv[6] ;

  printf("reading svm from %s...\n", svm_name) ;
  svm = SVMread(svm_name) ;
  if (!svm)
    ErrorExit(ERROR_NOFILE, "%s: could not read classifier from %s",
              Progname, svm_name) ;
  if (log_fname != NULL)
    printf("logging results to %s, true_class = %s\n",
           log_fname, true_class > 0 ? svm->class1_name : svm->class2_name) ;

  sprintf(fname, "%s/%s/surf/%s.%s", subjects_dir,subject_name,hemi,surf_name);
  printf("reading surface from %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s for %s",
              Progname, fname, subject_name) ;
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

  if (MRISreadCurvature(mris, input_name) != NO_ERROR)
    ErrorExit(ERROR_BADPARM, "%s: could not read curvature from %s", input_name) ;

    if (nannotations > 0)
    {
      int vno, a, found ;
      VERTEX *v ;
      
      if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, 
                  "%s: could not read annot file %s for subject %s",
                  Progname, annot_name, subject_name) ;
      for (a = 0 ; a < nannotations ; a++)
      {
        int index ;
        
        CTABfindName(mris->ct, anames[a], &index) ;
        CTABannotationAtIndex(mris->ct, index, &annotations[a]) ;
        printf("mapping annot %s to %d\n",
               anames[a], annotations[a]) ;
      }
      // rip all vertices that don't have one of the specified annotations
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
          continue ;
        found = 0 ;
        for (a = 0 ; a < nannotations ; a++)
          if (v->annotation == annotations[a])
            found = 1 ;
        if (found == 0)
          v->ripflag = 1 ;
      }
    }
  if (navgs > 0)
    MRISaverageCurvatures(mris, navgs) ;

  mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
  MRISfree(&mris) ;

  /* read in output surface */
  sprintf(fname, "%s/%s/surf/%s.%s", subjects_dir,output_subject_name,hemi,surf_name);
  printf("reading output surface from %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s for %s",
              Progname, fname, output_subject_name) ;
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISfromParameterization(mrisp, mris, 0) ;

  if (label_name) {
    area = LabelRead(output_subject_name, label_name) ;
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label %s", Progname,
                label_name) ;
    MRISmaskNotLabel(mris, area) ;
  } else
    area = NULL ;
  if (mris->nvertices != svm->ninputs)
    ErrorExit(ERROR_BADPARM, "%s: svm input (%d) does not match # of "
              "surface vertices (%d)",
              Progname, svm->ninputs, mris->nvertices);

  inputs = (float *)calloc(mris->nvertices, sizeof(float)) ;
  if (!inputs)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d input vector",
              Progname, mris->nvertices) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
    inputs[vno] = mris->vertices[vno].curv ;
  classification = SVMclassify(svm, inputs) ;
  printf("classification %f, class = %s",classification,
         classification > 0 ? svm->class1_name : svm->class2_name) ;
  if (true_class != 0)
    printf(", %s", true_class*classification>0 ? "CORRECT" : "INCORRECT") ;
  printf("\n") ;

  if (log_fname) {
    FILE *fp ;
    fp = fopen(log_fname, "a") ;
    if (!fp)
      ErrorExit(ERROR_BADPARM, "%s: could not open log file %s", log_fname) ;
    fprintf(fp, "%-30.30s %s %d %f %f\n",
            subject_name, hemi, (true_class*classification)>0, classification,
            true_class) ;
    fclose(fp) ;
  }
  free(inputs) ;
  MRISfree(&mris) ;
  SVMfree(&svm) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("classification took %d minutes and %d seconds.\n", minutes, seconds) ;
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
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  } else if (!stricmp(option, "aname")) {
    annot_name = argv[2] ;
    printf("using %s as name of annotation file\n", annot_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "annot")) {
    if (nannotations >= MAX_ANNOTATIONS)
      ErrorExit(ERROR_NOMEMORY, "%s: too many annotations specified (%d)",
                Progname, nannotations) ;
    anames[nannotations] = argv[2] ;
    printf("using annotation %s\n", anames[nannotations]) ;
    nannotations++;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'A':
      navgs = atoi(argv[2]) ;
      printf("averaging input values %d times before classifying...\n", navgs) ;
      nargs = 1 ;
      break ;
    case 'C':
      true_class = atof(argv[2]) ;
      if (true_class>1.0f)
        true_class = -1.0f ;
      log_fname = argv[3] ;
      nargs = 2 ;
      break ;
    case 'L':
      label_name = argv[2] ;
      printf("masking label %s\n", label_name) ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
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
  printf("usage: %s [options] <subject> <hemi> <canon surf name> <input curv name> "
         "<output subject> <classifier>\n", Progname) ;
  exit(code) ;
}

