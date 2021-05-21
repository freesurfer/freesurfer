/**
 * @brief main program for labeling a surface based on a previously trained random forest on the surface
 *
 * program for labeling a surfaced based on a set of features across subjects on the surface and that had
 * been used to train a Random Forest classifier (see mris_rf_train)
 */
/*
 * Original Author: Bruce Fischl
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
#include "rforest.h"
#include "rfutils.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;
static const char *hemi = "lh" ;
static const char *surf_name = "white" ;

static char sdir[STRLEN] ;


#define MAX_OVERLAYS 100

static int ndilates = 3 ;
static int nbhd_size = 0 ;
static MRI *mri_overlays[MAX_OVERLAYS] ;
static int noverlays = 0 ;
static char *overlay_names[MAX_OVERLAYS] ;

static const char *cortex_label_name = "cortex" ;
static const char *label_name = "FCD" ;

int
main(int argc, char *argv[]) {
  char          **av, *cp, fname[STRLEN], *subject, *out_fname ;
  int           ac, nargs, msec, minutes, seconds, i, nfeatures, vno ;
  Timer start ;
  LABEL         *cortex_label, *training_label ;
  RANDOM_FOREST *rf ;
  double        *feature, pval ;
  int           classnum ;
  MRI_SURFACE   *mris ;
  MRI           *mri_labels ;
  VERTEX        *v ;

  nargs = handleVersionOption(argc, argv, "mris_rf_label");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ; setRandomSeed(0L);
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;


  if (strlen(sdir) == 0) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in env or cmd line",Progname) ;
    strcpy(sdir, cp) ;
  }

  subject = argv[1] ;
  rf = RFread(argv[2]) ;
  noverlays = rf->nfeatures ;

  int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", sdir, subject, hemi, surf_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s",Progname, fname) ;
    
  MRIScomputeMetricProperties(mris) ;
  req = snprintf(fname, STRLEN, "%s/%s/label/%s.%s", sdir, subject, hemi, cortex_label_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  cortex_label = LabelRead(NULL, fname) ;
  if (cortex_label == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read cortex label %s\n", Progname, fname) ;
  LabelRipRestOfSurface(cortex_label, mris) ;

  req = snprintf(fname, STRLEN, "%s/%s/label/%s.%s", sdir, subject, hemi, label_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (FileExists(fname))
  {
    training_label = LabelRead(NULL, fname) ;
    if (training_label == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read training label %s\n", Progname, fname) ;
    
    MRISclearMarks(mris) ;
    LabelFillUnassignedVertices(mris, training_label, CURRENT_VERTICES);
    LabelDilate(training_label, mris, ndilates, CURRENT_VERTICES) ;
    LabelMark(training_label, mris) ;
  }
  else 
    training_label = NULL ;

  for (i = 0 ; i < noverlays ; i++)
  {
    int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", sdir, subject, hemi, rf->feature_names[i]) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    mri_overlays[i] = MRIread(fname) ;
    if (mri_overlays[i] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read overlay %s (%d)\n", Progname, fname, i) ;
  }
  mri_labels = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 2) ;
  nfeatures = noverlays*(nbhd_size+1) ;
  feature = (double *)calloc(nfeatures, sizeof(double)) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ; 
    if (v->ripflag)
      continue ;
    for (i = 0 ; i < noverlays ; i++)
      feature[i] = MRIgetVoxVal(mri_overlays[i], vno, 0, 0, 0) ;
    if (vno == Gdiag_no)
    {
      printf("classifying vertex %d\n", vno) ;
      for (i = 0 ; i < noverlays ; i++)
	printf("\t%s: %2.1f\n", rf->feature_names[i], feature[i]) ;
      Gdiag |= DIAG_VERBOSE ;
      DiagBreak() ;
    }

    classnum = RFclassify(rf, feature, &pval, v->marked);
    if (vno == Gdiag_no)
    {
      printf("\tclass = %d, pval = %2.2f\n", classnum, pval) ;
      Gdiag &= ~DIAG_VERBOSE ;
      DiagBreak() ;
    }
    if (classnum == 1)
      DiagBreak() ;
    MRIsetVoxVal(mri_labels, vno, 0, 0, 0, classnum) ;
    if (classnum == 0)
      pval = 1-pval ;
    MRIsetVoxVal(mri_labels, vno, 0, 0, 1, pval) ;
  }

  LabelFree(&cortex_label) ;
  if (training_label)
    LabelFree(&training_label) ;

  out_fname = argv[argc-1] ;
  printf("writing output to %s\n", out_fname) ;
  MRIwrite(mri_labels, out_fname) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;

  printf("random forest training took %d minutes and %d seconds.\n", minutes, seconds) ;

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
  } else if (!stricmp(option, "hemi")) {
    hemi = argv[2] ;
    printf("processing hemisphere = %s\n", hemi) ;
    nargs = 1 ;
  } else if (!stricmp(option, "surf")) {
    surf_name = argv[2] ;
    printf("using surface %s as input\n", surf_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "overlay")) {
    if (noverlays >= MAX_OVERLAYS)
      ErrorExit(ERROR_NOMEMORY, "%s: too many overlays specified (%d max)\n", Progname, MAX_OVERLAYS) ;
    overlay_names[noverlays] = argv[2] ;
    printf("setting overlay[%d] = %s\n", noverlays,overlay_names[noverlays]) ;
    noverlays++ ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  case 'N':
    nbhd_size = atof(argv[2]) ;
    nargs = 1 ;
    printf("using nbhd_size = %d\n", nbhd_size) ;
    if (nbhd_size > 0)
      ErrorExit(ERROR_UNSUPPORTED, "nsize>0 not supported yet", nbhd_size) ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
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
  printf("usage: %s [options]  <subject> <RF classifier> <output name>\n", Progname) ;
  printf(
    "\t--hemi <hemi>    - process <hemi> instead of lh\n"
    "\t--surf <surf>    - change default surface name from 'white' to <surf>\n"
  );
  exit(code) ;
}
