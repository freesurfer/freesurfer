/**
 * @brief main program for training a random forest on the surface
 *
 * program for building a set of features across subjects on the surface and using them to train
 * a Random Forest classifier
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
#ifdef HAVE_OPENMP
#include "romp_support.h"
#endif

static const char *class_names[] = 
{
  "Normal Cortex",
  "Dysplasia"
} ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;
static const char *hemi = "lh" ;
static const char *surf_name = "white" ;

static char sdir[STRLEN] ;


#define MAX_OVERLAYS 100
#define MAX_SUBJECTS 100

static int ndilates = 3 ;
static int nbhd_size = 0 ;
static MRI *mri_overlays[MAX_SUBJECTS][MAX_OVERLAYS] ;
static MRI_SURFACE *mris[MAX_SUBJECTS] ;
static int noverlays = 0 ;
static char *overlay_names[MAX_OVERLAYS] ;

static const char *cortex_label_name = "cortex" ;
static const char *label_name = "FCD" ;

static int ntrees = 40 ;
static int max_depth = 12 ;
static float training_fraction = .5 ;

#define ARGC_OFFSET 1
static int assemble_training_data_and_free_mris(MRI_SURFACE *mris[MAX_SUBJECTS], MRI *mri_overlays[MAX_SUBJECTS][MAX_OVERLAYS], int nsubjects, int noverlays, int **ptraining_classes, double ***ptraining_data, int *pntraining) ;

int
main(int argc, char *argv[]) {
  char          **av, *cp, fname[STRLEN], *subject, *out_fname ;
  int           ac, nargs, msec, minutes, n, seconds, nsubjects, i, sno, nfeatures, correct ;
  Timer start ;
  LABEL         *cortex_label, *training_label ;
  RANDOM_FOREST *rf ;
  double        **training_data ;
  int           *training_classes, ntraining, n_omp_threads;

  nargs = handleVersionOption(argc, argv, "mris_rf_train");
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

  if (noverlays == 0)
    ErrorExit(ERROR_NOFILE, "%s: must specify at least one training overlay with -overlay <overlay>",Progname);

#ifdef HAVE_OPENMP
  n_omp_threads = omp_get_max_threads();
  printf("\n== Number of threads available to %s for OpenMP = %d == \n",
         Progname, n_omp_threads);
#else
  n_omp_threads = 1;
#endif

  nsubjects = argc-(ARGC_OFFSET+1) ;
  printf("training random forest classifier using %d subjects\n", nsubjects) ;
  for (n = ARGC_OFFSET ; n < argc-1 ; n++)
  {
    subject = argv[n] ;
    sno = n-ARGC_OFFSET ;
    printf("processing subject %s: %d of %d\n", subject, sno+1, nsubjects) ;


    sprintf(fname, "%s/%s/label/lh.%s.label", sdir, subject, label_name) ;
    if (FileExists(fname) == 0)
    {
      sprintf(fname, "%s/%s/label/rh.%s.label", sdir, subject, label_name) ;
      if (FileExists(fname) == 0)
	ErrorExit(ERROR_NOFILE, "%s: subject %s has no training label (%s) for either hemisphere", Progname, subject,fname) ;
      hemi = "rh" ;
    }
    else
      hemi = "lh" ;

    training_label = LabelRead(NULL, fname) ;
    if (training_label == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read training label %s\n", Progname, fname) ;

    sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, hemi, surf_name) ;
    mris[sno] = MRISread(fname) ;
    if (!mris[sno])
      ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s",Progname, fname) ;
    
    MRIScomputeMetricProperties(mris[sno]) ;
    if (nbhd_size > mris[sno]->nsize)
      MRISsetNeighborhoodSizeAndDist(mris[sno], nbhd_size) ;
    sprintf(fname, "%s/%s/label/%s.%s", sdir, subject, hemi, cortex_label_name) ;
    cortex_label = LabelRead(NULL, fname) ;
    if (cortex_label == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read cortex label %s\n", Progname, fname) ;
    LabelRipRestOfSurface(cortex_label, mris[sno]) ;
    MRISreadCurvatureFile(mris[sno], "sulc") ;
    {
      int vno ;
      for (vno = 0 ; vno < mris[sno]->nvertices ; vno++)
	if (mris[sno]->vertices[vno].curv < 0)
	  mris[sno]->vertices[vno].ripflag = 1 ;
    }

    MRISclearMarks(mris[sno]) ;
    LabelFillUnassignedVertices(mris[sno], training_label, CURRENT_VERTICES);
    LabelDilate(training_label, mris[sno], ndilates, CURRENT_VERTICES) ;
    LabelMark(training_label, mris[sno]) ;
    LabelFree(&training_label) ;

    for (i = 0 ; i < noverlays ; i++)
    {
      sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, hemi, overlay_names[i]) ;
      mri_overlays[sno][i] = MRIread(fname) ;
      if (mri_overlays[sno][i] == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read overlay %s (%d)\n", Progname, fname, i) ;
    }
    LabelFree(&cortex_label) ;
  }

  nfeatures = noverlays*(nbhd_size+1) ;
  rf = RFalloc(ntrees, nfeatures, 2, max_depth, const_cast<char**>(class_names), 100) ;
  rf->feature_names = (char **)calloc(nfeatures, sizeof(char *)) ;
  for (i = 0 ; i < noverlays ; i++)
  {
    rf->feature_names[i] = (char *)calloc(strlen(overlay_names[i])+1, sizeof(char)) ;
    strcpy(rf->feature_names[i], overlay_names[i]) ;
  }
  assemble_training_data_and_free_mris(mris, mri_overlays, nsubjects, noverlays, &training_classes, &training_data, &ntraining) ;
  RFtrain(rf, 1.0, training_fraction, training_classes, training_data, ntraining);
  correct = RFcomputeOutOfBagCorrect(rf, training_classes, training_data, ntraining);
  printf("out of bag accuracy = %2.1f (%d of %d)\n", (float)correct*100.0f/ntraining, correct, ntraining) ;

  RFevaluateFeatures(rf, stdout) ;

  out_fname = argv[argc-1] ;
  printf("writing output to %s\n", out_fname) ;
  RFwrite(rf, out_fname) ;
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
  } else 
    switch (toupper(*option)) 
    {
    case '?':
    case 'U':
      usage_exit(0) ;
    break ;
    case 'L':
      label_name = argv[2] ;
      nargs = 1 ;
      printf("using label %s to define FCD\n", label_name) ;
      break ;
    case 'N':
      nbhd_size = atof(argv[2]) ;
      nargs = 1 ;
      printf("using nbhd_size = %d\n", nbhd_size) ;
      if (nbhd_size > 0)
	ErrorExit(ERROR_UNSUPPORTED, "nsize>0 not supported yet", nbhd_size) ;
      break ;
    case 'T':
      ntrees = atoi(argv[2]) ;
      nargs = 1 ;
      printf("training %d trees\n", ntrees) ;
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
  printf("usage: %s [options]  <subject 1> <subject 2> ... <output name>\n", Progname) ;
  printf(
    "\t--hemi <hemi>    - process <hemi> instead of lh\n"
    "\t--surf <surf>    - change default surface name from 'white' to <surf>\n"
  );
  exit(code) ;
}
static int
assemble_training_data_and_free_mris(MRI_SURFACE *mris[MAX_SUBJECTS], MRI *mri_overlays[MAX_SUBJECTS][MAX_OVERLAYS], 
				     int nsubjects, int noverlays, int **ptraining_classes, double ***ptraining_data, int *pntraining) 
{
  int    sno, i, ntraining, *training_classes, x, y, z, f, vno, tno ;
  double **training_data ;

  for (ntraining = sno = 0 ; sno < nsubjects ; sno++)
    ntraining += MRISvalidVertices(mris[sno]) ;

  printf("%2.1fM total training voxels found\n", (float)ntraining/(1024*1024.0)) ;
  training_classes = (int *)calloc(ntraining, sizeof(int)) ;
  if (training_classes == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not allocate %d-len training class vector", Progname, ntraining) ;
  training_data = (double **)calloc(ntraining, sizeof(double *)) ;
  if (training_data == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not allocate %d-len training data vector", Progname, ntraining) ;

  for (tno = sno = 0 ; sno < nsubjects ; sno++)
  {
    for (vno = x = 0 ; x < mri_overlays[sno][0]->width ; x++)
      for (y = 0 ; y < mri_overlays[sno][0]->height ; y++)
	for (z = 0 ; z < mri_overlays[sno][0]->depth ; z++)
	  for (f = 0 ; f < mri_overlays[sno][0]->nframes ; f++, vno++)
	  {
	    if (tno == Gdiag_no)
	      DiagBreak() ;
	    if (vno == Gdiag_no)
	      DiagBreak();
	    if (mris[sno]->vertices[vno].ripflag)  
	      continue ;   // not in cortex 
	    training_classes[tno] = mris[sno]->vertices[vno].marked ;
	    training_data[tno] = (double *)calloc(noverlays, sizeof(double)) ;
	    if (training_data[tno] == NULL)
	      ErrorExit(ERROR_NOFILE, "%s: could not allocate %d-len training data vector [%d]", Progname, noverlays, tno) ;
	    for (i = 0 ; i < noverlays ; i++)
	      training_data[tno][i] = MRIgetVoxVal(mri_overlays[sno][i], x, y, z, f) ;
	    tno++ ;
	  }
    
    for (i = 0 ; i < noverlays ; i++)
      MRIfree(&mri_overlays[sno][i]) ; 
    MRISfree(&mris[sno]) ;
  }

  *pntraining = ntraining ;
  *ptraining_classes = training_classes ;
  *ptraining_data = training_data ;

  return(NO_ERROR) ;
}
