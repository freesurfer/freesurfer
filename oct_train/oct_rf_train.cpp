/**
 * @file  oct_rf_train.c
 * @brief main program for training a cellular classifier from labeled OCT data
 *
 * program for trainin a random forest from a set of convolutional layers, and applying the
 * trained forest to an OCT image
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2014/09/26 15:27:08 $
 *    $Revision: 1.1 $
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
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "classify.h"
#include "rforest.h"
#include "voxlist.h"
#include "mrisegment.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

static double training_fraction = .25 ;
static double  feature_fraction = 1 ;
static int ntrees = 20 ;
static int max_tree_depth = 10 ;
static int dilate = 0 ;
static double pthresh = 0 ;
static int nthresh = 0 ;
static int gthresh = 0 ;
static char *class_names[] = 
{
  "Background",    // actually label 3, but will do replacement before calling traiining
  "Neuron",
  "Astrocyte",
  "Speckle",
  "Fiber"
} ;

#define OCT_BACKGROUND 0
#define OCT_NEURON     1
#define OCT_GLIA       2
#define OCT_SPECKLE    3
#define OCT_FIBER      4

static int  downsample = 2 ;
static int extract = 0 ;
static int whalf = 0 ;
static int min_training_samples = 20 ;

//static int  classifier_type = CLASSIFIER_RFOREST ;


static MRI *extract_subimage(MRI *mri_inputs)
{
  int x0, y0, w, h ;
  MRI *mri_tmp ;

  w = mri_inputs->width / extract ;
  h = mri_inputs->height / extract ;
  x0 = (mri_inputs->width - w)/2 ; y0 = (mri_inputs->height-h)/2 ;
  mri_tmp = MRIextract(mri_inputs, NULL, x0, y0, 0, w, h, 1) ;
  MRIfree(&mri_inputs) ; 
  return(mri_tmp) ;
}

int
main(int argc, char *argv[]) {
  char       **av, *int_fname, *label_fname, *out_fname ;
  int        x, y, z, f, ac, nargs, msec, minutes, seconds, nfeatures, nc, label, i, wsize ;
  Timer start ;
  MRI       *mri_inputs, *mri_labels, *mri_tmp, *mri_out ;
  float     min_label, max_label ;
  double    pval ;
  VOXEL_LIST *vl ;
  RANDOM_FOREST *rf ;
  double         **training_data, *feature ;
  int            *training_classes ;
  int xi, yi, xk, yk, fno ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: oct_rf_train.c,v 1.1 2014/09/26 15:27:08 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  setRandomSeed(-1L) ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;

  wsize = whalf*2 + 1 ;
  int_fname = argv[1] ; label_fname = argv[2] ; out_fname = argv[3];
  if (FileType(int_fname) == TEXT_FILE)
  {
    char line[STRLEN], *cp ;
    FILE *fp ;

    nfeatures = FileNumberOfEntries(int_fname) ;
    printf("%d intensity volumes specified in text file list %s\n", nfeatures, int_fname) ;
    fp = fopen(int_fname, "r") ;
    for (i = 0 ; i < nfeatures ; i++)
    {
      cp = fgetl(line, 199, fp) ;
      if (cp == NULL)
	ErrorExit(ERROR_BADPARM, "could not read feature volume name %d", i) ;
      mri_tmp = MRIread(cp) ;
      if (mri_tmp == NULL)
	ErrorExit(ERROR_NOFILE, "could not read feature volume %d: %s", i, cp) ;

      if (i == 0)
      {
	mri_inputs = MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth, MRI_FLOAT, nfeatures) ;
	MRIcopyHeader(mri_tmp, mri_inputs) ;
      }
      MRIcopyFrame(mri_tmp, mri_inputs, 0, i)  ;
      MRIfree(&mri_tmp) ;
    }
    fclose(fp) ;
  }
  else
  {
    printf("reading input intensities from %s\n", int_fname) ;
    mri_inputs = MRIread(int_fname) ;
    if (mri_inputs == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read inputs from %s", 
		Progname, int_fname) ;
  }

  if (extract > 0)
    mri_inputs = extract_subimage(mri_inputs) ;

  printf("reading input labels from %s\n", label_fname) ;
  mri_labels = MRIread(label_fname) ;
  if (mri_labels == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read inputs from %s", 
	      Progname, label_fname) ;
  if (stricmp(label_fname, int_fname) == 0)
    MRIsetValues(mri_labels, 1) ;
  if (extract > 0)
    mri_labels = extract_subimage(mri_labels) ;
  if (downsample > 0)
  {
    MRI *mri_tmp ;
    int  i ;

    for (i = 0 ; i < downsample ; i++)
    {
      mri_tmp = MRIdownsample2LabeledVolume(mri_labels, NULL) ;
      MRIfree(&mri_labels) ; mri_labels = mri_tmp ;
    }

    downsample = (int)pow(2.0, downsample) ;
    mri_tmp = MRIdownsampleN(mri_inputs, NULL,  downsample, downsample,1, 0) ;
    mri_inputs = mri_tmp ;
    MRIwrite(mri_inputs, "i.mgz") ; MRIwrite(mri_labels, "l.mgz") ;
  }

  if (dilate > 0)
  {
    MRI *mri_bin, *mri_dilate = NULL, *mri_ring, *mri_prev = NULL ;

    mri_bin = MRIcopy(mri_labels, NULL) ;
    MRIsetVoxelsWithValue(mri_bin, mri_bin, OCT_SPECKLE, 0) ;
    MRIsetVoxelsWithValue(mri_bin, mri_bin, OCT_FIBER, 0) ;
    MRIbinarize(mri_bin, mri_bin, 1, 0, 1) ;
    mri_dilate = mri_bin ;
    for (i = 0 ; i < dilate ; i++)
    {
      mri_prev = MRIcopy(mri_dilate, mri_prev) ;
      MRIdilate(mri_prev, mri_dilate) ;
    }
    mri_ring = MRIsubtract(mri_dilate, mri_prev, NULL) ;
    MRIsetVoxelsWithValue(mri_ring, mri_labels, 1, OCT_SPECKLE) ;
    MRIfree(&mri_bin) ; MRIfree(&mri_prev) ; MRIfree(&mri_ring) ;
  }

  MRInonzeroValRange(mri_labels, &min_label, &max_label) ;

  mri_tmp = MRIcloneDifferentType(mri_inputs, MRI_FLOAT) ;
//  MRIvalScale(mri_inputs, mri_tmp, 0, 1) ;
  MRIcopy(mri_inputs, mri_tmp) ;
  MRIfree(&mri_inputs) ; mri_inputs = mri_tmp ;
  vl = VLSTcreate(mri_labels, 1, 10000, NULL, 0, 0) ;

  nfeatures = mri_inputs->nframes*wsize*wsize ;
  rf = RFalloc(ntrees, nfeatures, max_label+1, max_tree_depth, class_names, 10) ;
  
  feature = (double *)calloc(nfeatures, sizeof(double)) ;
  training_classes = (int *)calloc(vl->nvox, sizeof(int)) ;
  training_data = (double **)calloc(vl->nvox, sizeof(double *)) ;
  for  (i = 0 ; i < vl->nvox ; i++)
  {
    training_data[i] = (double *)calloc(nfeatures, sizeof(double)) ;
    if (training_data[i] == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d of %d training vector\n", i, vl->nvox) ;
    x = vl->xi[i] ;  y = vl->yi[i] ;  z = vl->zi[i] ;
    label = MRIgetVoxVal(mri_labels, x, y, z, 0) ;
    if (label == OCT_SPECKLE) 
      label = 0 ;    // change background to class 0
    training_classes[i] = label ;
    for (fno = f = 0 ; f < mri_inputs->nframes ; f++)
    {
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
	xi = mri_inputs->xi[x+xk] ;
	for (yk = -whalf ; yk <= whalf ; yk++, fno++)
	{
	  yi = mri_inputs->yi[y+yk] ;
	  training_data[i][fno] = MRIgetVoxVal(mri_inputs, xi, yi, z, f) ;
	}
      }
    }
  }
  RFtrain(rf, feature_fraction, training_fraction, training_classes, training_data, vl->nvox);
  if (min_training_samples > 0)
    RFpruneTree(rf, min_training_samples) ;
  nc = RFcomputeOutOfBagCorrect(rf, training_classes, training_data, vl->nvox);
  printf("%d of %d correct (%2.1f%%)\n", nc, vl->nvox, 100.0*nc/vl->nvox) ;
  RFwrite(rf, out_fname) ;
  
  FileNameRemoveExtension(out_fname, out_fname) ; strcat(out_fname, ".mgz") ;
  mri_out = MRIallocSequence(mri_inputs->width, mri_inputs->height, mri_inputs->depth, MRI_FLOAT, 2) ; MRIcopyHeader(mri_inputs, mri_out) ;
  for (x = 0 ; x < mri_inputs->width ; x++)
  {
    for (y = 0 ; y < mri_inputs->height ; y++)
    {
      if ((MRIgetVoxVal(mri_inputs, x, y, 0, 0) < .0001)  /*| (MRIcountValInNbhd(mri_inputs, 7, x, y, 0, 0)  > 0)*/)
	continue ; // don't bother processing background
      
      if (x == Gx && y == Gy)
	DiagBreak() ;
      for (fno = f = 0 ; f < mri_inputs->nframes ; f++)
      {
	for (xk = -whalf ; xk <= whalf ; xk++)
	{
	  xi = mri_inputs->xi[x+xk] ;
	  for (yk = -whalf ; yk <= whalf ; yk++, fno++)
	  {
	    yi = mri_inputs->yi[y+yk] ;
	    feature[fno] = MRIgetVoxVal(mri_inputs, xi, yi, 0, f) ;
	  }
	}
      }
      label = RFclassify(rf, feature, &pval, -1) ;
      if (pthresh > 0)
      {
	if ((label != OCT_NEURON && label != OCT_GLIA) || pval < pthresh)
	  label = 0 ;
      }
      MRIsetVoxVal(mri_out, x, y, 0, 0, label) ;
      MRIsetVoxVal(mri_out, x, y, 0, 1, pval) ;
    }
  }
  if (nthresh > 0)
  {
    MRI_SEGMENTATION *mseg_neuron ;

    printf("applying neuron area threshold %d to segmentation output\n", nthresh) ;
    mseg_neuron = MRIsegment(mri_out, OCT_NEURON, OCT_NEURON) ;
    MRIeraseSmallSegments(mseg_neuron, mri_out, nthresh) ;
    MRIsegmentFree(&mseg_neuron) ;
  }
  if (gthresh > 0)
  {
    MRI_SEGMENTATION *mseg_glia ;

    printf("applying glia area threshold %d to segmentation output\n", gthresh) ;
    mseg_glia = MRIsegment(mri_out, OCT_GLIA, OCT_GLIA) ;
    MRIeraseSmallSegments(mseg_glia, mri_out, gthresh) ;
    MRIsegmentFree(&mseg_glia) ;
  }
  printf("writing outputs to %s\n", out_fname) ;
  MRIwrite(mri_out, out_fname) ;
  for  (i = 0 ; i < vl->nvox ; i++)
    free(training_data[i]) ;
  free(feature) ; free(training_data) ; free(training_classes) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ; minutes = seconds / 60 ; seconds = seconds % 60 ;
  fprintf(stderr, "OCT segmentation  took %d minutes and %d seconds.\n",  minutes, seconds) ;
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
  if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
  }
  else if (!stricmp(option, "trees"))
  {
    ntrees = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d trees in training random forest\n", ntrees) ;
  }
  else if (!stricmp(option, "dilate"))
  {
    dilate = atoi(argv[2]) ;
    nargs = 1 ;
    printf("dilating labels %d times and using outer ring as background\n", dilate) ;
  }
  else if (!stricmp(option, "depth"))
  {
    max_tree_depth = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d max tree depth in training random forest\n", max_tree_depth) ;
  }
  else if (!stricmp(option, "tfrac"))
  {
    training_fraction = atof(argv[2]) ;
    nargs = 1 ;
    printf("using %2.3f training fraction in training random forest\n", training_fraction) ;
  }
  else switch (toupper(*option)) {
    case 'D':
      downsample = atoi(argv[2]) ;
      nargs = 1 ;
      printf("downsampling inputs %d times before training\n", downsample) ;
      break ;
    case 'X':
      extract = atof(argv[2]) ;
      nargs = 1 ; 
      printf("setting extract to %d\n", extract) ;
      break ;
    case 'V':
      Gdiag_no = atof(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      usage_exit(0) ;
    break ;
    case 'W':
      whalf = atoi(argv[2]) ;
      printf("using whalf = %d, wsize = %d x %d\n", whalf, whalf*2+1, whalf*2+1) ;
      nargs  = 1 ;
      break ;
    case 'G':
      gthresh = atoi(argv[2]) ;
      printf("using glia segmentation area threshold %d\n", gthresh) ;
      nargs  = 1 ;
      break ;
    case 'N':
      nthresh = atoi(argv[2]) ;
      printf("using neuron segmentation area threshold %d\n", nthresh) ;
      nargs  = 1 ;
      break ;
    case 'P':
      pthresh = atof(argv[2]) ;
      printf("using segmentation p-value threshold %2.1f\n", pthresh) ;
      nargs  = 1 ;
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
  printf("usage: %s [options] <intensity volume> <label volume> <output RBM>\n",
         Progname) ;
  exit(code) ;
}



