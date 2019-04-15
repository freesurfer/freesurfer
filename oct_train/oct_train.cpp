/**
 * @file  oct_train.c
 * @brief main program for training a cellular classifier from labeled OCT data
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2014/04/17 13:30:03 $
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
#include "rbm.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;
static int train_xor(RBM_PARMS *parms)  ;


static int ngroups = 5 ;          // number of filters to learn
static int nhidden[MAX_LAYERS] ; // only for DBNs
static int ksizes[MAX_LAYERS] ;
static  int ksize = 11 ;          // size of filters
static int  downsample ;
static double momentum = .5 ;
static RBM_PARMS parms ;
static int extract = 0 ;
static int test_on_xor = 0 ;
static int rbm_input_type = RBM_TYPE_CONTINUOUS_INPUTS ;
static int force_dbn = 0 ;
static int force_cdbn = 0 ;

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

static char isXdigit(char c)
{
  return(c == '.' || isdigit(c)) ;
}
int
main(int argc, char *argv[]) {
  char       **av, *int_fname, *label_fname, *out_fname ;
  int        ac, nargs, msec, minutes, seconds ;
  Timer start ;
  MRI       *mri_inputs, *mri_labels, *mri_tmp  ;
  RBM       *rbm ;
  DBN       *dbn ;
  CDBN      *cdbn ;
  float     min_label, max_label ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: oct_train.c,v 1.1 2014/04/17 13:30:03 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  setRandomSeed(-1L) ;
  parms.nlayers = 1 ;
  parms.write_iterations = 100 ;
  {
    int l ;
    for (l = 0 ; l < MAX_LAYERS ; l++)
    {
      parms.training_rates[l] = .0001 ;
      parms.sparsity[l] = .003 ;
      parms.momentum[l] = momentum ;
      parms.weight_decays[l] = .00001 ;
      parms.l_sparsity[l] = 1 ;
    }
  }
  parms.mini_batch_size = 1000 ;
  parms.discriminative_training = 0 ;
  parms.label_trate = .0001 ;
  parms.l_label = 1 ;
  parms.held_out = 1000 ;
  parms.batches_per_step = 10 ;
  parms.nsteps = 5000 ;
  parms.learn_variance = 0 ;
  parms.nclasses = 3 ;
  parms.sparsity_decay = .9 ;
  parms.max_no_progress = 500 ;
  parms.ksize = ksize ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;

  int_fname = argv[1] ; label_fname = argv[2] ; out_fname = argv[argc-1];
  FileNameRemoveExtension(out_fname, parms.base_name) ;
  printf("reading input intensities from %s\n", int_fname) ;
  mri_inputs = MRIread(int_fname) ;
  if (mri_inputs == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read inputs from %s", 
	      Progname, int_fname) ;

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

  if (test_on_xor)
  {
    train_xor(&parms) ;
    exit(0) ;
  }
  MRInonzeroValRange(mri_labels, &min_label, &max_label) ;

  mri_tmp = MRIcloneDifferentType(mri_inputs, MRI_FLOAT) ;
//  MRIvalScale(mri_inputs, mri_tmp, 0, 1) ;
  MRIcopy(mri_inputs, mri_tmp) ;
  if (rbm_input_type == RBM_TYPE_BINARY_INPUTS)
  {
    char fname[STRLEN];
    MRIbinarize(mri_inputs, mri_inputs, 0.5, 0, 1) ;
    sprintf(fname, "%s.bin.mgz", parms.base_name) ;
    MRIwrite(mri_inputs, fname) ;
  }
    
  MRIfree(&mri_inputs) ; mri_inputs = mri_tmp ;

  if ((parms.nlayers == 1) && !force_dbn && !force_cdbn)
  {
    printf("training RBM\n") ;
    rbm = RBMalloc(rbm_input_type, parms.ksize*parms.ksize, ngroups,max_label+1, RBM_INPUT_IMAGE) ;
    RBMtrainFromImage(rbm, mri_inputs, mri_labels, &parms) ;
    RBMwriteNetwork(rbm,  -1, &parms, -1) ;
    {
      MRI *mri_tmp, *mri_seg ;
      char fname[STRLEN] ;
      
      mri_tmp = RBMreconstruct(rbm, mri_inputs, NULL, &mri_seg, &parms) ;
      sprintf(fname, "%s.recon.mgz", parms.base_name) ;
      printf("writing reconstructed image to %s\n", fname) ;
      MRIwrite(mri_tmp, fname) ;
      sprintf(fname, "%s.aseg.mgz", parms.base_name) ;
      printf("writing labeled image to %s\n", fname) ;
      MRIwrite(mri_seg, fname) ;
    }
    printf("writing trained RBM to %s\n", out_fname) ;
    RBMwrite(rbm, out_fname) ;
    RBMfree(&rbm) ;
  }
  else if (!force_cdbn)
  {
    printf("training DBN\n") ;
    dbn = DBNalloc(rbm_input_type, parms.nlayers, parms.ksize*parms.ksize, nhidden,max_label+1, RBM_INPUT_IMAGE) ;
    DBNtrainFromImage(dbn, mri_inputs, mri_labels, &parms) ;
    DBNwriteNetwork(dbn,  -1, &parms) ;
    {
      MRI *mri_tmp, *mri_seg ;
      char fname[STRLEN] ;
      
      mri_tmp = DBNreconstruct(dbn, mri_inputs, NULL, &mri_seg, &parms) ;
      sprintf(fname, "%s.recon.mgz", parms.base_name) ;
      printf("writing reconstructed image to %s\n", fname) ;
      MRIwrite(mri_tmp, fname) ;
      sprintf(fname, "%s.aseg.mgz", parms.base_name) ;
      printf("writing labeled image to %s\n", fname) ;
      MRIwrite(mri_seg, fname) ;
    }
    printf("writing trained RBM to %s\n", out_fname) ;
    DBNwrite(dbn, out_fname) ;
    DBNfree(&dbn) ;
  }
  else   // convolutional deep belief network
  {
    int layer ;

    printf("training CDBN\n") ;
    cdbn = CDBNalloc(rbm_input_type, parms.nlayers, ksizes, nhidden,max_label+1, mri_inputs) ;
    CDBNtrainFromImage(cdbn, mri_inputs, mri_labels, &parms) ;
    for (layer = 0 ; layer < cdbn->nlayers ; layer++)
      CDBNwriteNetwork(cdbn,  -1, &parms, layer) ;
    {
      MRI *mri_tmp, *mri_seg ;
      char fname[STRLEN] ;

      mri_tmp = CDBNcreateOutputs(cdbn,  &parms, mri_inputs, 0, cdbn->nlayers-1, &mri_seg) ;

      sprintf(fname, "%s.recon.mgz", parms.base_name) ;
      printf("writing reconstructed image to %s\n", fname) ;
      MRIwrite(mri_tmp, fname) ;
      sprintf(fname, "%s.aseg.mgz", parms.base_name) ;
      printf("writing labeled image to %s\n", fname) ;
      MRIwrite(mri_seg, fname) ;
    }
    printf("writing trained CDBN to %s\n", out_fname) ;
//    DBNwrite(dbn, out_fname) ;
//    DBNfree(&dbn) ;
  }

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "RBM training took %d minutes and %d seconds.\n", 
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
  int  nargs = 0, layer ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "wd"))
  {
    for (layer = 0 ; layer < MAX_LAYERS ; layer++)
    {
      if (isXdigit(*argv[2+layer]))
      {
	parms.weight_decays[layer] = atof(argv[2+layer]) ;
	printf("using layer %d weight decay %f\n", layer, parms.weight_decays[layer]) ;
	nargs++ ;
      }
      else
	break ;
    }
  }
  else if (!stricmp(option, "mb"))
  {
    parms.mini_batch_size = atoi(argv[2]) ;
    parms.batches_per_step = atoi(argv[3]) ;
    printf("using mini batch size=%d and batches/step = %d\n", parms.mini_batch_size, parms.batches_per_step) ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "nb") || !stricmp(option, "np"))
  {
    parms.max_no_progress = atoi(argv[2]) ;
    printf("set max no progress to %d\n", parms.max_no_progress) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "ho"))
  {
    parms.held_out = atoi(argv[2]) ;
    printf("holding out %d samples for testing\n", parms.held_out) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "sparsity"))
  {
    for (layer = 0 ; layer < MAX_LAYERS ; layer++)
    {
      if (isXdigit(*argv[2+2*layer]))
      {
	parms.sparsity[layer] = atof(argv[2+2*layer]) ;
	parms.l_sparsity[layer] = atof(argv[2+2*layer+1]) ;
	printf("setting sparsity target for layer %d to %f with weight %2.3f\n", layer, parms.sparsity[layer], parms.l_sparsity[layer]) ;
	nargs += 2 ;
      }
      else
	break ;
    }
  }
  else if (!stricmp(option, "variance"))
  {
    parms.variance = atof(argv[2]) ;
    parms.learn_variance = -1 ;
    printf("setting variance to %f (std=%f)\n", parms.variance, sqrt(parms.variance)) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "xor"))
  {
    test_on_xor = 1 ;
    printf("training on synthesized xor data\n") ;
  }
  else if (!stricmp(option, "layers") || !stricmp(option, "nlayers") || !stricmp(option, "dbn"))
  {
    int l ;
    force_dbn = 1 ;
    parms.nlayers = atoi(argv[2]) ;
    nargs = 1+parms.nlayers ;
    for (l = 0 ; l < parms.nlayers ; l++)
    {
      nhidden[l] = atoi(argv[3+l]) ;
      printf("hidden layer %d: %d units\n", l+1, nhidden[l]) ;
    }
      
    printf("training a %d-layer Deep Belief Network (DBN)\n", parms.nlayers) ;
  }
  else if (!stricmp(option, "cdbn"))
  {
    int l ;
    force_cdbn = 1 ;
    parms.nlayers = atoi(argv[2]) ;
    nargs = 1+2*parms.nlayers ;
    for (l = 0 ; l < parms.nlayers ; l++)
    {
      ksizes[l] =  atoi(argv[3+2*l]);
      nhidden[l] = atoi(argv[3+2*l+1]) ;
      printf("hidden layer %d: ksize: %d groups: %d\n", l+1, ksizes[l], nhidden[l]) ;
    }
      
    printf("training a %d-layer Deep Belief Network (DBN)\n", parms.nlayers) ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
  }
  else if (!stricmp(option, "debug"))
  {
    parms.debug = atoi(argv[2]) ;
    nargs = 1 ;
    printf("setting debug level to %d\n", parms.debug) ;
  }
  else if (!stricmp(option, "disc"))
  {
    parms.discriminative_training = atoi(argv[2]) ;
    nargs = 1 ;
    printf("setting discriminative traiing to %d\n", parms.discriminative_training) ;
  }
  else switch (toupper(*option)) {
    case 'W':
      parms.write_iterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("setting write iterations to  %d\n", parms.write_iterations) ;
      break ;
    case 'D':
      downsample = atoi(argv[2]) ;
      nargs = 1 ;
      printf("downsampling inputs %d times before training\n", downsample) ;
      break ;
    case 'B':
      rbm_input_type = RBM_TYPE_BINARY_INPUTS ;
      printf("using binary inputs\n") ;
      break ;
    case 'M':
      for (layer = 0 ; layer < MAX_LAYERS ; layer++)
      {
	if (isXdigit(*argv[2+layer]))
	{
	  parms.momentum[layer] = atof(argv[2+layer]) ;
	  printf("using layer %d momentum %f\n", layer, parms.momentum[layer]) ;
	  nargs++ ;
	}
      else
	break ;
      }
      break ;
    case 'X':
      extract = atof(argv[2]) ;
      nargs = 1 ; 
      printf("setting extract to %d\n", extract) ;
      break ;
    case 'G':
      ngroups = atoi(argv[2]) ;
      nargs = 1 ;
      printf("using %d groups in each hidden layer\n", ngroups) ;
      break ;
    case 'N':
      parms.nsteps = atoi(argv[2]) ;
      nargs = 1 ;
      printf("using %d training iterations\n", parms.nsteps) ;
      break ;
    case 'V':
      Gdiag_no = atof(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'K':
      parms.ksize = atoi(argv[2]) ;
      nargs = 1 ;
      printf("using kernel size %d in each hidden layer\n", parms.ksize) ;
      break ;
    case 'E':
      for (layer = 0 ; layer < MAX_LAYERS ; layer++)
      {
	if (isXdigit(*argv[2+layer]))
	{
	  parms.training_rates[layer] = atof(argv[2+layer]) ;
	  printf("using layer %d training rate %f\n", layer, parms.training_rates[layer]) ;
	  nargs++ ;
	}
      else
	break ;
      }
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
  printf("usage: %s [options] <intensity volume> <label volume> <output RBM>\n",
         Progname) ;
  exit(code) ;
}



#define NINPUTS 2
#define NTRAINING 10000
static int
train_xor(RBM_PARMS *parms) 
{
  RBM    *rbm ;
  MRI    *mri_inputs, *mri_labels ;
  int    x, val1, val2, n ;
  double r, visible[NINPUTS], rms ;
  VOXLIST *vl ;
  FILE  *fp ;

  fp = fopen("rbm.log", "w") ;
  parms->held_out = 1000 ;
  parms->mini_batch_size = 1 ;
  parms->batches_per_step = 1 ;
  parms->nsteps = 1000 ;
  
  mri_inputs = MRIallocSequence(NTRAINING, 1, 1, MRI_FLOAT,NINPUTS) ;
  mri_labels = MRIallocSequence(NTRAINING, 1, 1, MRI_UCHAR,1) ;
  for (x = 0 ; x < NTRAINING ; x++)
  {
    r = randomNumber(0, 1.0) ;
    if (r < .5)
      val1 = 0 ;
    else
      val1 = 1 ;
    MRIsetVoxVal(mri_inputs, x, 0, 0, 0, val1) ;
    if (NINPUTS == 2)
    {
      MRIsetVoxVal(mri_inputs, x, 0, 0, 1, !val1) ;
    }
    else
    {
      r = randomNumber(0, 1.0) ;
      if (r < .5)
	val2 = 0 ;
      else
	val2 = 1 ;
      MRIsetVoxVal(mri_inputs, x, 0, 0, 1,  val2) ;
      MRIsetVoxVal(mri_inputs, x, 0, 0, 2,  val1 ^ val2) ;
    }

//    MRIsetVoxVal(mri_inputs, x, 0, 0, 0, 0) ;
//    MRIsetVoxVal(mri_inputs, x, 0, 0, 1,  .5) ;
//    MRIsetVoxVal(mri_inputs, x, 0, 0, 2,  1) ;
  
//  MRIsetVoxVal(mri_inputs, x, 0, 0, 0, 0) ;
//    MRIsetVoxVal(mri_inputs, x, 0, 0, 1,  1) ;

    MRIsetVoxVal(mri_labels, x, 0, 0, 0, 1) ;
  }
  vl = VLSTcreate(mri_labels, 1, 255, NULL, 0, 0) ; 
  parms->weight_decays[0] = 0 ;
  vl->mri = mri_inputs ;
  rbm = RBMalloc(RBM_TYPE_BINARY_INPUTS, NINPUTS, ngroups, 0, RBM_INPUT_VALUE) ;

#if 0
  rbm->weights[0][0] = -3.323 ;
  rbm->weights[0][1] = 3.0157 ;
  rbm->weights[1][0] = 5.41 ;
  rbm->weights[1][1] = -5.557 ;

  rbm->visible_bias[0] = 0.3 ;
  rbm->visible_bias[1] = -.5 ;
  rbm->hidden_bias[0] = -.07 ;
  rbm->hidden_bias[1] = .326 ;
#else
  RBMtrainFromVoxlistImage(rbm, vl, parms) ;
  {
    char fname[STRLEN] ;
    MRI  *mri_tmp ;

    mri_tmp = RBMreconstruct(rbm, mri_inputs, NULL, NULL, parms) ;
    sprintf(fname, "%s.recon.mgz", parms->base_name) ;
    printf("writing reconstructed image to %s\n", fname) ;
    MRIwrite(mri_tmp, fname) ;
  }
  RBMwriteNetwork(rbm, NTRAINING+1, parms, -1) ;
#endif
  for (rms = 0.0, x = 0 ; x < NTRAINING ; x++)
  {
    visible[0] = randomNumber(0.0, 1.0) ;
    if (visible[0] < .5)
      visible[0] = 0 ;
    else
      visible[0] = 1 ;
    if (NINPUTS == 2)
      visible[1] = !nint(visible[0]) ;
    else
    {
      visible[1] = randomNumber(0.0, 1.0) ;
      if (visible[1] < .5)
	visible[1] = 0 ;
      else
	visible[1] = 1 ;
      visible[2] = nint(visible[0]) ^ nint(visible[1]) ;
    }
    RBMactivateForward(rbm, visible) ;
    RBMactivateBackward(rbm) ;
//    RBMprintNetworkActivations(rbm, stdout, 0, parms) ;
    for (n = 0 ; n < 0 ; n++)
    {
      RBMactivateForward(rbm, NULL) ;
      RBMactivateBackward(rbm) ;
//      RBMprintNetworkActivations(rbm, stdout, n+1, parms) ;
    }
    for (n = 0 ; n < rbm->nvisible ; n++)
      fprintf(fp, "%2.0f ", visible[n]) ;
    for (n = 0 ; n < rbm->nvisible ; n++)
      fprintf(fp, "%d ", nint(rbm->visible[n])) ;
    for (n = 0 ; n < rbm->nvisible ; n++)
      rms += SQR(visible[n] - rbm->visible[n]) ;

    fprintf(fp, "\n") ;
    fflush(fp) ;
    DiagBreak() ;
  }

  rms = sqrt(rms / (rbm->nvisible*NTRAINING)) ;
  printf("rms = %f\n", rms) ;
  fclose(fp) ;
  return(NO_ERROR) ;
}

