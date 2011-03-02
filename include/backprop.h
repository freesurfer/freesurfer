/**
 * @file  backprop.h
 * @brief back-propagation neural net
 *
 * backprop files can store multiple networks of varying size. In order to
 * accomodate this, I mimic the tiff file structure with a header indicating
 * the # of networks, and a pointer to the 1st one. Each one then starts with
 * a pointer to the next network, or NULL for the last network in the file.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.7 $
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

#ifndef BACKPROP_H
#define BACKPROP_H

#include "machine.h"

typedef struct
{
  long32    nnets ;         /* total # of nets in file */
  long32    magic ;         /* magic # for identification */
  long32    first ;         /* file offset of first net */
}
BPFILE_HEADER ;

#define BP_APPEND         0
#define BP_WRITE          1
#define BP_REWRITE        2
#define BP_MAGIC          4242L


/*
  for adaptive learning rate: If sse increases by more than ERROR_RATIO,
  then update the step size according to:

  learning_rate = learning_rate * TRATE_DECREASE

  otherwise, if the sse decresed, the increase the step size by:

  learning_rate = learning_rate * TRATE_INCREASE.
*/


#define BP_TRATE_DECREASE  0.7f
#define BP_TRATE_INCREASE  1.05f
#define BP_ERROR_RATIO     1.04f

typedef struct
{
  int     nunits ;        /* # of neurons in layer */
  int     ninputs ;       /* # of inputs to this layer */
  float   *x ;            /* neuron activations */
  float   *w ;            /* neuron weights */
  float   *biases ;       /* neuron biases (thresholds) */
  float   *db ;           /* last bias change for momentum */
  float   *deltas ;       /* for training */
  float   *dw ;           /* last weight change for momentum */
}
LAYER ;

typedef struct
{
  int     ninputs ;       /* # of inputs */
  int     noutputs ;      /* # of outputs (size of map field) */
  int     nhidden ;       /* # of hidden layers */
  int     learn ;         /* in learn mode */
  float   momentum ;      /* momentum coefficient */
  float   trate ;         /* learning rate */
  float   trate_up ;      /* adaptive step size increase */
  float   trate_down ;    /* adaptive step size decrease */
  float   error_ratio ;   /* ratio of new/old sse for step size modification */
  int     nepochs ;       /* # of epochs of training */
  int     ntrials ;       /* # of trials of training in this epoch */
  float   old_momentum ;  /* original momentum parameter */
  int     class ;         /* class of previous recognition */
  LAYER   hidden ;        /* hidden layer */
  LAYER   output ;        /* output layer */
  float   *mean_out ;     /* mean output value */
  float   *std_out ;      /* standard deviation of output values */
  float   *errors ;       /* output errors */
  float   sse ;           /* sum squared error */
  int     user_bytes ;    /* # of bytes in user memory */
  char    *user ;         /* ad hoc memory for user to use */
}
BACKPROP ;

BACKPROP *BackpropAlloc(int ninputs, int noutputs, int nhidden, float alpha,
                        float trate, float *mean_out, float *std_out) ;
BACKPROP *BackpropCopy(BACKPROP *bp_src, BACKPROP *bp_dst) ;
int      BackpropSetParms(BACKPROP *bp, float trate_up, float trate_down,
                          float error_rate) ;
int      BackpropFileNumberOfNets(char *fname) ;
BACKPROP *BackpropRead(char *fname, int netno) ;
int      BackpropWrite(BACKPROP *backprop, char *fname, int argc, char *argv[],
                       char *comments, int mode) ;
int      BackpropFree(BACKPROP **backprop) ;
int      BackpropProcess(BACKPROP *backprop, float *I) ;
float    BackpropLearn(BACKPROP *backprop, float *inputs, float *targets) ;
float    BackpropTrainEpoch(BACKPROP *bp, int ntrials,
                            int (*io_func)(float *inputs, float *targets,
                                           int index, int ninputs,
                                           int noutputs, void *user),
                            void *user) ;
float    BackpropError(BACKPROP *bp, float *targets) ;
int      BackpropErrorReset(BACKPROP *bp) ;
int      BackpropEpochComplete(BACKPROP *bp) ;

#endif
