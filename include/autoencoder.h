/**
 * @file  autoencoder.h
 * @brief header file for creating and training a stacked autoencoder for feature extraction.
 *
H.-C. Shin, M. R. Orton, D. J. Collins, S. J. Doran, and M. O. Leach,
"Stacked Autoencoders for
Unsupervised Feature Learning and Multiple Organ Detectionin a Pilot Study
Using 4D Patient Data,"
IEEE Transaction on Pattern Analysis and Machine Intelligence, 2012.

 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2013/11/03 19:56:13 $
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

extern int test ;
#ifndef _autoencoder_h
#define _autoencoder_h

#include "mri.h"
#include "matrix.h"

#define SIGMOID(x)    (1.0/(1.0+exp(-(x))))
#define D_SIGMOID(o)  ((o) * (1-(o)))

typedef struct _AE
{
  MATRIX  *m_input_to_hidden ;    // wts from input to hidden
  VECTOR  *v_output_bias ;
  VECTOR  *v_hidden_bias ;
  MATRIX  *m_hidden_to_output ;   // wts from hidden to output
  VECTOR  *v_hidden_net ;         // hidden unit 'activation' before nonlinearity
  VECTOR  *v_hidden ;             // hidden unit 'activation'
  VECTOR  *v_input ;              // input vector for this layer (unless previous is null, same as output of prev. layer
  VECTOR  *v_output ;             // output vector (same as input of next layer unless next is null)
  
  MATRIX  *m_saved_input_to_hidden ;    // wts from input to hidden
  VECTOR  *v_saved_output_bias ;
  VECTOR  *v_saved_hidden_bias ;
  MATRIX  *m_saved_hidden_to_output ;   // wts from hidden to output
  VECTOR  *v_saved_hidden_net ;         // hidden unit 'activation' before nonlinearity

  struct _AE *next ;
  struct _AE *prev ;
  VECTOR     *v_error ;                  // only needed for training
  MATRIX     *m_grad_input_to_hidden ;   // only needed for training
  MATRIX     *m_grad_hidden_to_output ;  // only needed for training
  MATRIX     *v_grad_hidden_bias ;       // only needed for training
  MATRIX     *v_grad_output_bias     ;   // only needed for training
  MATRIX     *m_previous_step_input_to_hidden ;   // only needed for training
  MATRIX     *m_previous_step_hidden_to_output ;  // only needed for training
  MATRIX     *v_previous_step_hidden_bias ;       // only needed for training
  MATRIX     *v_previous_step_output_bias     ;   // only needed for training
} AUTO_ENCODER, AE ;

typedef struct
{
  MRI    *mri_inputs ;
  int    nencoders ;
  int    whalf ;
  double scale ;  // the ratio of hidden to input nodes
  AE  *first ;
} STACKED_AUTO_ENCODER, SAE ;

/*
  double proposal_sigma = 5.0 ;   stddev of noise distribution
  double acceptance_sigma = .5 ;
*/

typedef struct
{
  double proposal_sigma ;    // stddev of noise distribution for boltzmann
  double acceptance_sigma ;
  double dt ;
  double orig_dt ;
  double tol ;
  double momentum ;
  int    integration_type ;
  char   *out_fname ;        // for diagnostics
  double temperature ;       // for annealing
  VECTOR *v_prev_grad_change_hidden_bias ;
  VECTOR *v_prev_grad_change_output_bias ;
  MATRIX *m_prev_grad_change_input_to_hidden ;
  MATRIX *m_prev_grad_change_hidden_to_output ;
  VECTOR *v_prev_grad_hidden_bias ;
  VECTOR *v_prev_grad_output_bias ;
  MATRIX *m_prev_grad_input_to_hidden ;
  MATRIX *m_prev_grad_hidden_to_output ;
  double norm_hidden_bias ;
  double norm_output_bias ;
  double norm_hidden_to_output ;
  double norm_input_to_hidden ;
  VECTOR *v_dir_hidden_bias ;
  VECTOR *v_dir_output_bias ;
  MATRIX *m_dir_input_to_hidden ;
  MATRIX *m_dir_hidden_to_output ;
} SAE_INTEGRATION_PARMS ;

#define INTEGRATE_GRADIENT_DESCENT   0
#define INTEGRATE_BOLTZMANN_MACHINE  1
#define INTEGRATE_ANNEALING          2
#define INTEGRATE_CONJUGATE_GRADIENT 3

SAE      *SAEalloc(int whalf, double scale) ;
void     SAEfree(SAE **psae) ;
AE       *SAEaddLayer(SAE *sae, float scale) ;
SAE      *SAEtrainLayer(SAE *sae, AE *layer, MRI *mri, double tol) ;
VECTOR   *SAEactivateLastHiddenLayer(SAE *sae, MRI *mri) ;
double   SAEtrainFromMRI(SAE *sae, MRI *mri, SAE_INTEGRATION_PARMS *parms) ;
MRI      *SAEvectorToMRI(VECTOR *v_input, int whalf, MRI *mri)  ;
VECTOR   *SAEactivateNetwork(SAE *sae) ;
double   SAEcomputeRMS(SAE *sae) ;
double   SAEcomputeTotalRMS(SAE *sae, MRI *mri) ;
int      SAEwrite(SAE *sae, char *fname) ;
SAE      *SAEread(char *fname) ;
MRI *    SAEinputWeightsToMRI(SAE *sae, MRI *mri)  ;
void     SAEdump(SAE *sae) ;
void     AEdump(AE *ae) ;
int      SAEfillInputVector(MRI *mri, int x, int y, int z, int whalf, VECTOR *v_input) ;

#endif
