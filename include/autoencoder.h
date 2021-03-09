/**
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

extern int test ;
#ifndef _autoencoder_h
#define _autoencoder_h

#include "mri.h"
#include "matrix.h"
#include "voxlist.h"
#include "transform.h"

#ifndef SIGMOID
#define SIGMOID(x)    (1.0/(1.0+exp(-(x))))
#endif
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
  struct _SAE *sae ;
  double     noise_fraction ;            // what fraction of inputs to zero out instead of propagating forward
  int        *zero_indices ;
  double     *saved_inputs ;             // to restore for training
  double     *average_act ;              // exponentially decaying average activation: p_t = p_t-1*.999 + .001*act
  double     sparsity_target ;           // [0,1] - average desired activation
  int        ksize ;                     // if it is a convolutional AE
  int        *hidden_class_labels ;      // if hidden nodes encode specific classes
} AUTO_ENCODER, AE ;

typedef struct _SAE
{
  MRI    *mri_inputs ;
  int    nencoders ;
  int    whalf ;
  double scale ;       // the ratio of hidden to input nodes
  int    type ;        // focused or not
  int    nlevels ;     // how many levels in the Gaussian pyramid for inputs
  AE     *first ;
  VOL_GEOM vg ;
} STACKED_AUTO_ENCODER, SAE ;

#define MAX_AE_LAYERS 50

typedef struct 
{ 
  int     nlayers ;
  AE      *aes[MAX_AE_LAYERS];            // array of layers
  SAE     *sae ;
  MRI     *mri_outputs[MAX_AE_LAYERS] ;   // outputs of the [ith] layer. Each frame is a "group"
} CONVOLUTIONAL_SAE, CSAE ;
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
  float   noise_fraction ;   // what fraction of inputs to zero out
  double  sparsity_trate ;
  double  sparsity_target ;
  double  weight_decays[MAX_AE_LAYERS] ;
  int     max_iter ;
  int     layer ;   // current layer being trained
  int     batches_per_step ;
  int     held_out ;  // how much data is held out for testing
  int     mini_batch_size ;
  char    base_name[STRLEN] ;
  int     write_iterations ;
  int     max_no_progress ;
  double  sparsity_decay ;  // l in active = l * previous_active_pval + (1-l) *  current_active_pval
  int     class_label ;     // true class label for use in training of nonzero
  double  class_weight ;
} SAE_INTEGRATION_PARMS ;

#define INTEGRATE_GRADIENT_DESCENT              0
#define INTEGRATE_BOLTZMANN_MACHINE             1
#define INTEGRATE_ANNEALING                     2
#define INTEGRATE_CONJUGATE_GRADIENT            3
#define INTEGRATE_STOCHASTIC_GRADIENT_DESCENT   4

#define NORMAL_AUTOENCODER    0x0001
#define FOCUSED_AUTOENCODER   0x0002
#define AUTOENCODER_2D        0x0008   // not 3D


SAE      *SAEalloc(int whalf, int nlevels, int type, double scale) ;
void     SAEfree(SAE **psae) ;
AE       *SAEaddLayer(SAE *sae, float scale) ;
SAE      *SAEtrainLayer(SAE *sae, AE *layer, MRI **mri, double tol) ;
double   SAEtrainFromVoxlist(SAE *sae, VOXEL_LIST *vl, MRI **mri_pyramid, SAE_INTEGRATION_PARMS *parms) ;
VECTOR   *SAEactivateLastHiddenLayer(SAE *sae, MRI *mri) ;
double   SAEtrainFromMRI(SAE *sae, MRI **mri, SAE_INTEGRATION_PARMS *parms) ;
MRI      *SAEvectorToMRI(VECTOR *v_input, int nlevels, int whalf, MRI *mri)  ;
VECTOR   *SAEactivateNetwork(SAE *sae) ;
double   SAEcomputeRMS(SAE *sae) ;
double   SAEcomputeTotalRMS(SAE *sae, MRI **mri) ;
int      SAEwrite(SAE *sae, char *fname) ;
SAE      *SAEread(char *fname) ;
MRI      *SAEinputWeightsToMRI(SAE *sae, MRI *mri)  ;
MRI      *CSAElayerWeightsToMRI(CSAE *csae, int layer);
MRI      *SAElayerWeightsToMRI(SAE *sae, int layer);
void     SAEdump(SAE *sae) ;
void     AEdump(AE *ae) ;
VECTOR   *SAEfillInputVector(MRI **mri, int nlevels, int x, int y, int z, int whalf, VECTOR *v_input) ;
AE       *SAEfindLastLayer(SAE *sae, AE *ae) ;

// convoutional autoencoder code
int        CSAEfillInputs(CSAE *csae, MRI *mri_inputs, VECTOR *v_visible, int x0, int y0, int z0, int ksize) ;
CSAE       *CSAEalloc(int type, int nlayers, int *ksizes, int *ngroups, MRI *mri_inputs) ;
AE         *CSAEaddLayer(CSAE *csae, int ksize, int nhidden) ;
int        CSAEfillInputse(CSAE *csae, MRI *mri_inputs, VECTOR *v_visible, int x0, int y0, int z0, int ksize);
double     CSAEtrainLayerFromVoxlist(CSAE *csae, int layer, VOXEL_LIST *vl, MRI **mri_pyramid, SAE_INTEGRATION_PARMS *parms) ;
int        CSAEwrite(CSAE *csae, char *fname) ;
MRI        *CSAEcreateOutputs(CSAE *csae, MRI *mri_inputs, int first_layer, int last_layer)  ;
double     CSAEcomputeTotalRMS(CSAE *csae, int layer, MRI **mri) ;
double     CSAEcomputeVoxlistRMS(CSAE *csae, SAE_INTEGRATION_PARMS *parms, int layer, MRI **mri, VOXEL_LIST *vl, int *indices, int start_index, int num_indices, int *always, int *never) ; 
CSAE       *CSAEread(char *fname) ;

#endif
