/**
 * @brief utilities for restricted Boltzmann Machines (RBMs) and Deep Belief
 *        Networks (DBNs)
 *
 * Reference:
 * Hinton, G. E., Osindero, S., and Teh, Y. W. (2006a). 
 * A fast learning algorithm for deep belief nets.
 * Neural Computation, 18(7):1527{1554.
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

#ifndef RBM_H
#define RBM_H

#include "mri.h"
#include "voxlist.h"

#define MAX_DBN_LAYERS 100
#define MAX_RBM_LABELS 100

// hidden layer is always binary
#define RBM_TYPE_BINARY_INPUTS     1
#define RBM_TYPE_CONTINUOUS_INPUTS 2
#define RBM_INPUT_IMAGE            1
#define RBM_INPUT_VALUE            2

#ifndef SIGMOID  // also defined in rbm.h
#define SIGMOID(x)        ((1) / (1 + exp(-(x))))
#endif
#define  GAUSSIAN(x,u,s)  ((1.0/(s*sqrt(2*M_PI)))*exp(-(((x)-u)*((x)-u))/(2.0*s*s)))

#define MAX_LAYERS 100

typedef struct
{
  int    label ;
  double pval ;
} LABEL_PVAL ;

typedef struct
{
  int    type ;
  int    input_type ;
  int    nvisible ;
  int    nhidden ;         
  int    nlabels ;
  int    ksize ;
  MRI    *mri_inputs ;    // for diagnostics

  // nlabels of these
  double *labels ;           // p-values of labels
  double *lact ;             // total input to this label node
  double *label_bias ;
  int    *label_states ;     // 1-of-C vector of class states
  int    *sorted_labels; 
  LABEL_PVAL *label_pvals ;  // list of sorted indices of labels (1st has highest pval, etc...)
  double **label_weights ;

  // nhidden of these
  double *hidden_bias ;
  double *act ;           // weighted summed inputs
  double *hidden ;        // output of sigmoid of summed inputs
  double *hidden_state ;  // binarized version of the previous
  double *active_pvals ;  // probability that each hidden node is active
  double *active ;        // the % of time this node was active

  // nvisible of these
  double *visible ;
  double *visible_bias ;
  double *variance ;
  double **weights ;         // nvisible x nhidden
} RBM, RESTRICTED_BOLTZMANN_MACHINE ;

typedef struct
{
  double  weight_decays[MAX_LAYERS] ;
  double  training_rates[MAX_LAYERS] ;
  double  momentum[MAX_LAYERS] ;
  double  sparsity[MAX_LAYERS] ;        // what fraction of the time we want hidden nodes active
  double  l_sparsity[MAX_LAYERS] ;      // weight of sparsity term
#if 0
  double  **dU_last ;
  double  *dlabel_bias ;
  double  **dw_last ;
  double  *dhidden_bias_last ;
  double  *dvisible_bias_last ;
#endif
  int     write_iterations ;
  char    base_name[STRLEN] ;
  int     mini_batch_size ;
  int     batches_per_step ;
  int     ksize ;
  int     nsteps ;
  int     held_out ;  // how much data is held out for testing
  int     debug ;
  int     max_no_progress ;
  int     learn_variance ;  // 1 --> learn the variance, 0 --> it is fixed
  int     nclasses ;        // divide the estimated variance by this
  double  variance ;        // user-specified variance
  int     nlayers ;
  double  sparsity_decay ;  // l in active = l * previous_active_pval + (1-l) *  current_active_pval
  int     discriminative_training ;  // do discriminative training after normal CD is done
  double  l_label ;          // weight of label term
  double  label_trate ;     // training_rate of label term
} RBM_PARMS ;


typedef struct 
{ 
  int     nlayers ;
  RBM     **rbms ;
} DEEP_BELIEF_NETWORK, DBN ;

typedef struct 
{ 
  int     nlayers ;
  RBM     **rbms ;
  MRI     **mri_outputs ;   // outputs of the [ith] layer. Each frame is a "group"
} CONVOLUTIONAL_DEEP_BELIEF_NETWORK, CDBN ;


RBM     *RBMcopy(RBM *rbm_src, RBM *rbm_dst);
RBM     *RBMalloc(int type, int nvisible, int nhidden, int nlabels, int input_type) ;
int     RBMfree(RBM **prbm) ;
RBM     *RBMread(char *fname) ; 
int     RBMwrite(RBM *rbm, char *fname) ; 
int     RBMactivateBackward(RBM *rbm) ;
int     RBMactivateForward(RBM *rbm, double *visible) ;
int     RBMtrainFromImage(RBM *rbm, MRI *mri_inputs, MRI *mri_labels, RBM_PARMS *parms);
int     RBMsampleHidden(RBM *rbm) ;
double  RBMvoxlistRMS(RBM *rbm, VOXLIST *vl, RBM_PARMS *parms, int *indices, int index, int num) ;
double  RBMvoxlistHistoRMS(RBM *rbm, VOXLIST *vl, RBM_PARMS *parms, int *indices, int index, int num) ;
int     RBMcomputeGradients(RBM *rbm, VOXLIST *vl, double **dw, 
			    double *dvisible_bias, double *dhidden_bias, double *dvariance,
			    double *dlabel_bias, double **dlabel_weights,
			    RBM_PARMS *parms, int *indices, int index)  ;
MRI     *RBMreconstruct(RBM *rbm, MRI *mri_inputs, MRI *mri_reconstructed, MRI **pmri_labeled, RBM_PARMS *parms) ;
int     RBMtrainFromVoxlistImage(RBM *rbm, VOXLIST *vl, RBM_PARMS *parms)  ;
int     RBMtrainFromVoxlist(RBM *rbm, VOXLIST *vl, RBM_PARMS *parms)  ;
int     RBMwriteNetwork(RBM *rbm, int n, RBM_PARMS *parms, int layer)  ;
int     RBMfillVisible(RBM *rbm, MRI *mri_inputs, double *visible, int x0, int y0, int z0, int f0, int ksize) ;
int     RBMprintNetworkActivations(RBM *rbm, FILE *fp, int n, RBM_PARMS *parms)  ;
int     RBMcountHiddenActive(RBM *rbm) ;
int     RBMsetLabel(RBM *rbm, int label) ;
int     RBMmostLikelyLabel(RBM *rbm) ;
double  RBMfreeEnergy(RBM *rbm, double *visible) ;

// Deep belief networks
int     DBNtrainFromImage(DBN *dbn, MRI *mri_inputs, MRI *mri_labels, RBM_PARMS *parms);
DBN     *DBNalloc(int type, int nlayers, int nvisible, int *nhidden, int nlabels, int input_type)  ;
int     DBNfree(DBN **pdbn)  ;
int     DBNtrainFromVoxlistImage(DBN *dbn, VOXLIST *vl, RBM_PARMS *parms)  ;
int     DBNwriteNetwork(DBN *dbn, int n, RBM_PARMS *parms)  ;
MRI     *DBNreconstruct(DBN *dbn, MRI *mri_inputs, MRI *mri_reconstructed, MRI **pmri_labeled, RBM_PARMS *parms) ;
int     DBNcomputeGradients(DBN *dbn, int layer, VOXLIST *vl, double **dw, 
			    double *dvisible_bias, double *dhidden_bias, double *dvariance,
			    double *dlabel_bias, double **dlabel_weights,
			    RBM_PARMS *parms, int *indices, int index)  ;
double  DBNvoxlistRMS(DBN *dbn, int nlayers, VOXLIST *vl, RBM_PARMS *parms, int *indices, int index, int num);
DBN     *DBNread(char *fname)  ;
int     DBNwrite(DBN *dbn, char *fname) ;
int     DBNactivateForward(DBN *dbn, double *visible, int nlayers) ;
int     DBNactivateBackward(DBN *dbn, int last_layer, int first_layer) ;

// convolutional deep belief networks
int     CDBNtrainFromImage(CDBN *cdbn, MRI *mri_inputs, MRI *mri_labels, RBM_PARMS *parms);
CDBN    *CDBNalloc(int type, int nlayers, int *ksizes, int *ngroups, int nlabels, MRI *mri_inputs) ;
int     CDBNtrainFromVoxlistImage(CDBN *cdbn, VOXLIST *vl, RBM_PARMS *parms, MRI *mri_inputs, MRI *mri_labels)  ;
double  CDBNvoxlistRMS(CDBN *cdbn, int nlayers, VOXLIST *vl, RBM_PARMS *parms, int *indices, int index, int num, double *plabel_rms);
int     CDBNcomputeDiscriminativeGradients(CDBN *cdbn, int layer, VOXLIST *vl, double **dw, 
					   double *dhidden_bias, double *dvariance,
					   double *dlabel_bias, double **dlabel_weights,
					   RBM_PARMS *parms, int *indices, int index) ;
int     CDBNcomputeLabelGradients(CDBN *cdbn, int layer, VOXLIST *vl, double **dw, 
				  double *dhidden_bias, double *dvariance,
				  double *dlabel_bias, double **dlabel_weights,
				  RBM_PARMS *parms, int *indices, int index) ;
int     CDBNcomputeGradients(CDBN *cdbn, int layer, VOXLIST *vl, double **dw, 
			    double *dvisible_bias, double *dhidden_bias, double *dvariance,
			    double *dlabel_bias, double **dlabel_weights,
			    RBM_PARMS *parms, int *indices, int index)  ;
int     CDBNactivateBackward(CDBN *cdbn, int last_layer, int first_layer) ;
int     CDBNactivateForward(CDBN *cdbn, double *visible, int nlayers) ;
int     CDBNfillVisible(CDBN *cdbn, MRI *mri_inputs, double *visible, int x0, int y0, int z0, int ksize) ;
MRI     *CDBNcreateOutputs(CDBN *cdbn, RBM_PARMS *parms, MRI *mri_inputs, int first_layer, int last_layer, MRI **pmri_labeled) ;
int     CDBNwriteNetwork(CDBN *cdbn, int n, RBM_PARMS *parms, int layer) ;
MRI     *cdbn_layer_weights(CDBN *cdbn, int layer) ;


#if 0
int     DBNwriteNetwork(DBN *dbn, int n, RBM_PARMS *parms)  ;
MRI     *DBNreconstruct(DBN *dbn, MRI *mri_inputs, MRI *mri_reconstructed, MRI **pmri_labeled, RBM_PARMS *parms) ;
#endif

#endif
