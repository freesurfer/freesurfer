/**
 * @file  rfa.h
 * @brief utilities for whole-brain segmentation with Random Forests
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2012/04/11 00:56:39 $
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


#ifndef RFA_H
#define RFA_H

#if defined(__cplusplus)
extern "C" {
#endif


#include "mri.h"
#include "transform.h"
#include "colortab.h"

#include "affine.h"
#include "rforest.h"


#define MAX_MRI_INPUTS  1000

typedef struct
{
  RANDOM_FOREST  *rf ;
  float          *label_priors ;
  unsigned short *labels ;
  short          nlabels;
  short          max_labels ;
} RF_NODE ;

typedef struct
{
  RF_NODE   ***nodes ;
  int       ninputs ;

  // direction cosine info
  VOL_GEOM      vg ;
  MRI          *mri_node__;       // these three MRI are helper to get
  MRI          *mri_tal__;
  MATRIX       *node_i_to_r__;
  MATRIX       *node_r_to_i__;
  MATRIX       *tal_i_to_r__;
  MATRIX       *tal_r_to_i__;
  MATRIX       *tmp__;
  int          total_training ;
  COLOR_TABLE  *ct ;
  int          wsize ;
  int          nfeatures ;   // the size of the input features
  int          width ;
  int          height ;
  int          depth ;
  float        spacing ;
} RF_CLASSIFIER_ARRAY, RFA ;


typedef struct
{
  int        max_depth ;
  int        ntrees ;
  int        nvols ;    // how many different types of input volumes are there
  float      spacing ;  // how densely to space the nodes
  int        width ;
  int        height ;
  int        depth ;
  int        wsize ;        // window size for inputs
  int        use_gradient ; // use spatial gradients in RF inputs
  int        training_size ;  // how often to train new trees 
  int        training_index ;
  double     training_fraction ;
  double     feature_fraction ;
  MRI        *mri_inputs[MAX_MRI_INPUTS] ;
  MRI        *mri_segs[MAX_MRI_INPUTS] ;
  TRANSFORM  *transforms[MAX_MRI_INPUTS] ;
} RFA_PARMS ;

  RFA  *RFAalloc(RFA_PARMS *parms, int alloc_trees) ;
int  RFAtrain(RFA *rfa, MRI *mri_inputs, MRI *mri_labels, TRANSFORM *transform) ;
int  RFAcompleteTraining(RFA *rfa, RFA_PARMS *parms) ;
int  RFAupdateTraining(RFA *rfa, RFA_PARMS *parms) ;
int  RFAaddInput(RFA *rfa, MRI *mri_seg, MRI *mri_inputs, TRANSFORM *transform, RFA_PARMS *parms) ;
int  RFAwrite(RFA *rfa, char *fname) ;
RFA  *RFAread(char *fname) ;
MRI  *RFAlabel(MRI *mri_inputs, RFA *rfa, MRI *mri_labeled, TRANSFORM *transform) ;
int  RFAsourceVoxelToAtlas( const RFA *rfa, MRI *mri, TRANSFORM *transform,
			    int xv, int yv, int zv,
			    double *px, double *py, double *pz ) ;

int  extract_feature(MRI *mri_in, int wsize, int x, int y, int z, double *feature, int xatlas, int yatlas, int zatlas) ;

  
#if defined(__cplusplus)
};
#endif



#endif
