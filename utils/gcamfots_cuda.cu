/**
 * @file  gcamfots.cu
 * @brief Implementation of gcamFindOptimalTimeStep for the GPU
 *
 * Reference:
  * "Whole Brain Segmentation: Automated Labeling of Neuroanatomical
  * Structures in the Human Brain", Fischl et al.
  * (2002) Neuron, 33:341-355.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/04/01 18:37:26 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2010,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"
#include "gcamorphenergy.hpp"


template<typename T>
float FindOptimalTimestep( GPU::Classes::GCAmorphGPU& gcam,
			   GCA_MORPH_PARMS *parms,
			   const GPU::Classes::MRIframeGPU<T>& mri ) {

  float min_dt = 0;



  return( min_dt );
}





// ============================================

template<typename T>
float FindOptimalTimestepDispatch( GCA_MORPH *gcam,
				   GCA_MORPH_PARMS *parms,
				   MRI *mri ) {
  
  GPU::Classes::GCAmorphGPU myGCAM;
  myGCAM.SendAll( gcam );

  GPU::Classes::MRIframeGPU<T> myMRI;
  myMRI.Allocate( mri );
  myMRI.Send( mri, 0 );
  myMRI.AllocateArray();
  myMRI.SendArray();


  float dt = FindOptimalTimestep( myGCAM, parms, myMRI );

  myGCAM.RecvAll( gcam );

  return( dt );
}

// ============================================


float gcamFindOptimalTimestepGPU( GCA_MORPH *gcam,
				  GCA_MORPH_PARMS *parms,
				  MRI *mri ) {
  
  float dt;

  switch( mri->type ) {

  case MRI_UCHAR:
    dt = FindOptimalTimestepDispatch<unsigned char>( gcam, parms, mri );
    break;

  default:
    std::cerr << __FUNCTION__
	      << ": Unrecognised MRI type" << std::endl;
    exit( EXIT_FAILURE );
  }

  return( dt );
  
}
  
