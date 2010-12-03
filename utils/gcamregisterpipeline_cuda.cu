/**
 * @file  gcamregisterpipeline_cuda.cu
 * @brief Implementation of GCAMregisterPipeline for the GPU
 *
 * Reference:
  * "Whole Brain Segmentation: Automated Labeling of Neuroanatomical
  * Structures in the Human Brain", Fischl et al.
  * (2002) Neuron, 33:341-355.
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/12/03 18:48:41 $
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


#ifdef GCAMORPH_ON_GPU

#include "macros.h"
#include "error.h"

#include "gcamorph.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"


#include "gcamorphenergy.hpp"
#include "gcamregisterlevel_cuda.hpp"


// ========================================================================

template<typename T, typename U>
void RegisterPipeline( GPU::Classes::GCAmorphGPU& gcam,
		       const GPU::Classes::MRIframeGPU<T>& mri,
		       const GPU::Classes::MRIframeGPU<U>& mri_smooth,
		       GCA_MORPH_PARMS *parms,
		       double *last_rms,
		       int *level_steps,
		       int i ) {
  GPU::Algorithms::GCAmorphEnergy gcamEnergy;

  *last_rms = gcamEnergy.ComputeRMS( gcam, mri, parms );
  if( i==0 ) {
    parms->start_rms = *last_rms;
  }
  *level_steps = parms->start_t;
  RegisterLevel( gcam, mri, mri_smooth, parms);
}


// ================

template<typename T, typename U>
void
gcamRPfinalDispatch( GCA_MORPH *gcam,
		     MRI *mri,
		     MRI *mri_smooth,
		     GCA_MORPH_PARMS *parms,
		     double *last_rms,
		     int *level_steps,
		     int i ) {
  GPU::Classes::GCAmorphGPU myGCAM;
  GPU::Classes::MRIframeGPU<T> myMRI;
  GPU::Classes::MRIframeGPU<U> myMRIsmooth;

  // Handle the MRIs
  myMRI.Allocate( mri );
  myMRI.Send( mri, 0 );
  myMRI.AllocateArray();
  myMRI.SendArray();

  myMRIsmooth.Allocate( mri_smooth );
  myMRIsmooth.Send( mri_smooth, 0 );
  myMRIsmooth.AllocateArray();
  myMRIsmooth.SendArray();

  // Put the GCAM on the GPU
  myGCAM.CheckIntegrity(); // Shouldn't be necessary....
  myGCAM.SendAll( gcam );

  // Run the computation
  RegisterPipeline( myGCAM, myMRI, myMRIsmooth, parms,
		    last_rms, level_steps, i);

  // Retrieve results
  myGCAM.RecvAll( gcam );
}


// -----------


template<typename T>
void
gcamRPsmoothDispatch(  GCA_MORPH *gcam,
		       MRI *mri,
		       MRI *mri_smooth,
		       GCA_MORPH_PARMS *parms,
		       double *last_rms,
		       int *level_steps,
		       int i ) {
  switch( mri_smooth->type ) {

  case MRI_UCHAR:
    gcamRPfinalDispatch<T,unsigned char>( gcam, mri, mri_smooth, parms,
					  last_rms, level_steps, i );
    break;
    default:
    std::cerr << __FUNCTION__
	      << ": Unrecognised type for mri_smooth "
	      << mri_smooth->type << std::endl;
    abort();
  }

}


// -------------------

void GCAMregisterPipelineGPU( GCA_MORPH *gcam,
			      MRI *mri,
			      MRI *mri_smooth,
			      GCA_MORPH_PARMS *parms,
			      double *last_rms,
			      int *level_steps,
			      int i ) {

  switch( mri->type ) {

  case MRI_UCHAR:
    gcamRPsmoothDispatch<unsigned char>( gcam, mri, mri_smooth, parms,
					 last_rms, level_steps, i );
    break;

  case MRI_FLOAT:
    gcamRPsmoothDispatch<float>( gcam, mri, mri_smooth, parms,
				 last_rms, level_steps, i );
    break;
    
  default:
    std::cerr << __FUNCTION__
	      << ": Unrecognised type for mri "
	      << mri->type << std::endl;
    abort();
  }

}


#endif
