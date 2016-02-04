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
 *    $Author: zkaufman $
 *    $Date: 2016/02/04 20:23:05 $
 *    $Revision: 1.5 $
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
float RegisterPipeline( GPU::Classes::GCAmorphGPU& gcam,
			const GPU::Classes::MRIframeGPU<T>& mri,
			const GPU::Classes::MRIframeGPU<U>& mri_smooth,
			GCA_MORPH_PARMS *parms,
			double *last_rms,
			int *level_steps,
			int i )
{
  GPU::Algorithms::GCAmorphEnergy gcamEnergy;

  *last_rms = gcamEnergy.ComputeRMS( gcam, mri, parms );
  if( i==0 )
  {
    parms->start_rms = *last_rms;
  }
  *level_steps = parms->start_t;
  RegisterLevel( gcam, mri, mri_smooth, parms);

  return gcamEnergy.ComputeRMS( gcam, mri, parms );
}


// ================

template<typename T, typename U>
float
gcamRPfinalDispatch( GCA_MORPH *gcam,
                     MRI *mri,
                     MRI *mri_smooth,
                     GCA_MORPH_PARMS *parms,
                     double *last_rms,
                     int *level_steps,
                     int i )
{
  GPU::Classes::GCAmorphGPU myGCAM;
  GPU::Classes::MRIframeGPU<T> myMRI;
  GPU::Classes::MRIframeGPU<U> myMRIsmooth;
  float result;

  // Handle the MRIs
  myMRI.Allocate( mri );
  myMRI.Send( mri, 0 );

  myMRIsmooth.Allocate( mri_smooth );
  myMRIsmooth.Send( mri_smooth, 0 );

  // Put the GCAM on the GPU
  myGCAM.CheckIntegrity(); // Shouldn't be necessary....
  myGCAM.SendAll( gcam );

  // Run the computation
  result = RegisterPipeline( myGCAM, myMRI, myMRIsmooth, parms,
			     last_rms, level_steps, i);

  // Retrieve results
  myGCAM.RecvAll( gcam );

  return result;
}


// -----------


template<typename T>
float
gcamRPsmoothDispatch(  GCA_MORPH *gcam,
                       MRI *mri,
                       MRI *mri_smooth,
                       GCA_MORPH_PARMS *parms,
                       double *last_rms,
                       int *level_steps,
                       int i )
{
  float result;

  switch( mri_smooth->type )
  {

  case MRI_UCHAR:
    result = gcamRPfinalDispatch<T,unsigned char>( gcam, mri, mri_smooth, parms,
						   last_rms, level_steps, i );
    break;
  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised type for mri_smooth "
              << mri_smooth->type << std::endl;
    abort();
  }

  return result;
}


// -------------------

float GCAMregisterPipelineAndComputeRMSGPU( GCA_MORPH *gcam,
					    MRI *mri,
					    MRI *mri_smooth,
					    GCA_MORPH_PARMS *parms,
					    double *last_rms,
					    int *level_steps,
					    int i )
{
  float result;

  switch( mri->type )
  {

  case MRI_UCHAR:
    result = gcamRPsmoothDispatch<unsigned char>( gcam, mri, mri_smooth, parms,
						  last_rms, level_steps, i );
    break;

  case MRI_FLOAT:
    result = gcamRPsmoothDispatch<float>( gcam, mri, mri_smooth, parms,
					  last_rms, level_steps, i );
    break;
    
  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised type for mri "
              << mri->type << std::endl;
    abort();
  }

  return result;
}


#endif
