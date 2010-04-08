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
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/04/08 13:35:33 $
 *    $Revision: 1.2 $
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

#include "macros.h"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"
#include "gcamorphenergy.hpp"


template<typename T>
float FindOptimalTimestep( GPU::Classes::GCAmorphGPU& gcam,
			   GCA_MORPH_PARMS *parms,
			   const GPU::Classes::MRIframeGPU<T>& mri ) {

  GPU::Algorithms::GCAmorphEnergy gcamEnergy;

  float min_dt = 0;
  float start_rms, min_rms, rms;
  float orig_dt, start_dt, max_dt;
  int bad, prev_neg;
  int suppressed = 0;

  gcam.ClearMomentum();

  start_rms = min_rms = gcamEnergy.ComputeRMS( gcam, mri, parms );

  // Find right order of magnitude for time step
  orig_dt = parms->dt;
  if( DZERO( orig_dt ) ) {
    orig_dt = 1e-6;
    std::cerr << "orig_dt = 0 in FindOptimalTimestep" << std::endl;
  }

  // Define a pretty broad search range initially
  start_dt = (sqrt(parms->navgs)+1)*orig_dt / (10*16.0*16.0);
  max_dt =   (sqrt(parms->navgs)+1)*orig_dt * (10*16.0*16.0);


  int i = 0;
  do {

    if( i++ > 0 ) {
      start_dt *= 0.01;
    }

    if( DEQUAL( start_dt, 0 ) ) {
      break;
    }

    bad = 0;
    prev_neg = 0;

    for( parms->dt = start_dt; parms->dt <= max_dt; parms->dt *= 4 ) {
      
      gcam.ApplyGradient( parms );
      rms = gcamEnergy.ComputeRMS( gcam, mri, parms );
      gcam.UndoGradient();

      if( (gcam.neg > 0) && 0 ) {
	std::cerr << __FUNCTION__
		  << ": How did you get here?" << std::endl;
	std::cerr << "Exiting on line " << __LINE__ << std::endl;
	exit( EXIT_FAILURE );
      } else {
	prev_neg = suppressed = 0;
      }


      if( (gcam.neg>0) && DEQUAL(parms->dt, start_dt) ) {
        start_dt *= 0.1;
        parms->dt /= 4;
        continue ;
      }

      if( (gcam.neg>0) && parms->noneg == True ) {
	break;
      }

      if( rms < min_rms ) {
	min_rms = rms;
	min_dt = parms->dt;
      }

      if( DEQUAL( min_dt, max_dt ) ) {
	max_dt *= 10;
      }
    }

    if( i > 10 ) {
      break;
    }

    max_dt = start_dt*10 ; /* will only iterate if min was at start_dt */
    

  } while( DEQUAL( min_dt, start_dt ) );


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
  
