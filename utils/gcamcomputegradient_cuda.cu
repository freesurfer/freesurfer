/**
 * @file  gcamcomputegradient_cuda.cu
 * @brief Implementation of gcamComputeGradient for the GPU
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
 *    $Date: 2010/10/29 19:15:08 $
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


#ifdef GCAMORPH_ON_GPU

#include "macros.h"

#include "gcamorph.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"
#include "gcamorphtermgpu.hpp"


// ========================================================================




template<typename T, typename U>
int ComputeGradient( GPU::Classes::GCAmorphGPU& gcam,
		     GPU::Classes::MRIframeGPU<T>& mri,
		     GPU::Classes::MRIframeGPU<U>& mri_smooth,
		     GCA_MORPH_PARMS *parms ) {

  static GPU::Algorithms::GCAmorphTerm myTerms;
  int invalid;

  // Start up
  gcam.ClearGradient();
  gcam.ComputeMetricProperties( invalid );

  SetGinvalid( invalid );

  // Map term unimplemented
  if( !DZERO(parms->l_map) ) {
    std::cerr << __FUNCTION__
	      << ": MapTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }
  
  // Do label term
  myTerms.LabelTerm( gcam, mri, parms->l_label, parms->label_dist );

  // Area intensity term unimplemented
  if( !DZERO(parms->l_area_intensity) ) {
    std::cerr << __FUNCTION__
	      << ": AreaIntensity not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }
  
  // Binary term unimplemented
  if( !DZERO(parms->l_binary) ) {
    std::cerr << __FUNCTION__
	      << ": BinaryTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }
  
  // Expansion term unimplemented
  if( !DZERO(parms->l_expansion) ) {
    std::cerr << __FUNCTION__
	      << ": ExpansionTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }
  

  // Likelihood term unimplemented
  if( !DZERO(parms->l_likelihood) ) {
    std::cerr << __FUNCTION__
	      << ": LikelihoodTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }
  
  // Distance transform term unimplemented
  if( !DZERO(parms->l_dtrans) ) {
    std::cerr << __FUNCTION__
	      << ": DistanceTransformTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }
  
  // Do the log likelihood term
  myTerms.LogLikelihood( gcam, mri, mri_smooth, parms->l_log_likelihood );

  
  // Multiscale term unimplemented
  if( !DZERO(parms->l_multiscale) ) {
    std::cerr << __FUNCTION__
	      << ": MultiscaleTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }


  // Distance term unimplemented
  if( !DZERO(parms->l_distance) ) {
    std::cerr << __FUNCTION__
	      << ": DistanceTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }

  // AreaSmoothness term unimplemented
  if( !DZERO(parms->l_area_smoothness) ) {
    std::cerr << __FUNCTION__
	      << ": AreaSmoothnessTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }
  
  // Area term unimplemented
  if( !DZERO(parms->l_area) ) {
    std::cerr << __FUNCTION__
	      << ": AreaTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Do the smoothness term
  myTerms.Smoothness( gcam, parms->l_smoothness );

  // LSmoothness term unimplemented
  if( !DZERO(parms->l_lsmoothness) ) {
    std::cerr << __FUNCTION__
	      << ": LSmoothnessTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }
  

  // Spring term unimplemented
  if( !DZERO(parms->l_spring) ) {
    std::cerr << __FUNCTION__
	      << ": SpringTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }
  

  // Do the Jacobian term
  myTerms.Jacobian( gcam, parms->l_jacobian, parms->ratio_thresh );
  // gcamLimitGradientMagnitude is null operation

  //gcamSmoothGradient(gcam, parms->navgs) ;
  gcam.SmoothGradient( parms->navgs );

  return( NO_ERROR );
}
















#endif
