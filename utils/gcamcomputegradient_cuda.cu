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
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.10 $
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

#include "gcamorph.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"
#include "gcamorphtermgpu.hpp"

#include "gcamcomputegradient_cuda.hpp"


// ========================================================================




template<typename T, typename U>
int ComputeGradient( GPU::Classes::GCAmorphGPU& gcam,
                     const GPU::Classes::MRIframeGPU<T>& mri,
                     const GPU::Classes::MRIframeGPU<U>& mri_smooth,
                     GCA_MORPH_PARMS *parms )
{

  static GPU::Algorithms::GCAmorphTerm myTerms;
  int invalid;

  // Start up
  gcam.ClearGradient();
  gcam.ComputeMetricProperties( invalid );

  SetGinvalid( invalid );

  // Map term unimplemented
  if( !DZERO(parms->l_map) )
  {
    std::cerr << __FUNCTION__
              << ": MapTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Do label term
  myTerms.LabelTerm( gcam, mri, parms->l_label, parms->label_dist );

  // Area intensity term unimplemented
  if( !DZERO(parms->l_area_intensity) )
  {
    std::cerr << __FUNCTION__
              << ": AreaIntensity not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Binary term unimplemented
  if( !DZERO(parms->l_binary) )
  {
    std::cerr << __FUNCTION__
              << ": BinaryTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Expansion term unimplemented
  if( !DZERO(parms->l_expansion) )
  {
    std::cerr << __FUNCTION__
              << ": ExpansionTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }


  // Likelihood term unimplemented
  if( !DZERO(parms->l_likelihood) )
  {
    std::cerr << __FUNCTION__
              << ": LikelihoodTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Distance transform term unimplemented
  if( !DZERO(parms->l_dtrans) )
  {
    std::cerr << __FUNCTION__
              << ": DistanceTransformTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Do the log likelihood term
  myTerms.LogLikelihood( gcam, mri, mri_smooth, parms->l_log_likelihood );


  // Multiscale term unimplemented
  if( !DZERO(parms->l_multiscale) )
  {
    std::cerr << __FUNCTION__
              << ": MultiscaleTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }


  // Distance term unimplemented
  if( !DZERO(parms->l_distance) )
  {
    std::cerr << __FUNCTION__
              << ": DistanceTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }

  // AreaSmoothness term unimplemented
  if( !DZERO(parms->l_area_smoothness) )
  {
    std::cerr << __FUNCTION__
              << ": AreaSmoothnessTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Area term unimplemented
  if( !DZERO(parms->l_area) )
  {
    std::cerr << __FUNCTION__
              << ": AreaTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Do the smoothness term
  myTerms.Smoothness( gcam, parms->l_smoothness );

  // LSmoothness term unimplemented
  if( !DZERO(parms->l_lsmoothness) )
  {
    std::cerr << __FUNCTION__
              << ": LSmoothnessTerm not implemented on GPU" << std::endl;
    exit( EXIT_FAILURE );
  }


  // Spring term unimplemented
  if( !DZERO(parms->l_spring) )
  {
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






template<typename T, typename U>
int
gcamCGFinalDispatch( GCA_MORPH *gcam,
                     const MRI *mri,
                     const MRI *mri_smooth,
                     GCA_MORPH_PARMS *parms )
{

  GPU::Classes::GCAmorphGPU myGCAM;
  GPU::Classes::MRIframeGPU<T> myMRI;
  GPU::Classes::MRIframeGPU<U> myMRIsmooth;

  // Handle the MRIs
  myMRI.Allocate( mri );
  myMRI.Send( mri, 0 );

  myMRIsmooth.Allocate( mri_smooth );
  myMRIsmooth.Send( mri_smooth, 0 );

  // Put the GCAM on the GPU
  myGCAM.CheckIntegrity(); // Shouldn't be necessary....
  myGCAM.SendAll( gcam );

  // Run the computation
  ComputeGradient( myGCAM, myMRI, myMRIsmooth, parms );

  // Retrieve results
  myGCAM.RecvAll( gcam );

  return( NO_ERROR );
}



template<typename T>
int
gcamCGSmoothDispatch( GCA_MORPH *gcam,
                      const MRI *mri,
                      const MRI *mri_smooth,
                      GCA_MORPH_PARMS *parms )
{

  switch( mri_smooth->type )
  {

  case MRI_UCHAR:
    gcamCGFinalDispatch<T,unsigned char>( gcam, mri, mri_smooth, parms );
    break;

  case MRI_FLOAT:
    gcamCGFinalDispatch<T,float>( gcam, mri, mri_smooth, parms );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised type for mri_smooth "
              << mri_smooth->type << std::endl;
    exit( EXIT_FAILURE );
  }

  return( NO_ERROR );
}



int gcamComputeGradientGPU( GCA_MORPH *gcam,
                            const MRI *mri,
                            const MRI *mri_smooth,
                            GCA_MORPH_PARMS *parms )
{

  switch( mri->type )
  {

  case MRI_UCHAR:
    gcamCGSmoothDispatch<unsigned char>( gcam, mri, mri_smooth, parms );
    break;

  case MRI_FLOAT:
    gcamCGSmoothDispatch<float>( gcam, mri, mri_smooth, parms );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised type for mri "
              << mri->type << std::endl;
    exit( EXIT_FAILURE );
  }

  return( NO_ERROR );
}



#endif
