/**
 * @file  gcamorphenergy.hpp
 * @brief Holds routines to compute GCAmorph energies on the GPU (header)
 *
 * This file holds a variety of routines which compute energies for
 * mri_ca_register
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.7 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef GCA_MORPH_ENERGY_HPP
#define GCA_MORPH_ENERGY_HPP

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"

#include "ctfactory.hpp"

namespace GPU
{
namespace Algorithms
{

//! Class to hold GCAMorph energy computations
class GCAmorphEnergy
{
public:

  //! Constructor (stream defaults to zero)
  GCAmorphEnergy( const cudaStream_t s = 0 ) : stream(s) {};

  //! Destructor
  ~GCAmorphEnergy( void ) {};

  //! Routine to print timing information to std.out
  static void ShowTimings( void );


  // --------------------------------------------------

  //! Implementation of gcamLogLikelihoodEnergy for the GPU
  template<typename T>
  float LogLikelihoodEnergy( const GPU::Classes::GCAmorphGPU& gcam,
                             const GPU::Classes::MRIframeGPU<T>& mri ) const;

  //! Dispatch wrapper for LogLikelihoodEnergy
  template<typename T>
  float LLEdispatch( const GCA_MORPH *gcam,
                     const MRI* mri ) const;

  // --------------------------------------------------

  //! Implementation of gcamComputeJacobianEnergy for the GPU
  float ComputeJacobianEnergy( const GPU::Classes::GCAmorphGPU& gcam,
                               const float mriThick ) const;


  // --------------------------------------------------

  //! Implementation of gcamLabelEnergy for the GPU
  float LabelEnergy( const GPU::Classes::GCAmorphGPU& gcam ) const;

  // --------------------------------------------------

  //! Implementation of gcamSmoothnessEnergy for the GPU
  float SmoothnessEnergy( GPU::Classes::GCAmorphGPU& gcam ) const;

  // --------------------------------------------------

  //! Implementation of gcamComputeSSE for the GPU
  template<typename T>
  float ComputeSSE( GPU::Classes::GCAmorphGPU& gcam,
                    const GPU::Classes::MRIframeGPU<T>& mri,
                    GCA_MORPH_PARMS *parms ) const;

  // --------------------------------------------------

  template<typename T>
  float ComputeRMS( GPU::Classes::GCAmorphGPU& gcam,
                    const GPU::Classes::MRIframeGPU<T>& mri,
                    GCA_MORPH_PARMS *parms ) const
  {

    float sse = this->ComputeSSE( gcam, mri, parms );

    const dim3 dims = gcam.d_rx.GetDims();
    float nVoxels = dims.x;
    nVoxels *= dims.y;
    nVoxels *= dims.z;

    float rms = sqrtf( sse/nVoxels );

    return( rms );
  }



  template<typename T>
  float RMSdispatch( GPU::Classes::GCAmorphGPU& gcam,
                     const MRI *mri,
                     GCA_MORPH_PARMS *parms ) const;


  // ######################################################
private:
  //! Stream to use for operations
  cudaStream_t stream;

  //! Timer for LLEdispatch
  static SciGPU::Utilities::Chronometer tLLEdispatch;

  //! Timer for log likelihood energy
  static SciGPU::Utilities::Chronometer tLLEtot;
  //! Timer for LLE 'good' assesments
  static SciGPU::Utilities::Chronometer tLLEgood;
  //! Timer for LLE calculation
  static SciGPU::Utilities::Chronometer tLLEcompute;

  //! Timer for Jacobian Energy
  static SciGPU::Utilities::Chronometer tJacobTot;
  //! Timer for Jacobian Energy kernel
  static SciGPU::Utilities::Chronometer tJacobCompute;

  //! Mutable for Label Energy
  static SciGPU::Utilities::Chronometer tLabelTot;
  //! Timer for Label Energy kernel
  static SciGPU::Utilities::Chronometer tLabelCompute;

  //! Timer for smoothness energy
  static SciGPU::Utilities::Chronometer tSmoothTot;
  //! Timer for subtractions in smoothness energy
  static SciGPU::Utilities::Chronometer tSmoothSubtract;
  //! Timer for smoothness computation itself
  static SciGPU::Utilities::Chronometer tSmoothCompute;

  // --------------------------------------------------------


  //! Templated texture binding wrapper
  template<typename T>
  GPU::Classes::CTfactory* BindMRI( const GPU::Classes::MRIframeGPU<T>& mri ) const;



};


}
}


#endif

#endif
