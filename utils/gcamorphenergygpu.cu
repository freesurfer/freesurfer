/**
 * @file  gcamorphenergygpu.cu
 * @brief Holds routines to compute GCAmorph energies on the GPU
 *
 * This file holds a variety of routines which compute energies for
 * mri_ca_register
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/03/23 13:34:01 $
 *    $Revision: 1.16 $
 *
 * Copyright (C) 2002-2008,
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

#include <thrust/device_new_allocator.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>

#include "cma.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"

// Stolen from gcamorph.c
#define MIN_STD 2
#define MIN_VAR (MIN_STD*MIN_STD)


//! Texture reference for an unsigned char mri
texture<unsigned char, 3, cudaReadModeNormalizedFloat> dt_mri_uchar;


// ==============================================================

namespace GPU {
  namespace Algorithms {


    const unsigned int kGCAmorphLLEkernelSize = 16;
    const unsigned int kGCAmorphJacobEnergyKernelSize = 16;

    //! Templated texture fetch
    template<typename T>
    __device__ float FetchMRIVoxel( const float3 r ) {
      /*!
	Does look ups into the textures.
	Recall that the textures are configured to return
	normalised floats (except for the float texture!)
	in order to enable linear filtering.
	This in turn requires a conversion from the normalised
	value back to the true value.
	Unspecialised version writes junk
      */
      return( 100000 );
    }

    template<>
    __device__ float FetchMRIVoxel<unsigned char>( const float3 r ) {

      float texVal;
      texVal = tex3D( dt_mri_uchar, r.x+0.5f, r.y+0.5f, r.z+0.5f );
      texVal *= UCHAR_MAX;

      return( texVal );
    }


    __device__
    float GCAmorphDist( const float mean, const float variance,
			const float val ) {
      float v;
      v = val - mean;
      v = v*v;
      v /= variance;

      return( v );
    }

    __device__
    bool IsUnknown( const int label ) {
      bool res;

      res = (label==Unknown);
      res = res || (label==255);
      res = res || (label==Bright_Unknown);
      res = res || (label==Dark_Unknown);
      
      return(res);
    }


    //! Kernel to compute whether a voxel should be included in LLE calculation
    __global__
    void ComputeGoodLLE( const GPU::Classes::VolumeArgGPU<char> invalid,
			 const GPU::Classes::VolumeArgGPU<int> label,
			 const GPU::Classes::VolumeArgGPU<int> status,
			 GPU::Classes::VolumeArgGPU<char> good ) {

      const unsigned int bx = ( blockIdx.x * blockDim.x );
      const unsigned int by = ( blockIdx.y * blockDim.y );
      const unsigned int ix = threadIdx.x + bx;
      const unsigned int iy = threadIdx.y + by;
      /*!
	Computes whether each voxel should be included in the LLE
	calculation.
	The 'border' loops should really use shared memory, but
	right now, it's pointless worrying about this.
	Getting data onto the GPU takes two orders of magnitude
	more time
      */

      // Loop over z slices
      for( unsigned int iz = 0; iz< good.dims.z; iz++ ) {

	// Only compute if ix, iy & iz are inside the bounding box
	if( good.InVolume(ix,iy,iz) ) {
	  good(ix,iy,iz) = 0;
	  
	  // Is it valid?
	  if( invalid(ix,iy,iz) == GCAM_POSITION_INVALID ) {
	    // If not, go to next z slice
	    continue;
	  }
	  
	  // What's the status?
	  if( status(ix,iy,iz) &
	      (GCAM_IGNORE_LIKELIHOOD|GCAM_NEVER_USE_LIKELIHOOD) ) {
	    // Go to next z slice
	    continue;
	  }
	  
	  // Don't use unknowns unless they border known
	  if( IS_UNKNOWN(label(ix,iy,iz)) ) {
	    unsigned int diffLabels = 0;
	    const int myLabel = label(ix,iy,iz);

	    for( unsigned int k=max(0,iz-1);
		 k<=min(invalid.dims.z-1,iz+1);
		 k++ ) {
	      for( unsigned int j=max(0,iy-1);
		   j<=min(invalid.dims.y-1,iy+1);
		   j++ ) {
		for( unsigned int i=max(0,ix-1);
		     i<=min(invalid.dims.x-1,ix+1);
		     i++ ) {
		  if( label(i,j,k) != myLabel ) {
		    diffLabels++;
		  }
		}
	      }
	    }

	    if( diffLabels == 0 ) {
	      // Go to next z slice
	      continue;
	    }
	  }

	  // If we get to here, it's OK
	  good(ix,iy,iz) = 1;
	}
      }
    }



    template<typename T>
    __global__
    void ComputeLLE( const GPU::Classes::VolumeArgGPU<float> rx,
		     const GPU::Classes::VolumeArgGPU<float> ry,
		     const GPU::Classes::VolumeArgGPU<float> rz,
		     const GPU::Classes::VolumeArgGPU<char> good,
		     const GPU::Classes::VolumeArgGPU<float> mean,
		     const GPU::Classes::VolumeArgGPU<float> variance,
		     float* energies ) {
      
      const unsigned int bx = ( blockIdx.x * blockDim.x );
      const unsigned int by = ( blockIdx.y * blockDim.y );
      const unsigned int ix = threadIdx.x + bx;
      const unsigned int iy = threadIdx.y + by;

      float myEnergy;


      // Loop over z slices
      for( unsigned int iz = 0; iz < rx.dims.z; iz++ ) {

	// Only compute if ix, iy & iz are inside the bounding box
	if( rx.InVolume(ix,iy,iz) ) {

	  const unsigned int iLoc = rx.Index1D( ix, iy, iz );
	  
	  // See if we want to do this pixel
	  if( good(ix,iy,iz) == 0 ) {
	    energies[iLoc] = 0;
	    continue;
	  }

	  float3 r = make_float3( rx( ix, iy, iz ),
				  ry( ix, iy, iz ),
				  rz( ix, iy, iz ) );

	  // Get the MRI value
	  // For some reason, bounds checking causes lots of zeros
	  float mriVal = FetchMRIVoxel<T>( r );
	  
	  // Compute contribution to the energy
	  if( variance(ix,iy,iz) >= 0 ) {
	    // We have a valid variance
	    myEnergy = GCAmorphDist( mean(ix,iy,iz),
				     variance(ix,iy,iz),
				     mriVal );
	    myEnergy += logf( variance(ix,iy,iz) );
	  } else {
	    myEnergy = mriVal*mriVal / MIN_VAR;
	  }

	  energies[iLoc] = myEnergy;
	}
      }
    }

    // --------------------------------------------------
    
    //! Implements FZERO macro
    __device__ bool NearZero( const float f ) {

      bool res = false;
      if( fabsf(f) < 0.0000001f ) {
	res = true;
      }

      return( res );
    }


    //! Calculates delta based on exponent
    __device__ float JacobDelta( const float exponent ) {
      /*!
	Calculates the value of delta for ComputeJacobEnergy.
	This states
	\f[
	\delta = \log \left ( 1 + e^{x} \right )
	\f]
	Uses Taylor expansions to avoid hideous accuracy losses.
	Note \f$2^{12} \approx e^{8.3}\f$ is used to assess
	changeover points for single precision.
      */
      float delta;

      const float kMaxExp = 200;

      if( exponent > kMaxExp ) {
	delta = 0;
      }
      else if( exponent > 9 ) {
	delta = exponent;
      } else if ( exponent < -9 ) {
	delta = exp( exponent );
      } else {
	delta = log( 1 + exp( exponent ) );
      }

      return( delta );
    }

    //! Kernel to implement loops of gcamComputeJacobianEnergy
    __global__
    void ComputeJacobEnergy( const GPU::Classes::VolumeArgGPU<char> invalid,
			     const GPU::Classes::VolumeArgGPU<float> origArea1,
			     const GPU::Classes::VolumeArgGPU<float> origArea2,
			     const GPU::Classes::VolumeArgGPU<float> area1,
			     const GPU::Classes::VolumeArgGPU<float> area2,
			     const float exp_k,
			     const float thick,
			     float *energies ) {

      const unsigned int bx = ( blockIdx.x * blockDim.x );
      const unsigned int by = ( blockIdx.y * blockDim.y );
      const unsigned int ix = threadIdx.x + bx;
      const unsigned int iy = threadIdx.y + by;

      float ratio, exponent, delta;

      // Loop over z slices
      for( unsigned int iz = 0; iz < invalid.dims.z; iz++ ) {

	// Only compute if ix, iy & iz are inside the bounding box
	if( invalid.InVolume(ix,iy,iz) ) {

	  const unsigned int iLoc = invalid.Index1D( ix, iy, iz );

	  // Check for validity
	  if( invalid(ix,iy,iz) != GCAM_VALID ) {
	    energies[iLoc] = 0;
	    continue;
	  }

	  float myEnergy = 0;

	  // Look at area1
	  if( !NearZero( origArea1(ix,iy,iz) ) ) {
	    ratio = area1(ix,iy,iz) / origArea1(ix,iy,iz);
	    exponent = -exp_k * ratio;
	    
	    delta = JacobDelta( exponent );

	    myEnergy += ( delta * thick );
	  }

	  if( !NearZero( origArea2(ix,iy,iz) ) ) {
	    ratio = area2(ix,iy,iz) / origArea2(ix,iy,iz);
	    exponent = -exp_k * ratio;
	    
	    delta = JacobDelta( exponent );
	    
	    //energies[iLoc] = exponent;
	    myEnergy += ( delta * thick );
	  }

	  energies[iLoc] = myEnergy;
	  
	}
	  
      }
    }

    // ##############################################################

    //! Class to hold GCAMorph energy computations
    class GCAmorphEnergy {
    public:

      //! Constructor
      GCAmorphEnergy( const cudaStream_t s = 0 ) : stream(s),
						   tLLEtot(),
						   tLLEgood(),
						   tLLEcompute() {}

      //! Destructor
      ~GCAmorphEnergy( void ) {
#ifdef CUDA_SHOW_TIMINGS
	std::cout << "==================================" << std::endl;
	std::cout << "GPU GCAmorph energy timers" << std::endl;
	std::cout << "--------------------------" << std::endl;
#ifndef CUDA_FORCE_SYNC
	std::cout << "WARNING: CUDA_FORCE_SYNC not #defined" << std::endl;
	std::cout << "Timings may not be accurate" << std::endl;
#endif
	std::cout << std::endl;

	std::cout << "LLEdispatch:" << std::endl;
	std::cout << "Total         : " << this->tLLEdispatch << std::endl;
	std::cout << std::endl;

	std::cout << "LogLikelihoodEnergy:" << std::endl;
	std::cout << "    Find good : " << this->tLLEgood << std::endl;
	std::cout << "      Compute : " << this->tLLEcompute << std::endl;
	std::cout << "Total         : " << this->tLLEtot << std::endl;
	std::cout << std::endl;


	std::cout << "==================================" << std::endl;
#endif
      }
      
      // --------------------------------------------------

      //! Implementation of gcamLogLikelihoodEnergy for the GPU
      template<typename T>
      float LogLikelihoodEnergy( const GPU::Classes::GCAmorphGPU& gcam,
				 const GPU::Classes::MRIframeGPU<T>& mri ) const {
	/*!
	  This the the host side function for
	  gcamLogLikelihoodEnergy on the GPU.
	  Note that a GCAmorphGPU implicitly only has one input for
	  each location.
	  This means that each covariance is just a variance,
	  and negative values flag
	*/

	tLLEtot.Start();

	// Make sure the GCAM is sane
	gcam.CheckIntegrity();

	const dim3 gcamDims = gcam.d_rx.GetDims();
	const unsigned int nVoxels = gcamDims.x * gcamDims.y * gcamDims.z;
	// Create a 'flag' array
	GPU::Classes::VolumeGPU<char> d_good;
	d_good.Allocate( gcamDims );

	// Allocate thrust arrays
	thrust::device_ptr<float> d_energies;
	d_energies = thrust::device_new<float>( nVoxels );

	// Get the MRI into a texture
	this->BindMRI( mri );


	// Run the computation
	dim3 grid, threads;
	threads.x = threads.y = kGCAmorphLLEkernelSize;
	threads.z = 1;

	grid = gcam.d_rx.CoverBlocks( kGCAmorphLLEkernelSize );
	grid.z = 1;

	tLLEgood.Start();
	ComputeGoodLLE<<<grid,threads,0,this->stream>>>( gcam.d_invalid,
							 gcam.d_label,
							 gcam.d_status,
							 d_good );
	CUDA_CHECK_ERROR( "ComputeGood kernel failed!\n" );
	tLLEgood.Stop();


	tLLEcompute.Start();
	ComputeLLE<T><<<grid,threads,0,this->stream>>>
	  ( gcam.d_rx, gcam.d_ry, gcam.d_rz,
	    d_good,
	    gcam.d_mean, gcam.d_variance,
	    thrust::raw_pointer_cast( d_energies ) );
	CUDA_CHECK_ERROR( "ComputeLLE kernel failed!\n" );
	tLLEcompute.Stop();

	CUDA_SAFE_CALL( cudaStreamSynchronize( this->stream ) );

	// Release the MRI texture
	this->UnbindMRI<T>();

	// Get the sum of the energies
	float energy = thrust::reduce( d_energies, d_energies+nVoxels );

	// Release thrust arrays
	thrust::device_delete( d_energies );

	tLLEtot.Stop();

	return( energy );
      }


      //! Dispatch wrapper for LogLikelihoodEnergy
      template<typename T>
      float LLEdispatch( const GCA_MORPH *gcam,
			 const MRI* mri ) const {
	
	float energy;

	tLLEdispatch.Start();

	GPU::Classes::GCAmorphGPU myGCAM;
	myGCAM.SendAll( gcam );

	GPU::Classes::MRIframeGPU<T> myMRI;
	myMRI.Allocate( mri );
	myMRI.Send( mri, 0 );
	myMRI.AllocateArray();
	myMRI.SendArray( this->stream );
	
	energy = this->LogLikelihoodEnergy( myGCAM, myMRI );

	tLLEdispatch.Stop();

	return( energy );

      }

      // --------------------------------------------------

      //! Implementation of gcamComputeJacobianEnergy for the GPU
      float ComputeJacobianEnergy( const GPU::Classes::GCAmorphGPU& gcam,
				   const MRI* mri ) const {
	
	

	const float thick = ( mri ? mri->thick : 1.0 );
	
	// Make sure the GCAM is sane
	gcam.CheckIntegrity();

	const dim3 gcamDims = gcam.d_rx.GetDims();
	const unsigned int nVoxels = gcamDims.x * gcamDims.y * gcamDims.z;

	// Allocate thrust arrays
	thrust::device_ptr<float> d_energies;
	d_energies = thrust::device_new<float>( nVoxels );


	// Run the computation
	dim3 grid, threads;
	threads.x = threads.y = kGCAmorphJacobEnergyKernelSize;
	threads.z = 1;

	grid = gcam.d_rx.CoverBlocks( kGCAmorphJacobEnergyKernelSize );
	grid.z = 1;

	ComputeJacobEnergy<<<grid,threads>>>
	  ( gcam.d_invalid,
	    gcam.d_origArea1, gcam.d_origArea2,
	    gcam.d_area1, gcam.d_area2,
	    gcam.exp_k, thick,
	    thrust::raw_pointer_cast( d_energies ) );
	CUDA_CHECK_ERROR( "ComputeJacobEnergy kernel failed!\n" );
#if 0
	for( unsigned int i=0; i<nVoxels; i++ ) {
	  if( d_energies[i] != 0 )
	    {
	    std::cout << i << " "
		      << std::scientific
		      << d_energies[i] << std::endl;
	  }
	}
#endif

	// Get the sum of the energies
	float jEnergy = thrust::reduce( d_energies, d_energies+nVoxels );

	// Release thrust arrays
	thrust::device_delete( d_energies );

	return( jEnergy );
      }


      // ------------------------------------------------
    private:
      //! Stream to use for operations
      cudaStream_t stream;

      //! Timer for LLEdispatch
      mutable SciGPU::Utilities::Chronometer tLLEdispatch;

      //! Timer for log likelihood energy
      mutable SciGPU::Utilities::Chronometer tLLEtot;
      //! Timer for LLE 'good' assesments
      mutable SciGPU::Utilities::Chronometer tLLEgood;
      //! Timer for LLE calculation
      mutable SciGPU::Utilities::Chronometer tLLEcompute;


      //! Templated texture binding wrapper
      template<typename T>
      void BindMRI( const GPU::Classes::MRIframeGPU<T>& mri ) const {
	std::cerr << __PRETTY_FUNCTION__
		  << ": Unrecognised MRI type" << std::endl;
	exit( EXIT_FAILURE );
      }

      //! Templated texture unbinding
      template<typename T>
      void UnbindMRI( void ) const {
	std::cerr << __PRETTY_FUNCTION__
		  << ": Unrecognised MRI type" << std::endl;
	exit( EXIT_FAILURE );
      }
    };


    template<>
    void GCAmorphEnergy::BindMRI<unsigned char>( const GPU::Classes::MRIframeGPU<unsigned char>& mri ) const {

      dt_mri_uchar.normalized = false;
      dt_mri_uchar.addressMode[0] = cudaAddressModeClamp;
      dt_mri_uchar.addressMode[1] = cudaAddressModeClamp;
      dt_mri_uchar.addressMode[2] = cudaAddressModeClamp;
      dt_mri_uchar.filterMode = cudaFilterModeLinear;
      
      CUDA_SAFE_CALL( cudaBindTextureToArray( dt_mri_uchar,
					      mri.GetArray() ) );
    }

    template<>
    void GCAmorphEnergy::UnbindMRI<unsigned char>( void ) const {
      CUDA_SAFE_CALL( cudaUnbindTexture( dt_mri_uchar ) );
    }
    
  }
}










static GPU::Algorithms::GCAmorphEnergy myEnergy;


//! Wrapper around GPU class for LLE
float gcamLogLikelihoodEnergyGPU( const GCA_MORPH *gcam,
				  const MRI* mri ) {
  
  float energy;

  switch( mri->type ) {
  
  case MRI_UCHAR:
    energy = myEnergy.LLEdispatch<unsigned char>( gcam, mri );
    break;


  default:
    std::cerr << __FUNCTION__
	      << ": Unrecognised MRI type" << std::endl;
    exit( EXIT_FAILURE );
  }

  return( energy );

}


//! Wrapper around GPU class for JacobianEnergy
float gcamJacobianEnergyGPU( const GCA_MORPH *gcam,
			     const MRI* mri ) {
  
  float energy;
  
  GPU::Classes::GCAmorphGPU myGCAM;
  myGCAM.SendAll( gcam );

  energy = myEnergy.ComputeJacobianEnergy( myGCAM, mri );

  return( energy );
}
