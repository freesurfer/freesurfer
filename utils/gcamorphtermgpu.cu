/**
 * @file  gcamorphtermgpu.cu
 * @brief Holds routines to compute GCAmorph terms on the GPU
 *
 * This file holds a variety of routines which compute terms for
 * mri_ca_register
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/05/14 17:50:13 $
 *    $Revision: 1.2 $
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

#include "macros.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"

#include "gcamorphtermgpu.hpp"




//! Texture for vx in Smoothness
texture<float, 3, cudaReadModeElementType> dt_smooth_vx;
//! Texture for vy in Smoothness
texture<float, 3, cudaReadModeElementType> dt_smooth_vy;
//! Texture for vz in Smoothness
texture<float, 3, cudaReadModeElementType> dt_smooth_vz;
//! Texure for 'invalid' in Smoothness
texture<char, 3, cudaReadModeElementType> dt_smooth_invalid;


// ==============================================================


namespace GPU {
  namespace Algorithms {

    const unsigned int kGCAmorphSmoothTermKernelSize = 16;


    // ##############################################################

    //! Kernel to perform \f$c = a - b\f$ on volumes
    template<typename T>
    __global__
    void SubtractVolumes( const GPU::Classes::VolumeArgGPU<T> a,
			  const GPU::Classes::VolumeArgGPU<T> b,
			  GPU::Classes::VolumeArgGPU<T> c ) {

      const unsigned int bx = ( blockIdx.x * blockDim.x );
      const unsigned int by = ( blockIdx.y * blockDim.y );
      const unsigned int ix = threadIdx.x + bx;
      const unsigned int iy = threadIdx.y + by;

      for( unsigned int iz = 0; iz < a.dims.z; iz++ ) {
	if( a.InVolume(ix,iy,iz) ) {

	  c(ix,iy,iz) = a(ix,iy,iz) - b(ix,iy,iz);
	}
      }

    }


    __device__ float Fetchvx( const int ix,
			      const int iy,
			      const int iz ) {
      return( tex3D( dt_smooth_vx, ix+0.5f, iy+0.5f, iz+0.5f ) );
    }

    __device__ float Fetchvy( const int ix,
			      const int iy,
			      const int iz ) {
      return( tex3D( dt_smooth_vy, ix+0.5f, iy+0.5f, iz+0.5f ) );
    }
    
    __device__ float Fetchvz( const int ix,
			      const int iy,
			      const int iz ) {
      return( tex3D( dt_smooth_vz, ix+0.5f, iy+0.5f, iz+0.5f ) );
    }



    //! Kernel to compute the smoothness term
    __global__
    void SmoothnessTermKernel( GPU::Classes::VolumeArgGPU<float> dx,
			       GPU::Classes::VolumeArgGPU<float> dy,
			       GPU::Classes::VolumeArgGPU<float> dz,
			       const float l_smoothness ) {
      /*!
	Main computation of the smoothness term.
      */

      const unsigned int bx = ( blockIdx.x * blockDim.x );
      const unsigned int by = ( blockIdx.y * blockDim.y );
      const unsigned int ix = threadIdx.x + bx;
      const unsigned int iy = threadIdx.y + by;

      for( unsigned int iz = 0; iz < dx.dims.z; iz++ ) {
	// Only compute if ix, iy & iz are inside the bounding box
	if( dx.InVolume(ix,iy,iz) ) {

	  // Only calculate if we're not invalid
	  const char myValid = tex3D( dt_smooth_invalid,
				      ix+0.5f, iy+0.5f, iz+0.5f );
	  if( myValid != GCAM_POSITION_INVALID ) {
	    
	    float vx = Fetchvx(ix,iy,iz);
	    float vy = Fetchvy(ix,iy,iz);
	    float vz = Fetchvz(ix,iy,iz);

	    float ddx, ddy, ddz;
	    ddx = ddy = ddz = 0;

	    unsigned int num = 0;

	    // Re-order loop nest in gcamSmoothnessEnergy
	    for( int zk=-1; zk<=1; zk++ ) {
	      int zn = iz + zk;
	      
	      for( int yk=-1; yk<=1; yk++ ) {
		int yn = iy + yk;
		
		for( int xk=-1; xk<=1; xk++ ) {
		  int xn = ix + xk;
		
		  // Don't include self
		  if( (!xk) && (!yk) && (!zk) ) {
		    continue;
		  }
		  
		  // Don't include invalid neighbours
		  
		  // Use texture to get boundary conditions
		  const char neighbourInvalid = tex3D( dt_smooth_invalid,
						       xn+0.5f,
						       yn+0.5f,
						       zn+0.5f );
		  if( neighbourInvalid == GCAM_POSITION_INVALID ) {
		    continue;
		  }


		  ddx += Fetchvx(xn,yn,zn) - vx;
		  ddy += Fetchvy(xn,yn,zn) - vy;
		  ddz += Fetchvz(xn,yn,zn) - vz;

		  num++;

		}
	      }
	    }

	    if( num > 0 ) {
	      // Rescale
	      ddx *= l_smoothness / num;
	      ddy *= l_smoothness / num;
	      ddz *= l_smoothness / num;

	      // Update the GCAmorph
	      dx(ix,iy,iz) += ddx;
	      dy(ix,iy,iz) += ddy;
	      dz(ix,iy,iz) += ddz;
	    }

	  } // End if( node valid )


	} // End if( inVolume )

      } // End loop over z
    }



    void GCAmorphTerm::Smoothness( GPU::Classes::GCAmorphGPU& gcam,
				   const float l_smoothness ) const {
      /*!
	Computes ther smoothness term of the given gcam, following
	the routine gcamSmoothnessTerm in gcamorph.c.
      */

      if( DZERO(l_smoothness) ) {
	return;
      }

      GCAmorphTerm::tSmoothTot.Start();

      // Make sure GCAM is sane
      gcam.CheckIntegrity();
      
      const dim3 gcamDims = gcam.d_rx.GetDims();
      
      // Compute vx, vy, and vz (see gcamSmoothnessEnergy)
      GPU::Classes::VolumeGPU<float> vx, vy, vz;
      
      vx.Allocate( gcamDims );
      vy.Allocate( gcamDims );
      vz.Allocate( gcamDims );
      vx.AllocateArray();
      vy.AllocateArray();
      vz.AllocateArray();
      
      dim3 grid, threads;
      threads.x = threads.y = kGCAmorphSmoothTermKernelSize;
      threads.z = 1;
      
      grid = gcam.d_rx.CoverBlocks( kGCAmorphSmoothTermKernelSize );
      grid.z = 1;
      
      GCAmorphTerm::tSmoothSubtract.Start();
      SubtractVolumes<float>
	<<<grid,threads>>>
	( gcam.d_rx, gcam.d_origx, vx );
      CUDA_CHECK_ERROR( "SubtractVolumes failed for x!" );
      
      SubtractVolumes<float>
	<<<grid,threads>>>
	  ( gcam.d_ry, gcam.d_origy, vy );
      CUDA_CHECK_ERROR( "SubtractVolumes failed for y!" );
      
      SubtractVolumes<float>
	<<<grid,threads>>>
	( gcam.d_rz, gcam.d_origz, vz );
      CUDA_CHECK_ERROR( "SubtractVolumes failed for z!" );
      GCAmorphTerm::tSmoothSubtract.Stop();
      
      // Get vx, vy and vz into CUDA arrays
      vx.SendArray();
      vy.SendArray();
      vz.SendArray();
      
      // Bind vx, vy and vz to their textures
      dt_smooth_vx.normalized = false;
      dt_smooth_vx.addressMode[0] = cudaAddressModeClamp;
      dt_smooth_vx.addressMode[1] = cudaAddressModeClamp;
      dt_smooth_vx.addressMode[2] = cudaAddressModeClamp;
      dt_smooth_vx.filterMode = cudaFilterModePoint;
      
      dt_smooth_vy.normalized = false;
      dt_smooth_vy.addressMode[0] = cudaAddressModeClamp;
      dt_smooth_vy.addressMode[1] = cudaAddressModeClamp;
      dt_smooth_vy.addressMode[2] = cudaAddressModeClamp;
      dt_smooth_vy.filterMode = cudaFilterModePoint;
      
      dt_smooth_vz.normalized = false;
      dt_smooth_vz.addressMode[0] = cudaAddressModeClamp;
      dt_smooth_vz.addressMode[1] = cudaAddressModeClamp;
      dt_smooth_vz.addressMode[2] = cudaAddressModeClamp;
      dt_smooth_vz.filterMode = cudaFilterModePoint;
      
      CUDA_SAFE_CALL( cudaBindTextureToArray( dt_smooth_vx,
					      vx.GetArray() ) );
      CUDA_SAFE_CALL( cudaBindTextureToArray( dt_smooth_vy,
					      vy.GetArray() ) );
      CUDA_SAFE_CALL( cudaBindTextureToArray( dt_smooth_vz,
					      vz.GetArray() ) );
      
      // Also have to get the 'invalid' field to its texture
      dt_smooth_invalid.normalized = false;
      dt_smooth_invalid.addressMode[0] = cudaAddressModeClamp;
      dt_smooth_invalid.addressMode[1] = cudaAddressModeClamp;
      dt_smooth_invalid.addressMode[2] = cudaAddressModeClamp;
      dt_smooth_invalid.filterMode = cudaFilterModePoint;
      
      gcam.d_invalid.AllocateArray();
      gcam.d_invalid.SendArray();
      CUDA_SAFE_CALL( cudaBindTextureToArray( dt_smooth_invalid,
					      gcam.d_invalid.GetArray() ) );
      

      // Run the main kernel
      GCAmorphTerm::tSmoothCompute.Start();
      SmoothnessTermKernel<<<grid,threads>>>
	( gcam.d_dx, gcam.d_dy, gcam.d_dz, l_smoothness );
      CUDA_CHECK_ERROR( "SmoothnessTermKernelFailed!" );
      GCAmorphTerm::tSmoothCompute.Stop();



      // Unbind textures
      CUDA_SAFE_CALL( cudaUnbindTexture( dt_smooth_vx ) );
      CUDA_SAFE_CALL( cudaUnbindTexture( dt_smooth_vy ) );
      CUDA_SAFE_CALL( cudaUnbindTexture( dt_smooth_vz ) );
      CUDA_SAFE_CALL( cudaUnbindTexture( dt_smooth_invalid ) );


      GCAmorphTerm::tSmoothTot.Stop();

    }




    // ##############################################################
    
    void GCAmorphTerm::ShowTimings( void ) {

#ifdef CUDA_SHOW_TIMINGS
      std::cout << "==================================" << std::endl;
      std::cout << "GPU GCAmorph term timers" << std::endl;
      std::cout << "------------------------" << std::endl;
#ifndef CUDA_FORCE_SYNC
      std::cout << "WARNING: CUDA_FORCE_SYNC not #defined" << std::endl;
      std::cout << "Timings may not be accurate" << std::endl;
#endif
      std::cout << std::endl;

      std::cout << "Smoothness:" << std::endl;
      std::cout << "  Subtraction : " << GCAmorphTerm::tSmoothSubtract
		<< std::endl;
      std::cout << "      Compute : " << GCAmorphTerm::tSmoothCompute
		<< std::endl;
      std::cout << "Total         : " << GCAmorphTerm::tSmoothTot
		<< std::endl;
      std::cout << std::endl;


     std::cout << "==================================" << std::endl;
#endif
    }


    // ##############################################################

    // Declare static members

    SciGPU::Utilities::Chronometer GCAmorphTerm::tSmoothTot;
    SciGPU::Utilities::Chronometer GCAmorphTerm::tSmoothSubtract;
    SciGPU::Utilities::Chronometer GCAmorphTerm::tSmoothCompute;
    

  }
}





// --------------------------------------------
static GPU::Algorithms::GCAmorphTerm myTerms;


//! Wrapper around GPU class for smoothness term

void gcamSmoothnessTermGPU( GCA_MORPH *gcam,
			    const float l_smoothness ) {
  
  GPU::Classes::GCAmorphGPU myGCAM;
  
  myGCAM.SendAll( gcam );
  
  myTerms.Smoothness( myGCAM, l_smoothness );
  
  myGCAM.RecvAll( gcam );
  
}
