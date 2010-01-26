/**
 * @file  mrivol2vol_cuda.cu
 * @brief Holds MRI volume to volume transform for CUDA
 *
 * Implements MRIvol2vol for the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/01/26 16:54:37 $
 *    $Revision: 1.7 $
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
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>


extern "C" {
#include "mri.h"
#include "matrix.h"
}


#include "chronometer.hpp"
#include "cudacheck.h"


#include "mriframegpu.hpp"
#include "affinegpu.hpp"

#include "mrivol2vol_cuda.h"

// ==============================================================

//! Texture reference for an unsigned char source
texture<unsigned char, 3, cudaReadModeNormalizedFloat> dt_src_uchar;
//! Texture reference for a float source
texture<float, 3, cudaReadModeElementType> dt_src_float;


namespace GPU {
  namespace Algorithms {

    const unsigned int kVol2VolBlockSize = 16;

    // ---------------------------------------

    //! Templated texture binding wrapper
    template<typename T>
    void BindSourceTexture( const GPU::Classes::MRIframeGPU<T>& src,
			    const int InterpMode ) {
      /*!
	Will bind the CUDA array of the source frame to the appropriate
	texture, with correct filtering set
	Unspecialised version aborts
      */
      std::cerr << __PRETTY_FUNCTION__ << ": Unrecognised type" << std::endl;
    }


    template<>
    void BindSourceTexture<unsigned char>( const GPU::Classes::MRIframeGPU<unsigned char>& src,
					   const int InterpMode ) {
      /*!
	Binds the unsigned char texture
      */

      dt_src_uchar.normalized = false;
      dt_src_uchar.addressMode[0] = cudaAddressModeClamp;
      dt_src_uchar.addressMode[1] = cudaAddressModeClamp;
      dt_src_uchar.addressMode[2] = cudaAddressModeClamp;

      // Set the interpolation
      switch( InterpMode ) {
      case SAMPLE_NEAREST:
	dt_src_uchar.filterMode = cudaFilterModePoint;
	break;

      case SAMPLE_TRILINEAR:
	dt_src_uchar.filterMode = cudaFilterModeLinear;
	break;

      default:
	std::cerr << __FUNCTION__ << ": Unrecognised InterpMode "
		  << InterpMode << std::endl;
	exit( EXIT_FAILURE );
      }

      CUDA_SAFE_CALL( cudaBindTextureToArray( dt_src_uchar, src.GetArray() ) );

    }

    // ---------------------------------------

    //! Templated texture unbinding
    template<typename T>
    void UnbindSourceTexture( void ) {
      /*!
	Will unbind the appropriate texture from the source frame.
	Unspecialised version aborts
      */
      std::cerr << __PRETTY_FUNCTION__ << ": Unrecognised type" << std::endl;
    }

    template<> void UnbindSourceTexture<unsigned char>( void ) {
      CUDA_SAFE_CALL( cudaUnbindTexture( dt_src_uchar ) );
    }

    // ---------------------------------------

    //! Templated texture fetch
    template<typename T>
    __device__ float FetchSrcVoxel( const float3 r ) {
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
    __device__ float FetchSrcVoxel<unsigned char>( const float3 r ) {

      float texVal;
      
      // Do the fetch, remembering texel centre offset
      texVal = tex3D( dt_src_uchar, r.x+0.5f, r.y+0.5f, r.z+0.5f );

      // We configured to return a normalised float, so get back in range
      texVal *= UCHAR_MAX;

      return( texVal );
    }
      

    //! Kernel to perform the affine transformation
    template<typename T>
    __global__ void MRIVol2VolKernel( GPU::Classes::MRIframeOnGPU<T> dst,
				      GPU::Classes::AffineTransformation afTrans ) {
      /*!
	Kernel to perform the affine transformation supplied on the
	MRI frame bound to the textures.
	Stores the result in the dst frame.
	Each thread produces one element of the result.
	The z co-ordinate will be in blockIdx.y
	The (x,y) location of the patch will be derived from blockIdx.x.
	Since the slices array is padded out to multiples of 16 by
	the other MRI slices routines (the host routine checks this), we
	don't have to worry about nasty edge cases
      */

      // Extract some co-ordinates
      const unsigned int by = blockIdx.x / (dst.dims.x/kVol2VolBlockSize);
      const unsigned int bx = blockIdx.x % (dst.dims.x/kVol2VolBlockSize);
      const unsigned int iz = blockIdx.y;
      const unsigned int tx = threadIdx.x;
      const unsigned int ty = threadIdx.y;

      const unsigned int ixStart = bx * kVol2VolBlockSize;
      const unsigned int iyStart = by * kVol2VolBlockSize;
      const unsigned int myx = ixStart + tx;
      const unsigned int myy = iyStart + ty;

      float3 rOut = afTrans.transform( make_float3( myx, myy, iz ) );

      dst( myx, myy, iz ) = dst.ConvertFloat( FetchSrcVoxel<T>( rOut ) );
    }

    //! Dispatch function for the transformation kernel
    template<typename T>
    void MRIVol2VolGPU( const GPU::Classes::MRIframeGPU<T> &src,
			GPU::Classes::MRIframeGPU<T> &dst,
			const GPU::Classes::AffineTransformation& transform,
			const int InterpMode,
			const cudaStream_t myStream = 0 ){
      /*!
	This is the dispatch function for the vol2vol kernel.
	The src MRI must already be in its CUDA array, and the
	dst MRI must be padded
      */

      // Run through a couple of sanity checks
      if( !src.HasArray() ) {
	std::cerr << __FUNCTION__ <<
	  ": Source array must be in CUDA array" << std::endl;
	exit( EXIT_FAILURE );
      }

      if( !dst.CheckPadding( kVol2VolBlockSize ) ) {
	std::cerr << __FUNCTION__ <<
	  ": Destination array must be padded" << std::endl;
	exit( EXIT_FAILURE );
      }
      
      // Bind the texture and array
      BindSourceTexture( src, InterpMode );


      const dim3 dstDims = dst.GetGPUDims();
      dim3 grid, threads;
      GPU::Classes::MRIframeOnGPU<T> dstGPU(dst);

      // Main dispatch
      switch( InterpMode ) {

      case SAMPLE_NEAREST:
      case SAMPLE_TRILINEAR:
	grid.x = (dstDims.x/kVol2VolBlockSize) * (dstDims.y/kVol2VolBlockSize);
	grid.y = dstDims.z;
	grid.z = 1;
	threads.x = threads.y = kVol2VolBlockSize;
	threads.z = 1;
	MRIVol2VolKernel<<<grid,threads>>>( dstGPU, transform );
	CUDA_CHECK_ERROR( "MRIVol2VolKernel call failed!\n" );
	break;

      case SAMPLE_SINC:
	std::cerr << __FUNCTION__ << ": Sinc interpolation not yet supported"
		  << std::endl;
	exit( EXIT_FAILURE );
      }


      // Unbind the texture
      UnbindSourceTexture<T>();

    }  

  }
}


int MRIvol2vol_cuda( const MRI* src, MRI* targ, 
		     const MATRIX* transformMatrix,
		     const int InterpMode,
		     const float param ) {
  /*!
    Reimplementation of MRIvol2vol for the GPU.
    Is meant to be called from within MRIvol2vol,
    once that routine has performed necessary checks
  */

  std::cout << __FUNCTION__ << ": Begin" << std::endl;

  // Need to track 'actual' and GPU dimensions
  // Recall that GPU dimensions will be padded
  // to multiples of 16
  
  
  GPU::Classes::AffineTransformation myTransform( transformMatrix );

  // Get the source MRI onto the GPU

  // Promote source MRI data to float on the GPU

  // Release original src data on GPU

  // Allocate destination float array on GPU

  // Allocate CUDA array to hold one frame of data


  // Loop over frames
  for( unsigned int iFrame=0; iFrame < src->nframes; iFrame++ ) {
    GPU::Classes::MRIframeGPU<unsigned char> src, dst;
    
    GPU::Algorithms::MRIVol2VolGPU( src, dst, myTransform, InterpMode );
    // Copy src frame into CUDA array

    // Bind the array to a texture

    // Run the kernel on this frame

    // Unbind texture and array
  }

  // Release CUDA array

  // Release src float data

  // Convert destination data back to desired data type

  // Get destination data back from the GPU

  return( 0 );
}
    
