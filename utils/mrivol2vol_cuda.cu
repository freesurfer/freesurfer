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
 *    $Date: 2010/01/27 16:02:34 $
 *    $Revision: 1.10 $
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

    // Some timers
    SciGPU::Utilities::Chronometer tVol2VolMem, tVol2VolMemHost;
    SciGPU::Utilities::Chronometer tVol2VolMRISendFrame, tVol2VolMRIRecvFrame;
    SciGPU::Utilities::Chronometer tVol2VolMRISendArray, tVol2VolCompute;
    SciGPU::Utilities::Chronometer tVol2VolTotal;

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

    
    template<>
    void BindSourceTexture<float>( const GPU::Classes::MRIframeGPU<float>& src,
				   const int InterpMode ) {
      /*!
	Binds the unsigned char texture
      */

      dt_src_float.normalized = false;
      dt_src_float.addressMode[0] = cudaAddressModeClamp;
      dt_src_float.addressMode[1] = cudaAddressModeClamp;
      dt_src_float.addressMode[2] = cudaAddressModeClamp;

      // Set the interpolation
      switch( InterpMode ) {
      case SAMPLE_NEAREST:
	dt_src_float.filterMode = cudaFilterModePoint;
	break;

      case SAMPLE_TRILINEAR:
	dt_src_float.filterMode = cudaFilterModeLinear;
	break;

      default:
	std::cerr << __FUNCTION__ << ": Unrecognised InterpMode "
		  << InterpMode << std::endl;
	exit( EXIT_FAILURE );
      }

      CUDA_SAFE_CALL( cudaBindTextureToArray( dt_src_float, src.GetArray() ) );

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

    template<> void UnbindSourceTexture<float>( void ) {
      CUDA_SAFE_CALL( cudaUnbindTexture( dt_src_float ) );
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
    

    template<>
    __device__ float FetchSrcVoxel<float>( const float3 r ) {
      // Do the fetch, remembering texel centre offset
      float texVal;
      
      texVal = tex3D( dt_src_float, r.x+0.5f, r.y+0.5f, r.z+0.5f );

      // No need to convert - float texture read mode was ElementType
      return( texVal );
    }

    // ---------------------------------------


    //! Kernel to perform the affine transformation
    template<typename T, typename U>
    __global__ void MRIVol2VolKernel( GPU::Classes::MRIframeOnGPU<U> dst,
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
    template<typename T, typename U>
    void MRIVol2VolGPU( const GPU::Classes::MRIframeGPU<T> &src,
			GPU::Classes::MRIframeGPU<U> &dst,
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
      GPU::Classes::MRIframeOnGPU<U> dstGPU(dst);

      // Main dispatch
      switch( InterpMode ) {

      case SAMPLE_NEAREST:
      case SAMPLE_TRILINEAR:
	grid.x = (dstDims.x/kVol2VolBlockSize) * (dstDims.y/kVol2VolBlockSize);
	grid.y = dstDims.z;
	grid.z = 1;
	threads.x = threads.y = kVol2VolBlockSize;
	threads.z = 1;
	MRIVol2VolKernel<T,U><<<grid,threads>>>( dstGPU, transform );
	CUDA_CHECK_ERROR_ASYNC( "MRIVol2VolKernel call failed!\n" );
	break;

      case SAMPLE_SINC:
	std::cerr << __FUNCTION__ << ": Sinc interpolation not yet supported"
		  << std::endl;
	exit( EXIT_FAILURE );
      }


      // Unbind the texture
      UnbindSourceTexture<T>();
    }  


    //! Routine to use GPU for vol2vol on all MRI frames
    template<typename T, typename U>
    void MRIVol2VolAllFramesGPU( const MRI* src, MRI* targ,
				 const MATRIX* transformMatrix,
				 const int InterpMode,
				 const float param ) {
      
      tVol2VolTotal.Start();

      // Get hold of the affine transformation
      GPU::Classes::AffineTransformation myTransform( transformMatrix );

      GPU::Classes::MRIframeGPU<T> srcGPU;
      GPU::Classes::MRIframeGPU<U> dstGPU;

      char* h_workspace;
      size_t srcWorkSize, dstWorkSize;

      // Allocate GPU arrays
      tVol2VolMem.Start();
      srcGPU.Allocate( src );
      srcGPU.AllocateArray();
      dstGPU.Allocate( targ, kVol2VolBlockSize );
      tVol2VolMem.Stop();

      // Sanity check
      srcGPU.VerifyMRI( src );
      dstGPU.VerifyMRI( targ );

      // Allocate workspace array
      tVol2VolMemHost.Start();
      srcWorkSize = srcGPU.GetBufferSize();
      dstWorkSize = dstGPU.GetBufferSize();
      
      if( srcWorkSize > dstWorkSize ) {
	CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_workspace,
				       srcWorkSize,
				       cudaHostAllocDefault ) );
      } else {
	CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_workspace,
				       dstWorkSize,
				       cudaHostAllocDefault ) );
      }
      tVol2VolMemHost.Stop();

      // Loop over frames
      for( unsigned int iFrame=0; iFrame < src->nframes; iFrame++ ) {
	// Get the next source frame
	tVol2VolMRISendFrame.Start();
	srcGPU.Send( src, iFrame, h_workspace );
	tVol2VolMRISendFrame.Stop();
	
	// Put it into a CUDA array
	tVol2VolMRISendArray.Start();
	srcGPU.SendArray();
	tVol2VolMRISendArray.Stop();

	// Run the convolution
	tVol2VolCompute.Start();
	MRIVol2VolGPU( srcGPU, dstGPU, myTransform, InterpMode );
	tVol2VolCompute.Stop();

	// Get the results back
	tVol2VolMRIRecvFrame.Start();
	dstGPU.Recv( targ, iFrame, h_workspace );
	tVol2VolMRIRecvFrame.Stop();
      }

      CUDA_CHECK_ERROR( "MRI Vol2Vol failed on GPU" );

      // Release workspace array
      tVol2VolMemHost.Start();
      CUDA_SAFE_CALL( cudaFreeHost( h_workspace ) );
      tVol2VolMemHost.Stop();

      // No need to release MRIframeGPU types - destructors will handle it
      tVol2VolTotal.Stop();
    }

    //! Dispatch routine to add templating
    template<typename T>
    void MRIVol2VolAllFramesDstDispatch( const MRI* src, MRI* targ,
					 const MATRIX* transformMatrix,
					 const int InterpMode,
					 const float param ) {
      /*!
	A thin wrapper to produce the correct template
	based on the type of the target MRI
      */

      switch( targ->type ) {
      case MRI_UCHAR:
	MRIVol2VolAllFramesGPU<T,unsigned char>( src, targ,
						 transformMatrix,
						 InterpMode,
						 param );
	break;

      case MRI_FLOAT:
	MRIVol2VolAllFramesGPU<T,float>( src, targ,
					 transformMatrix,
					 InterpMode,
					 param );
	break;

      default:
	std::cerr << __FUNCTION__ << ": Unrecognised target type "
		  << targ->type << std::endl;
	exit( EXIT_FAILURE );
      }
    }

  }
}


// ======================================================

int MRIvol2vol_cuda( const MRI* src, MRI* targ, 
		     const MATRIX* transformMatrix,
		     const int InterpMode,
		     const float param ) {
  /*!
    Reimplementation of MRIvol2vol for the GPU.
    Is meant to be called from within MRIvol2vol,
    once that routine has performed necessary checks
  */

  switch( src->type ) {
  case MRI_UCHAR:
    GPU::Algorithms::MRIVol2VolAllFramesDstDispatch<unsigned char>( src, targ,
								    transformMatrix,
								    InterpMode,
								    param );
    break;

  case MRI_FLOAT:
    GPU::Algorithms::MRIVol2VolAllFramesDstDispatch<float>( src, targ,
							    transformMatrix,
							    InterpMode,
							    param );
    break;

  default:
    std::cerr << __FUNCTION__ << ": Unrecognised data type "
	      << src->type << std::endl;
    exit( EXIT_FAILURE );
  }
  


  return( 0 );
}
    



// ======================================================


// ======================================================

//! Stream insertion operator for timer
static std::ostream& operator<<( std::ostream& os,
				 const SciGPU::Utilities::Chronometer& timer ) {
  
  os << std::setw(9) << std::setprecision(6) << timer.GetAverageTime() << " ms (avg) ";
  os << std::setw(9) << std::setprecision(6) << timer.GetTime() << " ms (tot)";

  return( os );
}


void MRIvol2volShowTimers( void ) {
  /*!
    Pretty prints timers to std.out
  */
  

  std::cout << "=============================================" << std::endl;
  std::cout << "GPU Vol2Vol timers" << std::endl;
  std::cout << "------------------" << std::endl;
#ifndef CUDA_FORCE_SYNC
  std::cout << "WARNING: CUDA_FORCE_SYNC not #defined" << std::endl;
  std::cout << "Timings may not be accurate" << std::endl;
#endif
  std::cout << std::endl;

  std::cout << "MRIVol2VolAllFramesGPU" << std::endl;
  std::cout << "Host Memory : " << GPU::Algorithms::tVol2VolMemHost << std::endl;
  std::cout << "GPU Memory  : " << GPU::Algorithms::tVol2VolMem << std::endl;
  std::cout << "Send Frame  : " << GPU::Algorithms::tVol2VolMRISendFrame << std::endl;
  std::cout << "Send Array  : " << GPU::Algorithms::tVol2VolMRISendArray << std::endl;
  std::cout << "Compute     : " << GPU::Algorithms::tVol2VolCompute << std::endl;
  std::cout << "Recv Frame  : " << GPU::Algorithms::tVol2VolMRIRecvFrame << std::endl;
  std::cout << "---------------------" << std::endl;
  std::cout << "Total : " << GPU::Algorithms::tVol2VolTotal << std::endl;

  std::cout << "=============================================" << std::endl;
}
