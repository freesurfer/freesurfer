/**
 * @file  mrimean_cuda.cu
 * @brief Holds MRI mean function for the GPU
 *
 * Implements MRI mean function on the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/02/02 21:03:03 $
 *    $Revision: 1.5 $
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
}

//#include "cuPrintf.cu"

#include "chronometer.hpp"
#include "cudacheck.h"


#include "mriframegpu.hpp"


#include "mrimean_cuda.h"






namespace GPU {
  namespace Algorithms {

    const unsigned int kMRImeanBlockSize = 16;

    //! Min function for ints
    __device__ int min( const int& a, const int& b ) {
      if( a < b ) {
	return( a );
      } else { 
	return( b );
      }
    }

    //! Max function for ints
    __device__ int max( const int& a, const int& b ) {
      if( a > b ) {
	return( a );
      } else { 
	return( b );
      }
    }

    //! Function to round up a value to multiple of given number
    template<typename T>
    __device__ T Roundup( const T val, const unsigned int num ) {
      
      float tmp;
      
      tmp = static_cast<float>( val ) / num;

      T tmpVal = static_cast<T>( ceilf( tmp ) );

      return( num * tmpVal );
    }

    //! Kernel to perform the mean filter
    template<typename T, typename U>
    __global__ void MRImeanKernel( const GPU::Classes::MRIframeOnGPU<T> src,
				   GPU::Classes::MRIframeOnGPU<U> dst,
				   const dim3 actualDims,
				   const unsigned int wSize ) {
      /*!
	Kernel to do the mean filtering.
	Each block will do one 16x16 (x,y) patch of the dst MRI.
	The z co-ordinate will be in blockIdx.y.
	The (x,y) location of the patch will be derived from blockIdx.x.
	The data dimensions on the GPU will be padded to multiples
	of 16 by the calling routine, but we still need the actual
	dimensions on the CPU
      */

      const unsigned int by = blockIdx.x / (src.dims.x/kMRImeanBlockSize);
      const unsigned int bx = blockIdx.x % (src.dims.x/kMRImeanBlockSize);
      const int iz = blockIdx.y;
      const unsigned int tx = threadIdx.x;
      const unsigned int ty = threadIdx.y;
      
      const unsigned int ixStart = bx * kMRImeanBlockSize;
      const unsigned int iyStart = by * kMRImeanBlockSize;
      const int ix = ixStart + tx;
      const int iy = iyStart + ty;

      // Accumulator
      float myVal = 0;

      // Declared int to eliminate negative number problems
      const int wHalf = wSize / 2;

      // Calculate voxels which will contribute, clamping to edges
      const unsigned int myxMin = max( 0           , ix - wHalf );
      const unsigned int myxMax = min( actualDims.x, ix + wHalf );
      const unsigned int myyMin = max( 0           , iy - wHalf );
      const unsigned int myyMax = min( actualDims.y, iy + wHalf );
      const unsigned int myzMin = max( 0           , iz - wHalf );
      const unsigned int myzMax = min( actualDims.z, iz + wHalf );

      // Again, declare int to remove need for some casts
      const int patchSize = Roundup( max(wHalf,1), kMRImeanBlockSize );

      // Loop over z levels (iz is the same for all threads in the block)
      for( unsigned int currZ = max( 0, iz - wHalf );
	   currZ <= min( actualDims.z-1, iz + wHalf );
	   currZ++ ) {

	__shared__ float currPatch[kMRImeanBlockSize][kMRImeanBlockSize];

	// Compute loop limits
	// Note that some of the arguments are declared integer to avoid need for casts
	const int xFirst = max( 0, ixStart - patchSize );
	const int xLast = min( src.dims.x - patchSize, ixStart + patchSize );

	const int yFirst = max( 0, iyStart - patchSize );
	const int yLast = min( src.dims.y - patchSize, iyStart + patchSize );

	// Loop over patches
	for( int yBegin = yFirst; yBegin <= yLast; yBegin += kMRImeanBlockSize ) {
	  for( int xBegin = xFirst; xBegin <= xLast; xBegin += kMRImeanBlockSize ) {
	    // Load up the patch
	    currPatch[ty][tx] = src( xBegin+tx, yBegin+ty, currZ );
	    __syncthreads();

	    // Accumulate within the patch
	    for( int curry = 0; curry < kMRImeanBlockSize; curry++ ) {
	      for( int currx = 0; currx < kMRImeanBlockSize; currx++ ) {
		
		int actx = xBegin+currx;
		int acty = yBegin+curry;

		if( (actx>=myxMin) && (actx<=myxMax) &&
		    (acty>=myyMin) && (acty<=myyMax) ) {
		  myVal += currPatch[curry][currx];
		}

	      }
	    }
	    __syncthreads();
	  }
	}


      }

      const unsigned long myVolume = ( myxMax - myxMin + 1 ) *
	(myyMax - myyMin + 1 ) *
	(myzMax - myzMin + 1 );

      dst( ix, iy, iz ) = dst.ConvertFloat( myVal / myVolume );
    }



    //! Runs the mean filtering kernel
    template<typename T, typename U>
    void MRImeanGPU( const GPU::Classes::MRIframeGPU<T> &src,
		     GPU::Classes::MRIframeGPU<U> &dst,
		     const unsigned int wSize,  
		     const cudaStream_t myStream = 0 ) {
      /*!
	Runs the mean filtering kernel on the GPU.
	Assumes most things are properly set already
      */
      const dim3 srcDims = src.GetGPUDims();
      const dim3 dstDims = dst.GetGPUDims();
      const dim3 srcCPUdims = src.GetCPUDims();

      // Check dimensions
      if( (srcDims.x != dstDims.x) &&
	  (srcDims.y != dstDims.y) &&
	  (srcDims.z != dstDims.z) ) {
	std::cerr << __FUNCTION__ << ": Dimension mismatch"
		  << std::endl;
	exit( EXIT_FAILURE );
      }
      
      // Check padding (only need to do one, given above check)
      if( !dst.CheckPadding( kMRImeanBlockSize ) ) {
	std::cerr << __FUNCTION__
		  << ": Arrays on GPU must be padded to multiples"
		  << " of kMRImeanBlockSize" << std::endl;
	exit( EXIT_FAILURE );
      }

      dim3 grid, threads;
      GPU::Classes::MRIframeOnGPU<T> srcGPU(src);
      GPU::Classes::MRIframeOnGPU<U> dstGPU(dst);

      grid.x = (dstDims.x/kMRImeanBlockSize) * (dstDims.y/kMRImeanBlockSize);
      grid.y = dstDims.z;
      grid.z = 1;
      threads.x = threads.y = kMRImeanBlockSize;
      threads.z = 1;

      MRImeanKernel<T,U><<<grid,threads,0,myStream>>>( srcGPU, dstGPU, srcCPUdims, wSize );
      CUDA_CHECK_ERROR_ASYNC( "MRImeanKernel failed!" );
    }
    


    //! Dispatch wrapper for MRI mean on the GPU
    template<typename T, typename U>
    void MRImeanDispatch( const MRI* src, MRI* dst,
			  const unsigned int wSize,
			  const int srcFrame, const int dstFrame ) {
      /*!
	Templated dispatch routine for MRI mean function on the
	GPU.
	Transfers data to the GPU, and retrieves the results
      */

      GPU::Classes::MRIframeGPU<T> srcGPU;
      GPU::Classes::MRIframeGPU<U> dstGPU;
      
      char* h_workspace;
      size_t srcWorkSize, dstWorkSize;
      
      // Allocate the GPU arrays
      srcGPU.Allocate( src, kMRImeanBlockSize );
      dstGPU.Allocate( dst, kMRImeanBlockSize );


      // Put in some sanity checks
      srcGPU.VerifyMRI( src );
      dstGPU.VerifyMRI( dst );
      
      // Allocate the workspace array
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

      // Send the source data
      srcGPU.Send( src, srcFrame, h_workspace );

      // Run the filter
      MRImeanGPU( srcGPU, dstGPU, wSize );

      // Get the results
      dstGPU.Recv( dst, dstFrame, h_workspace );

      // Release workspace array
      CUDA_SAFE_CALL( cudaFreeHost( h_workspace ) );

      CUDA_CHECK_ERROR( "Mean filtering failure" );
    }



    //! Wrapper for MRI mean
    template<typename T>
    void MRImeanDispatchWrap( const MRI* src, MRI* dst,
			      const unsigned int wSize,
			      const int srcFrame, const int dstFrame ) {
      switch( dst->type ) {
      case MRI_UCHAR:
	MRImeanDispatch<T,unsigned char>( src, dst, wSize, srcFrame, dstFrame );
	break;

      case MRI_SHORT:
	MRImeanDispatch<T,short>( src, dst, wSize, srcFrame, dstFrame );
	break;

      case MRI_FLOAT:
	MRImeanDispatch<T,float>( src, dst, wSize, srcFrame, dstFrame );
	break;

      default:
	std::cerr << __FUNCTION__
		  << ": Unrecognised destination MRI type "
		  << dst->type
		  << std::endl;
	exit( EXIT_FAILURE );
      }
    }

  }
}


MRI* MRImean_cuda( const MRI* src, MRI* dst,
		   const unsigned int wSize ) {
  /*!
    Wrapper around the GPU routine, to be called from the
    original MRImean routine.
    Note that the frames default to zero, per the original
    MRImean routine.
  */

  switch( src->type ) {
  case MRI_UCHAR:
    GPU::Algorithms::MRImeanDispatchWrap<unsigned char>( src, dst, wSize, 0, 0 );
    break;

  case MRI_SHORT:
    GPU::Algorithms::MRImeanDispatchWrap<short>( src, dst, wSize, 0, 0 );
    break;

  case MRI_FLOAT:
    GPU::Algorithms::MRImeanDispatchWrap<float>( src, dst, wSize, 0, 0 );
    break;

  default:
    std::cerr << __FUNCTION__
	      << ": Unrecognised source MRI type "
	      << src->type
	      << std::endl;
    exit( EXIT_FAILURE );
  }

  return( dst );
}
