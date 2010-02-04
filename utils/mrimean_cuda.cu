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
 *    $Date: 2010/02/04 18:05:23 $
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

    
    //! Kernel to compute x direction means
    template<typename T>
    __global__ void MRImeanXkernel( const GPU::Classes::MRIframeOnGPU<T> src,
				    GPU::Classes::MRIframeOnGPU<float> dst,
				    const dim3 actualDims,
				    const unsigned int wSize ) {
      /*!
	Kernel to compute means in the x direction, based on
	the given window size.
	Basically, does a 1D convolution, but with different
	boundary conditions to MRIConvolveKernelX.
	Also, since this is meant to be part of a pipeline,
	the destination type must be float
      */
      const unsigned int by = blockIdx.x / (src.dims.x/kMRImeanBlockSize);
      const unsigned int bx = blockIdx.x % (src.dims.x/kMRImeanBlockSize);
      const unsigned int tx = threadIdx.x;
      const unsigned int ty = threadIdx.y;
      
      const int ixStart = bx * kMRImeanBlockSize;
      const int iyStart = by * kMRImeanBlockSize;

      const int ix = ixStart + tx;
      const int iy = iyStart + ty;
      const int iz = blockIdx.y;

      const int wHalf = wSize/2;

      // Calculate voxels which will contribute, clamping to edges
      const unsigned int myxMin = max( 0           , ix - wHalf );
      const unsigned int myxMax = min( actualDims.x, ix + wHalf );

      // Again, declare int to remove need for some casts
      const int patchSize = Roundup( max(wHalf,1), kMRImeanBlockSize );

      // Accumulator
      float myVal = 0;

      __shared__ float currPatch[kMRImeanBlockSize][kMRImeanBlockSize];

      // Calculate patch limits (note integer declarations avoid -ve trouble)
      const int xFirst = max( 0, ixStart - patchSize );
      const int xLast = min( src.dims.x - patchSize, ixStart + patchSize );

      for( int xBegin = xFirst; xBegin <= xLast; xBegin += kMRImeanBlockSize ) {
	// Load the patch
	currPatch[ty][tx] = src( xBegin+tx, iy, iz );
	__syncthreads();

	// Accumulate desired values
	for( unsigned int i=0; i<kMRImeanBlockSize; i++ ) {
	  int actx = xBegin + i;

	  if( (actx>=myxMin) && (actx<=myxMax) ) {
	    myVal += currPatch[ty][i];
	  }

	}

	__syncthreads();
      }

      // Save result
      dst(ix,iy,iz) = dst.ConvertFloat( myVal );
    }

    
    //! Kernel to compute y direction means
    template<typename T>
    __global__ void MRImeanYkernel( const GPU::Classes::MRIframeOnGPU<T> src,
				    GPU::Classes::MRIframeOnGPU<float> dst,
				    const dim3 actualDims,
				    const unsigned int wSize ) {
      /*!
	Kernel to compute means in the y direction, based on
	the given window size.
	Basically, does a 1D convolution, but with different
	boundary conditions to MRIConvolveKernelY.
	Also, since this is meant to be part of a pipeline,
	the destination type must be float
      */
      const unsigned int by = blockIdx.x / (src.dims.x/kMRImeanBlockSize);
      const unsigned int bx = blockIdx.x % (src.dims.x/kMRImeanBlockSize);
      const unsigned int tx = threadIdx.x;
      const unsigned int ty = threadIdx.y;
      
      const int ixStart = bx * kMRImeanBlockSize;
      const int iyStart = by * kMRImeanBlockSize;

      const int ix = ixStart + tx;
      const int iy = iyStart + ty;
      const int iz = blockIdx.y;

      const int wHalf = wSize/2;

      // Calculate voxels which will contribute, clamping to edges
      const unsigned int myyMin = max( 0           , iy - wHalf );
      const unsigned int myyMax = min( actualDims.x, iy + wHalf );

      // Again, declare int to remove need for some casts
      const int patchSize = Roundup( max(wHalf,1), kMRImeanBlockSize );

      // Accumulator
      float myVal = 0;

      __shared__ float currPatch[kMRImeanBlockSize][kMRImeanBlockSize];

      // Calculate patch limits (note integer declarations avoid -ve trouble)
      const int yFirst = max( 0, iyStart - patchSize );
      const int yLast = min( src.dims.y - patchSize, iyStart + patchSize );

      for( int yBegin = yFirst; yBegin <= yLast; yBegin += kMRImeanBlockSize ) {
	// Load the patch
	currPatch[ty][tx] = src( ix, yBegin+ty, iz );
	__syncthreads();

	// Accumulate desired values
	for( unsigned int i=0; i<kMRImeanBlockSize; i++ ) {
	  int acty = yBegin + i;

	  if( (acty>=myyMin) && (acty<=myyMax) ) {
	    myVal += currPatch[i][tx];
	  }

	}

	__syncthreads();
      }

      // Save result
      dst(ix,iy,iz) = dst.ConvertFloat( myVal );
    }


    //! Kernel to compute z direction means
    template<typename T>
    __global__ void MRImeanZkernel( const GPU::Classes::MRIframeOnGPU<T> src,
				    GPU::Classes::MRIframeOnGPU<float> dst,
				    const dim3 actualDims,
				    const unsigned int wSize ) {
      /*!
	Kernel to compute means in the x direction, based on
	the given window size.
	Basically, does a 1D convolution, but with different
	boundary conditions to MRIConvolveKernelZ.
	Also, since this is meant to be part of a pipeline,
	the destination type must be float
      */
      const unsigned int bz = blockIdx.x / (src.dims.x/kMRImeanBlockSize);
      const unsigned int bx = blockIdx.x % (src.dims.x/kMRImeanBlockSize);
      const unsigned int tx = threadIdx.x;
      // Note... tz is threadIdx.y
      const unsigned int tz = threadIdx.y;
      
      const int ixStart = bx * kMRImeanBlockSize;
      const int izStart = bz * kMRImeanBlockSize;

      const int ix = ixStart + tx;
      const int iy = blockIdx.y;
      const int iz = izStart + tz;

      const int wHalf = wSize/2;

      // Calculate voxels which will contribute, clamping to edges
      const unsigned int myzMin = max( 0           , iz - wHalf );
      const unsigned int myzMax = min( actualDims.z, iz + wHalf );

      // Again, declare int to remove need for some casts
      const int patchSize = Roundup( max(wHalf,1), kMRImeanBlockSize );

      // Accumulator
      float myVal = 0;

      __shared__ float currPatch[kMRImeanBlockSize][kMRImeanBlockSize];

      // Calculate patch limits (note integer declarations avoid -ve trouble)
      const int zFirst = max( 0, izStart - patchSize );
      const int zLast = min( src.dims.z - patchSize, izStart + patchSize );

      for( int zBegin = zFirst; zBegin <= zLast; zBegin += kMRImeanBlockSize ) {
	// Load the patch
	currPatch[tz][tx] = src( ix, iy, zBegin+tz );
	__syncthreads();

	// Accumulate desired values
	for( unsigned int i=0; i<kMRImeanBlockSize; i++ ) {
	  int actz = zBegin + i;

	  if( (actz>=myzMin) && (actz<=myzMax) ) {
	    myVal += currPatch[i][tx];
	  }

	}

	__syncthreads();
      }

      // Save result
      dst(ix,iy,iz) = dst.ConvertFloat( myVal );
    }


    //! Kernel to normalise means
    template<typename U>
    __global__ void MRImeanNormal( const GPU::Classes::MRIframeOnGPU<float> src,
				   GPU::Classes::MRIframeOnGPU<U> dst,
				   const dim3 actualDims,
				   const unsigned int wSize ) {
      /*!
	Kernel to normalise the means computed by the earlier
	stages.
	As such, the input type must be a float
      */
      const unsigned int by = blockIdx.x / (src.dims.x/kMRImeanBlockSize);
      const unsigned int bx = blockIdx.x % (src.dims.x/kMRImeanBlockSize);
      const unsigned int tx = threadIdx.x;
      const unsigned int ty = threadIdx.y;
      
      const int ixStart = bx * kMRImeanBlockSize;
      const int iyStart = by * kMRImeanBlockSize;

      const int ix = ixStart + tx;
      const int iy = iyStart + ty;
      const int iz = blockIdx.y;

      const int wHalf = wSize/2;

      // Calculate voxels which contributed, clamping to edges
      const unsigned int myxMin = max( 0           , ix - wHalf );
      const unsigned int myxMax = min( actualDims.x, ix + wHalf );
      const unsigned int myyMin = max( 0           , iy - wHalf );
      const unsigned int myyMax = min( actualDims.y, iy + wHalf );
      const unsigned int myzMin = max( 0           , iz - wHalf );
      const unsigned int myzMax = min( actualDims.z, iz + wHalf );


      const unsigned long myVolume = ( myxMax - myxMin + 1 ) *
	(myyMax - myyMin + 1 ) *
	(myzMax - myzMin + 1 );

      dst( ix, iy, iz ) = dst.ConvertFloat( src( ix, iy, iz ) / myVolume );
    }



    //! Runs the mean filtering kernel
    template<typename T, typename U>
    void MRIDoMeanGPU( const GPU::Classes::MRIframeGPU<T> &src,
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

      // We need intermediates which are floats
      GPU::Classes::MRIframeGPU<float> f1, f2;

      // Get correctly sized intermediates
      f1.Allocate( src );
      f2.Allocate( src );

      // Create the GPU kernel objects
      GPU::Classes::MRIframeOnGPU<T> srcGPU(src);
      GPU::Classes::MRIframeOnGPU<float> f1GPU( f1 );
      GPU::Classes::MRIframeOnGPU<float> f2GPU( f2 );
      GPU::Classes::MRIframeOnGPU<U> dstGPU(dst);

      // Do the three convolutions. Recall objects have same dims
      dim3 grid, threads;

      grid.x = (dstDims.x/kMRImeanBlockSize) * (dstDims.y/kMRImeanBlockSize);
      grid.y = dstDims.z;
      grid.z = 1;
      threads.x = threads.y = kMRImeanBlockSize;
      threads.z = 1;

      // Do the X direction
      MRImeanXkernel<<<grid,threads,0,myStream>>>( srcGPU, f1GPU, srcCPUdims, wSize );
      // Do the Y direction
      MRImeanYkernel<<<grid,threads,0,myStream>>>( f1GPU, f2GPU, srcCPUdims, wSize );

      // Slight change for Z direction
      grid.x = (dstDims.x/kMRImeanBlockSize) * (dstDims.z/kMRImeanBlockSize);
      grid.y = dstDims.y;
      MRImeanZkernel<<<grid,threads,0,myStream>>>( f2GPU, f1GPU, srcCPUdims, wSize );

      // Normalise
      grid.x = (dstDims.x/kMRImeanBlockSize) * (dstDims.y/kMRImeanBlockSize);
      grid.y = dstDims.z;
      MRImeanNormal<<<grid,threads,0,myStream>>>( f1GPU, dstGPU, srcCPUdims, wSize );


      CUDA_CHECK_ERROR_ASYNC( "MRImeanKernel failed!" );
    }
    


    // ######################################################

    //! Wrapper class for the MRI mean algorithm
    class MRImean {
    private:
      cudaStream_t stream;
      mutable char* h_workspace;
      mutable size_t workSize;
      
      // =======================
      
      template<typename T>
      void DispatchWrap( const MRI* src, MRI* dst,
			 const unsigned int wSize,
			 const int srcFrame, const int dstFrame ) const {
	switch( dst->type ) {
	case MRI_UCHAR:
	  this->MeanDispatch<T,unsigned char>( src, dst, wSize, srcFrame, dstFrame );
	  break;
	  
	case MRI_SHORT:
	  this->MeanDispatch<T,short>( src, dst, wSize, srcFrame, dstFrame );
	  break;
	  
	case MRI_FLOAT:
	  this->MeanDispatch<T,float>( src, dst, wSize, srcFrame, dstFrame );
	  break;
	  
	default:
	  std::cerr << __FUNCTION__
		    << ": Unrecognised destination MRI type "
		    << dst->type
		    << std::endl;
	  exit( EXIT_FAILURE );
	}
      }
      
      // =========================

      //! Ensures internal pinned memory buffer is at least of size nBytes
      void Allocate( const size_t nBytes ) const {
	if( this->workSize < nBytes ) {
	  this->Release();

	  CUDA_SAFE_CALL( cudaHostAlloc( (void**)&(this->h_workspace),
					 nBytes,
					 cudaHostAllocDefault ) );
	  this->workSize = nBytes;
	}
      }
	  

      //! Releases internal pinned memory buffer
      void Release( void ) const {
	if( h_workspace != NULL ) {
	  CUDA_SAFE_CALL( cudaFreeHost( h_workspace ) );
	  h_workspace = NULL;
	  workSize = 0;
	}
      }


    public:
      //! Default constructor
      MRImean( void ) : stream(0),
			h_workspace(NULL),
			workSize(0) {}
      

      //! Dispatch for data on the CPU of unknown type
      void DoMean( const MRI* src, MRI* dst,
		   const unsigned int wSize,
		   const unsigned int srcFrame = 0,
		   const unsigned int dstFrame = 0 ) const {
	switch( src->type ) {
	case MRI_UCHAR:
	  this->DispatchWrap<unsigned char>( src, dst, wSize, srcFrame, dstFrame );
	  break;
	  
	case MRI_SHORT:
	  this->DispatchWrap<short>( src, dst, wSize, srcFrame, dstFrame );
	  break;
	  
	case MRI_FLOAT:
	  this->DispatchWrap<float>( src, dst, wSize, srcFrame, dstFrame );
	  break;
	  
	default:
	  std::cerr << __FUNCTION__
		    << ": Unrecognised source MRI type "
		    << src->type
		    << std::endl;
	  exit( EXIT_FAILURE );
	}
      }
      

      //! Templated dispatch for known data types
      template<typename T, typename U>
      void MeanDispatch( const MRI* src, MRI* dst,
			 const unsigned int wSize,
			 const int srcFrame, const int dstFrame ) const {
	/*!
	  Templated dispatch routine for MRI mean function on the
	  GPU.
	  Transfers data to the GPU, and retrieves the results
	*/

	GPU::Classes::MRIframeGPU<T> srcGPU;
	GPU::Classes::MRIframeGPU<U> dstGPU;

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
	  this->Allocate( srcWorkSize );
	} else {
	  this->Allocate( dstWorkSize );
	}

	// Send the source data
	srcGPU.Send( src, srcFrame, this->h_workspace, this->stream );

	// Run the filter
	MRIDoMeanGPU( srcGPU, dstGPU, wSize, this->stream );

	// Get the results
	dstGPU.Recv( dst, dstFrame, this->h_workspace, this->stream );

	CUDA_CHECK_ERROR( "Mean filtering failure" );
      }

    };

  }
}


static GPU::Algorithms::MRImean myMean;


MRI* MRImean_cuda( const MRI* src, MRI* dst,
		   const unsigned int wSize ) {
  /*!
    Wrapper around the GPU routine, to be called from the
    original MRImean routine.
    Note that the frames default to zero, per the original
    MRImean routine.
  */

  myMean.DoMean( src, dst, wSize );

  return( dst );
}
