/**
 * @file  mriconvolve_cuda.cu
 * @brief Holds MRI convolution functions for the GPU
 *
 * Implements MRI convolutions on the GPU. These routines will be hooked
 * into the higher level routines in the main library.
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/02/02 15:09:58 $
 *    $Revision: 1.8 $
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


#include "mriconvolve_cuda.h"


//#define SHOW_TIMINGS


/*
  I'd really like the following to be within the namespace, but
  I can't charm nvcc into allowing this. We obviously
  need some more C++ available for CUDA. And probably a real
  linker
*/
    
//! Device constant indicating the number of values in the convolution kernel
__constant__ unsigned int dc_mriconv1d_kernel_nVals;

//! Texture reference to the convolution kernel
texture<float, 1, cudaReadModeElementType> dtl_mriconv1d_kernel;


namespace GPU {
  namespace Algorithms {

    const unsigned int kConv1dBlockSize = 16;
    
    // Some timers
    SciGPU::Utilities::Chronometer tMRIconv1dMem, tMRIconv1dMemHost;
    SciGPU::Utilities::Chronometer tMRIconv1dSend, tMRIconv1dRecv;
    SciGPU::Utilities::Chronometer tMRIconv1dCompute;
    SciGPU::Utilities::Chronometer tMRIconv1dTotal;
    
    // =================================================
    
    //! Array to contain the convolution kernel on the device
    static float* d_mriconv1d_kernel;

    // =================================================
    
    //! Prepares convolution kernel for use
    void MRIconv1d_SendKernel( const float* kernel,
			       const unsigned int nVals ) {
      /*!
	This routine is responsible for preparing the dtl_mriconv1d_kernel
	texture for use.
	The array on the device is padded with an extra zero on each end,
	to save us some explicit boundary condition checks (the texture
	units will handle them).
	@param[in] kernel Array containing the kernel values
	@param[in] nVals The number of values in the kernel
      */

      // Allocate and zero GPU memory
      CUDA_SAFE_CALL( cudaMalloc( (void**)&d_mriconv1d_kernel,
				  (2+nVals)*sizeof(float) ) );
      CUDA_SAFE_CALL( cudaMemset( d_mriconv1d_kernel,
				  0,
				  (2+nVals)*sizeof(float) ) );

      // Copy the convolution kernel to the GPU
      // Note the extra offset
      CUDA_SAFE_CALL( cudaMemcpy( &(d_mriconv1d_kernel[1]),
				  kernel,
				  nVals*sizeof(float),
				  cudaMemcpyHostToDevice ) );

      // Copy the size of the texture to device constant memory
      CUDA_SAFE_CALL( cudaMemcpyToSymbol( "dc_mriconv1d_kernel_nVals",
					  &nVals,
					  sizeof(unsigned int) ) );

      // ------------------
      // Set up the texture
      
      cudaChannelFormatDesc cd_kernel = cudaCreateChannelDesc<float>();
      
      // Describe the addressing modes
      dtl_mriconv1d_kernel.normalized = false;
      dtl_mriconv1d_kernel.addressMode[0] = cudaAddressModeClamp;
      dtl_mriconv1d_kernel.filterMode = cudaFilterModePoint;
      
      // Bind the texture together
      CUDA_SAFE_CALL( cudaBindTexture( 0,
				       dtl_mriconv1d_kernel,
				       d_mriconv1d_kernel,
				       cd_kernel,
				       (2+nVals)*sizeof(float) ) );
    }
    
    
    // ----------------

    //! Releases convolution kernel after use
    void MRIconv1d_ReleaseKernel( void ) {
      /*!
	Releases everything on the device associated with the convolution kernel
      */
      CUDA_SAFE_CALL( cudaUnbindTexture( dtl_mriconv1d_kernel ) );
      CUDA_SAFE_CALL( cudaFree( d_mriconv1d_kernel ) );
      d_mriconv1d_kernel = NULL;
    }


    // ----------------

    //! Device function to perform convolution kernel look ups
    __device__ float GetMRIconvkernel( const int i ) {
      /*!
	Performs the texture look up of the convolution kernel.
	The texture will be cached, and the extra zeros padded
	mean we don't need to worry about look ups which are
	out of range.
      */
      
      // Note the additional offset for the zero padding
      return( tex1Dfetch( dtl_mriconv1d_kernel, i+1.5f ) );
    }


    // ==================================================
    
    //! Function to round up the convolution kernel size
    __device__ unsigned int MRIRoundUpConvolutionKernelSize( void ) {
      /*!
	To avoid edge cases, we want to make everything into
	units of kConv1dBlockSize.
	This performs the task for the convolution kernel
	Note that we have an extra divide by two since we're only
	really interested in the half-width of the convolution kernel
      */
      float temp;

      temp = static_cast<float>(dc_mriconv1d_kernel_nVals/2) / kConv1dBlockSize ;

      unsigned int n = static_cast<unsigned int>(ceilf( temp ) );

      return( n * kConv1dBlockSize );
    }




    //! Kernel to convolve in the X direction
    template<typename T, typename U>
    __global__ void MRIConvolveKernelX( const GPU::Classes::MRIframeOnGPU<T> src,
					GPU::Classes::MRIframeOnGPU<U> dst ) {
      /*!
	Kernel to do a convolution in the x direction, using the
	convolution kernel set up in the texture.
	Each block will do one 16x16 (x,y) patch of the MRI.
	The z co-ordinate will be in blockIdx.y
	The (x,y) location of the patch will be derived from blockIdx.x.
	Since the slices array is padded out to multiples of 16 by
	the other MRI slices routines (the host routine checks this), we
	don't have to worry about nasty edge cases
      */
      // Extract some co-ordinates
      const unsigned int by = blockIdx.x / (src.dims.x/kConv1dBlockSize);
      const unsigned int bx = blockIdx.x % (src.dims.x/kConv1dBlockSize);
      const unsigned int iz = blockIdx.y;
      const unsigned int tx = threadIdx.x;
      const unsigned int ty = threadIdx.y;
      
      const unsigned int ixStart = bx * kConv1dBlockSize;
      const unsigned int iyStart = by * kConv1dBlockSize;
      const unsigned int myx = ixStart + tx;
      const unsigned int myy = iyStart + ty;
      // Round up the convolution kernel size
      const unsigned int iConvKernelSize = MRIRoundUpConvolutionKernelSize();
      
      // Accumulator value
      float myVoxel = 0;
      
      // Declared as float because the GPU likes 4 byte datatypes
      __shared__ float srcCache[kConv1dBlockSize][kConv1dBlockSize];
      
      // Loop in the X direction, loading sub-blocks of the src array into the cache
      // Note that xBegin can be negative, so we must cast ixStart and iConvKernelSize
      // to signed integers
      int xFirst = static_cast<int>(ixStart) - static_cast<int>(iConvKernelSize);
      int xLast = ixStart + iConvKernelSize;
      
      for( int xBegin = xFirst; xBegin <= xLast; xBegin += kConv1dBlockSize ) {
	// Load the cache
	srcCache[ty][tx] = src( src.ClampCoord(xBegin+tx,src.dims.x), myy, iz );
	__syncthreads();
	
	// Accumulate
	for( unsigned int i=0; i<kConv1dBlockSize; i++ ) {
	  int convLoc;
	  convLoc = (xBegin-ixStart) - static_cast<int>(tx);
	  convLoc += i;
	  convLoc += (dc_mriconv1d_kernel_nVals/2);
	  myVoxel += srcCache[ty][i] * GetMRIconvkernel( convLoc );
	}
	__syncthreads();
      }
      
      dst(myx, myy, iz) = dst.ConvertFloat( myVoxel );
      
    }
    


    //! Kernel to convolve in the Y direction
    template<typename T, typename U>
    __global__ void MRIConvolveKernelY( const GPU::Classes::MRIframeOnGPU<T> src,
				     GPU::Classes::MRIframeOnGPU<U> dst ) {
      /*!
	Kernel to do a convolution in the y direction, using the
	convolution kernel set up in the texture.
	Each block will do one 16x16 (x,y) patch of the MRI.
	The z co-ordinate will be in blockIdx.y
	The (x,y) location of the patch will be derived from blockIdx.x.
	Since the slices array is padded out to multiples of 16 by
	the other MRI slices routines (the host routine checks this), we
	don't have to worry about nasty edge cases
      */
      // Extract some co-ordinates
      const unsigned int by = blockIdx.x / (src.dims.x/kConv1dBlockSize);
      const unsigned int bx = blockIdx.x % (src.dims.x/kConv1dBlockSize);
      const unsigned int iz = blockIdx.y;
      const unsigned int tx = threadIdx.x;
      const unsigned int ty = threadIdx.y;
      
      const unsigned int ixStart = bx * kConv1dBlockSize;
      const unsigned int iyStart = by * kConv1dBlockSize;
      const unsigned int myx = ixStart + tx;
      const unsigned int myy = iyStart + ty;
      // Round up the convolution kernel size
      const unsigned int iConvKernelSize = MRIRoundUpConvolutionKernelSize();
      
      // Accumulator value
      float myVoxel = 0;
      
      // Declared as float because the GPU likes 4 byte datatypes
      __shared__ float srcCache[kConv1dBlockSize][kConv1dBlockSize];
      
      // Loop in the Y direction, loading sub-blocks of the src array into the cache
      // Note that yBegin can be negative, so we must cast iyStart and iConvKernelSize
      // to signed integers
      int yFirst = static_cast<int>(iyStart) - static_cast<int>(iConvKernelSize);
      int yLast = iyStart + iConvKernelSize;
      
      for( int yBegin = yFirst; yBegin <= yLast; yBegin += kConv1dBlockSize ) {
	// Load the cache
	srcCache[ty][tx] = src( myx, src.ClampCoord(yBegin+ty,src.dims.y), iz );
	__syncthreads();
	
	// Accumulate
	for( unsigned int i=0; i<kConv1dBlockSize; i++ ) {
	  int convLoc;
	  convLoc = (yBegin-iyStart) - static_cast<int>(ty);
	  convLoc += i;
	  convLoc += (dc_mriconv1d_kernel_nVals/2);
	  myVoxel += srcCache[i][tx] * GetMRIconvkernel( convLoc );
	}
	__syncthreads();
      }
      
      dst(myx, myy, iz) = dst.ConvertFloat( myVoxel );
      
    }
    



    //! Kernel to convolve in the Z direction
    template<typename T, typename U>
    __global__ void MRIConvolveKernelZ( const GPU::Classes::MRIframeOnGPU<T> src,
					GPU::Classes::MRIframeOnGPU<U> dst ) {
      /*!
	Kernel to do a convolution in the z direction, using the
	convolution kernel set up in the texture.
	Each block will do one 16x16 (x,z) patch of the MRI.
	The y co-ordinate will be in blockIdx.y
	The (x,z) location of the patch will be derived from blockIdx.x.
	Since the slices array is padded out to multiples of 16 by
	the other MRI slices routines (the host routine checks this), we
	don't have to worry about nasty edge cases
      */
      // Extract some co-ordinates
      const unsigned int bz = blockIdx.x / (src.dims.x/kConv1dBlockSize);
      const unsigned int bx = blockIdx.x % (src.dims.x/kConv1dBlockSize);
      const unsigned int iy = blockIdx.y;
      const unsigned int tx = threadIdx.x;
      // Note that we assign y thread index to tz, for naming ease
      const unsigned int tz = threadIdx.y;
      
      const unsigned int ixStart = bx * kConv1dBlockSize;
      const unsigned int izStart = bz * kConv1dBlockSize;
      const unsigned int myx = ixStart + tx;
      const unsigned int myz = izStart + tz;
      // Round up the convolution kernel size
      const unsigned int iConvKernelSize = MRIRoundUpConvolutionKernelSize();
      
      // Accumulator value
      float myVoxel = 0;

      // Declared as float because the GPU likes 4 byte datatypes
      __shared__ float srcCache[kConv1dBlockSize][kConv1dBlockSize];
      
      // Loop in the z direction, loading sub-blocks of the src array into the cache
      // Note that zBegin can be negative, so we must cast izStart and iConvKernelSize
      // to signed integers
      int zFirst = static_cast<int>(izStart) - static_cast<int>(iConvKernelSize);
      int zLast = izStart + iConvKernelSize;
      
      for( int zBegin = zFirst; zBegin <= zLast; zBegin += kConv1dBlockSize ) {
	// Load the cache
	srcCache[tz][tx] = src( myx, iy, src.ClampCoord(zBegin+tz,src.dims.z) );
	__syncthreads();
	
	// Accumulate
	for( unsigned int i=0; i<kConv1dBlockSize; i++ ) {
	  int convLoc;
	  convLoc = (zBegin-izStart) - static_cast<int>(tz);
	  convLoc += i;
	  convLoc += (dc_mriconv1d_kernel_nVals/2);
	  myVoxel += srcCache[i][tx] * GetMRIconvkernel( convLoc );
	}
	__syncthreads();
      }
      
      dst(myx, iy, myz) = dst.ConvertFloat( myVoxel );
  
    }






    //! Runs the 1D convolution kernel on the GPU
    template<typename T, typename U>
    void MRIConvolve1dGPU( const GPU::Classes::MRIframeGPU<T> &src,
			   GPU::Classes::MRIframeGPU<U> &dst,
			   const int axis,  
			   const cudaStream_t myStream = 0 ) {
      
      /*!
	Runs the 1D convolution kernel on the GPU.
	Prior to calling this routine, MRIconv1d_SendKernel must
	be called to set up the convolution kernel
	@param[in] d_src The set of source frames
	@param[out] d_dst The set of destination frames
	@param[in] axis Which axis to do the convolution along
	@param[in] srcDims Dimensions of the source frames
	@param[in] dstDims Dimensions of the destinaiton frames
	@param[in] srcFrame Which frame of the source to use
	@param[in] dstFrame Which frame of the destination to use
	@param[in] myStream CUDA stream which should be used (Defaults to stream 0)
      */
      

      const dim3 srcDims = src.GetGPUDims();
      const dim3 dstDims = dst.GetGPUDims();

      // Check dimensions
      if( (srcDims.x != dstDims.x) &&
	  (srcDims.y != dstDims.y) &&
	  (srcDims.z != dstDims.z) ) {
	std::cerr << __FUNCTION__ << ": Dimension mismatch" << std::endl;
	exit( EXIT_FAILURE );
      }
      
      // Check padding (only need to do one, given above check)
      if( !dst.CheckPadding( kConv1dBlockSize ) ) {
	std::cerr << __FUNCTION__ <<
	  ": Arrays on GPU must be padded to multiples of kConv1dBlockSize" <<
	  std::endl;
	exit( EXIT_FAILURE );
      }
      
      
      dim3 grid, threads;
      GPU::Classes::MRIframeOnGPU<T> srcGPU(src);
      GPU::Classes::MRIframeOnGPU<U> dstGPU(dst);
      
      //cudaPrintfInit( 1024L*1024L*512 );
      
      switch( axis ) {
      case MRI_WIDTH:
	grid.x = (srcDims.x/kConv1dBlockSize) * (srcDims.y/kConv1dBlockSize);
	grid.y = srcDims.z;
	grid.z = 1;
	threads.x = threads.y = kConv1dBlockSize;
	threads.z = 1;
	
	MRIConvolveKernelX<T,U><<<grid,threads,0,myStream>>>( srcGPU, dstGPU );
	CUDA_CHECK_ERROR_ASYNC( "MRIconvolveKernelX failed!" );
	break;
	
      case MRI_HEIGHT:
	grid.x = (srcDims.x/kConv1dBlockSize) * (srcDims.y/kConv1dBlockSize);
	grid.y = srcDims.z;
	grid.z = 1;
	threads.x = threads.y = kConv1dBlockSize;
	threads.z = 1;
	
	MRIConvolveKernelY<T,U><<<grid,threads,0,myStream>>>( srcGPU, dstGPU );
	CUDA_CHECK_ERROR_ASYNC( "MRIconvolveKernelY failed!" );
	break;
	
      case MRI_DEPTH:
	// Slight change, since we do (x,z) patches
	grid.x = (srcDims.x/kConv1dBlockSize) * (srcDims.z/kConv1dBlockSize);
	grid.y = srcDims.y;
	threads.x = threads.y = kConv1dBlockSize;
	threads.z = 1;
	
	MRIConvolveKernelZ<T,U><<<grid,threads,0,myStream>>>( srcGPU, dstGPU );
	CUDA_CHECK_ERROR_ASYNC( "MRIconvolveKernelZ failed!" );
	break;
	
      default:
	std::cerr << __FUNCTION__ << ": Incompatible universe detected." << std::endl;
	std::cerr << "GPU functions are only tested ";
	std::cerr << "in a universe with three spatial dimensions" << std::endl;
	std::cerr << "Please adjust your reality accordingly, ";
	std::cerr << "and try again" << std::endl;
	std::cerr << "MRIframeGPU version was:\n" << src.VersionString() << std::endl;
	exit( EXIT_FAILURE );
      }
      
      //cudaPrintfDisplay( stdout, true );
      
      //cudaPrintfEnd();
    }


    
    //! Dispatch routine with transfers
    template<typename T, typename U>
    void MRIConv1dDispatch( const MRI* src, MRI* dst,
			    const int axis,
			    const int srcFrame, const int dstFrame ) {
      /*!
	This is a dispatch routine for the 1D convolution on the GPU.
	It transfers the data to the GPU, runs the convolution, and retrieves
	the results
	Things are written this way to avoid nastily nested switch statements.
      */

      tMRIconv1dTotal.Start();

      GPU::Classes::MRIframeGPU<T> srcGPU;
      GPU::Classes::MRIframeGPU<U> dstGPU;
      
      char* h_workspace;
      size_t srcWorkSize, dstWorkSize;
      
      // Allocate the GPU arrays
      tMRIconv1dMem.Start();
      srcGPU.Allocate( src, kConv1dBlockSize );
      dstGPU.Allocate( dst, kConv1dBlockSize );
      tMRIconv1dMem.Stop();
      
      // Put in some sanity checks
      srcGPU.VerifyMRI( src );
      dstGPU.VerifyMRI( dst );
      
      // Allocate the workspace array
      tMRIconv1dMemHost.Start();
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
      tMRIconv1dMemHost.Stop();

      // Send the source data
      tMRIconv1dSend.Start();
      srcGPU.Send( src, srcFrame, h_workspace );
      tMRIconv1dSend.Stop();
      
      // Run the convolution
      tMRIconv1dCompute.Start();
      MRIConvolve1dGPU( srcGPU, dstGPU, axis );
      tMRIconv1dCompute.Stop();
  
      // Retrieve the answers
      tMRIconv1dRecv.Start();
      dstGPU.Recv( dst, dstFrame, h_workspace );
      tMRIconv1dRecv.Stop();
      
      tMRIconv1dMemHost.Start();
      CUDA_SAFE_CALL( cudaFreeHost( h_workspace ) );
      tMRIconv1dMemHost.Stop();
      
      CUDA_CHECK_ERROR( "1D Convolution failure" );
      
      tMRIconv1dTotal.Stop();
      
      // No need to release - destructor will do that automatically
    }


    //! Wrapper for Gaussian convolution
    template<typename T>
    void MRIConvGaussianDispatch( const MRI* src, MRI* dst ) {
      /*!
	Function to run the 3D gaussian convolution on the GPU
	without pulling intermediate results back to the host.
	Assumes that src and dst types are the same.
	Assumes that the texture is already set up on the GPU
      */

      GPU::Classes::MRIframeGPU<T> frame1, frame2;

      char* h_workspace;
      size_t workSize;

      // Do some allocation
      frame1.Allocate( src, kConv1dBlockSize );
      frame2.Allocate( src, kConv1dBlockSize );

      // Verify (note frame2 verified from dst)
      frame1.VerifyMRI( src );
      frame2.VerifyMRI( dst );

      // Allocate workspace
      workSize = frame1.GetBufferSize();
      CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_workspace,
				     workSize,
				     cudaHostAllocDefault ) );

      // Loop over frame
      for( unsigned int iFrame=0; iFrame < src->nframes; iFrame++ ) {
	frame1.Send( src, iFrame, h_workspace );

	MRIConvolve1dGPU( frame1, frame2, MRI_WIDTH );
	MRIConvolve1dGPU( frame2, frame1, MRI_HEIGHT );
	MRIConvolve1dGPU( frame1, frame2, MRI_DEPTH );

	frame2.Recv( dst, iFrame, h_workspace );
      }

      CUDA_SAFE_CALL( cudaFreeHost( h_workspace ) );

      CUDA_CHECK_ERROR( "Gaussian convolution failure" );

    }

    //! Dispatch wrapper
    template<typename T>
    void MRIConv1dDispatchWrap( const MRI* src, MRI* dst,
				const int axis,
				const int srcFrame, const int dstFrame ) {
      switch( dst->type ) {
      case MRI_UCHAR:
	MRIConv1dDispatch<T,unsigned char>( src, dst, axis, srcFrame, dstFrame );
	break;
	
      case MRI_SHORT:
	MRIConv1dDispatch<T,short>( src, dst, axis, srcFrame, dstFrame );
	break;
	
      case MRI_FLOAT:
	MRIConv1dDispatch<T,float>( src, dst, axis, srcFrame, dstFrame );
	break;

      default:
	std::cerr << __FUNCTION__ << ": Unrecognised destination MRI type " << dst->type << std::endl;
	exit( EXIT_FAILURE );
      }
    }
    
  }
}

// =================================================


MRI* MRIconvolve1d_cuda( const MRI* src, MRI* dst,
			 const float *kernel,
			 const unsigned int kernelLength,
			 const int axis,
			 const int srcFrame, const int dstFrame ) {
  /*!
    Reimplementation of MRIconvolve1d for the GPU.
    This is the 'drop in' replacement, and so has to do a lot of
    data transfers.
    As such, I don't expect it to be fast
  */


  // Get the convolution kernel to the GPU
  GPU::Algorithms::MRIconv1d_SendKernel( kernel, kernelLength );

  switch( src->type ) {
  case MRI_UCHAR:
    GPU::Algorithms::MRIConv1dDispatchWrap<unsigned char>( src, dst, axis, srcFrame, dstFrame );
    break;

  case MRI_SHORT:
    GPU::Algorithms::MRIConv1dDispatchWrap<short>( src, dst, axis, srcFrame, dstFrame );
    break;

  case MRI_FLOAT:
    GPU::Algorithms::MRIConv1dDispatchWrap<float>( src, dst, axis, srcFrame, dstFrame );
    break;
    
  default:
    std::cerr << __FUNCTION__ << ": Unrecognised source MRI type " << src->type << std::endl;
    exit( EXIT_FAILURE );
  }
  
  // Release the convolution kernel
  GPU::Algorithms::MRIconv1d_ReleaseKernel();


  return( dst );
}

// =================================================

MRI* MRIconvolveGaussian_cuda( const MRI* src, MRI* dst,
			       const float *kernel,
			       const unsigned int kernelLength ) {
  /*!
    Implementation of MRIconvolveGaussian for the GPU.
    Designed to be called form that routine
  */

  // Send the convolution kernel
  GPU::Algorithms::MRIconv1d_SendKernel( kernel, kernelLength );

  switch( src->type ) {
  case MRI_UCHAR:
    GPU::Algorithms::MRIConvGaussianDispatch<unsigned char>( src, dst );
    break;

  default:
    std::cerr << __FUNCTION__ << ": Unrecognised source MRI type " << src->type << std::endl;
    exit( EXIT_FAILURE );
  }


  // Release the kernel
  GPU::Algorithms::MRIconv1d_ReleaseKernel();


  return( dst );
}




// ======================================================

//! Stream insertion operator for timer
static std::ostream& operator<<( std::ostream& os,
				 const SciGPU::Utilities::Chronometer& timer ) {
  
  os << std::setw(9) << std::setprecision(6) << timer.GetAverageTime() << " ms (avg) ";
  os << std::setw(9) << std::setprecision(6) << timer.GetTime() << " ms (tot)";

  return( os );
}

void MRIconvShowTimers( void ) {
  /*!
    Pretty prints timers to std.out
  */

  std::cout << "=============================================" << std::endl;
  std::cout << "GPU convolution timers" << std::endl;
  std::cout << "----------------------" << std::endl;
#ifndef CUDA_FORCE_SYNC
  std::cout << "WARNING: CUDA_FORCE_SYNC not #defined" << std::endl;
  std::cout << "Timings may not be accurate" << std::endl;
#endif
  std::cout << std::endl;

  std::cout << "MRIConv1dDispatch" << std::endl;
  std::cout << "Host Memory   : " << GPU::Algorithms::tMRIconv1dMemHost << std::endl;
  std::cout << "GPU Memory    : " << GPU::Algorithms::tMRIconv1dMem << std::endl;
  std::cout << "Send          : " << GPU::Algorithms::tMRIconv1dSend << std::endl;
  std::cout << "Compute       : " << GPU::Algorithms::tMRIconv1dCompute << std::endl;
  std::cout << "Receive       : " << GPU::Algorithms::tMRIconv1dRecv << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Total : " << GPU::Algorithms::tMRIconv1dTotal << std::endl;
  std::cout << std::endl;

  std::cout << "=============================================" << std::endl;
}
