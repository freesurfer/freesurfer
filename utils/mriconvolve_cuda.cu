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
 *    $Author: nicks $
 *    $Date: 2012/12/15 21:47:47 $
 *    $Revision: 1.33 $
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

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>


#include "chronometer.hpp"
#include "cudacheck.h"


#include "mriframegpu.hpp"

#include "mriconvolve_cuda.hpp"

#include "mriconvolve_cuda.h"


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


namespace GPU
{
//! Namespace to hold algorithms which operate on GPU data
namespace Algorithms
{

//! Break frames into tiles of this size
const unsigned int kConv1dBlockSize = 16;

// =================================================


//! Device function to perform convolution kernel look ups
__device__ float GetMRIconvkernel( const int i )
{
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
__device__ unsigned int MRIRoundUpConvolutionKernelSize( void )
{
  /*!
  To avoid edge cases, we want to make everything into
  units of \c kConv1dBlockSize.
  This performs the task for the convolution kernel
  Note that we have an extra divide by two since we're only
  really interested in the half-width of the convolution kernel
  */


  unsigned int n = NextMultiple( dc_mriconv1d_kernel_nVals/2,
                                 kConv1dBlockSize );

  return( n );
}




//! Kernel to convolve in the X direction
template<typename T, typename U>
__global__ void MRIConvolveKernelX( const GPU::Classes::MRIframeOnGPU<T> src,
                                    GPU::Classes::MRIframeOnGPU<U> dst,
                                    const dim3 coverGrid )
{
  /*!
  Kernel to do a convolution in the x direction, using the
  convolution kernel set up in the texture.
  Each block will do one 16x16 (x,y) patch of the MRI.
  The z co-ordinate will be in \c blockIdx.y.
  The (x,y) location of the patch will be derived from \c blockIdx.x
  and the \a coverGrid argument.
  Edge cases are handled through boundary checks.
  */
  // Extract some co-ordinates
  const unsigned int by = blockIdx.x / coverGrid.x;
  const unsigned int bx = blockIdx.x % coverGrid.x;
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

  for( int xBegin = xFirst; xBegin <= xLast; xBegin += kConv1dBlockSize )
  {
    // Load the cache
    srcCache[ty][tx] = src( src.ClampCoord(xBegin+tx,src.dims.x), myy, iz );
    __syncthreads();

    // Accumulate
    for( unsigned int i=0; i<kConv1dBlockSize; i++ )
    {
      int convLoc;
      convLoc = (xBegin-ixStart) - static_cast<int>(tx);
      convLoc += i;
      convLoc += (dc_mriconv1d_kernel_nVals/2);
      myVoxel += srcCache[ty][i] * GetMRIconvkernel( convLoc );
    }
    __syncthreads();
  }

  if( dst.InVolume( myx, myy, iz ) )
  {
    dst(myx, myy, iz) = dst.ConvertFloat( myVoxel );
  }

}



//! Kernel to convolve in the Y direction
template<typename T, typename U>
__global__ void MRIConvolveKernelY( const GPU::Classes::MRIframeOnGPU<T> src,
                                    GPU::Classes::MRIframeOnGPU<U> dst,
                                    const dim3 coverGrid )
{
  /*!
  Kernel to do a convolution in the y direction, using the
  convolution kernel set up in the texture.
  Each block will do one 16x16 (x,y) patch of the MRI.
  The z co-ordinate will be in \c blockIdx.y.
  The (x,y) location of the patch will be derived from \c blockIdx.x
  and \a coverGrid.
  */
  // Extract some co-ordinates
  const unsigned int by = blockIdx.x / coverGrid.x;
  const unsigned int bx = blockIdx.x % coverGrid.x;
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

  for( int yBegin = yFirst; yBegin <= yLast; yBegin += kConv1dBlockSize )
  {
    // Load the cache
    srcCache[ty][tx] = src( myx, src.ClampCoord(yBegin+ty,src.dims.y), iz );
    __syncthreads();

    // Accumulate
    for( unsigned int i=0; i<kConv1dBlockSize; i++ )
    {
      int convLoc;
      convLoc = (yBegin-iyStart) - static_cast<int>(ty);
      convLoc += i;
      convLoc += (dc_mriconv1d_kernel_nVals/2);
      myVoxel += srcCache[i][tx] * GetMRIconvkernel( convLoc );
    }
    __syncthreads();
  }

  if( dst.InVolume( myx, myy, iz ) )
  {
    dst(myx, myy, iz) = dst.ConvertFloat( myVoxel );
  }

}




//! Kernel to convolve in the Z direction
template<typename T, typename U>
__global__ void MRIConvolveKernelZ( const GPU::Classes::MRIframeOnGPU<T> src,
                                    GPU::Classes::MRIframeOnGPU<U> dst,
                                    const dim3 coverGrid )
{
  /*!
  Kernel to do a convolution in the z direction, using the
  convolution kernel set up in the texture.
  Each block will do one 16x16 (x,z) patch of the MRI.
  The y co-ordinate will be in \c blockIdx.y
  The (x,z) location of the patch will be derived from \c blockIdx.x
  and \a coverGrid.
  */
  // Extract some co-ordinates
  const unsigned int bz = blockIdx.x / coverGrid.x;
  const unsigned int bx = blockIdx.x % coverGrid.x;
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

  for( int zBegin = zFirst; zBegin <= zLast; zBegin += kConv1dBlockSize )
  {
    // Load the cache
    srcCache[tz][tx] = src( myx, iy, src.ClampCoord(zBegin+tz,src.dims.z) );
    __syncthreads();

    // Accumulate
    for( unsigned int i=0; i<kConv1dBlockSize; i++ )
    {
      int convLoc;
      convLoc = (zBegin-izStart) - static_cast<int>(tz);
      convLoc += i;
      convLoc += (dc_mriconv1d_kernel_nVals/2);
      myVoxel += srcCache[i][tx] * GetMRIconvkernel( convLoc );
    }
    __syncthreads();
  }

  if( dst.InVolume( myx, iy, myz ) )
  {
    dst(myx, iy, myz) = dst.ConvertFloat( myVoxel );
  }

}





MRIconvolve::~MRIconvolve( void )
{
  this->ReleaseKernel();
  this->ReleaseWorkspace();
}


void MRIconvolve::ShowTimings( void )
{

#ifdef CUDA_SHOW_TIMINGS
  std::cout << "==================================" << std::endl;
  std::cout << "GPU convolution timers" << std::endl;
  std::cout << "----------------------" << std::endl;
#ifndef CUDA_FORCE_SYNC
  std::cout << "WARNING: CUDA_FORCE_SYNC not #defined" << std::endl;
  std::cout << "Timings may not be accurate" << std::endl;
#endif
  std::cout << std::endl;

  std::cout << "General:" << std::endl;
  std::cout << "Host Memory   : " << MRIconvolve::tHostMem << std::endl;
  std::cout << "GPU Memory    : " << MRIconvolve::tMem << std::endl;
  std::cout << "Send          : " << MRIconvolve::tSend << std::endl;
  std::cout << "Receive       : " << MRIconvolve::tRecv << std::endl;
  std::cout << "Conv. kernel  : " << MRIconvolve::tTextureKernel << std::endl;
  std::cout << std::endl;

  std::cout << "Convolve 1D:" << std::endl;
  std::cout << "      Compute : " << MRIconvolve::tCompute1d << std::endl;
  std::cout << "Total         : " << MRIconvolve::tTotal1d << std::endl;
  std::cout << std::endl;

  std::cout << "Convolve 3D:" << std::endl;
  std::cout << "      Compute : " << MRIconvolve::tCompute3d << std::endl;
  std::cout << "Total         : " << MRIconvolve::tTotal3d << std::endl;
  std::cout << std::endl;


  std::cout << "==================================" << std::endl;
#endif

}

// ---------------------------------------

void MRIconvolve::Convolve1D( const MRI* src, MRI* dst,
                              const int srcFrame, const int dstFrame,
                              const float *kernel,
                              const unsigned int kernelLength,
                              const int axis )
{
  /*!
  If the data types of the source and destination MRI
  structures is not known, then this dispatch routine
  should be used for 1D convolutions.
  It puts the convolution kernel on the device,
  and then handles the templated dispatching of the
  appropriate convolution call.
  Note that the destination MRI must already be allocated.
  If called from the main Freesurfer routine, this will
  be the case.
  @param[in] src The source MRI
  @param[out] dst The destination MRI
  @param[in] srcFrame The frame of the source to be convolved
  @param[in] dstFrame The frame of the destination in which the result should be stored
  @param[in] kernel The convolution kernel to use
  @param[in] kernelLength The size of the convolution kernel
  @param[in] axis The axis along which the convolution should be made
  */
  MRIconvolve::tTotal1d.Start();

  this->BindKernel( kernel, kernelLength );

  switch( src->type )
  {
  case MRI_UCHAR:
    this->DispatchWrap1D<unsigned char>( src, dst,
                                         srcFrame, dstFrame,
                                         axis );
    break;

  case MRI_SHORT:
    this->DispatchWrap1D<short>( src, dst,
                                 srcFrame, dstFrame,
                                 axis );
    break;

  case MRI_FLOAT:
    this->DispatchWrap1D<float>( src, dst,
                                 srcFrame, dstFrame,
                                 axis );
    break;

  default:
    std::cerr << __FUNCTION__
              <<": Unrecognised source MRI type " << src->type
              << std::endl;
    exit( EXIT_FAILURE );
  }

  this->UnbindKernel();

  MRIconvolve::tTotal1d.Stop();
}

// --


template<typename T, typename U>
void MRIconvolve::Dispatch1D( const MRI* src, MRI* dst,
                              const int srcFrame, const int dstFrame,
                              const int axis ) const
{
  /*!
  This is a dispatch routine for the 1D convolution on
  the GPU.
  It transfers the data to the GPU, runs the convolution,
  and retrieves the results
  Things are written this way to avoid nastily nested
  switch statements.
  It assumes that the convolution kernel texture is already
  bound
  */

  GPU::Classes::MRIframeGPU<T> srcGPU;
  GPU::Classes::MRIframeGPU<U> dstGPU;

  size_t srcWorkSize, dstWorkSize;

  // Allocate the GPU arrays
  MRIconvolve::tMem.Start();
  srcGPU.Allocate( src );
  dstGPU.Allocate( dst );
  MRIconvolve::tMem.Stop();

  // Put in some sanity checks
  srcGPU.VerifyMRI( src );
  dstGPU.VerifyMRI( dst );

  // Allocate the workspace array
  srcWorkSize = srcGPU.BufferSize();
  dstWorkSize = dstGPU.BufferSize();

  if( srcWorkSize > dstWorkSize )
  {
    this->AllocateWorkspace( srcWorkSize );
  }
  else
  {
    this->AllocateWorkspace( dstWorkSize );
  }

  // Send the source data
  MRIconvolve::tSend.Start();
  srcGPU.SendFrame( src, srcFrame, this->h_workspace, this->stream );
  MRIconvolve::tSend.Stop();

  // Run the convolution
  MRIconvolve::tCompute1d.Start();
  this->RunGPU1D( srcGPU, dstGPU, axis );
  MRIconvolve::tCompute1d.Stop();

  // Retrieve the answers
  MRIconvolve::tRecv.Start();
  dstGPU.RecvFrame( dst, dstFrame, this->h_workspace, this->stream );
  MRIconvolve::tRecv.Stop();

  CUDA_CHECK_ERROR_ASYNC( "1D Convolution failure" );

  // No need to release - destructor will do that automatically
}

// --


template<typename T, typename U>
void MRIconvolve::RunGPU1D( const GPU::Classes::MRIframeGPU<T> &src,
                            GPU::Classes::MRIframeGPU<U> &dst,
                            const int axis ) const
{

  /*!
  Runs the 1D convolution kernel on the GPU.
  Prior to calling this routine, convolution
  kernel must be set up on the device, and
  the source MRI transferred to the GPU.
  */


  const dim3 srcDims = src.GetDims();
  const dim3 dstDims = dst.GetDims();

  // Check dimensions
  if( srcDims != dstDims )
  {
    std::cerr << __FUNCTION__
              << ": Dimension mismatch" << std::endl;
    exit( EXIT_FAILURE );
  }


  dim3 grid, threads;


  threads.x = threads.y = kConv1dBlockSize;
  threads.z = 1;

  const dim3 coverGrid = dst.CoverBlocks( kConv1dBlockSize );

  switch( axis )
  {
  case MRI_WIDTH:
    grid.x = coverGrid.x * coverGrid.y;
    grid.y = srcDims.z;
    grid.z = 1;

    MRIConvolveKernelX
    <T,U>
    <<<grid,threads,0,this->stream>>>
    ( src, dst, coverGrid );
    CUDA_CHECK_ERROR_ASYNC( "MRIconvolveKernelX failed!" );
    break;

  case MRI_HEIGHT:
    grid.x = coverGrid.x * coverGrid.y;
    grid.y = srcDims.z;
    grid.z = 1;

    MRIConvolveKernelY
    <T,U>
    <<<grid,threads,0,this->stream>>>
    ( src, dst, coverGrid );
    CUDA_CHECK_ERROR_ASYNC( "MRIconvolveKernelY failed!" );
    break;

  case MRI_DEPTH:
    // Slight change, since we do (x,z) patches
    grid.x = coverGrid.x * coverGrid.z;
    grid.y = srcDims.y;

    MRIConvolveKernelZ
    <T,U>
    <<<grid,threads,0,this->stream>>>
    ( src, dst, coverGrid );
    CUDA_CHECK_ERROR_ASYNC( "MRIconvolveKernelZ failed!" );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Incompatible universe detected." << std::endl;
    std::cerr << "GPU functions are only tested ";
    std::cerr << "in a universe with three spatial dimensions" << std::endl;
    std::cerr << "Please adjust your reality accordingly, ";
    std::cerr << "and try again" << std::endl;
    std::cerr << "MRIframeGPU version was:\n" << src.VersionString() << std::endl;
    exit( EXIT_FAILURE );
  }
}





// ---------------------------------------------------

void MRIconvolve::Convolve3D( const MRI* src, MRI* dst,
                              const int srcFrame, const int dstFrame,
                              const float *kernel,
                              const unsigned int kernelLength )
{
  /*!
  Performs a 3D convolution with a single 1D kernel of
  \a srcFrame of the MRI \a src, storing the result in
  \a dstFrame of the MRI \a dst.
  The convolution kernel is stored in the array \a kernel.
  @param[in] src The source MRI
  @param[out] dst The destination MRI
  @param[in] srcFrame The frame of the source to be convolved
  @param[in] dstFrame The frame of the destination in which the result should be stored
  @param[in] kernel The convolution kernel to use
  @param[in] kernelLength The size of the convolution kernel
  */

  MRIconvolve::tTotal3d.Start();

  // Send the convolution kernel
  this->BindKernel( kernel, kernelLength );

  switch( src->type )
  {
  case MRI_UCHAR:
    this->Conv3DSingleKernelDispatch<unsigned char>( src,
        dst,
        srcFrame,
        dstFrame );
    break;

  case MRI_SHORT:
    this->Conv3DSingleKernelDispatch<short>( src,
        dst,
        srcFrame,
        dstFrame );
    break;

  case MRI_FLOAT:
    this->Conv3DSingleKernelDispatch<float>( src,
        dst,
        srcFrame,
        dstFrame );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised source MRI type " << src->type
              << std::endl;
    exit( EXIT_FAILURE );
  }

  this->UnbindKernel();

  MRIconvolve::tTotal3d.Stop();
}

// --



template<typename T>
void MRIconvolve::Conv3DSingleKernelDispatch( const MRI* src,
    MRI* dst,
    const int srcFrame,
    const int dstFrame ) const
{
  /*!
  Function to run the 3D convolution on the GPU when
  all three directions use the same convolution kernel.
  This routine does not pull intermediate results
  back to the host.
  Assumes that src and dst types are the same.
  Assumes that the texture is already set up on the GPU.
  */

  GPU::Classes::MRIframeGPU<T> frame1, frame2;
  size_t workSize;

  // Do some allocation
  MRIconvolve::tMem.Start();
  frame1.Allocate( src );
  frame2.Allocate( src );
  MRIconvolve::tMem.Stop();

  // Verify (note frame2 verified from dst)
  frame1.VerifyMRI( src );
  frame2.VerifyMRI( dst );

  // Allocate workspace
  workSize = frame1.BufferSize();
  this->AllocateWorkspace( workSize );

  MRIconvolve::tSend.Start();
  frame1.SendFrame( src, srcFrame, this->h_workspace, this->stream );
  MRIconvolve::tSend.Stop();

  MRIconvolve::tCompute3d.Start();
  this->RunGPU1D( frame1, frame2, MRI_WIDTH );
  this->RunGPU1D( frame2, frame1, MRI_HEIGHT );
  this->RunGPU1D( frame1, frame2, MRI_DEPTH );
  MRIconvolve::tCompute3d.Stop();

  MRIconvolve::tRecv.Start();
  frame2.RecvFrame( dst, dstFrame, this->h_workspace, this->stream );
  MRIconvolve::tRecv.Stop();

  CUDA_CHECK_ERROR_ASYNC( "3D convolution  with single kernel failure" );
}


//! Dispatch for 3D convolution of unknown MRI with separate kernels
void MRIconvolve::Convolve3D( const MRI* src,
                              MRI* dst,
                              const int srcFrame,
                              const int dstFrame,
                              const float* xKernel,
                              const unsigned int xL,
                              const float* yKernel,
                              const unsigned int yL,
                              const float* zKernel,
                              const unsigned int zL )
{
  /*!
  Performs a 3D convolution of
  \a srcFrame of the MRI \a src, storing the result in
  \a dstFrame of the MRI \a dst.
  A separate convolution kernel is used for each axis.

  @param[in] src The source MRI
  @param[out] dst The destination MRI
  @param[in] srcFrame The frame of the source to be convolved
  @param[in] dstFrame The frame of the destination in which the result should be stored
  @param[in] xKernel The convolution kernel to use in the x direction
  @param[in] xL The size of the x convolution kernel
  @param[in] yKernel The convolution kernel to use in the y direciton
  @param[in] yL The size of the y convolution kernel
  @param[in] zKernel The z convolution kernel
  @param[in] zL The size of the z convolution kernel
  */
  switch( src->type )
  {
  case MRI_UCHAR:
    this->Convolve3DMultiKernelDispatch<unsigned char>( src,
        dst,
        srcFrame,
        dstFrame,
        xKernel,
        xL,
        yKernel,
        yL,
        zKernel,
        zL );
    break;

  case MRI_SHORT:
    this->Convolve3DMultiKernelDispatch<short>( src, dst,
        srcFrame, dstFrame,
        xKernel, xL,
        yKernel, yL,
        zKernel, zL );
    break;

  case MRI_FLOAT:
    this->Convolve3DMultiKernelDispatch<float>( src, dst,
        srcFrame, dstFrame,
        xKernel, xL,
        yKernel, yL,
        zKernel, zL );
    break;


  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised source MRI type " << src->type
              << std::endl;
    exit( EXIT_FAILURE );
  }
}



template<typename T>
void MRIconvolve::Convolve3DMultiKernelDispatch( const MRI* src,
    MRI* dst,
    const int srcFrame,
    const int dstFrame,
    const float* xKernel,
    const unsigned int xL,
    const float* yKernel,
    const unsigned int yL,
    const float* zKernel,
    const unsigned int zL )
{

  GPU::Classes::MRIframeGPU<T> frame1, frame2;
  size_t workSize;

  MRIconvolve::tTotal3d.Start();

  // Do some allocation
  MRIconvolve::tMem.Start();
  frame1.Allocate( src );
  frame2.Allocate( src );
  MRIconvolve::tMem.Stop();

  // Verify (note frame2 verified from dst)
  frame1.VerifyMRI( src );
  frame2.VerifyMRI( dst );

  // Allocate workspace
  workSize = frame1.BufferSize();
  this->AllocateWorkspace( workSize );

  // Send the data
  MRIconvolve::tSend.Start();
  frame1.Send( src, srcFrame, this->h_workspace, this->stream );
  MRIconvolve::tSend.Stop();

  // Convolve along x
  MRIconvolve::tCompute3d.Start();
  this->BindKernel( xKernel, xL );
  this->RunGPU1D( frame1, frame2, MRI_WIDTH );

  // And y
  this->BindKernel( yKernel, yL );
  this->RunGPU1D( frame2, frame1, MRI_HEIGHT );

  // And z
  this->BindKernel( zKernel, zL );
  this->RunGPU1D( frame1, frame2, MRI_DEPTH );
  MRIconvolve::tCompute3d.Stop();

  // Get the results
  MRIconvolve::tRecv.Start();
  frame2.Recv( dst, dstFrame, this->h_workspace, this->stream );
  MRIconvolve::tRecv.Stop();

  CUDA_CHECK_ERROR_ASYNC( "3D convolution  with multi-kernel failure" );
  MRIconvolve::tTotal3d.Stop();
}

// ---------------------------------------

void MRIconvolve::BindKernel( const float *kernel,
                              const unsigned int nVals )
{
  /*!
  This routine is responsible for preparing the dtl_mriconv1d_kernel
  texture for use.
  The array on the device is padded with an extra zero on each end,
  to save us some explicit boundary condition checks (the texture
  units will handle them).
  There's no need for a matching release, since rebinding
  a texture implicitly unbinds the previous, and the destructor
  handles the memory release
  @param[in] kernel Array containing the kernel values
  @param[in] nVals The number of values in the kernel
  */
  MRIconvolve::tTextureKernel.Start();

  // Allocate and zero GPU memory
  this->AllocateKernel( 2 + nVals );
  CUDA_SAFE_CALL( cudaMemset( this->d_kernel,
                              0,
                              (2+nVals)*sizeof(float) ) );

  // Copy the convolution kernel to the GPU
  // Note the extra offset
  CUDA_SAFE_CALL( cudaMemcpy( &(this->d_kernel[1]),
                              kernel,
                              nVals*sizeof(float),
                              cudaMemcpyHostToDevice ) );

  // Copy the size of the texture to device constant memory
  CUDA_SAFE_CALL( cudaMemcpyToSymbol( dc_mriconv1d_kernel_nVals,
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
                                   this->d_kernel,
                                   cd_kernel,
                                   (2+nVals)*sizeof(float) ) );

  MRIconvolve::tTextureKernel.Stop();
}


void MRIconvolve::UnbindKernel( void )
{
  CUDA_SAFE_CALL( cudaUnbindTexture( dtl_mriconv1d_kernel ) );
}

// ==============================================================


void MRIconvolve::AllocateKernel( const unsigned int nFloats )
{
  /*!
  Allocates space on the device to hold the
  convolution kernel, in necessary.
  @param[in] nFloats Number of floats required
  */
  if( this->kernelAllocSize < nFloats )
  {
    this->ReleaseKernel();

    CUDA_SAFE_CALL( cudaMalloc( (void**)&(this->d_kernel),
                                nFloats * sizeof(float) ) );
    this->kernelAllocSize = nFloats;
  }
}



void MRIconvolve::ReleaseKernel( void )
{
  if( this->d_kernel != NULL )
  {
    cudaFree( this->d_kernel );
    this->kernelAllocSize = 0;
  }
}

// ----------------------------------------------



void MRIconvolve::AllocateWorkspace( const size_t nBytes ) const
{
  if( this->workSize < nBytes )
  {
    this->ReleaseWorkspace();

    MRIconvolve::tHostMem.Start();
    CUDA_SAFE_CALL( cudaHostAlloc( (void**)&(this->h_workspace),
                                   nBytes,
                                   cudaHostAllocDefault ) );
    this->workSize = nBytes;
    MRIconvolve::tHostMem.Stop();
  }
}




void MRIconvolve::ReleaseWorkspace( void ) const
{
  if( h_workspace != NULL )
  {
    MRIconvolve::tHostMem.Start();
    CUDA_SAFE_CALL( cudaFreeHost( h_workspace ) );
    h_workspace = NULL;
    workSize = 0;
    MRIconvolve::tHostMem.Stop();
  }
}

// ----------------------------------------------
//! Dispatch wrapper for 1D convolutions based on dst type
template<typename T>
void MRIconvolve::DispatchWrap1D( const MRI* src, MRI* dst,
                                  const int srcFrame, const int dstFrame,
                                  const int axis ) const
{
  switch( dst->type )
  {
  case MRI_UCHAR:
    this->Dispatch1D<T,unsigned char>( src, dst,
                                       srcFrame, dstFrame,
                                       axis );
    break;

  case MRI_SHORT:
    this->Dispatch1D<T,short>( src, dst,
                               srcFrame, dstFrame,
                               axis );
    break;

  case MRI_FLOAT:
    this->Dispatch1D<T,float>( src, dst,
                               srcFrame, dstFrame,
                               axis );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised destination MRI type " << dst->type
              << std::endl;
    exit( EXIT_FAILURE );
  }
}




// Define all the static members
SciGPU::Utilities::Chronometer MRIconvolve::tMem;
SciGPU::Utilities::Chronometer MRIconvolve::tHostMem;
SciGPU::Utilities::Chronometer MRIconvolve::tSend, MRIconvolve::tRecv;
SciGPU::Utilities::Chronometer MRIconvolve::tTextureKernel;
SciGPU::Utilities::Chronometer MRIconvolve::tCompute1d, MRIconvolve::tTotal1d;
SciGPU::Utilities::Chronometer MRIconvolve::tCompute3d, MRIconvolve::tTotal3d;

}
}



// =================================================

static GPU::Algorithms::MRIconvolve myConvolve;

//! Routine to be called from C code.
MRI* MRIconvolve1d_cuda( const MRI* src, MRI* dst,
                         const float *kernel,
                         const unsigned int kernelLength,
                         const int axis,
                         const int srcFrame, const int dstFrame )
{
  /*!
    Reimplementation of MRIconvolve1d for the GPU.
    This is called from MRIconvolve1d in the main
    Freesurfer suite, after that routine has verified the
    input and output MRI are allocated.
  */

  myConvolve.Convolve1D( src, dst,
                         srcFrame, dstFrame,
                         kernel, kernelLength,
                         axis );

  return( dst );
}

// =================================================

//! Routine to be called from C code
MRI* MRIconvolveGaussian_cuda( const MRI* src, MRI* dst,
                               const float *kernel,
                               const unsigned int kernelLength )
{
  /*!
    Implementation of MRIconvolveGaussian for the GPU.
    Designed to be called form that routine.
    This operates on every frame, as the calling routine,
    but this may not be what you want.
  */

  for( int iFrame=0; iFrame < src->nframes; iFrame++ )
  {
    myConvolve.Convolve3D( src, dst,
                           iFrame, iFrame,
                           kernel, kernelLength );
  }

  return( dst );
}

