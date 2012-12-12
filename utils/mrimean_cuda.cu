/**
 * @file  mrimean_cuda.cu
 * @brief Holds MRI mean function for the GPU
 *
 * Implements MRI mean function on the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.22 $
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

#include "mrimean_cuda.hpp"

#include "mrimean_cuda.h"






namespace GPU
{
namespace Algorithms
{

const unsigned int kMRImeanBlockSize = 16;

//! Min function for ints
__device__ int min( const int& a, const int& b )
{
  if( a < b )
  {
    return( a );
  }
  else
  {
    return( b );
  }
}

//! Max function for ints
__device__ int max( const int& a, const int& b )
{
  if( a > b )
  {
    return( a );
  }
  else
  {
    return( b );
  }
}


//! Kernel to compute x direction means
template<typename T>
__global__ void MRImeanX( const GPU::Classes::MRIframeOnGPU<T> src,
                          GPU::Classes::MRIframeOnGPU<float> dst,
                          const dim3 coverGrid,
                          const unsigned int wSize )
{
  /*!
  Kernel to compute means in the x direction, based on
  the given window size.
  Basically, does a 1D convolution, but with different
  boundary conditions to MRIConvolveKernelX.
  Also, since this is meant to be part of a pipeline,
  the destination type must be float
  */
  const unsigned int by = blockIdx.x / coverGrid.x;
  const unsigned int bx = blockIdx.x % coverGrid.x;
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
  const unsigned int myxMax = min( dst.dims.x-1, ix + wHalf );

  // Again, declare int to remove need for some casts
  const int patchSize = NextMultiple( max(wHalf,1), kMRImeanBlockSize );

  // Accumulator
  float myVal = 0;

  __shared__ float currPatch[kMRImeanBlockSize][kMRImeanBlockSize];

  // Calculate patch limits (note integer declarations avoid -ve trouble)
  const int xDimRound = NextMultiple( src.dims.x, kMRImeanBlockSize );
  const int xFirst = max( 0, ixStart - patchSize );
  const int xLast  = min( xDimRound - kMRImeanBlockSize,
                          ixStart + patchSize );

  for( int xBegin = xFirst;
       xBegin <= xLast;
       xBegin += kMRImeanBlockSize )
  {
    // Load the patch
    currPatch[ty][tx] = src( xBegin+tx, iy, iz );
    __syncthreads();

    // Accumulate desired values
    for( unsigned int i=0; i<kMRImeanBlockSize; i++ )
    {
      int actx = xBegin + i;

      if( (actx>=myxMin) && (actx<=myxMax) )
      {
        myVal += currPatch[ty][i];
      }

    }

    __syncthreads();
  }

  // Save result
  if( dst.InVolume( ix, iy, iz ) )
  {
    dst(ix,iy,iz) = dst.ConvertFloat( myVal );
  }
}


//! Kernel to compute y direction means
template<typename T>
__global__ void MRImeanY( const GPU::Classes::MRIframeOnGPU<T> src,
                          GPU::Classes::MRIframeOnGPU<float> dst,
                          const dim3 coverGrid,
                          const unsigned int wSize )
{
  /*!
  Kernel to compute means in the y direction, based on
  the given window size.
  Basically, does a 1D convolution, but with different
  boundary conditions to MRIConvolveKernelY.
  Also, since this is meant to be part of a pipeline,
  the destination type must be float
  */
  const unsigned int by = blockIdx.x / coverGrid.x;
  const unsigned int bx = blockIdx.x % coverGrid.x;
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
  const unsigned int myyMax = min( dst.dims.y-1, iy + wHalf );

  // Again, declare int to remove need for some casts
  const int patchSize = NextMultiple( max(wHalf,1), kMRImeanBlockSize );

  // Accumulator
  float myVal = 0;

  __shared__ float currPatch[kMRImeanBlockSize][kMRImeanBlockSize];

  // Calculate patch limits (note integer declarations avoid -ve trouble)
  const int yDimRound = NextMultiple( src.dims.y, kMRImeanBlockSize );

  const int yFirst = max( 0, iyStart - patchSize );
  const int yLast  = min( yDimRound - kMRImeanBlockSize,
                          iyStart + patchSize );

  for( int yBegin = yFirst;
       yBegin <= yLast;
       yBegin += kMRImeanBlockSize )
  {
    // Load the patch
    currPatch[ty][tx] = src( ix, yBegin+ty, iz );
    __syncthreads();

    // Accumulate desired values
    for( unsigned int i=0; i<kMRImeanBlockSize; i++ )
    {
      int acty = yBegin + i;

      if( (acty>=myyMin) && (acty<=myyMax) )
      {
        myVal += currPatch[i][tx];
      }

    }

    __syncthreads();
  }

  // Save result
  if( dst.InVolume( ix, iy, iz ) )
  {
    dst(ix,iy,iz) = dst.ConvertFloat( myVal );
  }
}


//! Kernel to compute z direction means
template<typename T>
__global__ void MRImeanZ( const GPU::Classes::MRIframeOnGPU<T> src,
                          GPU::Classes::MRIframeOnGPU<float> dst,
                          const dim3 coverGrid,
                          const unsigned int wSize )
{
  /*!
  Kernel to compute means in the x direction, based on
  the given window size.
  Basically, does a 1D convolution, but with different
  boundary conditions to MRIConvolveKernelZ.
  Also, since this is meant to be part of a pipeline,
  the destination type must be float
  */
  const unsigned int bz = blockIdx.x / coverGrid.x;
  const unsigned int bx = blockIdx.x % coverGrid.x;
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
  const unsigned int myzMax = min( dst.dims.z-1, iz + wHalf );

  // Again, declare int to remove need for some casts
  const int patchSize = NextMultiple( max(wHalf,1), kMRImeanBlockSize );

  // Accumulator
  float myVal = 0;

  __shared__ float currPatch[kMRImeanBlockSize][kMRImeanBlockSize];

  // Calculate patch limits (note integer declarations avoid -ve trouble)
  const int zDimRound = NextMultiple( src.dims.z, kMRImeanBlockSize );

  const int zFirst = max( 0, izStart - patchSize );
  const int zLast  = min( zDimRound - kMRImeanBlockSize,
                          izStart + patchSize );

  for( int zBegin = zFirst;
       zBegin <= zLast;
       zBegin += kMRImeanBlockSize )
  {
    // Load the patch
    currPatch[tz][tx] = src( ix, iy, zBegin+tz );
    __syncthreads();

    // Accumulate desired values
    for( unsigned int i=0; i<kMRImeanBlockSize; i++ )
    {
      int actz = zBegin + i;

      if( (actz>=myzMin) && (actz<=myzMax) )
      {
        myVal += currPatch[i][tx];
      }

    }

    __syncthreads();
  }

  // Save result
  if( dst.InVolume( ix, iy, iz ) )
  {
    dst(ix,iy,iz) = dst.ConvertFloat( myVal );
  }
}


//! Kernel to normalise means
template<typename U>
__global__ void MRImeanNormal( const GPU::Classes::MRIframeOnGPU<float> src,
                               GPU::Classes::MRIframeOnGPU<U> dst,
                               const dim3 coverGrid,
                               const unsigned int wSize )
{
  /*!
  Kernel to normalise the means computed by the earlier
  stages.
  As such, the input type must be a float
  */
  const unsigned int by = blockIdx.x / coverGrid.x;
  const unsigned int bx = blockIdx.x % coverGrid.x;
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
  const unsigned int myxMax = min( dst.dims.x-1, ix + wHalf );
  const unsigned int myyMin = max( 0           , iy - wHalf );
  const unsigned int myyMax = min( dst.dims.y-1, iy + wHalf );
  const unsigned int myzMin = max( 0           , iz - wHalf );
  const unsigned int myzMax = min( dst.dims.z-1, iz + wHalf );


  const unsigned long myVolume = ( myxMax - myxMin + 1 ) *
                                 (myyMax - myyMin + 1 ) *
                                 (myzMax - myzMin + 1 );

  if( dst.InVolume( ix, iy, iz ) )
  {
    dst( ix, iy, iz ) = dst.ConvertFloat( src( ix, iy, iz ) / myVolume );
  }
}






// ######################################################

// Define the static timers
SciGPU::Utilities::Chronometer MRImean::tMem, MRImean::tHostMem;
SciGPU::Utilities::Chronometer MRImean::tSend, MRImean::tRecv;
SciGPU::Utilities::Chronometer MRImean::tCompute;
SciGPU::Utilities::Chronometer MRImean::tTotal;

// =======================

template<typename T>
void MRImean::DispatchWrap( const MRI* src, MRI* dst,
                            const unsigned int wSize,
                            const int srcFrame,
                            const int dstFrame ) const
{
  switch( dst->type )
  {
  case MRI_UCHAR:
    this->MeanDispatch<T,unsigned char>( src, dst, wSize,
                                         srcFrame, dstFrame );
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


void MRImean::Allocate( const size_t nBytes ) const
{
  if( this->workSize < nBytes )
  {
    this->Release();

    MRImean::tHostMem.Start();
    CUDA_SAFE_CALL( cudaHostAlloc( (void**)&(this->h_workspace),
                                   nBytes,
                                   cudaHostAllocDefault ) );
    this->workSize = nBytes;
    MRImean::tHostMem.Stop();
  }
}


//! Releases internal pinned memory buffer
void MRImean::Release( void ) const
{
  if( h_workspace != NULL )
  {
    MRImean::tHostMem.Start();
    cudaFreeHost( h_workspace );
    h_workspace = NULL;
    workSize = 0;
    MRImean::tHostMem.Stop();
  }
}


// ========================================================



MRImean::~MRImean( void )
{
  this->Release();
}


void MRImean::ShowTimings( void )
{
#ifdef CUDA_SHOW_TIMINGS
  std::cout << "==================================" << std::endl;
  std::cout << "GPU Mean timers" << std::endl;
  std::cout << "---------------" << std::endl;
#ifndef CUDA_FORCE_SYNC
  std::cout << "WARNING: CUDA_FORCE_SYNC not #defined" << std::endl;
  std::cout << "Timings may not be accurate" << std::endl;
#endif
  std::cout << std::endl;

  std::cout << "Host Memory : " << MRImean::tHostMem << std::endl;
  std::cout << "GPU Memory  : " << MRImean::tMem << std::endl;
  std::cout << "Send        : " << MRImean::tSend << std::endl;
  std::cout << "Compute     : " << MRImean::tCompute << std::endl;
  std::cout << "Receive     : " << MRImean::tRecv << std::endl;
  std::cout << "--------------" << std::endl;
  std::cout << "Total : " << MRImean::tTotal << std::endl;
  std::cout << "==================================" << std::endl;
#endif
}




void MRImean::DoMean( const MRI* src, MRI* dst,
                      const unsigned int wSize,
                      const unsigned int srcFrame,
                      const unsigned int dstFrame ) const
{
  switch( src->type )
  {
  case MRI_UCHAR:
    this->DispatchWrap<unsigned char>( src, dst, wSize,
                                       srcFrame, dstFrame );
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
void MRImean::MeanDispatch( const MRI* src, MRI* dst,
                            const unsigned int wSize,
                            const int srcFrame,
                            const int dstFrame ) const
{
  /*!
  Templated dispatch routine for MRI mean function on the
  GPU.
  Transfers data to the GPU, and retrieves the results
  */
  MRImean::tTotal.Start();

  GPU::Classes::MRIframeGPU<T> srcGPU;
  GPU::Classes::MRIframeGPU<U> dstGPU;

  size_t srcWorkSize, dstWorkSize;

  // Allocate the GPU arrays
  MRImean::tMem.Start();
  srcGPU.Allocate( src );
  dstGPU.Allocate( dst );
  MRImean::tMem.Stop();

  // Put in some sanity checks
  srcGPU.VerifyMRI( src );
  dstGPU.VerifyMRI( dst );

  // Allocate the workspace array
  srcWorkSize = srcGPU.BufferSize();
  dstWorkSize = dstGPU.BufferSize();

  if( srcWorkSize > dstWorkSize )
  {
    this->Allocate( srcWorkSize );
  }
  else
  {
    this->Allocate( dstWorkSize );
  }

  // Send the source data
  MRImean::tSend.Start();
  srcGPU.SendFrame( src, srcFrame, this->h_workspace, this->stream );
  MRImean::tSend.Stop();

  // Run the filter
  this->RunGPU( srcGPU, dstGPU, wSize );

  // Get the results
  MRImean::tRecv.Start();
  dstGPU.RecvFrame( dst, dstFrame, this->h_workspace, this->stream );
  MRImean::tRecv.Stop();

  CUDA_CHECK_ERROR( "Mean filtering failure" );

  MRImean::tTotal.Stop();
}



//! Runs the mean filtering kernel
template<typename T, typename U>
void MRImean::RunGPU( const GPU::Classes::MRIframeGPU<T> &src,
                      GPU::Classes::MRIframeGPU<U> &dst,
                      const unsigned int wSize ) const
{
  /*!
  Runs the mean filtering kernel on the GPU.
  Assumes most things are properly set already
  */
  const dim3 srcDims = src.GetDims();
  const dim3 dstDims = dst.GetDims();

  // Check dimensions
  if( srcDims != dstDims )
  {
    std::cerr << __FUNCTION__ << ": Dimension mismatch"
              << std::endl;
    exit( EXIT_FAILURE );
  }


  // We need intermediates which are floats
  GPU::Classes::MRIframeGPU<float> f1, f2;

  // Get correctly sized intermediates
  MRImean::tMem.Start();
  f1.Allocate( src );
  f2.Allocate( src );
  MRImean::tMem.Stop();


  // Do the three convolutions. Recall objects have same dims
  dim3 grid, threads;


  threads.x = threads.y = kMRImeanBlockSize;
  threads.z = 1;

  const dim3 coverGrid = dst.CoverBlocks( kMRImeanBlockSize );

  grid.x = coverGrid.x * coverGrid.y;
  grid.y = dstDims.z;
  grid.z = 1;

  MRImean::tCompute.Start();

  // Do the X direction
  MRImeanX<T>
  <<<grid,threads,0,this->stream>>>
  ( src, f1, coverGrid, wSize );
  CUDA_CHECK_ERROR_ASYNC( "MRImeanX kernel failed" );

  // Do the Y direction
  MRImeanY<float>
  <<<grid,threads,0,this->stream>>>
  ( f1, f2, coverGrid, wSize );
  CUDA_CHECK_ERROR_ASYNC( "MRImeanY kernel failed" );

  // Slight change for Z direction
  grid.x = coverGrid.x * coverGrid.z;
  grid.y = dstDims.y;
  MRImeanZ<float>
  <<<grid,threads,0,this->stream>>>
  ( f2, f1, coverGrid, wSize );
  CUDA_CHECK_ERROR_ASYNC( "MRImeanZ kernel failed" );

  // Normalise
  grid.x = coverGrid.x * coverGrid.y;
  grid.y = dstDims.z;
  MRImeanNormal<U>
  <<<grid,threads,0,this->stream>>>
  ( f1, dst, coverGrid, wSize );
  CUDA_CHECK_ERROR_ASYNC( "MRImeanNormal failed!" );

  MRImean::tCompute.Stop();
}

}
}


static GPU::Algorithms::MRImean myMean;


MRI* MRImean_cuda( const MRI* src, MRI* dst,
                   const unsigned int wSize )
{
  /*!
    Wrapper around the GPU routine, to be called from the
    original MRImean routine.
    Note that the frames default to zero, per the original
    MRImean routine.
  */

  myMean.DoMean( src, dst, wSize );

  return( dst );
}
