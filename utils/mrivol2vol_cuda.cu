/**
 * @file  mrivol2vol_cuda.cu
 * @brief Holds MRI volume to volume transform for CUDA
 *
 * Implements MRIvol2vol for the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.24 $
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

#include <memory>

#include "chronometer.hpp"
#include "cudacheck.h"


#include "mriframegpu.hpp"
#include "affinegpu.hpp"


#include "mrivol2vol_cuda.h"
#include "mrivol2vol_cuda.hpp"

// ==============================================================

//! Texture reference for an unsigned char source
texture<unsigned char, 3, cudaReadModeNormalizedFloat> dt_src_uchar;

//! Texture reference for a short source
texture<short, 3, cudaReadModeNormalizedFloat> dt_src_short;

//! Texture reference for a float source
texture<float, 3, cudaReadModeElementType> dt_src_float;


namespace GPU
{
namespace Algorithms
{

const unsigned int kVol2VolBlockSize = 16;

// ---------------------------------------


template<>
GPU::Classes::CTfactory* MRIvol2vol::BindSourceTexture<unsigned char>( const GPU::Classes::MRIframeGPU<unsigned char>& src,
    const int InterpMode )
{
  /*!
  Binds the unsigned char texture
  */
  GPU::Classes::CTfactory *tmp;

  // Set the interpolation
  switch( InterpMode )
  {
  case SAMPLE_NEAREST:
    tmp = new GPU::Classes::CTfactory( src, dt_src_uchar,
                                       cudaFilterModePoint );
    break;

  case SAMPLE_TRILINEAR:
    tmp = new GPU::Classes::CTfactory( src, dt_src_uchar,
                                       cudaFilterModeLinear );
    break;

  default:
    std::cerr << __FUNCTION__ << ": Unrecognised InterpMode "
              << InterpMode << std::endl;
    exit( EXIT_FAILURE );
  }


  return( tmp );
}


template<>
GPU::Classes::CTfactory* MRIvol2vol::BindSourceTexture<float>( const GPU::Classes::MRIframeGPU<float>& src,
    const int InterpMode )
{
  /*!
  Binds the float texture
  */

  GPU::Classes::CTfactory *tmp;

  // Set the interpolation
  switch( InterpMode )
  {
  case SAMPLE_NEAREST:
    tmp = new GPU::Classes::CTfactory( src, dt_src_float,
                                       cudaFilterModePoint );
    break;

  case SAMPLE_TRILINEAR:
    tmp = new GPU::Classes::CTfactory( src, dt_src_float,
                                       cudaFilterModeLinear );
    break;

  default:
    std::cerr << __FUNCTION__ << ": Unrecognised InterpMode "
              << InterpMode << std::endl;
    exit( EXIT_FAILURE );
  }


  return( tmp );
}

template<>
GPU::Classes::CTfactory* MRIvol2vol::BindSourceTexture<short>( const GPU::Classes::MRIframeGPU<short>& src,
    const int InterpMode )
{
  /*!
  Binds the short texture
  */
  GPU::Classes::CTfactory *tmp;

  // Set the interpolation
  switch( InterpMode )
  {
  case SAMPLE_NEAREST:
    tmp = new GPU::Classes::CTfactory( src, dt_src_short,
                                       cudaFilterModePoint );
    break;

  case SAMPLE_TRILINEAR:
    tmp = new GPU::Classes::CTfactory( src, dt_src_short,
                                       cudaFilterModeLinear );
    break;

  default:
    std::cerr << __FUNCTION__ << ": Unrecognised InterpMode "
              << InterpMode << std::endl;
    exit( EXIT_FAILURE );
  }


  return( tmp );

}


// ---------------------------------------



template<> void MRIvol2vol::UnbindSourceTexture<unsigned char>( void )
{
  CUDA_SAFE_CALL( cudaUnbindTexture( dt_src_uchar ) );
}

template<> void MRIvol2vol::UnbindSourceTexture<float>( void )
{
  CUDA_SAFE_CALL( cudaUnbindTexture( dt_src_float ) );
}

template<> void MRIvol2vol::UnbindSourceTexture<short>( void )
{
  CUDA_SAFE_CALL( cudaUnbindTexture( dt_src_short ) );
}

// ---------------------------------------

//! Templated texture fetch
template<typename T>
__device__ float FetchSrcVoxel( const float3 r )
{
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
__device__ float FetchSrcVoxel<unsigned char>( const float3 r )
{

  float texVal;

  // Do the fetch, remembering texel centre offset
  texVal = tex3D( dt_src_uchar, r.x+0.5f, r.y+0.5f, r.z+0.5f );

  // We configured to return a normalised float, so get back in range
  texVal *= UCHAR_MAX;

  return( texVal );
}

template<>
__device__ float FetchSrcVoxel<short>( const float3 r )
{

  float texVal;

  // Do the fetch, remembering texel centre offset
  texVal = tex3D( dt_src_short, r.x+0.5f, r.y+0.5f, r.z+0.5f );

  // We configured to return a normalised float, so get back in range
  texVal *= SHRT_MAX;

  return( texVal );
}


template<>
__device__ float FetchSrcVoxel<float>( const float3 r )
{
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
                                  GPU::Classes::AffineTransformation afTrans,
                                  const dim3 coverGrid )
{
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
  const unsigned int by = blockIdx.x / coverGrid.x;
  const unsigned int bx = blockIdx.x % coverGrid.x;
  const unsigned int iz = blockIdx.y;
  const unsigned int tx = threadIdx.x;
  const unsigned int ty = threadIdx.y;

  const unsigned int ixStart = bx * kVol2VolBlockSize;
  const unsigned int iyStart = by * kVol2VolBlockSize;
  const unsigned int myx = ixStart + tx;
  const unsigned int myy = iyStart + ty;

  float3 rOut = afTrans.transform( make_float3( myx, myy, iz ) );

  float res;

  const float tol = 0.5f;

  if( dst.InFuzzyVolume( rOut, tol ) )
  {
    res = FetchSrcVoxel<T>( rOut );
  }
  else
  {
    res = 0;
  }

  if( dst.InVolume( myx, myy, iz ) )
  {
    dst( myx, myy, iz ) = dst.ConvertFloat( res );
  }
}


template<typename T, typename U>
void MRIvol2vol::transformGPU( const GPU::Classes::MRIframeGPU<T> &src,
                               GPU::Classes::MRIframeGPU<U> &dst,
                               const GPU::Classes::AffineTransformation& transform,
                               const int InterpMode,
                               const cudaStream_t myStream )
{
  /*!
  This is the dispatch function for the vol2vol kernel.
  The src MRI must already be in its CUDA array, and the
  dst MRI must be padded
  */

  // Run through a couple of sanity checks



  // Bind the texture and array
  std::auto_ptr<GPU::Classes::CTfactory> srcCT( this->BindSourceTexture( src, InterpMode ) );


  const dim3 dstDims = dst.GetDims();
  dim3 grid, threads;

  threads.x = threads.y = kVol2VolBlockSize;
  threads.z = 1;

  const dim3 coverGrid = dst.CoverBlocks( kVol2VolBlockSize );

  grid.x = coverGrid.x * coverGrid.y;
  grid.y = dstDims.z;
  grid.z = 1;

  // Main dispatch
  switch( InterpMode )
  {

  case SAMPLE_NEAREST:
  case SAMPLE_TRILINEAR:
    MRIVol2VolKernel<T,U>
    <<<grid,threads>>>( dst, transform, coverGrid );
    CUDA_CHECK_ERROR_ASYNC( "MRIVol2VolKernel call failed!\n" );
    break;

  case SAMPLE_SINC:
    std::cerr << __FUNCTION__ << ": Sinc interpolation not yet supported"
              << std::endl;
    exit( EXIT_FAILURE );
  }


  // Unbind the texture
  this->UnbindSourceTexture<T>();
}


//! Routine to use GPU for vol2vol on all MRI frames
template<typename T, typename U>
void MRIvol2vol::AllFramesGPU( const MRI* src, MRI* targ,
                               const MATRIX* transformMatrix,
                               const int InterpMode,
                               const float param )
{

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
  dstGPU.Allocate( targ );
  tVol2VolMem.Stop();

  // Sanity check
  srcGPU.VerifyMRI( src );
  dstGPU.VerifyMRI( targ );

  // Allocate workspace array
  tVol2VolMemHost.Start();
  srcWorkSize = srcGPU.BufferSize();
  dstWorkSize = dstGPU.BufferSize();

  if( srcWorkSize > dstWorkSize )
  {
    CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_workspace,
                                   srcWorkSize,
                                   cudaHostAllocDefault ) );
  }
  else
  {
    CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_workspace,
                                   dstWorkSize,
                                   cudaHostAllocDefault ) );
  }
  tVol2VolMemHost.Stop();

  // Loop over frames
  for( unsigned int iFrame=0; iFrame < src->nframes; iFrame++ )
  {
    // Get the next source frame
    tVol2VolMRISendFrame.Start();
    srcGPU.SendFrame( src, iFrame, h_workspace );
    tVol2VolMRISendFrame.Stop();

    // Run the convolution
    tVol2VolCompute.Start();
    this->transformGPU( srcGPU, dstGPU, myTransform, InterpMode );
    tVol2VolCompute.Stop();

    // Get the results back
    tVol2VolMRIRecvFrame.Start();
    dstGPU.RecvFrame( targ, iFrame, h_workspace );
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
void MRIvol2vol::AllFramesDstDispatch( const MRI* src, MRI* targ,
                                       const MATRIX* transformMatrix,
                                       const int InterpMode,
                                       const float param )
{
  /*!
  A thin wrapper to produce the correct template
  based on the type of the target MRI
  */

  switch( targ->type )
  {
  case MRI_UCHAR:
    this->AllFramesGPU<T,unsigned char>( src, targ,
                                         transformMatrix,
                                         InterpMode,
                                         param );
    break;

  case MRI_FLOAT:
    this->AllFramesGPU<T,float>( src, targ,
                                 transformMatrix,
                                 InterpMode,
                                 param );
    break;

  case MRI_SHORT:
    this->AllFramesGPU<T,short>( src, targ,
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

//! Dispatch routine on src MRI
void MRIvol2vol::transform( const MRI* src, MRI* targ,
                            const MATRIX* transformMatrix,
                            const int InterpMode,
                            const float param )
{

  switch( src->type )
  {
  case MRI_UCHAR:
    this->AllFramesDstDispatch<unsigned char>( src, targ,
        transformMatrix,
        InterpMode,
        param );
    break;

  case MRI_FLOAT:
    this->AllFramesDstDispatch<float>( src, targ,
                                       transformMatrix,
                                       InterpMode,
                                       param );
    break;

  case MRI_SHORT:
    this->AllFramesDstDispatch<short>( src, targ,
                                       transformMatrix,
                                       InterpMode,
                                       param );
    break;


  default:
    std::cerr << __FUNCTION__ << ": Unrecognised src data type "
              << src->type << std::endl;
    exit( EXIT_FAILURE );
  }
}


// ===================================================================

// Define the static members
SciGPU::Utilities::Chronometer MRIvol2vol::tVol2VolMem;
SciGPU::Utilities::Chronometer MRIvol2vol::tVol2VolMemHost;
SciGPU::Utilities::Chronometer MRIvol2vol::tVol2VolMRISendFrame;
SciGPU::Utilities::Chronometer MRIvol2vol::tVol2VolMRIRecvFrame;
SciGPU::Utilities::Chronometer MRIvol2vol::tVol2VolCompute;
SciGPU::Utilities::Chronometer MRIvol2vol::tVol2VolTotal;


void MRIvol2vol::ShowTimings( void )
{
  std::cout << "=============================================" << std::endl;
  std::cout << "GPU Vol2Vol timers" << std::endl;
  std::cout << "------------------" << std::endl;
#ifndef CUDA_FORCE_SYNC
  std::cout << "WARNING: CUDA_FORCE_SYNC not #defined" << std::endl;
  std::cout << "Timings may not be accurate" << std::endl;
#endif
  std::cout << std::endl;

  std::cout << "MRIVol2VolAllFramesGPU" << std::endl;
  std::cout << "Host Memory : " << MRIvol2vol::tVol2VolMemHost << std::endl;
  std::cout << "GPU Memory  : " << MRIvol2vol::tVol2VolMem << std::endl;
  std::cout << "Send Frame  : " << MRIvol2vol::tVol2VolMRISendFrame << std::endl;
  std::cout << "Compute     : " << MRIvol2vol::tVol2VolCompute << std::endl;
  std::cout << "Recv Frame  : " << MRIvol2vol::tVol2VolMRIRecvFrame << std::endl;
  std::cout << "---------------------" << std::endl;
  std::cout << "Total : " << MRIvol2vol::tVol2VolTotal << std::endl;

  std::cout << "=============================================" << std::endl;
}
}
}


// ======================================================

static GPU::Algorithms::MRIvol2vol myvol2vol;


int MRIvol2vol_cuda( const MRI* src, MRI* targ,
                     const MATRIX* transformMatrix,
                     const int InterpMode,
                     const float param )
{
  /*!
    Reimplementation of MRIvol2vol for the GPU.
    Is meant to be called from within MRIvol2vol,
    once that routine has performed necessary checks
  */

  myvol2vol.transform( src, targ, transformMatrix, InterpMode, param );


  return( 0 );
}



