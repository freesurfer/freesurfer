/**
 * @file  volumegpu_impl.hpp
 * @brief Holds implementations for volume data on the GPU
 *
 * In order to keep the volumegpu.hpp file reasonably 'clean'
 * the implementations of the longer routines are kept
 * in this file
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.10 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <cuda_runtime.h>


namespace GPU
{
namespace Classes
{

// ================================================================
// Memory Management


template<typename T>
void VolumeGPU<T>::Allocate( const dim3 myDims )
{
  /*!
  Allocates GPU memory to hold a volume of the
  given size.
  This memory may include padding for GPU
  memory alignment requirements
  @param[in] myDims The dimensions required
  */

  // Check if we can re-use current memory
  if( myDims == this->dims )
  {
    return;
  }

  // Release existing memory
  this->Release();

  // Set the dimensions
  this->dims = myDims;


  // Make the extent
  cudaExtent tmpExtent = ExtentFromDims( this->dims );
  tmpExtent.width *= sizeof(T);

  // Allocate the memory
  CUDA_SAFE_CALL( cudaMalloc3D( &(this->d_data), tmpExtent ) );
}



// ---------------------------------

template<typename T>
void VolumeGPU<T>::Release( void )
{
  if( this->d_data.ptr != NULL )
  {
    cudaFree( d_data.ptr );
  }
  this->dims = make_uint3(0,0,0);
  this->d_data = make_cudaPitchedPtr(NULL,0,0,0);
}


// ---------------------------------

template<typename T>
size_t VolumeGPU<T>::BufferSize( void ) const
{
  /*!
  Prior to copying data to the GPU, the host must
  organise the volume into contiguous pinned memory.
  This routine supplies the number of bytes required
  by the current class.
  */
  size_t nElements;

  nElements = this->dims.x * this->dims.y * this->dims.z;

  return( nElements * sizeof(T) );
}


// ---------------------------------

template<typename T>
T* VolumeGPU<T>::AllocateHostBuffer( void ) const
{
  /*!
  Allocates a pinned memory buffer on the
  host of sufficient size to hold the volume
  in contiguous memory.
  It will abort if the allocation fails,
  so it should not be necessary to check
  the return value.
  Note that there is no provision within this
  class for releasing the memory so allocated.
  It is simply a convenience wrapper around
  cudaHostAlloc.
  */

  T* h_buf = NULL;

  CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_buf,
                                 this->BufferSize(),
                                 cudaHostAllocDefault ) );

  return( h_buf );
}


// ---------------------------------

#ifdef __CUDACC__
template<typename T>
__global__
void ZeroKernel( VolumeArgGPU<T> arr )
{
  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  if( arr.InVolume(ix,iy,0) )
  {
    for( unsigned int iz=0; iz<arr.dims.z; iz++ )
    {
      arr(ix,iy,iz) = 0;
    }
  }

}

template<typename T>
void VolumeGPU<T>::Zero( void )
{

  const unsigned int kKernelSize = 16;

  dim3 grid, threads;
  threads.x = threads.y = kKernelSize;
  threads.z = 1;

  grid = this->CoverBlocks( kKernelSize );
  grid.z = 1;

  VolumeArgGPU<T> vag( this->dims,
                       this->d_data.ptr,
                       this->d_data.pitch );

  ZeroKernel<<<grid,threads>>>( vag );
  CUDA_CHECK_ERROR( "ZeroKernel failed!\n" );
}
#endif

// ---------------------------------

template<typename T>
void VolumeGPU<T>::Copy( const VolumeGPU<T>& src,
                         const cudaStream_t myStream )
{
  this->Allocate( src.dims );

  cudaMemcpy3DParms copyParams = {0};

  copyParams.srcPtr = src.d_data;
  copyParams.dstPtr = this->d_data;
  copyParams.extent = ExtentFromDims( this->dims );
  copyParams.extent.width *= sizeof(T);
  copyParams.kind = cudaMemcpyDeviceToDevice;

  CUDA_SAFE_CALL( cudaMemcpy3DAsync( &copyParams, myStream ) );
}


// ================================================================
// Data Transfer

template<typename T>
void VolumeGPU<T>::SendBuffer( const T* const h_buffer,
                               const cudaStream_t stream )
{
  /*!
  Sends the given contiguous buffer on the
  host to the GPU.
  The data in the buffer must be ordered by
  the Index1D method.
  The copy is done asynchronously in the optionally
  defined stream.
  If no stream is given, it defaults to stream 0.
  Accordingly, this routine may return before the
  copy is complete - the user is responsible for
  stream management.
  @param[in] h_buffer The buffer of length BufferSize() to be send
  @param[in] stream The CUDA stream in which the copy should occur
  */

  cudaMemcpy3DParms copyParams = {0};

  copyParams.srcPtr = make_cudaPitchedPtr( (void*)h_buffer,
                      this->dims.x*sizeof(T),
                      this->dims.x,
                      this->dims.y );
  copyParams.dstPtr = this->d_data;
  copyParams.extent = ExtentFromDims( this->dims );
  copyParams.extent.width *= sizeof(T);
  copyParams.kind = cudaMemcpyHostToDevice;

  CUDA_SAFE_CALL_ASYNC( cudaMemcpy3DAsync( &copyParams, stream ) );
}


// ---------------------------------

template<typename T>
void VolumeGPU<T>::RecvBuffer( T* const h_buffer,
                               const cudaStream_t stream ) const
{
  /*!
  Retrieves data from the GPU into the
  contiguous buffer on the host.
  The copy is done asynchronously, so the
  user is responsible for stream management.
  If no stream is given, it defaults to stream 0.
  @param[in] h_buffer The buffer of length BufferSize() to be received
  @param[in] stream The CUDA stream in which the copy should occur
  */

  cudaMemcpy3DParms cpyPrms = {0};
  cpyPrms.srcPtr = this->d_data;
  cpyPrms.dstPtr = make_cudaPitchedPtr( (void*)h_buffer,
                                        this->dims.x*sizeof(T),
                                        this->dims.x,
                                        this->dims.y );
  cpyPrms.extent = ExtentFromDims( this->dims );
  cpyPrms.extent.width *= sizeof(T);
  cpyPrms.kind = cudaMemcpyDeviceToHost;

  CUDA_SAFE_CALL_ASYNC( cudaMemcpy3DAsync( &cpyPrms, stream ) );
}


// ================================================================
// Other routines

template<typename T>
dim3 VolumeGPU<T>::CoverBlocks( const unsigned int threadCount ) const
{
  /*!
  Method to return the number of blocks which
  will be required to cover the volume of cubic blocks
  of size threadCount.
  This is provided as an easy way of configuring
  kernel calls.
  */
  dim3 grid;

  grid.x = DivRoundUp( this->dims.x, threadCount );
  grid.y = DivRoundUp( this->dims.y, threadCount );
  grid.z = DivRoundUp( this->dims.z, threadCount );

  return( grid );
}


// ---------------------------------

template<typename T>
unsigned int VolumeGPU<T>::Index1D( const unsigned int ix,
                                    const unsigned int iy,
                                    const unsigned int iz ) const
{
  /*!
  Computes the offset into an unpadded contiguous
  array on the host of the given element location.
  This is for use by send and receive methods in
  derived classes
  */
  return( ix + ( this->dims.x * ( iy + ( this->dims.y * iz ) ) ) );
}


// ---------------------------------

template<typename T>
dim3 VolumeGPU<T>::Index3D( unsigned int idx ) const
{
  /*!
  Computes the x, y and z location of a given
  1D index.
  It is the inverse of the Index1D method
  Note that idx is not declared const.
  We rely on it being passed by value.
  */

  dim3 loc;

  loc.x = idx % this->dims.x;
  idx /= this->dims.x;

  loc.y = idx % this->dims.y;
  idx /= this->dims.y;

  loc.z = idx;

  return( loc );
}

}
}
