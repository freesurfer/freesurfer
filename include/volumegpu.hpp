/**
 * @file  volumegpu.hpp
 * @brief Holds templated datatype for volume data on the GPU
 *
 * Holds a templated datatype for volume data on the GPU.
 * This should be specialised by other datatypes as required
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.29 $
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

#include <cstdlib>

#include <iostream>

#include "cudacheck.h"
#include "cudatypeutils.hpp"

#ifndef VOLUME_GPU_CUDA_H
#define VOLUME_GPU_CUDA_H

//! Namespace to hold everything related to the GPU
namespace GPU
{

//! Namespace to hold the datatypes used by the GPU
namespace Classes
{

//! Templated ancilliary class for use in kernel calls
/*!
  This is an ancilliary class which should be used in
  kernel arguments when manipulating a VolumeGPU class.
  It exists because CUDA kernels are invoked with
  copies of the actual arguments given.
  The VolumeGPU class has the copy constructor declared
  private, since casually copying them around is
  going to kill performance, and if you just copy pointers,
  then when the kernel exits and the destructor is called
  on the copy, all of your memory just vanished.
  The solution is to declare kernel arguments as this type:
  \code
  template<typename T, typename U>
  __global__ void MyKernel( const VolumeArgGPU<T> a,
                            VolumeArgGPU<U> b ) {
     ...
  }
  \endcode
  Appropriate conversion operators are provided so that
  the kernel can be invoked using VolumeArg types.
  However, if templating then the invocation must be
  explicit, in order to enable the compiler to work out
  what you really want to happen
  \code
  VolumeGPU<float> a;
  VolumeGPU<short> b;
  ...
  MyKernel<float,short><<<grid,threads>>>( a, b );
  \endcode
*/
template<typename T>
class VolumeArgGPU
{
public:

  //! Size of the volume
  const dim3 dims;


  // --------------------------------------
  // Constructors

  //! Default constructor
  VolumeArgGPU( void ) : dims(make_uint3(0,0,0)),
    pitchedPtr(NULL),
    dataPitch(0) {};

  //! Constructor from inputs
  VolumeArgGPU( const dim3 myDims,
                void* const myPitchedPtr,
                const size_t myPitch ) : dims(myDims),
    pitchedPtr(myPitchedPtr),
    dataPitch(myPitch) {};

  // --------------------------------------
  // Subscripting operators

  //! Unsafe RHS subscripting operator
  __device__ T operator()( const unsigned int ix,
                           const unsigned int iy,
                           const unsigned int iz ) const
  {
    const char* data = reinterpret_cast<const char*>(this->pitchedPtr);
    // Rows are pitch apart
    const size_t pitch = this->dataPitch;
    // Slices are slicePitch apart
    const size_t slicePitch = pitch * this->dims.y;

    const char* slice = data + ( iz * slicePitch );
    const char* row = slice + ( iy * pitch );

    return( reinterpret_cast<const T*>(row)[ix] );
  }



  //! Unsafe LHS subscripting operator
  __device__ T& operator() ( const unsigned int ix,
                             const unsigned int iy,
                             const unsigned int iz )
  {
    char* data = reinterpret_cast<char*>(this->pitchedPtr);
    const size_t pitch = this->dataPitch;
    const size_t slicePitch = pitch * this->dims.y;

    char* slice = data + ( iz * slicePitch );
    char* row = slice + ( iy * pitch );

    return( reinterpret_cast<T*>(row)[ix] );
  }



  //! Function to return the 1D index of a location
  __device__ unsigned int Index1D ( const unsigned int ix,
                                    const unsigned int iy,
                                    const unsigned int iz ) const
  {
    /*!
      Provides the 1D index of the given location, if the
      current object were to be packed into linear memory
      (i.e. without the padding implied by the pitch).
    */
    unsigned int loc;

    loc = ix + ( this->dims.x * ( iy + ( this->dims.y * iz ) ) );

    return( loc );
  }


  // --------------------------------------

  //! Checks if given co-ordinate is in the volume
  __device__ bool InVolume( const unsigned int ix,
                            const unsigned int iy,
                            const unsigned int iz ) const
  {
    bool res = ( ix < (this->dims.x) );
    res = res && ( iy < (this->dims.y) );
    res = res && ( iz < (this->dims.z) );

    return( res );
  }

  //! Checks is given location is in volume within given tolerance
  __device__ bool InFuzzyVolume( const float3& r,
                                 const float tol ) const
  {
    /*!
      Sometimes we want to know if a point is almost within
      a volume.
      This routine performs that check.
      @param[in] r The location in question
      @param[in] tol The distance outside the volume still considered 'inside'
    */
    bool res = ( (r.x>-tol) && (r.x<(this->dims.x+tol-1)) );
    res = res && (r.y>-tol) && (r.y<(this->dims.y+tol-1));
    res = res && (r.z>-tol) && (r.z<(this->dims.z+tol-1));

    return( res );
  }

private:

  //! Pointer to the allocated memory
  void* const pitchedPtr;
  //! Pitch of the allocated memory
  const size_t dataPitch;
};


// Forward declaration
class CTfactory;

//! Templated class to hold volume data on the GPU
/*!
  This templated class provides a container for volume data.
  It provides a simple interface to padded arrays on the
  device.
  Although it should not (and cannot) be used directly as
  a kernel argument, the VolumeArgGPU is provided for
  this purpose.
  @see VolumeArgGPU
*/
template<typename T>
class VolumeGPU
{
public:

  friend class CTfactory;

  // -------------------------------------------
  // Constructors and destructor

  //! Default constructor
  VolumeGPU( void ) : dims(make_uint3(0,0,0)),
    d_data(make_cudaPitchedPtr(NULL,0,0,0)) {};

  //! Destructor
  ~VolumeGPU( void )
  {
    this->Release();
  }

  // -------------------------------------------
  // Conversion operator

  //! Converts a GPU volume to a kernel argument
  operator VolumeArgGPU<T>( void ) const
  {
    VolumeArgGPU<T> vag( this->dims,
                         this->d_data.ptr,
                         this->d_data.pitch );

    return( vag );
  }

  // -------------------------------------------
  // Data accessors

  //! Return the allocated dimensions
  const dim3 GetDims( void ) const
  {
    return( this->dims );
  }

  //! Return information about the file version
  const char* VersionString( void ) const
  {
    return "$Id: volumegpu.hpp,v 1.29 2012/12/12 21:18:23 nicks Exp $";
  }

  // -------------------------------------------
  // Memory management

  //! Allocates storage of the given dimensions
  void Allocate( const dim3 myDims );


  //! Allocates storage matching dimensions of potentially different type
  template<typename U>
  void Allocate( const VolumeGPU<U>& src )
  {
    /*!
      This routine allocates storage of identical dimensions to the
      given argument, but of potentially different type.
      It is useful if intermediate results have to be in
      higher precision
    */
    this->Allocate( src.GetDims() );
  }


  //! Releases all data on GPU
  void Release( void );


  //! Supplies the size of buffer required on the host
  size_t BufferSize( void ) const;


  //! Allocates a pinned memory buffer on the host
  T* AllocateHostBuffer( void ) const;


  //! Zeros contents of the volume
  void Zero( void );

  //! Copies a volume
  void Copy( const VolumeGPU<T>& src,
             const cudaStream_t myStream = 0 );


  // -------------------------------------------
  // Data transfer

  //! Send a buffer to the GPU
  void SendBuffer( const T* const h_buffer,
                   const cudaStream_t stream = 0 );



  //! Receive a buffer from the GPU
  void RecvBuffer( T* const h_buffer,
                   const cudaStream_t stream = 0 ) const;


  // --------------------------------------
  // Other routines

  //! Method to return number of 'covering' blocks
  dim3 CoverBlocks( const unsigned int threadCount ) const;

  // --------------------------------------

  //! Provides 1D indexing in to contiguous host memory
  unsigned int Index1D( const unsigned int ix,
                        const unsigned int iy,
                        const unsigned int iz ) const;


  //! Converts a 1D index to a dim3 location
  dim3 Index3D( unsigned int idx ) const;

  // ============================================
protected:

  //! Dimensions of the volume
  dim3 dims;
  //! Pointer to the allocated device memory
  cudaPitchedPtr d_data;


  // ============================================

private:
  // -------------------------------------------
  // Prevent copying

  //! Hidden copy constructor
  VolumeGPU( const VolumeGPU& src ) : dims(make_uint3(0,0,0)),
    d_data(make_cudaPitchedPtr(NULL,0,0,0))
  {
    std::cerr << __FUNCTION__
              << ": Please don't use copy constructor"
              << std::endl;
    abort();
  }

  //! Hidden assignment operator
  VolumeGPU& operator=( const VolumeGPU& src )
  {
    std::cerr << __FUNCTION__
              << ": Please don't use assignment operator"
              << std::endl;
    abort();
  }

};


}

}


#include "volumegpu_impl.hpp"

#endif
