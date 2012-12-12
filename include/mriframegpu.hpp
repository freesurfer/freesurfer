/**
 * @file  mriframegpu.hpp
 * @brief Holds MRI frame template for the GPU
 *
 * Holds an MRI frame template type for the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.49 $
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

#ifndef MRI_FRAME_CUDA_H
#define MRI_FRAME_CUDA_H

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>

#include <cuda_runtime.h>

#include "mri.h"

#include "volumegpu.hpp"
#include "affinegpu.hpp"

#include "cudatypeutils.hpp"
#include "cudacheck.h"

// ==================================================================

namespace GPU
{

namespace Classes
{


//! Ancillary class for use in kernel calls
/*!
  This is an auxillary class, for use in actual
  kernel calls.
  Host side functions should simply use the MRIframeGPU
  class.
  The kernel invocation has to be on a call-by-value
  copy, so we want to make sure that the destructor
  doesn't go zapping memory allocations prematurely.
  Most of its functionality is provided by the base
  class VolumeArgGPU.

  To use in a kernel, simply declare the kernel
  arguments to be of this type:
  \code
  template<typename T, typename U>
  __global__ void MyKernel( const MRIframeOnGPU<T> a,
                            MRIframeOnGPU<U> b ) {
...
  }
  \endcode
  Since appropriate conversion operators are provided,
  the kernel can be invoked simply with MRIframeGPU
  arguments:
  \code
  MRIframeGPU<T> a;
  MRIframeGPU<U> b;
  MyKernel<T,U><<<grid,threads>>>( a, b );
  \endcode
  The necessary casts and copies will occur automatically.
  @see VolumeArgGPU
  @see MRIframeGPU
*/
template<typename T>
class MRIframeOnGPU : public VolumeArgGPU<T>
{
public:

  // --------------------------------------------------------
  // Constructors

  //! Default constructor
  MRIframeOnGPU<T>( void ) : VolumeArgGPU<T>() {};


  //! Constructor from MRIframeGPU
  MRIframeOnGPU<T>( const VolumeArgGPU<T>& src ) :
    VolumeArgGPU<T>(src) {};


  // --------------------------------------------------------

  //! Out of range value
  static const T kOutOfRangeVal = 0;

  // --------------------------------------------------------
  // Subscripting operators

  //! RHS Subscripting operator
  __device__ T operator() ( const unsigned int ix,
                            const unsigned int iy,
                            const unsigned int iz ) const
  {
    /*!
      This version of the subscripting operator is safe.
      It overrides the unsafe version in VolumeArgGPU
      It bounds checks the requested dimensions,
      and returns kOutofRangeVal if out of range
    */
    if( this->InVolume( ix, iy, iz ) )
    {
      // Use the unsafe operator() of the base class
      return( this->VolumeArgGPU<T>::operator()( ix, iy, iz ) );
    }
    else
    {
      return( this->kOutOfRangeVal );
    }
  }


  //! LHS subscripting operator
  __device__ T& operator() ( const unsigned int ix,
                             const unsigned int iy,
                             const unsigned int iz )
  {
    /*!
      This is a thin wrapper around the base class operator.
      For reasons which are probably lost in the depths of the
      C++ standard (to do with resolving overloaded and overridden
      methods), this is necessary to use subscripting on the LHS
      of an expression
    */
    return( this->VolumeArgGPU<T>::operator()( ix, iy, iz ) );
  }

  //! Utility function to convert float to the class' datatype
  inline
  __device__ T ConvertFloat( const float in ) const
  {
    /*!
      This template is specialised for each supported class.
      The unspecialised default will write a constant
      negative value
    */
    return( -1 );
  }

  // --------------------------------------------------------

  //! Clamps input integer into range
  __device__ unsigned int ClampCoord( const int i,
                                      const unsigned int iMax ) const
  {
    /*!
      Performs clamping on the input index.
      Negative values are set to 0, values greater than iMax-1 are
      set to iMax-1
    */
    if( i < 0 )
    {
      return 0;
    }
    else if( i > (iMax-1) )
    {
      return( iMax-1 );
    }
    else
    {
      return i;
    }
  }
};

// ================================================================


//! Templated class to hold an MRI frame on the GPU
/*!
  This class is provided to host code to hold one frame
  of MRI data on the GPU.
  The data are held in linear memory, which is padded to
  optimise CUDA memory accesses.
  The data can also be copied into a CUDA array, in order
  to enable texturing.
  Most of the device-side functionality is provided by
  the VolumeGPU base class.
  On the host side, this derived class provides easy methods
  of copying data to the GPU.
  This class is templated based on the stored data type.
  As a result, switch blocks (and specialised templates) are
  often necessay to handle incoming MRI structures.
  @see MRIframeOnGPU
*/
template<typename T>
class MRIframeGPU : public VolumeGPU<T>
{
public:

  // --------------------------------------------------------
  // Constructors & destructors

  //! Default constructor does nothing
  MRIframeGPU<T>( void ) : VolumeGPU<T>(),
    thick(0),
    sizes(make_float3(0,0,0)),
    d_i_to_r(),
    d_r_to_i() {};

  //! Conversion operator to MRIframeOnGPU
  operator MRIframeOnGPU<T>( void ) const
  {
    VolumeArgGPU<T> vag( *this );

    return( MRIframeOnGPU<T>( vag ) );
  }

  // --------------------------------------------------------
  // Public data (should probably be private)

  AffineTransformation d_i_to_r;
  AffineTransformation d_r_to_i;

  // --------------------------------------------------------
  // Data accessors

  //! Return information about the file version
  const char* VersionString( void ) const
  {
    return "$Id: mriframegpu.hpp,v 1.49 2012/12/12 21:18:23 nicks Exp $";
  }

  //! Return the 'thick' field
  float GetThickness( void ) const
  {
    return( this->thick );
  }

  //! Return the 'sizes' field
  float3 GetSizes( void ) const
  {
    return( this->sizes );
  }

  // --------------------------------------------------------
  // Memory manipulation

  // -----

  //! Extracts frame dimensions from a given MRI and allocates the memory
  void Allocate( const MRI* src )
  {
    /*!
      Extracts frame dimensions from the supplied MRI, and uses
      this to allocate the appropriate amount of device
      memory.
      It first performs a sanity check on the type field
      of the MRI structure
      @param[in] src The MRI to use as a template
    */

    // Sanity checks
    if( src->type != this->MRItype()  )
    {
      std::cerr << __PRETTY_FUNCTION__
                << ": MRI type mismatch against "
                << src->type << std::endl;
      abort();
    }


    dim3 myDims = make_uint3( src->width, src->height, src->depth );

    this->VolumeGPU<T>::Allocate( myDims );
  }

  // -----

  //! Allocates matching storage of possibly different type
  template<typename U>
  void Allocate( const MRIframeGPU<U>& src )
  {
    /*!
      Some algorithms require intermediate results to be
      stored in higher precision.
      These will require identically dimensioned MRIframeGPU
      classes, but of a different type.
      This method clones the dimensions of the given source
      frame, but allocates space for its own type.
    */
    this->VolumeGPU<T>::Allocate( src );
  }

  //! Pass through wrapper to base class
  void Allocate( const dim3 dims )
  {
    this->VolumeGPU<T>::Allocate( dims );
  }


  // --------------------------------------------------------
  // Data transfer

  //! Sends the given MRI frame and the scalar data
  void Send( const MRI* src,
             const unsigned int iFrame,
             void* const h_work = NULL,
             const cudaStream_t stream = 0 )
  {

    // Copy the scalars over
    this->thick = src->thick;

    this->sizes = make_float3( src->xsize,
                               src->ysize,
                               src->zsize );

    // Send the transforms
    this->d_r_to_i.SetTransform( src->r_to_i__ );
    this->d_i_to_r.SetTransform( src->i_to_r__ );

    this->SendFrame( src, iFrame, h_work, stream );
  }

  //! Send the given MRI frame to the GPU
  void SendFrame( const MRI* src,
                  const unsigned int iFrame,
                  void* const h_work = NULL,
                  const cudaStream_t stream = 0 )
  {
    /*!
      Sends the given MRI frame to the GPU.
      Optional arguments can be used to supply page-locked
      host memory for the transfer, and a stream in which
      to perform the transfer.
      If supplied, the array h_work must be at least
      GetBufferSize() bytes long.
      Furthermore, the calling routine is responsible
      for synchronisation
      @param[in] src The MRI we wish to send to the GPU
      @param[in] iFrame Which frame of the MRI we wish to send
      @param h_work Optional pagelocked host memory
      @param[in] stream Optional CUDA stream for the copy
    */

    T* h_data;

    // Start with some sanity checks
    this->VerifyMRI( src );

    if( iFrame >= static_cast<unsigned int>(src->nframes) )
    {
      std:: cerr << __FUNCTION__
                 << ": Bad frame requested " << iFrame
                 << std::endl;
      abort();
    }

    // See if we need to allocate workspace
    const size_t bSize = this->BufferSize();
    // See if we were supplied with workspace
    if( h_work != NULL )
    {
      h_data = reinterpret_cast<T*>(h_work);
    }
    else
    {
      // Allocate contiguous host memory
      CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_data,
                                     bSize,
                                     cudaHostAllocDefault ) );
    }


    // Extract the data
    this->ExhumeFrame( src, h_data, iFrame );

    // Do the copy
    this->SendBuffer( h_data, stream );

    // Release host memory if needed
    if( h_work == NULL )
    {
      CUDA_SAFE_CALL( cudaStreamSynchronize( stream ) );
      CUDA_SAFE_CALL( cudaFreeHost( h_data ) );
    }
  }

  // -----

  //! Receives the given MRI frame and the scalar data
  void Recv( MRI* dst,
             const unsigned int iFrame,
             void* const h_work = NULL,
             const cudaStream_t stream = 0 ) const
  {

    // Retrieve the scalars
    dst->thick = this->thick;
    dst->xsize = this->sizes.x;
    dst->ysize = this->sizes.y;
    dst->zsize = this->sizes.z;

    this->RecvFrame( dst, iFrame, h_work, stream );
  }

  //! Receives the given MRI frame from the GPU
  void RecvFrame( MRI* dst,
                  const unsigned int iFrame,
                  void* const h_work = NULL,
                  const cudaStream_t stream = 0 ) const
  {
    /*!
      Retrieves the given MRI frame from the GPU.
      Optional arguments can be used to supply page-locked
      host memory for the transfer, and a stream in which
      to perform the transfer.
      If supplied, the array h_work must be at least
      GetBufferSize() bytes long
    */

    T* h_data;

    // Start with sanity checks
    this->VerifyMRI( dst );

    if( iFrame >= static_cast<unsigned int>(dst->nframes) )
    {
      std:: cerr << __FUNCTION__
                 << ": Bad frame requested " << iFrame
                 << std::endl;
      abort();
    }



    // See about buffers
    const size_t bSize = this->BufferSize();

    // Allocate contiguous host memory if needed
    if( h_work != NULL )
    {
      h_data = reinterpret_cast<T*>(h_work);
    }
    else
    {
      CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_data,
                                     bSize,
                                     cudaHostAllocDefault ) );
    }

    // Retrieve from GPU
    this->RecvBuffer( h_data, stream );

    CUDA_SAFE_CALL( cudaStreamSynchronize( stream ) );

    // Retrieve from contiguous RAM
    this->InhumeFrame( dst, h_data, iFrame );

    // Release host memory
    if( h_work == NULL )
    {
      CUDA_SAFE_CALL( cudaFreeHost( h_data ) );
    }
  }


  // -----

  // ----------------------------------------------------------------------
  //! Method to sanity check MRI
  void VerifyMRI( const MRI* mri ) const
  {
    /*!
      A very simple routine to check that both the type
      and dimensions of the current object match those
      of the given MRI structure.
    */

    if( mri->type != this->MRItype()  )
    {
      std::cerr << __PRETTY_FUNCTION__
                << ": MRI type mismatch against "
                << mri->type << std::endl;
      abort();
    }

    if( !this->CheckDims( mri ) )
    {
      std::cerr << __PRETTY_FUNCTION__
                << ": Size mismatch"
                << std::endl;
      abort();
    }
  }

  //! Method to return MRI type (has specialisations below)
  int MRItype( void ) const
  {
    /*!
      Returns the type of data held in a form that the
      rest of the Freesurfer package will recognise.
      This is implemented as a set of specialised templates.
    */
    return(-1);
  }


private:
  // --------------------------------------------------------------------

  //! Mirror of the 'thick' field of an MRI
  float thick;

  //! Mirror of the {x|y|z}size fields (the voxel sizes)
  float3 sizes;

  //! Function to sanity check dimensions
  bool CheckDims( const MRI* mri ) const
  {

    const dim3 mriDims = make_uint3( mri->width,
                                     mri->height,
                                     mri->depth );

    return( mriDims == this->dims );
  }


  // ----------------------------------------------------------------------
  // Functions for converting MRI frames data to and from contiguous memory

  //! Copies a single frame into contiguous memory
  void ExhumeFrame( const MRI* src,
                    T *h_slab,
                    const unsigned int iFrame ) const
  {
    /*!
      Copies a single MRI frame into contiguous memory on the host.
      Assumes that all the memory has been allocated and GPU dimensions set,
      so everything is not fanatically verified.
      The name of this method should not be taken as a value judgement
      of the CPU-side datastructure.
      @param[in] src The source MRI
      @param[out] h_slab Pointer to the destination array (must be already allocated)
      @param[in] iFrame Which frame to grab
    */

    // Start with a few sanity checks
    if( iFrame >= static_cast<unsigned int>(src->nframes) )
    {
      std::cerr << __FUNCTION__
                << ": iFrame out of range"
                << std::endl;
      abort();
    }

    if( h_slab == NULL )
    {
      std::cerr << __FUNCTION__
                << ": h_slab unallocated"
                << std::endl;
      abort();
    }

    this->VerifyMRI( src );


    // Main extraction loop
    for( unsigned int iz=0; iz<this->dims.z; iz++ )
    {
      for( unsigned int iy=0; iy<this->dims.y; iy++ )
      {
        unsigned int iStart = this->Index1D( 0, iy, iz );

        this->ExhumeRow( src, &( h_slab[iStart] ), iy, iz, iFrame );
      }
    }
  }

  // -----

  //! Copies contiguous memory to an MRI frame
  void InhumeFrame( MRI* dst,
                    const T *h_slab,
                    const unsigned int iFrame ) const
  {
    /*!
      Copies a block of contiguous host memory on the host
      into an MRI frame.
      Assumes that everything is all properly allocated,
      so things are not fanatically verified
    */

    // Start with a few sanity checks
    if( iFrame >= static_cast<unsigned int>(dst->nframes) )
    {
      std::cerr << __FUNCTION__
                << ": iFrame out of range"
                << std::endl;
      abort();
    }

    if( h_slab == NULL )
    {
      std::cerr << __FUNCTION__
                << ": h_slab unallocated"
                << std::endl;
      abort();
    }

    this->VerifyMRI( dst );

    // Main extraction loop
    for( unsigned int iz=0; iz<this->dims.z; iz++ )
    {
      for( unsigned int iy=0; iy<this->dims.y; iy++ )
      {
        unsigned int iStart = this->Index1D( 0, iy, iz );

        this->InhumeRow( dst, &( h_slab[iStart] ), iy, iz, iFrame );
      }
    }
  }

  // -----

  //! Wrapper around memcpy for MRI->contiguous transfers
  void ExhumeRow( const MRI* src,
                  T *h_slab,
                  const unsigned int iy,
                  const unsigned int iz,
                  const unsigned int iFrame ) const
  {
    /*!
      Specialised method to copy a row of an MRI frame
      into contiguous memory.
      The row of x voxels are specified by iy, iz and iFrame.
      This default method aborts the program
    */

    std::cerr << __PRETTY_FUNCTION__
              << ": Unrecognised data type "
              << src->type << std::endl;
    abort();
  }

  // -----

  //! Wrapper around memcpy for contiguous->MRI transfers
  void InhumeRow( MRI* dst,
                  const T *h_slab,
                  const unsigned int iy,
                  const unsigned int iz,
                  const unsigned int iFrame ) const
  {
    /*!
      Specialised method to copy a chunk of contiguous
      memory into a row of an MRI frame.
      The row of x voxels are specified by iy, iz and iFrame.
      This default method aborts the program
    */

    std::cerr << __PRETTY_FUNCTION__
              << ": Unrecognised data type "
              << dst->type << std::endl;
    abort();
  }


};


// Declarations of specialised methods

template<> int MRIframeGPU<unsigned char>::MRItype( void ) const;
template<> int MRIframeGPU<short>::MRItype( void ) const;
template<> int MRIframeGPU<float>::MRItype( void ) const;
template<> int MRIframeGPU<int>::MRItype( void ) const;



template<>
void MRIframeGPU<unsigned char>::ExhumeRow( const MRI* src,
    unsigned char* h_slab,
    const unsigned int iy,
    const unsigned int iz,
    const unsigned int iFrame ) const;
template<>
void MRIframeGPU<short>::ExhumeRow( const MRI* src,
                                    short* h_slab,
                                    const unsigned int iy,
                                    const unsigned int iz,
                                    const unsigned int iFrame ) const;
template<>
void MRIframeGPU<float>::ExhumeRow( const MRI* src,
                                    float* h_slab,
                                    const unsigned int iy,
                                    const unsigned int iz,
                                    const unsigned int iFrame ) const;

template<>
void MRIframeGPU<int>::ExhumeRow( const MRI* src,
                                  int* h_slab,
                                  const unsigned int iy,
                                  const unsigned int iz,
                                  const unsigned int iFrame ) const;



template<>
void MRIframeGPU<unsigned char>::InhumeRow( MRI* dst,
    const unsigned char* h_slab,
    const unsigned int iy,
    const unsigned int iz,
    const unsigned int iFrame ) const;
template<>
void MRIframeGPU<short>::InhumeRow( MRI* dst,
                                    const short* h_slab,
                                    const unsigned int iy,
                                    const unsigned int iz,
                                    const unsigned int iFrame ) const;

template<>
void MRIframeGPU<float>::InhumeRow( MRI* dst,
                                    const float* h_slab,
                                    const unsigned int iy,
                                    const unsigned int iz,
                                    const unsigned int iFrame ) const;

template<>
void MRIframeGPU<int>::InhumeRow( MRI* dst,
                                  const int* h_slab,
                                  const unsigned int iy,
                                  const unsigned int iz,
                                  const unsigned int iFrame ) const;

// ======================================================================



template<> __device__
inline
unsigned char MRIframeOnGPU<unsigned char>::ConvertFloat( const float in ) const
{
  // Copy 'nint' from utils.c
  return( static_cast<unsigned char>( in+0.5f ) );
}

template<> __device__
inline
short MRIframeOnGPU<short>::ConvertFloat( const float in ) const
{
  // Copy 'nint' from utils.c
  return( static_cast<short>( in+0.5f ) );
}

template<> __device__
inline
float MRIframeOnGPU<float>::ConvertFloat( const float in ) const
{
  return( in );
}

}
}


#endif
