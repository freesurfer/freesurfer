/**
 * @file  volumecpu.hpp
 * @brief Holds templated datatype for volume data on the CPU
 *
 * Holds a templated datatype for volume data on the CPU.
 * It is suspiciously similar to VolumeGPU - indeed it exists
 * because I can't come up with a simple method of implementing
 * a portion of mri_ca_register on the GPU (gcamLabelTerm, to be exact)
 * and I don't want a full repacking penalty.
 * You may well complain that this is horribly inelegant....
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.8 $
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
#include <cstring>
#include <iostream>

#include <cuda_runtime.h>


#include "volumegpu.hpp"

#ifndef VOLUME_CPU_HPP
#define VOLUME_CPU_HPP

namespace Freesurfer
{

//! Templated "mutator" class for the CPU
template<typename T>
class VolumeArgCPU
{
public:

  //! Volume sizes
  const unsigned int nx, ny, nz;

  // ------------------------------------------

  //! Default constructor
  VolumeArgCPU( void ) : nx(0), ny(0), nz(0),
    data(NULL) {};

  //! Constructor from inputs
  VolumeArgCPU( const unsigned int _nx,
                const unsigned int _ny,
                const unsigned int _nz,
                T* const _data ) : nx(_nx), ny(_ny), nz(_nz),
    data(_data) {};

  //! Empty destructor
  ~VolumeArgCPU( void ) {};

  // ---------------------------------------------------
  // Subscripting operators (unsafe)

  //! RHS subscripting operator
  inline const T& operator()( const unsigned int ix,
                              const unsigned int iy,
                              const unsigned int iz ) const
  {
    return( data[ix + ( nx * ( iy + (ny*iz) ) )] );
  }

  //! LHS subscripting operator
  inline T& operator()( const unsigned int ix,
                        const unsigned int iy,
                        const unsigned int iz )
  {
    return( data[ix + ( nx * ( iy + (ny*iz) ) )] );
  }

private:
  //! Data pointer
  T* const data;

};



//! Templated management class for the CPU
template<typename T>
class VolumeCPU
{
public:

  //! Default constructor
  VolumeCPU( void ) : nx(0), ny(0), nz(0),
    data(NULL) {};

  //! Destructor
  ~VolumeCPU( void )
  {
    this->Release();
  }

  //! Conversion operator
  operator VolumeArgCPU<T>( void ) const
  {
    VolumeArgCPU<T> vac( this->nx, this->ny, this->nz,
                         this->data );

    return( vac );
  }

  // --------------------------------
  // Memory management

  //! Allocates storage of given dimensions
  void Allocate( const unsigned int xn,
                 const unsigned int yn,
                 const unsigned int zn )
  {
    // Check for prior allocation
    if( this->data != NULL )
    {
      if( (xn == this->nx) &&
          (yn == this->ny) &&
          (zn == this->nz) )
      {
        // Same size
        return;
      }
      else
      {
        this->Release();
      }
    }

    this->nx = xn;
    this->ny = yn;
    this->nz = zn;

    CUDA_SAFE_CALL( cudaHostAlloc( &(this->data),
                                   this->BufferSize(),
                                   cudaHostAllocDefault ) );
  }

  //! Releases memory
  void Release( void )
  {
    if( this->data != NULL )
    {
      this->nx = 0;
      this->ny = 0;
      this->nz = 0;
      CUDA_SAFE_CALL( cudaFreeHost( this->data ) );
      this->data = NULL;
    }
  }

  //! Size of the buffer
  size_t BufferSize( void ) const
  {
    size_t mem;
    mem = this->nx * this->ny * this->nz * sizeof(T);
    return( mem );
  }

  //! Accessor for dimensions
  void GetDims( unsigned int& x, unsigned int& y, unsigned int& z ) const
  {
    x = this->nx;
    y = this->ny;
    z = this->nz;
  }


  //! Check for dimension match
  template<typename U>
  bool MatchDims( const VolumeCPU<U>& src ) const
  {
    bool match;

    unsigned int srcx, srcy, srcz;
    src.GetDims( srcx, srcy, srcz );

    match = ( this->nx == srcx );
    match = match && ( this->ny == srcy );
    match = match && ( this->nz == srcz );

    return( match );
  }


  //! Copy method
  void Copy( const VolumeCPU<T>& src )
  {

    // Check for NULL source
    if( src.data == NULL )
    {
      std::cerr << __FUNCTION__
                << ": Can't copy from NULL source"
                << std::endl;
      abort();
    }

    // Check for self-copy
    if( &src == this )
    {
      return;
    }

    this->Allocate( src.nx, src.ny, src.nz );

    memcpy( this->data, src.data, src.BufferSize() );
  }


  // -------------------------------
  // Data transfer

  //! Pulls from VolumeGPU
  void PullGPU( const GPU::Classes::VolumeGPU<T>& src )
  {
    src.RecvBuffer( this->data );
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
  }

  //! Returns data to a VolumeGPU
  void PushGPU( GPU::Classes::VolumeGPU<T>& dst ) const
  {
    dst.SendBuffer( this->data );
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
  }

  // ---------------------------------
  // Value manipulation

  //! Replace all occurences of one value by another
  void ReplaceValue( const T inVal, const T outVal )
  {


    VolumeArgCPU<T> me( *this );

    for( unsigned int iz=0; iz<this->nz; iz++ )
    {
      for( unsigned int iy=0; iy<this->ny; iy++ )
      {
        T *rowPtr = &me(0,iy,iz);

        for( unsigned int ix=0; ix<this->nx; ix++ )
        {
          if( inVal == *rowPtr )
          {
            *rowPtr = outVal;
          }
          rowPtr++;
        }
      }
    }
  }

private:
  unsigned int nx, ny, nz;
  T* data;

  //! Hidden copy constructor
  VolumeCPU( const VolumeCPU& src ) : nx(0), ny(0), nz(0),
    data(NULL)
  {
    std::cerr << __FUNCTION__
              << ": Please don't use copy constructor"
              << std::endl;
    abort();
  }

  //! Hidden assignment operator
  VolumeCPU& operator=( const VolumeCPU& src )
  {
    std::cerr << __FUNCTION__
              << ": Please don't use assignment operator"
              << std::endl;
    abort();
  }

};


}


#endif
