/**
 * @file  vecvolgpu.hpp
 * @brief Holds datatype for float vector volume data on the GPU
 *
 * Holds a datatype for vector volume data on the GPU.
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.6 $
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

#include <iostream>

#ifndef VEC_VOL_GPU_CUDA_HPP
#define VEC_VOL_GPU_CUDA_HPP

namespace GPU
{

namespace Classes
{

class VecVolArgGPU
{
public:
  //! Size of the volume
  const dim3 dims;


  // ----------------------------
  // Constructor from inputs
  VecVolArgGPU( const dim3 myDims,
                float* const d_x,
                float* const d_y,
                float* const d_z,
                const size_t mynx ) : dims(myDims),
    x(d_x),
    y(d_y),
    z(d_z),
    nxRow(mynx) {};

  // ----------------------------

  //! Unsafe accessor
  __device__
  float3 operator()( const unsigned int ix,
                     const unsigned int iy,
                     const unsigned int iz )
  {
    const size_t loc = this->PitchIndex1D( ix, iy, iz );
    float3 res;
    res.x = this->x[loc];
    res.y = this->y[loc];
    res.z = this->z[loc];

    return( res );
  }

  //! Unsafe accessor (const)
  __device__
  float3 operator()( const unsigned int ix,
                     const unsigned int iy,
                     const unsigned int iz ) const
  {
    const size_t loc = this->PitchIndex1D( ix, iy, iz );
    float3 res;
    res.x = this->x[loc];
    res.y = this->y[loc];
    res.z = this->z[loc];

    return( res );
  }


  //! Unsafe mutator
  __device__
  void Set( const float3 val,
            const unsigned int ix,
            const unsigned int iy,
            const unsigned int iz )
  {
    const size_t loc = this->PitchIndex1D( ix, iy, iz );
    this->x[loc] = val.x;
    this->y[loc] = val.y;
    this->z[loc] = val.z;
  }

private:
  //! Pointer to x components
  float* const x;
  //! Pointer to y components
  float* const y;
  //! Pointer to z components
  float* const z;
  //! Number of elements allocated per row (pitch/sizeof(float))
  const unsigned int nxRow;


  //! Function to return the 1D index into the pitched memory
  __device__
  size_t PitchIndex1D( const unsigned int ix,
                       const unsigned int iy,
                       const unsigned int iz ) const
  {
    size_t loc;
    loc = ix + ( this->nxRow * ( iy + ( this->dims.y * iz ) ) );
    return( loc );
  }
};


// =============================================================

class VecVolGPU
{
public:
  // --------------------------------------
  // Constructors & destructors

  //! Default constructor
  VecVolGPU( void ) : dims(make_uint3(0,0,0)),
    d_x(make_cudaPitchedPtr(NULL,0,0,0)),
    d_y(make_cudaPitchedPtr(NULL,0,0,0)),
    d_z(make_cudaPitchedPtr(NULL,0,0,0)) {};


  //! Destructor
  ~VecVolGPU( void )
  {
    this->Release();
  }

  // --------------------------------------------
  // Conversion operator

  //! Converts into the kernel argument
  operator VecVolArgGPU( void ) const
  {

    // First do a simple sanity check
    if( (this->d_x.pitch % sizeof(float)) != 0 )
    {
      std::cerr << __FUNCTION__
                << ": Length of row not multiple of float"
                << std::endl;
      abort();
    }

    size_t nxVals = this->d_x.pitch / sizeof(float);

    VecVolArgGPU vvag( this->dims,
                       static_cast<float *>(this->d_x.ptr),
                       static_cast<float *>(this->d_y.ptr),
                       static_cast<float *>(this->d_z.ptr),
                       nxVals );

    return( vvag );
  }


  // -------------------------------------------
  // Data accessors

  //! Return the allocated dimensions
  const dim3 GetDims( void ) const
  {
    return( this->dims );
  }


  // ------------------------------------
  // Memory management

  //! Allocates storage of given dimensions
  void Allocate( const dim3 myDims );

  //! Releases storage
  void Release( void );

  //! Supplies the size of buffer required on the host
  size_t BufferSize( void ) const;

  //! Allocates pinned memory buffers on the host
  void AllocateHostBuffers( float** h_x,
                            float** h_y,
                            float** h_z ) const;


  // -----------------------------------
  // Data transfer

  //! Send buffers to the GPU
  void SendBuffers( const float* const h_x,
                    const float* const h_y,
                    const float* const h_z );

  //! Retrieve buffers from the GPU
  void RecvBuffers( float* const h_x,
                    float* const h_y,
                    float* const h_z ) const;



  //! Copy data from another VecVolGPU
  void Copy( const VecVolGPU& src );


protected:
  //! Dimensions of the volume
  dim3 dims;
  //! Pointer to x data on the device
  cudaPitchedPtr d_x;
  //! Pointer to y data on the device
  cudaPitchedPtr d_y;
  //! Pointer to z data on the device
  cudaPitchedPtr d_z;
private:

  // ---------------------------------------
  // Prevent copying

  //! Hide copy constructor
  VecVolGPU( const VecVolGPU& src ) : dims(),
    d_x(),
    d_y(),
    d_z()
  {
    std::cerr << __FUNCTION__
              << ": Please do not copy"
              << std::endl;
    abort();
    dim3 tmp = src.GetDims();
    tmp.x = 1;
  }

  //! Hide assignment operator
  VecVolGPU& operator=( const VecVolGPU& src )
  {
    std::cerr << __FUNCTION__
              << ": Please do not copy"
              << std::endl;
    abort();
    dim3 tmp = src.GetDims();
    tmp.x = 1;
  }
};

}
}


#endif
