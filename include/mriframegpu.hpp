/**
 * @file  mriframegpu.hpp
 * @brief Holds MRI frame template for the GPU
 *
 * Holds an MRI frame template type for the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/01/20 18:16:16 $
 *    $Revision: 1.6 $
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
 *
 */

#ifndef MRI_FRAME_CUDA_H
#define MRI_FRAME_CUDA_H

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>

#include <cuda_runtime.h>

extern "C" {
#include "mri.h"
}


#include "cudacheck.h"

// ==================================================================


//! Templated function to convert type to MRI->type
template<typename T>
int GetAsMRItype( const T tmp ) {
  return( -1 );
}

template<> int GetAsMRItype<unsigned char>( const unsigned char tmp );
template<> int GetAsMRItype<short>( const short tmp );
template<> int GetAsMRItype<float>( const float tmp );

// ------


// ------

//! Templated memory copy for a row of MRI frame data
template<typename T>
void CopyMRIrowToContiguous( const MRI *src, T* h_slab,
			     const unsigned int iy,
			     const unsigned int iz,
			     const unsigned int iFrame ) {
  std::cerr << __PRETTY_FUNCTION__ << ": Unrecognised type" << std::endl;
  exit( EXIT_FAILURE );
}


template<>
void CopyMRIrowToContiguous<unsigned char>( const MRI* src,
					    unsigned char* h_slab,
					    const unsigned int iy,
					    const unsigned int iz,
					    const unsigned int iFrame );

template<>
void CopyMRIrowToContiguous<short>( const MRI* src,
				    short* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame );


template<>
void CopyMRIrowToContiguous<float>( const MRI* src,
				    float* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame );


//! Templated memory copy for a row of MRI frame data
template<typename T>
void CopyMRIcontiguousToRow( MRI *dst, const
			     T* h_slab,
			     const unsigned int iy,
			     const unsigned int iz,
			     const unsigned int iFrame ) {
  std::cerr << __PRETTY_FUNCTION__ << ": Unrecognised type" << std::endl;
  exit( EXIT_FAILURE );
}


template<>
void CopyMRIcontiguousToRow<unsigned char>( MRI* dst,
					    const unsigned char* h_slab,
					    const unsigned int iy,
					    const unsigned int iz,
					    const unsigned int iFrame );

template<>
void CopyMRIcontiguousToRow<short>( MRI* dst,
				    const short* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame );

template<>
void CopyMRIcontiguousToRow<float>( MRI* dst,
				    const float* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame );


// ================================================================


//! Templated class to hold an MRI frame on the GPU
template<typename T>
class MRIframeGPU {
public:

  //! Original data size
  dim3 cpuDims;
  //! Padded data size
  dim3 gpuDims;
  //! Extent of allocated 3D array
  cudaExtent extent;
  //! Pointer to the allocated memory
  cudaPitchedPtr d_data;
  
  // --------------------------------------------------------
  // Constructors & destructors

  //! Default constructor
  MRIframeGPU( void ) : cpuDims(make_uint3(0,0,0)),
			gpuDims(make_uint3(0,0,0)),
			extent(make_cudaExtent(0,0,0)),
			d_data(make_cudaPitchedPtr(NULL,0,0,0)) {};

 
  //! Destructor
  ~MRIframeGPU( void ) {
    this->Release();
  }

  
  
  // --------------------------------------------------------
  // Memory manipulation
  
  //! Supplies the size of buffer required on the host
  size_t GetBufferSize( void ) const {
    unsigned int nElements;

    nElements = this->gpuDims.x * this->gpuDims.y * this->gpuDims.z;

    return( nElements * sizeof(T) );
  }

  //! Extracts frame dimensions from a given MRI and allocates the memory
  void Allocate( const MRI* src, const unsigned int padSize = 1 ) {
    /*!
      Fills out the cpuDims, gpuDims and extent data members.
      Uses this information to call cudaMalloc3D
      @params[in] src The MRI to use as a template
      @params[in] padSize Each dimension will be padded to be a multiple of this
    */
    
    // Sanity check
    T tmp = 0;
    if( src->type != GetAsMRItype(tmp)  ) {
      std::cerr << __PRETTY_FUNCTION__ << ": MRI type mismatch against " <<
	src->type << std::endl;
      exit( EXIT_FAILURE );
    }

    if( padSize == 0 ) {
      std::cerr << __FUNCTION__ << ": Must have non-zero padSize" << std::endl;
      exit( EXIT_FAILURE );
    }


    // Get rid of old data
    this->Release();

    // Make a note of the dimensions on the CPU
    this->cpuDims.x = src->width;
    this->cpuDims.y = src->height;
    this->cpuDims.z = src->depth;

    // Generate the dimensions on the GPU
    this->gpuDims.x = this->RoundToBlock( cpuDims.x, padSize );
    this->gpuDims.y = this->RoundToBlock( cpuDims.y, padSize );
    this->gpuDims.z = this->RoundToBlock( cpuDims.z, padSize );

    // Make the extent
    this->extent = make_cudaExtent( this->gpuDims.x * sizeof(T),
				    gpuDims.y,
				    gpuDims.z );

    // Allocate the memory
    CUDA_SAFE_CALL( cudaMalloc3D( &(this->d_data), this->extent ) );
  }


  //! Releases memory associated with class instance
  void Release( void ) {
    if( d_data.ptr != NULL ) {
      CUDA_SAFE_CALL( cudaFree( d_data.ptr ) );
    }
    this->cpuDims = make_uint3(0,0,0);
    this->gpuDims = make_uint3(0,0,0);
    this->extent = make_cudaExtent(0,0,0);
    this->d_data = make_cudaPitchedPtr(NULL,0,0,0);
  }


  // --------------------------------------------------------
  // Data transfer

  //! Send the given MRI frame to the GPU
  void Send( const MRI* src, const unsigned int iFrame ) {
    /*!
      Sends the given MRI frame to the GPU.
      For now, this allocates its own memory, and does
      the copy synchronously
    */

    T* h_data;

    // Start with some sanity checks
    this->VerifyMRI( src );

    if( iFrame >= src->nframes ) {
      std:: cerr << __FUNCTION__ << ": Bad frame requested " << iFrame << std::endl;
      exit( EXIT_FAILURE );
    }
    
    const size_t bSize = this->GetBufferSize();

    // Allocate contiguous host memory
    CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_data,
				   bSize,
				   cudaHostAllocDefault ) );
    memset( h_data, 0 , bSize );

    // Extract the data
    this->ExhumeFrame( src, h_data, iFrame );

    // Do the copy
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr = make_cudaPitchedPtr( (void*)h_data,
					     gpuDims.x*sizeof(T),
					     gpuDims.x,
					     gpuDims.y );
    copyParams.dstPtr = this->d_data;
    copyParams.extent = this->extent;
    copyParams.kind = cudaMemcpyHostToDevice;

    CUDA_SAFE_CALL( cudaMemcpy3D( &copyParams ) );

    // Release host memory
    CUDA_SAFE_CALL( cudaFreeHost( h_data ) );
  }



  //! Receives the given MRI frame from the GPU
  void Recv( MRI* dst, const unsigned int iFrame ) const {
    /*!
      Retrieves the given MRI frame from the GPU.
      For now, this allocates its own memory and does the
      copy synchronously
    */

    T* h_data;

    // Start with sanity checks
    this->VerifyMRI( dst );

    if( iFrame >= dst->nframes ) {
      std:: cerr << __FUNCTION__ << ": Bad frame requested " << iFrame << std::endl;
      exit( EXIT_FAILURE );
    }

    const size_t bSize = this->GetBufferSize();

    // Allocate contiguous host memory
    CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_data,
				   bSize,
				   cudaHostAllocDefault ) );

    // Retrieve from GPU
    cudaMemcpy3DParms cpyPrms = {0};
    cpyPrms.srcPtr = this->d_data;
    cpyPrms.dstPtr = make_cudaPitchedPtr( (void*)h_data,
					  gpuDims.x*sizeof(T),
					  gpuDims.x,
					  gpuDims.y );
    cpyPrms.extent = this->extent;
    cpyPrms.kind = cudaMemcpyDeviceToHost;

    CUDA_SAFE_CALL( cudaMemcpy3D( &cpyPrms ) );

    // Retrieve from contiguous RAM
    this->InhumeFrame( dst, h_data, iFrame );

    // Release host memory
    CUDA_SAFE_CALL( cudaFreeHost( h_data ) );
  }


  // ----------------------------------------------------------------------
  //! Method to sanity check MRI
  void VerifyMRI( const MRI* mri ) const {

    T tmp;
    tmp = 0;
    if( mri->type != GetAsMRItype(tmp)  ) {
      std::cerr << __PRETTY_FUNCTION__ << ": MRI type mismatch against " <<
	mri->type << std::endl;
      exit( EXIT_FAILURE );
      // Shut the compiler up
      std::cout << tmp << std::endl;
    }

    if( !this->CheckDims( mri ) ) {
      std::cerr << __PRETTY_FUNCTION__ << ": Size mismatch" << std::endl;
      exit( EXIT_FAILURE );
    }
  }

 
private:

  // ----------------------------------------------------------------------
  // Prevent copying

  //! Copy constructor - don't use
  MRIframeGPU( const MRIframeGPU& src ) : cpuDims(make_uint3(0,0,0)),
					  gpuDims(make_uint3(0,0,0)),
					  extent(make_cudaExtent(0,0,0)),
					  d_data(make_cudaPitchedPtr(NULL,0,0,0)) {
    std::cerr << __PRETTY_FUNCTION__ << ": Please don't use copy constructor" << std::endl;
    exit( EXIT_FAILURE );
  }

  //! Assignment operator - don't use
  MRIframeGPU& operator=( const MRIframeGPU &src ) {
    std::cerr << __PRETTY_FUNCTION__ << ": Please don't use assignment operator" << std::endl;
    exit( EXIT_FAILURE );
  }


  // ----------------------------------------------------------------------

  //! Rounds an array length up to multiple of given padding size
  unsigned int RoundToBlock( const unsigned int arrayLength,
			     const unsigned int padSize ) const {
    /*!
      Rounds the given array length up to the next multiple of
      padSize
    */
    unsigned int nBlocks;
    
    float nBfloat;

    nBfloat = static_cast<float>(arrayLength) / padSize;

    nBlocks = static_cast<unsigned int>( ceilf( nBfloat ) );

    return( nBlocks * padSize );
  }

  // ----------------------------------------------------------------------
  
  //! Function to sanity check dimensions
  bool CheckDims( const MRI* mri ) const {
    
    bool goodDims;

    goodDims = ( mri->width == this->cpuDims.x );
    goodDims = goodDims && ( mri->height == this->cpuDims.y );
    goodDims = goodDims && ( mri->depth == this->cpuDims.z );

    return( goodDims );
  }


  // ----------------------------------------------------------------------
  // Functions for converting MRI frames data to and from contiguous memory
  
  //! Copies a single frame into contiguous memory
  void ExhumeFrame( const MRI* src,
		    T *h_slab,
		    const unsigned int iFrame ) const {
    /*!
      Copies a single MRI frame into contiguous memory on the host.
      Assumes that all the memory has been allocated and GPU dimensions set,
      so everything is not fanatically verified
      @param[in] src The source MRI
      @param[out] h_slab Pointer to the destination array (must be already allocated)
      @param[in] iFrame Which frame to grab
    */
  
    // Start with a few sanity checks
    if( iFrame >= src->nframes ) {
      std::cerr << __FUNCTION__ << ": iFrame out of range" << std::endl;
      exit( EXIT_FAILURE );
    }
    
    if( h_slab == NULL ) {
      std::cerr << __FUNCTION__ << ": h_slab unallocated" << std::endl;
      exit( EXIT_FAILURE );
    }

    this->VerifyMRI( src );
      

    // Main extraction loop
    for( unsigned int iz=0; iz<this->cpuDims.z; iz++ ) {
      for( unsigned int iy=0; iy<this->cpuDims.y; iy++ ) {
	unsigned int iStart = this->Index1D( 0, iy, iz );

	CopyMRIrowToContiguous( src, &( h_slab[iStart] ), iy, iz, iFrame );
      }
    } 
  }


  //! Copies contiguous memory to an MRI frame
  void InhumeFrame( MRI* dst,
		    const T *h_slab,
		    const unsigned int iFrame ) const {
    /*!
      Copies a block of contiguous host memory on the host into an MRI frame.
      Assumes that everything is all properly allocated, so things are not fanatically
      verified
    */
    
    // Start with a few sanity checks
    if( iFrame >= dst->nframes ) {
      std::cerr << __FUNCTION__ << ": iFrame out of range" << std::endl;
      exit( EXIT_FAILURE );
    }
    
    if( h_slab == NULL ) {
      std::cerr << __FUNCTION__ << ": h_slab unallocated" << std::endl;
      exit( EXIT_FAILURE );
    }

    this->VerifyMRI( dst );

    // Main extraction loop
    for( unsigned int iz=0; iz<this->cpuDims.z; iz++ ) {
      for( unsigned int iy=0; iy<this->cpuDims.y; iy++ ) {
	unsigned int iStart = this->Index1D( 0, iy, iz );
	
	CopyMRIcontiguousToRow( dst, &( h_slab[iStart] ), iy, iz, iFrame );
      }
    }
  }


  // ----------------------------------------------------------------------

  // Index into host memory defined by gpuDims
  unsigned int Index1D( const unsigned int ix,
			const unsigned int iy,
			const unsigned int iz ) const {
    return( ix + ( this->gpuDims.x * ( iy + ( this->gpuDims.y * iz ) ) ) );
  }


    

};




 
//! Ancillary class for use in kernel calls
/*!
  This is an auxillary class, for use in actual
  kernel calls.
  The kernel invocation has to be on a call-by-value
  copy, so we want to make sure that the destructor
  doesn't go zapping memory allocations prematurely
*/
template<typename T>
class MRIframeOnGPU {
public:
  //! Padded data size
  dim3 dims;
  //! Extent of allocated 3D array
  cudaExtent extent;
  //! Pointer to the allocated memory
  cudaPitchedPtr data;
  // --------------------------------------------------------
  // Constructors
  
  //! Default constructor
  MRIframeOnGPU( void ) : dims(make_uint3(0,0,0)),
		  extent(make_cudaExtent(0,0,0)),
		  data(make_cudaPitchedPtr(NULL,0,0,0)) {};
  
  //! Constructor from MRIframeGPU
  MRIframeOnGPU( const MRIframeGPU<T>& src ) : dims(src.gpuDims),
				       extent(src.extent),
				       data(src.d_data) {};
  
  
  // --------------------------------------------------------
  // Subscripting operators
  
  __device__ T operator() ( const unsigned int ix,
			    const unsigned int iy,
			    const unsigned int iz ) const {
    const char* data = reinterpret_cast<const char*>(this->data.ptr);
    // Rows are pitch apart
    size_t pitch = this->data.pitch;
    // Slices are slicePitch apart
    size_t slicePitch = pitch * this->extent.height;
    
    const char* slice = data + ( iz * slicePitch );
    const char* row = slice + ( iy * pitch );
    
      return( reinterpret_cast<const T*>(row)[ix] );
  }
  
  
  __device__ T& operator() ( const unsigned int ix,
			     const unsigned int iy,
			     const unsigned int iz ) {
    char* data = reinterpret_cast<char*>(this->data.ptr);
    size_t pitch = this->data.pitch;
    size_t slicePitch = pitch * this->extent.height;
    
    char* slice = data + ( iz * slicePitch );
    char* row = slice + ( iy * pitch );
    
    return( reinterpret_cast<T*>(row)[ix] );
  }
  
    
  //! Clamps input integer into range
  __device__ unsigned int ClampCoord( const int i,
				      const unsigned int iMax ) const {
    /*!
      Performs clamping on the input index.
      Negative values are set to 0, values greater than iMax-1 are
      set to iMax-1
    */
    if( i < 0 ) {
      return 0;
    } else if( i > (iMax-1) ) {
      return( iMax-1 );
    } else {
	return i;
    }
  }
};




#endif
