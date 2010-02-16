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
 *    $Author: rge21 $
 *    $Date: 2010/02/16 16:58:24 $
 *    $Revision: 1.5 $
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

#include <cstdlib>

#include <iostream>

#include "cudacheck.h"
#include "cudatypeutils.hpp"

#ifndef VOLUME_GPU_CUDA_H
#define VOLUME_GPU_CUDA_H

namespace GPU {

  namespace Classes {

    //! Templated ancilliary class for use in kernel calls
    template<typename T>
    class VolumeArgGPU {
    public:

      //! Out of range value
      static const T kOutOfRangeVal = 0;

      //! Padded data size
      dim3 dims;
      //! Extent of allocated 3D array
      cudaExtent extent;
      //! Pointer to the allocated memory
      cudaPitchedPtr data;

      // --------------------------------------
      // Constructors
      VolumeArgGPU( void ) : dims(make_uint3(0,0,0)),
			     extent(make_cudaExtent(0,0,0)),
			     data(make_cudaPitchedPtr(NULL,0,0,0)) {};
      
      // --------------------------------------
      // Subscripting operators
      
      //! Unsafe RHS subscripting routine
      __device__ T UnsafeAt( const unsigned int ix,
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

      //! RHS Subscripting operator
      __device__ T operator() ( const unsigned int ix,
				const unsigned int iy,
				const unsigned int iz ) const {
	/*!
	  This version of the subscripting operator is safe.
	  It bounds checks the requested dimensions,
	  and returns kOutofRangeVal if out of range
	*/
	if( this->InVolume( ix, iy, iz ) ) {
	  return( this->UnsafeAt( ix, iy, iz ) );
	} else {
	  return( this->kOutOfRangeVal );
	}
      }
      
      
      //! LHS subscripting operator
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


      // --------------------------------------
      
      //! Checks if given co-ordinate is in the volume
      __device__ bool InVolume( const unsigned int ix,
				const unsigned int iy,
				const unsigned int iz ) const {
	bool res = ( ix < (this->dims.x) );
	res = res && ( iy < (this->dims.y) );
	res = res && ( iz < (this->dims.z) );

	return( res );
      }

      //! Checks is given location is in volume within given tolerance
      __device__ bool InFuzzyVolume( const float3& r,
				     const float tol ) const {
	bool res = ( (r.x>-tol) && (r.x<(this->dims.x+tol)) );
	res = res && (r.y>-tol) && (r.y<(this->dims.y+tol));
	res = res && (r.z>-tol) && (r.z<(this->dims.z+tol));

	return( res );
      }
    };



    //! Templated class to hold volume data on the GPU
    template<typename T>
    class VolumeGPU {
    public:

      // -------------------------------------------
      // Constructors and destructor
      
      //! Default constructor
      VolumeGPU( void ) : dims(make_uint3(0,0,0)),
			  extent(make_cudaExtent(0,0,0)),
			  d_data(make_cudaPitchedPtr(NULL,0,0,0)),
			  dca_data(NULL) {};

      //! Destructor
      ~VolumeGPU( void ) {
	this->Release();
      }

      // -------------------------------------------
      // Conversion operator

      //! Converts volume to a kernel argument
      operator VolumeArgGPU<T>( void ) {
	VolumeArgGPU<T> vag;

	vag.dims = this->dims;
	vag.extent = this->extent;
	vag.data = this->d_data;
	
	return( vag );
      }

      // -------------------------------------------
      // Data accessors
      
      //! Return the allocated dimensions
      const dim3 GetDims( void ) const {
	return( this->dims );
      }

      //! Return information about the file version
      const char* VersionString( void ) const {
	return "$Id: volumegpu.hpp,v 1.5 2010/02/16 16:58:24 rge21 Exp $";
      }
      
      //! Return pointer to the cudaArray
      const cudaArray* GetArray( void ) const {
	return( this->dca_data );
      }

      //! Checks to see if CUDA array is allocated
      bool HasArray( void ) const {
	return( this->dca_data != NULL );
      }

      // -------------------------------------------
      // Memory management

      //! Allocates storage of the given dimensions
      void Allocate( const dim3 myDims ) {
	/*!
	  Allocates GPU memory to hold a volume of the
	  given size.
	  This memory may include padding for GPU
	  memory alignment requirements
	  @params[in] myDims The dimensions required
	*/
	
	// Check if we can re-use current memory
	if( myDims == this->dims ) {
	  return;
	}

	// Release existing memory
	this->Release();

	// Set the dimensions
	this->dims = myDims;


	// Make the extent
	this->extent = make_cudaExtent( this->dims.x * sizeof(T),
					this->dims.y,
					this->dims.z );
	
	// Allocate the memory
	CUDA_SAFE_CALL( cudaMalloc3D( &(this->d_data), this->extent ) );
      }

      //! Allocates storage matching dimensions of potentially different type
      template<typename U>
      void Allocate( const VolumeGPU<U>& src ) {

	this->Allocate( src.GetDims() );
      }


      //! Allocates the CUDA array member
      void AllocateArray( void ) {
	/*!
	  Allocates the CUDA array member based on the existing
	  data stored in the class.
	  Checks to see that the sizes are non-zero by
	  examining the d_data pointer
	*/

	// Check for initialisation
	if( this->d_data.ptr == NULL ) {
	  std::cerr << __FUNCTION__
		    << ": d_data is NULL!"
		    << std::endl;
	  exit( EXIT_FAILURE );
	}

	// See if array is already allocated
	if( this->HasArray() ) {
	  /*
	    We don't need to do anything.
	    The current dimensions can only be set
	    via the Allocate method, and this will
	    release everything is there's a mismatch
	  */
	  return;
	}

	cudaChannelFormatDesc cd = cudaCreateChannelDesc<T>();

	CUDA_SAFE_CALL( cudaMalloc3DArray( &(this->dca_data),
					   &cd,
					   this->extent ) );
      }


      //! Releases all data on GPU
      void Release( void ) {
	this->ReleaseArray();
	if( this->d_data.ptr != NULL ) {
	  CUDA_SAFE_CALL( cudaFree( d_data.ptr ) );
	}
	this->dims = make_uint3(0,0,0);
	this->extent = make_cudaExtent(0,0,0);
	this->d_data = make_cudaPitchedPtr(NULL,0,0,0);
      }


      //! Releases CUDA array member
      void ReleaseArray( void ) {
	if( this->dca_data != NULL ) {
	  CUDA_SAFE_CALL( cudaFreeArray( this->dca_data ) );
	  this->dca_data = NULL;
	}
      }

      //! Supplies the size of buffer required on the host
      size_t BufferSize( void ) const {
	unsigned int nElements;

	nElements = this->dims.x * this->dims.y * this->dims.z;

	return( nElements * sizeof(T) );
      }

      // -------------------------------------------
      // Data transfer

      //! Send a buffer to the GPU
      void SendBuffer( const T* const h_buffer,
		       const cudaStream_t stream = 0 ) {
	/*!
	  Sends the given contiguous buffer on the
	  host to the GPU.
	  The data in the buffer must be ordered by
	  the Index1D method.
	  The copy is done asynchronously in the optionally
	  defined stream.
	  Accordingly, this routine may return before the
	  copy is complete - the user is responsible for
	  stream management.
	*/

	cudaMemcpy3DParms copyParams = {0};

	copyParams.srcPtr = make_cudaPitchedPtr( (void*)h_buffer,
						 this->dims.x*sizeof(T),
						 this->dims.x,
						 this->dims.y );
	copyParams.dstPtr = this->d_data;
	copyParams.extent = this->extent;
	copyParams.kind = cudaMemcpyHostToDevice;

	CUDA_SAFE_CALL_ASYNC( cudaMemcpy3DAsync( &copyParams, stream ) );
      }



      //! Receive a buffer from the GPU
      void RecvBuffer( T* const h_buffer,
		       const cudaStream_t stream = 0 ) const {
	/*!
	  Retrieves data from the GPU into the
	  contiguous buffer on the host.
	  The copy is done asynchronously, so the
	  user is responsible for stream management
	*/
	
	cudaMemcpy3DParms cpyPrms = {0};
	cpyPrms.srcPtr = this->d_data;
	cpyPrms.dstPtr = make_cudaPitchedPtr( (void*)h_buffer,
					      this->dims.x*sizeof(T),
					      this->dims.x,
					      this->dims.y );
	cpyPrms.extent = this->extent;
	cpyPrms.kind = cudaMemcpyDeviceToHost;

	CUDA_SAFE_CALL_ASYNC( cudaMemcpy3DAsync( &cpyPrms, stream ) );
      }



      //! Copies data into a CUDA array for texturing
      void SendArray( const cudaStream_t myStream = 0 ) {
	/*!
	  Method to copy the data currently on the GPU
	  into a CUDA array, which can then be used in a
	  texture.
	  Copy is done asynchronously within the (optionally)
	  given stream
	*/
	if( this->d_data.ptr == NULL ) {
	  std::cerr << __FUNCTION__
		    << ": GPU data not available!"
		    << std::endl;
	}
	
	if( this->dca_data == NULL ) {
	  std::cerr << __FUNCTION__
		    << ": CUDA array not allocated!"
		    << std::endl;
	  exit( EXIT_FAILURE );
	}
	
	cudaMemcpy3DParms cp = {0};
	
	cp.srcPtr = this->d_data;
	cp.dstArray = this->dca_data;
	cp.extent = this->extent;
	cp.kind = cudaMemcpyDeviceToDevice;
	
	CUDA_SAFE_CALL_ASYNC( cudaMemcpy3DAsync( &cp, myStream ) );
      }
      
      // --------------------------------------
      //! Method to return number of 'covering' blocks
      dim3 CoverBlocks( const unsigned int threadCount ) const {
	/*!
	  Method to return the number of blocks which
	  will be required to cover the volume of cubic blocks
	  of size threadCount
	*/
	dim3 grid;

	grid.x = DivRoundUp( this->dims.x, threadCount );
	grid.y = DivRoundUp( this->dims.y, threadCount );
	grid.z = DivRoundUp( this->dims.z, threadCount );

	return( grid );
      }

      // --------------------------------------

      //! Provides 1D indexing in to contiguous host memory
      unsigned int Index1D( const unsigned int ix,
			    const unsigned int iy,
			    const unsigned int iz ) const {
	/*!
	  Computes the offset into an unpadded contiguous
	  array on the host of the given element location.
	  This is for use by send and receive methods in
	  derived classes
	*/
	return( ix + ( this->dims.x * ( iy + ( this->dims.y * iz ) ) ) );
      }

      // ============================================
    protected:
   
      //! Dimensions of the volume
      dim3 dims;
      
      //! Extent of the allocated 3D array
      cudaExtent extent;
      //! Pointer to the allocated device memory
      cudaPitchedPtr d_data;
      //! CUDA array pointer for texturing
      cudaArray *dca_data;


      // ============================================
   
    private:
      // -------------------------------------------
      // Prevent copying

      //! Hidden copy constructor
      VolumeGPU( const VolumeGPU& src ) : dims(make_uint3(0,0,0)),
					  extent(make_cudaExtent(0,0,0)),
					  d_data(make_cudaPitchedPtr(NULL,0,0,0)),
					  dca_data(NULL) {
	std::cerr << __FUNCTION__
		  << ": Please don't use copy constructor"
		  << std::endl;
	exit( EXIT_FAILURE );
      }

      //! Hidden assignment operator
      VolumeGPU& operator=( const VolumeGPU& src ) {
	std::cerr << __FUNCTION__
		  << ": Please don't use assignment operator"
		  << std::endl;
	exit( EXIT_FAILURE );
      }

    };


  }

}

#endif
