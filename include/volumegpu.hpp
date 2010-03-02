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
 *    $Date: 2010/03/02 20:09:33 $
 *    $Revision: 1.18 $
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

//! Namespace to hold everything related to the GPU
namespace GPU {

  //! Namespace to hold the datatypes used by the GPU
  namespace Classes {

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
    class VolumeArgGPU {
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
			       const unsigned int iz ) const {
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
				 const unsigned int iz ) {
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
					const unsigned int iz ) const {
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
				const unsigned int iz ) const {
	bool res = ( ix < (this->dims.x) );
	res = res && ( iy < (this->dims.y) );
	res = res && ( iz < (this->dims.z) );

	return( res );
      }

      //! Checks is given location is in volume within given tolerance
      __device__ bool InFuzzyVolume( const float3& r,
				     const float tol ) const {
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



    //! Templated class to hold volume data on the GPU
    /*!
      This templated class provides a container for volume data.
      It provides a simple interface to padded arrays on the
      device, as well as the ability to copy data into a CUDA
      array (for texturing).
      Although it should not (and cannot) be used directly as
      a kernel argument, the VolumeArgGPU is provided for
      this purpose.
      @see VolumeArgGPU
    */
    template<typename T>
    class VolumeGPU {
    public:

      // -------------------------------------------
      // Constructors and destructor
      
      //! Default constructor
      VolumeGPU( void ) : dims(make_uint3(0,0,0)),
			  d_data(make_cudaPitchedPtr(NULL,0,0,0)),
			  dca_data(NULL) {};

      //! Destructor
      ~VolumeGPU( void ) {
	this->Release();
      }

      // -------------------------------------------
      // Conversion operator

      //! Converts a GPU volume to a kernel argument
      operator VolumeArgGPU<T>( void ) const {
	VolumeArgGPU<T> vag( this->dims, this->d_data.ptr, this->d_data.pitch );

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
	return "$Id: volumegpu.hpp,v 1.18 2010/03/02 20:09:33 rge21 Exp $";
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
	  @param[in] myDims The dimensions required
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
	cudaExtent tmpExtent = ExtentFromDims( this->dims );
	tmpExtent.width *= sizeof(T);

	// Allocate the memory
	CUDA_SAFE_CALL( cudaMalloc3D( &(this->d_data), tmpExtent ) );
      }

      //! Allocates storage matching dimensions of potentially different type
      template<typename U>
      void Allocate( const VolumeGPU<U>& src ) {
	/*!
	  This routine allocates storage of identical dimensions to the
	  given argument, but of potentially different type.
	  It is useful if intermediate results have to be in 
	  higher precision
	*/
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

	cudaExtent tmp = ExtentFromDims( this->dims );

	CUDA_SAFE_CALL( cudaMalloc3DArray( &(this->dca_data),
					   &cd,
					   tmp ) );
      }


      //! Releases all data on GPU
      void Release( void ) {
	this->ReleaseArray();
	if( this->d_data.ptr != NULL ) {
	  CUDA_SAFE_CALL( cudaFree( d_data.ptr ) );
	}
	this->dims = make_uint3(0,0,0);
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
	/*!
	  Prior to copying data to the GPU, the host must
	  organise the volume into contiguous pinned memory.
	  This routine supplies the number of bytes required
	  by the current class.
	*/
	unsigned int nElements;

	nElements = this->dims.x * this->dims.y * this->dims.z;

	return( nElements * sizeof(T) );
      }


      //! Allocates a pinned memory buffer on the host
      T* AllocateHostBuffer( void ) const {
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



      //! Receive a buffer from the GPU
      void RecvBuffer( T* const h_buffer,
		       const cudaStream_t stream = 0 ) const {
	/*!
	  Retrieves data from the GPU into the
	  contiguous buffer on the host.
	  The copy is done asynchronously, so the
	  user is responsible for stream management
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



      //! Copies data into a CUDA array for texturing
      void SendArray( const cudaStream_t myStream = 0 ) {
	/*!
	  Method to copy the data currently on the GPU
	  into a CUDA array, which can then be used in a
	  texture.
	  Copy is done asynchronously within the (optionally)
	  given stream
	  @param[in] myStream The CUDA stream in which the copy should occur
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
	cp.extent = ExtentFromDims( this->dims );
	cp.kind = cudaMemcpyDeviceToDevice;
	
	CUDA_SAFE_CALL_ASYNC( cudaMemcpy3DAsync( &cp, myStream ) );
      }
      
      // --------------------------------------
      //! Method to return number of 'covering' blocks
      dim3 CoverBlocks( const unsigned int threadCount ) const {
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
