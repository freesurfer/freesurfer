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
 *    $Author: rge21 $
 *    $Date: 2010/11/11 15:00:47 $
 *    $Revision: 1.4 $
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

#include <cuda_runtime.h>


namespace GPU {
  namespace Classes {

    // ================================================================
    // Memory Management


    template<typename T>
    void VolumeGPU<T>::Allocate( const dim3 myDims ) {
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


    // ---------------------------------

    template<typename T>
    void VolumeGPU<T>::AllocateArray( void ) {
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
    
    // ---------------------------------

    template<typename T>
    void VolumeGPU<T>::Release( void ) {
      this->ReleaseArray();
      if( this->d_data.ptr != NULL ) {
	CUDA_SAFE_CALL( cudaFree( d_data.ptr ) );
      }
      this->dims = make_uint3(0,0,0);
      this->d_data = make_cudaPitchedPtr(NULL,0,0,0);
    }


    // ---------------------------------

    template<typename T>
    void VolumeGPU<T>::ReleaseArray( void ) {
      if( this->dca_data != NULL ) {
	CUDA_SAFE_CALL( cudaFreeArray( this->dca_data ) );
	this->dca_data = NULL;
      }
    }

    // ---------------------------------

    template<typename T>
    size_t VolumeGPU<T>::BufferSize( void ) const {
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


    // ---------------------------------
    
    template<typename T>
    T* VolumeGPU<T>::AllocateHostBuffer( void ) const {
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

    template<typename T>
    void VolumeGPU<T>::Zero( void ) {
      
      // Have to set up an extent too
      cudaExtent extent = ExtentFromDims( this->dims );
      extent.width *= sizeof(T);

      CUDA_SAFE_CALL( cudaMemset3D( this->d_data,
				    0,
				    extent ) );
    }

    // ---------------------------------

    template<typename T>
    void VolumeGPU<T>::Copy( const VolumeGPU<T>& src,
			     const cudaStream_t myStream ) {
      this->Allocate( src.dims );

      cudaMemcpy3DParms copyParams = {0};

      copyParams.srcPtr = src->d_data;
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
				   const cudaStream_t stream ) {
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
				   const cudaStream_t stream ) const {
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
    

    // ---------------------------------
    
    template<typename T>
    void VolumeGPU<T>::SendArray( const cudaStream_t myStream ) {
      /*!
	Method to copy the data currently on the GPU
	into a CUDA array, which can then be used in a
	texture.
	Copy is done asynchronously within the (optionally)
	given stream.
	If no stream is given, it defaults to stream 0.
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


    // ================================================================
    // Other routines

    template<typename T>
    dim3 VolumeGPU<T>::CoverBlocks( const unsigned int threadCount ) const {
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
					const unsigned int iz ) const {
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
    dim3 VolumeGPU<T>::Index3D( unsigned int idx ) const {
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
