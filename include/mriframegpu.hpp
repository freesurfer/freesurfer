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
 *    $Date: 2010/02/04 18:25:57 $
 *    $Revision: 1.22 $
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

namespace GPU {

  namespace Classes {

    // Forward declaration of the class which will be passed to kernels
    template<typename T> class MRIframeOnGPU;

    // ================================================================


    //! Templated class to hold an MRI frame on the GPU
    template<typename T>
    class MRIframeGPU {
    public:
      // -------------------------------------------------------

      // Declare the related 'kernel' class a friend for private access
      friend class MRIframeOnGPU<T>;
  
      // --------------------------------------------------------
      // Constructors & destructors

      //! Default constructor
      MRIframeGPU( void ) : cpuDims(make_uint3(0,0,0)),
			    gpuDims(make_uint3(0,0,0)),
			    extent(make_cudaExtent(0,0,0)),
			    d_data(make_cudaPitchedPtr(NULL,0,0,0)),
			    dca_data(NULL) {};

 
      //! Destructor
      ~MRIframeGPU( void ) {
	this->ReleaseArray();
	this->Release();
      }

      // --------------------------------------------------------
      // Data accessors

      //! Return the dimensions on the GPU of this MRI frame
      dim3 GetGPUDims( void ) const {
	return( this->gpuDims );
      }

      //! Return the dimensions on the CPU of this MRI frame
      dim3 GetCPUDims( void ) const {
	return( this->cpuDims );
      }

      //! Return information about the file version
      const char* VersionString( void ) const {
	return "$Id: mriframegpu.hpp,v 1.22 2010/02/04 18:25:57 rge21 Exp $";
      }
      
      //! Return pointer to the cudaArray
      cudaArray* GetArray( void ) const {
	return( this->dca_data );
      }

  
      // --------------------------------------------------------
      // Memory manipulation
      
      //! Supplies the size of buffer required on the host
      size_t GetBufferSize( void ) const {
	unsigned int nElements;
	
	nElements = this->gpuDims.x * this->gpuDims.y * this->gpuDims.z;

	return( nElements * sizeof(T) );
      }

      // -----

      //! Extracts frame dimensions from a given MRI and allocates the memory
      void Allocate( const MRI* src,
		     const unsigned int padSize = 1 ) {
	/*!
	  Fills out the cpuDims, gpuDims and extent data members.
	  Uses this information to call cudaMalloc3D
	  @params[in] src The MRI to use as a template
	  @params[in] padSize Each dimension will be padded to be a multiple of this
	*/
	
	// Sanity checks
	if( src->type != this->MRItype()  ) {
	  std::cerr << __PRETTY_FUNCTION__
		    << ": MRI type mismatch against "
		    << src->type << std::endl;
	  exit( EXIT_FAILURE );
	}

	if( padSize == 0 ) {
	  std::cerr << __FUNCTION__
		    << ": Must have non-zero padSize"
		    << std::endl;
	  exit( EXIT_FAILURE );
	}
	
	
	dim3 myDims = make_uint3( src->width, src->height, src->depth );

	this->Allocate( myDims, padSize );
      }

      // -----

      //! Allocates storage based on input dimensions and padding
      void Allocate( const dim3 cpuDims,
		     const unsigned int padSize = 1 ) {
	/*!
	  Given a set of dimensions on the CPU, sets up cpuDims, gpuDims
	  and extent data members based on the padSize argument.
	  Uses this to do a cudaMalloc3D
	*/

	if( padSize == 0 ) {
	  std::cerr << __FUNCTION__
		    << ": Must have non-zero padSize"
		    << std::endl;
	  exit( EXIT_FAILURE );
	}
	
	
	// Get rid of old data
	this->ReleaseArray();
	this->Release();
	
	// Make a note of the dimensions on the CPU
	this->cpuDims = cpuDims;
	
	// Generate the dimensions on the GPU
	this->gpuDims.x = this->RoundToBlock( this->cpuDims.x, padSize );
	this->gpuDims.y = this->RoundToBlock( this->cpuDims.y, padSize );
	this->gpuDims.z = this->RoundToBlock( this->cpuDims.z, padSize );
	
	// Do the actual allocation
	this->AllocateFromDims();
	
	
      }

      // -----

      //! Allocates storage to match dimensions of given MRIframe
      template<typename U>
      void Allocate( const MRIframeGPU<U>& src ) {
	/*!
	  Copies the dimensions from the given source
	  MRIframeGPU (which may be of different datatype)
	  and then allocates memory for the current object
	*/
	
	// Get rid of the old
	this->ReleaseArray();
	this->Release();

	// Get the new
	this->cpuDims = src.GetCPUDims();
	this->gpuDims = src.GetGPUDims();

	// Do the allocation
	this->AllocateFromDims();
      }

      // -----
      
      //! Releases memory associated with class instance
      void Release( void ) {
	if( this->d_data.ptr != NULL ) {
	  CUDA_SAFE_CALL( cudaFree( d_data.ptr ) );
	}
	this->cpuDims = make_uint3(0,0,0);
	this->gpuDims = make_uint3(0,0,0);
	this->extent = make_cudaExtent(0,0,0);
	this->d_data = make_cudaPitchedPtr(NULL,0,0,0);
      }

      // -----

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

	// Shouldn't need to release old memory - happens in Allocate

	cudaChannelFormatDesc cd = cudaCreateChannelDesc<T>();

	CUDA_SAFE_CALL( cudaMalloc3DArray( &(this->dca_data),
					   &cd,
					   this->extent ) );
      }
      
      // -----

      //! Releases the CUDA array member
      void ReleaseArray( void ) {
	if( this->dca_data != NULL ) {
	  CUDA_SAFE_CALL( cudaFreeArray( this->dca_data ) );
	  this->dca_data = NULL;
	}
      }

      // --------------------------------------------------------
      // Data transfer
      
      //! Send the given MRI frame to the GPU
      void Send( const MRI* src,
		 const unsigned int iFrame,
		 void* const h_work = NULL,
		 const cudaStream_t stream = 0 ) {
	/*!
	  Sends the given MRI frame to the GPU.
	  Optional arguments can be used to supply page-locked
	  host memory for the transfer, and a stream in which
	  to perform the transfer.
	  If supplied, the array h_work must be at least
	  this->GetBufferSize() bytes long.
	  Furthermore, the calling routine is responsible
	  for synchronisation
	*/
	
	T* h_data;
	
	// Start with some sanity checks
	this->VerifyMRI( src );
	
	if( iFrame >= static_cast<unsigned int>(src->nframes) ) {
	  std:: cerr << __FUNCTION__
		     << ": Bad frame requested " << iFrame
		     << std::endl;
	  exit( EXIT_FAILURE );
	}
	
	const size_t bSize = this->GetBufferSize();
	// See if we were supplied with workspace
	if( h_work != NULL ) {
	  h_data = reinterpret_cast<T*>(h_work);
	} else {
	  // Allocate contiguous host memory
	  CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_data,
					 bSize,
					 cudaHostAllocDefault ) );
	}
	
	// Zero the memory if needed
	if( this->IsPadded() ) {
	  memset( h_data, 0 , bSize );
	}
	
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
	
	CUDA_SAFE_CALL_ASYNC( cudaMemcpy3DAsync( &copyParams, stream ) );
	
	// Release host memory if needed
	if( h_work == NULL ) {
	  CUDA_SAFE_CALL( cudaStreamSynchronize( stream ) );
	  CUDA_SAFE_CALL( cudaFreeHost( h_data ) );
	}
      }
      
      // -----
      
      
      //! Receives the given MRI frame from the GPU
      void Recv( MRI* dst,
		 const unsigned int iFrame,
		 void* const h_work = NULL,
		 const cudaStream_t stream = 0 ) const {
	/*!
	  Retrieves the given MRI frame from the GPU.
	  Optional arguments can be used to supply page-locked
	  host memory for the transfer, and a stream in which
	  to perform the transfer.
	  If supplied, the array h_work must be at least
	  this->GetBufferSize() bytes long
	*/
	
	T* h_data;
	
	// Start with sanity checks
	this->VerifyMRI( dst );
	
	if( iFrame >= static_cast<unsigned int>(dst->nframes) ) {
	  std:: cerr << __FUNCTION__
		     << ": Bad frame requested " << iFrame
		     << std::endl;
	  exit( EXIT_FAILURE );
	}
	
	const size_t bSize = this->GetBufferSize();
	
	// Allocate contiguous host memory if needed
	if( h_work != NULL ) {
	  h_data = reinterpret_cast<T*>(h_work);
	} else {
	  CUDA_SAFE_CALL( cudaHostAlloc( (void**)&h_data,
					 bSize,
					 cudaHostAllocDefault ) );
	}
	
	// Retrieve from GPU
	cudaMemcpy3DParms cpyPrms = {0};
	cpyPrms.srcPtr = this->d_data;
	cpyPrms.dstPtr = make_cudaPitchedPtr( (void*)h_data,
					      gpuDims.x*sizeof(T),
					      gpuDims.x,
					      gpuDims.y );
	cpyPrms.extent = this->extent;
	cpyPrms.kind = cudaMemcpyDeviceToHost;

	CUDA_SAFE_CALL_ASYNC( cudaMemcpy3DAsync( &cpyPrms, stream ) );
	CUDA_SAFE_CALL( cudaStreamSynchronize( stream ) );
	
	// Retrieve from contiguous RAM
	this->InhumeFrame( dst, h_data, iFrame );
	
	// Release host memory
	if( h_work == NULL ) {
	  CUDA_SAFE_CALL( cudaFreeHost( h_data ) );
	}
      }


      // -----

      //! Copies frame data to a CUDA array for texturing
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
      
      // ----------------------------------------------------------------------
      //! Method to sanity check MRI
      void VerifyMRI( const MRI* mri ) const {
	
	if( mri->type != this->MRItype()  ) {
	  std::cerr << __PRETTY_FUNCTION__
		    << ": MRI type mismatch against "
		    << mri->type << std::endl;
	  exit( EXIT_FAILURE );
	}
	
	if( !this->CheckDims( mri ) ) {
	  std::cerr << __PRETTY_FUNCTION__
		    << ": Size mismatch"
		    << std::endl;
	  exit( EXIT_FAILURE );
	}
      }
      
      //! Method to return MRI type (has specialisations below)
      int MRItype( void ) const {
	return(-1);
      }
      
      //! Method to check padding size of an MRI
      bool CheckPadding( const unsigned int padSize ) const {
	/*!
	  Some routines may require that the data on the GPU
	  are sized to a multiple of padSize.
	  This routine checks the size of gpuDims.
	*/
	
	bool res;
	
	res = ( (this->gpuDims.x % padSize) == 0 );
	res = res && ( (this->gpuDims.y % padSize) == 0 );
	res = res && ( (this->gpuDims.z % padSize) == 0 );
	
	return( res );
      }

      //! Method to check if the CUDA array is allocated
      bool HasArray( void ) const {
	return( this->dca_data != NULL );
      }
  
 
    private:
      // --------------------------------------------------------------------
      // Data members

      //! Original data size
      dim3 cpuDims;
      //! Data size on GPU (may be padded)
      dim3 gpuDims;
      //! Extent of allocated 3D array
      cudaExtent extent;
      //! Pointer to the allocated memory
      cudaPitchedPtr d_data;
      //! CUDA array for texturing
      cudaArray *dca_data;

      // ----------------------------------------------------------------------
      // Prevent copying
      
      //! Copy constructor - don't use
      MRIframeGPU( const MRIframeGPU& src ) : cpuDims(make_uint3(0,0,0)),
					      gpuDims(make_uint3(0,0,0)),
					      extent(make_cudaExtent(0,0,0)),
					      d_data(make_cudaPitchedPtr(NULL,0,0,0)) {
	std::cerr << __PRETTY_FUNCTION__
		  << ": Please don't use copy constructor"
		  << std::endl;
	exit( EXIT_FAILURE );
      }
      
      //! Assignment operator - don't use
      MRIframeGPU& operator=( const MRIframeGPU &src ) {
	std::cerr << __PRETTY_FUNCTION__
		  << ": Please don't use assignment operator"
		  << std::endl;
	exit( EXIT_FAILURE );
      }

      
      // ----------------------------------------------------------------------

      //! Allocates extent and memory space from object's gpuDims member
      void AllocateFromDims( void ) {
	/*!
	  This routine allocates GPU memory to hold an MRI frame, and
	  also sets up the extent data member.
	  It assumes that the gpuDims data member has been correctly
	  set prior to the call.
	*/
	// Make the extent
	this->extent = make_cudaExtent( this->gpuDims.x * sizeof(T),
					this->gpuDims.y,
					this->gpuDims.z );
	
	// Allocate the memory
	CUDA_SAFE_CALL( cudaMalloc3D( &(this->d_data), this->extent ) );
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


      //! Checks to see if we're padding on the GPU
      bool IsPadded( void ) const {
	/*!
	  It's sometimes useful to know if the sizes
	  on the CPU and GPU match exactly.
	  This method is used for checking for
	  this property.
	*/
	bool res;
	
	// Check for equality, return inverse
	res = (this->cpuDims.x == this->gpuDims.x );
	res = res && (this->cpuDims.y == this->gpuDims.y );
	res = res && (this->cpuDims.z == this->gpuDims.z );
	
	return( !res );
      }
      
      // ----------------------------------------------------------------------
      
      //! Function to sanity check dimensions
      bool CheckDims( const MRI* mri ) const {
	
	bool goodDims;
	
	goodDims = ( static_cast<unsigned int>(mri->width) == this->cpuDims.x );
	goodDims = goodDims &&
	  ( static_cast<unsigned int>(mri->height) == this->cpuDims.y );
	goodDims = goodDims &&
	  ( static_cast<unsigned int>(mri->depth) == this->cpuDims.z );
	
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
	  so everything is not fanatically verified.
	  The name of this method should not be taken as a value judgement
	  of the CPU-side datastructure.
	  @param[in] src The source MRI
	  @param[out] h_slab Pointer to the destination array (must be already allocated)
	  @param[in] iFrame Which frame to grab
	*/
	
	// Start with a few sanity checks
	if( iFrame >= static_cast<unsigned int>(src->nframes) ) {
	  std::cerr << __FUNCTION__
		    << ": iFrame out of range"
		    << std::endl;
	  exit( EXIT_FAILURE );
	}
	
	if( h_slab == NULL ) {
	  std::cerr << __FUNCTION__
		    << ": h_slab unallocated"
		    << std::endl;
	  exit( EXIT_FAILURE );
	}
	
	this->VerifyMRI( src );
	
	
	// Main extraction loop
	for( unsigned int iz=0; iz<this->cpuDims.z; iz++ ) {
	  for( unsigned int iy=0; iy<this->cpuDims.y; iy++ ) {
	    unsigned int iStart = this->Index1D( 0, iy, iz );
	    
	    this->ExhumeRow( src, &( h_slab[iStart] ), iy, iz, iFrame );
	  }
	} 
      }
      
      // -----

      //! Copies contiguous memory to an MRI frame
      void InhumeFrame( MRI* dst,
			const T *h_slab,
			const unsigned int iFrame ) const {
	/*!
	  Copies a block of contiguous host memory on the host
	  into an MRI frame.
	  Assumes that everything is all properly allocated,
	  so things are not fanatically verified
	*/
	
	// Start with a few sanity checks
	if( iFrame >= static_cast<unsigned int>(dst->nframes) ) {
	  std::cerr << __FUNCTION__
		    << ": iFrame out of range"
		    << std::endl;
	  exit( EXIT_FAILURE );
	}
	
	if( h_slab == NULL ) {
	  std::cerr << __FUNCTION__
		    << ": h_slab unallocated"
		    << std::endl;
	  exit( EXIT_FAILURE );
	}
	
	this->VerifyMRI( dst );
	
	// Main extraction loop
	for( unsigned int iz=0; iz<this->cpuDims.z; iz++ ) {
	  for( unsigned int iy=0; iy<this->cpuDims.y; iy++ ) {
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
		      const unsigned int iFrame ) const {
	/*!
	  Specialised method to copy a row of an MRI frame
	  into contiguous memory.
	  The row of x voxels are specified by iy, iz and iFrame.
	  This default method aborts the program
	*/
	
	std::cerr << __PRETTY_FUNCTION__
		  << ": Unrecognised data type "
		  << src->type << std::endl;
	exit( EXIT_FAILURE );
      }

      // -----

      //! Wrapper around memcpy for contiguous->MRI transfers
      void InhumeRow( MRI* dst,
		      const T *h_slab,
		      const unsigned int iy,
		      const unsigned int iz,
		      const unsigned int iFrame ) const {
	/*!
	  Specialised method to copy a chunk of contiguous
	  memory into a row of an MRI frame.
	  The row of x voxels are specified by iy, iz and iFrame.
	  This default method aborts the program
	*/
	
	std::cerr << __PRETTY_FUNCTION__
		  << ": Unrecognised data type "
		  << dst->type << std::endl;
	exit( EXIT_FAILURE );
      }
      
      
      // ----------------------------------------------------------------------
      
      // Index into host memory defined by gpuDims
      unsigned int Index1D( const unsigned int ix,
			    const unsigned int iy,
			    const unsigned int iz ) const {
	return( ix + ( this->gpuDims.x * ( iy + ( this->gpuDims.y * iz ) ) ) );
      }

    };


    // Declarations of specialised methods
    
    template<> int MRIframeGPU<unsigned char>::MRItype( void ) const;
    template<> int MRIframeGPU<short>::MRItype( void ) const;
    template<> int MRIframeGPU<float>::MRItype( void ) const;
    


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
    
    // ===================================================================================

    
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
      

      //! Utility function to convert float to the class' datatype
      __device__ T ConvertFloat( const float in ) const {
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

    template<> __device__
    unsigned char MRIframeOnGPU<unsigned char>::ConvertFloat( const float in ) const {
      return( static_cast<unsigned char>( rintf( in ) ) );
    }

    template<> __device__
    short MRIframeOnGPU<short>::ConvertFloat( const float in ) const {
      return( static_cast<short>( rintf( in ) ) );
    }

    template<> __device__
    float MRIframeOnGPU<float>::ConvertFloat( const float in ) const {
      return( in );
    }
    
  }
}


#endif
