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
 *    $Date: 2010/02/16 16:22:02 $
 *    $Revision: 1.30 $
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

#include "mri.h"

#include "volumegpu.hpp"

#include "cudatypeutils.hpp"
#include "cudacheck.h"

// ==================================================================

namespace GPU {

  namespace Classes {

    // Forward declaration of the class which will be passed to kernels
    template<typename T> class MRIframeOnGPU;

    // ================================================================


    //! Templated class to hold an MRI frame on the GPU
    template<typename T>
    class MRIframeGPU : public VolumeGPU<T> {
    public:
      // -------------------------------------------------------

      // Declare the related 'kernel' class a friend for private access
      friend class MRIframeOnGPU<T>;
  
      // --------------------------------------------------------
      // Constructors & destructors
      MRIframeGPU<T>( void ) : VolumeGPU<T>() {};

      // --------------------------------------------------------
      // Data accessors

      //! Return information about the file version
      const char* VersionString( void ) const {
	return "$Id: mriframegpu.hpp,v 1.30 2010/02/16 16:22:02 rge21 Exp $";
      }
      
      // --------------------------------------------------------
      // Memory manipulation
      
      // -----

      //! Extracts frame dimensions from a given MRI and allocates the memory
      void Allocate( const MRI* src ) {
	/*!
	  Fills out the cpuDims, gpuDims and extent data members.
	  Uses this information to call cudaMalloc3D
	  @params[in] src The MRI to use as a template
	*/
	
	// Sanity checks
	if( src->type != this->MRItype()  ) {
	  std::cerr << __PRETTY_FUNCTION__
		    << ": MRI type mismatch against "
		    << src->type << std::endl;
	  exit( EXIT_FAILURE );
	}

	
	dim3 myDims = make_uint3( src->width, src->height, src->depth );

	static_cast<VolumeGPU<T>*>(this)->Allocate( myDims );
      }

      // -----

      //! Allocates matching storage of possibly different type
      template<typename U>
      void Allocate( const MRIframeGPU<U>& src ) {
	static_cast<VolumeGPU<T>*>(this)->Allocate( static_cast<const VolumeGPU<U>& >(src) );
      }

      // -----


      // -----
     

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
	
	const size_t bSize = this->BufferSize();
	// See if we were supplied with workspace
	if( h_work != NULL ) {
	  h_data = reinterpret_cast<T*>(h_work);
	} else {
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
	
	const size_t bSize = this->BufferSize();
	
	// Allocate contiguous host memory if needed
	if( h_work != NULL ) {
	  h_data = reinterpret_cast<T*>(h_work);
	} else {
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
      

      //! Method to check if the CUDA array is allocated
      bool HasArray( void ) const {
	return( this->dca_data != NULL );
      }

      // --------------------------------------------------------------------
      //! Method to return number of 'covering' blocks
      dim3 CoverBlocks( const dim3 threadBlock ) const {
	/*!
	  Method to return the number of blocks required to
	  'cover' the frame with threads, given the dimensions
	  of the thread block.
	*/
	dim3 grid;

	grid.x = this->DivRoundUp( this->dims.x, threadBlock.x );
	grid.y = this->DivRoundUp( this->dims.y, threadBlock.y );
	grid.z = this->DivRoundUp( this->dims.z, threadBlock.z );

	return( grid );
      }

      //! Returns a / b, rounded up to next integer
      unsigned int DivRoundUp( const unsigned int a,
			       const unsigned int b ) const {
	float tmp;

	tmp = static_cast<float>(a) / b;

	return( static_cast<unsigned int>( ceilf( tmp ) ) );
      }
 
    private:
      // --------------------------------------------------------------------
      

      // ----------------------------------------------------------------------
      // Prevent copying
  
      
      // ----------------------------------------------------------------------

      // ----------------------------------------------------------------------
      
      //! Function to sanity check dimensions
      bool CheckDims( const MRI* mri ) const {

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
	for( unsigned int iz=0; iz<this->dims.z; iz++ ) {
	  for( unsigned int iy=0; iy<this->dims.y; iy++ ) {
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
	for( unsigned int iz=0; iz<this->dims.z; iz++ ) {
	  for( unsigned int iy=0; iy<this->dims.y; iy++ ) {
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
	return( ix + ( this->dims.x * ( iy + ( this->dims.y * iz ) ) ) );
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
    
    // ======================================================================

    
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
      //! Out of range value
      static const T kOutOfRangeVal = 0;

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
      MRIframeOnGPU( const MRIframeGPU<T>& src ) : dims(src.dims),
						   extent(src.extent),
						   data(src.d_data) {};
      
      
      // --------------------------------------------------------
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
