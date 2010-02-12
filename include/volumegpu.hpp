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
 *    $Date: 2010/02/12 20:33:19 $
 *    $Revision: 1.1 $
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
      // Memory management

      //! Allocates storage of the given dimensions
      void Allocate( const dim3& myDims ) {
	
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

	this->ReleaseArray();

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
