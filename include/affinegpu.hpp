/**
 * @file  affinegpu.hpp
 * @brief Holds affine transformation class for the GPU
 *
 * Holds an affine transformation type for the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/02/05 17:30:19 $
 *    $Revision: 1.2 $
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

#ifndef AFFINE_TRANSFORMATION_GPU_H
#define AFFINE_TRANSFORMATION_GPU_H

extern "C" {
#include "matrix.h"
}

// ==================================================================

namespace GPU {

  namespace Classes {
    
    //! Class to hold an affine transformation
    class AffineTransformation {
    public:
      //! Elements in an affine vector
      static const unsigned int kVectorSize = 4;
      //! Size of the affine matrix
      static const unsigned int kMatrixSize = kVectorSize*kVectorSize;

      //! Default Constructor
      AffineTransformation( void );

      //! Constructor from MATRIX
      AffineTransformation( const MATRIX* src );

      //! RHS subscripting operator (no bounds check)
      __host__ __device__ float operator() ( const unsigned int i,
					     const unsigned int j ) const {
	return( this->matrix[j+(i*kVectorSize)] );
      }

      //! LHS subscripting operator (no bounds check)
      __host__ __device__ float& operator() ( const unsigned int i,
					      const unsigned int j ) {
	return( this->matrix[j+(i*kVectorSize)] );
      }

      //! Function to transform a location
      __device__ float3 transform( const float3& rIn ) const {
	/*!
	  Performs the transformation described by
	  the class to the given input location
	*/
	float r1[kVectorSize], r2[kVectorSize];

	// Set up in put affine vector
	r1[0] = rIn.x;
	r1[1] = rIn.y;
	r1[2] = rIn.z;
	r1[3] = 1;

	// Do the multiplication
	for( unsigned int i=0; i<kVectorSize; i++ ) {
	  r2[i] = 0;
	  #pragma unroll 4
	  for( unsigned int j=0; j<kVectorSize; j++ ) {
	    r2[i] += this->operator() ( i, j ) * r1[j];
	  }
	}

	// Construct float3 return value
	return( make_float3( r2[0], r2[1], r2[2] ) );
      }

    private:
      //! The matrix itself, stored row major
      float matrix[kMatrixSize];
    };
  }
}


#endif
