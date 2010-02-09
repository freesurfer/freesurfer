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
 *    $Date: 2010/02/09 18:29:38 $
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

      //! Accessor for the actual pointer
      __device__ const float* GetPointer( void ) const {
	return( &(this->matrix[0]) );
      }

    private:
      //! The matrix itself, stored row major
      float matrix[kMatrixSize];
    };



    //! Class to hold an affine transformation in shared memory
    class AffineTransShared {
    public:
      //! Elements in an affine vector
      static const unsigned int kVectorSize = AffineTransformation::kVectorSize;
      //! Size of the affine matrix
      static const unsigned int kMatrixSize = kVectorSize*kVectorSize;

      //! Constructor from pointer
      __device__ AffineTransShared( float *ptr ) : matrix(ptr) {};

      //! RHS subscripting operator
      __device__ float operator() ( const unsigned int i,
				    const unsigned int j ) const {
	return( this->matrix[j+(i*kVectorSize)] );
      }

      //! LHS subscripting operator
      __device__ float& operator() ( const unsigned int i,
				     const unsigned int j ) {
	return( this->matrix[j+(i*kVectorSize)] );
      }


      //! Routine to compute i and j for this thread
      __device__ void Getij( unsigned int& i, unsigned int& j ) const {
	i = threadIdx.x / kVectorSize;
	j = threadIdx.x % kVectorSize;
      }

      //! Routine to set up identity
      __device__ void SetIdentity( void ) {
	// Only done by first 16 threads
	if( threadIdx.x < kMatrixSize ) {
	 
	  unsigned int i, j;
	  this->Getij( i, j );

	  if( i == j ) {
	    matrix[threadIdx.x] = 1;
	  } else {
	    matrix[threadIdx.x] = 0;
	  }

	}

      }

      //! Routine to set up translation
      __device__ void SetTranslation( const float3& trans ) {
	/*!
	  This routine sets the affine matrix locations
	  corresponding to translation to the given
	  translation vector.
	  It leaves everything else untouched, so if you
	  want a pure transformation matrix, call
	  SetIdentity first.
	*/
	// Only done by first thread
	if( threadIdx.x == 0 ) {
	  (*this)(0,3) = trans.x;
	  (*this)(1,3) = trans.y;
	  (*this)(2,3) = trans.z;
	  (*this)(3,3) = 1;
	}
      }

      //! Routine to set up scaling
      __device__ void SetScaling( const float3& scale ) {
	/*!
	  This routine sets the affine matrix locations
	  corresponding to a scaling to the given
	  vector.
	  It leaves everything else untouched, so if you
	  want a pure scaling matrix, call
	  SetIdentity first.
	*/
	if( threadIdx.x == 0 ) {
	  (*this)(0,0) = scale.x;
	  (*this)(1,1) = scale.y;
	  (*this)(2,2) = scale.z;
	}
      }

      //! Routine to set up x rotation
      __device__void SetXRotation( const float theta ) {
	/*!
	  This routine sets the affine matrix locations
	  corresponding to a rotation about the X axis
	  of the given angle.
	  It leaves everything else untouched, so if you
	  want a pure rotation matrix, call
	  SetIdentity first.
	*/
	if( threadIdx.x == 0 ) {
	  float s, c;
	  sincosf( theta, &s, &c );
	  
	  (*this)(1,1) =  c;
	  (*this)(1,2) =  s;
	  (*this)(2,1) = -s;
	  (*this)(2,2) =  c;
	}
      }


      //! Routine to set up y rotation
      __device__void SetYRotation( const float theta ) {
	/*!
	  This routine sets the affine matrix locations
	  corresponding to a rotation about the Y axis
	  of the given angle.
	  It leaves everything else untouched, so if you
	  want a pure rotation matrix, call
	  SetIdentity first.
	*/
	if( threadIdx.x == 0 ) {
	  float s, c;
	  sincosf( theta, &s, &c );
	  
	  (*this)(0,0) =  c;
	  (*this)(0,2) = -s;
	  (*this)(2,0) =  s;
	  (*this)(2,2) =  c;
	}
      }


      //! Routine to set up z rotation
      __device__void SetZRotation( const float theta ) {
	/*!
	  This routine sets the affine matrix locations
	  corresponding to a rotation about the Z axis
	  of the given angle.
	  It leaves everything else untouched, so if you
	  want a pure rotation matrix, call
	  SetIdentity first.
	*/
	if( threadIdx.x == 0 ) {
	  float s, c;
	  sincosf( theta, &s, &c );
	  
	  (*this)(0,0) =  c;
	  (*this)(0,1) =  s;
	  (*this)(1,0) = -s;
	  (*this)(1,1) =  c;
	}
      }


      //! Routine to multiply two matrices
      __device__ void Multiply( const AffineTransShared& a,
			       const AffineTransShared& b ) {
	/*!
	  Performs the matrix-multiply a*b and stores the
	  result in the current object.
	  Multiplication will be performed by the first
	  warp of the block, and will have bank conflicts
	*/

	if( threadIdx.x < kMatrixSize ) {
	  unsigned int i, j;
	  this->Getij( i, j );

	  float myVal = 0;

	  #pragma unroll 4
	  for( unsigned int k=0; k<kVectorSize; k++ ) {
	    myVal += a(i,k) * b(k,j);
	  }

	  (*this)(i,j) = myVal;
	  
	}
      }


      //! Routine to transform a float3 vector
      __device__ float3 transform( const float3& a ) const {
	float r1[kVectorSize], r2[kVectorSize];
	
	// Set up in put affine vector
	r1[0] = a.x;
	r1[1] = a.y;
	r1[2] = a.z;
	r1[3] = 1;

	// Do the multiplication
	for( unsigned int i=0; i<kVectorSize; i++ ) {
	  r2[i] = 0;
	  #pragma unroll 4
	  for( unsigned int j=0; j<kVectorSize; j++ ) {
	    r2[i] += (*this)( i, j ) * r1[j];
	  }
	}

	// Construct float3 return value
	return( make_float3( r2[0], r2[1], r2[2] ) );
      }
	

    private:
      //! Default constructor (don't use!)
      __device__ AffineTransShared( void ) : matrix(NULL) {};

      //! Pointer to relevant memory
      float *matrix;
    };
      
  }
}


#endif
