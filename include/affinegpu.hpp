/**
 * @file  affinegpu.hpp
 * @brief Holds affine transformation class for the GPU
 *
 * Holds an affine transformation type for the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.11 $
 *
 * Copyright (C) 2002-2010,
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

#include "affine.h"

// ==================================================================

namespace GPU
{

namespace Classes
{


//! Class to hold an affine transformation
class AffineTransformation
{
public:
  //! Elements in an affine vector
  static const unsigned int kVectorSize = 4;
  //! Size of the affine matrix
  static const unsigned int kMatrixSize = kVectorSize*kVectorSize;

  //! Default Constructor
  AffineTransformation( void );

  //! Constructor from MATRIX
  AffineTransformation( const MATRIX* src );


  //! Set transform from MATRIX
  void SetTransform( const MATRIX* src );

  //! Set transform from AffineMatrix
  void SetTransform( const AffineMatrix* src );

  //! RHS subscripting operator (no bounds check)
  __host__ __device__ float operator() ( const unsigned int i,
                                         const unsigned int j ) const
  {
    return( this->matrix[j+(i*kVectorSize)] );
  }

  //! LHS subscripting operator (no bounds check)
  __host__ __device__ float& operator() ( const unsigned int i,
                                          const unsigned int j )
  {
    return( this->matrix[j+(i*kVectorSize)] );
  }

  //! Function to transform a location
  __device__ float3 transform( const float3& rIn ) const
  {
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
#pragma unroll 4
    for( unsigned int i=0; i<kVectorSize; i++ )
    {
      r2[i] = 0;
#pragma unroll 4
      for( unsigned int j=0; j<kVectorSize; j++ )
      {
        r2[i] += this->operator() ( i, j ) * r1[j];
      }
    }

    // Construct float3 return value
    return( make_float3( r2[0], r2[1], r2[2] ) );
  }

  //! Accessor for the actual pointer
  __device__ const float* GetPointer( void ) const
  {
    return( &(this->matrix[0]) );
  }

private:
  //! The matrix itself, stored row major
  float matrix[kMatrixSize];
};



//! Class to hold an affine transformation in shared memory
class AffineTransShared
{
public:
  //! Elements in an affine vector
  static const unsigned int kVectorSize = AffineTransformation::kVectorSize;
  //! Size of the affine matrix
  static const unsigned int kMatrixSize = kVectorSize*kVectorSize;

  //! Constructor from pointer
  __device__ AffineTransShared( float *ptr ) : matrix(ptr) {};

  //! RHS subscripting operator
  __device__ float operator() ( const unsigned int i,
                                const unsigned int j ) const
  {
    return( this->matrix[j+(i*kVectorSize)] );
  }

  //! LHS subscripting operator
  __device__ float& operator() ( const unsigned int i,
                                 const unsigned int j )
  {
    return( this->matrix[j+(i*kVectorSize)] );
  }


  //! Routine to compute i and j for this thread
  __device__ void Getij( unsigned int& i, unsigned int& j ) const
  {
    i = threadIdx.x / kVectorSize;
    j = threadIdx.x % kVectorSize;
  }

  //! Routine to set up identity
  __device__ void SetIdentity( void )
  {
    // Only done by first 16 threads
    if( threadIdx.x < kMatrixSize )
    {

      unsigned int i, j;
      this->Getij( i, j );

      if( i == j )
      {
        matrix[threadIdx.x] = 1;
      }
      else
      {
        matrix[threadIdx.x] = 0;
      }

    }

  }

  //! Routine to set up translation
  __device__ void SetTranslation( const float3& trans )
  {
    /*!
      This routine sets the affine matrix locations
      corresponding to translation to the given
      translation vector.
      It leaves everything else untouched, so if you
      want a pure transformation matrix, call
      SetIdentity first.
    */
    // Only done by first thread
    if( threadIdx.x == 0 )
    {
      (*this)(0,3) = trans.x;
      (*this)(1,3) = trans.y;
      (*this)(2,3) = trans.z;
      (*this)(3,3) = 1;
    }
  }

  //! Routine to set up scaling
  __device__ void SetScaling( const float3& scale )
  {
    /*!
      This routine sets the affine matrix locations
      corresponding to a scaling to the given
      vector.
      It leaves everything else untouched, so if you
      want a pure scaling matrix, call
      SetIdentity first.
    */
    if( threadIdx.x == 0 )
    {
      (*this)(0,0) = scale.x;
      (*this)(1,1) = scale.y;
      (*this)(2,2) = scale.z;
    }
  }

  //! Routine to set up x rotation
  __device__ void SetXRotation( const float theta )
  {
    /*!
      This routine sets the affine matrix locations
      corresponding to a rotation about the X axis
      of the given angle.
      It leaves everything else untouched, so if you
      want a pure rotation matrix, call
      SetIdentity first.
    */
    if( threadIdx.x == 0 )
    {
      float s, c;
      sincosf( theta, &s, &c );

      (*this)(1,1) =  c;
      (*this)(1,2) =  s;
      (*this)(2,1) = -s;
      (*this)(2,2) =  c;
    }
  }


  //! Routine to set up y rotation
  __device__ void SetYRotation( const float theta )
  {
    /*!
      This routine sets the affine matrix locations
      corresponding to a rotation about the Y axis
      of the given angle.
      It leaves everything else untouched, so if you
      want a pure rotation matrix, call
      SetIdentity first.
    */
    if( threadIdx.x == 0 )
    {
      float s, c;
      sincosf( theta, &s, &c );

      (*this)(0,0) =  c;
      (*this)(0,2) = -s;
      (*this)(2,0) =  s;
      (*this)(2,2) =  c;
    }
  }


  //! Routine to set up z rotation
  __device__ void SetZRotation( const float theta )
  {
    /*!
      This routine sets the affine matrix locations
      corresponding to a rotation about the Z axis
      of the given angle.
      It leaves everything else untouched, so if you
      want a pure rotation matrix, call
      SetIdentity first.
    */
    if( threadIdx.x == 0 )
    {
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
                            const AffineTransShared& b )
  {
    /*!
      Performs the matrix-multiply a*b and stores the
      result in the current object.
      Multiplication will be performed by the first
      warp of the block, and will have bank conflicts
    */

    if( threadIdx.x < kMatrixSize )
    {
      unsigned int i, j;
      this->Getij( i, j );

      float myVal = 0;

#pragma unroll 4
      for( unsigned int k=0; k<kVectorSize; k++ )
      {
        myVal += a(i,k) * b(k,j);
      }

      (*this)(i,j) = myVal;

    }
  }


  //! Routine to transform a float3 vector
  __device__ float3 transform( const float3& a ) const
  {
    float r1[kVectorSize], r2[kVectorSize];

    // Set up in put affine vector
    r1[0] = a.x;
    r1[1] = a.y;
    r1[2] = a.z;
    r1[3] = 1;

    // Do the multiplication
#pragma unroll 4
    for( unsigned int i=0; i<kVectorSize; i++ )
    {
      r2[i] = 0;
#pragma unroll 4
      for( unsigned int j=0; j<kVectorSize; j++ )
      {
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


//! Class to hold a 2x2 matrix
class Matrix2x2
{
public:
  static const unsigned int kVectorSize = 2;
  static const unsigned int kMatrixSize = kVectorSize*kVectorSize;


  //! RHS subscripting operator
  __device__ float operator() ( const unsigned int i,
                                const unsigned int j ) const
  {
    return( this->matrix[j+(i*kVectorSize)] );
  }

  //! LHS subscripting operator
  __device__ float& operator() ( const unsigned int i,
                                 const unsigned int j )
  {
    return( this->matrix[j+(i*kVectorSize)] );
  }


  //! Assignment operator
  __device__ Matrix2x2& operator=( const Matrix2x2& src )
  {
#pragma unroll 2
    for( unsigned int i=0; i<kVectorSize; i++ )
    {
#pragma unroll 2
      for( unsigned int j=0; j<kVectorSize; j++ )
      {
        (*this)(i,j) = src(i,j);
      }
    }

    return *this;
  }


  //! Self addition of another 2x2 matrix
  __device__ Matrix2x2& operator+=( const Matrix2x2& b )
  {
#pragma unroll 2
    for( unsigned int i=0; i<kVectorSize; i++ )
    {
#pragma unroll 2
      for( unsigned int j=0; j<kVectorSize; j++ )
      {
        (*this)(i,j) += b(i,j);
      }
    }

    return *this;
  }

  //! Addition of 2x2 matrices
  __device__ Matrix2x2 operator+( const Matrix2x2& b ) const
  {
    return (Matrix2x2( *this ) += b );
  }


  //! Self subtraction of another 2x2 matrix
  __device__ Matrix2x2& operator-=( const Matrix2x2& b )
  {
#pragma unroll 2
    for( unsigned int i=0; i<kVectorSize; i++ )
    {
#pragma unroll 2
      for( unsigned int j=0; j<kVectorSize; j++ )
      {
        (*this)(i,j) -= b(i,j);
      }
    }

    return *this;
  }

  //! Subtraction of 2x2 matrices
  __device__ Matrix2x2 operator-( const Matrix2x2& b ) const
  {
    return (Matrix2x2( *this ) -= b );
  }

  //! Multiplication of 2x2 matrices
  __device__ Matrix2x2 operator*( const Matrix2x2& b ) const
  {
    Matrix2x2 res;

#pragma unroll 2
    for( unsigned int i=0; i<kVectorSize; i++ )
    {
#pragma unroll 2
      for( unsigned int j=0; j<kVectorSize; j++ )
      {
        res(i,j) = 0;
#pragma unroll 2
        for( unsigned int k=0; k<kVectorSize; k++ )
        {
          res(i,j) += (*this)(i,k) * b(k,j);
        }
      }
    }

    return( res );
  }

  //! Determinant calculation
  __device__ float det( void ) const
  {
    float res;
    res = (*this)(0,0) * (*this)(1,1);
    res -= (*this)(0,1) * (*this)(1,0);
    return( res );
  }

  //! Inversion in place
  __device__ Matrix2x2 inverse( void ) const
  {
    Matrix2x2 res;

    // Get determinant
    float d = this->det();

    // Swap leading diagonal (and use determinant)
    res(0,0) = (*this)(1,1) / d;
    res(1,1) = (*this)(0,0) / d;

    // Negate counter-diagonal (and use determinant)
    res(1,0) = - (*this)(0,1) / d;
    res(0,1) = - (*this)(1,0) / d;

    return( res );
  }

private:
  //! The matrix itself
  float matrix[kMatrixSize];
};
}
}


#endif
