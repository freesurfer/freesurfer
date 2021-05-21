/**
 * @brief Affine transformations
 *
 */
/*
 * Original Authors: Richard Edgar
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 */

#ifndef AFFINE_HPP
#define AFFINE_HPP

#include <xmmintrin.h>
#include <cstring>
#include <iostream>
#include <iomanip>

#include "matrix.h"

#include "affine.h"


namespace Freesurfer
{

const unsigned int kAffineVectorSize = 4;
const unsigned int kAffineMatrixSize = kAffineVectorSize * kAffineVectorSize;

// ===============================================================

//! Affine vector class
template<typename T>
class AffineVector
{

public:
  // Indexing enum
  enum ComponentIndices { iX=0, iY=1, iZ=2, iW=3 };

  // ----------------------------------------------------
  // No constructors or destructors


  // ----------------------------------------------------
  // Get and set

  template<typename U>
  void Set( const U& x, const U& y, const U& z )
  {
    this->vec[iX] = x;
    this->vec[iY] = y;
    this->vec[iZ] = z;
    this->vec[iW] = 1;
  }

  template<typename U>
  void Get( U& x, U& y, U& z ) const
  {
    x = this->vec[iX];
    y = this->vec[iY];
    z = this->vec[iZ];
  }

  void GetFloor( int& x, int& y, int& z ) const
  {
    x = static_cast<int>( floor( this->vec[iX] ) );
    y = static_cast<int>( floor( this->vec[iY] ) );
    z = static_cast<int>( floor( this->vec[iZ] ) );
  }


  // ----------------------------------------------------
  // I/O

  void Print( std::ostream& os = std::cout ) const
  {
    os << "("
       << std::setw(9) << std::setprecision(3) << this->vec[iX]
       << ","
       << std::setw(9) << std::setprecision(3) << this->vec[iY]
       << ","
       << std::setw(9) << std::setprecision(3) << this->vec[iZ]
       << " )";
  }

  // ----------------------------------------------------

  //! The actual vector, aligned on 16 byte boundary
  T vec[kAffineVectorSize] __attribute__ ((aligned (16)));


};


// ================================================================

//! Affine matrix class
template<typename T>
class AffineMatrix
{

public:
  // ----------------------------------------------------
  // No constructors or destructors


  // ----------------------------------------------------
  // Indexing operators

  T operator()( const unsigned int i, const unsigned int j ) const
  {
    return( this->mat[i+j*kAffineVectorSize] );
  }

  T& operator()( const unsigned int i, const unsigned int j )
  {
    return( this->mat[i+j*kAffineVectorSize] );
  }

  // ----------------------------------------------------
  // Setting

  void Set( const MATRIX* src )
  {

    for( unsigned int i=0; i<kAffineVectorSize; i++ )
    {
      for( unsigned int j=0; j<kAffineVectorSize; j++ )
      {
        (*this)(i,j) = src->data[j+(i*kAffineVectorSize)];
      }
    }
  }

  void Set( const ::AffineMatrix* src )
  {
    /*
    Note the descoping on AffineMatrix, since this
    is for getting data from the C version of the
    structure.
    Can't use memcpy because of possible type
    conversion
    */
    for( unsigned int i=0; i<kAffineMatrixSize; i++ )
    {
      this->mat[i] = src->mat[i];
    }
  }


  void Clear( void )
  {
    memset( &(this->mat[0]), 0, kAffineMatrixSize*sizeof(T) );
  }

  void Identity( void )
  {
    this->Clear();
    this->mat[0]  = 1;
    this->mat[5]  = 1;
    this->mat[10] = 1;
    this->mat[15] = 1;
  }

  // ---------------------------------------------------
  // Mathematical operations

  //! Matrix-Vector multiplication
  AffineVector<T> operator*( const AffineVector<T>& v ) const
  {
    AffineVector<T> res;

    for( unsigned int i=0; i<kAffineVectorSize; i++ )
    {
      T tmp = 0;
      for( unsigned int j=0; j<kAffineVectorSize; j++ )
      {
        tmp += (*this)(i,j) * v.vec[j];
      }
      res.vec[i] = tmp;
    }

    return( res );
  }

  //! Matrix-Matrix multiplication
  AffineMatrix<T> operator*( const AffineMatrix<T>& m ) const
  {
    AffineMatrix<T> res;

    for( unsigned int i=0; i<kAffineVectorSize; i++ )
    {
      for( unsigned int j=0; j<kAffineVectorSize; j++ )
      {
        T tmp = 0;
        for( unsigned int k=0; k<kAffineVectorSize; k++ )
        {
          tmp += (*this)(i,k) * m(k,j);
        }
        res(i,j) = tmp;
      }
    }

    return( res );
  }

  // ----------------------------------------------------
  // I/O

  void Print( std::ostream& os = std::cout ) const
  {

    for( unsigned int i=0; i<kAffineVectorSize; i++ )
    {
      for( unsigned int j=0; j<kAffineVectorSize; j++ )
      {
        os << std::setw(9) << std::setprecision(3) << (*this)( i, j );
      }
      os << std::endl;
    }

  }


  //! The actual matrix, column major and 16 byte aligned
  T mat[kAffineMatrixSize] __attribute__ ((aligned (16)));

};


//! Specialise matrix-vector for float and use SSE
template<>
AffineVector<float>
AffineMatrix<float>::operator*( const AffineVector<float>& v ) const
{

  __m128 resLine, col, val;

  // Use x component of input vector
  val = _mm_set1_ps( v.vec[0] );
  col = _mm_load_ps( &(this->mat[0]) );
  resLine = _mm_mul_ps( col, val );

  // Use y component of input vector
  val = _mm_set1_ps( v.vec[1] );
  col = _mm_load_ps( &(this->mat[4] ) );
  resLine = _mm_add_ps( _mm_mul_ps( col, val), resLine );

  // Use z component of input vector
  val = _mm_set1_ps( v.vec[2] );
  col = _mm_load_ps( &(this->mat[8] ) );
  resLine = _mm_add_ps( _mm_mul_ps( col, val), resLine );

  // And the w
  val = _mm_set1_ps( v.vec[3] );
  col = _mm_load_ps( &(this->mat[12] ) );
  resLine = _mm_add_ps( _mm_mul_ps( col, val), resLine );

  AffineVector<float> res;
  _mm_store_ps( &(res.vec[0]), resLine );

  return( res );
}

//! Specialise matrix-matrix for float and use SSE
template<>
AffineMatrix<float>
AffineMatrix<float>::operator*( const AffineMatrix<float>& m ) const
{

  __m128 resLine, col, val;
  AffineMatrix<float> res;
  const float *mLoc = &m.mat[0];
  const float *thisLoc = &(this->mat[0]);
  float *resLoc = &(res.mat[0]);

  // Loop over columns
  for( unsigned int j=0; j<kAffineVectorSize; j++ )
  {
    val = _mm_set1_ps( *mLoc );
    col = _mm_load_ps( thisLoc );
    resLine = _mm_mul_ps( col, val );
    mLoc += 1;
    thisLoc += kAffineVectorSize;

    val = _mm_set1_ps( *mLoc );
    col = _mm_load_ps( thisLoc );
    resLine = _mm_add_ps( _mm_mul_ps( col, val), resLine );
    mLoc += 1;
    thisLoc += kAffineVectorSize;

    val = _mm_set1_ps( *mLoc );
    col = _mm_load_ps( thisLoc );
    resLine = _mm_add_ps( _mm_mul_ps( col, val), resLine );
    mLoc += 1;
    thisLoc += kAffineVectorSize;

    val = _mm_set1_ps( *mLoc );
    col = _mm_load_ps( thisLoc );
    resLine = _mm_add_ps( _mm_mul_ps( col, val), resLine );
    mLoc += 1;
    thisLoc = &(this->mat[0]); // Reset for next iteration

    _mm_store_ps( resLoc, resLine );
    resLoc += kAffineVectorSize;
  }
  return( res );

}

}


#endif
