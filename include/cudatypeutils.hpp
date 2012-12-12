/**
 * @file  cudatypeutils.hpp
 * @brief Holds utility routines for the native CUDA types
 *
 * Holds things like mathematical operators and stream operators
 * for the native CUDA types
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.8 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef CUDA_TYPE_UTILS_H
#define CUDA_TYPE_UTILS_H

#include <cmath>
#include <iostream>

// =======================================================
// Stream operators

//! Stream insertion for float3
std::ostream& operator<<( std::ostream& os, const float3& me );
//! Stream insertion for dim3
std::ostream& operator<<( std::ostream& os, const dim3& me );


// =======================================================
// dim3 operators

//! Returns a / b, rounded up to next integer
template<typename T>
__device__ __host__
unsigned int DivRoundUp( const T a,
                         const unsigned int b )
{
  /*!
    This isn't strictly a routine just for CUDA types,
    but it's the sort of code which comes up repeatedly
    when configuring kernel launches
  */
  float tmp;

  tmp = static_cast<float>(a) / b;

  return( static_cast<unsigned int>( ceilf( tmp ) ) );
}

//! Returns the next multiple of b which is greater than a
template<typename T>
__device__ __host__
unsigned int NextMultiple( const T a,
                           const unsigned int b )
{
  /*!
    This isn't strictly a routine for CUDA types
    but it's the sort of thing which comes up
    repeatedly in CUDA code
  */

  unsigned int tmp = DivRoundUp( a, b );

  return( b * tmp );
}


//! Equality for dim3
static inline
__device__ __host__
bool operator==( const dim3& a, const dim3& b )
{

  bool res = ( a.x == b.x );
  res = res && ( a.y == b.y );
  res = res && ( a.z == b.z );

  return( res );
}

//! Inequality for dim3
static inline
__device__ __host__
bool operator!=( const dim3& a, const dim3& b )
{
  return( !(a==b) );
}

// =======================================================
// float3 operators

//! Unary negation for float3
static inline
__device__ __host__
float3 operator-( const float3& a )
{
  float3 b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;

  return( b );
}


//! Self addition for float3
static inline
__device__ __host__
void operator+=( float3& a, const float3 b )
{

  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
}

//! Addition for float3
static inline
__device__ __host__
const float3 operator+( const float3& a, const float3& b )
{
  float3 res = a;

  res += b;

  return( res );
}



//! Self subtraction for float3
static inline
__device__ __host__
void operator-=( float3& a, const float3 b )
{
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
}

//! Subtraction for float3
static inline
__device__ __host__
const float3 operator-( const float3& a, const float3& b )
{
  float3 res = a;

  res -= b;

  return( res );
}



//! Scalar self-multiplication for float3
static inline
__device__ __host__
void operator*=( float3& a, const float b )
{
  a.x *= b;
  a.y *= b;
  a.z *= b;
}

//! Scalar multiplication for float3 (on right)
static inline
__device__ __host__
const float3 operator*( const float3& a, const float b )
{
  float3 res = a;

  res *= b;

  return( res );
}

//! Scalar multiplication for float3 (on left)
static inline
__device__ __host__
const float3 operator*( const float b, const float3& a )
{

  return( a * b );
}



//! Scalar self-division for float3
static inline
__device__ __host__
void operator/=( float3& a, const float b )
{

  a *= (1/b);
}

//! Scalar division for float3
static inline
__device__ __host__
const float3 operator/( const float3& a, const float b )
{
  float3 res = a;

  res /= b;

  return( res );
}



//! Dot product for float3
static inline
__device__ __host__
float dot( const float3& a, const float3& b )
{
  float res;

  res = (a.x*b.x) + (a.y*b.y) + (a.z*b.z);

  return( res );
}

//! Modulus squared for float3
static inline
__device__ __host__
float modsq( const float3& a )
{
  return( dot( a, a ) );
}

//! Modulus for float3
static inline
__device__ __host__
float mod( const float3& a )
{
  return( sqrtf( modsq( a ) ) );
}



//! Cross product for float3
static inline
__device__ __host__
const float3 cross( const float3& a, const float3& b )
{
  float3 res;

  res.x = (a.y*b.z) - (a.z*b.y);
  res.y = (a.z*b.x) - (a.x*b.z);
  res.z = (a.x*b.y) - (a.y*b.x);

  return( res );
}




//! Scalar triple product \f$ \mathbf{a} \times \mathbf{b} \cdot \mathbf{c} \f$
static inline
__device__ __host__
float stp( const float3& a, const float3& b, const float3& c )
{
  return( dot( cross( a, b ), c ) );
}



// =======================================================
// cudaExtent operators


//! Copy a dim3 to a cudaExtent
static
__host__
cudaExtent ExtentFromDims( const dim3 dims )
{
  cudaExtent e;

  e.width = dims.x;
  e.height = dims.y;
  e.depth = dims.z;

  return( e );
}



#endif
