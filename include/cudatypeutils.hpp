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
 *    $Author: rge21 $
 *    $Date: 2010/02/26 16:34:50 $
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
			 const unsigned int b )  {
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
			   const unsigned int b ) {
  /*!
    This isn't strictly a routine for CUDA types
    but it's the sort of thing which comes up
    repeatedly in CUDA code
  */

  unsigned int tmp = DivRoundUp( a, b );

  return( b * tmp );
}


//! Equality for dim3
static
__device__ __host__
bool operator==( const dim3& a, const dim3& b ) {

  bool res = ( a.x == b.x );
  res = res && ( a.y == b.y );
  res = res && ( a.z == b.z );

  return( res );
}

//! Inequality for dim3
static
__device__ __host__
bool operator!=( const dim3& a, const dim3& b ) {
  return( !(a==b) );
}

// =======================================================
// float3 operators

//! Unary negation for float3
static
__device__ __host__
float3 operator-( const float3& a ) {
  float3 b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;

  return( b );
}

// =======================================================
// cudaExtent operators


//! Copy a dim3 to a cudaExtent
static
__host__
cudaExtent ExtentFromDims( const dim3 dims ) {
  cudaExtent e;

  e.width = dims.x;
  e.height = dims.y;
  e.depth = dims.z;

  return( e );
}



#endif
