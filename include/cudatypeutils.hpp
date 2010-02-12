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
 *    $Date: 2010/02/12 15:50:09 $
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

#ifndef CUDA_TYPE_UTILS_H
#define CUDA_TYPE_UTILS_H

#include <iostream>

// =======================================================
// Stream operators

//! Stream insertion for float3
std::ostream& operator<<( std::ostream& os, const float3& me );
//! Stream insertion for dim3
std::ostream& operator<<( std::ostream& os, const dim3& me );


// =======================================================
// dim3 operators

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


#endif
