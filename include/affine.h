/**
 * @file  affine.hpp
 * @brief Affine transformations
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/12/03 16:34:00 $
 *    $Revision: 1.2 $
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
 */

#ifndef AFFINE_H
#define AFFINE_H

#include <xmmintrin.h>

#include "matrix.h"


/*
  Apparently 'const' in C isn't quite the same as 'const' in C++
  So we have enums here
*/
enum { kAffineVectorSize = 4 };
enum { kAffineMatrixSize = 16 };

// =====================================================

typedef struct _av {
  float vec[kAffineVectorSize] __attribute__ ((aligned (16)));
} AffineVector;


inline
void SetAffineVector( AffineVector* av,
		      const float x,
		      const float y,
		      const float z ) {
  av->vec[0] = x;
  av->vec[1] = y;
  av->vec[2] = z;
  av->vec[3] = 1;
}


inline
void GetAffineVector( const AffineVector* av,
		      float* x,
		      float* y,
		      float* z ) {
  *x = av->vec[0];
  *y = av->vec[1];
  *z = av->vec[2];
}

inline
void GetFloorAffineVector( const AffineVector* av,
			   int *x,
			   int *y,
			   int *z ) {
  *x = (int)( floorf( av->vec[0] ) );
  *y = (int)( floorf( av->vec[1] ) );
  *z = (int)( floorf( av->vec[2] ) );
}




// =====================================================


typedef struct _am {
  float mat[kAffineMatrixSize] __attribute__ ((aligned (16)));
} AffineMatrix;



inline
void SetAffineMatrix( AffineMatrix* am,
		      const MATRIX* src ) {
  
  unsigned int i, j;

  if( (src->rows!=(int)kAffineVectorSize) ||
      (src->cols!=(int)kAffineVectorSize) ) {
    fprintf( stderr, "%s: Bad matrix sizes\n", __FUNCTION__ );
    exit( EXIT_FAILURE );
  }

  for( i=0; i<kAffineVectorSize; i++ ) {
    for( j=0; j<kAffineVectorSize; j++ ) {
      am->mat[i+(kAffineVectorSize*j)] =
	src->data[j+(i*kAffineVectorSize)];
    }
  }
}



inline
void AffineMV( AffineVector* y,
	       const AffineMatrix* A,
	       const AffineVector* x ) {
  
  __m128 resLine, col, val;

  val = _mm_set1_ps( x->vec[0] );
  col = _mm_load_ps( &(A->mat[0] ) );
  resLine = _mm_mul_ps( col, val );

  val = _mm_set1_ps( x->vec[1] );
  col = _mm_load_ps( &(A->mat[4] ) );
  resLine = _mm_add_ps( _mm_mul_ps( col, val ), resLine );

  val = _mm_set1_ps( x->vec[2] );
  col = _mm_load_ps( &(A->mat[8] ) );
  resLine = _mm_add_ps( _mm_mul_ps( col, val ), resLine );

  val = _mm_set1_ps( x->vec[3] );
  col = _mm_load_ps( &(A->mat[12] ) );
  resLine = _mm_add_ps( _mm_mul_ps( col, val ), resLine );

  _mm_store_ps( &(y->vec[0]), resLine );
}



inline
void AffineMM( AffineMatrix* C,
	       const AffineMatrix* A,
	       const AffineMatrix* B ) {
  unsigned int j;
  __m128 resLine, col, val;

  const float *Bloc = &(B->mat[0]);
  const float *Aloc = &(A->mat[0]);
  float *cLoc = &(C->mat[0]);

  for( j=0; j<kAffineVectorSize; j++ ) {
    val = _mm_set1_ps( *Bloc );
    col = _mm_load_ps( Aloc );
    resLine = _mm_mul_ps( col, val );
    Bloc += 1;
    Aloc += kAffineVectorSize;
    
    val = _mm_set1_ps( *Bloc );
    col = _mm_load_ps( Aloc );
    resLine = _mm_add_ps( _mm_mul_ps( col, val), resLine );
    Bloc += 1;
    Aloc += kAffineVectorSize;
    
    val = _mm_set1_ps( *Bloc );
    col = _mm_load_ps( Aloc );
    resLine = _mm_add_ps( _mm_mul_ps( col, val), resLine );
    Bloc += 1;
    Aloc += kAffineVectorSize;
    
    val = _mm_set1_ps( *Bloc );
    col = _mm_load_ps( Aloc );
    resLine = _mm_add_ps( _mm_mul_ps( col, val), resLine );
    Bloc += 1;
    Aloc = &(A->mat[0]); // Reset for next iteration
      
    _mm_store_ps( cLoc, resLine );
    cLoc += kAffineVectorSize;
    
  }
}


#endif
