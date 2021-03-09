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
 *
 */

#ifndef AFFINE_H
#define AFFINE_H

#if (__GNUC__ > 3) && !defined(HAVE_MCHECK)     // mcheck does not understand _mm_alloc et. al.
#define AFFINE_MATRIX_USE_SSE
#endif

#ifdef AFFINE_MATRIX_USE_SSE
#include <xmmintrin.h>
#endif

#include "matrix.h"

/*
  Apparently 'const' in C isn't quite the same as 'const' in C++
  So we have enums here
*/
enum { kAffineVectorSize = 4 };
enum { kAffineMatrixSize = 16 };

typedef struct _av {
  float vec[kAffineVectorSize] __attribute__ ((aligned (16)));
} AffineVector;


inline static
void SetAffineVector( AffineVector* av,
		      const float x,
		      const float y,
		      const float z ) {
  av->vec[0] = x;
  av->vec[1] = y;
  av->vec[2] = z;
  av->vec[3] = 1;
}


inline static
void GetAffineVector( const AffineVector* av,
		      float* x,
		      float* y,
		      float* z ) {
  *x = av->vec[0];
  *y = av->vec[1];
  *z = av->vec[2];
}

inline static
void GetFloorAffineVector( const AffineVector* av,
			   int *x,
			   int *y,
			   int *z ) {
  *x = (int)( floorf( av->vec[0] ) );
  *y = (int)( floorf( av->vec[1] ) );
  *z = (int)( floorf( av->vec[2] ) );
}


typedef struct _am {
  float mat[kAffineMatrixSize] __attribute__ ((aligned (16)));
} AffineMatrix;


inline static
void SetAffineMatrix( AffineMatrix* am,
		      const MATRIX* src ) {
  
  unsigned int i, j;

  if( (src->rows!=(int)kAffineVectorSize) ||
      (src->cols!=(int)kAffineVectorSize) ) {
    fprintf( stderr, "%s: Bad matrix size\n", __FUNCTION__ );
    abort();
  }

  for( i=0; i<kAffineVectorSize; i++ ) {
    for( j=0; j<kAffineVectorSize; j++ ) {
      am->mat[i+(kAffineVectorSize*j)] =
	src->data[j+(i*kAffineVectorSize)];
    }
  }
}


inline static
void GetAffineMatrix( MATRIX* dst,
		      const AffineMatrix* am ) {
  unsigned int i, j;

  if( (dst->rows!=(int)kAffineVectorSize) ||
      (dst->cols!=(int)kAffineVectorSize) ) {
    fprintf( stderr, "%s: Bad matrix size\n", __FUNCTION__ );
    abort();
  }

  for( i=0; i<kAffineVectorSize; i++ ) {
    for( j=0; j<kAffineVectorSize; j++ ) {
      dst->data[j+(i*kAffineVectorSize)] =
	am->mat[i+(kAffineVectorSize*j)];
    }
  }

}

inline static
void AffineMatrixFree( AffineMatrix **am ) {
  if( *am != NULL ) {
    free( *am );
    *am = NULL;
  }
}

inline static
void AffineMatrixAlloc( AffineMatrix **am ) {

  void* tmp = NULL;     // gcc complains otherwise
  posix_memalign(&tmp, 16, sizeof(AffineMatrix));

  if( tmp == NULL ) {
    fprintf( stderr, "%s: FAILED\n", __FUNCTION__ );
    abort();
  }
  
  AffineMatrixFree( am );

  *am = (AffineMatrix *)tmp;
}


inline static
AffineMatrix* AffineMatrixCopy( const AffineMatrix *src,
				AffineMatrix* dst ) {
  if( dst == NULL ) {
    AffineMatrixAlloc( &dst );
  }

  memcpy( dst->mat, src->mat, kAffineMatrixSize*sizeof(float) );

  return( dst );
}

inline static
void AffineMV( AffineVector * const y,
	       const AffineMatrix * const A,
	       const AffineVector * const x ) {
  
  const float * const Amat = &(A->mat[0]);
  const float * const xVec = &(x->vec[0]);
  float *yVec = &(y->vec[0]);

#ifdef AFFINE_MATRIX_USE_SSE
  __m128 resLine, col, val;

  val = _mm_set1_ps( xVec[0] );
  col = _mm_load_ps( &(Amat[0]) );
  resLine = _mm_mul_ps( col, val );

  val = _mm_set1_ps( xVec[1] );
  col = _mm_load_ps( &(Amat[4]) );
  resLine = _mm_add_ps( _mm_mul_ps( col, val ), resLine );

  val = _mm_set1_ps( xVec[2] );
  col = _mm_load_ps( &(Amat[8]) );
  resLine = _mm_add_ps( _mm_mul_ps( col, val ), resLine );

  val = _mm_set1_ps( xVec[3] );
  col = _mm_load_ps( &(Amat[12]) );
  resLine = _mm_add_ps( _mm_mul_ps( col, val ), resLine );

  _mm_store_ps( yVec, resLine );
#else
  unsigned int i, j;

  for( i=0; i<kAffineVectorSize; i++ ) {
    yVec[i] = 0;
  }

  for( j=0; j<kAffineVectorSize; j++ ) {
    for( i=0; i<kAffineVectorSize; i++ ) {
      yVec[i] += Amat[i+(j*kAffineVectorSize)] * xVec[j];
    }
  }
#endif
}



inline static
void AffineMM( AffineMatrix * const C,
	       const AffineMatrix * const A,
	       const AffineMatrix * const B ) {
#ifdef AFFINE_MATRIX_USE_SSE
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
#else
  unsigned int i, j, k;
  
  const float * const Amat = &(A->mat[0]);
  const float * const Bmat = &(B->mat[0]);
  float * const Cmat = &(C->mat[0]);

  for( i=0; i<kAffineMatrixSize; i++ ) {
    Cmat[i] = 0;
  }


  for( j=0; j<kAffineVectorSize; j++ ) {
    for( k=0; k<kAffineVectorSize; k++ ) {
      for( i=0; i<kAffineVectorSize; i++ ) {
	Cmat[i+(j*kAffineVectorSize)] +=
	  Amat[i+(k*kAffineVectorSize)] * Bmat[k+(j*kAffineVectorSize)];
      }
    }
  }
#endif
}

#endif
