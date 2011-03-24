/*! \file
  File to hold some useful CUDA macros
  Largely taken from the SDK
*/

#ifndef CUDACHECK_H
#define CUDACHECK_H

#if defined(__cplusplus)
#include <cstdio>
#include <cstdlib>
#else
#include <stdio.h>
#include <stdlib.h>
#endif


//! Macro to wrap up a call to the CUDA API. Forces synchronisation
#define CUDA_SAFE_CALL( call ) do {                                     \
    cudaError err = call;						\
    if( cudaSuccess != err ) {						\
      fprintf( stderr, "CUDA Error in file '%s' on line %i : %s.\n",	\
	       __FILE__, __LINE__, cudaGetErrorString( err ) );		\
      abort();                                                          \
    }									\
    err = cudaThreadSynchronize();					\
    if( cudaSuccess != err ) {						\
      fprintf( stderr, "CUDA Error in file '%s' on line %i : %s.\n",	\
	       __FILE__, __LINE__, cudaGetErrorString( err ) );		\
      abort();                                                          \
    }									\
  } while( 0 );

//! Synchronising macro to check for CUDA errors
#define CUDA_CHECK_ERROR( errorMessage ) do {	\
    cudaError err = cudaGetLastError();	\
    if( cudaSuccess != err) {						\
      fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",	\
	      errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) ); \
      abort();                                                          \
    }									\
    err = cudaThreadSynchronize();					\
    if( cudaSuccess != err) {						\
      fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",	\
	      errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) ); \
      abort();                                                          \
    }									\
  } while( 0 );


#if CUDA_FORCE_SYNC
#define CUDA_SAFE_CALL_ASYNC( call ) CUDA_SAFE_CALL( call );
#define CUDA_CHECK_ERROR_ASYNC( call ) CUDA_CHECK_ERROR( call );

#else
#define CUDA_SAFE_CALL_ASYNC( call ) do { \
    cudaError err = call;						\
    if( cudaSuccess != err ) {						\
      fprintf( stderr, "CUDA Error in file '%s' on line %i : %s.\n",	\
	       __FILE__, __LINE__, cudaGetErrorString( err ) );		\
      fprintf( stderr, "Check made asynchronously. May be reporting earlier error\n" ); \
      abort();                                                          \
    } } while( 0 );

#define CUDA_CHECK_ERROR_ASYNC( errorMessage ) do {	\
    cudaError err = cudaGetLastError();					\
    if( cudaSuccess != err) {						\
      fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",	\
	      errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) ); \
      fprintf( stderr, \
	       "Check made asynchronously. May be reporting earlier error\n" ); \
      abort();                                                          \
    }									\
  } while( 0 );


#endif


#endif
