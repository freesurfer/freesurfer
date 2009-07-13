/*! \file
  Header file to hold CUDA error checking macros
  Original Author : Richard Edgar
  $Id: cudacheck.h,v 1.1 2009/07/13 21:25:33 krish Exp $
*/

//! Wrapper for CUDA library routines
#define CUDA_SAFE_CALL( call) do {					\
    cudaError err;       						\
    err = call;								\
    if( cudaSuccess != err) {						\
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",	\
	      __FILE__, __LINE__, cudaGetErrorString( err) );		\
      exit(EXIT_FAILURE);						\
    }									\
    err = cudaThreadSynchronize();					\
    if( cudaSuccess != err) {						\
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",	\
	      __FILE__, __LINE__, cudaGetErrorString( err) );		\
      exit(EXIT_FAILURE);						\
    } } while (0)

//! Check for CUDA kernel execution error
#define CUDA_CHECK_KERNEL_ERROR(errorMessage) do {			\
    cudaError_t err = cudaGetLastError();				\
    if( cudaSuccess != err) {						\
      fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",	\
	      errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) ); \
      exit(EXIT_FAILURE);						\
    }									\
    err = cudaThreadSynchronize();					\
    if( cudaSuccess != err) {						\
      fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",	\
	      errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) ); \
      exit(EXIT_FAILURE);						\
    } } while (0)



#ifdef __DEVICE_EMULATION__

//! Do a call once per block
#define PER_BLOCK( call ) do {			\
    if( threadIdx.x == 0 ) { \
      call;		     \
    }			     \
  } while( 0 ); __syncthreads();


//! Do a call once per grid
#define PER_GRID( call, blockID ) do {		\
    if( blockIdx.x == blockID ) { \
      PER_BLOCK( call ); \
    }			 \
  } while( 0 ); __syncthreads();

#else

#define PER_BLOCK( call );
#define PER_GRID( call, blockID );

#endif
