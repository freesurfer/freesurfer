#pragma once

#include "cudaexception.hpp"

#define CUDA_SAFE_CALL( call ) do {		\
    cudaError err = call;			\
    if( cudaSuccess != err ) {			\
      throw kvl::cuda::CUDAException(err);	\
    }						\
    err = cudaDeviceSynchronize();		\
    if( cudaSuccess != err ) {			\
      throw kvl::cuda::CUDAException(err);	\
    }						\
  } while( 0 );
