#pragma once

#include <stdexcept>

#include <cuda_runtime.h>

namespace kvl {
  namespace cuda {
    class CUDAException : public std::runtime_error {
    public:
      cudaError errorCode;

      CUDAException(const cudaError error) : errorCode(error),
					     runtime_error(cudaGetErrorString(errorCode)) {}
      
    };
  }
}
