#include <cuda_runtime.h>

#include "cudacheck.hpp"

#include "cudacontroller.hpp"

namespace kvl {
  namespace cuda {
    void InitialiseCUDA(const int deviceID) {
      CUDA_SAFE_CALL( cudaSetDevice(deviceID) );
      int* d_tmp;
      CUDA_SAFE_CALL( cudaMalloc( &d_tmp, 1 ) );
      CUDA_SAFE_CALL( cudaFree( d_tmp ) );
    }

    void FinalizeCUDA() {
    }
  }
}
