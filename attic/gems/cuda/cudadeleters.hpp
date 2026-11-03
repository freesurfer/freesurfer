#pragma once

#include <memory>

#include "cudacheck.hpp"


namespace kvl {
  namespace cuda {
    
    class CudaDeviceDeleter {
    public:
      void operator()(cudaPitchedPtr* d_ptr) {
	if( d_ptr != NULL ) {
	  if( d_ptr->ptr != NULL ) {
	    CUDA_SAFE_CALL(cudaFree(d_ptr->ptr));
	  }
	  delete d_ptr;
	}
      }
    };
  }
}
