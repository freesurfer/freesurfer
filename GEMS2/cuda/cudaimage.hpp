#pragma once

#include <stdexcept>
#include <memory>
#include <string>
#include <iostream>
#include <algorithm>
#include <limits>

#include <cuda_runtime.h>

#include "dimensioncuda.hpp"

#include "cudadeleters.hpp"

namespace kvl {
  namespace cuda {
    template<typename ElementType,unsigned char nDims,typename IndexType = size_t>
    class CudaImage {
    public:
      typedef Dimension<nDims,IndexType> DimensionType;
      
      void Recv( std::vector<ElementType>& dst, DimensionType& d ) const {
	throw std::runtime_error("Recv not implemented");
      }

      void Send( const std::vector<ElementType>& src, const DimensionType& d ) {
	throw std::runtime_error("Send not implemented");
      }

    private:
      DimensionType dims;
      std::unique_ptr<cudaPitchedPtr,CudaDeviceDeleter> d_elements;
    };    
  }
}
