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
	throw std::runtime_error("CudaImage::Recv not implemented");
      }

      void Send( const std::vector<ElementType>& source, const DimensionType& d ) {
	this->SetDimensions(d);

	// Setup the extent
	cudaExtent extent = this->GetCudaExtent();
	
	cudaPitchedPtr src = this->GetCudaPitchedPtr(source);
	
	// Set up the copy params
	cudaMemcpy3DParms cpyPrms = cudaMemcpy3DParms();
	cpyPrms.srcPtr = src;
	cpyPrms.dstPtr = *(this->d_elements);
	cpyPrms.extent = extent;
	cpyPrms.kind = cudaMemcpyHostToDevice;
	
	// Do the copy
	CUDA_SAFE_CALL( cudaMemcpy3D(&cpyPrms) );
      }
      
      void SetDimensions(const DimensionType& srcDims) {
	this->dims = srcDims;
	
	this->SetDimensionsCommon();
      }

      cudaExtent GetCudaExtent() const {
	// Creates the CudaExtent from the dims
	cudaExtent res;
	res.width = this->dims[nDims-1] * sizeof(ElementType);
	res.height = 1;
	res.depth = 1;
	
	if( nDims >= 2 ) {
	  res.height = this->dims[nDims-2];
	  
	  for( int i=0; i<nDims-2; i++ ) {
	    res.depth *= this->dims[i];
	  }
	}
	
	return res;
      }

      cudaPitchedPtr GetCudaPitchedPtr( const std::vector<ElementType>& target ) const {
	cudaPitchedPtr res;
      
	auto tmpExtent = this->GetCudaExtent();

	/*
	  This is not nice, but unfortunately the cudaMemcpy3D API
	  doesn't allow us to do it in a nicer way. If there were
	  a const_cudaPitchedPtr, the cast wouldn't be necessary
	*/
	res.ptr = const_cast<ElementType*>(&(target.at(0)));
	
	res.pitch = tmpExtent.width;
	res.xsize = this->dims[nDims-1];
	res.ysize = tmpExtent.height;
	
	return res;
      }

      void SetMemory( const int value ) {
	CUDA_SAFE_CALL( cudaMemset3D( *(this->d_elements), value, this->GetCudaExtent() ) ); 
      }

    private:
      DimensionType dims;
      std::unique_ptr<cudaPitchedPtr,CudaDeviceDeleter> d_elements;

      void SetDimensionsCommon() {
	// Assumes that this->dims has been set
	cudaExtent tmpExtent = this->GetCudaExtent();

	this->d_elements = std::unique_ptr<cudaPitchedPtr,CudaDeviceDeleter>(new cudaPitchedPtr);

	CUDA_SAFE_CALL( cudaMalloc3D(this->d_elements.get(), tmpExtent) );
      }
    };    
  }
}
