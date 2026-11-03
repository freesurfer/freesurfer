#pragma once

#include <stdexcept>
#include <memory>
#include <string>
#include <iostream>
#include <algorithm>
#include <limits>
#include <vector>

#include <cuda_runtime.h>

#include "dimensioncuda.hpp"

#include "cudadeleters.hpp"

namespace kvl {
  namespace cuda {
    template<typename ElementType,unsigned char nDims,typename IndexType = size_t>
    class Image_GPU {
    public:
      IndexType dims[nDims];
      void* pitchedPtr;
      size_t dataPitch;

      template<typename NextType, typename... Values>
      __device__
      ElementType operator()( const NextType iVal, const Values... vals ) const {
	static_assert( std::numeric_limits<IndexType>::is_integer, "Must have integral IndexType" );
	IndexType location[nDims];
	this->copyToArray(location, iVal, vals...);
	
	IndexType iRow = 0;
	if( nDims>=2 ) {
	  for( unsigned char i=0; i<=nDims-2; i++ ) {
	    iRow = location[i] + (iRow * this->dims[i]);
	  }
	}

	const char* data = reinterpret_cast<const char*>(this->pitchedPtr);
	const char* row = data + (iRow*this->dataPitch);
	
	return( reinterpret_cast<const ElementType*>(row)[location[nDims-1]] );
      }

      template<typename NextType, typename... Values>
      __device__
      ElementType& operator()( const NextType iVal, const Values... vals ) {
	static_assert( std::numeric_limits<IndexType>::is_integer, "Must have integral IndexType" );
	IndexType location[nDims];
	this->copyToArray(location, iVal, vals...);

	IndexType iRow = 0;
	if( nDims>=2 ) {
	  for( unsigned char i=0; i<=nDims-2; i++ ) {
	    iRow = location[i] + (iRow * this->dims[i]);
	  }
	}
	
	char* data = reinterpret_cast<char*>(this->pitchedPtr);
	char* row = data + (iRow*this->dataPitch);
	
	return( reinterpret_cast<ElementType*>(row)[location[nDims-1]] );
      }

      
      template<typename NextType, typename... Values>
      __device__
      bool PointInRange(const NextType iVal, const Values... vals) const {
	static_assert((1+sizeof...(Values))==nDims,
		      "Must call with nDims arguments");
	IndexType location[nDims];
      
	this->copyToArray(location, iVal, vals...);
      
	return this->PointInRangeFromArray(location);
      }

      __device__
      bool PointInRangeFromArray(const IndexType location[nDims]) const {
	bool result = true;
	
	for( unsigned char i=0; i<nDims; i++ ) {
	  result = result && (location[i] >=0) && (location[i] < this->dims[i]);
	}

	return result;
      }

    private:
      // Base case of copyToArray
      __device__
      void copyToArray(IndexType*) const { }

      // Extracts variadic template arguments into a array of length nDims
      template<typename NextType, typename... Values>
      __device__
      void copyToArray(IndexType* nextLoc, const NextType iVal, const Values... vals) const {
	*nextLoc = iVal;
	this->copyToArray(++nextLoc,vals...);
      }
    };

    // --------------------------------------------

    template<typename ElementType,unsigned char nDims,typename IndexType = size_t>
    class CudaImage {
    public:
      typedef Dimension<nDims,IndexType> DimensionType;
      
      void Recv( std::vector<ElementType>& dest, DimensionType& dims ) const {
	dims = this->dims;

	// Set up the result
	dest.resize(dims.ElementCount());

	// Setup the extent
	auto extent = this->GetCudaExtent();
	
	auto dst = this->GetCudaPitchedPtr(dest);
      
	// Set up the copy params
	cudaMemcpy3DParms cpyPrms = cudaMemcpy3DParms();
	cpyPrms.srcPtr = *(this->d_elements);
	cpyPrms.dstPtr = dst;
	cpyPrms.extent = extent;
	cpyPrms.kind = cudaMemcpyDeviceToHost;
	
	// Do the copy
	CUDA_SAFE_CALL( cudaMemcpy3D(&cpyPrms) );
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
      
      DimensionType GetDimensions() const {
	return this->dims;
      }

      size_t ElementCount() const {
	return this->dims.ElementCount();
      }

      cudaExtent GetCudaExtent() const {
	// Creates the CudaExtent from the dims
	cudaExtent res;
	res.width = this->dims[nDims-1] * sizeof(ElementType);
	res.height = 1;
	res.depth = 1;
	
	if( nDims >= 2 ) {
	  res.height = this->dims[nDims-2];
	  
	  for( unsigned char i=0; i<nDims-2; i++ ) {
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

      Image_GPU<ElementType,nDims,IndexType> getArg() const {
	Image_GPU<ElementType,nDims,IndexType> gpuArg;

	for( unsigned char i=0; i<nDims; i++ ) {
	  gpuArg.dims[i] = this->dims[i];
	}
	gpuArg.pitchedPtr =  this->d_elements->ptr;
	gpuArg.dataPitch = this->d_elements->pitch;

	return gpuArg;
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
