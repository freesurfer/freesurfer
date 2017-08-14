#pragma once

#include <stdexcept>

#include "cudaimage.hpp"


template<typename ElementType,unsigned char nDims,typename IndexType>
void runPlusTest( kvl::cuda::CudaImage<ElementType,nDims,IndexType>& dst,
		  const kvl::cuda::CudaImage<ElementType,nDims,IndexType>& src,
		  const ElementType value ) {
  throw std::runtime_error("Not implemented");
}

template<>
void runPlusTest<unsigned char,1,size_t>( kvl::cuda::CudaImage<unsigned char,1,size_t>& dst,
					  const kvl::cuda::CudaImage<unsigned char,1,size_t>& src,
					  const unsigned char value );

template<>
void runPlusTest<int,1,size_t>( kvl::cuda::CudaImage<int,1,size_t>& dst,
				const kvl::cuda::CudaImage<int,1,size_t>& src,
				const int value );

template<>
void runPlusTest<float,1,size_t>( kvl::cuda::CudaImage<float,1,size_t>& dst,
				  const kvl::cuda::CudaImage<float,1,size_t>& src,
				  const float value );

template<>
void runPlusTest<double,1,size_t>( kvl::cuda::CudaImage<double,1,size_t>& dst,
				   const kvl::cuda::CudaImage<double,1,size_t>& src,
				   const double value );

template<>
void runPlusTest<unsigned char,2,size_t>( kvl::cuda::CudaImage<unsigned char,2,size_t>& dst,
					  const kvl::cuda::CudaImage<unsigned char,2,size_t>& src,
					  const unsigned char value );

template<>
void runPlusTest<int,2,size_t>( kvl::cuda::CudaImage<int,2,size_t>& dst,
				const kvl::cuda::CudaImage<int,2,size_t>& src,
				const int value );

template<>
void runPlusTest<float,2,size_t>( kvl::cuda::CudaImage<float,2,size_t>& dst,
				  const kvl::cuda::CudaImage<float,2,size_t>& src,
				  const float value );

template<>
void runPlusTest<double,2,size_t>( kvl::cuda::CudaImage<double,2,size_t>& dst,
				   const kvl::cuda::CudaImage<double,2,size_t>& src,
				   const double value );

template<>
void runPlusTest<unsigned char,3,size_t>( kvl::cuda::CudaImage<unsigned char,3,size_t>& dst,
					  const kvl::cuda::CudaImage<unsigned char,3,size_t>& src,
					  const unsigned char value );

template<>
void runPlusTest<int,3,size_t>( kvl::cuda::CudaImage<int,3,size_t>& dst,
				const kvl::cuda::CudaImage<int,3,size_t>& src,
				const int value );

template<>
void runPlusTest<float,3,size_t>( kvl::cuda::CudaImage<float,3,size_t>& dst,
				  const kvl::cuda::CudaImage<float,3,size_t>& src,
				  const float value );

template<>
void runPlusTest<double,3,size_t>( kvl::cuda::CudaImage<double,3,size_t>& dst,
				   const kvl::cuda::CudaImage<double,3,size_t>& src,
				   const double value );
