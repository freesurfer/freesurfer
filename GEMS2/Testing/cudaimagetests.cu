#include "cudaimagetests.hpp"

template<typename ElementType, typename IndexType>
__global__
void PlusKernel1D( kvl::cuda::Image_GPU<ElementType,1,IndexType> dst,
		   const kvl::cuda::Image_GPU<ElementType,1,IndexType> src,
		   const ElementType value ) {
  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int ix = threadIdx.x + bx;

  if( dst.PointInRange(ix) ) {
    dst(ix) = src(ix) + value;
  }
}

template<typename ElementType, typename IndexType>
void LaunchPlusKernel( kvl::cuda::CudaImage<ElementType,1,IndexType>& dst,
		       const kvl::cuda::CudaImage<ElementType,1,IndexType>& src,
		       const ElementType value ) {
  const unsigned int kKernelSize = 128;

  dst.SetDimensions(src.GetDimensions());

  dim3 grid, threads;
  threads.x = kKernelSize;
  threads.y = threads.z = 1;

  const IndexType nx = src.GetDimensions()[0];

  grid.x = nx / kKernelSize;
  if( (nx%kKernelSize) != 0 ) {
    grid.x++;
  }
  grid.y = grid.z = 1;

  auto err = cudaGetLastError();
  if( cudaSuccess != err ) {
    throw kvl::cuda::CUDAException(err);
  }
  PlusKernel1D<<<grid,threads>>>( dst.getArg(), src.getArg(), value );
  err = cudaThreadSynchronize(); if( cudaSuccess != err ) {
    throw kvl::cuda::CUDAException(err);
  }
}

// =====================================================

template<>
void runPlusTest<unsigned char,1,size_t>( kvl::cuda::CudaImage<unsigned char,1,size_t>& dst,
					  const kvl::cuda::CudaImage<unsigned char,1,size_t>& src,
					  const unsigned char value ) {
  LaunchPlusKernel(dst, src, value);
}

template<>
void runPlusTest<int,1,size_t>( kvl::cuda::CudaImage<int,1,size_t>& dst,
				const kvl::cuda::CudaImage<int,1,size_t>& src,
				const int value ) {
  LaunchPlusKernel(dst, src, value);
}

template<>
void runPlusTest<float,1,size_t>( kvl::cuda::CudaImage<float,1,size_t>& dst,
				  const kvl::cuda::CudaImage<float,1,size_t>& src,
				  const float value ) {
  LaunchPlusKernel(dst, src, value);
}

template<>
void runPlusTest<double,1,size_t>( kvl::cuda::CudaImage<double,1,size_t>& dst,
				   const kvl::cuda::CudaImage<double,1,size_t>& src,
				   const double value ) {
  LaunchPlusKernel(dst, src, value);
}
