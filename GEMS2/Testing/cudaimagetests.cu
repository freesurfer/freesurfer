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
__global__
void PlusKernel2D( kvl::cuda::Image_GPU<ElementType,2,IndexType> dst,
		   const kvl::cuda::Image_GPU<ElementType,2,IndexType> src,
		   const ElementType value ) {
  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  if( dst.PointInRange(iy,ix) ) {
    dst(iy,ix) = src(iy,ix) + value;
  }
}

template<typename ElementType, typename IndexType>
__global__
void PlusKernel3D( kvl::cuda::Image_GPU<ElementType,3,IndexType> dst,
		   const kvl::cuda::Image_GPU<ElementType,3,IndexType> src,
		   const ElementType value ) {
  const unsigned int bx = ( blockIdx.x * blockDim.x );
  const unsigned int by = ( blockIdx.y * blockDim.y );
  const unsigned int ix = threadIdx.x + bx;
  const unsigned int iy = threadIdx.y + by;

  if( dst.PointInRange(0,iy,ix) ) {
    for( unsigned int iz=0; iz<src.dims[0]; iz++ ) {
      dst(iz,iy,ix) = src(iz,iy,ix) + value;
    }
  }
}

// ===================================================

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
  err = cudaDeviceSynchronize(); if( cudaSuccess != err ) {
    throw kvl::cuda::CUDAException(err);
  }
}

template<typename ElementType, typename IndexType>
void LaunchPlusKernel( kvl::cuda::CudaImage<ElementType,2,IndexType>& dst,
		       const kvl::cuda::CudaImage<ElementType,2,IndexType>& src,
		       const ElementType value ) {
  const unsigned int kKernelWidth = 32;
  const unsigned int kKernelHeight = 8;

  dst.SetDimensions(src.GetDimensions());

  dim3 grid, threads;
  threads.x = kKernelWidth;
  threads.y = kKernelHeight;
  threads.z = 1;

  auto dims = src.GetDimensions();

  grid.x = dims[1] / kKernelWidth;
  if( (dims[1]%kKernelWidth) != 0 ) {
    grid.x++;
  }

  grid.y = dims[0] / kKernelHeight;
  if( (dims[0]%kKernelHeight) != 0 ) {
    grid.y++;
  }

  grid.z = 1;

  auto err = cudaGetLastError();
  if( cudaSuccess != err ) {
    throw kvl::cuda::CUDAException(err);
  }
  PlusKernel2D<<<grid,threads>>>( dst.getArg(), src.getArg(), value );
  err = cudaDeviceSynchronize(); if( cudaSuccess != err ) {
    throw kvl::cuda::CUDAException(err);
  }
}

template<typename ElementType, typename IndexType>
void LaunchPlusKernel( kvl::cuda::CudaImage<ElementType,3,IndexType>& dst,
		       const kvl::cuda::CudaImage<ElementType,3,IndexType>& src,
		       const ElementType value ) {
  const unsigned int kKernelWidth = 32;
  const unsigned int kKernelHeight = 8;

  dst.SetDimensions(src.GetDimensions());

  dim3 grid, threads;
  threads.x = kKernelWidth;
  threads.y = kKernelHeight;
  threads.z = 1;

  auto dims = src.GetDimensions();

  grid.x = dims[2] / kKernelWidth;
  if( (dims[2]%kKernelWidth) != 0 ) {
    grid.x++;
  }

  grid.y = dims[1] / kKernelHeight;
  if( (dims[1]%kKernelHeight) != 0 ) {
    grid.y++;
  }

  grid.z = 1;

  auto err = cudaGetLastError();
  if( cudaSuccess != err ) {
    throw kvl::cuda::CUDAException(err);
  }
  PlusKernel3D<<<grid,threads>>>( dst.getArg(), src.getArg(), value );
  err = cudaDeviceSynchronize(); if( cudaSuccess != err ) {
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

// -----------

template<>
void runPlusTest<unsigned char,2,size_t>( kvl::cuda::CudaImage<unsigned char,2,size_t>& dst,
					  const kvl::cuda::CudaImage<unsigned char,2,size_t>& src,
					  const unsigned char value ) {
  LaunchPlusKernel(dst, src, value);
}

template<>
void runPlusTest<int,2,size_t>( kvl::cuda::CudaImage<int,2,size_t>& dst,
				const kvl::cuda::CudaImage<int,2,size_t>& src,
				const int value ) {
  LaunchPlusKernel(dst, src, value);
}

template<>
void runPlusTest<float,2,size_t>( kvl::cuda::CudaImage<float,2,size_t>& dst,
				  const kvl::cuda::CudaImage<float,2,size_t>& src,
				  const float value ) {
  LaunchPlusKernel(dst, src, value);
}

template<>
void runPlusTest<double,2,size_t>( kvl::cuda::CudaImage<double,2,size_t>& dst,
				   const kvl::cuda::CudaImage<double,2,size_t>& src,
				   const double value ) {
  LaunchPlusKernel(dst, src, value);
}

// -----------

template<>
void runPlusTest<unsigned char,3,size_t>( kvl::cuda::CudaImage<unsigned char,3,size_t>& dst,
					  const kvl::cuda::CudaImage<unsigned char,3,size_t>& src,
					  const unsigned char value ) {
  LaunchPlusKernel(dst, src, value);
}

template<>
void runPlusTest<int,3,size_t>( kvl::cuda::CudaImage<int,3,size_t>& dst,
				const kvl::cuda::CudaImage<int,3,size_t>& src,
				const int value ) {
  LaunchPlusKernel(dst, src, value);
}

template<>
void runPlusTest<float,3,size_t>( kvl::cuda::CudaImage<float,3,size_t>& dst,
				  const kvl::cuda::CudaImage<float,3,size_t>& src,
				  const float value ) {
  LaunchPlusKernel(dst, src, value);
}

template<>
void runPlusTest<double,3,size_t>( kvl::cuda::CudaImage<double,3,size_t>& dst,
				   const kvl::cuda::CudaImage<double,3,size_t>& src,
				   const double value ) {
  LaunchPlusKernel(dst, src, value);
}
