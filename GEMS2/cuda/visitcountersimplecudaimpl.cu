#include <stdexcept>

#include "visitcountersimplecudaimpl.hpp"

template<typename T>
__global__
void SimpleVisitCounterKernel( kvl::cuda::Image_GPU<int,3,unsigned short> output,
			       const kvl::cuda::Image_GPU<T,3,size_t> tetrahedra ) {
  const size_t iTet = blockIdx.x + (gridDim.x * blockIdx.y);

}


namespace kvl {
  namespace cuda {
    template<typename T>
    void SimpleVisitCounter( CudaImage<int,3,unsigned short>& d_output,
			     const CudaImage<T,3,size_t>& d_tetrahedra ) {
      const unsigned int nBlockx = 1024;

      const unsigned int nThreadsx = 32;
      const unsigned int nThreadsy = 32;
      const unsigned int nThreadsz = 1;

      dim3 grid, threads;

      const size_t nTetrahedra = d_tetrahedra.GetDimensions()[0];
      
      if( nTetrahedra > nBlockx ) {
	grid.x = nBlockx;
	grid.y = nTetrahedra / grid.x;
	if( (grid.y * grid.x) < nTetrahedra ) {
	  grid.y++;
	}
      } else {
	grid.x = nTetrahedra;
      }

      threads.x = nThreadsx;
      threads.y = nThreadsy;
      threads.z = nThreadsz;
      
      // Run the kernel
      auto err = cudaGetLastError();
      if( cudaSuccess != err ) {
	throw CUDAException(err);
      }
      SimpleVisitCounterKernel<<<grid,threads>>>( d_output.getArg(), d_tetrahedra.getArg() );
      err = cudaThreadSynchronize();
      if( cudaSuccess != err ) {
	throw CUDAException(err);
      }
    }

    // -----------------------------------------------------------

    template<>
    void RunVisitCounterSimpleCUDA( CudaImage<int,3,unsigned short>& d_output,
				    const CudaImage<float,3,size_t>& d_tetrahedra ) {
      SimpleVisitCounter( d_output, d_tetrahedra );
      throw std::runtime_error("RunVisitCounterSimpleCUDA not yet implemented for float");
    }
  }
}
