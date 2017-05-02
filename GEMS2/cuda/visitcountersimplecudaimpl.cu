#include <stdexcept>

#include "visitcountersimplecudaimpl.hpp"

const unsigned int nDims = 3;
const unsigned int nVertices = 4;

template<typename T>
__global__
void SimpleVisitCounterKernel( kvl::cuda::Image_GPU<int,3,unsigned short> output,
			       const kvl::cuda::Image_GPU<T,3,size_t> tetrahedra ) {
  const size_t iTet = blockIdx.x + (gridDim.x * blockIdx.y);
  
  // Load the tetrahedron and determine bounding box
  __shared__ T tetrahedron[nVertices][nDims];
  __shared__ unsigned short min[nDims], max[nDims];
  if( (threadIdx.x < nDims) && (threadIdx.y==0) ) {
    for( unsigned int iVert=0; iVert<nVertices; iVert++ ) {
      tetrahedron[iVert][threadIdx.x] = tetrahedra(iTet,iVert,threadIdx.x);
    }

    // No need to sync since we're only using the first 3 threads
    T minVal = tetrahedron[0][threadIdx.x];
    T maxVal = tetrahedron[0][threadIdx.x];
    for( unsigned int iVert=1; iVert<nVertices; iVert++ ) {
      T nxt = tetrahedron[iVert][threadIdx.x];
      if( nxt < minVal ) {
	minVal = nxt;
      }
      if( nxt > maxVal ) {
	maxVal = nxt;
      }
    }
    
    // Indices are always positive
    if( minVal < 0 ) {
      minVal = 0;
    }

    min[threadIdx.x] = floor(minVal);
    max[threadIdx.x] = ceil(maxVal);
  }
  __syncthreads();
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
      throw std::runtime_error("SimpleVisitCounterKernel not yet implemented");
    }

    // -----------------------------------------------------------

    template<>
    void RunVisitCounterSimpleCUDA( CudaImage<int,3,unsigned short>& d_output,
				    const CudaImage<float,3,size_t>& d_tetrahedra ) {
      if( d_tetrahedra.GetDimensions()[1] != nVertices ) {
	throw std::runtime_error("Must have four vertices per tetrahedron!");
      }
      if( d_tetrahedra.GetDimensions()[2] != nDims ) {
	throw std::runtime_error("Only implemented for 3D space");
      }

      SimpleVisitCounter( d_output, d_tetrahedra );
    }
  }
}
