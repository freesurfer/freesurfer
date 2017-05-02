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

    // Add one to the max to make the loops
    // simpler later
    // It also avoids some pathological cases of
    // planar tetrahedra with all integral vertices
    min[threadIdx.x] = floor(minVal);
    max[threadIdx.x] = ceil(maxVal)+1;
  }
  __syncthreads();

  // Compute barycentric co-ordinate conversion matrix
  // TODO!!!!!!!

  // Figure out how to cover the bounding box with the current thread block
  // We assume that each thread block is strictly 2D

  unsigned short box[nDims];
  for( unsigned int i=0; i<nDims; i++ ) {
    box[i] = max[i] - min[i];
  }
  const unsigned short nBx = (box[0] / blockDim.x)+1;
  const unsigned short nBy = (box[1] / blockDim.y)+1;

  // Divide the bounding box into blocks equal to the blockDim
  for( unsigned short iBy=0; iBy<nBy; iBy++ ) {
    for( unsigned short iBx=0; iBx<nBx; iBx++ ) {
      const unsigned short ix = min[0] + (iBx*blockDim.x) + threadIdx.x;
      const unsigned short iy = min[1] + (iBy*blockDim.y) + threadIdx.y;

      // Could probably do this test a little better
      if( output.PointInRange(0,iy,ix) ) {

	for( unsigned short iz=min[2]; iz<max[2]; iz++ ) {
	  bool inside = true;
	  
	  // Figure out if point lies inside tetrahedron
	  // TODO!
	  
	  if( inside ) {
	    atomicAdd(&output(iz,iy,ix),1);
	  }
	}
      }
    }
  }
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
	grid.y = (nTetrahedra / grid.x)+1;
	if( (grid.y * grid.x) < nTetrahedra ) {
	  grid.y++;
	}
      } else {
	grid.x = nTetrahedra;
	grid.y = 1;
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
