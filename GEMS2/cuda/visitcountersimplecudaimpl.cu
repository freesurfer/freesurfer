#include <stdexcept>

#include "visitcountersimplecudaimpl.hpp"

const unsigned int nDims = 3;
const unsigned int nVertices = 4;

template<typename ArgType,typename InvertType>
class SimpleSharedTetrahedron {
public:
  __device__
  void LoadAndBoundingBox( const kvl::cuda::Image_GPU<ArgType,3,size_t> tetrahedra,
			   const size_t iTet,
			   ArgType tetrahedron[nVertices][nDims],
			   unsigned short min[nDims],
			   unsigned short max[nDims] ) {
    // We presume that the tetrahedron, min and max are in shared memory
    if( (threadIdx.x < nDims) && (threadIdx.y==0) ) {
      for( unsigned int iVert=0; iVert<nVertices; iVert++ ) {
	tetrahedron[iVert][threadIdx.x] = tetrahedra(iTet,iVert,threadIdx.x);
      }
      
      // No need to sync since we're only using the first 3 threads
      ArgType minVal = tetrahedron[0][threadIdx.x];
      ArgType maxVal = tetrahedron[0][threadIdx.x];
      for( unsigned int iVert=1; iVert<nVertices; iVert++ ) {
	ArgType nxt = tetrahedron[iVert][threadIdx.x];
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
  }
private:
};

template<typename T,typename Internal>
__global__
void SimpleVisitCounterKernel( kvl::cuda::Image_GPU<int,3,unsigned short> output,
			       const kvl::cuda::Image_GPU<T,3,size_t> tetrahedra ) {
  const size_t iTet = blockIdx.x + (gridDim.x * blockIdx.y);
  
  // Check if this block has an assigned tetrahedron
  if( iTet >= tetrahedra.dims[0] ) {
    return;
  }

  // Load the tetrahedron and determine bounding box
  __shared__ T tetrahedron[nVertices][nDims];
  __shared__ unsigned short min[nDims], max[nDims];
  __shared__ T M[nDims][nDims];
  SimpleSharedTetrahedron<T,Internal> tet;

  tet.LoadAndBoundingBox( tetrahedra, iTet, tetrahedron, min, max );

  // Compute barycentric co-ordinate conversion matrix
  // This is taken from kvlTetrahedronInteriorConstIterator.hxx

  // Do single threaded
  if( (threadIdx.x==0) && (threadIdx.y==0) ) {
    // Start by computing locations relative to the first vertex
    // of the tetrahedron
    const T a = tetrahedron[1][0] - tetrahedron[0][0];
    const T b = tetrahedron[2][0] - tetrahedron[0][0];
    const T c = tetrahedron[3][0] - tetrahedron[0][0];
    const T d = tetrahedron[1][1] - tetrahedron[0][1];
    const T e = tetrahedron[2][1] - tetrahedron[0][1];
    const T f = tetrahedron[3][1] - tetrahedron[0][1];
    const T g = tetrahedron[1][2] - tetrahedron[0][2];
    const T h = tetrahedron[2][2] - tetrahedron[0][2];
    const T i = tetrahedron[3][2] - tetrahedron[0][2];
    
    const T A = ( e * i - f * h );
    const T D = -( b * i - c * h );
    const T G = ( b * f - c * e );
    const T B = -(d * i - f * g );
    const T E = ( a * i - c * g );
    const T H = -( a * f - c * d );
    const T C = ( d * h - e * g );
    const T F = - (a * h - b * g );
    const T I = ( a * e - b * d );

    const T determinant = a * A + b * B + c * C;

    M[0][0] = A / determinant;
    M[1][0] = B / determinant;
    M[2][0] = C / determinant;
    M[0][1] = D / determinant;
    M[1][1] = E / determinant;
    M[2][1] = F / determinant;
    M[0][2] = G / determinant;
    M[1][2] = H / determinant;
    M[2][2] = I / determinant;
  }
  __syncthreads();

  // Figure out how to cover the bounding box with the current thread block
  // We assume that each thread block is strictly 2D

  // Divide the bounding box into blocks equal to the blockDim
  for( unsigned short iyStart=min[1]; iyStart<max[1]; iyStart += blockDim.y ) {
    for( unsigned short ixStart=min[0]; ixStart<max[0]; ixStart += blockDim.x ) {
      const unsigned short ix = ixStart + threadIdx.x;
      const unsigned short iy = iyStart + threadIdx.y;

      // Could probably do this test a little better
      if( output.PointInRange(0,iy,ix) ) {

	for( unsigned short iz=min[2]; iz<max[2]; iz++ ) {
	  bool inside = true;
	  
	  // Figure out if point lies inside tetrahedron
	  T r[nDims];
	  r[0] = ix - tetrahedron[0][0];
	  r[1] = iy - tetrahedron[0][1];
	  r[2] = iz - tetrahedron[0][2];

	  T p[nDims];
	  for( unsigned int i=0; i<nDims; i++ ) {
	    p[i] = 0;
	    for( unsigned int j=0; j<nDims; j++ ) {
	      p[i] += M[i][j] * r[j];
	    }
	  }

	  // p now contains three of the barycentric co-ordinates
	  // We have the additional constraint that all four barycentric
	  // co-ordinates must sum to 1
	  // The point is inside the tetrahedron if all barycentric
	  // co-ordinates lie between 0 and 1
	  for( unsigned int i=0; i<nDims; i++ ) {
	    inside = inside && ( p[i] > 0 ) && (p[i] <= 1);
	  }

	  // Check the 4th (uncomputed) co-ordinate
	  inside = inside && ( (p[0]+p[1]+p[2]) <= 1 );
	  
	  // TODO Handle special cases (see IsOutsideTetrahedron method
	  // of kvlTetrahedronInteriorConstIterator.hxx
	  
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
    template<typename T,typename Internal>
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
      SimpleVisitCounterKernel<T,Internal><<<grid,threads>>>( d_output.getArg(), d_tetrahedra.getArg() );
      err = cudaDeviceSynchronize();
      if( cudaSuccess != err ) {
	throw CUDAException(err);
      }
    }

    // -----------------------------------------------------------

    template<>
    void RunVisitCounterSimpleCUDA<float,float>( CudaImage<int,3,unsigned short>& d_output,
						 const CudaImage<float,3,size_t>& d_tetrahedra ) {
      if( d_tetrahedra.GetDimensions()[1] != nVertices ) {
	throw std::runtime_error("Must have four vertices per tetrahedron!");
      }
      if( d_tetrahedra.GetDimensions()[2] != nDims ) {
	throw std::runtime_error("Only implemented for 3D space");
      }

      SimpleVisitCounter<float,float>( d_output, d_tetrahedra );
    }

    template<>
    void RunVisitCounterSimpleCUDA<double,double>( CudaImage<int,3,unsigned short>& d_output,
						   const CudaImage<double,3,size_t>& d_tetrahedra ) {
      if( d_tetrahedra.GetDimensions()[1] != nVertices ) {
	throw std::runtime_error("Must have four vertices per tetrahedron!");
      }
      if( d_tetrahedra.GetDimensions()[2] != nDims ) {
	throw std::runtime_error("Only implemented for 3D space");
      }

      SimpleVisitCounter<double,double>( d_output, d_tetrahedra );
    }
  }
}
