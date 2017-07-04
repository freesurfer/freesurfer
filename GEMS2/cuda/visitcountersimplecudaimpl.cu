#include <stdexcept>

#include "visitcountersimplecudaimpl.hpp"

#include "visitcounteraction.hpp"
#include "simplesharedtetrahedroninterior.hpp"

const unsigned int nDims = 3;
const unsigned int nVertices = 4;

// ---------

template<typename ArgType>
class SimpleMeshSupply {
public:
  __device__
  SimpleMeshSupply( const kvl::cuda::Image_GPU<ArgType,3,size_t>& tetrahedra ) : tetInfo(tetrahedra) {}

  __device__
  ArgType GetVertexCoordinate( size_t iTet, size_t iVert, size_t iDim ) const {
    return this->tetInfo(iTet,iVert,threadIdx.x);
  }

  __device__
  size_t GetTetrahedraCount() const {
    return this->tetInfo.dims[0];
  }

private:
  const kvl::cuda::Image_GPU<ArgType,3,size_t>& tetInfo;
};

// ---------

template<typename T,typename Internal>
__global__
void SimpleVisitCounterKernel( kvl::cuda::Image_GPU<int,3,unsigned short> output,
			       const kvl::cuda::Image_GPU<T,3,size_t> tetrahedra ) {
  const size_t iTet = blockIdx.x + (gridDim.x * blockIdx.y);
  
  SimpleMeshSupply<T> mesh(tetrahedra);

  // Check if this block has an assigned tetrahedron
  if( iTet >= mesh.GetTetrahedraCount() ) {
    return;
  }

  // Load the tetrahedron and determine bounding box
  __shared__ T tetrahedron[nVertices][nDims];
  __shared__ unsigned short min[nDims], max[nDims];
  __shared__ T M[nDims][nDims];
  SimpleSharedTetrahedron<T,Internal> tet(tetrahedron, M);

  tet.LoadAndBoundingBox( mesh, iTet, min, max );

  tet.ComputeBarycentricTransform();

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
	  bool inside = tet.PointInside(iz,iy,ix);
	  
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

      const size_t nTetrahedra = d_tetrahedra.GetDimensions()[0];

      const unsigned int nThreadsx = GetBlockSize( d_output.ElementCount(), nTetrahedra );
      const unsigned int nThreadsy = GetBlockSize( d_output.ElementCount(), nTetrahedra );
      const unsigned int nThreadsz = 1;

      dim3 grid, threads;

      
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

    template<>
    void RunVisitCounterSimpleCUDA<float,double>( CudaImage<int,3,unsigned short>& d_output,
						  const CudaImage<float,3,size_t>& d_tetrahedra ) {
      if( d_tetrahedra.GetDimensions()[1] != nVertices ) {
	throw std::runtime_error("Must have four vertices per tetrahedron!");
      }
      if( d_tetrahedra.GetDimensions()[2] != nDims ) {
	throw std::runtime_error("Only implemented for 3D space");
      }

      SimpleVisitCounter<float,double>( d_output, d_tetrahedra );
    }
  }
}
