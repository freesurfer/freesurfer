#include "visitcountertetrahedralmeshcudaimpl.hpp"

#include "cudautils.hpp"

#include "simplesharedtetrahedron.hpp"

const unsigned int nDims = 3;
const unsigned int nVertices = 4;

class VisitCounterAction {
public:
  VisitCounterAction(kvl::cuda::Image_GPU<int,3,unsigned short> target) : output(target) {} 

  template<typename T, typename Internal, typename IndexType>
  __device__
  void operator()( const SimpleSharedTetrahedron<T,Internal>& tet,
		   const IndexType iz,
		   const IndexType iy,
		   const IndexType ix ) {
    atomicAdd(&(this->output(iz,iy,ix)),1);
  }

private:
  kvl::cuda::Image_GPU<int,3,unsigned short> output;
};

// Largely copied from visitcountersimplecudaimpl.cu
// TODO Refactor common code
template<typename T,typename Internal,typename MeshSupplier, typename VertexAction>
__global__
void TetrahedronInteriorKernel( kvl::cuda::Image_GPU<int,3,unsigned short> output,
				const MeshSupplier mesh,
				VertexAction action ) {
  const size_t iTet = blockIdx.x + (gridDim.x * blockIdx.y);
  
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
	  bool inside = tet.PointInside(ix,iy,iz);
	  
	  if( inside ) {
	    action(tet, iz, iy, ix );
	  }
	}
      }
    }
  }
}



namespace kvl {
  namespace cuda {
    
    void RunVisitCounterTetrahedralMeshCUDA( CudaImage<int,3,unsigned short>& d_output,
					     const CudaTetrahedralMesh<double,unsigned long>& ctm ) {
      const unsigned int nBlockx = 1024;

      const size_t nTetrahedra = ctm.GetTetrahedraCount();

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
      
      VisitCounterAction vca(d_output.getArg());

      // Run the kernel
      auto err = cudaGetLastError();
      if( cudaSuccess != err ) {
	throw CUDAException(err);
      }
      TetrahedronInteriorKernel<double,double,kvl::cuda::TetrahedralMesh_GPU<double,unsigned long>,VisitCounterAction><<<grid,threads>>>( d_output.getArg(), ctm.getArg(), vca );
      err = cudaDeviceSynchronize();
      if( cudaSuccess != err ) {
	throw CUDAException(err);
      }
    } 
  }
}
