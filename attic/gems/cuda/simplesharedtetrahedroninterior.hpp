#pragma once

#include "simplesharedtetrahedron.hpp"
#include "cudautils.hpp"

template<typename MeshSupplier, typename VertexAction, typename IndexType, typename T, typename Internal>
__global__
void SimpleSharedTetrahedronInteriorKernel( const IndexType nz,
					    const IndexType ny,
					    const IndexType nx,
					    const MeshSupplier mesh,
					    VertexAction action ) {
  const unsigned int nDims = 3;
  const unsigned int nVertices = 4;
  const size_t iTet = blockIdx.x + (gridDim.x * blockIdx.y);
  
  // Check if this block has an assigned tetrahedron
  if( iTet >= mesh.GetTetrahedraCount() ) {
    return;
  }

  // Load the tetrahedron and determine bounding box
  __shared__ T tetrahedron[nVertices][nDims];
  __shared__ IndexType min[nDims], max[nDims];
  __shared__ T M[nDims][nDims];
  SimpleSharedTetrahedron<MeshSupplier,T,Internal> tet(mesh, tetrahedron, M);

  tet.LoadAndBoundingBox( iTet, min, max );

  tet.ComputeBarycentricTransform();

  // Figure out how to cover the bounding box with the current thread block
  // We assume that each thread block is strictly 2D

  // Divide the bounding box into blocks equal to the blockDim
  for( IndexType iyStart=min[1]; iyStart<max[1]; iyStart += blockDim.y ) {
    for( IndexType ixStart=min[0]; ixStart<max[0]; ixStart += blockDim.x ) {
      const IndexType ix = ixStart + threadIdx.x;
      const IndexType iy = iyStart + threadIdx.y;

      // Could probably do this test a little better
      if( (iy<ny) && (ix<nx) ) {

	for( IndexType iz=min[2]; iz<max[2]; iz++ ) {
	  if( iz<nz ) {
	    bool inside = tet.PointInside(iz,iy,ix);
	    
	    if( inside ) {
	      action(tet, iz, iy, ix );
	    }
	  }
	}
      }
    }
  }
}

namespace kvl {
  namespace cuda {
    template<typename MeshSupplier, typename VertexAction, typename IndexType, typename Internal = double>
    void RunSimpleSharedTetrahedron( const IndexType nz,
				     const IndexType ny,
				     const IndexType nx,
				     const MeshSupplier& mesh,
				     VertexAction& va ) {
      const unsigned int nBlockx = 1024;
      
      const size_t nTetrahedra = mesh.GetTetrahedraCount();
      
      const size_t volume = static_cast<size_t>(nz) * static_cast<size_t>(ny) * static_cast<size_t>(nx);
      const unsigned int nThreadsx = GetBlockSize( volume, nTetrahedra );
      const unsigned int nThreadsy = GetBlockSize( volume, nTetrahedra );
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
      
      typedef typename MeshSupplier::GPUType MeshArg;
      typedef typename MeshSupplier::CoordType MeshCoordType;

      SimpleSharedTetrahedronInteriorKernel<MeshArg,VertexAction,IndexType,MeshCoordType,Internal><<<grid,threads>>>(nz, ny, nx, mesh.getArg(), va );
      err = cudaDeviceSynchronize();
      if( cudaSuccess != err ) {
	throw CUDAException(err);
      }
    }
  }
}
