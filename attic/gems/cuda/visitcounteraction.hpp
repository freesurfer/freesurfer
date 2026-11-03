#pragma once

#include "simplesharedtetrahedron.hpp"

class VisitCounterAction {
public:
  VisitCounterAction(kvl::cuda::Image_GPU<int,3,unsigned short> target) : output(target) {} 

  template<typename MeshSupplier, typename T, typename Internal, typename IndexType>
  __device__
  void operator()( const SimpleSharedTetrahedron<MeshSupplier,T,Internal>& tet,
		   const IndexType iz,
		   const IndexType iy,
		   const IndexType ix ) {
    atomicAdd(&(this->output(iz,iy,ix)),1);
  }

private:
  kvl::cuda::Image_GPU<int,3,unsigned short> output;
};
