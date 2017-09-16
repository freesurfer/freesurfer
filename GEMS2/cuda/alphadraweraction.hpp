# pragma once

#include "simplesharedtetrahedron.hpp"

template<typename AlphasType>
class AlphaDrawerAction {
  AlphaDrawerAction(kvl::cuda::Image_GPU<AlphasType,3,unsigned short> target) : output(target) {} 

  template<typename MeshSupplier, typename T, typename Internal, typename IndexType>
  __device__
  void operator()( const SimpleSharedTetrahedron<MeshSupplier,T,Internal>& tet,
		   const IndexType iz,
		   const IndexType iy,
		   const IndexType ix ) {
    const int nVertices = 4;

    __shared__ AlphasType alphas[nVertices];

    tet->LoadAlphas( &(alphas[0]), this->iAlpha );
    this->output(iz,iy,ix)) = tet.BarycentricInterpolation( alphas, iz, iy. ix );;
  }

private:
  kvl::cuda::Image_GPU<AlphasType,3,unsigned short> output;
  unsigned short iAlpha;
}
