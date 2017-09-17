# pragma once

#include "simplesharedtetrahedron.hpp"

template<typename AlphasType>
class AlphaDrawerAction {
public:
  AlphaDrawerAction(kvl::cuda::Image_GPU<AlphasType,3,unsigned short> target,
		    const unsigned int iA) : output(target), iAlpha(iA) {} 

  template<typename MeshSupplier, typename T, typename Internal, typename IndexType>
  __device__
  void operator()( const SimpleSharedTetrahedron<MeshSupplier,T,Internal>& tet,
		   const IndexType iz,
		   const IndexType iy,
		   const IndexType ix ) {
    const int nVertices = 4;

    AlphasType alphas[nVertices];

    tet.LoadAlphas( &(alphas[0]), this->iAlpha );
    AlphasType result = tet.BarycentricInterpolation( &(alphas[0]), iz, iy, ix );
    this->output(iz,iy,ix) = result;
  }

private:
  kvl::cuda::Image_GPU<AlphasType,3,unsigned short> output;
  const unsigned int iAlpha;
};
