#pragma once

#include "kvlAtlasMesh.h"

namespace kvl {
  namespace Testing {

    const int nDims = 3;
    const int nVertices = 4;

    template<typename ImageType>
    typename ImageType::Pointer CreateImageCube( const int sideLength, const int value ) {
      const int nx = sideLength;
      const int ny = sideLength;
      const int nz = sideLength;
      
      typename ImageType::RegionType region;
      typename ImageType::IndexType start;
      typename ImageType::SizeType size;
      start[0] = start[1] = start[2] = 0;
      size[0] = nx;
      size[1] = ny;
      size[2] = nz;
      
      region.SetSize(size);
      region.SetIndex(start);
      
      typename ImageType::Pointer image = ImageType::New();
      image->SetRegions(region);
      image->Allocate();
      
      for( int k=0; k<nz; k++ ) {
	for( int j=0; j<ny; j++ ) {
	  for( int i=0; i<nx; i++ ) {
	    typename ImageType::IndexType idx;
	    idx[0] = i;
	    idx[1] = j;
	    idx[2] = k;
	    
	    image->SetPixel(idx, value);
	  }
	}
      }
      
      return image;
    }
    
    kvl::AtlasMesh::Pointer CreateSingleTetrahedronMesh( const float vertices[nVertices][nDims],
							 const unsigned int nAlphas );
  }
}
