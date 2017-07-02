#pragma once

#include "kvlAtlasMesh.h"

#include "cudaimage.hpp"

namespace kvl {
  namespace cuda {
    template<typename CoordinateType, typename MeshIndexType>
    class CudaTetrahedralMesh {
    public:
      const unsigned char nDims = 3;
      const unsigned char nVertices = 4;

      void Send( kvl::AtlasMesh::ConstPointer mesh ) {
	auto tetIds = this->GetTetrahedronIds( mesh );

	// Sanity check things
	std::set<size_t> activePointIndices;
	for( size_t iTet=0; iTet<tetIds.size(); iTet++ ) {
	  kvl::AtlasMesh::CellAutoPointer cell;
	  mesh->GetCell( tetIds.at(iTet), cell );
	  
	  size_t pointCount = 0;
	  for( auto pit = cell->PointIdsBegin();
	       pit != cell->PointIdsEnd();
	       ++pit ) {
	    activePointIndices.insert(*pit);
	    pointCount++;
	  }
    
	  if( pointCount != nVertices ) {
	    throw std::runtime_error("Found tetrahedron with wrong vertex count");
	  }
	}

	if( activePointIndices.size() != mesh->GetPoints()->size() ) {
	  throw std::runtime_error("Point index remapping not supported!");
	}

	
      }

    private:
      CudaImage<CoordinateType,2,MeshIndexType> d_vertices;
      CudaImage<MeshIndexType,2,MeshIndexType> d_vertexMap;

      std::vector<kvl::AtlasMesh::CellIdentifier> GetTetrahedronIds( kvl::AtlasMesh::ConstPointer mesh ) const {
	std::vector<kvl::AtlasMesh::CellIdentifier> ids;
	for( auto cellIt = mesh->GetCells()->Begin();
	     cellIt != mesh->GetCells()->End();
	     ++cellIt ) {
	  if( cellIt.Value()->GetType() == kvl::AtlasMesh::CellType::TETRAHEDRON_CELL ) {
	    ids.push_back( cellIt.Index() );
	  }
	}
	return ids;
      }

      void SendVertices( kvl::AtlasMesh::ConstPointer mesh ) {
	Dimension<2,MeshIndexType> vertexDims;
	vertexDims[0] = mesh->GetPointData()->size();
	vertexDims[1] = nDims;

	std::vector<CoordinateType> vertices;
	vertices.resize(vertexDims.GetElementCount());

	for( auto pointIt = mesh->GetPoints()->Begin();
	     pointIt != mesh->GetPoints()->End();
	     ++pointIt ) {
	  for( unsigned int i=0; i<nDims; i++ ) {
	    vertices.at(vertexDims.GetLinearIndex(pointIt.Index(),i)) = pointIt.Value()[i];
	  }
	}

	this->d_vertices.Send( vertices, vertexDims );
      }
    };
  }
}
