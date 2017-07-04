#pragma once

#include "kvlAtlasMesh.h"

#include "cudaimage.hpp"

namespace kvl {
  namespace cuda {
    template<typename CoordinateType, typename MeshIndexType>
    class TetrahedralMesh_GPU {
    public:
      Image_GPU<CoordinateType,2,MeshIndexType> vertices;
      Image_GPU<MeshIndexType,2,MeshIndexType> vertexMap;

      __device__
      CoordinateType GetVertexCoordinate( MeshIndexType iTet, MeshIndexType iVert, unsigned char iDim ) const {
	MeshIndexType iVertexId = this->vertexMap(iTet,iVert);

	return this->vertices(iVertexId,iDim);
      }

      __device__
      MeshIndexType GetTetrahedraCount() const {
	return this->vertexMap.dims[0];
      }
    };

    template<typename CoordinateType, typename MeshIndexType>
    class CudaTetrahedralMesh {
    public:
      const unsigned char nDims = 3;
      const unsigned char nVertices = 4;

      typedef TetrahedralMesh_GPU<CoordinateType,MeshIndexType> GPUType;
      typedef CoordinateType CoordType;

      void Send( kvl::AtlasMesh::ConstPointer mesh ) {
	auto tetIds = this->GetTetrahedronIds( mesh );

	// Sanity check size of MeshIndexType
	if( tetIds.size() > std::numeric_limits<MeshIndexType>::max() ) {
	  throw std::out_of_range("Too many tetrahedra for MeshIndexType");
	}

	// Sanity check the points used to define the tetrahedra
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

	// Transfer data to the GPU
	this->SendVertices(mesh);
	this->SendVertexMap(mesh, tetIds);
      }

      TetrahedralMesh_GPU<CoordinateType,MeshIndexType> getArg() const {
	TetrahedralMesh_GPU<CoordinateType,MeshIndexType> gpuArg;

	gpuArg.vertices = this->d_vertices.getArg();
	gpuArg.vertexMap = this->d_vertexMap.getArg();

	return gpuArg;
      }

      MeshIndexType GetTetrahedraCount() const {
	return this->d_vertexMap.GetDimensions()[0];
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
	vertexDims[0] = mesh->GetPoints()->size();
	vertexDims[1] = nDims;

	std::vector<CoordinateType> vertices;
	vertices.resize(vertexDims.ElementCount());

	for( auto pointIt = mesh->GetPoints()->Begin();
	     pointIt != mesh->GetPoints()->End();
	     ++pointIt ) {
	  for( MeshIndexType i=0; i<nDims; i++ ) {
	    size_t idx = vertexDims.GetLinearIndex(static_cast<MeshIndexType>(pointIt.Index()),i);
	    vertices.at(idx) = pointIt.Value()[i];
	  }
	}

	this->d_vertices.Send( vertices, vertexDims );
      }

      void SendVertexMap( kvl::AtlasMesh::ConstPointer mesh, const std::vector<kvl::AtlasMesh::CellIdentifier>& ids ) {
	Dimension<2,MeshIndexType> mapDims;
	mapDims[0] = ids.size();
	mapDims[1] = nVertices;

	std::vector<MeshIndexType> vertexMap;
	vertexMap.resize(mapDims.ElementCount());

	for( MeshIndexType iTet=0; iTet<ids.size(); iTet++ ) {
	  AtlasMesh::CellAutoPointer cell;
	  mesh->GetCell( ids.at(iTet), cell );

	  MeshIndexType pIdx = 0;
	  for( auto pit = cell->PointIdsBegin();
	       pit != cell->PointIdsEnd();
	       ++pit ) {
	     size_t idx = mapDims.GetLinearIndex(iTet,pIdx);
	     vertexMap.at(idx) = *pit;
	     pIdx++;
	  }
	}

	this->d_vertexMap.Send( vertexMap, mapDims );
      }
    };
  }
}
