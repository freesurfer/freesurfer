#include "imageutils.hpp"

#include <itkAutomaticTopologyMeshSource.h>

typedef itk::AutomaticTopologyMeshSource<kvl::AtlasMesh> MeshSource;
typedef MeshSource::IdentifierType  IdentifierType;

namespace kvl {
  namespace Testing {
    kvl::AtlasMesh::Pointer CreateSingleTetrahedronMesh( const float vertices[nVertices][nDims],
							 const unsigned int nAlphas ) {      
      MeshSource::Pointer meshSource = MeshSource::New();
      
      const IdentifierType  id0 = meshSource->AddPoint( vertices[0] );
      const IdentifierType  id1 = meshSource->AddPoint( vertices[1] );
      const IdentifierType  id2 = meshSource->AddPoint( vertices[2] );
      const IdentifierType  id3 = meshSource->AddPoint( vertices[3] );
      meshSource->AddTetrahedron( id0, id1, id2, id3 );
      
      auto mesh = meshSource->GetOutput();
      
      kvl::PointParameters universalParams;
      universalParams.m_Alphas = kvl::AtlasAlphasType(nAlphas);
      for( unsigned int i=0; i<nAlphas; i++ ) {
	universalParams.m_Alphas[i] = i;
      }

      for( unsigned int i=0; i<mesh->GetNumberOfPoints(); i++ ) {
	mesh->SetPointData(i, universalParams );
      }
      
      return mesh;
    }
  }
}
