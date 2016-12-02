#include "kvlAtlasMeshAlphaDrawerGPU.h"

#include "kvlTetrahedronInteriorIterator.h"


namespace kvl
{

//
//
//
AtlasMeshAlphaDrawerGPU
::AtlasMeshAlphaDrawerGPU()
{
  m_LabelNumber = 0;
  m_Image = 0; 
}  



//
//
//
AtlasMeshAlphaDrawerGPU
::~AtlasMeshAlphaDrawerGPU()
{

  
}  

 
//
//
//
bool
AtlasMeshAlphaDrawerGPU
::RasterizeTetrahedron( const AtlasMesh* mesh, 
                        AtlasMesh::CellIdentifier tetrahedronId,
                        int threadNumber )
{
  // Retrieve everything we need to know 
  AtlasMesh::CellAutoPointer  cell;
  mesh->GetCell( tetrahedronId, cell );

  AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
  const AtlasMesh::PointIdentifier  id0 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id1 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id2 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id3 = *pit;
  
  AtlasMesh::PointType p0;
  AtlasMesh::PointType p1;
  AtlasMesh::PointType p2;
  AtlasMesh::PointType p3;
  mesh->GetPoint( id0, &p0 );
  mesh->GetPoint( id1, &p1 );
  mesh->GetPoint( id2, &p2 );
  mesh->GetPoint( id3, &p3 );
  
  const float alphaInVertex0 = ( mesh->GetPointData()->ElementAt( id0 ).m_Alphas )[ m_LabelNumber ];
  const float alphaInVertex1 = ( mesh->GetPointData()->ElementAt( id1 ).m_Alphas )[ m_LabelNumber ];
  const float alphaInVertex2 = ( mesh->GetPointData()->ElementAt( id2 ).m_Alphas )[ m_LabelNumber ];
  const float alphaInVertex3 = ( mesh->GetPointData()->ElementAt( id3 ).m_Alphas )[ m_LabelNumber ];

  
  // Loop over all voxels within the tetrahedron and do The Right Thing  
  TetrahedronInteriorIterator< ImageType::PixelType >  it( m_Image, p0, p1, p2, p3 );
  it.AddExtraLoading( alphaInVertex0, alphaInVertex1, alphaInVertex2, alphaInVertex3 );
  for ( ; !it.IsAtEnd(); ++it )
    {
    it.Value() = it.GetExtraLoadingInterpolatedValue( 0 );
    }
    
  return true;
}


  
} // End namespace kvl