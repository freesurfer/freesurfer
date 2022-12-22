#include "kvlAtlasMeshMultiAlphaDrawer.h"

#include "kvlTetrahedronInteriorIterator.h"


namespace kvl
{

//
//
//
AtlasMeshMultiAlphaDrawer
::AtlasMeshMultiAlphaDrawer()
{
  m_Image = 0; 
}  



//
//
//
AtlasMeshMultiAlphaDrawer
::~AtlasMeshMultiAlphaDrawer()
{

  
}  



//
//
//
void
AtlasMeshMultiAlphaDrawer
::Rasterize( const AtlasMesh* mesh )
{
  // Fill image with empty result (needed because area outside of the mesh will never be visited)
  const int  numberOfClasses = mesh->GetPointData()->Begin().Value().m_Alphas.Size();
  AtlasAlphasType  emptyEntry( numberOfClasses );
  emptyEntry.Fill( 0.0f );
  m_Image->FillBuffer( emptyEntry );
  
  //
  Superclass::Rasterize( mesh );
  
}  

//
//
//
bool
AtlasMeshMultiAlphaDrawer
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
  
  // Loop over all voxels within the tetrahedron and do The Right Thing  
  TetrahedronInteriorIterator< ImageType::PixelType >  it( m_Image, p0, p1, p2, p3 );
  const AtlasAlphasType&  alphasInVertex0 = mesh->GetPointData()->ElementAt( id0 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex1 = mesh->GetPointData()->ElementAt( id1 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex2 = mesh->GetPointData()->ElementAt( id2 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex3 = mesh->GetPointData()->ElementAt( id3 ).m_Alphas;
  const int  numberOfClasses = alphasInVertex0.Size();
  for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
    if (alphasInVertex0[ classNumber ] != 0 || alphasInVertex1[ classNumber ] != 0 || alphasInVertex2[ classNumber ] != 0 || alphasInVertex3[ classNumber ] != 0) 
      it.AddExtraLoading( alphasInVertex0[ classNumber ], 
                          alphasInVertex1[ classNumber ], 
                          alphasInVertex2[ classNumber ], 
                          alphasInVertex3[ classNumber ] );
    }

  for ( ; !it.IsAtEnd(); ++it )
    {
      it.Value() = AtlasAlphasType( numberOfClasses );
      int classIdx = 0;
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      if (alphasInVertex0[ classNumber ] != 0 || alphasInVertex1[ classNumber ] != 0 || alphasInVertex2[ classNumber ] != 0 || alphasInVertex3[ classNumber ] != 0)
	{
        it.Value()[ classNumber ] = it.GetExtraLoadingInterpolatedValue( classIdx );
        classIdx++;
	}
      else
        it.Value()[ classNumber ] = 0;
      }
    }
    
  return true;
}


  
} // End namespace kvl
