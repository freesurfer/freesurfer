#include "kvlAtlasMeshSummaryDrawer.h"

#include "kvlTetrahedronInteriorIterator.h"


namespace kvl
{

//
//
//
AtlasMeshSummaryDrawer
::AtlasMeshSummaryDrawer()
{
  m_Image = 0; 
  m_CompressionLookupTable = 0;
}  



//
//
//
AtlasMeshSummaryDrawer
::~AtlasMeshSummaryDrawer()
{

  
}  

 
//
//
//
bool
AtlasMeshSummaryDrawer
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

  const AtlasAlphasType&  alphasInVertex0 = mesh->GetPointData()->ElementAt( id0 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex1 = mesh->GetPointData()->ElementAt( id1 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex2 = mesh->GetPointData()->ElementAt( id2 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex3 = mesh->GetPointData()->ElementAt( id3 ).m_Alphas;


  // Let's compute the resulting RGBA colors in each of the vertices
  std::vector< double >  colorInVertex0( 4, 0.0 );
  std::vector< double >  colorInVertex1( 4, 0.0 );
  std::vector< double >  colorInVertex2( 4, 0.0 );
  std::vector< double >  colorInVertex3( 4, 0.0 );
  const int  numberOfClasses = m_CompressionLookupTable->GetNumberOfClasses();
  for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
    const double  r = m_CompressionLookupTable->GetColor( classNumber )[ 0 ];
    const double  g = m_CompressionLookupTable->GetColor( classNumber )[ 1 ];
    const double  b = m_CompressionLookupTable->GetColor( classNumber )[ 2 ];
    const double  a = m_CompressionLookupTable->GetColor( classNumber )[ 3 ];
    //std::cout << r << " " << g << " " << b << " " << a << std::endl;
      
    const double  weight0 = alphasInVertex0[ classNumber ];  
    colorInVertex0[ 0 ] += ( a / 255 ) * r * weight0;
    colorInVertex0[ 1 ] += ( a / 255 ) * g * weight0;
    colorInVertex0[ 2 ] += ( a / 255 ) * b * weight0;
    colorInVertex0[ 3 ] += a * weight0;
    
    const double  weight1 = alphasInVertex1[ classNumber ];  
    colorInVertex1[ 0 ] += ( a / 255 ) * r * weight1;
    colorInVertex1[ 1 ] += ( a / 255 ) * g * weight1;
    colorInVertex1[ 2 ] += ( a / 255 ) * b * weight1;
    colorInVertex1[ 3 ] += a * weight1;

    const double  weight2 = alphasInVertex2[ classNumber ];  
    colorInVertex2[ 0 ] += ( a / 255 ) * r * weight2;
    colorInVertex2[ 1 ] += ( a / 255 ) * g * weight2;
    colorInVertex2[ 2 ] += ( a / 255 ) * b * weight2;
    colorInVertex2[ 3 ] += a * weight2;

    const double  weight3 = alphasInVertex3[ classNumber ];  
    colorInVertex3[ 0 ] += ( a / 255 ) * r * weight3;
    colorInVertex3[ 1 ] += ( a / 255 ) * g * weight3;
    colorInVertex3[ 2 ] += ( a / 255 ) * b * weight3;
    colorInVertex3[ 3 ] += a * weight3;
    }  
  
  
  // Loop over all voxels within the tetrahedron and do The Right Thing  
  TetrahedronInteriorIterator< ImageType::PixelType >  it( m_Image, p0, p1, p2, p3 );
  for ( int colorNumber = 0; colorNumber < 4; colorNumber++ )
    {
    it.AddExtraLoading( colorInVertex0[ colorNumber ], 
                        colorInVertex1[ colorNumber ], 
                        colorInVertex2[ colorNumber ], 
                        colorInVertex3[ colorNumber ] ); 
    }
  for ( ; !it.IsAtEnd(); ++it )
    {
    for ( int colorNumber = 0; colorNumber < 4; colorNumber++ )
      {
      //std::cout << it.GetExtraLoadingInterpolatedValue( colorNumber ) << " " << std::flush;
      it.Value()[ colorNumber ] = 
             static_cast< unsigned char >( it.GetExtraLoadingInterpolatedValue( colorNumber ) + 0.5 );
      }
    //std::cout << std::endl;  
  
    }
    
  return true;
}


  
} // End namespace kvl
