#include "kvlAtlasMeshLabelImageStatisticsCollector.h"

#include "kvlTetrahedronInteriorConstIterator.h"



namespace kvl
{

//
//
//
AtlasMeshLabelImageStatisticsCollector
::AtlasMeshLabelImageStatisticsCollector()
{

  m_LabelImage = 0;
  m_CompressionLookupTable = 0;

}


//
//
//
AtlasMeshLabelImageStatisticsCollector
::~AtlasMeshLabelImageStatisticsCollector()
{
}



//
//
//
void 
AtlasMeshLabelImageStatisticsCollector
::SetLabelImage( const LabelImageType* labelImage,
                 const CompressionLookupTable* lookupTable )
{
  m_LabelImage = labelImage;
  m_CompressionLookupTable = lookupTable;
}



//
//
//
void
AtlasMeshLabelImageStatisticsCollector
::GetContributionOfTetrahedron( const AtlasMesh::PointType& p0,
                                const AtlasMesh::PointType& p1,
                                const AtlasMesh::PointType& p2,
                                const AtlasMesh::PointType& p3,
                                const AtlasAlphasType&  alphasInVertex0,
                                const AtlasAlphasType&  alphasInVertex1,
                                const AtlasAlphasType&  alphasInVertex2,
                                const AtlasAlphasType&  alphasInVertex3,
                                double&  minLogLikelihood,
                                AtlasAlphasType&  statisticsInVertex0,
                                AtlasAlphasType&  statisticsInVertex1,
                                AtlasAlphasType&  statisticsInVertex2,
                                AtlasAlphasType&  statisticsInVertex3 )
{           
  
  // We start with an empty slate
  minLogLikelihood = 0.0;
  const int  numberOfClasses = m_CompressionLookupTable->GetNumberOfClasses();
  statisticsInVertex0 = AtlasAlphasType( numberOfClasses );
  statisticsInVertex0.Fill( 0.0f );
  statisticsInVertex1 = AtlasAlphasType( numberOfClasses );
  statisticsInVertex1.Fill( 0.0f );
  statisticsInVertex2 = AtlasAlphasType( numberOfClasses );
  statisticsInVertex2.Fill( 0.0f );
  statisticsInVertex3 = AtlasAlphasType( numberOfClasses );
  statisticsInVertex3.Fill( 0.0f );

  
  // Loop over all voxels within the tetrahedron and do The Right Thing  
  TetrahedronInteriorConstIterator< LabelImageType::PixelType >  it( m_LabelImage, p0, p1, p2, p3 );
  for ( ; !it.IsAtEnd(); ++it )
    {
    double  denominator = 0.0;
    const std::vector< int >&  classNumbers = m_CompressionLookupTable->GetClassNumbers( it.Value() );  
    for ( std::vector< int >::const_iterator classIt = classNumbers.begin(); 
          classIt != classNumbers.end();
          ++classIt )
      {
      double  W0 = alphasInVertex0[ *classIt ] * it.GetPi0();
      double  W1 = alphasInVertex1[ *classIt ] * it.GetPi1();
      double  W2 = alphasInVertex2[ *classIt ] * it.GetPi2();
      double  W3 = alphasInVertex3[ *classIt ] * it.GetPi3();
  
      denominator += ( W0 + W1 + W2 + W3 + 1e-15 );
      }
    minLogLikelihood -= log( denominator );
  
    for ( std::vector< int >::const_iterator classIt = classNumbers.begin(); 
          classIt != classNumbers.end();
          ++classIt )
      {
      double  W0 = alphasInVertex0[ *classIt ] * it.GetPi0() / denominator;
      double  W1 = alphasInVertex1[ *classIt ] * it.GetPi1() / denominator;
      double  W2 = alphasInVertex2[ *classIt ] * it.GetPi2() / denominator; 
      double  W3 = alphasInVertex3[ *classIt ] * it.GetPi3() / denominator;
    
      statisticsInVertex0[ *classIt ] += W0;
      statisticsInVertex1[ *classIt ] += W1;
      statisticsInVertex2[ *classIt ] += W2;
      statisticsInVertex3[ *classIt ] += W3;
      }
      
    } // End loop over all pixels within tetrahedron  

}






} // end namespace kvl
