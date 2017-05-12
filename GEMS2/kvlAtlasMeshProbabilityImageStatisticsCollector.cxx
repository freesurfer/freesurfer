#include "kvlAtlasMeshProbabilityImageStatisticsCollector.h"

#include "kvlTetrahedronInteriorConstIterator.h"



namespace kvl
{

//
//
//
AtlasMeshProbabilityImageStatisticsCollector
::AtlasMeshProbabilityImageStatisticsCollector()
{

  m_ProbabilityImage = 0;

}


//
//
//
AtlasMeshProbabilityImageStatisticsCollector
::~AtlasMeshProbabilityImageStatisticsCollector()
{
}




//
//
//
void
AtlasMeshProbabilityImageStatisticsCollector
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
  const int  numberOfClasses = alphasInVertex0.Size();
  statisticsInVertex0 = AtlasAlphasType( numberOfClasses );
  statisticsInVertex0.Fill( 0.0f );
  statisticsInVertex1 = AtlasAlphasType( numberOfClasses );
  statisticsInVertex1.Fill( 0.0f );
  statisticsInVertex2 = AtlasAlphasType( numberOfClasses );
  statisticsInVertex2.Fill( 0.0f );
  statisticsInVertex3 = AtlasAlphasType( numberOfClasses );
  statisticsInVertex3.Fill( 0.0f );

  
  // Loop over all voxels within the tetrahedron and do The Right Thing  
  TetrahedronInteriorConstIterator< ProbabilityImageType::PixelType >  it( m_ProbabilityImage, p0, p1, p2, p3 );
  for ( ; !it.IsAtEnd(); ++it )
    {
      
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      // Get the weight of this class's contribution
      const double  weight = it.Value()[ classNumber ];

      //
      double  W0 = alphasInVertex0[ classNumber ] * it.GetPi0();
      double  W1 = alphasInVertex1[ classNumber ] * it.GetPi1();
      double  W2 = alphasInVertex2[ classNumber ] * it.GetPi2();
      double  W3 = alphasInVertex3[ classNumber ] * it.GetPi3();
  
      //
      const double  denominator = ( W0 + W1 + W2 + W3 + 1e-15 );
      minLogLikelihood -= weight * log( denominator );
 
      // Normalize to obtain W0, W1, W2, and W3
      W0 /= denominator;
      W1 /= denominator;
      W2 /= denominator;
      W3 /= denominator;
      
      // Update the histogram entries in the vertices accordingly
      statisticsInVertex0[ classNumber ] += W0 * weight;
      statisticsInVertex1[ classNumber ] += W1 * weight;
      statisticsInVertex2[ classNumber ] += W2 * weight;
      statisticsInVertex3[ classNumber ] += W3 * weight;
      } // End loop over all classes
  
    } // End loop over all pixels within tetrahedron  

}






} // end namespace kvl
