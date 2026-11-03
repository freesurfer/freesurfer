#include "kvlAtlasMeshToLabelImageCostAndGradientCalculator.h"

#include <itkMath.h>
#include "vnl/vnl_matrix_fixed.h"
#include "kvlTetrahedronInteriorConstIterator.h"



namespace kvl
{

//
//
//
AtlasMeshToLabelImageCostAndGradientCalculator
::AtlasMeshToLabelImageCostAndGradientCalculator()
{

  m_LabelImage = 0;
  m_CompressionLookupTable = 0;

}


//
//
//
AtlasMeshToLabelImageCostAndGradientCalculator
::~AtlasMeshToLabelImageCostAndGradientCalculator()
{
}



//
//
//
void 
AtlasMeshToLabelImageCostAndGradientCalculator
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
AtlasMeshToLabelImageCostAndGradientCalculator
::AddDataContributionOfTetrahedron( const AtlasMesh::PointType& p0,
                                    const AtlasMesh::PointType& p1,
                                    const AtlasMesh::PointType& p2,
                                    const AtlasMesh::PointType& p3,
                                    const AtlasAlphasType&  alphasInVertex0,
                                    const AtlasAlphasType&  alphasInVertex1,
                                    const AtlasAlphasType&  alphasInVertex2,
                                    const AtlasAlphasType&  alphasInVertex3,
                                    ThreadAccumDataType&  priorPlusDataCost,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex0,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex1,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex2,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex3 )
{

  // Set up voxel iterator
  TetrahedronInteriorConstIterator< LabelImageType::PixelType >  it( m_LabelImage, p0, p1, p2, p3 );

  // Precompute some intermediate results we're going to use over and over
  const int  numberOfClasses = m_CompressionLookupTable->GetNumberOfClasses();
  std::vector< double >  xGradientBasisHelper( numberOfClasses );
  std::vector< double >  yGradientBasisHelper( numberOfClasses );
  std::vector< double >  zGradientBasisHelper( numberOfClasses );
  for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
    const double  alpha0 = alphasInVertex0[ classNumber ];  
    const double  alpha1 = alphasInVertex1[ classNumber ];  
    const double  alpha2 = alphasInVertex2[ classNumber ];  
    const double  alpha3 = alphasInVertex3[ classNumber ];  
        
    xGradientBasisHelper[ classNumber ] =
                        alpha0 * it.GetExtraLoadingNextRowAddition( -4 ) +
                        alpha1 * it.GetExtraLoadingNextRowAddition( -3 ) +
                        alpha2 * it.GetExtraLoadingNextRowAddition( -2 ) +
                        alpha3 * it.GetExtraLoadingNextRowAddition( -1 );
      
    yGradientBasisHelper[ classNumber ] = 
                        alpha0 * it.GetExtraLoadingNextColumnAddition( -4 ) +
                        alpha1 * it.GetExtraLoadingNextColumnAddition( -3 ) +
                        alpha2 * it.GetExtraLoadingNextColumnAddition( -2 ) +
                        alpha3 * it.GetExtraLoadingNextColumnAddition( -1 );
          
    zGradientBasisHelper[ classNumber ] = 
                        alpha0 * it.GetExtraLoadingNextSliceAddition( -4 ) +
                        alpha1 * it.GetExtraLoadingNextSliceAddition( -3 ) +
                        alpha2 * it.GetExtraLoadingNextSliceAddition( -2 ) +
                        alpha3 * it.GetExtraLoadingNextSliceAddition( -1 );
    }
  
  
  // Loop over all voxels within the tetrahedron and do The Right Thing  
  for ( ; !it.IsAtEnd(); ++it )
    {
    double  alpha0 = 0.0;
    double  alpha1 = 0.0;
    double  alpha2 = 0.0;
    double  alpha3 = 0.0;
    double  xGradientBasis = 0.0;
    double  yGradientBasis = 0.0;
    double  zGradientBasis = 0.0;
    const std::vector< int >&  classNumbers = m_CompressionLookupTable->GetClassNumbers( it.Value() );  
    for ( std::vector< int >::const_iterator classIt = classNumbers.begin(); 
          classIt != classNumbers.end();
          ++classIt )
      {
      alpha0 += alphasInVertex0[ *classIt ];
      alpha1 += alphasInVertex1[ *classIt ];
      alpha2 += alphasInVertex2[ *classIt ];
      alpha3 += alphasInVertex3[ *classIt ];
        
      xGradientBasis += xGradientBasisHelper[ *classIt ];
      
      yGradientBasis += yGradientBasisHelper[ *classIt ];
                        
      zGradientBasis += zGradientBasisHelper[ *classIt ];
                        
      } // End loop over all classes associated with this label                  

    //  Add contribution to log-likelihood
    const double  likelihood = alpha0 * it.GetPi0() + 
                               alpha1 * it.GetPi1() +
                               alpha2 * it.GetPi2() +
                               alpha3 * it.GetPi3() + 1e-15; 
    priorPlusDataCost -= log( likelihood );


    //
    xGradientBasis /= likelihood;
    yGradientBasis /= likelihood;
    zGradientBasis /= likelihood;

    // Add contribution to gradient in vertex 0
    gradientInVertex0[ 0 ] += xGradientBasis * it.GetPi0();
    gradientInVertex0[ 1 ] += yGradientBasis * it.GetPi0();
    gradientInVertex0[ 2 ] += zGradientBasis * it.GetPi0();

    // Add contribution to gradient in vertex 1
    gradientInVertex1[ 0 ] += xGradientBasis * it.GetPi1();
    gradientInVertex1[ 1 ] += yGradientBasis * it.GetPi1();
    gradientInVertex1[ 2 ] += zGradientBasis * it.GetPi1();
    
    // Add contribution to gradient in vertex 2
    gradientInVertex2[ 0 ] += xGradientBasis * it.GetPi2();
    gradientInVertex2[ 1 ] += yGradientBasis * it.GetPi2();
    gradientInVertex2[ 2 ] += zGradientBasis * it.GetPi2();
    
    // Add contribution to gradient in vertex 3
    gradientInVertex3[ 0 ] += xGradientBasis * it.GetPi3();
    gradientInVertex3[ 1 ] += yGradientBasis * it.GetPi3();
    gradientInVertex3[ 2 ] += zGradientBasis * it.GetPi3();
    
    
    } // End loop over all voxels within the tetrahedron
   
}



} // end namespace kvl
