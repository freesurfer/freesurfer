#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"

#include <itkMath.h>
#include "vnl/vnl_matrix_fixed.h"
#include "kvlTetrahedronInteriorConstIterator.h"



namespace kvl
{

//
//
//
AtlasMeshToIntensityImageCostAndGradientCalculator
::AtlasMeshToIntensityImageCostAndGradientCalculator()
{

  m_LikelihoodFilter = LikelihoodFilterType::New();  
  
}


//
//
//
AtlasMeshToIntensityImageCostAndGradientCalculator
::~AtlasMeshToIntensityImageCostAndGradientCalculator()
{
}



//
//
//
void 
AtlasMeshToIntensityImageCostAndGradientCalculator
::SetImages( const std::vector< ImageType::ConstPointer >& images )
{
  // 
  for ( int contrastNumber = 0; contrastNumber < images.size(); contrastNumber++ )
    {
    m_LikelihoodFilter->SetInput( contrastNumber, images[ contrastNumber ] );
    }

}


//
//
//
void 
AtlasMeshToIntensityImageCostAndGradientCalculator
::SetParameters( const std::vector< vnl_vector< double > >& means,
                 const std::vector< vnl_matrix< double > >& variances,
                 const std::vector< double >&  mixtureWeights,
                 const std::vector< int >&  numberOfGaussiansPerClass )
{
  m_LikelihoodFilter->SetParameters( means, variances, mixtureWeights, numberOfGaussiansPerClass );
}


//
//
//
void
AtlasMeshToIntensityImageCostAndGradientCalculator
::Rasterize( const AtlasMesh* mesh )
{

  // Make sure the likelihoods are up-to-date
  //m_LikelihoodFilter->SetNumberOfThreads( 1 );
  m_LikelihoodFilter->Update();
  
  // Now rasterize
  Superclass::Rasterize( mesh );
  
}    
  
    
  
 

//
//
//
void 
AtlasMeshToIntensityImageCostAndGradientCalculator
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
  
  // Loop over all voxels within the tetrahedron and do The Right Thing  
  const int  numberOfClasses = alphasInVertex0.Size();
  TetrahedronInteriorConstIterator< LikelihoodFilterType::OutputPixelType >  it( m_LikelihoodFilter->GetOutput(), p0, p1, p2, p3 );
  for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
    it.AddExtraLoading( alphasInVertex0[ classNumber ], 
                        alphasInVertex1[ classNumber ], 
                        alphasInVertex2[ classNumber ], 
                        alphasInVertex3[ classNumber ] );
    }  
    
  for ( ; !it.IsAtEnd(); ++it )
    {
    // Skip voxels for which nothing is known
    if ( it.Value().Size() == 0 )
      {
      //std::cout << "Skipping: " << it.Value().Size() << std::endl;
      continue;
      }
      
    //
    double likelihood = 0.0;
    double  xGradientBasis = 0.0;
    double  yGradientBasis = 0.0;
    double  zGradientBasis = 0.0;
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      // Get the Gaussian mixture model likelihood of this class at the intensity of this pixel
      const double mixture = it.Value()[ classNumber ];
        
      // Add contribution of the likelihood
      likelihood += mixture * it.GetExtraLoadingInterpolatedValue( classNumber );
      
      //
      xGradientBasis += mixture * it.GetExtraLoadingNextRowAddition( classNumber );
      yGradientBasis += mixture * it.GetExtraLoadingNextColumnAddition( classNumber );
      zGradientBasis += mixture * it.GetExtraLoadingNextSliceAddition( classNumber );
      } // End loop over all classes
      
      
    //  Add contribution to log-likelihood
    likelihood = likelihood + 1e-15; //dont want to divide by zero
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
