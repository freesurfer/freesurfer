#include "kvlMutualInformationCostAndGradientCalculator.h"

#include <itkMath.h>
#include "kvlTetrahedronInteriorConstIterator.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"


namespace kvl
{

//
//
//
MutualInformationCostAndGradientCalculator
::MutualInformationCostAndGradientCalculator()
{
  
  m_Histogrammer = Histogrammer::New();
  
}


//
//
//
MutualInformationCostAndGradientCalculator
::~MutualInformationCostAndGradientCalculator()
{
}




//
//
//
void 
MutualInformationCostAndGradientCalculator
::SetImage( const ImageType* image )
{
  m_Histogrammer->SetImage( image );
  
}
 


//
//
//
void
MutualInformationCostAndGradientCalculator
::Rasterize( const AtlasMesh* mesh )
{

  // Estimate parameters of generative model using EM
  const int  numberOfBins = 64;
  const int  numberOfClasses = mesh->GetPointData()->Begin().Value().m_Alphas.Size();
  Histogrammer::ConditionalIntensityDistributionType  
            uniformDistribution( numberOfBins, 
                                 1.0 / static_cast< double >( numberOfBins ) );
  std::vector< Histogrammer::ConditionalIntensityDistributionType >
            conditionalIntensityDistributions( numberOfClasses, uniformDistribution );
  const int  maximumNumberOfIterations = 10;   
  double  minLogLikelihood = itk::NumericTraits< double >::max();
  for ( int iterationNumber = 0; iterationNumber < maximumNumberOfIterations; iterationNumber++ )
    {
    // E-step
    m_Histogrammer->SetConditionalIntensityDistributions( conditionalIntensityDistributions );          
    m_Histogrammer->Rasterize( mesh );
  

    // M-step
    const Histogrammer::HistogramType&  histogram = m_Histogrammer->GetHistogram();
    double  numberOfVoxels = 0.0;
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      double  numberOfVoxelsInThisClass = 1e-15;
      for ( int binNumber = 0; binNumber < numberOfBins; binNumber++ )
        {
        numberOfVoxelsInThisClass += histogram[ classNumber ][ binNumber ];  
        }  
      for ( int binNumber = 0; binNumber < numberOfBins; binNumber++ )
        {
        conditionalIntensityDistributions[ classNumber ][ binNumber ] =  
              histogram[ classNumber ][ binNumber ] / numberOfVoxelsInThisClass;  
        }  
       
      numberOfVoxels += numberOfVoxelsInThisClass;
      } // End loop over classes
      
    // Check convergence
    const double  previousMinLogLikelihood = minLogLikelihood;
    minLogLikelihood = m_Histogrammer->GetMinLogLikelihood();
    std::cout << "minLogLikelihood: " << minLogLikelihood << std::endl;
    std::cout << "numberOfVoxels: " << numberOfVoxels << std::endl;
    const double  changeInCostPerVoxel = ( previousMinLogLikelihood - minLogLikelihood ) 
                                         / static_cast< double >( numberOfVoxels );
    std::cout << "changeInCostPerVoxel: " << changeInCostPerVoxel << std::endl;                                     
    if ( changeInCostPerVoxel < 1e-3 )
      {
      break;
      }  
      
    } // End loop over EM iterations  

    
    
  //
  if ( 0 )
    {
    const Histogrammer::HistogramType&  histogram = m_Histogrammer->GetHistogram();
    std::cout << "histogram = [ ..." << std::endl;
    for ( int binNumber = 0; binNumber < numberOfBins; binNumber++ )
      {
        
      for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
        {
        std::cout << histogram[ classNumber ][ binNumber ] << " ";  
        }
      if ( binNumber == ( numberOfBins-1 ) )
        {
        std::cout << "]" << std::endl;
        }
      else
        {
        std::cout << "; ..." << std::endl;
        }
      }
    }  

    
  // Now rasterize
  Superclass::Rasterize( mesh );
  
}    
  
    
 

//
//
//
void 
MutualInformationCostAndGradientCalculator
::AddDataContributionOfTetrahedron( const AtlasMesh::PointType& p0,
                                    const AtlasMesh::PointType& p1,
                                    const AtlasMesh::PointType& p2,
                                    const AtlasMesh::PointType& p3,
                                    const AtlasAlphasType&  alphasInVertex0,
                                    const AtlasAlphasType&  alphasInVertex1,
                                    const AtlasAlphasType&  alphasInVertex2,
                                    const AtlasAlphasType&  alphasInVertex3,
                                    double&  priorPlusDataCost,
                                    AtlasPositionGradientType&  gradientInVertex0,
                                    AtlasPositionGradientType&  gradientInVertex1,
                                    AtlasPositionGradientType&  gradientInVertex2,
                                    AtlasPositionGradientType&  gradientInVertex3 )
{
#if 0  
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
      // Get the Gaussian likelihood of this class at the intensity of this pixel
      const double gauss = it.Value()[ classNumber ];
        
      // Add contribution of the likelihood
      likelihood += gauss * it.GetExtraLoadingInterpolatedValue( classNumber );
      
      //
      xGradientBasis += gauss * it.GetExtraLoadingNextRowAddition( classNumber );
      yGradientBasis += gauss * it.GetExtraLoadingNextColumnAddition( classNumber );
      zGradientBasis += gauss * it.GetExtraLoadingNextSliceAddition( classNumber );
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

#endif

}



} // end namespace kvl
