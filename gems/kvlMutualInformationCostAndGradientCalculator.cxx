#include "kvlMutualInformationCostAndGradientCalculator.h"

#include <itkMath.h>
#include "kvlTetrahedronInteriorConstIterator.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include <fstream>

namespace kvl
{

//
//
//
MutualInformationCostAndGradientCalculator
::MutualInformationCostAndGradientCalculator()
{
  
  m_Histogrammer = Histogrammer::New();
  m_NumberOfVoxels = 0.0;
  
  //this->SetNumberOfThreads( 1 );
  
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
  m_Histogrammer->SetConditionalIntensityDistributions( conditionalIntensityDistributions );          
  const int  maximumNumberOfIterations = 10; // 10;   
  double  minLogLikelihood = itk::NumericTraits< double >::max();
  for ( int iterationNumber = 0; iterationNumber < maximumNumberOfIterations; iterationNumber++ )
    {
    // E-step
    m_Histogrammer->Rasterize( mesh );
  

    // M-step
    const Histogrammer::HistogramType&  histogram = m_Histogrammer->GetHistogram();
    m_NumberOfVoxels = 0.0;
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
       
      m_NumberOfVoxels += numberOfVoxelsInThisClass;
      } // End loop over classes
    m_Histogrammer->SetConditionalIntensityDistributions( conditionalIntensityDistributions );          

      
    // Check convergence
    const double  previousMinLogLikelihood = minLogLikelihood;
    minLogLikelihood = m_Histogrammer->GetMinLogLikelihood();
    //std::cout << "minLogLikelihood: " << minLogLikelihood << std::endl;
    //std::cout << "m_NumberOfVoxels: " << m_NumberOfVoxels << std::endl;
    const double  changeInCostPerVoxel = ( previousMinLogLikelihood - minLogLikelihood ) 
                                         / static_cast< double >( m_NumberOfVoxels );
    //std::cout << "changeInCostPerVoxel: " << changeInCostPerVoxel << std::endl;                                     
    if ( changeInCostPerVoxel < 1e-3 )
      {
      break;
      }  
      
    } // End loop over EM iterations  

    
    
  //
  if ( 0 )
    {
    const Histogrammer::HistogramType&  histogram = m_Histogrammer->GetHistogram();
  
    std::ofstream  out( "loadHistogram.m" );
    if ( out.bad() )
      {
      std::cerr << "Can't open file loadHistogram.m for writing." << std::endl;
      }

    out << "histogram = [ ..." << std::endl;
    for ( int binNumber = 0; binNumber < numberOfBins; binNumber++ )
      {
        
      for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
        {
        out << histogram[ classNumber ][ binNumber ] << " ";  
        }
      if ( binNumber == ( numberOfBins-1 ) )
        {
        out << "]" << std::endl;
        }
      else
        {
        out << "; ..." << std::endl;
        }
      }
    }  

    
  // Now rasterize to get approximate gradients (but actual data cost is computed separately)
  Superclass::Rasterize( mesh );
  
  
  // Compute Mutual Information, and add it to cost from prior
  if ( m_MinLogLikelihoodTimesPrior == itk::NumericTraits< double >::max() )
    {
    return;
    }
  
  const Histogrammer::HistogramType&  histogram = m_Histogrammer->GetHistogram();
  double  negativeMutualInformation = 0.0;
  std::vector< double >  marginalIntensityDistribution( numberOfBins, 0.0 );  
  for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
    double   marginalProbabilityOfClass = 0.0;
    for ( int binNumber = 0; binNumber < numberOfBins; binNumber++ )
      {
      const double  jointProbability = histogram[ classNumber ][ binNumber ] / m_NumberOfVoxels + 1e-15;  
      negativeMutualInformation -= jointProbability * log( jointProbability ); 
      marginalProbabilityOfClass += jointProbability;
      marginalIntensityDistribution[ binNumber ] += jointProbability;
      }
    negativeMutualInformation += marginalProbabilityOfClass * log( marginalProbabilityOfClass );  
    }
  for ( int binNumber = 0; binNumber < numberOfBins; binNumber++ )
    {
    const double  marginalProbabilityOfIntensity = marginalIntensityDistribution[ binNumber ];
    negativeMutualInformation += marginalProbabilityOfIntensity * log( marginalProbabilityOfIntensity );
    }
  //std::cout << "negativeMutualInformation: " << negativeMutualInformation << std::endl;
  //std::cout << "m_NumberOfVoxels: " << m_NumberOfVoxels << std::endl;
  //std::cout << "priorCost: " << m_MinLogLikelihoodTimesPrior << std::endl;
  //std::cout << "dataCost: " << m_NumberOfVoxels * negativeMutualInformation << std::endl;

  m_MinLogLikelihoodTimesPrior += negativeMutualInformation;
    
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
                                    ThreadAccumDataType&  priorPlusDataCost,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex0,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex1,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex2,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex3 )
{

  // Loop over all voxels within the tetrahedron and do The Right Thing  
  const int  numberOfClasses = alphasInVertex0.Size();
  TetrahedronInteriorConstIterator< Histogrammer::BinnedImageType::PixelType >  
                                                it( m_Histogrammer->GetBinnedImage(), p0, p1, p2, p3 );
  for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
    if (alphasInVertex0[ classNumber ] != 0 || alphasInVertex1[ classNumber ] != 0 || alphasInVertex2[ classNumber ] != 0 || alphasInVertex3[ classNumber ] != 0)
      it.AddExtraLoading( alphasInVertex0[ classNumber ], 
                          alphasInVertex1[ classNumber ], 
                          alphasInVertex2[ classNumber ], 
                          alphasInVertex3[ classNumber ] );
    }  
    
  for ( ; !it.IsAtEnd(); ++it )
    {
    const int  binNumber = it.Value();
      
     // Skip sentinel values
    if ( binNumber < 0 )
      {
      continue;  
      }  
     
      
    //
    double likelihood = 1e-15;
    double  xGradientBasis = 0.0;
    double  yGradientBasis = 0.0;
    double  zGradientBasis = 0.0;
#if 1    
    int classIdx = 0;
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      if (alphasInVertex0[ classNumber ] != 0 || alphasInVertex1[ classNumber ] != 0 || alphasInVertex2[ classNumber ] != 0 || alphasInVertex3[ classNumber ] != 0)
	{
        // Get the class-conditional likelihood of this class at the intensity of this pixel
        const double classConditionalLikelihood = 
                       m_Histogrammer->GetConditionalIntensityDistributions()[ classNumber ][ binNumber ];
        
        // Add contribution of the likelihood
        likelihood += classConditionalLikelihood * it.GetExtraLoadingInterpolatedValue( classIdx );
      
        //
        xGradientBasis += classConditionalLikelihood * it.GetExtraLoadingNextRowAddition( classIdx );
        yGradientBasis += classConditionalLikelihood * it.GetExtraLoadingNextColumnAddition( classIdx );
        zGradientBasis += classConditionalLikelihood * it.GetExtraLoadingNextSliceAddition( classIdx );

        classIdx++;
	}
      } // End loop over all classes
      
      
    //
    xGradientBasis /= likelihood;
    yGradientBasis /= likelihood;
    zGradientBasis /= likelihood;
#else
    
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      //std::cout << "classNumber: " << classNumber << std::endl;  
      //std::cout << "binNumber: " << binNumber << std::endl;  
        
      // Get the class-conditional likelihood of this class at the intensity of this pixel
      const double classConditionalLikelihood = 
                     m_Histogrammer->GetConditionalIntensityDistributions()[ classNumber ][ binNumber ];
      //std::cout << "classConditionalLikelihood: " << classConditionalLikelihood << std::endl;  
      const double  logClassConditionalLikelihood = log( classConditionalLikelihood + 1e-15 );               
      //std::cout << "logClassConditionalLikelihood: " << logClassConditionalLikelihood << std::endl;  
        
      // Add contribution of the likelihood
      //likelihood += classConditionalLikelihood * it.GetExtraLoadingInterpolatedValue( classNumber );
      
      //
      xGradientBasis += logClassConditionalLikelihood * it.GetExtraLoadingNextRowAddition( classNumber );
      yGradientBasis += logClassConditionalLikelihood * it.GetExtraLoadingNextColumnAddition( classNumber );
      zGradientBasis += logClassConditionalLikelihood * it.GetExtraLoadingNextSliceAddition( classNumber );
      } // End loop over all classes

#endif    
    
    //
    xGradientBasis /= m_NumberOfVoxels;
    yGradientBasis /= m_NumberOfVoxels;
    zGradientBasis /= m_NumberOfVoxels;

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
