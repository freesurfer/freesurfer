#include "kvlMutualInformationCostAndGradientCalculator.h"

#include <itkMath.h>
#include "kvlTetrahedronInteriorConstIterator.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include "itkImageFileWriter.h"


namespace kvl
{

//
//
//
MutualInformationCostAndGradientCalculator
::MutualInformationCostAndGradientCalculator()
{
  
  m_BinnedImage = 0;
  
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
::ComputeRobustRange( const ImageType* image, double& robustMin, double& robustMax )
{
  
  // Compute min and max intensity
  ImageType::PixelType  min =  itk::NumericTraits< ImageType::PixelType >::max();
  ImageType::PixelType  max =  itk::NumericTraits< ImageType::PixelType >::min();
  for ( itk::ImageRegionConstIterator< ImageType >  it( image, image->GetBufferedRegion() );
        !it.IsAtEnd(); ++it )
    {
    // Skip zeroes
    if ( it.Value() == 0 )
      {
      continue;  
      }  
      
    if ( it.Value() < min )
      {
      min = it.Value();
      }

    if ( it.Value() > max )
      {
      max = it.Value();
      }
    }
    
    
  // Build histogram
  const int  numberOfBins = 4000;
  const double  slope = static_cast< double >( numberOfBins - 1 ) / 
                        ( static_cast< double >( max ) - static_cast< double >( min ) );
  const double  offset = static_cast< double >( min );                      
  std::vector< double >  histogram( numberOfBins, 0.0 );
  double  sumOfHistogram = 0.0;
  for ( itk::ImageRegionConstIterator< ImageType >  it( image, image->GetBufferedRegion() );
        !it.IsAtEnd(); ++it )
    {
    // Skip zeroes
    if ( it.Value() == 0 )
      {
      continue;  
      }  
      
    const int  binNumber = itk::Math::Round< int >( slope * ( static_cast< double >( it.Value() ) - offset ) );
    histogram[ binNumber ]++;
    sumOfHistogram++;
    }
    
  // Compute cumulative histogram
  std::vector< double >  cumulativeHistogram( numberOfBins, 0.0 );
  cumulativeHistogram[ 0 ] = histogram[ 0 ] / sumOfHistogram;
  for ( int binNumber = 1; binNumber < numberOfBins; binNumber++ )
    {
    cumulativeHistogram[ binNumber ] = cumulativeHistogram[ binNumber-1 ] + histogram[ binNumber ] / sumOfHistogram;
    }  
  
  // Determine bin number of robust min and max
  int  robustMinBinNumber = 0;
  int  robustMaxBinNumber = 0;
  for ( int binNumber = ( numberOfBins-1 ); binNumber >= 0; binNumber-- )
    {
    if ( cumulativeHistogram[ binNumber ] > 0.9995 )
      {
      robustMaxBinNumber = binNumber;
      }
      
    if ( cumulativeHistogram[ binNumber ] > 0.0005 )
      {
      robustMinBinNumber = binNumber;
      }
      
    }
    
  // Determine robust min and max
  robustMin = static_cast< double >( robustMinBinNumber ) / slope + offset;
  robustMax = static_cast< double >( robustMaxBinNumber ) / slope + offset;
  
  //
  std::cout << "min: " << min << std::endl;
  std::cout << "max: " << max << std::endl;
  std::cout << "robustMin: " << robustMin << std::endl;
  std::cout << "robustMax: " << robustMax << std::endl;

  for ( int binNumber = 0; binNumber < 3; binNumber++ )
    {
    std::cout << "histogram[ " << binNumber << " ]: " << histogram[ binNumber ] << std::endl;
    }
  for ( int binNumber = (numberOfBins-3); binNumber < numberOfBins; binNumber++ )
    {
    std::cout << "histogram[ " << binNumber << " ]: " << histogram[ binNumber ] << std::endl;
    }

  for ( int binNumber = 0; binNumber < 3; binNumber++ )
    {
    std::cout << "cumulativeHistogram[ " << binNumber << " ]: " << cumulativeHistogram[ binNumber ] << std::endl;
    }
  for ( int binNumber = (numberOfBins-3); binNumber < numberOfBins; binNumber++ )
    {
    std::cout << "cumulativeHistogram[ " << binNumber << " ]: " << cumulativeHistogram[ binNumber ] << std::endl;
    }
  
  
  
}


//
//
//
void 
MutualInformationCostAndGradientCalculator
::SetImage( const ImageType* image )
{
  // Compute robust range
  double  robustMin = 0.0;
  double  robustMax = 0.0;
  this->ComputeRobustRange( image, robustMin, robustMax );
    
  // Create an image with binned intensities
  const int  numberOfBins = 32;
  const double  slope = static_cast< double >( numberOfBins - 1 ) / 
                        ( robustMax - robustMin );
  const double  offset = robustMin;
  m_BinnedImage = BinnedImageType::New();
  m_BinnedImage->SetRegions( image->GetBufferedRegion() );
  m_BinnedImage->Allocate();
  
  itk::ImageRegionConstIterator< ImageType >  it( image, image->GetBufferedRegion() );
  itk::ImageRegionIterator< BinnedImageType >  binnedIt( m_BinnedImage, image->GetBufferedRegion() );
  for ( ; !it.IsAtEnd(); ++it, ++binnedIt )
    {
    int  binNumber = itk::Math::Round< int >( slope * ( static_cast< double >( it.Value() ) - offset ) );

    // Zeroes or intensities outside of robust range get a sentinel value
    if ( ( it.Value() == 0 ) || 
         ( it.Value() < robustMin ) ||
         ( it.Value() > robustMax ) )
      {
      binNumber = -1;
      }
      
    binnedIt.Value() = binNumber;
    }
     
  //
  typedef itk::ImageFileWriter< BinnedImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetInput( m_BinnedImage );
  writer->SetFileName( "binnedImage.nii" );
  writer->Write();
  
         
}
 


//
//
//
void
MutualInformationCostAndGradientCalculator
::Rasterize( const AtlasMesh* mesh )
{

  // Estimate parameters of generative model using EM
  
  
  
  

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
