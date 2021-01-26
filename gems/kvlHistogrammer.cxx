#include "kvlHistogrammer.h"

#include <itkMath.h>
#include "kvlTetrahedronInteriorConstIterator.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"


namespace kvl
{

//
//
//
Histogrammer
::Histogrammer()
{

  m_MinLogLikelihood = 0;
  m_Image = 0;
  m_BinnedImage = 0;
  m_NumberOfBins = 0;
  
}


//
//
//
Histogrammer
::~Histogrammer()
{
}




//
//
//
void 
Histogrammer
::ComputeRobustRange( const ImageType* image, double& robustMin, double& robustMax )
{
  
  // Compute min and max intensity
  ImageType::PixelType  min =  itk::NumericTraits< ImageType::PixelType >::max();
  ImageType::PixelType  max =  itk::NumericTraits< ImageType::PixelType >::min();
  for ( itk::ImageRegionConstIterator< ImageType >  it( image, image->GetBufferedRegion() );
        !it.IsAtEnd(); ++it )
    {
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
Histogrammer
::UpdateBinnedImage()
{
  //
  m_NumberOfBins = m_ConditionalIntensityDistributions[ 0 ].size();
  
  
  // Compute robust range
  double  robustMin = 0.0;
  double  robustMax = 0.0;
  this->ComputeRobustRange( m_Image, robustMin, robustMax );
    
  // Create an image with binned intensities
  const double  slope = static_cast< double >( m_NumberOfBins - 1 ) / 
                        ( robustMax - robustMin );
  const double  offset = robustMin;
  m_BinnedImage = BinnedImageType::New();
  m_BinnedImage->SetRegions( m_Image->GetBufferedRegion() );
  m_BinnedImage->Allocate();
  
  itk::ImageRegionConstIterator< ImageType >  it( m_Image, m_Image->GetBufferedRegion() );
  itk::ImageRegionIterator< BinnedImageType >  binnedIt( m_BinnedImage, m_Image->GetBufferedRegion() );
  for ( ; !it.IsAtEnd(); ++it, ++binnedIt )
    {
    int  binNumber = itk::Math::Round< int >( slope * ( static_cast< double >( it.Value() ) - offset ) );

    // Intensities outside of robust range get a sentinel value
    if ( ( it.Value() < robustMin ) ||
         ( it.Value() > robustMax ) )
      {
      binNumber = -1;
      }
      
    binnedIt.Value() = binNumber;
    }
     
  
}
 


//
//
//
void
Histogrammer
::Rasterize( const AtlasMesh* mesh )
{

  // Make sure the provided input is fool proof
  const int  numberOfClasses = m_ConditionalIntensityDistributions.size();
  if ( numberOfClasses == 0 )
    {
    itkExceptionMacro( << "m_ConditionalIntensityDistributions need to be set" );  
    }
  if ( mesh->GetPointData()->Begin().Value().m_Alphas.Size() != numberOfClasses )
    {
    itkExceptionMacro( << "number of classes in m_ConditionalIntensityDistributions and provided mesh don't match" )  
    }  

  
  // Make sure we have pre-computed a binned image
  if ( !m_BinnedImage )
    {
    this->UpdateBinnedImage();
    }  
  
  // Initialize from a clean slate
  HistogramType  emptyHistogram;
  for ( int classNumber = 0; classNumber < numberOfClasses; ++classNumber )
    {
    emptyHistogram.push_back( std::vector< double >( m_NumberOfBins, 0.0 ) );  
    }
  m_Histogram = emptyHistogram;
  m_MinLogLikelihood = 0.0;
  m_ThreadSpecificHistograms.clear();
  m_ThreadSpecificMinLogLikelihoods.clear();

  // Initialize thread-specific histogram from a clean slate
  HistogramThreadAccumType  emptyThreadHistogram;
  for ( int classNumber = 0; classNumber < numberOfClasses; ++classNumber )
    {
    emptyThreadHistogram.push_back( std::vector< ThreadAccumDataType >( m_NumberOfBins, 0.0 ) );  
    }

  // For each thread, create an empty histogram and cost so that
  // different threads never interfere with one another
  for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
    {
    // Initialize cost to zero for this thread
    m_ThreadSpecificMinLogLikelihoods.push_back( 0.0 );  
      
    // Initialize to zero-filled histogram
    m_ThreadSpecificHistograms.push_back( emptyThreadHistogram );
    } // End loop over threads
    
  // Now rasterize
  Superclass::Rasterize( mesh );

  // Collect the results of all threads and copy to final MinLogLikelihood
  for ( int i = 1 ; i < m_ThreadSpecificMinLogLikelihoods.size() ; i++ )
    {
    m_ThreadSpecificMinLogLikelihoods[0] += m_ThreadSpecificMinLogLikelihoods[i];
    }
  m_MinLogLikelihood = m_ThreadSpecificMinLogLikelihoods[0];

  // Collect the results of all the histogram threads
  for ( int i = 1 ; i < m_ThreadSpecificHistograms.size() ; i++ )
    {
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      for ( int binNumber = 0; binNumber < m_NumberOfBins; binNumber++ )
        {
        m_ThreadSpecificHistograms[ 0 ][ classNumber ][ binNumber ] += 
            m_ThreadSpecificHistograms[ i ][ classNumber ][ binNumber ];
        }  
      }
    }

  // Copy to final histogram
  for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
    for ( int binNumber = 0; binNumber < m_NumberOfBins; binNumber++ )
      {
      m_Histogram[ classNumber ][ binNumber ] =
          m_ThreadSpecificHistograms[ 0 ][ classNumber ][ binNumber ]; 
      }  
    }
    
}    
  
    
  

//
//
//
bool
Histogrammer
::RasterizeTetrahedron( const AtlasMesh* mesh, 
                        AtlasMesh::CellIdentifier tetrahedronId,
                        int threadNumber )
{
  // Retrieve necessary info about tetrahedron
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
  TetrahedronInteriorConstIterator< BinnedImageType::PixelType >  it( m_BinnedImage, p0, p1, p2, p3 );
  const AtlasAlphasType&  alphasInVertex0 = mesh->GetPointData()->ElementAt( id0 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex1 = mesh->GetPointData()->ElementAt( id1 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex2 = mesh->GetPointData()->ElementAt( id2 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex3 = mesh->GetPointData()->ElementAt( id3 ).m_Alphas;
  const int  numberOfClasses = m_ConditionalIntensityDistributions.size();
  for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
    it.AddExtraLoading( alphasInVertex0[ classNumber ], 
                        alphasInVertex1[ classNumber ], 
                        alphasInVertex2[ classNumber ], 
                        alphasInVertex3[ classNumber ] );
    }
  std::vector< double >  unnormalizedPosterior( numberOfClasses, 0.0 );
  for ( ; !it.IsAtEnd(); ++it )
    {
    //
    const int  binNumber = it.Value();

    // Skip sentinel values
    if ( binNumber < 0 )
      {
      continue;  
      }  
      
    //
    double  denominator = 1e-15;
    double  weightToDistribute = 1e-15;
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      const double  tmp =  m_ConditionalIntensityDistributions[ classNumber ][ binNumber ] // Likelihood   
                           *
                           it.GetExtraLoadingInterpolatedValue( classNumber ); // Prior
      unnormalizedPosterior[ classNumber ] = tmp;
      denominator += tmp;
      weightToDistribute += it.GetExtraLoadingInterpolatedValue( classNumber ); // Prior
      }
    m_ThreadSpecificMinLogLikelihoods[ threadNumber ] -= weightToDistribute * log( denominator );
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      m_ThreadSpecificHistograms[ threadNumber ][ classNumber ][ binNumber ] 
           += weightToDistribute * unnormalizedPosterior[ classNumber ] / denominator;  
      }
    
    } // End loop over all pixels within tetrahedron  

  return true;
}






} // end namespace kvl
