#include "kvlAtlasMeshCollectionModelLikelihoodCalculator.h"

#include "kvlAtlasMeshLabelStatisticsCollector.h"
#include "itkTimeProbe.h"



namespace kvl
{


//
//
//
AtlasMeshCollectionModelLikelihoodCalculator
::AtlasMeshCollectionModelLikelihoodCalculator()
{
  m_MeshCollection = 0;
  m_mapCompToComp = 0;
}
  
  

//
//
//
AtlasMeshCollectionModelLikelihoodCalculator
::~AtlasMeshCollectionModelLikelihoodCalculator()
{

}




//
//
//
void
AtlasMeshCollectionModelLikelihoodCalculator
::SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages )
{
  if ( labelImages.size() == 0 )
    return;

  m_LabelImages = labelImages;
  m_NumberOfLabelImages = m_LabelImages.size();

}





//
//
//
void
AtlasMeshCollectionModelLikelihoodCalculator
::GetDataCostAndAlphasCost( float& dataCost, float& alphasCost, bool verbose ) const
{
  // Sanity check on input
  if ( ( !m_MeshCollection ) || ( m_LabelImages.size() == 0 ) )
    {
    itkExceptionMacro( << "No mesh collection or label images set." );
    }


  // Allocate a container to hold the total amount of pixels assigned to each vertex, 
  // summed over all label images
  itk::TimeProbe  probe;
  if ( verbose )
    {
    probe.Start();
    // std::cout << "AtlasMeshCollectionModelLikelihoodCalculator: allocating pixel count container..." << std::flush;
    }

  typedef itk::MapContainer< AtlasMesh::PointIdentifier , float >   PixelCountContainerType;
  PixelCountContainerType::Pointer  pixelCounts = PixelCountContainerType::New();
  for ( AtlasMesh::PointDataContainer::ConstIterator  pointParamIt = m_MeshCollection->GetPointParameters()->Begin();
        pointParamIt != m_MeshCollection->GetPointParameters()->End(); ++pointParamIt )
    {
    pixelCounts->InsertElement( pointParamIt.Index(), 0.0f );
    }

  if ( verbose )
    {
    probe.Stop();
    // std::cout << "done (" << probe.GetMeanTime() << " seconds)" << std::endl;
    }


  // Loop over all label images, adding each image's contribution to the data cost and to 
  // the total amount of pixels assigned to each vertex
  dataCost = 0;
  typedef AtlasMeshLabelStatisticsCollector::StatisticsContainerType StatisticsContainerType;
  AtlasMeshLabelStatisticsCollector::Pointer  statisticsCollector = AtlasMeshLabelStatisticsCollector::New();
  for ( int  labelImageNumber = 0; labelImageNumber < m_NumberOfLabelImages; labelImageNumber++ )
    {
    // Calculate statistics for this label image
    statisticsCollector->SetLabelImage( m_LabelImages[ labelImageNumber ] );
    statisticsCollector->SetMapCompToComp( m_mapCompToComp );
    if ( verbose )
      {
      probe = itk::TimeProbe();
      probe.Start();
      // std::cout << "AtlasMeshCollectionModelLikelihoodCalculator: calculating statistics for label image " << labelImageNumber << "..." << std::flush;
      }
  
    statisticsCollector->Rasterize( m_MeshCollection->GetMesh( labelImageNumber ) );

    if ( verbose )
      {
      probe.Stop();
      // std::cout << "done (" << probe.GetMeanTime() << " seconds)" << std::endl;
      // std::cout << "AtlasMeshCollectionModelLikelihoodCalculator: number of threads used was " << statisticsCollector->GetNumberOfThreads() << std::endl;

      probe = itk::TimeProbe();
      probe.Start();
      // std::cout << "AtlasMeshCollectionModelLikelihoodCalculator: addition likelihood contribution from label image " << labelImageNumber << "..." << std::flush;
      }

    // Add contribution to data cost of this label image 
    dataCost += statisticsCollector->GetMinLogLikelihood();

    if ( verbose )
      {
      probe.Stop();
      // std::cout << "done (" << probe.GetMeanTime() << " seconds)" << std::endl;

      probe = itk::TimeProbe();
      probe.Start();
      // std::cout << "AtlasMeshCollectionModelLikelihoodCalculator: retreiving statistics from label image " << labelImageNumber << "..." << std::flush;
      }
    
    // Add pixel counts
    StatisticsContainerType::ConstIterator  statIt = statisticsCollector->GetLabelStatistics()->Begin();
    if ( verbose )
      {
      probe.Stop();
      // std::cout << "done (" << probe.GetMeanTime() << " seconds)" << std::endl;
      
      probe = itk::TimeProbe();
      probe.Start();
      // std::cout << "AtlasMeshCollectionModelLikelihoodCalculator: adding statistics from label image " << labelImageNumber << "..." << std::flush;
     }

    PixelCountContainerType::Iterator  pixelCountIt = pixelCounts->Begin();
    for ( ; statIt != statisticsCollector->GetLabelStatistics()->End(); ++statIt, ++pixelCountIt )
      {
      pixelCountIt.Value() += statIt.Value().sum();
      }
      
    if ( verbose )
      {
      probe.Stop();
      // std::cout << "done (" << probe.GetMeanTime() << " seconds)" << std::endl;
      }

      
    } // End loop over all images


  // Loop over all vertices, each time adding each contribution to the alphas cost
  if ( verbose )
    {
    probe = itk::TimeProbe();
    probe.Start();
    // std::cout << "AtlasMeshCollectionModelLikelihoodCalculator: computing alphas cost..." << std::flush;
    }

  alphasCost = 0;
  unsigned int  numberOfLabels = m_MeshCollection->GetPointParameters()->Begin().Value().m_Alphas.Size();
#if 1
  // We don't "allow" more than three labels simultaneously to have non-zero alpha. The alpha cost is therefore
  // the cost of sending three alpha values, PLUS the cost of sending which combination we are sending. The last
  // cost is simply -log( 1/C ), with C the number of combinations of three classes in the total amount of classes:
  // C = numberOfLabels! / 3! / ( numberOfLabels - 3)!
  float  combinationSelectionCost = 0.0f;
  if ( numberOfLabels > 3 )
    {
    for ( unsigned int i = numberOfLabels-2; i <= numberOfLabels; i++ )
      {
      combinationSelectionCost += log( i );
      }
    for ( int i = 1; i <= 3; i++ )
      {
      combinationSelectionCost -= log( i );
      }
    numberOfLabels = 3;
    }
  //// std::cout << "combinationSelectionCost: " << combinationSelectionCost << std::endl;
#endif
  PixelCountContainerType::ConstIterator  pixelCountIt = pixelCounts->Begin();
  AtlasMesh::PointDataContainer::ConstIterator  pointParamIt = m_MeshCollection->GetPointParameters()->Begin();
  for ( ; pixelCountIt != pixelCounts->End(); ++pixelCountIt, ++pointParamIt )
    {
#if 0
    if ( !pointParamIt.Value().m_CanChangeAlphas )
      {
      continue;
      }
#endif   

    const float  totalNumberOfPixels = pixelCountIt.Value();
    for ( unsigned int i = 1; i < numberOfLabels; i++ )
      {
      alphasCost += log( totalNumberOfPixels + i );
      alphasCost -= log( i );
      }
#if 1
    alphasCost += combinationSelectionCost;
#endif

    } // End loop over all points

  if ( verbose )
    {
    probe.Stop();
    // std::cout << "done (" << probe.GetMeanTime() << " seconds)" << std::endl;
    }

}




} // end namespace kvl
