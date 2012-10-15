/**
 * @file  kvlAtlasMeshCollectionModelLikelihoodCalculator.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#include "kvlAtlasMeshCollectionModelLikelihoodCalculator.h"

#include "kvlAtlasMeshLabelStatisticsCollector.h"



namespace kvl
{


//
//
//
AtlasMeshCollectionModelLikelihoodCalculator
::AtlasMeshCollectionModelLikelihoodCalculator()
{
  m_MeshCollection = 0;

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
  {
    return;
  }

  m_LabelImages = labelImages;
  m_NumberOfLabelImages = m_LabelImages.size();

}





//
//
//
void
AtlasMeshCollectionModelLikelihoodCalculator
::GetDataCostAndAlphasCost( float& dataCost, float& alphasCost ) const
{
  // Sanity check on input
  if ( ( !m_MeshCollection ) || ( m_LabelImages.size() == 0 ) )
  {
    itkExceptionMacro( << "No mesh collection or label images set." );
  }


  // Allocate a container to hold the total amount of pixels assigned to each vertex,
  // summed over all label images
  typedef itk::MapContainer< AtlasMesh::PointIdentifier , float >   PixelCountContainerType;
  PixelCountContainerType::Pointer  pixelCounts = PixelCountContainerType::New();
  for ( AtlasMesh::PointDataContainer::ConstIterator  pointParamIt = m_MeshCollection->GetPointParameters()->Begin();
        pointParamIt != m_MeshCollection->GetPointParameters()->End(); ++pointParamIt )
  {
    pixelCounts->InsertElement( pointParamIt.Index(), 0.0f );
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
    statisticsCollector->Rasterize( m_MeshCollection->GetMesh( labelImageNumber ) );

    // Add contribution to data cost of this label image
    dataCost += statisticsCollector->GetMinLogLikelihood();

    // Add pixel counts
    StatisticsContainerType::ConstIterator  statIt = statisticsCollector->GetLabelStatistics()->Begin();
    PixelCountContainerType::Iterator  pixelCountIt = pixelCounts->Begin();
    for ( ; statIt != statisticsCollector->GetLabelStatistics()->End(); ++statIt, ++pixelCountIt )
    {
      pixelCountIt.Value() += statIt.Value().sum();
    }

  } // End loop over all images


  // Loop over all vertices, each time adding each contribution to the alphas cost
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
  //std::cout << "combinationSelectionCost: " << combinationSelectionCost << std::endl;
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


}




} // end namespace kvl
