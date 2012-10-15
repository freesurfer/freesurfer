/**
 * @file  kvlAtlasParameterEstimator.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
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
#include "kvlAtlasParameterEstimator.h"

#include "itkImageRegionConstIterator.h"
#include "kvlAtlasMeshLabelStatisticsCollector.h"
#include "kvlAtlasMeshPositionGradientCalculator.h"
#include "kvlAtlasMeshMinLogLikelihoodCalculator.h"
#include "kvlAtlasMeshSmoother.h"


namespace kvl
{


//
//
//
AtlasParameterEstimator
::AtlasParameterEstimator()
{
  m_MeshCollection = 0;
  m_NumberOfClasses = 0;
  m_IterationNumber = 0;
  m_MaximumNumberOfIterations = 300;
  m_LabelImageNumber = 0;
  m_NumberOfLabelImages = 0;
  m_AlphasEstimationIterationNumber = 0;
  m_AlphasEstimationMaximumNumberOfIterations = 20;
  m_PositionEstimationIterationNumber = 0;
  m_PositionEstimationMaximumNumberOfIterations = 200;
  m_PositionGradientDescentStepSize = 1.0f;
  m_PositionEstimationIterationEventResolution = 10;
  m_CurrentMinLogLikelihoodTimesPrior = 0.0f;
  m_AlphaEstimationStopCriterion = 0.0005f;
  m_PositionEstimationStopCriterion = 0.00005f;
  m_AlphasSmoothingFactor = 0.0f;
  m_StopCriterion = 0.001f;

  m_UseGaussians = false;
  m_IgnoreLastLabelImage = false;
}



//
//
//
AtlasParameterEstimator
::~AtlasParameterEstimator()
{

}



//
//
//
void
AtlasParameterEstimator
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{



}



//
//
//
const AtlasParameterEstimator::LabelImageType*
AtlasParameterEstimator
::GetLabelImage( unsigned int labelImageNumber ) const
{
  // Sanity check
  if ( labelImageNumber >= m_LabelImages.size() )
  {
    return 0;
  }

  return m_LabelImages[ labelImageNumber ];
}




//
//
//
void
AtlasParameterEstimator
::SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages )
{
  if ( labelImages.size() == 0 )
  {
    return;
  }

  m_LabelImages = labelImages;
  m_NumberOfLabelImages = m_LabelImages.size();

  // Calculate the number of classes present in the label images
  LabelImageType::PixelType  maximum = 0;
  for ( unsigned int labelImageNumber = 0; labelImageNumber < m_NumberOfLabelImages; labelImageNumber++ )
  {

    // Get maximum intensity value
    typedef itk::ImageRegionConstIterator< LabelImageType >   IteratorType;
    IteratorType  it( labelImages[ labelImageNumber ],
                      labelImages[ labelImageNumber ]->GetLargestPossibleRegion() );
    it = it.Begin();
    for (; !it.IsAtEnd(); ++it)
    {
      if ( it.Get() > maximum )
      {
        maximum = it.Get();
      }
    }

  }
  m_NumberOfClasses = static_cast< unsigned int >( maximum ) + 1;
  //std::cout << "Calculated number of classes to be: " << m_NumberOfClasses << std::endl;


}



//
//
//
void
AtlasParameterEstimator
::Estimate()
{
  // Sanity checking
  if ( m_MeshCollection->GetNumberOfMeshes() == 0 )
  {
    itkExceptionMacro( "No initial meshes set!" );
  }
  if ( m_LabelImages.size() == 0 )
  {
    itkExceptionMacro( "No label images set!" );
  }
  if ( m_MeshCollection->GetNumberOfMeshes() != m_LabelImages.size() )
  {
    itkExceptionMacro( "Number of meshes must match number of label images!" );
  }


  // Main loop
  //const float  stopCriterion = 0.001 /* m_PositionEstimationStopCriterion * m_LabelImages.size() */;
  float  previousMinLogLikelihoodTimesPrior = itk::NumericTraits< float >::max();
  m_CurrentMinLogLikelihoodTimesPrior = previousMinLogLikelihoodTimesPrior / 2;
  m_IterationNumber = 0;
  this->InvokeEvent( itk::StartEvent() );
  while ( ( ( ( previousMinLogLikelihoodTimesPrior - m_CurrentMinLogLikelihoodTimesPrior ) /
              fabsf( m_CurrentMinLogLikelihoodTimesPrior ) ) > m_StopCriterion ) &&
          ( m_IterationNumber < m_MaximumNumberOfIterations ) )
  {
    // Estimate alphas
    this->EstimateAlphas();

    // Smooth results
    this->SmoothAlphas();

    // Register. This will also provide us with the minLogLikelihoodTimesPrior of
    // the current parameter set
    previousMinLogLikelihoodTimesPrior = m_CurrentMinLogLikelihoodTimesPrior;
    m_CurrentMinLogLikelihoodTimesPrior = this->EstimatePositions();

    // Prepare for the next iteration
    //std::cout << "================ IterationNumber: " << m_IterationNumber
    //          << " CurrentMinLogLikelihoodTimesPrior: " << m_CurrentMinLogLikelihoodTimesPrior
    //          << "================" << std::endl;
    m_IterationNumber++;
    this->InvokeEvent( itk::IterationEvent() );
  }

  // Estimate final alphas
  this->EstimateAlphas();

  this->InvokeEvent( itk::EndEvent() );


}



//
//
//
void
AtlasParameterEstimator
::EstimateAlphas()
{

  // Nothing to estimate if only one image which doesn't contribute
  if ( m_IgnoreLastLabelImage && ( m_NumberOfLabelImages == 1 ) )
  {
    return;
  }



  // If we're smoothing, we have to start with flat alpha entries
  if ( m_AlphasSmoothingFactor != 0 )
  {
    m_MeshCollection->FlattenAlphas();
  }


  // Allocate a container to sum the statistics over all label images
  typedef AtlasMeshLabelStatisticsCollector::StatisticsContainerType StatisticsContainerType;
  StatisticsContainerType::Pointer  pooledStatistics = StatisticsContainerType::New();
  AtlasMesh::PointDataContainer::ConstIterator  pointParamIt = m_MeshCollection->GetPointParameters()->Begin();
  while ( pointParamIt != m_MeshCollection->GetPointParameters()->End() )
  {
    pooledStatistics->CreateIndex( pointParamIt.Index() );

    ++pointParamIt;
  }


  float  previousCost = itk::NumericTraits< float >::max();
  float  currentCost = previousCost / 2;
  m_AlphasEstimationIterationNumber = 0;
  this->InvokeEvent( AlphasEstimationStartEvent() );

  while ( ( ( ( previousCost - currentCost ) / fabsf( currentCost ) ) > m_AlphaEstimationStopCriterion ) &&
          ( m_AlphasEstimationIterationNumber < m_AlphasEstimationMaximumNumberOfIterations ) )
  {
    previousCost = currentCost;

    // Initialize pooled statistics to zero
    unsigned int  numberOfLabeles = m_MeshCollection->GetPointParameters()->Begin().Value().m_Alphas.Size();
    StatisticsContainerType::Element  zeroEntry( numberOfLabeles );
    zeroEntry.Fill( 0.0f );
    StatisticsContainerType::Iterator  pooledIt = pooledStatistics->Begin();
    while ( pooledIt != pooledStatistics->End() )
    {
      pooledIt.Value() = zeroEntry;
      pooledIt++;
    }


    // Loop over all label images, retrieve their statistics, and add those to the pooled statistics
    currentCost = 0;
    int  numberOfLabelImagesToVisit = m_NumberOfLabelImages;
    if ( m_IgnoreLastLabelImage )
    {
      std::cout << "Ignoring last image in the alphas estimation" << std::endl;
      numberOfLabelImagesToVisit--;
    }
    AtlasMeshLabelStatisticsCollector::Pointer  statisticsCollector = AtlasMeshLabelStatisticsCollector::New();
    for ( int  labelImageNumber = 0; labelImageNumber < numberOfLabelImagesToVisit; labelImageNumber++ )
    {
      // Calculate statistics for this label image
      statisticsCollector->SetLabelImage( m_LabelImages[ labelImageNumber ] );
      statisticsCollector->Rasterize( m_MeshCollection->GetMesh( labelImageNumber ) );

      // Add statistics to pooled statistics
      pooledIt = pooledStatistics->Begin();
      StatisticsContainerType::ConstIterator  statIt = statisticsCollector->GetLabelStatistics()->Begin();
      while ( pooledIt != pooledStatistics->End() )
      {
        pooledIt.Value() += statIt.Value();

        pooledIt++;
        statIt++;
      }

      // Also add contribution of cost function
      currentCost += statisticsCollector->GetMinLogLikelihood();
    }


    if ( !m_UseGaussians )
    {
      // Estimated alphas are simply normalized pooled statistics
      //std::cout << "Not using gaussians" << std::endl;

      pooledIt = pooledStatistics->Begin();
      PointDataContainerType::Iterator  pointParamIt = m_MeshCollection->GetPointParameters()->Begin();
      while ( pooledIt != pooledStatistics->End() )
      {
        if ( pointParamIt.Value().m_CanChangeAlphas )
        {
          pointParamIt.Value().m_Alphas = pooledIt.Value();
          pointParamIt.Value().m_Alphas /= pooledIt.Value().sum();
        }

        pooledIt++;
        pointParamIt++;
      }

    }
    else
    {
      // Estimated alphas are assumed Gaussian distributed. Calculate the mean and variance, and impose
      // this on the alphas
      std::cout << "Using gaussians" << std::endl;

      const int  numberOfLabels = m_MeshCollection->GetPointParameters()->Begin().Value().m_Alphas.Size();
      AtlasAlphasType  tmp1( numberOfLabels );
      AtlasAlphasType  tmp2( numberOfLabels );
      for ( int labelNumber = 0; labelNumber < numberOfLabels; labelNumber++ )
      {
        tmp1[ labelNumber ] = labelNumber;
        tmp2[ labelNumber ] = labelNumber * labelNumber;
      }

#if 1
      pooledIt = pooledStatistics->Begin();
      PointDataContainerType::Iterator  pointParamIt = m_MeshCollection->GetPointParameters()->Begin();
      while ( pooledIt != pooledStatistics->End() )
      {
        // Calculate the mean and variance
        if ( pointParamIt.Value().m_CanChangeAlphas )
        {
          const float  numberOfSamples = pooledIt.Value().sum();
          const float  mean = dot_product( pooledIt.Value(), tmp1 ) / numberOfSamples;
          //const float  variance = dot_product( pooledIt.Value(), tmp2 ) / numberOfSamples - mean * mean;
          const float  variance = 2.0f;
          //std::cout << "For this node I have mean " << mean << " and sigma " << sqrt( variance )
          //          << " (numberOfSamples: " << numberOfSamples << ")" << std::endl;

          // Now inforce a Gaussian on the labels
          for ( int labelNumber = 0; labelNumber < numberOfLabels; labelNumber++ )
          {
            pointParamIt.Value().m_Alphas[ labelNumber ] = exp( -pow( labelNumber - mean, 2 ) / 2.0f / variance ) /
                sqrt( 2 * 3.14 * variance );
          }

        }


        pooledIt++;
        pointParamIt++;
      }
#else

      // First calculate the variance by looping over all nodes
      pooledIt = pooledStatistics->Begin();
      PointDataContainerType::Iterator  pointParamIt = m_MeshCollection->GetPointParameters()->Begin();
      float  varianceNumerator = 0.0f;
      float  varianceDenominator = 0.0f;
      for ( ; pooledIt != pooledStatistics->End(); ++pooledIt, ++pointParamIt )
      {
        // Calculate the mean and variance
        if ( pointParamIt.Value().m_CanChangeAlphas )
        {
          const float  numberOfSamples = pooledIt.Value().sum();
          const float  mean = dot_product( pooledIt.Value(), tmp1 ) / numberOfSamples;
          varianceNumerator += dot_product( pooledIt.Value(), tmp2 ) - mean * mean * numberOfSamples;
          varianceDenominator += numberOfSamples;
        }

      }
      const float  variance = varianceNumerator / varianceDenominator;
      std::cout << "Calculated variance to be " << sqrt( variance ) << "^2" << std::endl;
      std::cout << "     varianceNumerator: " << varianceNumerator << std::endl;
      std::cout << "     varianceDenominator: " << varianceDenominator << std::endl;


      // Now loop again over all nodes, and inforce a Gaussian on the labels with the
      // variance with just estimated
      for ( pointParamIt = m_MeshCollection->GetPointParameters()->Begin(), pooledIt = pooledStatistics->Begin();
            pointParamIt != m_MeshCollection->GetPointParameters()->End(); ++pointParamIt, ++pooledIt )
      {
        // Calculate the mean and variance
        if ( pointParamIt.Value().m_CanChangeAlphas )
        {
          const float  numberOfSamples = pooledIt.Value().sum();
          const float  mean = dot_product( pooledIt.Value(), tmp1 ) / numberOfSamples;

          for ( int labelNumber = 0; labelNumber < numberOfLabels; labelNumber++ )
          {
            pointParamIt.Value().m_Alphas[ labelNumber ] = exp( -pow( labelNumber - mean, 2 ) / 2.0f / variance ) /
                sqrt( 2 * 3.14 * variance );
          }

        }

      }
      std::cout << "Successfully inforced Gaussian distribution" << std::endl;

#endif

    }

    // Prepare for next iteration
    m_AlphasEstimationIterationNumber++;
    this->InvokeEvent( AlphasEstimationIterationEvent() );

  } // End iterative loop

  this->InvokeEvent( AlphasEstimationEndEvent() );


}





//
//
//
void
AtlasParameterEstimator
::SmoothAlphas()
{

  if ( m_AlphasSmoothingFactor != 0 )
  {
    // Retrieve sigma of the smoothing kernel
    const float  maximalSigma = 0.05 * ( m_LabelImages[ 0 ]->GetLargestPossibleRegion().GetSize()[ 0 ] );
    const float  sigma = maximalSigma * ( 1 - log( 2 - m_AlphasSmoothingFactor ) / log( 2 ) );

    // Smooth
    AtlasMeshSmoother::Pointer  smoother = AtlasMeshSmoother::New();
    smoother->SetMeshCollection( m_MeshCollection );
    smoother->SetSigma( sigma );
    m_MeshCollection = smoother->GetSmoothedMeshCollection();

#if 1
    {
      std::cout << "sigma: " << sigma << std::endl;
      std::ostringstream  fileNameStream;
      fileNameStream << "smoothed_" << sigma << ".txt";
      m_MeshCollection->Write( fileNameStream.str().c_str() );
    }
#endif
  }

}




//
//
//
float
AtlasParameterEstimator
::EstimatePositions()
{

  float  totalMinLogLikelihoodTimesPrior = 0.0f;

  for ( m_LabelImageNumber = 0; m_LabelImageNumber < m_NumberOfLabelImages; m_LabelImageNumber++ )
  {
    //
    if ( m_IgnoreLastLabelImage && ( m_LabelImageNumber == ( m_NumberOfLabelImages - 1 ) ) )
    {
      std::cout << "Doing the registration of the last image, but it doesn't count towards the cost function" << std::endl;
      this->EstimatePosition( m_LabelImageNumber );
    }
    else
    {
      totalMinLogLikelihoodTimesPrior += this->EstimatePosition( m_LabelImageNumber );
    }
  }

  return totalMinLogLikelihoodTimesPrior;
}



//
//
//
float
AtlasParameterEstimator
::EstimatePosition( unsigned int  labelImageNumber )
{

  // Let's check if we are expected to just stay put
  if ( m_MeshCollection->GetK() >= 1000 )
  {
#if 0
    // Just copy the reference position
    PointsContainerType::ConstIterator  refPosIt = m_MeshCollection->GetReferencePosition()->Begin();
    PointsContainerType::Iterator  posIt = m_MeshCollection->GetPositions()[ labelImageNumber ]->Begin();
    for ( ; refPosIt != m_MeshCollection->GetReferencePosition()->End(); ++refPosIt, ++posIt )
    {
      posIt.Value() = refPosIt.Value();
    }
#endif

    // Now return the cost at this position
    float  cost;
    this->CalculateCurrentPositionGradient( labelImageNumber, cost );
    return cost;
  }


  float  previousCost = itk::NumericTraits< float >::max();
  float  currentCost = previousCost / 2;
  AtlasPositionGradientContainerType::Pointer  previousGradient = 0;
  AtlasPositionGradientContainerType::Pointer  currentGradient = 0;
  this->InvokeEvent( PositionEstimationStartEvent() );

  for ( m_PositionEstimationIterationNumber = 0;
        m_PositionEstimationIterationNumber < m_PositionEstimationMaximumNumberOfIterations;
        m_PositionEstimationIterationNumber++ )
  {

    // Get current gradient and current cost
    previousCost = currentCost;
    previousGradient = currentGradient;
    currentGradient = this->CalculateCurrentPositionGradient( labelImageNumber, currentCost );

    //std::cout << "Cost going from " << previousCost << " to " << currentCost << std::endl;

    // Compare the current cost with the previous cost, and decide what to do
    if ( ( currentCost > previousCost ) && ( m_PositionEstimationIterationNumber != 0 ) )
    {
      // Last step was wrong. Undo last move and stop
      //std::cout << "!!!!! Going uphill. Undoing last move and stopping." << std::endl;
#if 1
      PointsContainerType::Iterator  posIt = m_MeshCollection->GetPositions()[ labelImageNumber ]->Begin();
      AtlasPositionGradientContainerType::ConstIterator  gradIt = previousGradient->Begin();
      while ( posIt != m_MeshCollection->GetPositions()[ labelImageNumber ]->End() )
      {
        posIt.Value() += gradIt.Value();

        ++posIt;
        ++gradIt;
      }

      currentCost = previousCost;
      break;
#endif
    }
    else if ( ( ( previousCost - currentCost ) / fabsf( currentCost ) ) < m_PositionEstimationStopCriterion )
    {
      // Converged. Stop
      //std::cout << "!!!!! Converged. Stopping." << std::endl;
      break;
    }
    else
    {
      //std::cout << "!!!!! Applying gradient " << std::endl;

      // Calculate the maximal gradient magnitude of the current gradient

      AtlasMesh::PointDataContainer::ConstIterator  pointParamIt =
        m_MeshCollection->GetPointParameters()->Begin();
      AtlasPositionGradientContainerType::Iterator  gradIt = currentGradient->Begin();
      float  maximumGradientMagnitude = 0.0f;
      while ( gradIt != currentGradient->End() )
      {
        // Reset gradient components to zero for points that are immobile
        if ( pointParamIt.Value().m_CanMoveX == false )
        {
          gradIt.Value()[ 0 ] = 0;
        }

        if ( pointParamIt.Value().m_CanMoveY == false )
        {
          gradIt.Value()[ 1 ] = 0;
        }

        if ( pointParamIt.Value().m_CanMoveZ == false )
        {
          gradIt.Value()[ 2 ] = 0;
        }

        // Calculate the magnitude and compare with the maximal value so far
        float  magnitude = gradIt.Value().GetNorm();
        if ( magnitude > maximumGradientMagnitude )
        {
          maximumGradientMagnitude = magnitude;
        }

        ++pointParamIt;
        ++gradIt;
      }

      if ( maximumGradientMagnitude > 0 )
      {
        // Scale current gradient, and apply
        PointsContainerType::Iterator  posIt = m_MeshCollection->GetPositions()[ labelImageNumber ]->Begin();
        gradIt = currentGradient->Begin();
        while ( posIt != m_MeshCollection->GetPositions()[ labelImageNumber ]->End() )
        {
          //std::cout << "Applying gradient to point " << posIt.Index() << std::endl;
          //std::cout << "      position: " << posIt.Value() << std::endl;
          gradIt.Value() *= ( m_PositionGradientDescentStepSize / maximumGradientMagnitude );
          //std::cout << "      gradient: " << gradIt.Value() << std::endl;
          posIt.Value() -= gradIt.Value();

          ++posIt;
          ++gradIt;
        }
      }
      else
      {
        //std::cout << "        (Actually, gradient is 0 everywhere)" << std::endl;
      }

    }


    if ( !( m_PositionEstimationIterationNumber % m_PositionEstimationIterationEventResolution ) )
    {
      this->InvokeEvent( PositionEstimationIterationEvent() );
    }

  }


  this->InvokeEvent( PositionEstimationEndEvent() );

  return currentCost;

}







//
//
//
AtlasPositionGradientContainerType::Pointer
AtlasParameterEstimator
::CalculateCurrentPositionGradient( unsigned int labelImageNumber, float& minLogLikelihoodTimesPrior ) const
{
  // Sanity checking
  if ( ( labelImageNumber >= m_MeshCollection->GetNumberOfMeshes() ) ||
       ( labelImageNumber >= m_LabelImages.size() ) )
  {
    return 0;
  }


  AtlasMeshPositionGradientCalculator::Pointer  gradientCalculator = AtlasMeshPositionGradientCalculator::New();
  gradientCalculator->SetLabelImage( m_LabelImages[ labelImageNumber ] );
  gradientCalculator->Rasterize( m_MeshCollection->GetMesh( labelImageNumber ) );
  minLogLikelihoodTimesPrior = gradientCalculator->GetMinLogLikelihoodTimesPrior();
  if ( std::isnan( minLogLikelihoodTimesPrior ) || std::isinf( minLogLikelihoodTimesPrior ) )
  {
    minLogLikelihoodTimesPrior = itk::NumericTraits< float >::max();
  }

  return gradientCalculator->GetPositionGradient();

}



//
//
//
AtlasPositionGradientContainerType::Pointer
AtlasParameterEstimator
::GetCurrentPositionGradient( unsigned int labelImageNumber ) const
{
  float  dummy;

  return this->CalculateCurrentPositionGradient( labelImageNumber, dummy );
}



//
//
//
float
AtlasParameterEstimator
::GetCurrentMinLogLikelihood() const
{

  return this->GetMinLogLikelihood( m_MeshCollection );

}



//
//
//
float
AtlasParameterEstimator
::GetMinLogLikelihood( const AtlasMeshCollection* meshCollection ) const
{
  float  minLogLikelihood = 0.0f;

  AtlasMeshMinLogLikelihoodCalculator::Pointer  calculator = AtlasMeshMinLogLikelihoodCalculator::New();
  for ( unsigned int  labelImageNumber = 0; labelImageNumber < m_NumberOfLabelImages; labelImageNumber++ )
  {
    // Calculate minLogLikelihood for this label image
    calculator->SetLabelImage( m_LabelImages[ labelImageNumber ] );
    calculator->Rasterize( meshCollection->GetMesh( labelImageNumber ) );
    minLogLikelihood += calculator->GetMinLogLikelihood();
  }

  return minLogLikelihood;

}



} // end namespace kvl

