/**
 * @file  kvlAtlasMeshSegmenter.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
 *    $Revision: 1.5 $
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
#include "kvlAtlasMeshSegmenter.h"

#include "kvlAtlasMeshToIntensityImageGradientCalculator.h"
#include "itkCastImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkIntensityWindowingImageFilter.h"
#include "kvlAtlasMeshVisitCounter.h"

#include "kvlAtlasMeshToProbabilityImageGradientCalculator.h"
#include "kvlAtlasMeshDeformationLevenbergMarquardt.h"

#if 1
#include "kvlAtlasMeshCollection.h"
#endif

namespace kvl
{

//
//
//
AtlasMeshSegmenter
::AtlasMeshSegmenter()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_PartialVolumeUpsamplingFactors[ i ] = 1;
  }
  m_UsePartialVoluming = false;

  m_Image = 0;
  m_InternalImage = 0;
  m_Mesh = 0;
  m_EMSegmenter = EMSegmenter::New();

  m_PositionUpdatingIterationNumber = 0;
  m_PositionUpdatingMaximumNumberOfIterations = 11 /* 2000 */ /* itk::NumericTraits< unsigned int >::max() */;
  m_PositionUpdatingIterationEventResolution = 10;

  m_IterationNumber = 0;
  m_MaximumNumberOfIterations = 30;
  m_CoregisterToPosteriorProbabilities = false;
  m_StopCriterion = 1e-5;
  m_MaximalDeformationStopCriterion = 0.05;

  m_MeshToImageTransform = TransformType::New();

  m_PartialVolumeGradientCalculator = 0;
  m_PosteriorProbabilityImage = 0;

  m_DontDeform = false;

}




//
//
//
AtlasMeshSegmenter
::~AtlasMeshSegmenter()
{

}



//
//
//
void
AtlasMeshSegmenter
::Segment( bool useAffine )
{
  // Sanity check on input
  if ( !m_Image || !m_Mesh )
  {
    itkExceptionMacro( << "Insufficient inputs set!" );
  }

  // Now loop over the two steps
  this->InvokeEvent( itk::StartEvent() );

  float  currentCost = itk::NumericTraits< float >::max();
  for ( m_IterationNumber = 0; m_IterationNumber < m_MaximumNumberOfIterations; m_IterationNumber++ )
  {
    // Update the model parameters
    const float  previousCost = currentCost;
    currentCost = this->UpdateModelParameters();
    currentCost += this->CalculateCurrentWarpCost();

    // Evaluate convergence.
    std::cout << "============================================" << std::endl;
    std::cout << "Cost for iteration " << m_IterationNumber << " was " << currentCost << std::endl;
    std::cout << "============================================" << std::endl;

    if ( ( ( ( previousCost - currentCost ) / fabsf( currentCost ) ) < m_StopCriterion ) || m_DontDeform )
    {
      std::cout << "Converged!" << std::endl;
      break;
    }

    // Make sure to invalidate the posterior image, and update the position
    m_PosteriorProbabilityImage = 0;
    if ( !this->UpdatePosition( useAffine ) )
    {
      std::cout << "Converged!" << std::endl;
      break;
    }


    this->InvokeEvent( itk::IterationEvent() );
  }

  // Retrieve the posterior
  m_Posteriors.clear();
  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
  {
    m_Posteriors.push_back( m_EMSegmenter->GetPosterior( classNumber ) );
  }

  this->InvokeEvent( itk::EndEvent() );
}



//
//
//
void
AtlasMeshSegmenter
::SetPartialVolumeUpsamplingFactors( const int* partialVolumeUpsamplingFactors )
{

  m_UsePartialVoluming = false;
  for ( int i = 0; i < 3; i++ )
  {
    m_PartialVolumeUpsamplingFactors[ i ] = partialVolumeUpsamplingFactors[ i ];

    if ( m_PartialVolumeUpsamplingFactors[ i ] > 1 )
    {
      m_UsePartialVoluming = true;
    }
  }

  m_EMSegmenter->SetPartialVolumeUpsamplingFactors( m_PartialVolumeUpsamplingFactors );

}


//
//
//
void
AtlasMeshSegmenter
::SetImage( const ImageType*  image )
{

  // Calculate internal image
  std::cout << "TODO: we're using internally an image of unsigned char pixel type"
            " to calculate the position gradients" << std::endl;


  // Calculate minimum and maximum
  typedef itk::MinimumMaximumImageCalculator< ImageType >  RangeCalculatorType;
  RangeCalculatorType::Pointer  rangeCalculator = RangeCalculatorType::New();
  rangeCalculator->SetImage( image );
  rangeCalculator->Compute();
  std::cout << "minimum: " << rangeCalculator->GetMinimum() << std::endl;
  std::cout << "maximum: " << rangeCalculator->GetMaximum() << std::endl;


  // Scale and clip intensities to be between 0 and 255
  typedef itk::IntensityWindowingImageFilter< ImageType, InternalImageType >   WindowerType;
  WindowerType::Pointer  windower = WindowerType::New();
  windower->SetInput( image );
  windower->SetWindowMinimum( rangeCalculator->GetMinimum() );
  windower->SetWindowMaximum( rangeCalculator->GetMaximum() );
  windower->SetOutputMinimum( 0 );
  windower->SetOutputMaximum( 255 );
  windower->Update();
  m_InternalImage = windower->GetOutput();



#if 0
  // Cast internal image back to original pixel type
  typedef itk::CastImageFilter< InternalImageType, ImageType >  BackCasterType;
  BackCasterType::Pointer  backCaster = BackCasterType::New();
  backCaster->SetInput( m_InternalImage );
  backCaster->Update();
  m_Image = backCaster->GetOutput();
#else
  m_Image = image;
#endif

  m_EMSegmenter->SetImage( m_Image.GetPointer() );

}



//
//
//
void
AtlasMeshSegmenter
::SetMesh( AtlasMesh* mesh,
           const std::vector< unsigned int >&  lookupTable,
           const std::vector< unsigned int >&  independentParametersLookupTable )
{

  m_EMSegmenter->SetAtlasMesh( mesh, lookupTable, independentParametersLookupTable );

  m_NumberOfClasses = m_EMSegmenter->GetNumberOfClasses();
  m_Mesh = m_EMSegmenter->GetAtlasMesh();

  std::cout << "OK" << std::endl;
}





//
//
//
bool
AtlasMeshSegmenter
::UpdatePosition( bool useAffine )
{

  this->InvokeEvent( PositionUpdatingStartEvent() );

  std::cout << "m_PositionUpdatingMaximumNumberOfIterations: " << m_PositionUpdatingMaximumNumberOfIterations << std::endl;

  // Keep track if we've moved at all
  bool  haveMoved = false;

  // Convert means and variances to correct format
  itk::Array< float >   means( m_NumberOfClasses );
  itk::Array< float >   variances( m_NumberOfClasses );
  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
  {
    means[ classNumber ] = m_Means[ classNumber ];
    variances[ classNumber ] = m_Variances[ classNumber ];
  }

  // Some initialization
  AtlasMesh::Pointer  meshToRasterize = AtlasMesh::New();
  meshToRasterize->SetPoints( m_Mesh->GetPoints() );
  meshToRasterize->SetCells( m_Mesh->GetCells() );
  meshToRasterize->SetPointData( m_Mesh->GetPointData() );
  meshToRasterize->SetCellData( m_Mesh->GetCellData() );
  meshToRasterize = m_EMSegmenter->GetSuperResolutionMesh( meshToRasterize );

#if 0
  {
    AtlasMeshCollection::Pointer  collection = AtlasMeshCollection::New();
    collection->GenerateFromSingleMesh( meshToRasterize, 1, 1000.0f );
    collection->Write( "meshToRasterize" );
    //exit( -1 );
  }
#endif

  bool  stopping = false;
  float  lambda = 1.0f; // Marquardt's lambda
  const float  lambdaIncreaseFactor = 2.0f; // Factor by which lambda is increased if failure detected
  const float  lambdaDecreaseFactor = 1.1f; // Factor by which lambda is decreased if success detected
  AtlasMeshDeformationLevenbergMarquardt::Pointer  levenbergMarquardt =
    AtlasMeshDeformationLevenbergMarquardt::New();
  InternalImageType::ConstPointer  dummyTemplateImage = m_InternalImage;
  if ( m_UsePartialVoluming )
  {
    typedef itk::CastImageFilter< EMSegmenter::ImageType, InternalImageType >  CasterType;
    CasterType::Pointer  caster = CasterType::New();
    caster->SetInput( m_EMSegmenter->GetSuperResolutionImage() );
    caster->Update();
    dummyTemplateImage = caster->GetOutput();
  }
  levenbergMarquardt->SetLabelImage( dummyTemplateImage );

  if ( m_CoregisterToPosteriorProbabilities )
  {
    levenbergMarquardt->SetProbabilityImage( this->GetPosteriorProbabilityImage() );
    levenbergMarquardt->SetUseProbabilityImage( true );
  }
  else
  {
    if ( !m_UsePartialVoluming )
    {
      levenbergMarquardt->SetImage( m_EMSegmenter->GetBiasCorrectedImage() );
      levenbergMarquardt->SetMeans( means );
      levenbergMarquardt->SetVariances( variances );
    }
    else
    {
      //levenbergMarquardt->SetLikelihoods( m_EMSegmenter->GetSuperResolutionLikelihoods() );
      itkExceptionMacro( "That's not implemented yet!" );
    }
    levenbergMarquardt->SetUseProbabilityImage( false );
  }

  //
  TransformType::Pointer  meshToImageTransform = TransformType::New();
  meshToImageTransform->SetParameters( m_MeshToImageTransform->GetParameters() );
  std::cout << "m_MeshToImageTransform: " << std::endl;
  m_MeshToImageTransform->Print( std::cout );
  if ( m_UsePartialVoluming )
  {
    TransformType::OutputVectorType  extraScaling;
    for ( int i = 0; i < 3; i++ )
    {
      extraScaling[ i ] = m_PartialVolumeUpsamplingFactors[ i ];
    }
    meshToImageTransform->Scale( extraScaling );
  }
  std::cout << "meshToImageTransform: " << std::endl;
  meshToImageTransform->Print( std::cout );
  levenbergMarquardt->SetMeshToImageTransform( meshToImageTransform );


  // Evaluate cost at current position
  std::cout << "Rasterizing Levenberg-Marquardt..." << std::endl;
  levenbergMarquardt->Rasterize( meshToRasterize );
  std::cout << "...done!" << std::endl;
  double  currentCost = levenbergMarquardt->GetMinLogLikelihoodTimesPrior();


  for ( m_PositionUpdatingIterationNumber = 0;
        m_PositionUpdatingIterationNumber < m_PositionUpdatingMaximumNumberOfIterations;
        m_PositionUpdatingIterationNumber++ )
  {

    //
    while ( true )
    {
      // Save current position
      AtlasMesh::PointsContainer::Pointer  currentPosition = m_EMSegmenter->GetSuperResolutionMesh( m_Mesh )->GetPoints();


      // Calculate the trial step with the current lambda
      std::cout << "Calculating Levenberg-Marquardt step with lambda " << lambda << "..." << std::endl;
      AtlasPositionGradientContainerType::Pointer  trialStep = levenbergMarquardt->GetStep( lambda, false );
      std::cout << "...done!" << std::endl;


      // Add the step to the current position to obtain the trial position
      float  maximalDeformation = 0.0f;
      AtlasMesh::PointsContainer::Pointer  trialPosition = AtlasMesh::PointsContainer::New();
      AtlasMesh::PointsContainer::ConstIterator  posIt = currentPosition->Begin();
      AtlasPositionGradientContainerType::ConstIterator  stepIt = trialStep->Begin();
      while ( stepIt != trialStep->End() )
      {
        AtlasPositionGradientType  trialStep = stepIt.Value();

        if ( trialStep.GetNorm() > maximalDeformation )
        {
          maximalDeformation = trialStep.GetNorm();
        }

        trialPosition->InsertElement( posIt.Index(), posIt.Value() + trialStep );

        // Move on to next vertex
        ++posIt;
        ++stepIt;
      }

      std::cout << "maximalDeformation: " << maximalDeformation << std::endl;


      // Evaluate the trial position
      AtlasMeshDeformationLevenbergMarquardt::Pointer  trialLevenbergMarquardt =
        AtlasMeshDeformationLevenbergMarquardt::New();
      trialLevenbergMarquardt->SetLabelImage( dummyTemplateImage );
      if ( m_CoregisterToPosteriorProbabilities )
      {
        trialLevenbergMarquardt->SetProbabilityImage( this->GetPosteriorProbabilityImage() );
        trialLevenbergMarquardt->SetUseProbabilityImage( true );
      }
      else
      {
        if ( !m_UsePartialVoluming )
        {
          trialLevenbergMarquardt->SetImage( m_EMSegmenter->GetBiasCorrectedImage() );
          trialLevenbergMarquardt->SetMeans( means );
          trialLevenbergMarquardt->SetVariances( variances );
        }
        else
        {
          //levenbergMarquardt->SetLikelihoods( m_EMSegmenter->GetSuperResolutionLikelihoods() );
          itkExceptionMacro( "That's not implemented yet!" );
        }
        trialLevenbergMarquardt->SetUseProbabilityImage( false );

      }
      trialLevenbergMarquardt->SetMeshToImageTransform( meshToImageTransform );

      meshToRasterize->SetPoints( trialPosition );
      m_Mesh->SetPoints( m_EMSegmenter->GetNormalResolutionMesh( meshToRasterize )->GetPoints() );
      std::cout << "Rasterizing Levenberg-Marquardt..." << std::endl;
      trialLevenbergMarquardt->Rasterize( meshToRasterize );
      std::cout << "...done!" << std::endl;
      double  trialCost = trialLevenbergMarquardt->GetMinLogLikelihoodTimesPrior();


      // Do the Right Thing
      std::cout << "=======================" << std::endl;
      if ( trialCost < currentCost )
      {
        // Good move. Keep the move and decrease Marquardt's lambda
        std::cout << "Good trial move (maximal deformation: " << maximalDeformation
                  << ", lambda: " << lambda << ")" << std::endl;
        std::cout << "        Cost going from " << currentCost << " to " << trialCost << std::endl;
        std::cout << "        Decreasing lambda from " << lambda;
        lambda /= lambdaDecreaseFactor;
        std::cout << " to " << lambda << std::endl;

        currentCost = trialCost;
        levenbergMarquardt = trialLevenbergMarquardt;

        haveMoved = true;

        // This was a good step, but if we're not making much progress anymore, simply stop
        if ( maximalDeformation < m_MaximalDeformationStopCriterion )
        {
          std::cout << "        maximalDeformation is too small; stopping" << std::endl;
          stopping = true;
        }

        break;
      }
      else
      {
        // Bad move. Restore previous position
        std::cout << "Bad trial move (maximal deformation: " << maximalDeformation
                  << ", lambda: " << lambda << ")" << std::endl;
        std::cout << "        Cost going from " << currentCost << " to " << trialCost << std::endl;

        meshToRasterize->SetPoints( currentPosition );
        m_Mesh->SetPoints( m_EMSegmenter->GetNormalResolutionMesh( meshToRasterize )->GetPoints() );

        // Increase lambda and try again - unless the proposed deformation was already so small you want to stop altogether
        if ( maximalDeformation > 0.01 /* 0.01 */ )
        {
          std::cout << "        Increasing lambda from " << lambda;
          lambda *= lambdaIncreaseFactor;
          std::cout << " to " << lambda << " and trying again." << std::endl;
        }
        else
        {
          std::cout << "        proposed deformation is already too small; bailing out" << std::endl;
          stopping = true;
          break;
        }

      } // End test if good or bad trial move


    } // End trying loop until we find a good step or we're converged

    if ( stopping )
    {
      break;
    }

    if ( !( m_PositionUpdatingIterationNumber % m_PositionUpdatingIterationEventResolution ) )
    {
      this->InvokeEvent( PositionUpdatingIterationEvent() );
    }

  } // End loop over all positive moves


  this->InvokeEvent( PositionUpdatingEndEvent() );


  return haveMoved;

}




//
//
//
float
AtlasMeshSegmenter
::UpdateModelParameters()
{
  const float  cost = m_EMSegmenter->Segment();

  m_Means = m_EMSegmenter->GetMeans();
  m_Variances = m_EMSegmenter->GetVariances();

  return cost;
}



//
//
//
void
AtlasMeshSegmenter
::MakeAffine( AtlasPositionGradientContainerType* gradient ) const
{
  //
  std::cout << "====== Making gradient affine =========== " << std::endl;


  // The parameterization of the affine transform is simply the deformation (displacement) of
  // the vectors o = (0, 0, 0)^T, e1 = (1, 0, 0)^T, e2 = (0,1,0)^T, and e3 = (0,0,1)^T. In order
  // to calculate the gradient w.r.t. this parameterization, we can simply apply the chain rule,
  // i.e. by adding all contributions of the mesh nodes to the total gradient
  itk::Vector< float, 3 >  gradientOfE1( 0.0f );
  itk::Vector< float, 3 >  gradientOfE2( 0.0f );
  itk::Vector< float, 3 >  gradientOfE3( 0.0f );
  itk::Vector< float, 3 >  gradientOfOrigin( 0.0f );
  AtlasMesh::PointsContainer::ConstIterator  posIt = m_Mesh->GetPoints()->Begin();
  AtlasPositionGradientContainerType::ConstIterator  gradIt = gradient->Begin();
  for( ; posIt != m_Mesh->GetPoints()->End(); ++posIt, ++gradIt )
  {
    gradientOfE1 += posIt.Value()[ 0 ] * gradIt.Value();
    gradientOfE2 += posIt.Value()[ 1 ] * gradIt.Value();
    gradientOfE3 += posIt.Value()[ 2 ] * gradIt.Value();
    gradientOfOrigin += gradIt.Value();
  }

  std::cout << "   gradientOfE1: " << gradientOfE1 << std::endl;
  std::cout << "   gradientOfE2: " << gradientOfE2 << std::endl;
  std::cout << "   gradientOfE3: " << gradientOfE3 << std::endl;

  // Now map the gradient in the parameterization back onto a gradient in the mesh nodes
  AtlasPositionGradientContainerType::Iterator  grad2It = gradient->Begin();
  for( posIt = m_Mesh->GetPoints()->Begin();
       posIt != m_Mesh->GetPoints()->End(); ++posIt, ++grad2It )
  {
    grad2It.Value() = gradientOfE1 * posIt.Value()[ 0 ] +
                      gradientOfE2 * posIt.Value()[ 1 ] +
                      gradientOfE3 * posIt.Value()[ 2 ] +
                      gradientOfOrigin;
  }



}




//
//
//
AtlasMeshSegmenter::SummaryImageType::Pointer
AtlasMeshSegmenter
::GetSummaryImage( bool  usePrior, bool  useCrisp, const std::vector< float >&  means ) const
{

  // Create an empty image to hold the reconstruction
  SummaryImageType::Pointer  summary = SummaryImageType::New();
  summary->SetRegions( m_EMSegmenter->GetPosterior( 0 )->GetLargestPossibleRegion() );
  summary->SetSpacing( m_EMSegmenter->GetPosterior( 0 )->GetSpacing() );
  summary->SetOrigin( m_EMSegmenter->GetPosterior( 0 )->GetOrigin() );
  summary->Allocate();
  summary->FillBuffer( 0 );

  // Create an empty image to hold the maximum probability in each voxel
  SummaryImageType::Pointer  maximumProbability = SummaryImageType::New();
  maximumProbability->SetRegions( m_EMSegmenter->GetPosterior( 0 )->GetLargestPossibleRegion() );
  maximumProbability->SetSpacing( m_EMSegmenter->GetPosterior( 0 )->GetSpacing() );
  maximumProbability->SetOrigin( m_EMSegmenter->GetPosterior( 0 )->GetOrigin() );
  maximumProbability->Allocate();
  maximumProbability->FillBuffer( 0 );

  // Loop over all the posterior images
  for ( unsigned int  classNumber = 0; classNumber < this->GetNumberOfClasses(); classNumber++ )
  {
    // Decide what mean to use
    float  mean;
    if ( means.size() != 0 )
    {
      mean = means[ classNumber ];
    }
    else
    {
      //mean = classNumber;
      mean = m_EMSegmenter->GetLookupTable()[ classNumber ];
    }

    // Loop over all voxels
    EMSegmenter::ClassificationImageType::ConstPointer  classification = 0;
    if ( usePrior )
    {
      classification =  m_EMSegmenter->GetPrior( classNumber );
    }
    else
    {
      classification =  m_EMSegmenter->GetPosterior( classNumber );
    }
    itk::ImageRegionConstIterator< EMSegmenter::ClassificationImageType >
    classificationIt( classification,
                      classification->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< SummaryImageType >  maximumIt( maximumProbability,
        maximumProbability->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< SummaryImageType >  summaryIt( summary,
        summary->GetLargestPossibleRegion() );

    for ( ; !classificationIt.IsAtEnd(); ++classificationIt, ++maximumIt, ++summaryIt )
    {
      if ( useCrisp )
      {
        // Only override voxel value if this is the MAP so far
        if ( classificationIt.Value() > maximumIt.Value() )
        {
          maximumIt.Value() = classificationIt.Value();
          summaryIt.Value() = mean;
        }

      }
      else
      {
        summaryIt.Value() += mean * classificationIt.Value();
      }


    } // End loop over all voxels

  } // End loop over all labels


  return summary;


}



//
//
//
AtlasMeshSegmenter::ProbabilityImageType::Pointer
AtlasMeshSegmenter
::GetPosteriorProbabilityImage() const
{

  if ( !m_PosteriorProbabilityImage )
  {
    if ( !m_UsePartialVoluming )
    {
      // Combine the posteriors, which are now contained in different images,
      // into one image
      std::cout << "Constructing posterior probability image...";
      m_PosteriorProbabilityImage = ProbabilityImageType::New();
      m_PosteriorProbabilityImage->SetRegions( m_EMSegmenter->GetPosterior( 0 )->GetBufferedRegion() );
      m_PosteriorProbabilityImage->Allocate();

      AtlasAlphasType  emptyAlphas( m_NumberOfClasses );
      emptyAlphas.Fill( 0 );
      m_PosteriorProbabilityImage->FillBuffer( emptyAlphas );

      for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
      {
        // Loop over all voxels
        itk::ImageRegionConstIterator< ClassificationImageType >  sourceIt( m_EMSegmenter->GetPosterior( classNumber ),
            m_EMSegmenter->GetPosterior( classNumber )->GetBufferedRegion() );
        itk::ImageRegionIterator< ProbabilityImageType >  targetIt( m_PosteriorProbabilityImage,
            m_PosteriorProbabilityImage->GetBufferedRegion() );
        for ( ; !sourceIt.IsAtEnd(); ++sourceIt, ++targetIt )
        {
          targetIt.Value()[ classNumber ] = sourceIt.Value();
        }

      }

      std::cout << " done!" << std::endl;
    }
    else
    {
      std::cout << "AtlasMeshSegmenter: retrieving super resolution posterior..." << std::endl;
      m_PosteriorProbabilityImage = m_EMSegmenter->GetSuperResolutionPosteriors();
      std::cout << "...done!" << std::endl;
    }

  } // End test if m_PosteriorProbabilityImage already exists

  return m_PosteriorProbabilityImage;

}



//
//
//
float
AtlasMeshSegmenter
::CalculateCurrentWarpCost() const
{

  // Rasterize image with only zeros using gradient calculator: since
  // zeros are skipped, the resulting cost will be the warp cost only :-)
  InternalImageType::Pointer  zeroImage = InternalImageType::New();
  zeroImage->SetRegions( m_InternalImage->GetLargestPossibleRegion() );
  zeroImage->Allocate();
  zeroImage->FillBuffer( 0 );


  // Convert means and variances to correct format
  itk::Array< float >   means( m_NumberOfClasses );
  itk::Array< float >   variances( m_NumberOfClasses );
  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
  {
    means[ classNumber ] = m_Means[ classNumber ];
    variances[ classNumber ] = m_Variances[ classNumber ];
  }

  // Rasterize using gradient calculator, just so we can evaluate the current cost
  AtlasMeshToIntensityImageGradientCalculator::Pointer  gradientCalculator =
    AtlasMeshToIntensityImageGradientCalculator::New();
  gradientCalculator->SetLabelImage( zeroImage );
  gradientCalculator->SetMeans( means );
  gradientCalculator->SetVariances( variances );
  //gradientCalculator->SetIgnoreDeformationPrior( useAffine );
  gradientCalculator->Rasterize( m_Mesh );

  return gradientCalculator->GetMinLogLikelihoodTimesPrior();

}


} // end namespace kvl
