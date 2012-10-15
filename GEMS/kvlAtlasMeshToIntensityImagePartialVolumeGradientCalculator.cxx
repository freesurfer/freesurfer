/**
 * @file  kvlAtlasMeshToIntensityImagePartialVolumeGradientCalculator.cxx
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
#include "kvlAtlasMeshToIntensityImagePartialVolumeGradientCalculator.h"


#include "kvlAtlasMeshMultiAlphaDrawer.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#if 0
#include "itkImageFileWriter.h"
#endif



namespace kvl
{

//
//
//
AtlasMeshToIntensityImagePartialVolumeGradientCalculator
::AtlasMeshToIntensityImagePartialVolumeGradientCalculator()
{
  int  partialVolumeUpsamplingFactors[ 3 ] = { 2, 2, 2 };
  this->SetPartialVolumeUpsamplingFactors( partialVolumeUpsamplingFactors );

  m_MinLogLikelihood = 0;

  m_RasterizedHighResolutionPriors = 0;
  m_SiblingPriorsImage = 0;
  m_SubvoxelPriorsImage = 0;
  m_LikelihoodImage = 0;

  m_GradientBasisWeights = 0;
  m_LikelihoodPointerImage = 0;

  m_HighResolutionImage = 0;
  m_Image = 0;

}


//
//
//
AtlasMeshToIntensityImagePartialVolumeGradientCalculator
::~AtlasMeshToIntensityImagePartialVolumeGradientCalculator()
{

}



//
//
//
void
AtlasMeshToIntensityImagePartialVolumeGradientCalculator
::SetLabelImage( const LabelImageType*  labelImage )
{

  if ( ( m_PartialVolumeUpsamplingFactors[ 0 ] != m_PartialVolumeUpsamplingFactors[ 1 ] ) ||
       ( m_PartialVolumeUpsamplingFactors[ 0 ] != m_PartialVolumeUpsamplingFactors[ 2 ] ) ||
       ( m_PartialVolumeUpsamplingFactors[ 1 ] != m_PartialVolumeUpsamplingFactors[ 2 ] ) )
  {
    itkExceptionMacro( "Anisotropic partial voluming isn't implemented: you should scale the different components of the gradient accordingly!!!" );
  }

  m_Image = labelImage;


  // Overall picture:
  //
  // The rasterizor needs to visit each of the subvoxels. In each subvoxel, it needs to know
  // the prior probability of the other subvoxels belonging to the same big "averaging" voxel,
  // and the intensity of that big voxel.
  //
  //
  // Practical implementation:
  //
  // We're simply gonna upsample the original image by repeating the intensity of the big voxel
  // to the subvoxels; this will allow easy access during the rasterization to the original big
  // voxel intensity for the subvoxels.
  // Also, we're gonna construct a high-res image that has a vector value in each subvoxel, holding
  // the prior probabilities of the remaining subvoxels. The content of this image will need to
  // be filled in prior to rasterizing a given mesh. Also, make sure that subvoxels that lie outside
  // of the mesh are always filled in as background
  //
  //

  // Create an empty high-resolution label image
  LabelImageType::SizeType  highResolutionSize;
  for ( int i = 0; i < 3; i++ )
  {
    highResolutionSize[ i ] = labelImage->GetBufferedRegion().GetSize()[ i ] * m_PartialVolumeUpsamplingFactors[ i ];
  }
  m_HighResolutionImage = LabelImageType::New();
  m_HighResolutionImage->SetRegions( highResolutionSize );
  m_HighResolutionImage->Allocate();
  m_HighResolutionImage->FillBuffer( 0 );

  // Fill in the partial volumed data
  for ( itk::ImageRegionConstIteratorWithIndex< LabelImageType >  it( labelImage,
        labelImage->GetBufferedRegion() );
        !it.IsAtEnd(); ++it )
  {
    //std::cout << "Filling in label in subvoxels of big voxel with index " << it.GetIndex() << std::endl;

    LabelImageType::IndexType  index;
    for ( int xStep = 0; xStep < m_PartialVolumeUpsamplingFactors[ 0 ]; xStep++ )
    {
      index[ 0 ] = it.GetIndex()[ 0 ] * m_PartialVolumeUpsamplingFactors[ 0 ] + xStep;
      for ( int yStep = 0; yStep < m_PartialVolumeUpsamplingFactors[ 1 ]; yStep++ )
      {
        index[ 1 ] = it.GetIndex()[ 1 ] * m_PartialVolumeUpsamplingFactors[ 1 ] + yStep;
        for ( int zStep = 0; zStep < m_PartialVolumeUpsamplingFactors[ 2 ]; zStep++ )
        {
          index[ 2 ] = it.GetIndex()[ 2 ] * m_PartialVolumeUpsamplingFactors[ 2 ] + zStep;

          //std::cout << "         Copying label to high-res voxel with index " << index << std::endl;
          m_HighResolutionImage->SetPixel( index, it.Value() );
        }
      }
    }

  }

#if 0
  // Write out the partial volumed image
  typedef itk::ImageFileWriter< LabelImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetFileName( "debugHighResolutionImage.mhd" );
  writer->SetInput( m_HighResolutionImage );
  writer->Update();
#endif


  // Allocate a high-resolution image that will hold the rasterized priors
  m_RasterizedHighResolutionPriors = RasterizedHighResolutionPriorsType::New();
  m_RasterizedHighResolutionPriors->SetRegions( highResolutionSize );
  m_RasterizedHighResolutionPriors->Allocate();


  // Allocate a high-resolution image that will hold the weights, and pass it
  // to the rasterizor
  m_GradientBasisWeights = GradientBasisWeightsType::New();
  m_GradientBasisWeights->SetRegions( highResolutionSize );
  m_GradientBasisWeights->Allocate();
  this->GetFragmentProcessor().SetGradientBasisWeights( m_GradientBasisWeights );


  // Allocate a high-resolution image that will hold a pointer of each subvoxel's big
  // "averaging" voxel likelihood, and pass it to the rasterizor
  m_LikelihoodPointerImage = LikelihoodPointerImageType::New();
  m_LikelihoodPointerImage->SetRegions( highResolutionSize );
  m_LikelihoodPointerImage->Allocate();
  this->GetFragmentProcessor().SetLikelihoodPointerImage( m_LikelihoodPointerImage );


  // Allocate a high-resolution image that, for each subvoxel, points to
  // the high-resolution prior probabilities in the subvoxels sharing the
  // same big, "averaging" voxel. At the same time, construct a low-resolution
  // image that holds pointers to the prior probabilities in each of its voxel's
  // subvoxels, and a high-resolution image to hold a pointer to the low-resolution
  // likelihood of its parent "averaging" voxel.
  m_SiblingPriorsImage = SiblingPriorsImageType::New();
  m_SiblingPriorsImage->SetRegions( highResolutionSize );
  m_SiblingPriorsImage->Allocate();

  m_SubvoxelPriorsImage = SubvoxelPriorsImageType::New();
  m_SubvoxelPriorsImage->SetRegions( labelImage->GetBufferedRegion().GetSize() );
  m_SubvoxelPriorsImage->Allocate();

  m_LikelihoodImage = LikelihoodImageType::New();
  m_LikelihoodImage->SetRegions( labelImage->GetBufferedRegion().GetSize() );
  m_LikelihoodImage->Allocate();

  itk::ImageRegionIteratorWithIndex< SubvoxelPriorsImageType >  it( m_SubvoxelPriorsImage,
      m_SubvoxelPriorsImage->GetBufferedRegion() );
  itk::ImageRegionConstIterator< LikelihoodImageType >  likelihoodIt( m_LikelihoodImage,
      m_LikelihoodImage->GetBufferedRegion() );
  for (; !it.IsAtEnd(); ++it, ++likelihoodIt )
  {
    //std::cout << "Filling in pointers to sibling priors in big voxel with index " << it.GetIndex() << std::endl;

    // Collect *all* pointers in this voxel's subvoxels
    std::vector< AtlasAlphasType* >  priorsOfSubvoxels( m_NumberOfSubvoxels );
    LabelImageType::IndexType  index;
    int  counter = -1;
    for ( int xStep = 0; xStep < m_PartialVolumeUpsamplingFactors[ 0 ]; xStep++ )
    {
      index[ 0 ] = it.GetIndex()[ 0 ] * m_PartialVolumeUpsamplingFactors[ 0 ] + xStep;
      for ( int yStep = 0; yStep < m_PartialVolumeUpsamplingFactors[ 1 ]; yStep++ )
      {
        index[ 1 ] = it.GetIndex()[ 1 ] * m_PartialVolumeUpsamplingFactors[ 1 ] + yStep;
        for ( int zStep = 0; zStep < m_PartialVolumeUpsamplingFactors[ 2 ]; zStep++ )
        {
          index[ 2 ] = it.GetIndex()[ 2 ] * m_PartialVolumeUpsamplingFactors[ 2 ] + zStep;
          counter++;

          priorsOfSubvoxels[ counter ] = &( m_RasterizedHighResolutionPriors->GetPixel( index ) );
          m_LikelihoodPointerImage->SetPixel( index, &( likelihoodIt.Value() ) );

#if 0
          {
            if ( ( index[ 0 ] == 3 ) && ( index[ 1 ] == 16 ) && ( index[ 2 ] == 16 ) )
            {
              std::cout << "   m_LikelihoodPointerImage voxel with index " << index << " made to point to address " <<  &( likelihoodIt.Value() ) << std::endl;
              std::cout << "   Correspoding low-res index is " << likelihoodIt.GetIndex() << std::endl;
              std::cout << "   Enter character to continue" << std::endl;
              char  dummy;
              std::cin >> dummy;
            }


          }
#endif

        }
      }
    }

    // Assign it in the low-res image m_SubvoxelPriorsImage
    it.Value() = priorsOfSubvoxels;

    // Now simply copy all pointers to each subvoxel, except don't copy the one affilidated
    // with itself
    counter = -1;
    for ( int xStep = 0; xStep < m_PartialVolumeUpsamplingFactors[ 0 ]; xStep++ )
    {
      index[ 0 ] = it.GetIndex()[ 0 ] * m_PartialVolumeUpsamplingFactors[ 0 ] + xStep;
      for ( int yStep = 0; yStep < m_PartialVolumeUpsamplingFactors[ 1 ]; yStep++ )
      {
        index[ 1 ] = it.GetIndex()[ 1 ] * m_PartialVolumeUpsamplingFactors[ 1 ] + yStep;
        for ( int zStep = 0; zStep < m_PartialVolumeUpsamplingFactors[ 2 ]; zStep++ )
        {
          index[ 2 ] = it.GetIndex()[ 2 ] * m_PartialVolumeUpsamplingFactors[ 2 ] + zStep;
          counter++;

          std::vector< AtlasAlphasType* >  siblingPriors( m_NumberOfSiblings );
          int i = 0;
          int j = 0;
          for ( ; i < m_NumberOfSubvoxels; i++ )
          {
            if ( i != counter )
            {
              siblingPriors[ j ] = priorsOfSubvoxels[ i ];
              j++;
            }

          }
          m_SiblingPriorsImage->SetPixel( index, siblingPriors );

        }
      }
    }

  } // End loop over all big, "averaging" voxels




  // Invoke superclass' implementation
  Superclass::SetLabelImage( m_HighResolutionImage );

}




//
//
//
void
AtlasMeshToIntensityImagePartialVolumeGradientCalculator
::Rasterize( const AtlasMesh* mesh )
{

  // Construct a mesh that lives in the high-res space
  float  translation[ 3 ];
  float  scaling[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    translation[ i ] = ( m_PartialVolumeUpsamplingFactors[ i ] - 1 ) / 2.0f;
    scaling[ i ] = m_PartialVolumeUpsamplingFactors[ i ];
  }

  AtlasMesh::PointsContainer::Pointer  highResolutionPoints = AtlasMesh::PointsContainer::New();
  for ( AtlasMesh::PointsContainer::ConstIterator it = mesh->GetPoints()->Begin();
        it != mesh->GetPoints()->End(); ++it )
  {
    AtlasMesh::PointType  highResolutionPoint;
    for ( int i = 0; i < 3; i++ )
    {
      highResolutionPoint[ i ] = scaling[ i ] * it.Value()[ i ] + translation[ i ];
    }

    highResolutionPoints->InsertElement( it.Index(), highResolutionPoint );
  }

  AtlasMesh::Pointer  highResolutionMesh = AtlasMesh::New();
  highResolutionMesh->SetPoints( highResolutionPoints );
  highResolutionMesh->SetCells( const_cast< AtlasMesh::CellsContainer* >( mesh->GetCells() ) );
  highResolutionMesh->SetPointData( const_cast< AtlasMesh::PointDataContainer* >( mesh->GetPointData() ) );
  highResolutionMesh->SetCellData( const_cast< AtlasMesh::CellDataContainer* >( mesh->GetCellData() ) );



  // Rasterize the mesh into the prior probabilities in each high-resolution voxel
  AtlasMeshMultiAlphaDrawer::Pointer  drawer = AtlasMeshMultiAlphaDrawer::New();
  drawer->SetLabelImage( this->GetLabelImage() );
  drawer->SetAlphasImage( m_RasterizedHighResolutionPriors );
  drawer->Rasterize( highResolutionMesh );


#if 0
  {
    typedef itk::Image< float, 3 >  DebugImageType;
    DebugImageType::Pointer  debugImage = DebugImageType::New();
    debugImage->SetRegions( m_RasterizedHighResolutionPriors->GetBufferedRegion() );
    debugImage->Allocate();

    itk::ImageRegionConstIterator< RasterizedHighResolutionPriorsType >  sourceIt( m_RasterizedHighResolutionPriors,
        m_RasterizedHighResolutionPriors->GetBufferedRegion() );
    itk::ImageRegionIterator< DebugImageType >  targetIt( debugImage, debugImage->GetBufferedRegion() );
    for ( ; !sourceIt.IsAtEnd(); ++sourceIt, ++targetIt )
    {
      targetIt.Value() = sourceIt.Value()[ 0 ];
    }

    typedef itk::ImageFileWriter< DebugImageType >  DebugWriterType;
    DebugWriterType::Pointer  debugWriter = DebugWriterType::New();
    debugWriter->SetInput( debugImage );
    debugWriter->SetFileName( "debugImage.mhd" );
    debugWriter->Update();

    std::cout << "Wrote debugImage.mhd" << std::endl;
  }
#endif



  // Compute the likelihoods in each big "averaging" voxel
  this->CalculateLikelihoods();


  // Using the pointers to each subvoxel's siblings, compute the statistics the rasterizor'
  // fragment processor will need
  int  numberOfClasses = mesh->GetPointData()->Begin().Value().m_Alphas.Size();
  itk::Array< float >  defaultGradientBasisWeight( numberOfClasses );
  defaultGradientBasisWeight.Fill( 0.0f );
  defaultGradientBasisWeight( 0 ) = 1.0f;
  m_GradientBasisWeights->FillBuffer( defaultGradientBasisWeight );
  itk::ImageRegionConstIterator< SiblingPriorsImageType >  siblingIt( m_SiblingPriorsImage,
      m_SiblingPriorsImage->GetBufferedRegion() );
  itk::ImageRegionConstIterator< LabelImageType >  intensityIt( m_HighResolutionImage,
      m_HighResolutionImage->GetBufferedRegion() );
  itk::ImageRegionIterator< GradientBasisWeightsType >  gradBWIt( m_GradientBasisWeights,
      m_GradientBasisWeights->GetBufferedRegion() );
  for ( ; !siblingIt.IsAtEnd(); ++siblingIt, ++intensityIt, ++gradBWIt )
  {
    gradBWIt.Value().Fill( 0 );

    // Contribution of pure tissue classes
    bool  dontBotherLookingFurther = false;
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
      // Calculate prior of that all siblings are this class
      float  priorOfSiblingCombination = 1;
      for ( int siblingNumber = 0; siblingNumber < m_NumberOfSiblings; siblingNumber++ )
      {
        priorOfSiblingCombination *= ( *( siblingIt.Value()[ siblingNumber ] ) )[ classNumber ];
      }

      // Calculate likelihood
      const float  mean = m_Means[ classNumber ];
      const float  variance = m_Variances[ classNumber ];
      const float  gauss = exp( -pow( intensityIt.Value() - mean, 2 ) / 2.0f / variance ) /
                           sqrt( 2 * 3.14 * variance );

      // Add contribution to corresponding class's gradient basis factor
      gradBWIt.Value()[ classNumber ] += gauss * priorOfSiblingCombination;

      // Check if we have certainty about this subvoxel's siblings all being this class. If so, then don't
      // waste time calculating zeros - especially since a large part of the high-res image will be background
      // and never be visited by the rasterizor since it falls outside of the mesh's boundary box.
      if ( priorOfSiblingCombination > 0.98 )
      {
        dontBotherLookingFurther = true;
        break;
      }
    }


    if ( dontBotherLookingFurther )
    {
      continue;
    }


    // Contribution of PV tissue classes
    int pvClassNumber = -1;
    const int  numberOfGaussiansPerPV = m_NumberOfSiblings;
    for ( int classNumber1 = 0; classNumber1 < numberOfClasses; classNumber1++ )
    {
      for ( int classNumber2 = classNumber1+1; classNumber2 < numberOfClasses; classNumber2++ )
      {
        pvClassNumber++;

        // Retrieve parameters of the pure tissue classes this PV class is mixing for
        const float  mean1 = m_Means[ classNumber1 ];
        const float  variance1 = m_Variances[ classNumber1 ];
        const float  mean2 = m_Means[ classNumber2 ];
        const float  variance2 = m_Variances[ classNumber2 ];


        itk::Array< float >  pureProbabilitiesOfClass1( m_NumberOfSiblings );
        itk::Array< float >  pureProbabilitiesOfClass2( m_NumberOfSiblings );
        for ( int siblingNumber = 0; siblingNumber < m_NumberOfSiblings; siblingNumber++ )
        {
          pureProbabilitiesOfClass1[ siblingNumber ] = ( *( siblingIt.Value()[ siblingNumber ] ) )[ classNumber1 ];
          pureProbabilitiesOfClass2[ siblingNumber ] = ( *( siblingIt.Value()[ siblingNumber ] ) )[ classNumber2 ];
        }
        itk::Array< float >  priorOfSiblingCombinations = this->GetMixingProbabilities( pureProbabilitiesOfClass1,
            pureProbabilitiesOfClass2 );
        for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussiansPerPV; gaussianNumber++ )
        {
          // Retrieve mean and variance for this sub-Gaussian
          const float  alpha = static_cast< float >( gaussianNumber + 1 ) / static_cast< float >( m_NumberOfSubvoxels );
          const float  mean = alpha * mean1 + ( 1 - alpha ) * mean2;
          const float  variance = alpha * variance1 + ( 1 - alpha ) * variance2;

          // Calculate likelihood
          const float  gauss = exp( -pow( intensityIt.Value() - mean, 2 ) / 2.0f / variance ) /
                               sqrt( 2 * 3.14 * variance );

          // Add contribution to class 1's gradient basis factor
          gradBWIt.Value()[ classNumber1 ] += gauss * priorOfSiblingCombinations[ gaussianNumber ];

          // Add contribution to class 2's gradient basis factor
          gradBWIt.Value()[ classNumber2 ] += gauss * priorOfSiblingCombinations[ gaussianNumber+1 ];

        } // End loop over all sub-gaussians of this PV class

      }
    } // End looping over PV tissue classes


  } // End looping over all sub-voxels


  // Invoke superclass' implementation
  Superclass::Rasterize( highResolutionMesh );

}



//
//
//
void
AtlasMeshToIntensityImagePartialVolumeGradientCalculator
::CalculateLikelihoods()
{

  //std::cout << "Calculating minLogLikelihoods" << std::endl;

  const int  numberOfClasses = ( *( m_RasterizedHighResolutionPriors->GetPixelContainer() ) )[0].Size();
  //std::cout << "      found that numberOfClasses = " << numberOfClasses << std::endl;

  //std::cout << "m_Means: " << m_Means << std::endl;
  //std::cout << "m_Variances: " << m_Variances << std::endl;


  // Loop over all voxels of the low-res intensity image, and calculate it's likelihood. Also
  // compute the entire minLogLikelihood at the same time.
  //
  // Calculating the likelihood in each low-res voxel is simple, because you know have an image
  // (m_SubvoxelPriorsImage) that points to the priors in each of the voxel's subvoxels.
  m_MinLogLikelihood = 0;
  itk::ImageRegionConstIterator< LabelImageType >  intensityIt( m_Image, m_Image->GetBufferedRegion() );
  itk::ImageRegionConstIterator< SubvoxelPriorsImageType >  subvoxelPriorIt( m_SubvoxelPriorsImage,
      m_SubvoxelPriorsImage->GetBufferedRegion() );
  itk::ImageRegionIterator< LikelihoodImageType >  likelihoodIt( m_LikelihoodImage,
      m_LikelihoodImage->GetBufferedRegion() );
  for ( ; !intensityIt.IsAtEnd(); ++intensityIt, ++subvoxelPriorIt, ++likelihoodIt )
  {
    float likelihood = 0;


    // Contribution of pure tissue classes
    bool  dontBotherLookingFurther = false;
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
      // Calculate prior that all subvoxels are this class
      float  priorOfSubvoxelCombination = 1;
      for ( int subvoxelNumber = 0; subvoxelNumber < m_NumberOfSubvoxels; subvoxelNumber++ )
      {
        priorOfSubvoxelCombination *= ( *( subvoxelPriorIt.Value()[ subvoxelNumber ] ) )[ classNumber ];
      }

      // Calculate likelihood
      const float  mean = m_Means[ classNumber ];
      const float  variance = m_Variances[ classNumber ];
      const float  gauss = exp( -pow( intensityIt.Value() - mean, 2 ) / 2.0f / variance ) /
                           sqrt( 2 * 3.14 * variance );

      // Add contribution to likelihood
      likelihood += gauss * priorOfSubvoxelCombination;

      // Check if we have certainty about this voxel's subvoxels all being this class. If so, then don't
      // waste time calculating zeros - especially since a large part of the high-res image will be background
      // and never be visited by the rasterizor since it falls outside of the mesh's boundary box.
      if ( priorOfSubvoxelCombination > 0.98 )
      {
        dontBotherLookingFurther = true;
        break;
      }
    }


    if ( !dontBotherLookingFurther )
    {

      // Contribution of PV tissue classes
      int pvClassNumber = -1;
      const int  numberOfGaussiansPerPV = m_NumberOfSiblings;
      for ( int classNumber1 = 0; classNumber1 < numberOfClasses; classNumber1++ )
      {
        for ( int classNumber2 = classNumber1+1; classNumber2 < numberOfClasses; classNumber2++ )
        {
          pvClassNumber++;

          // Retrieve parameters of the pure tissue classes this PV class is mixing for
          const float  mean1 = m_Means[ classNumber1 ];
          const float  variance1 = m_Variances[ classNumber1 ];
          const float  mean2 = m_Means[ classNumber2 ];
          const float  variance2 = m_Variances[ classNumber2 ];


          itk::Array< float >  pureProbabilitiesOfClass1( m_NumberOfSubvoxels );
          itk::Array< float >  pureProbabilitiesOfClass2( m_NumberOfSubvoxels );
          for ( int subvoxelNumber = 0; subvoxelNumber < m_NumberOfSubvoxels; subvoxelNumber++ )
          {
            pureProbabilitiesOfClass1[ subvoxelNumber ] = ( *( subvoxelPriorIt.Value()[ subvoxelNumber ] ) )[ classNumber1 ];
            pureProbabilitiesOfClass2[ subvoxelNumber ] = ( *( subvoxelPriorIt.Value()[ subvoxelNumber ] ) )[ classNumber2 ];
          }
          itk::Array< float >  priorOfSubvoxelCombinations = this->GetMixingProbabilities( pureProbabilitiesOfClass1,
              pureProbabilitiesOfClass2 );
          for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussiansPerPV; gaussianNumber++ )
          {
            // Retrieve mean and variance for this sub-Gaussian
            const float  alpha = static_cast< float >( gaussianNumber + 1 ) / static_cast< float >( m_NumberOfSubvoxels );
            const float  mean = alpha * mean1 + ( 1 - alpha ) * mean2;
            const float  variance = alpha * variance1 + ( 1 - alpha ) * variance2;

            // Calculate likelihood
            const float  gauss = exp( -pow( intensityIt.Value() - mean, 2 ) / 2.0f / variance ) /
                                 sqrt( 2 * 3.14 * variance );

            // Add contribution to likelihood
            likelihood += gauss * priorOfSubvoxelCombinations[ gaussianNumber + 1 ];

          } // End loop over all sub-gaussians of this PV class

        }
      } // End looping over PV tissue classes

    }

    likelihoodIt.Value() = ( likelihood + 1E-5 );
    m_MinLogLikelihood -= log( likelihood + 1E-15 );


#if 0
    std::cout << "    Just finished new voxel: " << std::endl;
    std::cout << "              intensity: " << static_cast< int >( intensityIt.Value() ) << std::endl;
    for ( int subvoxelNumber = 0; subvoxelNumber < m_NumberOfSubvoxels; subvoxelNumber++ )
    {
      std::cout << "             probabilities in subvoxel " << subvoxelNumber << ": "
                << *( subvoxelPriorIt.Value()[ subvoxelNumber ] ) << std::endl;
    }
    std::cout << "               likelihood: " << likelihood << "\n\n" << std::endl;
#endif


#if 0
    if ( ( likelihoodIt.GetIndex()[ 0 ] == 1 ) && ( likelihoodIt.GetIndex() [ 1 ] == 8 ) && ( likelihoodIt.GetIndex()[ 2 ] == 8 ) )
    {
      std::cout << "           This is likelihood voxel content with address " << &( likelihoodIt.Value() )
                << " and content " << likelihoodIt.Value() << std::endl;
      std::cout << "           low-res index is " << likelihoodIt.GetIndex() << std::endl;
      char  dummy;
      std::cin >> dummy;
    }

#endif

  }  // End loop over all low-resolution voxels


  //std::cout << "Done calculating likelihoods in the low-resolution image (m_MinLogLikelihood: "
  //          << m_MinLogLikelihood << ")" << std::endl;

}


//
//
//
itk::Array< float >
AtlasMeshToIntensityImagePartialVolumeGradientCalculator
::GetMixingProbabilities( const itk::Array< float >& pureProbabilitiesOfClass1,
                          const itk::Array< float >& pureProbabilitiesOfClass2 )
{

  /**
   Imagine a number of jars filled with white, black, and colored balls, where the jars' fraction
   of white balls is given by "pureProbabilities1" and of black balls by "pureProbabilitiesOfClass2".
   This function returns the probability of ending up with N white balls and (numberOfJars-N) in total
   if one would randomly pick one ball from each jar.
  */


  // This stuff is easily implemented using a recursive paradigm
  itk::Array< float >  mixingProbabilities( pureProbabilitiesOfClass1.Size() + 1 );
  mixingProbabilities.Fill( 0 );
  if ( pureProbabilitiesOfClass1.Size() == 1 )
  {
    // Trivial case
    mixingProbabilities[ 0 ] = pureProbabilitiesOfClass2[ 0 ];
    mixingProbabilities[ 1 ] = pureProbabilitiesOfClass1[ 0 ];
  }
  else
  {
    // Non-trivial case: we can end up with N white balls either by drawing black
    // from the first jar and N times white from the remaining jars OR drawing white
    // from the first jar and (N-1) times white from the remaining jars
    itk::Array< float >  pureProbabilitiesOfClass1OfOfRemainingJars( pureProbabilitiesOfClass1.Size()-1 );
    itk::Array< float >  pureProbabilitiesOfClass2OfOfRemainingJars( pureProbabilitiesOfClass1.Size()-1 );
    for ( unsigned int i = 1; i < pureProbabilitiesOfClass1.Size(); i++ )
    {
      pureProbabilitiesOfClass1OfOfRemainingJars[ i-1 ] = pureProbabilitiesOfClass1[ i ];
      pureProbabilitiesOfClass2OfOfRemainingJars[ i-1 ] = pureProbabilitiesOfClass2[ i ];
    }

    itk::Array< float >  otherMixingProbabilities = Self::GetMixingProbabilities( pureProbabilitiesOfClass1OfOfRemainingJars,
        pureProbabilitiesOfClass2OfOfRemainingJars );
    for ( unsigned int i = 0; i < ( mixingProbabilities.Size() - 1 ); i++ )
    {
      mixingProbabilities[ i ] = pureProbabilitiesOfClass2[ 0 ] * otherMixingProbabilities[ i ];
    }

    for ( unsigned int i = 1; i < mixingProbabilities.Size(); i++ )
    {
      mixingProbabilities[ i ] += pureProbabilitiesOfClass1[ 0 ] * otherMixingProbabilities[ i-1 ];
    }

  }


  return mixingProbabilities;

}



} // end namespace kvl
