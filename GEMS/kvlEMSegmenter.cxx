/**
 * @file  kvlEMSegmenter.cxx
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
#include "kvlEMSegmenter.h"

#include "kvlAtlasMeshAlphaDrawer.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "kvlAtlasMeshMultiAlphaDrawer.h"
#include "itkCastImageFilter.h"


namespace kvl
{

//
//
//
EMSegmenter
::EMSegmenter()
{
  m_Image = 0;
  m_BiasCorrectedImage = 0;
  m_BiasField = 0;
  m_AtlasMesh = 0;

  m_BiasFieldOrder = 4;

  m_ReestimatedRelativeWeightOfSharedClasses = true;

  m_IterationNumber = 0;
  m_MaximumNumberOfIterations = 100;
  m_StopCriterion = 1e-7;

  int  partialVolumeUpsamplingFactors[ 3 ] = { 1, 1, 1 };
  this->SetPartialVolumeUpsamplingFactors( partialVolumeUpsamplingFactors );

  m_SuperResolutionPriors = 0;
  m_SubvoxelPriorsImage = 0;


  m_BiasFieldResidueSmoother = 0;

}




//
//
//
EMSegmenter
::~EMSegmenter()
{

}


//
//
//
void
EMSegmenter
::SetImage( const ImageType*  image )
{

  m_Image = image;

  //
  m_Priors.clear();
  m_Posteriors.clear();
  //m_Means.clear();
  //m_Variances.clear();
  //m_PriorWeights.clear();

  // Initialize bias field corrected image to image itself
  typedef itk::CastImageFilter< ImageType, BiasCorrectedImageType >  CasterType;
  CasterType::Pointer  caster = CasterType::New();
  caster->SetInput( m_Image );
  caster->Update();
  m_BiasCorrectedImage = caster->GetOutput();

}



//
//
//
double
EMSegmenter
::Segment()
{
  // Sanity check on input
  if ( !m_Image || !m_AtlasMesh )
  {
    itkExceptionMacro( << "Insufficient inputs set!" );
  }

  std::cout << "Starting EM: " << std::endl;
  //std::cout << "         m_MaximumNumberOfIterations: " << m_MaximumNumberOfIterations << std::endl;


  // First initialize if not already done so
  this->Initialize();

  // Now loop over the EM-steps
  this->InvokeEvent( itk::StartEvent() );

  // If no previous estimates for the parameters are available,
  // obtain some initialization by usings the priors as "posteriors"
  if ( m_Means.size() == 0 )
  {
    this->UpdateMixtureModelParameters();

#if 1
    {
      std::cout << "m_Means: [ ";
      for ( unsigned int i = 0; i < m_Means.size(); i++ )
      {
        std::cout << " " << m_Means[ i ];
      }
      std::cout << " ]" << std::endl;
    }
#endif

    // Split the Gaussians that have the same priors
    this->SplitSharedClasses();
  }


  double  currentCost = 0;
  if ( m_NumberOfSubvoxels == 1 )
  {
    // Start with bias field order 0 - this will be gradually increased later on
    int biasFieldOrderUsed  = 0;
    currentCost = itk::NumericTraits< double >::max();
    for ( m_IterationNumber = 0; m_IterationNumber < m_MaximumNumberOfIterations; m_IterationNumber++ )
    {
      // E-step
      const double  previousCost = currentCost;
      currentCost = this->UpdatePosteriors();

      std::cout << "EM Iteration: " << m_IterationNumber << "  ->  " << currentCost << std::endl;


      if ( ( ( previousCost - currentCost ) / fabsf( currentCost ) ) < m_StopCriterion )
      {
        std::cout << "EM converged" << std::endl;
        //std::cout << "( ( " << previousCost << " - " << currentCost << " ) / " << " fabsf( " << currentCost << " ) ) = "
        //          << ( ( previousCost - currentCost ) / fabsf( currentCost ) ) << " < "
        //          << m_StopCriterion << std::endl;
        if ( biasFieldOrderUsed == m_BiasFieldOrder )
        {
          std::cout << "Maximal bias field order already reached; EM finished" << std::endl;
          break;
        }
        else
        {
          biasFieldOrderUsed++;
          std::cout << "Increasing bias field order to " << biasFieldOrderUsed << std::endl;
        }
      }

      // M-step
      this->UpdateModelParameters( biasFieldOrderUsed );
#if 1
      {
        std::cout << "m_Means: [ ";
        for ( unsigned int i = 0; i < m_Means.size(); i++ )
        {
          std::cout << " " << m_Means[ i ];
        }
        std::cout << " ]" << std::endl;
      }
#endif


      this->InvokeEvent( itk::IterationEvent() );
    } // End loop over iterations

  }
  else
  {
    currentCost = itk::NumericTraits< double >::max();
    for ( m_IterationNumber = 0; m_IterationNumber < m_MaximumNumberOfIterations; m_IterationNumber++ )
    {
      const double  previousCost = currentCost;
      currentCost = this->DoOneEMIterationWithPartialVoluming();

      if ( ( ( previousCost - currentCost ) / fabsf( currentCost ) ) < m_StopCriterion )
      {
        std::cout << "EM converged" << std::endl;
        break;
      }

      std::cout << "EM Iteration: " << m_IterationNumber << "  ->  " << currentCost << std::endl;
      this->InvokeEvent( itk::IterationEvent() );
    }
  }


  std::cout << "Stopped EM" << std::endl;
  std::cout << "         m_IterationNumber: " << m_IterationNumber << std::endl;
  std::cout << "         m_MaximumNumberOfIterations: " << m_MaximumNumberOfIterations << std::endl;
  this->InvokeEvent( itk::EndEvent() );

  return  currentCost;
}




//
//
//
AtlasMesh::Pointer
EMSegmenter
::GetDistributedMesh( const AtlasMesh* atlasMesh, const std::vector< unsigned int >& lookupTable ) const
{

  // Pre-calculate in how much classes each atlas prior will be divided
  std::vector< unsigned int >  numberOfFinalClassesPerOriginalClass;
  for ( unsigned int finalClassNumber = 0; finalClassNumber < lookupTable.size(); finalClassNumber++ )
  {
    // Make sure we have enough memory allocated
    if ( numberOfFinalClassesPerOriginalClass.size() < ( lookupTable[ finalClassNumber ] + 1 ) )
    {
      numberOfFinalClassesPerOriginalClass.resize( lookupTable[ finalClassNumber ] + 1, 0 );
    }

    numberOfFinalClassesPerOriginalClass[ lookupTable[ finalClassNumber ] ]++;
  }
  for ( unsigned int originalClassNumber = 0; originalClassNumber < numberOfFinalClassesPerOriginalClass.size(); originalClassNumber++ )
  {
    std::cout << "Original class " << originalClassNumber << " is distributed over "
              << numberOfFinalClassesPerOriginalClass[ originalClassNumber ] << " classes" << std::endl;

  }

  // Copy the original point parameters
  AtlasMesh::PointDataContainer::Pointer  distributedPointParameters = AtlasMesh::PointDataContainer::New();
  for ( AtlasMesh::PointDataContainer::ConstIterator  it = atlasMesh->GetPointData()->Begin();
        it != atlasMesh->GetPointData()->End(); ++it )
  {
    distributedPointParameters->InsertElement( it.Index(), it.Value() );
  }

  // Adjust the alphas in each vertex
  for ( AtlasMesh::PointDataContainer::Iterator  it = distributedPointParameters->Begin();
        it != distributedPointParameters->End(); ++it )
  {
    // Retrieve the original alphas
    AtlasAlphasType  originalAlphas = it.Value().m_Alphas;

    // Loop over all classes in the final configuration
    AtlasAlphasType  finalAlphas( lookupTable.size() );
    for ( unsigned int finalClassNumber = 0; finalClassNumber < lookupTable.size(); finalClassNumber++ )
    {
      const unsigned int  originalClassNumber = lookupTable[ finalClassNumber ];
      const unsigned int  distributionFactor = numberOfFinalClassesPerOriginalClass[ originalClassNumber ];
      finalAlphas[ finalClassNumber ] = originalAlphas[ originalClassNumber ] /
                                        static_cast< float >( distributionFactor );
    }

    // Set the alpha in this vertex to the final alphas
    it.Value().m_Alphas = finalAlphas;
  }

  // Construct the final atlas mesh from its components
  AtlasMesh::Pointer  distributedMesh = AtlasMesh::New();
  distributedMesh->SetPoints( const_cast< AtlasMesh::PointsContainer* >( atlasMesh->GetPoints() ) );
  distributedMesh->SetCells(
    const_cast< AtlasMesh::CellsContainer* >( atlasMesh->GetCells() ) );
  distributedMesh->SetPointData( distributedPointParameters );
  distributedMesh->SetCellData(
    const_cast< AtlasMesh::CellDataContainer* >( atlasMesh->GetCellData() ) );

  return distributedMesh.GetPointer();
}



//
//
//
void
EMSegmenter
::SplitSharedClasses()
{

  std::cout << "Splitting shared classes" << std::endl;

  // Build inverse lookup table
  std::vector< std::vector< unsigned int > >  inverseLookupTable;
  for ( unsigned int finalClassNumber = 0; finalClassNumber < m_LookupTable.size(); finalClassNumber++ )
  {
    // Make sure we have enough memory allocated
    if ( inverseLookupTable.size() < ( m_LookupTable[ finalClassNumber ] + 1 ) )
    {
      inverseLookupTable.resize( m_LookupTable[ finalClassNumber ] + 1 );
    }

    //
    inverseLookupTable[ m_LookupTable[ finalClassNumber ] ].push_back( finalClassNumber );
  }

  for ( unsigned int originalClassNumber = 0;
        originalClassNumber < inverseLookupTable.size();
        originalClassNumber++ )
  {
    const unsigned int  numberOfSharedClasses = inverseLookupTable[ originalClassNumber ].size();
    if ( numberOfSharedClasses < 2 )
    {
      continue;
    }

    const float  sharedMean = m_Means[ inverseLookupTable[ originalClassNumber ][ 0 ] ];
    const float  sharedVariance = m_Variances[ inverseLookupTable[ originalClassNumber ][ 0 ] ];
    const float  minimumOfRange = sharedMean - 2 * sqrt( sharedVariance );
    const float  maximumOfRange = sharedMean + 2 * sqrt( sharedVariance );
    const float  slope = ( maximumOfRange - minimumOfRange ) /
                         static_cast< float >( numberOfSharedClasses - 1 );
    const float  offset = minimumOfRange;

    std::cout << "Original prior " << originalClassNumber << " is distributed over classes: ";
    for ( unsigned int i = 0; i < numberOfSharedClasses; i++ )
    {
      const unsigned int finalClassNumber = inverseLookupTable[ originalClassNumber ][ i ];
      std::cout << finalClassNumber << " ";
      m_Means[ finalClassNumber ] = slope * i + offset;
    }
    std::cout << std::endl;
  }


}



//
//
//
void
EMSegmenter
::SetAtlasMesh( const AtlasMesh* atlasMesh,
                const std::vector< unsigned int >&  lookupTable,
                const std::vector< unsigned int >& independentParametersLookupTable )
{

  if ( lookupTable.size() == 0 )
  {
    m_LookupTable.clear();
    const unsigned int  numberOfClasses = atlasMesh->GetPointData()->Begin().Value().m_Alphas.Size();
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
      m_LookupTable.push_back( classNumber );
    }
  }
  else
  {
    m_LookupTable = lookupTable;
  }

  m_AtlasMesh = this->GetDistributedMesh( atlasMesh, m_LookupTable );

  m_NumberOfClasses = m_AtlasMesh->GetPointData()->Begin().Value().m_Alphas.Size();

  if ( independentParametersLookupTable.size() == 0 )
  {
    m_IndependentParametersLookupTable.clear();
    for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
    {
      m_IndependentParametersLookupTable.push_back( classNumber );
    }
  }
  else
  {
    m_IndependentParametersLookupTable = independentParametersLookupTable;
  }


  m_Priors.clear();
  m_Posteriors.clear();
  m_Means.clear();
  m_Variances.clear();
  m_PriorWeights.clear();
}



//
//
//
void
EMSegmenter
::Initialize() const
{

  if ( m_NumberOfSubvoxels == 1 )
  {
    // Create an empty image that serves as a template for the alpha rasterizor
    typedef AtlasMeshAlphaDrawer::LabelImageType  LabelImageType;
    LabelImageType::Pointer  templateImage = LabelImageType::New();
    templateImage->SetRegions( m_Image->GetLargestPossibleRegion().GetSize() );
    templateImage->Allocate();
    templateImage->FillBuffer( 0 );

    // Now rasterize the prior for each class one class at a time
    m_Priors.clear();
    for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
    {
      AtlasMeshAlphaDrawer::Pointer  alphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
      alphaDrawer->SetLabelImage( templateImage );
      alphaDrawer->SetLabelNumber( classNumber );
      alphaDrawer->Rasterize( m_AtlasMesh );

      AtlasMeshAlphaDrawer::AlphaImageType::Pointer  prior =
        const_cast< AtlasMeshAlphaDrawer::AlphaImageType* >( alphaDrawer->GetAlphaImage() );
#if 1
      // Voxels with intensity 0 are not considered: set prior to zero there
      itk::ImageRegionConstIterator< ImageType >  imageIt( m_Image, m_Image->GetLargestPossibleRegion() );
      itk::ImageRegionIterator< AtlasMeshAlphaDrawer::AlphaImageType >  priorIt( prior, prior->GetLargestPossibleRegion() );
      for ( ; !imageIt.IsAtEnd(); ++imageIt, ++priorIt )
      {
        // if ( classNumber == 0 )
        //   {
        //   if ( priorIt.Value() > ( 0.9 / 3.0 ) )
        //     {
        //     const_cast< ImageType::PixelType& >( imageIt.Value() ) = 0;
        //     }
        //   }

        if ( imageIt.Value() == 0 )
        {
          priorIt.Value() = 0;
        }
#if 1
        if ( priorIt.Value() < 0 )
        {
          priorIt.Value() = 0;
        }
#endif

      }
#endif

      m_Priors.push_back( prior );
    }


#if 0
    // Voxels with prior for background large are not considered: set priors to zero there
    for ( int classNumber = m_NumberOfClasses-1; classNumber > -1; classNumber-- )
    {
      itk::ImageRegionIterator< AtlasMeshAlphaDrawer::AlphaImageType >  priorIt( m_Priors[ classNumber ],
          m_Priors[ classNumber ]->GetLargestPossibleRegion() );
      itk::ImageRegionConstIterator< AtlasMeshAlphaDrawer::AlphaImageType >  firstPriorIt( m_Priors[ 0 ],
          m_Priors[ 0 ]->GetLargestPossibleRegion() );
      for ( ; !priorIt.IsAtEnd(); ++priorIt, ++firstPriorIt )
      {
        if ( firstPriorIt.Value() > 0.90 )
        {
          priorIt.Value() = 0;
        }
      }
    }
#endif


  }
  else
  {
    std::cout << "Initializing for partialVolumeUpsamplingFactors " << m_PartialVolumeUpsamplingFactors[ 0 ]
              << " " << m_PartialVolumeUpsamplingFactors[ 1 ] << " " << m_PartialVolumeUpsamplingFactors[ 2 ] << std::endl;

    // Allocate a high-resolution image that will hold the rasterized priors
    ImageType::SizeType  highResolutionSize;
    for ( int i = 0; i < 3; i++ )
    {
      highResolutionSize[ i ] = m_Image->GetBufferedRegion().GetSize()[ i ] * m_PartialVolumeUpsamplingFactors[ i ];
    }
    std::cout << "  highResolutionSize: " << highResolutionSize << std::endl;

    m_SuperResolutionPriors = SuperResolutionProbabilityImageType::New();
    m_SuperResolutionPriors->SetRegions( highResolutionSize );
    m_SuperResolutionPriors->Allocate();


    // Construct a mesh that lives in the high-res space
#if 0
    float  translation[ 3 ];
    float  scaling[ 3 ];
    for ( int i = 0; i < 3; i++ )
    {
      translation[ i ] = ( m_PartialVolumeUpsamplingFactors[ i ] - 1 ) / 2.0f;
      scaling[ i ] = m_PartialVolumeUpsamplingFactors[ i ];
    }

    AtlasMesh::PointsContainer::Pointer  highResolutionPoints = AtlasMesh::PointsContainer::New();
    for ( AtlasMesh::PointsContainer::ConstIterator it = m_AtlasMesh->GetPoints()->Begin();
          it != m_AtlasMesh->GetPoints()->End(); ++it )
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
    highResolutionMesh->SetCells( m_AtlasMesh->GetCells() );
    highResolutionMesh->SetPointData( const_cast< AtlasMesh::PointDataContainer* >( m_AtlasMesh->GetPointData() ) );
    highResolutionMesh->SetCellData( m_AtlasMesh->GetCellData() );
#else
    AtlasMesh::Pointer  highResolutionMesh = this->GetSuperResolutionMesh( m_AtlasMesh );
#endif


    // Create a high-resolution template image to be passed on to the rasterizor
    typedef AtlasMeshMultiAlphaDrawer::LabelImageType  LabelImageType;
    LabelImageType::Pointer  templateImage = LabelImageType::New();
    templateImage->SetRegions( m_SuperResolutionPriors->GetLargestPossibleRegion().GetSize() );
    templateImage->Allocate();
    templateImage->FillBuffer( 0 );


    // Rasterize the high-resoluation mesh into the prior probabilities in each high-resolution voxel
    AtlasMeshMultiAlphaDrawer::Pointer  drawer = AtlasMeshMultiAlphaDrawer::New();
    drawer->SetLabelImage( templateImage );
    drawer->SetAlphasImage( m_SuperResolutionPriors );
    drawer->Rasterize( highResolutionMesh );


    // Set prior probabilities in voxels with intensity 0 to 0 for all classes
    std::cout << "Setting priors for zero-intensity voxels to zero everywhere..." << std::endl;
    ImageType::ConstPointer  superResolutionImage = this->GetSuperResolutionImage().GetPointer();
    itk::ImageRegionConstIterator< ImageType >  intensityIt( superResolutionImage,
        superResolutionImage->GetBufferedRegion() );
    itk::ImageRegionIterator< SuperResolutionProbabilityImageType >  priorIt( m_SuperResolutionPriors,
        m_SuperResolutionPriors->GetBufferedRegion() );
    for ( ; !intensityIt.IsAtEnd(); ++intensityIt, ++priorIt )
    {
      if ( intensityIt.Value() == 0 )
      {
        priorIt.Value().Fill( 0.0f );
      }
    }


    // Allocate a low-resolution image that will hold pointers to the prior probabilities in
    // each of its voxel's subvoxels, and fill it in
    m_SubvoxelPriorsImage = SubvoxelPriorsImageType::New();
    m_SubvoxelPriorsImage->SetRegions( m_Image->GetBufferedRegion().GetSize() );
    m_SubvoxelPriorsImage->Allocate();


    for ( itk::ImageRegionIteratorWithIndex< SubvoxelPriorsImageType >  it( m_SubvoxelPriorsImage,
          m_SubvoxelPriorsImage->GetBufferedRegion() );
          !it.IsAtEnd(); ++it )
    {
      //std::cout << "Filling in pointers to sibling priors in big voxel with index " << it.GetIndex() << std::endl;

      // Collect *all* pointers in this voxel's subvoxels
      std::vector< AtlasAlphasType* >  priorsOfSubvoxels( m_NumberOfSubvoxels );
      ImageType::IndexType  index;
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

            priorsOfSubvoxels[ counter ] = &( m_SuperResolutionPriors->GetPixel( index ) );

          }
        }
      }

      // Assign it in the low-res image m_SubvoxelPriorsImage
      it.Value() = priorsOfSubvoxels;

    } // End loop over all big, "averaging" voxels


    // Now allocate and fill in the prior expectations; these are not actually useful in the EM algorithm but they will help
    // the user of this class in visualizing the algorithm's progression
    m_Priors.clear();
    for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
    {
      ClassificationImageType::Pointer  prior = ClassificationImageType::New();
      prior->SetRegions( m_Image->GetLargestPossibleRegion().GetSize() );
      prior->Allocate();
      prior->FillBuffer( 0 );

#if 0
      // Make sure spacing, origin, and directions are copied
      prior->SetSpacing( m_Image->GetSpacing() );
      prior->SetOrigin( m_Image->GetOrigin() );
      prior->SetDirection( m_Image->GetDirection() );
#endif

      m_Priors.push_back( prior.GetPointer() );
    }

    itk::ImageRegionConstIterator< SubvoxelPriorsImageType >  subvoxelPriorIt( m_SubvoxelPriorsImage,
        m_SubvoxelPriorsImage->GetBufferedRegion() );
    std::vector< itk::ImageRegionIterator< ClassificationImageType > >   priorIts;
    for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
    {
      itk::ImageRegionIterator< ClassificationImageType >  priorIt( m_Priors[ classNumber ],
          m_Priors[ classNumber ]->GetBufferedRegion() );
      priorIts.push_back( priorIt );
    }

    for ( ; !subvoxelPriorIt.IsAtEnd(); )
    {
      // Contribution of pure tissue classes
      bool  dontBotherLookingFurther = false;
      for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
      {
        // Calculate prior that all subvoxels are this class
        float  priorOfSubvoxelCombination = 1;
        for ( int subvoxelNumber = 0; subvoxelNumber < m_NumberOfSubvoxels; subvoxelNumber++ )
        {
          priorOfSubvoxelCombination *= ( *( subvoxelPriorIt.Value()[ subvoxelNumber ] ) )[ classNumber ];
        }

        // Add contribution to prior expectation
        ( priorIts[ classNumber ] ).Value() += priorOfSubvoxelCombination;

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
        const int  numberOfGaussiansPerPV = m_NumberOfSubvoxels-1;
        for ( unsigned int classNumber1 = 0; classNumber1 < m_NumberOfClasses; classNumber1++ )
        {
          for ( unsigned int classNumber2 = classNumber1+1; classNumber2 < m_NumberOfClasses; classNumber2++ )
          {
            pvClassNumber++;


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
              // Add contribution to prior expectations
              const float  alpha = static_cast< float >( gaussianNumber + 1 ) / static_cast< float >( m_NumberOfSubvoxels );
              ( priorIts[ classNumber1 ] ).Value() += alpha * priorOfSubvoxelCombinations[ gaussianNumber + 1 ];
              ( priorIts[ classNumber2 ] ).Value() += ( 1 - alpha ) * priorOfSubvoxelCombinations[ gaussianNumber + 1 ];

            } // End loop over all sub-gaussians of this PV class

          }
        } // End looping over PV tissue classes

      } // End test if we need to bother going into PV classes


      // Normalize the prior: for more than 2 classes, we're not allowing mixing of 3 or more
      float  sum = 1e-15f;
      for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
      {
        sum += ( priorIts[ classNumber ] ).Value();
      }
      for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
      {
        ( priorIts[ classNumber ] ).Value() /= sum;
      }


      // Advance to the next voxel
      ++subvoxelPriorIt;
      for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
      {
        ++( priorIts[ classNumber ] );
      }

    }  // End loop over all low-resolution voxels



  } // End test if we're using partial volume modeling or not






  // Initialize the posteriors with the priors
  m_Posteriors.clear();
  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
  {
    ClassificationImageType::ConstPointer  prior = m_Priors[ classNumber ].GetPointer();
    ClassificationImageType::Pointer  posterior = ClassificationImageType::New();
    posterior->SetRegions( prior->GetLargestPossibleRegion().GetSize() );
    posterior->Allocate();

    itk::ImageRegionConstIterator< ClassificationImageType >  priorIt( prior, prior->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< ClassificationImageType >  posteriorIt( posterior, posterior->GetLargestPossibleRegion() );
    for ( ; !priorIt.IsAtEnd(); ++priorIt, ++posteriorIt )
    {
      posteriorIt.Value() = priorIt.Value();
    }

    m_Posteriors.push_back( posterior );
  }

}



//
//
//
double
EMSegmenter
::UpdatePosteriors()
{

  // Loop over all classes, and update the unnormalized posteriors
  ClassificationImageType::Pointer  likelihoodImage = ClassificationImageType::New();
  likelihoodImage->SetRegions( m_Image->GetLargestPossibleRegion().GetSize() );
  likelihoodImage->Allocate();
  likelihoodImage->FillBuffer( 0 );

  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
  {
    const float  mean = m_Means[ classNumber ];
    const float  variance = m_Variances[ classNumber ];
    const float  weight = m_PriorWeights[ classNumber ];

    BiasCorrectedImageType::ConstPointer  imageToUse = m_BiasCorrectedImage.GetPointer();
#if 0
    // If this is the very first class, it isn't affected by the bias field...
    if ( classNumber == 0 )
    {
      typedef itk::CastImageFilter< ImageType, BiasCorrectedImageType >  CasterType;
      CasterType::Pointer  caster = CasterType::New();
      caster->SetInput( m_Image );
      caster->Update();
      imageToUse = caster->GetOutput();
    }
#endif

    itk::ImageRegionConstIterator< BiasCorrectedImageType >  imageIt( imageToUse,
        imageToUse->GetLargestPossibleRegion() );
    itk::ImageRegionConstIterator< ClassificationImageType >
    priorIt( m_Priors[ classNumber ], m_Priors[ classNumber ]->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< ClassificationImageType >
    posteriorIt( m_Posteriors[ classNumber ], m_Posteriors[ classNumber ]->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< ClassificationImageType >
    likelihoodIt( likelihoodImage, likelihoodImage->GetLargestPossibleRegion() );
    for ( ; !imageIt.IsAtEnd(); ++imageIt, ++priorIt, ++posteriorIt, ++likelihoodIt )
    {
      const float  unnormalizedPosterior = exp( -pow( imageIt.Value() - mean, 2 ) / variance / 2 ) /
                                           sqrt( 2 * 3.14 * variance ) *
                                           weight * priorIt.Value();
#if 0
      {
        if ( unnormalizedPosterior < 0 )
        {
          std::cout << "unnormalizedPosterior: " << unnormalizedPosterior << std::endl;
          std::cout << "    imageIt.Value(): " << imageIt.Value() << std::endl;
          std::cout << "    mean:" << mean << std::endl;
          std::cout << "    variance:" << variance << std::endl;
          std::cout << "    weight:" << weight << std::endl;
          std::cout << "    priorIt.Value():" << priorIt.Value() << std::endl;

          exit( -1 );
        }
      }
#endif
      posteriorIt.Value() = unnormalizedPosterior;
      likelihoodIt.Value() += unnormalizedPosterior;
    }

  }


  // Loop again over all classes, this time normalizing
  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
  {
    itk::ImageRegionConstIterator< ClassificationImageType >
    likelihoodIt( likelihoodImage, likelihoodImage->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< ClassificationImageType >
    posteriorIt( m_Posteriors[ classNumber ], m_Posteriors[ classNumber ]->GetLargestPossibleRegion() );
    for ( ; !likelihoodIt.IsAtEnd(); ++likelihoodIt, ++posteriorIt )
    {
      posteriorIt.Value() /= ( likelihoodIt.Value() + 1e-15 );
    }

  }


  // Calculate the likelihood
  double  minLogLikelihood = 0.0f;
  itk::ImageRegionConstIterator< ClassificationImageType >
  likelihoodIt( likelihoodImage, likelihoodImage->GetLargestPossibleRegion() );

  for ( ; !likelihoodIt.IsAtEnd(); ++likelihoodIt )
  {
    if ( !likelihoodIt.Value() )
    {
      // Skip voxels that are not being considered
      continue;
    }

    minLogLikelihood += -log( likelihoodIt.Value() + 1e-15 );
#if 0
    {
      if ( std::isnan( minLogLikelihood ) )
      {
        std::cout << "likelihoodIt.Value(): " << likelihoodIt.Value() << std::endl;
        exit( -1 );
      }
    }
#endif
  }


  return  minLogLikelihood;
}



//
//
//
void
EMSegmenter
::UpdateModelParameters( int biasFieldOrderUsed )
{
  this->UpdateMixtureModelParameters();
  this->UpdateBiasFieldModelParameters( biasFieldOrderUsed );

}

//
//
//
void
EMSegmenter
::UpdateBiasFieldModelParameters( int biasFieldOrderUsed )
{

  // if ( m_BiasFieldOrder == 0 )
  //   {
  //   // Nothing to do
  //   return;
  //   }

  std::cout << "Updating bias field model parameters" << std::endl;

  // Set up residue smoother
  if ( !m_BiasFieldResidueSmoother )
  {
    std::cout << "Setting up m_BiasFieldResidueSmoother" << std::endl;

    // Make a mask image
    typedef ImageSmoother::MaskImageType  MaskImageType;
    MaskImageType::Pointer  maskImage = MaskImageType::New();
    maskImage->SetRegions( m_Image->GetBufferedRegion() );
    maskImage->Allocate();

    itk::ImageRegionConstIterator< ImageType >  imageIt( m_Image, m_Image->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< MaskImageType >  maskIt( maskImage, maskImage->GetLargestPossibleRegion() );
    for ( ; !imageIt.IsAtEnd(); ++imageIt, ++maskIt )
    {
      if ( imageIt.Value() == 0 )
      {
        maskIt.Value() = false;
      }
      else
      {
        maskIt.Value() = true;
      }
    }

    // Create a smoother and pass the mask image to it
    m_BiasFieldResidueSmoother = ImageSmoother::New();
    m_BiasFieldResidueSmoother->SetMaskImage( maskImage );
    m_BiasFieldResidueSmoother->SetPolynomialOrder( m_BiasFieldOrder );
  }


  // Construct a residue and a weight image
  typedef ImageSmoother::ImageType  ResidueImageType;
  ResidueImageType::Pointer  residueImage = ResidueImageType::New();
  residueImage->SetRegions( m_Image->GetBufferedRegion() );
  residueImage->Allocate();
  residueImage->FillBuffer( 0 );

  typedef ImageSmoother::ImageType  WeightImageType;
  WeightImageType::Pointer  weightImage = WeightImageType::New();
  weightImage->SetRegions( m_Image->GetBufferedRegion() );
  weightImage->Allocate();
  weightImage->FillBuffer( 0 );

  // Loop over all classes, each time adding contribution to residue and weight.
#if 0
  // Remember that the very first class if modeled as not being affected by the bias field
  for ( unsigned int classNumber = 1; classNumber < m_NumberOfClasses; classNumber++ )
#else
  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
#endif
  {
    const float  mean = m_Means[ classNumber ];
    const float  variance = m_Variances[ classNumber ];

    // Loop over all voxels
    itk::ImageRegionConstIterator< ClassificationImageType >
    posteriorIt( m_Posteriors[ classNumber ],
                 m_Posteriors[ classNumber ]->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< ResidueImageType >  residueIt( residueImage,
        residueImage->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< WeightImageType >  weightIt( weightImage,
        weightImage->GetLargestPossibleRegion() );

    for ( ; !posteriorIt.IsAtEnd(); ++posteriorIt, ++residueIt, ++weightIt )
    {
      // Using notation of the paper
      const float  w_ij = posteriorIt.Value() / variance;

      residueIt.Value() -= w_ij * mean;
      weightIt.Value() += w_ij;
    }

  } // End loop over all classes


  // Loop over all voxels, complementing remaining residue terms
  itk::ImageRegionConstIterator< WeightImageType >  weightIt( weightImage,
      weightImage->GetLargestPossibleRegion() );
  itk::ImageRegionConstIterator< ImageType >  imageIt( m_Image,
      m_Image->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ResidueImageType >  residueIt( residueImage,
      residueImage->GetLargestPossibleRegion() );
  for ( ; !imageIt.IsAtEnd(); ++imageIt, ++weightIt, ++residueIt )
  {
    residueIt.Value() /= ( weightIt.Value() + 1e-15 );
    residueIt.Value() += imageIt.Value();
  }


  // Calculate bias field
  m_BiasFieldResidueSmoother->SetImage( residueImage );
  m_BiasFieldResidueSmoother->SetWeightImage( weightImage );
  m_BiasFieldResidueSmoother->SetPolynomialOrderUsed( biasFieldOrderUsed );
  m_BiasField = m_BiasFieldResidueSmoother->GetSmoothedImage().GetPointer();


  // Calculate corrected image
  imageIt.GoToBegin();
  itk::ImageRegionConstIterator< BiasFieldImageType >  biasIt( m_BiasField,
      m_BiasField->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< BiasCorrectedImageType >  correctedIt( m_BiasCorrectedImage,
      m_BiasCorrectedImage->GetLargestPossibleRegion() );
  for ( ; !imageIt.IsAtEnd(); ++imageIt, ++correctedIt, ++biasIt )
  {
    correctedIt.Value() = imageIt.Value() - biasIt.Value();
  }


}



//
//
//
void
EMSegmenter
::UpdateMixtureModelParameters()
{
  std::cout << "Updating mixture model parameters" << std::endl;

  // Loop over all classes, adding the contributions for the parameter calculations
  std::vector< float >  independentLinearSum;
  std::vector< float >  independentQuadraticSum;
  std::vector< float >  independentNumberOfPixels;
  std::vector< float >  unnormalizedWeights( m_NumberOfClasses );
  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
  {
    BiasCorrectedImageType::ConstPointer  imageToUse = m_BiasCorrectedImage.GetPointer();
#if 0
    // If this is the very first class, it isn't affected by the bias field...
    if ( classNumber == 0 )
    {
      typedef itk::CastImageFilter< ImageType, BiasCorrectedImageType >  CasterType;
      CasterType::Pointer  caster = CasterType::New();
      caster->SetInput( m_Image );
      caster->Update();
      imageToUse = caster->GetOutput();
    }
#endif

    itk::ImageRegionConstIterator< BiasCorrectedImageType >  imageIt( imageToUse,
        imageToUse->GetLargestPossibleRegion() );
    itk::ImageRegionConstIterator< ClassificationImageType >
    posteriorIt( m_Posteriors[ classNumber ], m_Posteriors[ classNumber ]->GetLargestPossibleRegion() );
    float  linearSum = 0.0f;
    float  quadraticSum = 0.0f;
    float  numberOfPixels = 0.0f;
    for ( ; !imageIt.IsAtEnd(); ++imageIt, ++posteriorIt )
    {
      linearSum += posteriorIt.Value() * imageIt.Value();
      quadraticSum += posteriorIt.Value() * pow( imageIt.Value() , 2 );
      numberOfPixels += posteriorIt.Value();
    }
    unnormalizedWeights[ classNumber ] = numberOfPixels;

    // Add contributions to the correct independent class; first make sure we have enough memory allocated
    const unsigned int  independentClassNumber = m_IndependentParametersLookupTable[ classNumber ];
    // std::cout << "Adding contributions of class " << classNumber <<
    //               " to independent parameter " << independentClassNumber << std::endl;
    // std::cout << "     linearSum: " << linearSum << std::endl;
    // std::cout << "     quadraticSum: " << quadraticSum << std::endl;
    // std::cout << "     numberOfPixels: " << numberOfPixels << std::endl;
    if ( independentClassNumber >= independentLinearSum.size() )
    {
      independentLinearSum.resize( independentClassNumber + 1, 0.0f );
      independentQuadraticSum.resize( independentClassNumber + 1, 0.0f );
      independentNumberOfPixels.resize( independentClassNumber + 1, 0.0f );
    }
    independentLinearSum[ independentClassNumber ] += linearSum;
    independentQuadraticSum[ independentClassNumber ] += quadraticSum;
    independentNumberOfPixels[ independentClassNumber ] += numberOfPixels;
  }

  // Calculate independent parameters
  std::vector< float >  independentMeans;
  std::vector< float >  independentVariances;
  float  globalVariance = 0.0f;
  float  globalNumberOfPixels = 0.0f;
  for ( unsigned int independentClassNumber = 0;
        independentClassNumber < independentLinearSum.size();
        independentClassNumber++ )
  {
    const float  linearSum = independentLinearSum[ independentClassNumber ];
    const float  quadraticSum = independentQuadraticSum[ independentClassNumber ];
    const float  numberOfPixels = independentNumberOfPixels[ independentClassNumber ];

    const float  mean = linearSum / numberOfPixels;
    float  variance = quadraticSum / numberOfPixels - pow( mean , 2 );
#if 1
    if ( variance < 1 )
    {
      variance = 1;
    }
#endif

    globalVariance += variance * numberOfPixels;
    globalNumberOfPixels += numberOfPixels;

    // std::cout << "   independent class " << independentClassNumber << ":" << std::endl;
    // std::cout << "          mean " << mean << std::endl;
    // std::cout << "          variance " << variance << std::endl;
    // std::cout << "                linearSum:  " << linearSum << std::endl;
    // std::cout << "                quadraticSum: " << quadraticSum << std::endl;
    // std::cout << "                numberOfPixels: " << numberOfPixels << std::endl;

    independentMeans.push_back( mean );
    independentVariances.push_back( variance );
  }
  globalVariance /= globalNumberOfPixels;


  // Map independent parameters back to original parameter space
  m_Means.clear();
  m_Variances.clear();
  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
  {
    const unsigned int  independentClassNumber = m_IndependentParametersLookupTable[ classNumber ];

    m_Means.push_back( independentMeans[ independentClassNumber ] );
#if 1
    m_Variances.push_back( independentVariances[ independentClassNumber ] );
#else
    m_Variances.push_back( globalVariance );
#endif

  }


  m_PriorWeights.assign( m_NumberOfClasses, 1.0f );
  if ( m_ReestimatedRelativeWeightOfSharedClasses )
  {
    // If a class is sharing the atlas prior with other classes, we allow the relative
    // weight if the different classes (in effect the weight of the mixture model for
    // that class in the atlas) to be re-estimated to something else than just 1/N, with
    // N the number of shared classes. Note that the prior is split as 1/N; therefore the
    // weights are N * relative fractions, to "undo" the 1/N effect.
    std::cout << "Readjusting weight of shared classes" << std::endl;

    // Build inverse lookup table: mapping from original atlas mesh (which we don't have anymore) to
    // the classes we do have
    std::vector< std::vector< unsigned int > >  inverseLookupTable;
    for ( unsigned int finalClassNumber = 0; finalClassNumber < m_LookupTable.size(); finalClassNumber++ )
    {
      // Make sure we have enough memory allocated
      if ( inverseLookupTable.size() < ( m_LookupTable[ finalClassNumber ] + 1 ) )
      {
        inverseLookupTable.resize( m_LookupTable[ finalClassNumber ] + 1 );
      }

      //
      inverseLookupTable[ m_LookupTable[ finalClassNumber ] ].push_back( finalClassNumber );
    }

    for ( unsigned int originalClassNumber = 0;
          originalClassNumber < inverseLookupTable.size();
          originalClassNumber++ )
    {
      const unsigned int  numberOfSharedClasses = inverseLookupTable[ originalClassNumber ].size();

      // If this class was not split, no need to do anything
      if ( numberOfSharedClasses < 2 )
      {
        continue;
      }

      // Count how many voxels were assigned to the originalClassNumber in total
      float  numberOfVoxelsAssignedToOriginalClass = 0;
      for ( unsigned int i = 0; i < numberOfSharedClasses; i++ )
      {
        const unsigned int finalClassNumber = inverseLookupTable[ originalClassNumber ][ i ];
        numberOfVoxelsAssignedToOriginalClass += unnormalizedWeights[ finalClassNumber ];
      }

      // Now update relative weights in the classes split from this originalClassNumber
      std::cout << "Original prior " << originalClassNumber << " is distributed over " << numberOfSharedClasses << " classes." << std::endl;
      for ( unsigned int i = 0; i < numberOfSharedClasses; i++ )
      {
        const unsigned int finalClassNumber = inverseLookupTable[ originalClassNumber ][ i ];
        m_PriorWeights[ finalClassNumber ] = ( unnormalizedWeights[ finalClassNumber ] / numberOfVoxelsAssignedToOriginalClass ) * numberOfSharedClasses;
        std::cout << "    Changed prior weight for " << finalClassNumber << " from 1 to " << m_PriorWeights[ finalClassNumber ] << std::endl;
      }
      std::cout << std::endl;
    } // End loop over all classes in the original, unsplit mesh

  }


}


//
//
//
itk::Array< float >
EMSegmenter
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





//
//
//
float
EMSegmenter
::DoOneEMIterationWithPartialVoluming()
{



  /**
   *
   * Perform one iteration of the EM algorithm when partial voluming is allowed.
   *
   * Because the number of classes (pure or mixed) can quickly become astronimous, it's not
   * feasible to actually first do the E-step, and store the classification, and then do the
   * the M-step explicitly. Instead, here we fold the two steps together, by looping over each
   * voxel, classifying it and adding it's contribution to the quantitaties required for the
   * theoretical M-step.
   *
   */

  // Allocate some variables that we'll be collecting while we go through all voxels and all
  // classes (pure and mixed). These variables will serve as sufficient statistics later on
  // to calculate the M-step
  float  minLogLikelihood = 0.0f;

  itk::Array< float >  totalWeights( m_NumberOfClasses );
  itk::Array< float >  means( m_NumberOfClasses );
  itk::Array< float >  variancesTerm1( m_NumberOfClasses );
  itk::Array< float >  variancesTerm2( m_NumberOfClasses );
  itk::Array< float >  variancesTerm3( m_NumberOfClasses );

  totalWeights.Fill( 0 );
  means.Fill( 0 );
  variancesTerm1.Fill( 0 );
  variancesTerm2.Fill( 0 );
  variancesTerm3.Fill( 0 );


  // Also allocate some temporary variables in which quantities of interest are first stored while
  // finishing the processing of one voxel. This is needed because we can't normalize each voxel's
  // classification before seeing all classes.
  itk::Array< float >  totalWeightsTmp( m_NumberOfClasses );
  itk::Array< float >  meansTmp( m_NumberOfClasses );
  itk::Array< float >  variancesTerm1Tmp( m_NumberOfClasses );
  itk::Array< float >  variancesTerm2Tmp( m_NumberOfClasses );
  itk::Array< float >  variancesTerm3Tmp( m_NumberOfClasses );


  // Loop over all voxels
  std::cout << "Partial volume EM iteration: Starting to loop over all voxels" << std::endl;
  itk::ImageRegionConstIterator< ImageType >  intensityIt( m_Image, m_Image->GetBufferedRegion() );
  itk::ImageRegionConstIterator< SubvoxelPriorsImageType >  subvoxelPriorIt( m_SubvoxelPriorsImage,
      m_SubvoxelPriorsImage->GetBufferedRegion() );
  std::vector< itk::ImageRegionIterator< ClassificationImageType > >   posteriorIts;
  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
  {
    itk::ImageRegionIterator< ClassificationImageType >  posteriorIt( m_Posteriors[ classNumber ],
        m_Posteriors[ classNumber ]->GetBufferedRegion() );
    posteriorIts.push_back( posteriorIt );
  }
  for ( ; !intensityIt.IsAtEnd(); ++intensityIt, ++subvoxelPriorIt )
  {

    const float  y = intensityIt.Value();

    if ( y == 0 )
    {
      // Skip zero voxels, but remember to manually forward the posteriorIts. Wow, this is
      // a real mess -> definitely should have the posteriors stored in one vector-valued
      // image! ...
      for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
      {
        ++( posteriorIts[ classNumber ] );
      }
      continue;
    }

    float likelihood = 0;

    totalWeightsTmp.Fill( 0 );
    meansTmp.Fill( 0 );
    variancesTerm1Tmp.Fill( 0 );
    variancesTerm2Tmp.Fill( 0 );
    variancesTerm3Tmp.Fill( 0 );


    // Contribution of pure tissue classes
    bool  dontBotherLookingFurther = false;
    for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
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
      const float  gauss = exp( -pow( y - mean, 2 ) / 2.0f / variance ) /
                           sqrt( 2 * 3.14 * variance );

      // Add contribution to quantities of interest
      const float  tau = variance * ( m_NumberOfSubvoxels - 1 ) / pow( m_NumberOfSubvoxels, 2 );
      const float  weight = gauss * priorOfSubvoxelCombination;
      meansTmp[ classNumber ] += y / m_NumberOfSubvoxels * weight;
      variancesTerm1Tmp[ classNumber ] += ( tau + pow( y / m_NumberOfSubvoxels, 2 ) ) * weight;
      variancesTerm2Tmp[ classNumber ] += y / m_NumberOfSubvoxels * weight;
      variancesTerm3Tmp[ classNumber ] += weight;
      totalWeightsTmp[ classNumber ] += weight;
      likelihood += weight;

      // Check if we have certainty about this voxel's subvoxels all being this class. If so, then don't
      // waste time calculating zeros - especially since a large part of the high-res image will be background
      // and never be visited by the rasterizor since it falls outside of the mesh's boundary box.
      if ( priorOfSubvoxelCombination > 0.98 )
      {
        //std::cout << "     NOT BOTHERING FURTHER" << std::endl;
        dontBotherLookingFurther = true;
        break;
      }

    }


    if ( !dontBotherLookingFurther )
    {

      // Contribution of PV tissue classes
      int pvClassNumber = -1;
      const int  numberOfGaussiansPerPV = m_NumberOfSubvoxels-1;
      for ( unsigned int classNumber1 = 0; classNumber1 < m_NumberOfClasses; classNumber1++ )
      {
        for ( unsigned int classNumber2 = classNumber1+1; classNumber2 < m_NumberOfClasses; classNumber2++ )
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
            const float  gauss = exp( -pow( y - mean, 2 ) / 2.0f / variance ) /
                                 sqrt( 2 * 3.14 * variance );

            // Add contribution to quantities of interest
            const float  weight = gauss * priorOfSubvoxelCombinations[ gaussianNumber + 1 ];

            // Contribution to class 1
            const float  tau1 = mean1 / m_NumberOfSubvoxels + variance1 / m_NumberOfSubvoxels / variance *
                                ( y - mean );
            const float  sigma1 = ( variance - variance1 / m_NumberOfSubvoxels ) / variance * variance1 / m_NumberOfSubvoxels;
            meansTmp[ classNumber1 ] += weight * alpha * tau1;
            variancesTerm1Tmp[ classNumber1 ] += weight * alpha * ( sigma1 + pow( tau1, 2 ) );
            variancesTerm2Tmp[ classNumber1 ] += weight * alpha * tau1;
            variancesTerm3Tmp[ classNumber1 ] += weight * alpha;
            totalWeightsTmp[ classNumber1 ] += weight * alpha;

            // Contribution to class 2
            const float  tau2 = mean2 / m_NumberOfSubvoxels + variance2 / m_NumberOfSubvoxels / variance *
                                ( y - mean );
            const float  sigma2 = ( variance - variance2 / m_NumberOfSubvoxels ) / variance * variance2 / m_NumberOfSubvoxels;
            meansTmp[ classNumber2 ] += weight * ( 1 - alpha ) * tau2;
            variancesTerm1Tmp[ classNumber2 ] += weight * ( 1 - alpha ) * ( sigma2 + pow( tau2, 2 ) );
            variancesTerm2Tmp[ classNumber2 ] += weight * ( 1 - alpha ) * tau2;
            variancesTerm3Tmp[ classNumber2 ] += weight * ( 1 - alpha );
            totalWeightsTmp[ classNumber2 ] += weight * ( 1 - alpha );

            likelihood += weight;
          } // End loop over all sub-gaussians of this PV class

        }
      } // End looping over PV tissue classes

    } // End test if we need to loop over non-pure classes or not

    likelihood += 1e-15;
    minLogLikelihood -= log( likelihood );


    // Now that we have calculated the likelihood of this voxel, we can normalize the contributions (in the previous
    // equations, "weight" should have been, but was not, normalized) and add them with those of previous voxels
    for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
    {
      means[ classNumber ] += meansTmp[ classNumber ] / likelihood;
      variancesTerm1[ classNumber ] += variancesTerm1Tmp[ classNumber ] / likelihood;
      variancesTerm2[ classNumber ] += variancesTerm2Tmp[ classNumber ] / likelihood;
      variancesTerm3[ classNumber ] += variancesTerm3Tmp[ classNumber ] / likelihood;
      totalWeights[ classNumber ] += totalWeightsTmp[ classNumber ] / likelihood;
    }


    // For visualization purposes for the user, also calculate posterior expected fraction in this voxel
    for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
    {
      ( posteriorIts[ classNumber ] ).Value() = totalWeightsTmp[ classNumber ] / likelihood;
      ++( posteriorIts[ classNumber ] );
    }

  }  // End loop over all low-resolution voxels


  // Now update the means and variances; keep in mind that some classes' parameters are coupled
  std::cout << "Partial volume EM iteration: finished looping over all voxels" << std::endl;

  // Add contributions of classes that are coupled together
  std::vector< float >   independentMeanszz;
  std::vector< float >   independentVariancesTerm1;
  std::vector< float >   independentVariancesTerm2;
  std::vector< float >   independentVariancesTerm3;
  std::vector< float >   independentTotalWeights;
  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
  {
    // Add contributions to the correct independent class; first make sure we have enough memory allocated
    const unsigned int  independentClassNumber = m_IndependentParametersLookupTable[ classNumber ];
    //std::cout << "Adding contributions of class " << classNumber <<
    //              " to independent parameter " << independentClassNumber << std::endl;

    if ( independentClassNumber >= independentMeanszz.size() )
    {
      independentMeanszz.resize( independentClassNumber + 1, 0.0f );
      independentVariancesTerm1.resize( independentClassNumber + 1, 0.0f );
      independentVariancesTerm2.resize( independentClassNumber + 1, 0.0f );
      independentVariancesTerm3.resize( independentClassNumber + 1, 0.0f );
      independentTotalWeights.resize( independentClassNumber + 1, 0.0f );
    }
    independentMeanszz[ independentClassNumber ] += means[ classNumber ];
    independentVariancesTerm1[ independentClassNumber ] += variancesTerm1[ classNumber ];
    independentVariancesTerm2[ independentClassNumber ] += variancesTerm2[ classNumber ];
    independentVariancesTerm3[ independentClassNumber ] += variancesTerm3[ classNumber ];
    independentTotalWeights[ independentClassNumber ] += totalWeights[ classNumber ];
  }


  // Calculate independent parameters
  std::vector< float >  independentMeans;
  std::vector< float >  independentVariances;
  for ( unsigned int independentClassNumber = 0;
        independentClassNumber < independentMeanszz.size();
        independentClassNumber++ )
  {
    const float  mean = m_NumberOfSubvoxels * independentMeanszz[ independentClassNumber ] /
                        independentTotalWeights[ independentClassNumber ];
    const float  variance = m_NumberOfSubvoxels * ( independentVariancesTerm1[ independentClassNumber ]
                            - 2 * mean / m_NumberOfSubvoxels *
                            independentVariancesTerm2[ independentClassNumber ]
                            + pow( mean / m_NumberOfSubvoxels, 2 )
                            * independentVariancesTerm3[ independentClassNumber ] )
                            / independentTotalWeights[ independentClassNumber ];

    independentMeans.push_back( mean );
    independentVariances.push_back( variance );
  }

  // Map independent parameters back to original parameter space
  m_Means.clear();
  m_Variances.clear();
  for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
  {
    const unsigned int  independentClassNumber = m_IndependentParametersLookupTable[ classNumber ];

    m_Means.push_back( independentMeans[ independentClassNumber ] );
    m_Variances.push_back( independentVariances[ independentClassNumber ] );
  }


  return minLogLikelihood;

}






//
//
//
EMSegmenter::BiasCorrectedImageType::Pointer
EMSegmenter
::BiasCorrect( const ImageType* image, int upsamplingFactor ) const
{

  // Construct upsampled bias field
  typedef ImageSmoother::MaskImageType  MaskImageType;
  MaskImageType::Pointer  maskImage = MaskImageType::New();
  maskImage->SetRegions( image->GetBufferedRegion() );
  maskImage->Allocate();

  itk::ImageRegionConstIterator< ImageType >  imageIt( image, image->GetBufferedRegion() );
  itk::ImageRegionIterator< MaskImageType >  maskIt( maskImage, maskImage->GetBufferedRegion() );
  for ( ; !imageIt.IsAtEnd(); ++imageIt, ++maskIt )
  {
    if ( imageIt.Value() == 0 )
    {
      maskIt.Value() = false;
    }
    else
    {
      maskIt.Value() = true;
    }

  }

  BiasFieldImageType::ConstPointer  biasField = m_BiasFieldResidueSmoother->GetSmoothedImage( maskImage, upsamplingFactor ).GetPointer();

#if 0
  // Compute the region that overlaps in both the image and the bias field
  ImageType::IndexType  overlappingRegionIndex = { 0, 0, 0 };
  ImageType::SizeType  overlappingRegionSize;
  for ( int i = 0; i < 3; i++ )
  {
    overlappingRegionSize[ i ] = image->GetBufferedRegion().GetSize()[ i ];
    if ( overlappingRegionSize[ i ] >  biasField->GetBufferedRegion().GetSize()[ i ] )
    {
      overlappingRegionSize[ i ] =  biasField->GetBufferedRegion().GetSize()[ i ];
    }

  }
  ImageType::RegionType  overlappingRegion( overlappingRegionIndex, overlappingRegionSize );
#endif

  // Create an empty bias corrected image
  BiasCorrectedImageType::Pointer  biasCorrectedImage = BiasCorrectedImageType::New();
  biasCorrectedImage->SetRegions( image->GetBufferedRegion() );
  biasCorrectedImage->Allocate();
  biasCorrectedImage->FillBuffer( 0 );

  // Calculate bias field corrected intensities in the overlapping area
  //itk::ImageRegionConstIterator< ImageType >  imageIt( image, image->GetBufferedRegion() );
  itk::ImageRegionConstIterator< BiasFieldImageType >  biasIt( biasField, biasField->GetBufferedRegion() );
  itk::ImageRegionIterator< BiasCorrectedImageType >  correctedIt( biasCorrectedImage, biasCorrectedImage->GetBufferedRegion() );
  for ( imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, ++correctedIt, ++biasIt )
  {
    if ( imageIt.Value() == 0 )
    {
      continue;
    }

    correctedIt.Value() = imageIt.Value() - biasIt.Value();
  }


  //
  return biasCorrectedImage;

}




//
//
//
EMSegmenter::ImageType::Pointer
EMSegmenter
::ConvertToNativeImageType( const BiasCorrectedImageType* image )
{

  // Convert it to native pixel type, simply by rounding. Make sure we don't exceed range of what can be represented
  ImageType::Pointer  nativeImage = ImageType::New();
  nativeImage->SetRegions( image->GetBufferedRegion() );
  nativeImage->Allocate();
  itk::ImageRegionConstIterator< BiasCorrectedImageType >  nonNativeIt( image,
      image->GetBufferedRegion() );
  itk::ImageRegionIterator< ImageType >  nativeIt( nativeImage,
      nativeImage->GetBufferedRegion() );
  for ( ; !nonNativeIt.IsAtEnd(); ++nonNativeIt, ++nativeIt )
  {
    if ( ( nonNativeIt.Value() + 0.5 ) > itk::NumericTraits< ImageType::PixelType >::max() )
    {
      nativeIt.Value() = itk::NumericTraits< ImageType::PixelType >::max();
    }
    else if ( ( nonNativeIt.Value() + 0.5 ) < itk::NumericTraits< ImageType::PixelType >::min() )
    {
      nativeIt.Value() = itk::NumericTraits< ImageType::PixelType >::min();
    }
    else
    {
      nativeIt.Value() = static_cast< ImageType::PixelType >( nonNativeIt.Value() + 0.5 );
    }
  }

  //
  return nativeImage;

}




//
//
//
void
EMSegmenter
::SetImageAndSegmentItWithCurrentParametersAndWithoutBiasFieldCorrection( const ImageType* image )
{
  std::cout << "-- Starting" << std::endl;
  this->SetImage( image );
  std::cout << "-- Initializing" << std::endl;
  this->Initialize();
  std::cout << "-- Doing the rest" << std::endl;
  if ( m_NumberOfSubvoxels == 1 )
  {
    std::cout << "-- Starting to update posteriors" << std::endl;
    this->UpdatePosteriors();
    std::cout << "-- done" << std::endl;
  }
  else
  {
    this->DoOneEMIterationWithPartialVoluming();
  }

}




//
//
//
EMSegmenter::SuperResolutionProbabilityImageType::Pointer
EMSegmenter
::GetSuperResolutionLikelihoods() const
{

  // Sanity check
  if ( !m_SuperResolutionPriors )
  {
    return 0;
  }

  // Allocate an image filled with zeroes for the likelihoods
  SuperResolutionProbabilityImageType::Pointer  superResolutionLikelihoods = SuperResolutionProbabilityImageType::New();
  superResolutionLikelihoods->SetRegions( m_SuperResolutionPriors->GetLargestPossibleRegion().GetSize() );
  superResolutionLikelihoods->Allocate();
  SuperResolutionProbabilityImageType::PixelType  allZeroes( m_NumberOfClasses );
  allZeroes.Fill( 0.0f );
  superResolutionLikelihoods->FillBuffer( allZeroes );


  // We want to know, in every subvoxel (i.e., voxel in the super resolution image), what the the likelihood
  // is of the big average voxel's intensity *given* the label in the subvoxel. This is given by summing over
  //  all configurations of the subvoxel's siblings, i.e. voxels in the superresolution image underlying the
  // same averge voxel. Summing over all these siblings is implemented more efficiently by grouping all the
  // siblings' configurations that result in the same fraction, and summing over all possible fractions instead.

  // Loop over all voxels
  std::cout << "Likelihood calculation at super-resolution: Starting to loop over all voxels" << std::endl;
  itk::ImageRegionConstIteratorWithIndex< ImageType >  intensityIt( m_Image, m_Image->GetBufferedRegion() );
  itk::ImageRegionConstIterator< SubvoxelPriorsImageType >  subvoxelPriorIt( m_SubvoxelPriorsImage,
      m_SubvoxelPriorsImage->GetBufferedRegion() );
  for ( ; !intensityIt.IsAtEnd(); ++intensityIt, ++subvoxelPriorIt )
  {

    const float  y = intensityIt.Value();


    // Contribution of pure tissue classes
    bool  dontBotherLookingFurther = false;
    for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
    {
      // Precalculate stuff of interest
      const float  mean = m_Means[ classNumber ];
      const float  variance = m_Variances[ classNumber ];
      const float  gauss = exp( -pow( y - mean, 2 ) / 2.0f / variance ) /
                           sqrt( 2 * 3.14 * variance );



      // Loop over all subvoxels underlying the current voxel
      ImageType::IndexType  index;  // Will hold the index in the superresolution image
      int  counter = -1;  // Will hold the index in the container holding the super resolution priors for the current voxel
      for ( int xStep = 0; xStep < m_PartialVolumeUpsamplingFactors[ 0 ]; xStep++ )
      {
        index[ 0 ] = intensityIt.GetIndex()[ 0 ] * m_PartialVolumeUpsamplingFactors[ 0 ] + xStep;
        for ( int yStep = 0; yStep < m_PartialVolumeUpsamplingFactors[ 1 ]; yStep++ )
        {
          index[ 1 ] = intensityIt.GetIndex()[ 1 ] * m_PartialVolumeUpsamplingFactors[ 1 ] + yStep;
          for ( int zStep = 0; zStep < m_PartialVolumeUpsamplingFactors[ 2 ]; zStep++ )
          {
            index[ 2 ] = intensityIt.GetIndex()[ 2 ] * m_PartialVolumeUpsamplingFactors[ 2 ] + zStep;
            counter++;

            // Calculate prior that all siblings are this class
            float  priorOfSiblingCombination = 1;
            for ( int subvoxelNumber = 0; subvoxelNumber < m_NumberOfSubvoxels; subvoxelNumber++ )
            {
              if ( subvoxelNumber == counter )
              {
                continue;
              }

              priorOfSiblingCombination *= ( *( subvoxelPriorIt.Value()[ subvoxelNumber ] ) )[ classNumber ];
            }


            // Add contribution to likelihood
            superResolutionLikelihoods->GetPixel( index )[ classNumber ] += gauss * priorOfSiblingCombination;
          }
        }
      } // End loop over all subvoxels underlying the current voxel



      // Check if we have certainty about this voxel's subvoxels all being this class. If so, then don't
      // waste time calculating zeros - especially since a large part of the high-res image will be background
      // and never be visited by the rasterizor since it falls outside of the mesh's boundary box.
      float  priorOfSubvoxelCombination = 1.0f;
      for ( int subvoxelNumber = 0; subvoxelNumber < m_NumberOfSubvoxels; subvoxelNumber++ )
      {
        priorOfSubvoxelCombination *= ( *( subvoxelPriorIt.Value()[ subvoxelNumber ] ) )[ classNumber ];
      }
      if ( priorOfSubvoxelCombination > 0.98 )
      {
        //std::cout << "     NOT BOTHERING FURTHER" << std::endl;
        dontBotherLookingFurther = true;
        break;
      }

    } // End loop over all pure classes


    if ( dontBotherLookingFurther )
    {
      continue;
    }


    // Contribution of PV tissue classes
    int pvClassNumber = -1;
    const int  numberOfGaussiansPerPV = m_NumberOfSubvoxels-1;
    for ( unsigned int classNumber1 = 0; classNumber1 < m_NumberOfClasses; classNumber1++ )
    {
      for ( unsigned int classNumber2 = classNumber1+1; classNumber2 < m_NumberOfClasses; classNumber2++ )
      {
        pvClassNumber++;

        // Retrieve parameters of the pure tissue classes this PV class is mixing for
        const float  mean1 = m_Means[ classNumber1 ];
        const float  variance1 = m_Variances[ classNumber1 ];
        const float  mean2 = m_Means[ classNumber2 ];
        const float  variance2 = m_Variances[ classNumber2 ];


        // Pre-calculate gauss value for each mixing fraction
        std::vector< float >  gausses( numberOfGaussiansPerPV );
        for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussiansPerPV; gaussianNumber++ )
        {
          // Retrieve mean and variance for this sub-Gaussian
          const float  alpha = static_cast< float >( gaussianNumber + 1 ) / static_cast< float >( m_NumberOfSubvoxels );
          const float  mean = alpha * mean1 + ( 1 - alpha ) * mean2;
          const float  variance = alpha * variance1 + ( 1 - alpha ) * variance2;

          // Calculate gauss
          const float  gauss = exp( -pow( y - mean, 2 ) / 2.0f / variance ) /
                               sqrt( 2 * 3.14 * variance );

          // Store for later usage
          gausses[ gaussianNumber ] = gauss;
        }


        // Loop over all subvoxels underlying the current voxel
        ImageType::IndexType  index;  // Will hold the index in the superresolution image
        int  counter = -1;  // Will hold the index in the container holding the super resolution priors for the current voxel
        for ( int xStep = 0; xStep < m_PartialVolumeUpsamplingFactors[ 0 ]; xStep++ )
        {
          index[ 0 ] = intensityIt.GetIndex()[ 0 ] * m_PartialVolumeUpsamplingFactors[ 0 ] + xStep;
          for ( int yStep = 0; yStep < m_PartialVolumeUpsamplingFactors[ 1 ]; yStep++ )
          {
            index[ 1 ] = intensityIt.GetIndex()[ 1 ] * m_PartialVolumeUpsamplingFactors[ 1 ] + yStep;
            for ( int zStep = 0; zStep < m_PartialVolumeUpsamplingFactors[ 2 ]; zStep++ )
            {
              index[ 2 ] = intensityIt.GetIndex()[ 2 ] * m_PartialVolumeUpsamplingFactors[ 2 ] + zStep;
              counter++;


              // Get prior probabilities for mixing fractions generated by this subvoxel's siblings
              itk::Array< float >  pureProbabilitiesOfClass1( m_NumberOfSubvoxels - 1 );
              itk::Array< float >  pureProbabilitiesOfClass2( m_NumberOfSubvoxels - 1 );

              for ( int subvoxelNumber = 0, siblingNumber = 0; subvoxelNumber < m_NumberOfSubvoxels; subvoxelNumber++ )
              {
                if ( subvoxelNumber == counter )
                {
                  continue;
                }

                pureProbabilitiesOfClass1[ siblingNumber ] = ( *( subvoxelPriorIt.Value()[ subvoxelNumber ] ) )[ classNumber1 ];
                pureProbabilitiesOfClass2[ siblingNumber ] = ( *( subvoxelPriorIt.Value()[ subvoxelNumber ] ) )[ classNumber2 ];
                siblingNumber++;
              }
              itk::Array< float >  priorOfSiblingCombinations = this->GetMixingProbabilities( pureProbabilitiesOfClass1,
                  pureProbabilitiesOfClass2 );



              // Loop over all sub-gaussians
              for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussiansPerPV; gaussianNumber++ )
              {
                const float  gauss = gausses[ gaussianNumber ];

                // Add contribution to likelihood
                superResolutionLikelihoods->GetPixel( index )[ classNumber1 ] += gauss * priorOfSiblingCombinations[ gaussianNumber ];
                superResolutionLikelihoods->GetPixel( index )[ classNumber2 ] += gauss * priorOfSiblingCombinations[ gaussianNumber + 1 ];

              } // End loop over all sub-gaussians


            }
          }
        } // End loop over all subvoxels underlying the current voxel



      }
    } // End looping over PV tissue classes


  } // End loop over all voxels


  return superResolutionLikelihoods;
}



//
//
//
EMSegmenter::SuperResolutionProbabilityImageType::Pointer
EMSegmenter
::GetSuperResolutionPosteriors() const
{

  // Sanity check
  if ( !m_SuperResolutionPriors )
  {
    return 0;
  }

  std::cout << "Entering main body of EMSegmenter::GetSuperResolutionPosteriors()" << std::endl;

  // Initialize posterior with likelihoods
  SuperResolutionProbabilityImageType::Pointer  superResolutionPosteriors = this->GetSuperResolutionLikelihoods();

  // Loop over all subvoxels (i.e., super resolution voxels), multipy likelihood by prior, and normalize
  itk::ImageRegionConstIterator< SuperResolutionProbabilityImageType >  priorIt( m_SuperResolutionPriors,
      m_SuperResolutionPriors->GetBufferedRegion() );
  itk::ImageRegionIterator< SuperResolutionProbabilityImageType >  posteriorIt( superResolutionPosteriors,
      superResolutionPosteriors->GetBufferedRegion() );
  for ( ; !priorIt.IsAtEnd(); ++priorIt, ++posteriorIt )
  {

    float  sum = 1e-15f;
    for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
    {
      posteriorIt.Value()[ classNumber ] *= priorIt.Value()[ classNumber ];
      sum += posteriorIt.Value()[ classNumber ];
    }
    for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
    {
      posteriorIt.Value()[ classNumber ] /= sum;
    }

  }


  // Return result
  return  superResolutionPosteriors;

}



//
//
//
EMSegmenter::ImageType::Pointer
EMSegmenter
::GetSuperResolutionImage() const
{

  // Create an empty super resolution image
  ImageType::SizeType  superResolutionSize;
  for ( int i = 0; i < 3; i++ )
  {
    superResolutionSize[ i ] = m_Image->GetBufferedRegion().GetSize()[ i ] * m_PartialVolumeUpsamplingFactors[ i ];
  }
  ImageType::Pointer  superResolutionImage = ImageType::New();
  superResolutionImage->SetRegions( superResolutionSize );
  superResolutionImage->Allocate();


  // Loop over voxels in the original image, and simply copy the intensity to the corresponding
  // subvoxels in the super resolution image
  for ( itk::ImageRegionConstIteratorWithIndex< ImageType >  it( m_Image,
        m_Image->GetBufferedRegion() );
        !it.IsAtEnd(); ++it )
  {
    //std::cout << "Filling in intensity in subvoxels of big voxel with index " << it.GetIndex() << std::endl;

    ImageType::IndexType  index; // Will hold the index in the super resolution image
    for ( int xStep = 0; xStep < m_PartialVolumeUpsamplingFactors[ 0 ]; xStep++ )
    {
      index[ 0 ] = it.GetIndex()[ 0 ] * m_PartialVolumeUpsamplingFactors[ 0 ] + xStep;
      for ( int yStep = 0; yStep < m_PartialVolumeUpsamplingFactors[ 1 ]; yStep++ )
      {
        index[ 1 ] = it.GetIndex()[ 1 ] * m_PartialVolumeUpsamplingFactors[ 1 ] + yStep;
        for ( int zStep = 0; zStep < m_PartialVolumeUpsamplingFactors[ 2 ]; zStep++ )
        {
          index[ 2 ] = it.GetIndex()[ 2 ] * m_PartialVolumeUpsamplingFactors[ 2 ] + zStep;

          //std::cout << "         Copying intensity to high-res voxel with index " << index << std::endl;
          superResolutionImage->SetPixel( index, it.Value() );
        }
      }
    }

  }


  // Return result
  return  superResolutionImage;

}




//
//
//
AtlasMesh::Pointer
EMSegmenter
::GetTransformedMesh( const AtlasMesh* mesh, const float* translation, const float* scaling )
{
  std::cout << "Transforming mesh" << std::endl;

  // Transform mesh points
  AtlasMesh::PointsContainer::Pointer  transformedPoints = AtlasMesh::PointsContainer::New();
  for ( AtlasMesh::PointsContainer::ConstIterator it = mesh->GetPoints()->Begin();
        it != mesh->GetPoints()->End(); ++it )
  {
    AtlasMesh::PointType  transformedPoint;
    for ( int i = 0; i < 3; i++ )
    {
      transformedPoint[ i ] = scaling[ i ] * it.Value()[ i ] + translation[ i ];
    }

    transformedPoints->InsertElement( it.Index(), transformedPoint );
  }


  // Also the reference position of the mesh has changed. Note, however, that we don't
  // actually have access to the original reference position, only some sufficient
  // statistics calculated from it. In particular, we have only the three first columns
  // of the matrix Z = inv( X ) = inv( [ p0 p1 p2 p3; 1 1 1 1 ] ). The transformed
  // position Xtrans is given by
  //
  //     Xtrans = diag( scaling[0] scaling[ 1 ] scaling[ 2 ] 1 ) * X + [ translation[ 0 ] translation[ 1 ] translation[ 2 ] 1 ]'
  //
  // but since we'll also end up calculating the upper 3x3 matrix of Ytrans * Ztrans
  // (with Y the equivalent of X but in the deformed mesh)
  // to evaluate the deformation penatly, we can safely drop the translation part.
  // In short, we will calculate the first 3 columns of the matrix
  //
  //    Ztrans = inv( diag( scaling[0] scaling[ 1 ] scaling[ 2 ] 1 ) * X )
  //
  //           = Z * diag( 1/scaling[0] 1/scaling[1] 1/scaling[2] 1 )
  //
  // which is given by multiplying each column i of Z with a factor 1/scaling[i]
  //
  AtlasMesh::CellDataContainer::Pointer   transformedReferenceTetrahedronInfos = AtlasMesh::CellDataContainer::New();
  for ( AtlasMesh::CellDataContainer::ConstIterator  it = mesh->GetCellData()->Begin();
        it != mesh->GetCellData()->End(); ++it )
  {
    ReferenceTetrahedronInfo  info;

    info.m_ReferenceVolumeTimesK = it.Value().m_ReferenceVolumeTimesK;

    info.m_Z11 = it.Value().m_Z11 / scaling[ 0 ];
    info.m_Z21 = it.Value().m_Z21 / scaling[ 0 ];
    info.m_Z31 = it.Value().m_Z31 / scaling[ 0 ];
    info.m_Z41 = it.Value().m_Z41 / scaling[ 0 ];

    info.m_Z12 = it.Value().m_Z12 / scaling[ 1 ];
    info.m_Z22 = it.Value().m_Z22 / scaling[ 1 ];
    info.m_Z32 = it.Value().m_Z32 / scaling[ 1 ];
    info.m_Z42 = it.Value().m_Z42 / scaling[ 1 ];

    info.m_Z13 = it.Value().m_Z13 / scaling[ 2 ];
    info.m_Z23 = it.Value().m_Z23 / scaling[ 2 ];
    info.m_Z33 = it.Value().m_Z33 / scaling[ 2 ];
    info.m_Z43 = it.Value().m_Z43 / scaling[ 2 ];


    transformedReferenceTetrahedronInfos->InsertElement( it.Index(), info );

#if 0
    {
      if ( it.Index() == 82444 )
      {
        std::cout << "info.m_ReferenceVolumeTimesK: " << info.m_ReferenceVolumeTimesK << std::endl;
        std::cout << std::endl;
        std::cout << "info.m_Z11: " << info.m_Z11 << std::endl;
        std::cout << "info.m_Z21: " << info.m_Z21 << std::endl;
        std::cout << "info.m_Z31: " << info.m_Z31 << std::endl;
        std::cout << "info.m_Z41: " << info.m_Z41 << std::endl;
        std::cout << std::endl;
        std::cout << "info.m_Z12: " << info.m_Z12 << std::endl;
        std::cout << "info.m_Z22: " << info.m_Z22 << std::endl;
        std::cout << "info.m_Z32: " << info.m_Z32 << std::endl;
        std::cout << "info.m_Z42: " << info.m_Z42 << std::endl;
        std::cout << std::endl;
        std::cout << "info.m_Z13: " << info.m_Z13 << std::endl;
        std::cout << "info.m_Z23: " << info.m_Z23 << std::endl;
        std::cout << "info.m_Z33: " << info.m_Z33 << std::endl;
        std::cout << "info.m_Z43: " << info.m_Z43 << std::endl;
      }
    }
#endif
  }


  // Put together the pieces of the mesh
  AtlasMesh::Pointer  transformedMesh = AtlasMesh::New();
  transformedMesh->SetPoints( transformedPoints );
  transformedMesh->SetCells( const_cast< AtlasMesh::CellsContainer* >( mesh->GetCells() ) );
  transformedMesh->SetPointData( const_cast< AtlasMesh::PointDataContainer* >( mesh->GetPointData() ) );
  transformedMesh->SetCellData( transformedReferenceTetrahedronInfos );


  // Return result
  std::cout << "Done transforming mesh" << std::endl;
  return  transformedMesh;

}



//
//
//
AtlasMesh::Pointer
EMSegmenter
::GetSuperResolutionMesh( const AtlasMesh* mesh ) const
{

  // Calculate the translation and scaling of mesh points in super resolution
  float  translation[ 3 ];
  float  scaling[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    translation[ i ] = ( m_PartialVolumeUpsamplingFactors[ i ] - 1 ) / 2.0f;
    scaling[ i ] = m_PartialVolumeUpsamplingFactors[ i ];
  }

  // Calculate result
  return Self::GetTransformedMesh( mesh, translation, scaling );

}



//
//
//
AtlasMesh::Pointer
EMSegmenter
::GetNormalResolutionMesh( const AtlasMesh* superResolutionMesh ) const
{
  // Calculate the translation and scaling of mesh points in normal resolution
  float  translation[ 3 ];
  float  scaling[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    scaling[ i ] = 1.0f / m_PartialVolumeUpsamplingFactors[ i ];

    translation[ i ] = -( ( m_PartialVolumeUpsamplingFactors[ i ] - 1 ) / 2.0f ) / m_PartialVolumeUpsamplingFactors[ i ];
  }

  // Calculate result
  return Self::GetTransformedMesh( superResolutionMesh, translation, scaling );

}



} // end namespace kvl
