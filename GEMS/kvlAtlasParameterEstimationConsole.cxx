/**
 * @file  kvlAtlasParameterEstimationConsole.cxx
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
#include "kvlAtlasParameterEstimationConsole.h"

#include "FL/Fl.H"
#include "FL/fl_ask.H"
#include "itkImageFileReader.h"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "kvlAtlasMeshSummaryDrawer.h"
#include "kvlCompressionLookupTable.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkCommand.h"
#include "kvlAtlasMeshCollectionModelLikelihoodCalculator.h"
#include "kvlAtlasMeshCollectionPositionCostCalculator2.h"
#include "itkCastImageFilter.h"




namespace kvl
{

//
//
//
AtlasParameterEstimationConsole
::AtlasParameterEstimationConsole()
{

  m_Estimator = AtlasParameterEstimator::New();

#if 0
  std::cout << "Setting maximum number of iterations to 10 for debugging purposes. Please remove this!" << std::endl;
  m_Estimator->SetMaximumNumberOfIterations( 10 );
#endif

  m_Interrupted = false;
  m_Stepping = false;

  // Add observers to the estimator
  typedef itk::MemberCommand< AtlasParameterEstimationConsole >   MemberCommandType;
  MemberCommandType::Pointer  command = MemberCommandType::New();
  command->SetCallbackFunction( this, &AtlasParameterEstimationConsole::HandleEstimatorEvent );
  m_Estimator->AddObserver( itk::StartEvent(), command );
  m_Estimator->AddObserver( itk::IterationEvent(), command );
  m_Estimator->AddObserver( itk::EndEvent(), command );
  m_Estimator->AddObserver( AlphasEstimationStartEvent(), command );
  m_Estimator->AddObserver( AlphasEstimationIterationEvent(), command );
  m_Estimator->AddObserver( AlphasEstimationEndEvent(), command );
  m_Estimator->AddObserver( PositionEstimationStartEvent(), command );
  m_Estimator->AddObserver( PositionEstimationIterationEvent(), command );
  m_Estimator->AddObserver( PositionEstimationEndEvent(), command );

  // Set up GUI to reflect estimator state
  m_PositionEstimationResolution->value( m_Estimator->GetPositionEstimationIterationEventResolution() );

  //
  m_NumberOfUpsamplingSteps->do_callback();
}



//
//
//
AtlasParameterEstimationConsole
::~AtlasParameterEstimationConsole()
{

}



//
//
//
void
AtlasParameterEstimationConsole
::SetLabelImages( const std::vector< std::string >&  fileNames, bool useGaussians, bool ignoreLastImage )
{
  // Sanity checking
  if ( fileNames.size() == 0 )
  {
    return;
  }




  // Read the images (in original ushort format)
  typedef CompressionLookupTable::ImageType  InputImageType;
  std::vector< InputImageType::ConstPointer >  originalImages;
  for ( std::vector< std::string >::const_iterator it = fileNames.begin();
        it != fileNames.end(); ++it )
  {
    // Read the input image
    typedef itk::ImageFileReader< InputImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( ( *it ).c_str() );
    reader->Update();
    InputImageType::ConstPointer  originalImage = reader->GetOutput();

    // Over-ride the spacing and origin since at this point we can't deal with that
    const double spacing[] = { 1, 1, 1 };
    const double origin[] = { 0, 0, 0 };
    const_cast< InputImageType* >( originalImage.GetPointer() )->SetSpacing( spacing );
    const_cast< InputImageType* >( originalImage.GetPointer() )->SetOrigin( origin );

    // Remember this image
    originalImages.push_back( originalImage );
  }


  // Build up the internal label images (in uchar - whereas the original images are in ushort)
  typedef CompressionLookupTable::CompressedImageType  OutputImageType;
  CompressionLookupTable::Pointer  compressor = CompressionLookupTable::New();
  std::vector< OutputImageType::ConstPointer >  labelImages;
  if ( !useGaussians )
  {
    // Build a lookup table that maps the original intensities onto uchar starting
    // at 0 and densely packed
    compressor->Construct( originalImages );
    compressor->Write( "compressionLookupTable.txt" );

    // Collect the label images resulting from pushing the original images through the
    // lookup table
    for ( std::vector< InputImageType::ConstPointer >::const_iterator it = originalImages.begin();
          it != originalImages.end(); ++it )
    {
      labelImages.push_back( ( compressor->CompressImage( ( *it ).GetPointer() ) ).GetPointer() );
    }
  }
  else
  {
    // Just convert the original images into uchar
    for ( std::vector< InputImageType::ConstPointer >::const_iterator it = originalImages.begin();
          it != originalImages.end(); ++it )
    {
      typedef itk::CastImageFilter< InputImageType, OutputImageType >   CasterType;
      CasterType::Pointer  caster = CasterType::New();
      caster->SetInput( ( *it ).GetPointer() );
      caster->Update();
      labelImages.push_back( caster->GetOutput() );
    }

  }


  // Give the label images to the estimator
  m_Estimator->SetLabelImages( labelImages );
  m_Estimator->SetUseGaussians( useGaussians );
  m_Estimator->SetIgnoreLastLabelImage( ignoreLastImage );

  // Construct the initial mesh
  this->InitializeMesh();

  // Change the GUI according to the new images
  for ( unsigned int i = 0; i < labelImages.size(); i++ )
  {
    std::ostringstream  labelStream;
    labelStream << "Image " << i;
    m_LabelImageNumber->add( labelStream.str().c_str() );
  }

  for ( unsigned int i = 0; i < m_Estimator->GetNumberOfClasses(); i++ )
  {
    std::ostringstream  labelStream;
    //if ( labelStringLookupTable.find( static_cast< OutputImageType::PixelType >( i ) ) != labelStringLookupTable.end() )
    if ( compressor->GetLabelStringLookupTable().find( i ) != compressor->GetLabelStringLookupTable().end() )
    {
      // Found an anatomical name in our lookup table. Provide that to the user
      labelStream << ( *( compressor->GetLabelStringLookupTable().find( i ) ) ).second;
    }
    else
    {
      // Didn't find an anatomical label; just provide a silly label name then
      labelStream << "Label " << i;
    }
    m_LabelNumber->add( labelStream.str().c_str() );
  }
  m_LabelNumber->value( 0 );

  // Set the GUI to reflect the first image
  m_LabelImageNumber->value( 0 );
  m_LabelImageNumber->do_callback();

}



//
//
//
void
AtlasParameterEstimationConsole
::InitializeMesh()
{
  // Don't build a mesh if no label images have been set
  if ( m_Estimator->GetNumberOfLabelImages() == 0 )
  {
    return;
  }

  // Construct mesh
  if ( m_UseExplicitStartCollection->value() )
  {
    // We are using an explicit start collection, read it from file and start the estimator with it
    AtlasMeshCollection::Pointer  explicitStartCollection = AtlasMeshCollection::New();
    if ( !explicitStartCollection->Read( m_ExplicitStartCollection->value() ) )
    {
      fl_alert( "Couldn't read explicit start mesh collection from file!" );
      return;
    }

    m_Estimator->SetInitialMeshCollection( explicitStartCollection );
  }
  else
  {

    // Retrieve mesh parameters from the GUI
    unsigned int  meshSize[ 3 ];
    meshSize[ 0 ] = static_cast< unsigned int >( m_Size0->value() );
    meshSize[ 1 ] = static_cast< unsigned int >( m_Size1->value() );
    meshSize[ 2 ] = static_cast< unsigned int >( m_Size2->value() );

    // Construct a mesh collection in Talairach space
    AtlasParameterEstimator::LabelImageType::SizeType  labelImageSize =
      m_Estimator->GetLabelImage( 0 )->GetLargestPossibleRegion().GetSize();
    unsigned int  domainSize[ 3 ];
    domainSize[ 0 ] = static_cast< unsigned int >( labelImageSize[ 0 ] );
    domainSize[ 1 ] = static_cast< unsigned int >( labelImageSize[ 1 ] );
    domainSize[ 2 ] = static_cast< unsigned int >( labelImageSize[ 2 ] );



    unsigned int  numberOfClasses = m_Estimator->GetNumberOfClasses();

    float  initialStiffness = m_InitialStiffness->value();

    unsigned int  numberOfMeshes = m_Estimator->GetNumberOfLabelImages();

    // Create a mesh collection accordingly, and initialize the estimator with it
    AtlasMeshCollection::Pointer  meshCollection = AtlasMeshCollection::New();
    meshCollection->Construct( meshSize, domainSize, initialStiffness,
                               numberOfClasses, numberOfMeshes );


#if 0
    typedef AtlasMeshCollection::TransformType  TransformType;
    TransformType::Pointer  transform = TransformType::New();
    TransformType::OutputVectorType  compensationScale;
    compensationScale[ 0 ] = 1.0f / domainSize[ 0 ];
    compensationScale[ 1 ] = 1.0f / domainSize[ 1 ];
    compensationScale[ 2 ] = 1.0f / domainSize[ 2 ];
    transform->Scale( compensationScale );

    // OK, so now we have our mesh defined as ( 0 0 0) -> ( 1 1 1 ). Now specify how it maps
    // into the label images' voxel grid
    TransformType::OutputVectorType  scale;
    scale[ 0 ] = 60; // Size in number of voxels
    scale[ 1 ] = 60; // Size in number of voxels
    scale[ 2 ] = 60; // Size in number of voxels
    transform->Scale( scale );
    TransformType::OutputVectorType  translation;
    translation[ 0 ] = 20;
    translation[ 1 ] = 40;
    translation[ 2 ] = 50;
    transform->Translate( translation );
    transform->Print( std::cout );

    // Tranform the mesh collection
    std::cout << "Transforming mesh collection" << std::endl;
    for ( int meshNumber = -1; meshNumber < static_cast< int >( numberOfMeshes ); meshNumber++ )
    {
      meshCollection->Transform( meshNumber, transform );
    }

    meshCollection->Write( "debug.txt" );

#endif

    // Now initialize the estimator with the mesh
    m_Estimator->SetInitialMeshCollection( meshCollection );

  }

}



//
//
//
void
AtlasParameterEstimationConsole
::Show()
{
  m_Window->show();
  Fl::run();

}


//
//
//
void
AtlasParameterEstimationConsole
::DisplayLabelImage( unsigned int labelImageNumber )
{
  // Show the label image
  //m_Estimator->GetLabelImage( labelImageNumber )->Print( std::cout );
  m_LabelImageViewer->SetImage( m_Estimator->GetLabelImage( labelImageNumber ) );
  m_LabelImageViewer->SetScaleToFillRange();

  if ( !m_Estimator->GetUseGaussians() )
  {
    // Show the alpha image
    AtlasMeshAlphaDrawer::Pointer  alphaDrawer = AtlasMeshAlphaDrawer::New();
    alphaDrawer->SetLabelImage( m_Estimator->GetLabelImage( labelImageNumber ) );
    alphaDrawer->SetLabelNumber( static_cast< unsigned char >( m_LabelNumber->value() ) );
    alphaDrawer->Rasterize( m_Estimator->GetCurrentMeshCollection()->GetMesh( labelImageNumber ) );

    typedef itk::IntensityWindowingImageFilter< AtlasMeshAlphaDrawer::AlphaImageType,
            ImageViewer::ImageType >   WindowerType;
    WindowerType::Pointer  windower = WindowerType::New();
    windower->SetInput( alphaDrawer->GetAlphaImage() );
    windower->SetWindowMinimum( 0 );
    windower->SetWindowMaximum( 1 );
    windower->SetOutputMinimum( 0 );
    windower->SetOutputMaximum( 255 );
    windower->Update();

    m_AlphaImageViewer->SetImage( windower->GetOutput() );
    m_AlphaImageViewer->SetScale( 1.0f );
    m_LabelImageViewer->SetOverlayImage( windower->GetOutput() );
    m_LabelImageViewer->SetOverlayScale( 1.0f );
  }
  else
  {
    // Show the reconstucted image
    AtlasMeshSummaryDrawer::Pointer  summaryDrawer = AtlasMeshSummaryDrawer::New();
    summaryDrawer->SetLabelImage( m_Estimator->GetLabelImage( labelImageNumber ) );
    summaryDrawer->Rasterize( m_Estimator->GetCurrentMeshCollection()->GetMesh( labelImageNumber ) );

    typedef itk::IntensityWindowingImageFilter< AtlasMeshSummaryDrawer::SummaryImageType,
            ImageViewer::ImageType >   WindowerType;
    WindowerType::Pointer  windower = WindowerType::New();
    windower->SetInput( summaryDrawer->GetSummaryImage() );
    windower->SetWindowMinimum( 0 );
    windower->SetWindowMaximum( m_Estimator->GetNumberOfClasses() );
    windower->SetOutputMinimum( 0 );
    windower->SetOutputMaximum( 255 );
    windower->Update();

    m_AlphaImageViewer->SetImage( windower->GetOutput() );
    m_AlphaImageViewer->SetScale( 1.0f );
  }


  // Show the mesh overlaid
  if ( m_ShowMesh->value() )
  {
    AtlasMesh::ConstPointer  mesh = m_Estimator->GetCurrentMeshCollection()->GetMesh( labelImageNumber );
    m_LabelImageViewer->SetMesh( mesh );
    m_AlphaImageViewer->SetMesh( mesh );
  }
  else
  {
    m_LabelImageViewer->SetMesh( 0 );
    m_AlphaImageViewer->SetMesh( 0 );
  }

  // Show the position gradient overlaid
  if ( m_ShowGradient->value() )
  {
    AtlasPositionGradientContainerType::Pointer  positionGradient =
      m_Estimator->GetCurrentPositionGradient( labelImageNumber );
    m_LabelImageViewer->SetPositionGradient( positionGradient );
    m_AlphaImageViewer->SetPositionGradient( positionGradient );
  }
  else
  {
    m_LabelImageViewer->SetPositionGradient( 0 );
    m_AlphaImageViewer->SetPositionGradient( 0 );
  }

  // Redraw
  m_LabelImageViewer->redraw();
  m_AlphaImageViewer->redraw();
  Fl::check();

}



//
//
//
void
AtlasParameterEstimationConsole
::Estimate()
{
  if ( m_UseExplicitStartCollection->value() )
  {
#if 1
    m_Estimator->Estimate();
#elif 1
    const float  originalK = m_Estimator->GetCurrentMeshCollection()->GetK();
    float  initialStiffness = originalK;
    if ( initialStiffness < 0.1f )
    {
      initialStiffness = 0.1f;
    }
    for ( float K = initialStiffness / 2; K > originalK; K /= 2 )
    {
      std::cout << "Dividing K by factor of 2 and estimating for K " << K << std::endl;
      AtlasMeshCollection::Pointer  meshCollection =
        const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );
      meshCollection->SetK( K );
      meshCollection->FlattenAlphas();
      m_Estimator->SetInitialMeshCollection( meshCollection );
      m_Estimator->Estimate();


#if 1
      {
        // Calculate the data cost and the alpha cost
        AtlasMeshCollectionModelLikelihoodCalculator::Pointer  dataAndAlphaCostCalculator =
          AtlasMeshCollectionModelLikelihoodCalculator::New();
        dataAndAlphaCostCalculator->SetMeshCollection( m_Estimator->GetCurrentMeshCollection() );
        dataAndAlphaCostCalculator->SetLabelImages( m_Estimator->GetLabelImages() );
        float  dataCost;
        float  alphasCost;
        dataAndAlphaCostCalculator->GetDataCostAndAlphasCost( dataCost, alphasCost );

        // Calculate the position cost
        AtlasMeshCollectionPositionCostCalculator::Pointer  positionCostCalculator =
          AtlasMeshCollectionPositionCostCalculator::New();
        positionCostCalculator->SetMeshCollection(
          const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() ) );
        positionCostCalculator->SetLabelImages( m_Estimator->GetLabelImages() );
        const float positionCost = positionCostCalculator->GetPositionCost();

        // Output total cost
        std::cout << "Total cost: " << std::endl;
        std::cout << "                 dataCost: " << dataCost << std::endl;
        std::cout << "               alphasCost: " << alphasCost << std::endl;
        std::cout << "             positionCost: " << positionCost << std::endl;
        std::cout << "  + ---------------------  : " << std::endl;
        std::cout << "                 " << dataCost + alphasCost + positionCost << std::endl;
        std::cout << std::endl;
      }
#endif


    }

    //
    std::cout << "Now that we're here, let's try estimating with K " << originalK << std::endl;
    AtlasMeshCollection::Pointer  meshCollection
    = const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );
    meshCollection->SetK( originalK );
    meshCollection->FlattenAlphas();
    m_Estimator->SetInitialMeshCollection( meshCollection );
    m_Estimator->Estimate();

#else
    AtlasMeshCollection::Pointer  meshCollection =
      const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );
    const float  originalK = meshCollection->GetK();
    meshCollection->SetK( 0.1 );
    m_Estimator->SetInitialMeshCollection( meshCollection );

    const int  numberOfSmoothingSteps = 10;
    for ( int smoothingStepNumber = 0; smoothingStepNumber < numberOfSmoothingSteps; smoothingStepNumber++ )
    {
      const float  alphasSmoothingFactor =
        static_cast< float >( numberOfSmoothingSteps - 1 - smoothingStepNumber ) /
        static_cast< float >( numberOfSmoothingSteps - 1 );
      std::cout << "Using K " << 0.1 << "smoothingFactor " << alphasSmoothingFactor << std::endl;
      m_Estimator->SetAlphasSmoothingFactor( alphasSmoothingFactor );
      m_Estimator->Estimate();
    }

    //
    for ( float K = 0.1 / 2; K > originalK; K /= 2 )
    {
      std::cout << "Dividing K by factor of 2 and estimating for K " << K << std::endl;
      meshCollection = const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );
      meshCollection->SetK( K );
      meshCollection->FlattenAlphas();
      m_Estimator->SetInitialMeshCollection( meshCollection );
      m_Estimator->Estimate();
    }

    //
    std::cout << "Now that we're here, let's try estimating with K " << originalK << std::endl;
    meshCollection = const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );
    meshCollection->SetK( originalK );
    meshCollection->FlattenAlphas();
    m_Estimator->SetInitialMeshCollection( meshCollection );
    m_Estimator->Estimate();
#endif
  }
  else
  {
    // Retrieve number of upsampling steps from the GUI
    const int  numberOfUpsamplingSteps = static_cast< const int >( m_NumberOfUpsamplingSteps->value() );
#if 1
    float  initialStiffness = m_InitialStiffness->value();
#else
    float  initialStiffness = m_InitialStiffness->value();
    if ( initialStiffness < 0.1f )
    {
      initialStiffness = 0.1f;
    }

#endif

    for ( int upsamplingStepNumber = 0; upsamplingStepNumber <= numberOfUpsamplingSteps; upsamplingStepNumber++ )
    {
#if 0
      m_Estimator->Estimate();
#else
      AtlasMeshCollection::Pointer  meshCollection =
        const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );
      meshCollection->SetK( initialStiffness );
      m_Estimator->SetInitialMeshCollection( meshCollection );
      m_Estimator->Estimate();
#endif


      // Write out what we have so far
      std::ostringstream  fileNameStream;
      fileNameStream << "estimatedAtUpsamplingStep" << upsamplingStepNumber << "Of"
                     << numberOfUpsamplingSteps << ".txt";
      m_Estimator->GetCurrentMeshCollection()->Write( fileNameStream.str().c_str() );


      // If this is not the final resolution yet, upsample mesh collection
      if ( upsamplingStepNumber != numberOfUpsamplingSteps )
      {
        AtlasMeshCollection::Pointer  upsampledMeshCollection =
          m_Estimator->GetCurrentMeshCollection()->GetUpsampled();
        m_Estimator->SetInitialMeshCollection( upsampledMeshCollection );
      }




    } // End loop over upsampling steps

#if 0

#else
    const float  originalK = m_InitialStiffness->value();
    for ( float K = initialStiffness / 2; K > originalK; K /= 2 )
    {
      std::cout << "Dividing K by factor of 2 and estimating for K " << K << std::endl;
      AtlasMeshCollection::Pointer  meshCollection =
        const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );
      meshCollection->SetK( K );
      meshCollection->FlattenAlphas();
      m_Estimator->SetInitialMeshCollection( meshCollection );
      m_Estimator->Estimate();

#if 1
      if ( true )
      {
        // Calculate the data cost and the alpha cost
        AtlasMeshCollectionModelLikelihoodCalculator::Pointer  dataAndAlphaCostCalculator =
          AtlasMeshCollectionModelLikelihoodCalculator::New();
        dataAndAlphaCostCalculator->SetMeshCollection( m_Estimator->GetCurrentMeshCollection() );
        dataAndAlphaCostCalculator->SetLabelImages( m_Estimator->GetLabelImages() );
        float  dataCost;
        float  alphasCost;
        dataAndAlphaCostCalculator->GetDataCostAndAlphasCost( dataCost, alphasCost );

        // Calculate the position cost
        AtlasMeshCollectionPositionCostCalculator::Pointer  positionCostCalculator =
          AtlasMeshCollectionPositionCostCalculator::New();
        positionCostCalculator->SetMeshCollection(
          const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() ) );
        positionCostCalculator->SetLabelImages( m_Estimator->GetLabelImages() );
        const float positionCost = positionCostCalculator->GetPositionCost();

        // Output total cost
        std::cout << "Total cost: " << std::endl;
        std::cout << "                 dataCost: " << dataCost << std::endl;
        std::cout << "               alphasCost: " << alphasCost << std::endl;
        std::cout << "             positionCost: " << positionCost << std::endl;
        std::cout << "  + ---------------------  : " << std::endl;
        std::cout << "                 " << dataCost + alphasCost + positionCost << std::endl;
        std::cout << std::endl;
      }
#endif


    }

    if (  initialStiffness != originalK )
    {
      //
      std::cout << "Now that we're here, let's try estimating with K " << originalK << std::endl;
      AtlasMeshCollection::Pointer  meshCollection
      = const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );
      meshCollection->SetK( originalK );
      meshCollection->FlattenAlphas();
      m_Estimator->SetInitialMeshCollection( meshCollection );
      m_Estimator->Estimate();
    }

#endif

  } // End test for explicit start collection

  // Write result out
  m_Estimator->GetCurrentMeshCollection()->Write( "estimated.txt" );

#if 0
  // Calculate the data cost and the alpha cost
  AtlasMeshCollectionModelLikelihoodCalculator::Pointer  dataAndAlphaCostCalculator =
    AtlasMeshCollectionModelLikelihoodCalculator::New();
  dataAndAlphaCostCalculator->SetMeshCollection( m_Estimator->GetCurrentMeshCollection() );
  dataAndAlphaCostCalculator->SetLabelImages( m_Estimator->GetLabelImages() );
  float  dataCost;
  float  alphasCost;
  dataAndAlphaCostCalculator->GetDataCostAndAlphasCost( dataCost, alphasCost );

  // Calculate the position cost
  AtlasMeshCollectionPositionCostCalculator::Pointer  positionCostCalculator =
    AtlasMeshCollectionPositionCostCalculator::New();
  positionCostCalculator->SetMeshCollection(
    const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() ) );
  positionCostCalculator->SetLabelImages( m_Estimator->GetLabelImages() );
  const float positionCost = positionCostCalculator->GetPositionCost();

  // Output total cost
  std::cout << "Total cost: " << std::endl;
  std::cout << "                 dataCost: " << dataCost << std::endl;
  std::cout << "               alphasCost: " << alphasCost << std::endl;
  std::cout << "             positionCost: " << positionCost << std::endl;
  std::cout << "  + ---------------------  : " << std::endl;
  std::cout << "                 " << dataCost + alphasCost + positionCost << std::endl;
  std::cout << std::endl;
#endif

}




//
//
//
void
AtlasParameterEstimationConsole
::HandleEstimatorEvent( itk::Object* object, const itk::EventObject & event )
{

  if ( typeid( event ) == typeid( itk::EndEvent ) )
  {
    m_TotalProgressLabel = "Inactive";
    m_TotalProgress->label( m_TotalProgressLabel.c_str() );
    m_TotalProgress->value( 0 );
    m_TotalProgress->redraw();
  }
  else if ( ( typeid( event ) == typeid( itk::StartEvent ) ) ||
            ( typeid( event ) == typeid( itk::IterationEvent ) ) )
  {
    // Show total progress
    float  totalProgress = static_cast< float >( m_Estimator->GetIterationNumber() + 1 ) /
                           static_cast< float >( m_Estimator->GetMaximumNumberOfIterations() );
    std::ostringstream  labelStream;
    labelStream << "Iteration " << m_Estimator->GetIterationNumber() + 1 << " / "
                << m_Estimator->GetMaximumNumberOfIterations()
                << "  (" << m_Estimator->GetCurrentMinLogLikelihoodTimesPrior() << ")";
    m_TotalProgressLabel = labelStream.str();
    m_TotalProgress->label( m_TotalProgressLabel.c_str() );
    m_TotalProgress->value( totalProgress * 100 );
    m_TotalProgress->redraw();

    // Show current result by pretending to alter the labelImageNumber GUI
    m_LabelImageNumber->do_callback();
  }
  else if ( typeid( event ) == typeid( AlphasEstimationEndEvent ) )
  {
    m_SubProgressLabel = "";
    m_SubProgress->label( m_SubProgressLabel.c_str() );
    m_SubProgress->value( 0 );
    m_SubProgress->redraw();
  }
  else if ( ( typeid( event ) == typeid( AlphasEstimationStartEvent ) ) ||
            ( typeid( event ) == typeid( AlphasEstimationIterationEvent ) ) )
  {
    // Show sub progress
    m_SubProgressLabel = "Estimating alphas...";
    m_SubProgress->label( m_SubProgressLabel.c_str() );
    float  subProgress =
      static_cast< float >( m_Estimator->GetAlphasEstimationIterationNumber() + 1 ) /
      static_cast< float >( m_Estimator->GetAlphasEstimationMaximumNumberOfIterations() );
    m_SubProgress->value( subProgress * 100 );
    m_SubProgress->redraw();

    // Only update the viewers if we are asked for it
    if ( m_ShowAlphasEstimationIterations->value() )
    {
      m_LabelImageNumber->do_callback();
    }
  }
  else if ( typeid( event ) == typeid( PositionEstimationEndEvent ) )
  {
    m_SubProgressLabel = "";
    m_SubProgress->label( m_SubProgressLabel.c_str() );
    m_SubProgress->value( 0 );
    m_SubProgress->redraw();
  }
  else if ( ( typeid( event ) == typeid( PositionEstimationStartEvent ) ) ||
            ( typeid( event ) == typeid( PositionEstimationIterationEvent ) ) )
  {
    // Show sub progress
    std::ostringstream  labelStream;
    labelStream << "Estimating position of label image " << m_Estimator->GetLabelImageNumber() + 1
                << " / " << m_Estimator->GetNumberOfLabelImages();
    m_SubProgressLabel = labelStream.str();
    m_SubProgress->label( m_SubProgressLabel.c_str() );
    float  subProgress =
      ( static_cast< float >( m_Estimator->GetLabelImageNumber() ) +
        static_cast< float >( m_Estimator->GetPositionEstimationIterationNumber() + 1 ) /
        static_cast< float >( m_Estimator->GetPositionEstimationMaximumNumberOfIterations() ) ) /
      static_cast< float >( m_Estimator->GetNumberOfLabelImages() );
    m_SubProgress->value( subProgress * 100 );
    m_SubProgress->redraw();

    // Only update the viewers if we are asked for it
    //if ( m_ShowPositionEstimationIterations->value() )
    //  {
    m_LabelImageNumber->value( m_Estimator->GetLabelImageNumber() );
    m_LabelImageNumber->do_callback();
    //  }
  }

  Fl::check();


  while ( m_Interrupted )
  {
    Fl::wait();
  }

  if ( m_Stepping )
  {
    this->Interrupt();
    m_Stepping = false;
  }


}





#if 0
//
//
//
static bool IsLeftRotatingTriangle( AtlasMesh::PointType p0,
                                    AtlasMesh::PointType p1,
                                    AtlasMesh::PointType p2 )
{

  AtlasMesh::PointType::VectorType  vec01 = p1 - p0;
  AtlasMesh::PointType::VectorType  vecOrth02;
  vecOrth02[0] = p2[1] - p0[1];
  vecOrth02[1] = -p2[0] + p0[0];
  if ( vec01 * vecOrth02 <= 0  )
  {
    // Invalid triangle.
    return false;
  }

  return true;

}


//
//
//
void
AtlasParameterEstimationConsole
::SelectTriangleContainingPoint( float x, float y )
{

  unsigned long  triangleId;
  if ( m_LabelImageViewer->GetTriangleContainingPoint( x, y, triangleId ) )
  {
    // Display the selected triangle in the viewers
    m_LabelImageViewer->SetTriangleIdToHighlight( triangleId );
    m_AlphaImageViewer->SetTriangleIdToHighlight( triangleId );

    // Retrieve the point id of the nearest vertex clicked
    unsigned long  nearestVertexPointId = m_LabelImageViewer->GetNearestVertexPointId( x, y );

    // Display the approximative Gaussian in the viewers
    m_LabelImageViewer->SetPointIdToShowGaussianApproximation( nearestVertexPointId );
    m_AlphaImageViewer->SetPointIdToShowGaussianApproximation( nearestVertexPointId );

    // Redraw
    m_LabelImageViewer->redraw();
    m_AlphaImageViewer->redraw();
    Fl::check();
  }

}

#endif



//
//
//
void
AtlasParameterEstimationConsole
::Interrupt()
{
  m_Interrupted = true;

}



//
//
//
void
AtlasParameterEstimationConsole
::Step()
{
  m_Stepping = true;
  this->Continue();

}



//
//
//
void
AtlasParameterEstimationConsole
::Continue()
{
  m_Interrupted = false;

}


//
//
//
void
AtlasParameterEstimationConsole
::SetPositionEstimationResolution( unsigned int positionEstimationResolution )
{
  m_Estimator->SetPositionEstimationIterationEventResolution( positionEstimationResolution );
}







} // end namespace kvl

