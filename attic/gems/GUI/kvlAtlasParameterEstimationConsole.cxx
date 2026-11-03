#include "kvlAtlasParameterEstimationConsole.h"

#include "FL/Fl.H"
#include "FL/fl_ask.H"
#include "itkImageFileReader.h"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkCommand.h"
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
   
  switch( m_Estimator->GetPositionOptimizer() ) 
    {
    case AtlasParameterEstimator::FIXED_STEP_GRADIENT_DESCENT: 
      {
      m_FixedStepGradientDescent->setonly();
      break;
      } 
    case AtlasParameterEstimator::GRADIENT_DESCENT: 
      {
      m_GradientDescent->setonly();
      break;
      } 
    case AtlasParameterEstimator::CONJUGATE_GRADIENT: 
      {
      m_ConjugateGradient->setonly();
      break;
      } 
    default:
      {
      m_LBFGS->setonly();
      break;
      }
    }

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
::SetLabelImages( const std::vector< std::string >&  fileNames )
{
  // Sanity checking
  if ( fileNames.size() == 0 )
    return;
    
  // Read the images
  typedef CompressionLookupTable::ImageType  InputImageType;
  std::vector< InputImageType::ConstPointer >  labelImages;
  for ( std::vector< std::string >::const_iterator it = fileNames.begin();
         it != fileNames.end(); ++it )
    {
    // Read the input image
    typedef itk::ImageFileReader< InputImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( ( *it ).c_str() );
    reader->Update();
    InputImageType::ConstPointer  labelImage = reader->GetOutput();

    // Over-ride the spacing and origin since at this point we can't deal with that
    const double spacing[] = { 1, 1, 1 };
    const double origin[] = { 0, 0, 0 };
    const_cast< InputImageType* >( labelImage.GetPointer() )->SetSpacing( spacing );
    const_cast< InputImageType* >( labelImage.GetPointer() )->SetOrigin( origin );

    // Remember this image
    labelImages.push_back( labelImage );
    }


  // Build a lookup table that maps the original intensities into classes
  CompressionLookupTable::Pointer  lookupTable = CompressionLookupTable::New();
  lookupTable->Construct( labelImages );
  lookupTable->Write( "compressionLookupTable.txt" );

  // Give the label images to the estimator
  m_Estimator->SetLabelImages( labelImages, lookupTable );
  
  // Construct the initial mesh
  this->InitializeMesh();

  // Change the GUI according to the new images
  for ( unsigned int i = 0; i < labelImages.size(); i++ )
    {
    std::ostringstream  labelStream;
    labelStream << "Image " << i;
    m_LabelImageNumber->add( labelStream.str().c_str() );
    }

  for ( int classNumber = 0; classNumber < lookupTable->GetNumberOfClasses(); classNumber++ )
    {
    m_LabelNumber->add( ( lookupTable->GetLabelName( classNumber ) ).c_str() );
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
    return;

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
    std::cout << "meshSize: " << meshSize[0] << " " << meshSize[1] << " " << meshSize[2] << std::endl;
    std::cout << "domainSize: " << domainSize[0] << " " << domainSize[1] << " " << domainSize[2] << std::endl;
    std::cout << "initialStiffness: " << initialStiffness << std::endl;
    std::cout << "numberOfClasses: " << numberOfClasses << std::endl;
    std::cout << "numberOfMeshes: " << numberOfMeshes << std::endl;
  

    meshCollection->Construct( meshSize, domainSize, initialStiffness, 
                              numberOfClasses, numberOfMeshes );

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
  typedef ImageViewer::ImageType  CompressedImageType;
  CompressedImageType::Pointer  compressedImage = CompressedImageType::New();
  compressedImage->SetRegions( m_Estimator->GetLabelImage( labelImageNumber )->GetLargestPossibleRegion() );
  compressedImage->Allocate();
  compressedImage->FillBuffer( 0 );
  typedef CompressionLookupTable::ImageType  ImageType;
  itk::ImageRegionConstIterator< ImageType >  origIt( m_Estimator->GetLabelImage( labelImageNumber ), 
                                                      m_Estimator->GetLabelImage( labelImageNumber )->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< CompressedImageType >  compIt( compressedImage, 
                                                           compressedImage->GetLargestPossibleRegion() );
  for ( ; !origIt.IsAtEnd(); ++origIt, ++compIt )
    {
    compIt.Set( static_cast< CompressedImageType::PixelType >( m_Estimator->GetCompressionLookupTable()->GetClassNumbers( origIt.Get() )[0] ) );  
    }
  m_LabelImageViewer->SetImage( compressedImage );
  m_LabelImageViewer->SetScaleToFillRange();

  // Show the alpha image
  AtlasMeshAlphaDrawer::Pointer  alphaDrawer = AtlasMeshAlphaDrawer::New();
  alphaDrawer->SetRegions( compressedImage->GetLargestPossibleRegion() );
  alphaDrawer->SetClassNumber( m_LabelNumber->value() );
  alphaDrawer->Rasterize( m_Estimator->GetCurrentMeshCollection()->GetMesh( labelImageNumber ) );

  typedef itk::IntensityWindowingImageFilter< AtlasMeshAlphaDrawer::ImageType, 
                                              ImageViewer::ImageType >   WindowerType;
  WindowerType::Pointer  windower = WindowerType::New();
  windower->SetInput( alphaDrawer->GetImage() );  
  windower->SetWindowMinimum( 0 );
  windower->SetWindowMaximum( 1 );
  windower->SetOutputMinimum( 0 );
  windower->SetOutputMaximum( 255 );
  windower->Update();

  m_AlphaImageViewer->SetImage( windower->GetOutput() );
  m_AlphaImageViewer->SetScale( 1.0f );
  m_LabelImageViewer->SetOverlayImage( windower->GetOutput() );
  m_LabelImageViewer->SetOverlayScale( 1.0f );

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
  // Inform the estimator about the user-selected deformation optimizer
  if ( m_FixedStepGradientDescent->value() )
    {
    m_Estimator->SetPositionOptimizer( AtlasParameterEstimator::FIXED_STEP_GRADIENT_DESCENT );  
    }
  else if ( m_GradientDescent->value() )
    {
    m_Estimator->SetPositionOptimizer( AtlasParameterEstimator::GRADIENT_DESCENT );        
    }
  else if ( m_ConjugateGradient->value() )
    {
    m_Estimator->SetPositionOptimizer( AtlasParameterEstimator::CONJUGATE_GRADIENT );        
    }
  else
    {
    m_Estimator->SetPositionOptimizer( AtlasParameterEstimator::LBFGS );        
    }
  
  //
  if ( m_UseExplicitStartCollection->value() )
    {
    m_Estimator->Estimate();
    }
  else
    {
    // Retrieve number of upsampling steps from the GUI
    const int  numberOfUpsamplingSteps = static_cast< const int >( m_NumberOfUpsamplingSteps->value() );
    float  initialStiffness = m_InitialStiffness->value();
    
    for ( int upsamplingStepNumber = 0; upsamplingStepNumber <= numberOfUpsamplingSteps; upsamplingStepNumber++ )
      {
      AtlasMeshCollection::Pointer  meshCollection = 
            const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );
      meshCollection->SetK( initialStiffness );
      m_Estimator->SetInitialMeshCollection( meshCollection );
      m_Estimator->Estimate();
      
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
    } // End test for explicit start collection

  // Write result out
  m_Estimator->GetCurrentMeshCollection()->Write( "estimated.txt" );
  
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
    Fl::wait();

  if ( m_Stepping )
    {
    this->Interrupt();
    m_Stepping = false;
    }


}
 
  




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

