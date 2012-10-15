/**
 * @file  kvlAtlasMeshSegmenterConsole.cxx
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
#include "kvlAtlasMeshSegmenterConsole.h"

#include "FL/Fl.H"
#include "FL/fl_ask.H"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkCommand.h"


namespace kvl
{

//
//
//
AtlasMeshSegmenterConsole
::AtlasMeshSegmenterConsole()
{

  m_Driver = AtlasMeshSegmentationDriver::New();

  m_Interrupted = false;
  m_Stepping = false;

  // m_GradientDescentStepSize->value( m_Driver->GetAtlasMeshSegmenter()->GetPositionGradientDescentStepSize() );
  // m_PositionUpdatingStopCriterion->value( m_Driver->GetAtlasMeshSegmenter()->GetPositionUpdatingStopCriterion() );

  // Add observers to the driver's segmenter's EMSegmenter
  typedef itk::MemberCommand< AtlasMeshSegmenterConsole >   MemberCommandType;
  MemberCommandType::Pointer  EMSegmenterCommand = MemberCommandType::New();
  EMSegmenterCommand->SetCallbackFunction( this, &AtlasMeshSegmenterConsole::HandleEMSegmenterEvent );
  m_Driver->AddEMSegmenterObserver( itk::StartEvent(), EMSegmenterCommand );
  m_Driver->AddEMSegmenterObserver( itk::IterationEvent(), EMSegmenterCommand );
  m_Driver->AddEMSegmenterObserver( itk::EndEvent(), EMSegmenterCommand );

  // Add observers to the driver's segmenter
  MemberCommandType::Pointer  segmenterCommand = MemberCommandType::New();
  segmenterCommand->SetCallbackFunction( this, &AtlasMeshSegmenterConsole::HandleSegmenterEvent );
  m_Driver->AddSegmenterObserver( itk::StartEvent(), segmenterCommand );
  m_Driver->AddSegmenterObserver( itk::IterationEvent(), segmenterCommand );
  m_Driver->AddSegmenterObserver( itk::EndEvent(), segmenterCommand );
  m_Driver->AddSegmenterObserver( PositionUpdatingStartEvent(), segmenterCommand );
  m_Driver->AddSegmenterObserver( PositionUpdatingIterationEvent(), segmenterCommand );
  m_Driver->AddSegmenterObserver( PositionUpdatingEndEvent(), segmenterCommand );

  // Add observers to the driver's segmenter
  MemberCommandType::Pointer  driverCommand = MemberCommandType::New();
  driverCommand->SetCallbackFunction( this, &AtlasMeshSegmenterConsole::HandleDriverEvent );
  m_Driver->AddObserver( MultiResolutionStartEvent(), driverCommand );
  m_Driver->AddObserver( MultiResolutionIterationEvent(), driverCommand );
  m_Driver->AddObserver( MultiResolutionEndEvent(), driverCommand );


}



//
//
//
AtlasMeshSegmenterConsole
::~AtlasMeshSegmenterConsole()
{

}



//
//
//
void
AtlasMeshSegmenterConsole
::SetUp( const std::string& setUpFileName )
{

  // Set up the driver
  m_Driver->SetUp( setUpFileName );

  // Change the GUI according to the input
  typedef AtlasMeshSegmenter::ImageType  ImageType;
  typedef itk::MinimumMaximumImageCalculator< ImageType >  RangeCalculatorType;
  RangeCalculatorType::Pointer  rangeCalculator = RangeCalculatorType::New();
  rangeCalculator->SetImage( m_Driver->GetAtlasMeshSegmenter()->GetImage() );
  rangeCalculator->Compute();
  typedef itk::IntensityWindowingImageFilter< ImageType, ImageViewer::ImageType >   WindowerType;
  WindowerType::Pointer  windower = WindowerType::New();
  windower->SetInput( m_Driver->GetAtlasMeshSegmenter()->GetImage() );
  windower->SetWindowMinimum( rangeCalculator->GetMinimum() );
  windower->SetWindowMaximum( rangeCalculator->GetMaximum() );
  windower->SetOutputMinimum( 0 );
  windower->SetOutputMaximum( 255 );
  windower->Update();
  m_ImageViewer->SetImage( windower->GetOutput() );
  m_ImageViewer->SetScale( 1.0f );


  // Fill in the GUI for the label names
  for ( std::vector< std::string >::const_iterator it = m_Driver->GetLabels().begin();
        it != m_Driver->GetLabels().end(); ++it )
  {
    m_LabelSelecter->add( it->c_str() );
  }


  // Set the GUI elements corresponding to the slice locations
  m_SagittalSliceNumber->maximum( m_ImageViewer->GetMaximumImageIndex()[ 0 ] );
  m_CoronalSliceNumber->maximum( m_ImageViewer->GetMaximumImageIndex()[ 1 ] );
  m_AxialSliceNumber->maximum( m_ImageViewer->GetMaximumImageIndex()[ 2 ] );

  unsigned int  sagittalSliceNumber;
  unsigned int  coronalSliceNumber;
  unsigned int  axialSliceNumber;
  m_ImageViewer->GetSliceLocation( sagittalSliceNumber, coronalSliceNumber, axialSliceNumber );
  m_SagittalSliceNumber->value( sagittalSliceNumber );
  m_CoronalSliceNumber->value( coronalSliceNumber );
  m_AxialSliceNumber->value( axialSliceNumber );


  m_LabelSelecter->value( 1 );
  m_LabelSelecter->do_callback();

}


//
//
//
void
AtlasMeshSegmenterConsole
::Show()
{
  m_Window->show();
  Fl::run();

}



//
//
//
void
AtlasMeshSegmenterConsole
::DisplayLabel( unsigned int labelNumber )
{
  std::cout << "DisplayLabel - begin" << std::endl;

  // Show the prior
  typedef itk::IntensityWindowingImageFilter< EMSegmenter::ClassificationImageType,
          ImageViewer::ImageType >   WindowerType;
  WindowerType::Pointer  priorWindower = WindowerType::New();
  priorWindower->SetInput( m_Driver->GetAtlasMeshSegmenter()->GetEMSegmenter()->GetPrior( labelNumber ) );
  priorWindower->SetWindowMinimum( 0 );
  priorWindower->SetWindowMaximum( 1 );
  priorWindower->SetOutputMinimum( 0 );
  priorWindower->SetOutputMaximum( 255 );
  priorWindower->Update();
  m_PriorViewer->SetImage( priorWindower->GetOutput() );
  m_PriorViewer->SetScale( 1.0f );

  // Show the posterior
  WindowerType::Pointer  posteriorWindower = WindowerType::New();
  posteriorWindower->SetInput( m_Driver->GetAtlasMeshSegmenter()->GetEMSegmenter()->GetPosterior( labelNumber ) );
  posteriorWindower->SetWindowMinimum( 0 );
  posteriorWindower->SetWindowMaximum( 1 );
  posteriorWindower->SetOutputMinimum( 0 );
  posteriorWindower->SetOutputMaximum( 255 );
  posteriorWindower->Update();
  m_PosteriorViewer->SetImage( posteriorWindower->GetOutput() );
  m_PosteriorViewer->SetScale( 1.0f );

  // Overlay the posterior on the gray scale image
  if ( m_ShowOverlay->value() )
  {
    if ( m_OverlayPrior->value() )
    {
      m_ImageViewer->SetOverlayImage( priorWindower->GetOutput() );
    }
    else
    {
      m_ImageViewer->SetOverlayImage( posteriorWindower->GetOutput() );
    }
    m_ImageViewer->SetOverlayScale( 1.0f );
    m_ImageViewer->SetOverlayAlpha( 0.7 );
  }
  else
  {
    //m_ImageViewer->SetOverlayImage( 0 );
    m_ImageViewer->SetOverlayAlpha( 0 );
  }

  // Show the reconstruction
  if ( m_Driver->GetAtlasMeshSegmenter()->GetMeans().size() != 0 )
  {
    // Get the reconstruction
    typedef AtlasMeshSegmentationDriver::SummaryImageType  SummaryImageType;
    SummaryImageType::Pointer  reconstructedImage =
      m_Driver->GetSummaryImage( m_OverlayPrior->value(), false, true, m_Driver->GetAtlasMeshSegmenter()->GetMeans() );

    // Display the reconstruction
    typedef itk::MinimumMaximumImageCalculator< SummaryImageType >  RangeCalculatorType;
    RangeCalculatorType::Pointer  rangeCalculator = RangeCalculatorType::New();
    rangeCalculator->SetImage( reconstructedImage );
    rangeCalculator->Compute();
    typedef itk::IntensityWindowingImageFilter<  SummaryImageType,
            ImageViewer::ImageType >   WindowerType;
    WindowerType::Pointer  windower = WindowerType::New();
    windower->SetInput( reconstructedImage );
    windower->SetWindowMinimum( rangeCalculator->GetMinimum() );
    windower->SetWindowMaximum( rangeCalculator->GetMaximum() );
    windower->SetOutputMinimum( 0 );
    windower->SetOutputMaximum( 255 );
    windower->Update();
    m_ReconstructionViewer->SetImage( windower->GetOutput() );
    m_ReconstructionViewer->SetScale( 1.0f );

    // Overlay the posterior on the reconstruction
    if ( m_ShowOverlay->value() )
    {
      if ( m_OverlayPrior->value() )
      {
        m_ReconstructionViewer->SetOverlayImage( priorWindower->GetOutput() );
      }
      else
      {
        m_ReconstructionViewer->SetOverlayImage( posteriorWindower->GetOutput() );
      }
      m_ReconstructionViewer->SetOverlayScale( 1.0f );
      m_ReconstructionViewer->SetOverlayAlpha( 0.7 );
    }
    else
    {
      //m_ReconstructionViewer->SetOverlayImage( 0 );
      m_ReconstructionViewer->SetOverlayAlpha( 0 );
    }

  }  // End test if means.size() != 0



  // Show the mesh overlaid
  if ( m_ShowMesh->value() )
  {
    m_ImageViewer->SetMesh( m_Driver->GetAtlasMeshSegmenter()->GetMesh() );
    m_PriorViewer->SetMesh( m_Driver->GetAtlasMeshSegmenter()->GetMesh() );
    m_PosteriorViewer->SetMesh( m_Driver->GetAtlasMeshSegmenter()->GetMesh() );
    m_ReconstructionViewer->SetMesh( m_Driver->GetAtlasMeshSegmenter()->GetMesh() );
  }
  else
  {
    m_ImageViewer->SetMesh( 0 );
    m_PriorViewer->SetMesh( 0 );
    m_PosteriorViewer->SetMesh( 0 );
    m_ReconstructionViewer->SetMesh( 0 );
  }

  // Show the position gradient overlaid
  if ( m_ShowGradient->value() )
  {
  }
  else
  {
    m_ImageViewer->SetPositionGradient( 0 );
    m_PriorViewer->SetPositionGradient( 0 );
    m_PosteriorViewer->SetPositionGradient( 0 );
    m_ReconstructionViewer->SetPositionGradient( 0 );
  }

  // Redraw
  // m_ImageViewer->redraw();
  // m_PriorViewer->redraw();
  // m_PosteriorViewer->redraw();
  // m_ReconstructionViewer->redraw();
  // Fl::check();
  m_AxialSliceNumber->do_callback();

  std::cout << "DisplayLabel - end" << std::endl;
}



//
//
//
void
AtlasMeshSegmenterConsole
::Segment( bool useAffine )
{

  // Let the beast go
  m_Driver->Segment( useAffine );

}



//
//
//
void
AtlasMeshSegmenterConsole
::HandleDriverEvent( itk::Object* object, const itk::EventObject & event )
{

  if ( typeid( event ) == typeid( MultiResolutionIterationEvent ) )
  {
    // Driver has put up a new image and/or a new mesh; adjust GUI
    typedef AtlasMeshSegmenter::ImageType  ImageType;
    typedef itk::MinimumMaximumImageCalculator< ImageType >  RangeCalculatorType;
    RangeCalculatorType::Pointer  rangeCalculator = RangeCalculatorType::New();
    rangeCalculator->SetImage( m_Driver->GetAtlasMeshSegmenter()->GetImage() );
    rangeCalculator->Compute();
    typedef itk::IntensityWindowingImageFilter< ImageType, ImageViewer::ImageType >   WindowerType;
    WindowerType::Pointer  windower = WindowerType::New();
    windower->SetInput( m_Driver->GetAtlasMeshSegmenter()->GetImage() );
    windower->SetWindowMinimum( rangeCalculator->GetMinimum() );
    windower->SetWindowMaximum( rangeCalculator->GetMaximum() );
    windower->SetOutputMinimum( 0 );
    windower->SetOutputMaximum( 255 );
    windower->Update();
    m_ImageViewer->SetImage( windower->GetOutput() );
    m_ImageViewer->SetScale( 1.0f );

  }



}



//
//
//
void
AtlasMeshSegmenterConsole
::HandleSegmenterEvent( itk::Object* object, const itk::EventObject & event )
{


  if ( typeid( event ) == typeid( itk::EndEvent ) )
  {
    m_ProgressLabel = "Inactive";
    m_Progress->label( m_ProgressLabel.c_str() );
    m_Progress->value( 0 );
    m_Progress->redraw();

    // Show current result by pretending to alter the labelImageNumber GUI
    m_LabelSelecter->do_callback();
  }
  else if ( typeid( event ) == typeid( itk::StartEvent ) )
  {
    // Show total progress
    float  totalProgress = 0.0f;
    std::ostringstream  labelStream;
    labelStream << "Busy";
    m_ProgressLabel = labelStream.str();
    m_Progress->label( m_ProgressLabel.c_str() );
    m_Progress->value( totalProgress * 100 );
    m_Progress->redraw();

    // Show current result by pretending to alter the labelImageNumber GUI
    m_LabelSelecter->do_callback();
  }
  else if ( typeid( event ) == typeid( itk::IterationEvent ) )
  {
    // Show total progress
    // float  totalProgress = 0.0f;
    // std::ostringstream  labelStream;
    // labelStream << "Busy";
    // m_ProgressLabel = labelStream.str();
    // m_Progress->label( m_ProgressLabel.c_str() );
    // m_Progress->value( totalProgress * 100 );
    // m_Progress->redraw();

    // Show current result by pretending to alter the labelImageNumber GUI
    m_LabelSelecter->do_callback();

#if 0

#if 0
    // Write the result out
    for ( unsigned int classNumber = 0; classNumber < m_Driver->GetAtlasMeshSegmenter()->GetNumberOfClasses(); classNumber++ )
    {
      std::ostringstream  fileNameStream;
      fileNameStream << "iteration" << m_Driver->GetAtlasMeshSegmenter()->GetIterationNumber() << "_prior" << classNumber << ".mhd";
      typedef itk::ImageFileWriter< EMSegmenter::ClassificationImageType >  WriterType;
      WriterType::Pointer  writer = WriterType::New();
      writer->SetInput( m_Driver->GetAtlasMeshSegmenter()->GetEMSegmenter()->GetPrior( classNumber ) );
      writer->SetFileName( fileNameStream.str().c_str() );
      writer->Write();
    }

    for ( unsigned int classNumber = 0; classNumber < m_Driver->GetAtlasMeshSegmenter()->GetNumberOfClasses(); classNumber++ )
    {
      std::ostringstream  fileNameStream;
      fileNameStream << "iteration" << m_Driver->GetAtlasMeshSegmenter()->GetIterationNumber() << "_posterior" << classNumber << ".mhd";
      typedef itk::ImageFileWriter< EMSegmenter::ClassificationImageType >  WriterType;
      WriterType::Pointer  writer = WriterType::New();
      writer->SetInput( m_Driver->GetAtlasMeshSegmenter()->GetEMSegmenter()->GetPosterior( classNumber ) );
      writer->SetFileName( fileNameStream.str().c_str() );
      writer->Write();
    }
#endif


    // Write out the mesh
    AtlasMeshCollection::Pointer  collection = AtlasMeshCollection::New();
    collection->GenerateFromSingleMesh( const_cast< AtlasMesh* >( m_Driver->GetAtlasMeshSegmenter()->GetMesh() ), 1, 1000.0f );
#if 0
    collection->SetReferencePosition( m_OriginalMesh->GetPoints() );
#endif
    std::ostringstream  fileNameStream;
    fileNameStream << "iteration" << m_Driver->GetAtlasMeshSegmenter()->GetIterationNumber() << "_mesh.txt";
    collection->Write( fileNameStream.str().c_str() );



    // Write out a prior summary as well
    fileNameStream.str( "" );
    fileNameStream << "iteration" << m_Driver->GetAtlasMeshSegmenter()->GetIterationNumber() << "_priorSummary.mhd";
    typedef itk::ImageFileWriter< SummaryImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( this->GetSummaryImage( true ) );
    writer->SetFileName( fileNameStream.str().c_str() );
    writer->Write();
    fileNameStream.str( "" );
    fileNameStream << "iteration" << m_Driver->GetAtlasMeshSegmenter()->GetIterationNumber() << "_priorSummaryCrisp.mhd";
    writer->SetInput( this->GetSummaryImage( true, true ) );
    writer->SetFileName( fileNameStream.str().c_str() );
    writer->Write();

    // Write out a posterior summary as well
    fileNameStream.str( "" );
    fileNameStream << "iteration" << m_Driver->GetAtlasMeshSegmenter()->GetIterationNumber() << "_posteriorSummary.mhd";
    //writer->SetInput( this->GetSummaryImage() );
    writer->SetInput( this->GetSummaryImage( false, false, true, m_Driver->GetAtlasMeshSegmenter()->GetMeans() ) );
    writer->SetFileName( fileNameStream.str().c_str() );
    writer->Write();
    fileNameStream.str( "" );
    fileNameStream << "iteration" << m_Driver->GetAtlasMeshSegmenter()->GetIterationNumber() << "_posteriorSummaryCrisp.mhd";
    writer->SetInput( this->GetSummaryImage( false, true ) );
    writer->SetFileName( fileNameStream.str().c_str() );
    writer->Write();

#if 0
    {
      typedef itk::Image< unsigned char, 3 >  DebugImageType;
      typedef itk::CastImageFilter< SummaryImageType, DebugImageType >  CasterType;
      CasterType::Pointer  caster = CasterType::New();
      caster->SetInput( this->GetSummaryImage( false, true ) );
      caster->Update();

      typedef itk::ImageFileWriter< DebugImageType >  DebugWriterType;
      DebugWriterType::Pointer  writer = DebugWriterType::New();
      writer->SetInput( caster->GetOutput() );
      writer->SetFileName( "DEBUGGGG.mhd" );
      writer->Update();
    }
#endif


#endif
    // The following only works if none of the classes has more than one Gaussian
    if ( m_Driver->GetAtlasMeshSegmenter()->GetMeans().size() ==
         static_cast< std::size_t >( m_LabelSelecter->size() ) )
    {
      for ( int i = 1; i <= m_LabelSelecter->size(); i++ )
      {
        std::cout << "Mean (+/- st.dev) for " << m_LabelSelecter->text( i ) << ": "
                  << m_Driver->GetAtlasMeshSegmenter()->GetMeans()[ i-1 ]
                  << " (" << sqrt( m_Driver->GetAtlasMeshSegmenter()->GetVariances()[ i-1 ] )
                  << ")" << std::endl;
      }
    }
    else
    {
      // std::cout << "Apparently someone has more than one Gaussian; not showing means estimates" << std::endl;
      // std::cout << "  m_Driver->GetAtlasMeshSegmenter()->GetMeans().size() : "
      //           << m_Driver->GetAtlasMeshSegmenter()->GetMeans().size() << std::endl;
      // std::cout << "  m_LabelSelecter->size() :" << m_LabelSelecter->size() << std::endl;
      for ( unsigned int i = 0; i < m_Driver->GetAtlasMeshSegmenter()->GetMeans().size(); i++ )
      {
        std::cout << "Mean (+/- st.dev) for " << i << ": "
                  << m_Driver->GetAtlasMeshSegmenter()->GetMeans()[ i ]
                  << " (" << sqrt( m_Driver->GetAtlasMeshSegmenter()->GetVariances()[ i ] )
                  << ")" << std::endl;
      }
    }


  }
#if 0
  else if ( typeid( event ) == typeid( kvl::PositionUpdatingStartEvent ) )
  {
    if ( m_Driver->GetAtlasMeshSegmenter()->GetIterationNumber() == 0 )
    {
      std::cout << "This is the start of the registration of the first iteration. " << std::endl;
      std::cout << "For testing purposes, we're gonna fill smth else in as the alphas of the mesh" << std::endl;
      AtlasMeshCollection::Pointer  collection = AtlasMeshCollection::New();
      collection->Read( "testje.txt.gz" );
      ( const_cast< AtlasMesh* >( m_Driver->GetAtlasMeshSegmenter()->GetMesh() ) )->SetPointData( collection->GetPointParameters() );
    }
  }
#endif
  else
  {
    // We have some PositionUpdatingIterationEvent; update the progress bar
    std::ostringstream  progressLabelStream;
    progressLabelStream << "updating position (iteration "
                        << m_Driver->GetAtlasMeshSegmenter()->GetPositionUpdatingIterationNumber() << ")";
    m_ProgressLabel = progressLabelStream.str();
    m_Progress->label( m_ProgressLabel.c_str() );
    const float  progress = static_cast< float >( m_Driver->GetAtlasMeshSegmenter()->GetPositionUpdatingIterationNumber() ) /
                            static_cast< float >( m_Driver->GetAtlasMeshSegmenter()->GetPositionUpdatingMaximumNumberOfIterations() );
    m_Progress->value( progress * 100 );
    m_Progress->redraw();

    // Show current result by pretending to alter the labelImageNumber GUI
    std::cout << "PositionUpdatingIterationEvent is causing the pics to be re-drawn." << std::endl;
    std::cout << "     note that the updated mesh position is *NOT* reflected in the prior, \n"
              << "     since those are only re-calculated just before the GMM is updated in \n"
              << "     the EMSegmenter. So if you want to watch the deformation progress, you'll \n"
              << "     have to switch on the mesh visualization" << std::endl;
    std::cout << "     Re-drawing starting..." << std::flush;
    m_LabelSelecter->do_callback();
    std::cout << "done!" << std::endl;

#if 1
    {
      std::ostringstream  fileNameStream;
      fileNameStream << "debugMesh_iteration" << m_Driver->GetAtlasMeshSegmenter()->GetIterationNumber()
                     << "_positionUpdatingIterationNumber" << m_Driver->GetAtlasMeshSegmenter()->GetPositionUpdatingIterationNumber()
                     << ".txt";

      AtlasMeshCollection::Pointer  collection = AtlasMeshCollection::New();
      collection->GenerateFromSingleMesh( const_cast< AtlasMesh* >( m_Driver->GetAtlasMeshSegmenter()->GetMesh() ), 1, 1000.0f );
      collection->Write( fileNameStream.str().c_str() );
    }
#endif


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

//
//
//
void
AtlasMeshSegmenterConsole
::HandleEMSegmenterEvent( itk::Object* object, const itk::EventObject & event )
{

  // Take care of progress bar
  std::ostringstream  progressLabelStream;
  progressLabelStream << "EM iteration " << m_Driver->GetAtlasMeshSegmenter()->GetEMSegmenter()->GetIterationNumber();
  m_ProgressLabel = progressLabelStream.str();
  m_Progress->label( m_ProgressLabel.c_str() );
  const float  progress = static_cast< float >( m_Driver->GetAtlasMeshSegmenter()->GetEMSegmenter()->GetIterationNumber() ) /
                          static_cast< float >( m_Driver->GetAtlasMeshSegmenter()->GetEMSegmenter()->GetMaximumNumberOfIterations() );
  m_Progress->value( progress * 100 );
  m_Progress->redraw();


  if ( typeid( event ) == typeid( itk::EndEvent ) )
  {
    // Show what we have in the viewers
    m_LabelSelecter->do_callback();
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
void
AtlasMeshSegmenterConsole
::SelectTriangleContainingPoint( float x, float y )
{

  unsigned long  triangleId;
  if ( m_ImageViewer->GetTriangleContainingPoint( x, y, triangleId ) )
  {
    // Set the viewers accordingly
    m_ImageViewer->SetTriangleIdToHighlight( triangleId );
    m_PriorViewer->SetTriangleIdToHighlight( triangleId );
    m_PosteriorViewer->SetTriangleIdToHighlight( triangleId );
    m_ReconstructionViewer->SetTriangleIdToHighlight( triangleId );

    // Retrieve the point id of the nearest vertex clicked
    unsigned long  nearestVertexPointId = m_ImageViewer->GetNearestVertexPointId( x, y );

    // Display the approximative Gaussian in the viewers
    m_ImageViewer->SetPointIdToShowGaussianApproximation( nearestVertexPointId );
    m_PriorViewer->SetPointIdToShowGaussianApproximation( nearestVertexPointId );
    m_PosteriorViewer->SetPointIdToShowGaussianApproximation( nearestVertexPointId );
    m_ReconstructionViewer->SetPointIdToShowGaussianApproximation( nearestVertexPointId );

    // Redraw
    m_ImageViewer->redraw();
    m_PriorViewer->redraw();
    m_PosteriorViewer->redraw();
    m_ReconstructionViewer->redraw();
    Fl::check();
  }

}


//
//
//
void
AtlasMeshSegmenterConsole
::SelectCurrentTriangleContainingPoint( float x, float y )
{

  unsigned long  triangleId;
  if ( m_ImageViewer->GetTriangleContainingPoint( x, y, triangleId ) )
  {
    // Set the viewers accordingly
    m_ImageViewer->SetTriangleIdToHighlight( triangleId );
    m_PriorViewer->SetTriangleIdToHighlight( triangleId );
    m_PosteriorViewer->SetTriangleIdToHighlight( triangleId );
    m_ReconstructionViewer->SetTriangleIdToHighlight( triangleId );

    // Redraw
    m_ImageViewer->redraw();
    m_PriorViewer->redraw();
    m_PosteriorViewer->redraw();
    m_ReconstructionViewer->redraw();
    Fl::check();
  }

}
#endif


//
//
//
void
AtlasMeshSegmenterConsole
::Interrupt()
{
  m_Interrupted = true;

}


//
//
//
void
AtlasMeshSegmenterConsole
::Step()
{
  m_Stepping = true;
  this->Continue();

}



//
//
//
void
AtlasMeshSegmenterConsole
::Continue()
{
  m_Interrupted = false;

}



//
//
//
void
AtlasMeshSegmenterConsole
::ApplyAffineParameters()
{
#if 0
  // Construct the required transform from the GUI elements
  typedef itk::Affine3DTransform< float >   TransformType;
  TransformType::ParametersType  parameters( 12 );
  parameters[ 0 ] = m_RotX->value() / 180 * 3.14159265358979;
  parameters[ 1 ] = m_RotY->value() / 180 * 3.14159265358979;
  parameters[ 2 ] = m_RotZ->value() / 180 * 3.14159265358979;
  parameters[ 3 ] = m_TransX->value();
  parameters[ 4 ] = m_TransY->value();
  parameters[ 5 ] = m_TransZ->value();
  parameters[ 6 ] = m_ScaleX->value();
  parameters[ 7 ] = m_ScaleY->value();
  parameters[ 8 ] = m_ScaleZ->value();
  parameters[ 9 ] = m_SkewX->value();
  parameters[ 10 ] = m_SkewY->value();
  parameters[ 11 ] = m_SkewZ->value();


  // Check if we really need to do anything
  bool  isIdentityTransform = true;
  for ( int i = 0; i < 12; i++ )
  {
    if ( ( i != 6 ) && ( i != 7 ) && ( i != 8 ) )
    {
      if ( parameters[ i] != 0 )
      {
        isIdentityTransform = false;
        break;
      }
    }
    else
    {
      if ( parameters[ i ] != 1 )
      {
        isIdentityTransform = false;
        break;
      }
    }
  }

  if ( isIdentityTransform )
  {
    return;
  }



  //
  if ( ( parameters[ 6 ] != 1 ) || ( parameters[ 7 ] != 1 ) || ( parameters[ 8 ] != 1 ) ||
       ( parameters[ 9 ] != 1 ) || ( parameters[ 10 ] != 1 ) || ( parameters[ 11 ] != 1 ) )
  {
    std::cout << "The deformation field model is going to want to undo the non-rigid components "
              << " of your transform!" << std::endl;
  }

  TransformType::Pointer  transform = TransformType::New();
  transform->SetParameters( parameters );
  transform->Print( std::cout );


  // Make sure the mesh doesn't fall outside the image area. If it does, the rasterizor
  // will draw on invalid memory locations, so let's explicitly forbid that.
  const float  minimalAllowedCoordinate[] = { 0.0f, 0.0f, 0.0f };
  float  maximalAllowedCoordinate[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    maximalAllowedCoordinate[ i ] = m_Driver->GetAtlasMeshSegmenter()->GetImage()->GetBufferedRegion().GetSize()[ i ] - 1;
  }
  m_Driver->GetAtlasMeshSegmenter()->GetImage()->Print( std::cout );
  std::vector< TransformType::InputPointType >  meshCornerPoints;
  TransformType::InputPointType  corner;
  corner[ 2 ] = 0;
  corner[ 0 ] = 0;
  corner[ 1 ] = 0;
  meshCornerPoints.push_back( corner );
  corner[ 0 ] = m_OriginalMeshSize[ 0 ]-1;
  corner[ 1 ] = 0;
  meshCornerPoints.push_back( corner );
  corner[ 0 ] = m_OriginalMeshSize[ 0 ]-1;
  corner[ 1 ] = m_OriginalMeshSize[ 1 ]-1;
  meshCornerPoints.push_back( corner );
  corner[ 0 ] = 0;
  corner[ 1 ] = m_OriginalMeshSize[ 1 ]-1;
  meshCornerPoints.push_back( corner );
  corner[ 2 ] = m_OriginalMeshSize[ 2 ]-1;
  corner[ 0 ] = 0;
  corner[ 1 ] = 0;
  meshCornerPoints.push_back( corner );
  corner[ 0 ] = m_OriginalMeshSize[ 0 ]-1;
  corner[ 1 ] = 0;
  meshCornerPoints.push_back( corner );
  corner[ 0 ] = m_OriginalMeshSize[ 0 ]-1;
  corner[ 1 ] = m_OriginalMeshSize[ 1 ]-1;
  meshCornerPoints.push_back( corner );
  corner[ 0 ] = 0;
  corner[ 1 ] = m_OriginalMeshSize[ 1 ]-1;
  meshCornerPoints.push_back( corner );

  for ( std::vector< TransformType::InputPointType >::const_iterator  it = meshCornerPoints.begin();
        it != meshCornerPoints.end(); ++it )
  {
    TransformType::OutputPointType  mappedCoordinate =  transform->TransformPoint( *it );
    for ( int i = 0; i < 3; i++ )
    {
      if ( ( mappedCoordinate[ i ] < 0  ) || ( mappedCoordinate[ i ] >  maximalAllowedCoordinate[ i ] ) )
      {
        std::cout << "Transform would map mesh outside of image: not allowed." << std::endl;
        this->ResetAffineParameters();
        return;
      }
    }
  }


  std::cout << "Applying affine parameters: " << parameters << std::endl;


  // Create a new mesh that is indentical to the m_OriginalMesh, except that
  // the positions are transformed.
  AtlasMesh::Pointer  mesh = AtlasMesh::New();
  mesh->SetPoints( m_OriginalMesh->GetPoints() );
  mesh->SetCells( m_OriginalMesh->GetCells() );
  mesh->SetPointData( m_OriginalMesh->GetPointData() );
  mesh->SetCellData( m_OriginalMesh->GetCellData() );

  // The positions are the original mesh positions pushed through the transform
  AtlasMesh::PointsContainer::Pointer  transformedPosition = AtlasMesh::PointsContainer::New();
  for ( AtlasMesh::PointsContainer::ConstIterator  it = m_OriginalMesh->GetPoints()->Begin();
        it != m_OriginalMesh->GetPoints()->End(); ++it )
  {
    transformedPosition->InsertElement( it.Index(), transform->TransformPoint( it.Value() ) );
  }
  mesh->SetPoints( transformedPosition );

#if 0
  kvl::AtlasMeshCollection::Pointer  debugCollection = kvl::AtlasMeshCollection::New();
  debugCollection->GenerateFromSingleMesh( mesh, 1, 1000.0f );
  debugCollection->Write( "debug.txt" );
#endif

  // Give the mesh to the segmenter
  m_Segmenter->SetMesh( mesh );

#endif
}


//
//
//
void
AtlasMeshSegmenterConsole
::ResetAffineParameters()
{

  m_RotX->value( 0 );
  m_RotY->value( 0 );
  m_RotZ->value( 0 );
  m_TransX->value( 0 );
  m_TransY->value( 0 );
  m_TransZ->value( 0 );
  m_ScaleX->value( 1 );
  m_ScaleY->value( 1 );
  m_ScaleZ->value( 1 );
  m_SkewX->value( 0 );
  m_SkewY->value( 0 );
  m_SkewZ->value( 0 );

  this->ApplyAffineParameters();

}



//
//
//
void
AtlasMeshSegmenterConsole
::SetGradientDescentStepSize( float stepSize )
{
  // std::cout << "Setting step size to " << stepSize << std::endl;
  // AtlasMeshSegmenter::Pointer  segmenter = const_cast< AtlasMeshSegmenter* >( m_Driver->GetAtlasMeshSegmenter() );
  // segmenter->SetPositionGradientDescentStepSize( stepSize );

}



//
//
//
void
AtlasMeshSegmenterConsole
::SetPositionUpdatingStopCriterion( float stopCriterion )
{
  // std::cout << "Setting position updating stop criterion to " << stopCriterion << std::endl;
  // AtlasMeshSegmenter::Pointer  segmenter = const_cast< AtlasMeshSegmenter* >( m_Driver->GetAtlasMeshSegmenter() );
  // segmenter->SetPositionUpdatingStopCriterion( stopCriterion );

}


//
//
//
void
AtlasMeshSegmenterConsole
::ShowSelectedView()
{
  int  quadrantNumber = 5;
  if ( m_ViewOne->value() )
  {
    quadrantNumber = 1;
  }
  else if ( m_ViewTwo->value() )
  {
    quadrantNumber = 2;
  }
  else if ( m_ViewThree->value() )
  {
    quadrantNumber = 3;
  }
  else if ( m_ViewFour->value() )
  {
    quadrantNumber = 4;
  }

  std::cout << "Showing quadrantNumber: " << quadrantNumber << std::endl;
  m_ImageViewer->LookAt( quadrantNumber );
  m_ReconstructionViewer->LookAt( quadrantNumber );
  m_PriorViewer->LookAt( quadrantNumber );
  m_PosteriorViewer->LookAt( quadrantNumber );

  m_ImageViewer->redraw();
  m_ReconstructionViewer->redraw();
  m_PriorViewer->redraw();
  m_PosteriorViewer->redraw();
  Fl::check();

}


void
AtlasMeshSegmenterConsole
::SetSliceLocation( unsigned int  sagittalSliceNumber,
                    unsigned int  coronalSliceNumber,
                    unsigned int  axialSliceNumber )
{

  m_ImageViewer->SetSliceLocation( sagittalSliceNumber, coronalSliceNumber, axialSliceNumber );
  m_PriorViewer->SetSliceLocation( sagittalSliceNumber, coronalSliceNumber, axialSliceNumber );
  m_PosteriorViewer->SetSliceLocation( sagittalSliceNumber, coronalSliceNumber, axialSliceNumber );
  m_ReconstructionViewer->SetSliceLocation( sagittalSliceNumber, coronalSliceNumber, axialSliceNumber );


  // Redraw
  m_ImageViewer->redraw();
  m_PriorViewer->redraw();
  m_PosteriorViewer->redraw();
  m_ReconstructionViewer->redraw();
  Fl::check();

}


} // end namespace kvl

