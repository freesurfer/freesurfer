/**
 * @file  kvlImageThresholdingConsole.cxx
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
#include "kvlImageThresholdingConsole.h"

#include "FL/Fl.H"
#include "FL/fl_ask.H"
#include "itkIntensityWindowingImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMaskImageFilter.h"



namespace kvl
{

//
//
//
ImageThresholdingConsole
::ImageThresholdingConsole()
{

  m_OverlayOpacity->value( m_ImageViewer->GetOverlayAlpha() );

  // Set up a special overlay lookup table, that is blue for background
  // and transparent for foreground
  vtkSmartPointer< vtkLookupTable >  overlayImageLookupTable = vtkSmartPointer< vtkLookupTable >::New();
  overlayImageLookupTable->SetTableRange( 0, 255 );
  overlayImageLookupTable->SetHueRange( 0.75, 0.75 );  // The actual color
  overlayImageLookupTable->SetSaturationRange( 1, 1 );
  overlayImageLookupTable->SetValueRange( 1, 1 );
  overlayImageLookupTable->SetAlphaRange( 1, 0 );
  overlayImageLookupTable->Build();
  m_ImageViewer->SetOverlayImageLookupTable( overlayImageLookupTable );


  // Set up the thresholder
  m_Thresholder = ThresholderType::New();

}


//
//
//
ImageThresholdingConsole
::~ImageThresholdingConsole()
{

}




//
//
//
void
ImageThresholdingConsole
::LoadImages( const std::vector< std::string >& fileNames )
{

  // Loop over all file names, read the image and add its label
  // to the GUI
  for ( std::vector< std::string >::const_iterator  it = fileNames.begin();
        it != fileNames.end(); ++it )
  {
    const std::string  fileName = *it;

    // Read the image
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( fileName.c_str() );
    reader->Update();

#if 0
    // Unset the origin and spacing because we can't deal with that right now
    const double spacing[] = { 1, 1, 1 };
    const double origin[] = { 0, 0, 0 };
    reader->GetOutput()->SetSpacing( spacing );
    reader->GetOutput()->SetOrigin( origin );
#endif

    // Remember the image
    m_Images.push_back( reader->GetOutput() );


    // Add file name to GUI elements
    const std::string  fileNameToDisplay = itksys::SystemTools::GetFilenameWithoutExtension( fileName );
    m_ImageToMask->add( fileNameToDisplay.c_str() );
    m_ImageToThreshold->add( fileNameToDisplay.c_str() );

  } // End loop over all image file names


  //
  m_ImageToMask->value( 1 );
  m_ImageToMask->do_callback();
  m_ImageToThreshold->value( 1 );
  m_ImageToThreshold->do_callback();

}




//
//
//
void
ImageThresholdingConsole
::Show()
{
  m_Window->show();
  Fl::run();

}



//
//
//
void
ImageThresholdingConsole
::Draw()
{

  // Redraw
  m_ImageViewer->redraw();
  Fl::check();

}



//
//
//
void
ImageThresholdingConsole
::SetOverlayOpacity( float overlayOpacity )
{
  m_ImageViewer->SetOverlayAlpha( overlayOpacity );
  this->Draw();
}



//
//
//
void
ImageThresholdingConsole
::SetSliceLocation( unsigned int  sagittalSliceNumber,
                    unsigned int  coronalSliceNumber,
                    unsigned int  axialSliceNumber )
{

  m_ImageViewer->SetSliceLocation( sagittalSliceNumber, coronalSliceNumber, axialSliceNumber );
  this->Draw();

}




//
//
//
void
ImageThresholdingConsole
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
  m_ImageViewer->redraw();
  Fl::check();

}



//
//
//
void
ImageThresholdingConsole
::SetImageToMask( int imageToMask )
{

  //
  ImageType::ConstPointer  image = m_Images[ imageToMask ];

  // Convert to uchar as required by the viewer
  typedef itk::MinimumMaximumImageCalculator< ImageType >  RangeCalculatorType;
  RangeCalculatorType::Pointer  rangeCalculator = RangeCalculatorType::New();
  rangeCalculator->SetImage( image );
  rangeCalculator->Compute();
  typedef itk::IntensityWindowingImageFilter< ImageType, ImageViewer::ImageType >   WindowerType;
  WindowerType::Pointer  windower = WindowerType::New();
  windower->SetInput( image );
  windower->SetWindowMinimum( rangeCalculator->GetMinimum() );
  windower->SetWindowMaximum( rangeCalculator->GetMaximum() );
  windower->SetOutputMinimum( 0 );
  windower->SetOutputMaximum( 255 );
  windower->Update();
  m_ImageViewer->SetImage( windower->GetOutput() );
  m_ImageViewer->SetScale( 1.0f );


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

  this->Draw();

}



//
//
//
void
ImageThresholdingConsole
::SetImageToThreshold( int imageToThreshold )
{

  //
  ImageType::ConstPointer  image = m_Images[ imageToThreshold ];

  // Calculate the min and max of the image, and update the GUI accordingly
  typedef itk::MinimumMaximumImageCalculator< ImageType >  RangeCalculatorType;
  RangeCalculatorType::Pointer  rangeCalculator = RangeCalculatorType::New();
  rangeCalculator->SetImage( image );
  rangeCalculator->Compute();

  std::cout << rangeCalculator->GetMinimum() << std::endl;
  std::cout << rangeCalculator->GetMaximum() << std::endl;

  m_LowerThreshold->minimum( rangeCalculator->GetMinimum() );
  m_LowerThreshold->maximum( rangeCalculator->GetMaximum() );
  m_LowerThreshold->value( rangeCalculator->GetMinimum() + 1 );

  m_UpperThreshold->minimum( rangeCalculator->GetMinimum() );
  m_UpperThreshold->maximum( rangeCalculator->GetMaximum() );
  m_UpperThreshold->value( rangeCalculator->GetMaximum() );

  // Connect the image to the threshold filter
  m_Thresholder->SetInput( image );

  m_ImageViewer->SetOverlayImage( m_Thresholder->GetOutput() );
  m_ImageViewer->SetOverlayScale( 1.0f );

  m_LowerThreshold->do_callback();

}



//
//
//
void
ImageThresholdingConsole
::SetThresholds( float lowerThreshold, float upperThreshold )
{

  m_Thresholder->SetLowerThreshold( static_cast< ThresholderType::InputPixelType >( lowerThreshold ) );
  m_Thresholder->SetUpperThreshold( static_cast< ThresholderType::InputPixelType >( upperThreshold ) );

  this->Draw();
}



//
//
//
void
ImageThresholdingConsole
::MaskImage()
{

  // Sanity check
  if ( m_ImageToMask->value() == 0 )
  {
    return;
  }

  // Get image and file name
  const int  imageNumber = ( m_ImageToMask->value() - 1 );
  ImageType::ConstPointer  image = m_Images[ imageNumber ];
  std::string  fileName = m_ImageToMask->text( imageNumber + 1 );
  fileName += "_masked.mgz";

  // Mask the image that's currently selected
  typedef itk::MaskImageFilter< ImageType, ThresholderType::OutputImageType, ImageType >  MaskerType;
  MaskerType::Pointer  masker = MaskerType::New();
  masker->SetInput1( image );
  masker->SetInput2( m_Thresholder->GetOutput() );


  // Write out
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetInput( masker->GetOutput() );
  writer->SetFileName( fileName.c_str() );
  writer->Update();

  std::cout << "Wrote masked image to file " << fileName << std::endl;

}


//
//
//
void
ImageThresholdingConsole
::WriteMask()
{

  // Write out
  const std::string  fileName = "mask.mgz";
  typedef itk::ImageFileWriter< ThresholderType::OutputImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetInput( m_Thresholder->GetOutput() );
  writer->SetFileName( fileName.c_str() );
  writer->Update();

  std::cout << "Wrote mask to file " << fileName << std::endl;

}



//
//
//
void
ImageThresholdingConsole
::ShowInverseMask( bool showInverseMask )
{

  if ( showInverseMask )
  {
    m_ImageViewer->GetOverlayImageLookupTable()->SetAlphaRange( 0, 1 );
  }
  else
  {
    m_ImageViewer->GetOverlayImageLookupTable()->SetAlphaRange( 1, 0 );
  }

  this->Draw();

}


} // end namespace kvl

