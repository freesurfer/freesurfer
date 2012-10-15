/**
 * @file  kvlSegmentationEvaluationConsole.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:40 $
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
#include "kvlSegmentationEvaluationConsole.h"

#include "FL/Fl.H"
#include "FL/fl_ask.H"
#include "itkIntensityWindowingImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageFileReader.h"
#include "kvlCompressionLookupTable.h"



namespace kvl
{

//
//
//
SegmentationEvaluationConsole
::SegmentationEvaluationConsole()
{

  m_OverlayOpacity->value( m_ImageViewer1->GetOverlayAlpha() );

}



//
//
//
SegmentationEvaluationConsole
::~SegmentationEvaluationConsole()
{

}




//
//
//
void
SegmentationEvaluationConsole
::LoadImage( const char*  fileName )
{
  // Read the image
  typedef itk::Image< unsigned short, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  // Unset the origin and spacing because we can't deal with that right now
  const double spacing[] = { 1, 1, 1 };
  const double origin[] = { 0, 0, 0 };
  reader->GetOutput()->SetSpacing( spacing );
  reader->GetOutput()->SetOrigin( origin );


  // Convert to uchar as required by the viewers
  typedef itk::MinimumMaximumImageCalculator< ImageType >  RangeCalculatorType;
  RangeCalculatorType::Pointer  rangeCalculator = RangeCalculatorType::New();
  rangeCalculator->SetImage( reader->GetOutput() );
  rangeCalculator->Compute();
  typedef itk::IntensityWindowingImageFilter< ImageType, ImageViewer::ImageType >   WindowerType;
  WindowerType::Pointer  windower = WindowerType::New();
  windower->SetInput( reader->GetOutput() );
  windower->SetWindowMinimum( rangeCalculator->GetMinimum() );
  windower->SetWindowMaximum( rangeCalculator->GetMaximum() );
  windower->SetOutputMinimum( 0 );
  windower->SetOutputMaximum( 255 );
  windower->Update();
  m_ImageViewer1->SetImage( windower->GetOutput() );
  m_ImageViewer1->SetScale( 1.0f );
  m_ImageViewer2->SetImage( windower->GetOutput() );
  m_ImageViewer2->SetScale( 1.0f );


  // Set the GUI elements corresponding to the slice locations
  m_SagittalSliceNumber->maximum( m_ImageViewer1->GetMaximumImageIndex()[ 0 ] );
  m_CoronalSliceNumber->maximum( m_ImageViewer1->GetMaximumImageIndex()[ 1 ] );
  m_AxialSliceNumber->maximum( m_ImageViewer1->GetMaximumImageIndex()[ 2 ] );

  unsigned int  sagittalSliceNumber;
  unsigned int  coronalSliceNumber;
  unsigned int  axialSliceNumber;
  m_ImageViewer1->GetSliceLocation( sagittalSliceNumber, coronalSliceNumber, axialSliceNumber );
  m_SagittalSliceNumber->value( sagittalSliceNumber );
  m_CoronalSliceNumber->value( coronalSliceNumber );
  m_AxialSliceNumber->value( axialSliceNumber );

}



//
//
//
void
SegmentationEvaluationConsole
::LoadOverlayImage1( const char*  fileName )
{
  // Read the image
  typedef itk::Image< unsigned char, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  // Unset the origin and spacing because we can't deal with that right now
  const double spacing[] = { 1, 1, 1 };
  const double origin[] = { 0, 0, 0 };
  reader->GetOutput()->SetSpacing( spacing );
  reader->GetOutput()->SetOrigin( origin );

  // Give it to the viewer
  m_ImageViewer1->SetOverlayImage( reader->GetOutput() );

}


//
//
//
void
SegmentationEvaluationConsole
::LoadOverlayImage2( const char*  fileName )
{
  // Read the image
  typedef itk::Image< unsigned char, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  // Unset the origin and spacing because we can't deal with that right now
  const double spacing[] = { 1, 1, 1 };
  const double origin[] = { 0, 0, 0 };
  reader->GetOutput()->SetSpacing( spacing );
  reader->GetOutput()->SetOrigin( origin );

  // Give it to the viewer
  m_ImageViewer2->SetOverlayImage( reader->GetOutput() );

}



//
//
//
void
SegmentationEvaluationConsole
::LoadOverlayLookupTable( const char*  fileName )
{
  // Get the overlay color lookup table
  kvl::CompressionLookupTable::Pointer  compressor = kvl::CompressionLookupTable::New();
  compressor->Read( fileName );
  kvl::CompressionLookupTable::ColorLookupTableType  colors = compressor->GetColorLookupTable();

  // Construct a VTK color lookup table
  vtkLookupTable*  lookupTable = m_ImageViewer1->GetOverlayImageLookupTable();
  m_ImageViewer2->SetOverlayImageLookupTable( lookupTable );
  lookupTable->Print( std::cout );
  kvl::CompressionLookupTable::ColorLookupTableType::const_iterator  it = colors.end();
  --it;
  int  numberOfLabels = ( *it ).first + 1;
  std::cout << "Having " << numberOfLabels << " in the overlay image lookup table" << std::endl;
  lookupTable->SetTableRange( 0, numberOfLabels - 1 );
  lookupTable->SetNumberOfTableValues( numberOfLabels );

  for ( kvl::CompressionLookupTable::ColorLookupTableType::const_iterator  it = colors.begin();
        it != colors.end(); ++it )
  {
    const int  labelNumber = ( *it ).first;
    const kvl::CompressionLookupTable::ColorType  color = ( *it ).second;

    lookupTable->SetTableValue( labelNumber,
                                color[ 0 ] / 255.0, color[ 1 ] / 255.0, color[ 2 ] / 255.0,
                                color[ 3 ] / 255.0 );
    std::cout << "Set table value " << labelNumber << " to ["
              << color[ 0 ] / 255.0 << "  "
              << color[ 1 ] / 255.0 << "  "
              << color[ 2 ] / 255.0 << "  "
              << color[ 3 ] / 255.0 << "]" << std::endl;

  }
  lookupTable->Print( std::cout );

  double   rgb[3];
  for ( int i = 0; i < lookupTable->GetNumberOfTableValues(); i++ )
  {
    lookupTable->GetColor( i, rgb );
    std::cout << i << " maps to [" << rgb[ 0 ] << "  " << rgb[ 1 ] << "  " << rgb[ 2 ] << "]" << std::endl;
  }

}



//
//
//
void
SegmentationEvaluationConsole
::Show()
{
  m_Window->show();
  Fl::run();

}



//
//
//
void
SegmentationEvaluationConsole
::Draw()
{


  // Redraw
  m_ImageViewer1->redraw();
  m_ImageViewer2->redraw();
  Fl::check();

}



//
//
//
void
SegmentationEvaluationConsole
::SetOverlayOpacity( float overlayOpacity )
{
  m_ImageViewer1->SetOverlayAlpha( overlayOpacity );
  m_ImageViewer2->SetOverlayAlpha( overlayOpacity );
  this->Draw();
}



//
//
//
void
SegmentationEvaluationConsole
::SetSliceLocation( unsigned int  sagittalSliceNumber,
                    unsigned int  coronalSliceNumber,
                    unsigned int  axialSliceNumber )
{

  m_ImageViewer1->SetSliceLocation( sagittalSliceNumber, coronalSliceNumber, axialSliceNumber );
  m_ImageViewer2->SetSliceLocation( sagittalSliceNumber, coronalSliceNumber, axialSliceNumber );
  this->Draw();

}




//
//
//
void
SegmentationEvaluationConsole
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
  m_ImageViewer1->LookAt( quadrantNumber );
  m_ImageViewer2->LookAt( quadrantNumber );
  this->Draw();
}



//
//
//
void
SegmentationEvaluationConsole
::Swap()
{
  ImageViewer::ImageBaseType::ConstPointer  tmp = m_ImageViewer1->GetOverlayImage();
  m_ImageViewer1->SetOverlayImage( m_ImageViewer2->GetOverlayImage() );
  m_ImageViewer2->SetOverlayImage( tmp );

  this->Draw();

}



} // end namespace kvl

