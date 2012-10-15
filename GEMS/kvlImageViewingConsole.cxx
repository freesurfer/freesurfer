/**
 * @file  kvlImageViewingConsole.cxx
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
#include "kvlImageViewingConsole.h"

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
ImageViewingConsole
::ImageViewingConsole()
{

  m_ImageViewer->SetOverlayAlpha( 0.33f );
  m_OverlayOpacity->value( m_ImageViewer->GetOverlayAlpha() );

}



//
//
//
ImageViewingConsole
::~ImageViewingConsole()
{

}




//
//
//
void
ImageViewingConsole
::LoadImage( const char*  fileName )
{
  // Read the image
  typedef itk::Image< unsigned short, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  // Print some info about the image
  std::cout << "Spacing: " << reader->GetOutput()->GetSpacing() << std::endl;
  std::cout << "Origin: " << reader->GetOutput()->GetOrigin() << std::endl;
  std::cout << "Direction: " << std::endl << reader->GetOutput()->GetDirection() << std::endl;

#if 0
  // Unset the origin and spacing because we can't deal with that right now
  const double spacing[] = { 1, 1, 1 };
  const double origin[] = { 0, 0, 0 };
  reader->GetOutput()->SetSpacing( spacing );
  reader->GetOutput()->SetOrigin( origin );
#endif

  // Convert to uchar as required by the viewer
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

}



//
//
//
void
ImageViewingConsole
::LoadOverlayImage( const char*  fileName )
{

  // Try to read image information
  itk::ImageIOBase::Pointer   imageIO =
    itk::ImageIOFactory::CreateImageIO( fileName, itk::ImageIOFactory::ReadMode );
  if ( imageIO.IsNull() )
  {
    std::cerr << "Can't find a suitable filter to import this type of image." << std::endl;
    return;
  }
  imageIO->SetFileName( fileName );
  imageIO->ReadImageInformation();


  // Read the image
  ImageViewer::ImageBaseType::Pointer  image = 0;
  std::cout << "PixelType: " << imageIO->GetPixelTypeAsString( imageIO->GetPixelType() ) << std::endl;
  if ( imageIO->GetPixelType() == itk::ImageIOBase::VECTOR )
  {
    std::cout << "This seems to be an RGBA file" << std::endl;

    typedef ImageViewer::RGBAImageType  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( fileName );
    reader->Update();
    image = reader->GetOutput();
  }
  else
  {
    typedef itk::Image< unsigned char, 3 >  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( fileName );
    reader->Update();
    image = reader->GetOutput();
  }

#if 0
  // Unset the origin and spacing because we can't deal with that right now
  const double spacing[] = { 1, 1, 1 };
  const double origin[] = { 0, 0, 0 };
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
#endif

  // Give it to the viewer
  m_ImageViewer->SetOverlayImage( image );

}


//
//
//
void
ImageViewingConsole
::LoadMesh( const char*  fileName, int meshNumber )
{
  // Read the mesh collection
  AtlasMeshCollection::Pointer  meshCollection = AtlasMeshCollection::New();
  if ( !meshCollection->Read( fileName ) )
  {
    std::cerr << "Couldn't read mesh collection from file " << fileName << std::endl;
    return;
  }

  // Give the correct mesh to the viewer
  if ( ( meshNumber < 0 ) || ( meshNumber >= static_cast< int >( meshCollection->GetNumberOfMeshes() ) ) )
  {
    m_ImageViewer->SetMesh( meshCollection->GetReferenceMesh() );
  }
  else
  {
    m_ImageViewer->SetMesh( meshCollection->GetMesh( meshNumber ) );
  }

}



//
//
//
void
ImageViewingConsole
::LoadOverlayLookupTable( const char*  fileName )
{
  // Get the overlay color lookup table
  kvl::CompressionLookupTable::Pointer  compressor = kvl::CompressionLookupTable::New();
  compressor->Read( fileName );
  kvl::CompressionLookupTable::ColorLookupTableType  colors = compressor->GetColorLookupTable();

  // Construct a VTK color lookup table
  vtkLookupTable*  lookupTable = m_ImageViewer->GetOverlayImageLookupTable();
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
ImageViewingConsole
::Show()
{
  m_Window->show();
  Fl::run();

}



//
//
//
void
ImageViewingConsole
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
ImageViewingConsole
::SetOverlayOpacity( float overlayOpacity )
{
  m_ImageViewer->SetOverlayAlpha( overlayOpacity );
  this->Draw();
}



//
//
//
void
ImageViewingConsole
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
ImageViewingConsole
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
ImageViewingConsole
::GetScreenShot()
{

  m_ImageViewer->WriteScreenShot( "screenshot.png" );

}




//
//
//
void
ImageViewingConsole
::GetScreenShotSeries( int directionNumber )
{

  unsigned int  sliceLocation[]
  = { static_cast< unsigned int >( m_SagittalSliceNumber->value() ),
      static_cast< unsigned int >( m_CoronalSliceNumber->value() ),
      static_cast< unsigned int >( m_AxialSliceNumber->value() )
    };
  std::string  orientationName;
  switch ( directionNumber )
  {
  case 0:
    orientationName = "Sagittal";
    break;
  case 1:
    orientationName = "Coronal";
    break;
  case 2:
    orientationName = "Axial";
    break;
  default:
    std::cerr << "Invalid directionNumber: " << directionNumber << std::endl;
    return;
  }

  sliceLocation[ directionNumber ] = 0;
  for ( int i = 0;
        i < m_ImageViewer->GetMaximumImageIndex()[ directionNumber ];
        i++, sliceLocation[ directionNumber ]++ )
  {
    this->SetSliceLocation( sliceLocation[ 0 ],
                            sliceLocation[ 1 ],
                            sliceLocation[ 2 ] );

    std::ostringstream  fileNameStream;
    fileNameStream << "screenShot" << orientationName << "Slice";
    if ( i < 100 )
    {
      fileNameStream << "0";
      if ( i < 10 )
      {
        fileNameStream << "0";
      }
    }
    fileNameStream << i << ".png";
    m_ImageViewer->WriteScreenShot( fileNameStream.str() );

  } // End loop over all slices

}


} // end namespace kvl

