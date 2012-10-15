/**
 * @file  kvlAtlasMeshViewingConsole.cxx
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
#include "kvlAtlasMeshViewingConsole.h"

#include "FL/Fl.H"
#include "FL/fl_ask.H"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "kvlAtlasMeshSummaryDrawer.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkMinimumMaximumImageCalculator.h"



namespace kvl
{

//
//
//
AtlasMeshViewingConsole
::AtlasMeshViewingConsole()
{

  m_MeshCollection = 0;
  m_BackgroundImage = 0;
  m_Compressor = 0;

}



//
//
//
AtlasMeshViewingConsole
::~AtlasMeshViewingConsole()
{

}




//
//
//
void
AtlasMeshViewingConsole
::LoadMeshCollection( const char*  fileName,
                      const int* templateSize,
                      const std::string& backgroundImageFileName )
{
  // Loaded mesh collection from file
  m_MeshCollection = AtlasMeshCollection::New();
  if ( !m_MeshCollection->Read( fileName ) )
  {
    m_MeshCollection = 0;
    std::ostringstream  alertStream;
    alertStream <<  "Couldn't read mesh collection from file!\n (" << fileName << ")";
    fl_alert( alertStream.str().c_str() );
    exit( -1 );
  }

  // Determine some stuff from the collection
  const unsigned int  numberOfMeshes = m_MeshCollection->GetNumberOfMeshes();
  const unsigned int  numberOfLabels = m_MeshCollection->GetPointParameters()->Begin().Value().m_Alphas.Size();
  ImageType::SizeType  size;
  if ( templateSize )
  {
    for ( int i = 0; i < 3; i++ )
    {
      size[ i ] = templateSize[ i ];
    }
  }
  else
  {
    size[ 0 ] = 0;
    size[ 1 ] = 0;
    size[ 2 ] = 0;
    for ( AtlasMesh::PointsContainer::ConstIterator  it = m_MeshCollection->GetReferencePosition()->Begin();
          it != m_MeshCollection->GetReferencePosition()->End(); ++it )
    {
      for ( int i = 0; i < 3; i++ )
      {
        if ( it.Value()[ i ] > size[ i ] )
        {
          size[ i ] = static_cast< ImageType::SizeType::SizeValueType >( it.Value()[ i ] );
        }
      }
    }
    size[ 0 ]++;
    size[ 1 ]++;
    size[ 2 ]++;
  }
  std::cout << "numberOfMeshes: " << numberOfMeshes << std::endl;
  std::cout << "numberOfLabels: " << numberOfLabels << std::endl;
  std::cout << "size: " << size << std::endl;

  // If background image is given, read it. Otherwise, create an empty image
  if ( backgroundImageFileName.size() )
  {
    std::cout << "Reading background image from file" << std::endl;
    m_BackgroundImage = AtlasMeshViewingConsole::ReadBackgroundImage( backgroundImageFileName );

  }
  else
  {
    std::cout << "Creating zero-filled background image" << std::endl;
    m_BackgroundImage = ImageType::New();
    m_BackgroundImage->SetRegions( size );
    m_BackgroundImage->Allocate();
    m_BackgroundImage->FillBuffer( 0 );

    // Override the default overlayImagelookupTable in the image viewer
    vtkSmartPointer< vtkLookupTable >  overlayImageLookupTable = vtkLookupTable::New();
    overlayImageLookupTable->SetTableRange( 0, 255 );
    overlayImageLookupTable->SetSaturationRange( 0, 0 );
    overlayImageLookupTable->SetHueRange( 0, 0 );
    overlayImageLookupTable->SetValueRange( 0, 1 );
    overlayImageLookupTable->SetAlphaRange( 0, 1 );
    overlayImageLookupTable->Build();
    m_ImageViewer->SetOverlayImageLookupTable( overlayImageLookupTable );

  }
  m_ImageViewer->SetImage( m_BackgroundImage );
  m_ImageViewer->SetScale( 1.0f );


  // Set up the GUI elements
  m_MeshNumber->add( "Reference mesh" );
  for ( unsigned int i = 0; i < numberOfMeshes; i++ )
  {
    std::ostringstream  labelStream;
    labelStream << "Mesh " << i;
    m_MeshNumber->add( labelStream.str().c_str() );
  }
  m_MeshNumber->value( 1 );


  // Try to read the label strings from a lookup table, and fill in the GUI accordingly
  CompressionLookupTable::Pointer  compressor = CompressionLookupTable::New();
  if ( compressor->Read( "compressionLookupTable.txt" ) )
  {
    std::cout << "Found file compressionLookupTable.txt; using it" << std::endl;

    for ( CompressionLookupTable::LabelStringLookupTableType::const_iterator it = compressor->GetLabelStringLookupTable().begin();
          it != compressor->GetLabelStringLookupTable().end(); ++it )
    {
      m_LabelNumber->add( it->second.c_str() );
    }
    m_Compressor = compressor;
  }
  else
  {
    std::cout << "Did not find file compressionLookupTable.txt; using default labels" << std::endl;

    for ( unsigned int i = 0; i < numberOfLabels; i++ )
    {
      std::ostringstream  labelStream;
      labelStream << "Label " << i;
      m_LabelNumber->add( labelStream.str().c_str() );
    }
  }
  m_LabelNumber->value( 0 );


  // Draw
  this->Draw();


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
AtlasMeshViewingConsole
::Show()
{
  m_Window->show();
  Fl::run();

}



//
//
//
void
AtlasMeshViewingConsole
::Draw()
{
  if ( !m_MeshCollection )
  {
    return;
  }

  // Retrieve the correct mesh
  const int  meshNumber = static_cast< int >( m_MeshNumber->value() ) - 1;
  AtlasMesh::ConstPointer  mesh;
  if ( meshNumber < 0 )
  {
    mesh = m_MeshCollection->GetReferenceMesh();
  }
  else
  {
    mesh = m_MeshCollection->GetMesh( meshNumber );
  }


  if ( !m_ShowSummary->value() )
  {
    // Show the alpha image
    AtlasMeshAlphaDrawer::Pointer  alphaDrawer = AtlasMeshAlphaDrawer::New();
    alphaDrawer->SetLabelImage( m_BackgroundImage );
    alphaDrawer->SetLabelNumber( static_cast< unsigned char >( m_LabelNumber->value() ) );
    alphaDrawer->Rasterize( mesh );

    typedef itk::IntensityWindowingImageFilter< AtlasMeshAlphaDrawer::AlphaImageType,
            ImageViewer::ImageType >   WindowerType;
    WindowerType::Pointer  windower = WindowerType::New();
    windower->SetInput( alphaDrawer->GetAlphaImage() );
    windower->SetWindowMinimum( 0 );
    windower->SetWindowMaximum( 1 );
    windower->SetOutputMinimum( 0 );
    windower->SetOutputMaximum( 255 );
    windower->Update();

    m_ImageViewer->SetOverlayImage( windower->GetOutput() );
    m_ImageViewer->SetOverlayScale( 1.0f );
  }
  else
  {
    if ( !m_Compressor )
    {
      // Show the reconstucted image
      AtlasMeshSummaryDrawer::Pointer  summaryDrawer = AtlasMeshSummaryDrawer::New();
      summaryDrawer->SetLabelImage( m_BackgroundImage );
      summaryDrawer->Rasterize( mesh );

      typedef itk::IntensityWindowingImageFilter< AtlasMeshSummaryDrawer::SummaryImageType,
              ImageViewer::ImageType >   WindowerType;
      WindowerType::Pointer  windower = WindowerType::New();
      windower->SetInput( summaryDrawer->GetSummaryImage() );
      windower->SetWindowMinimum( 0 );
      windower->SetWindowMaximum( mesh->GetPointData()->Begin().Value().m_Alphas.Size() );
      windower->SetOutputMinimum( 0 );
      windower->SetOutputMaximum( 255 );
      windower->Update();

      m_ImageViewer->SetOverlayImage( windower->GetOutput() );
      m_ImageViewer->SetOverlayScale( 1.0f );
    }
    else
    {
      std::cout << "We can actually fill in the colors for a nice summary image" << std::endl;

      // Generate an empty RGBA image of the float type
      typedef itk::RGBAPixel< float >  FloatRGBAPixelType;
      typedef itk::Image< FloatRGBAPixelType, 3 >  FloatRGBAImageType;
      FloatRGBAImageType::Pointer  floatRgbaImage = FloatRGBAImageType::New();
      floatRgbaImage->SetRegions( m_BackgroundImage->GetBufferedRegion() );
      floatRgbaImage->Allocate();
      FloatRGBAImageType::PixelType  initialPixelValue( 0.0f );
      //float  tmp[] = { 0.25f, 0.25f, 0.25f, 1.0f };
      //FloatRGBAImageType::PixelType  initialPixelValue( tmp );
      floatRgbaImage->FillBuffer( initialPixelValue );


      // Loop over all entries in the lookup table, rasterize the corresponding class in the
      // mesh, and add to each color channel. The color is weighted by the alpha channel, so
      // it's possible NOT to show some label(s) by editing the lookupTable text file
      for ( CompressionLookupTable::ColorLookupTableType::const_iterator it = m_Compressor->GetColorLookupTable().begin();
            it != m_Compressor->GetColorLookupTable().end(); ++it )
      {
        // Retrieve labelNumber and associated color
        const int  labelNumber = ( *it ).first;
        const kvl::CompressionLookupTable::ColorType  color = ( *it ).second;
        const float  red = color[ 0 ] / 255.0f;
        const float  green = color[ 1 ] / 255.0f;
        const float  blue = color[ 2 ] / 255.0f;
        const float  alpha = color[ 3 ] / 255.0f;
        std::cout << "Label " << labelNumber << " maps to ["
                  << red << "  "
                  << green << "  "
                  << blue << "  "
                  << alpha << "]" << std::endl;

        // Rasterize the prior
        kvl::AtlasMeshAlphaDrawer::Pointer  alphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
        alphaDrawer->SetLabelImage( m_BackgroundImage );
        alphaDrawer->SetLabelNumber( labelNumber );
        alphaDrawer->Rasterize( mesh );
        kvl::AtlasMeshAlphaDrawer::AlphaImageType::ConstPointer  priorImage = alphaDrawer->GetAlphaImage();


        // Now loop over all voxels, and do the Right Thing
        itk::ImageRegionConstIterator< kvl::AtlasMeshAlphaDrawer::AlphaImageType >  priorIt( priorImage,
            priorImage->GetBufferedRegion() );
        itk::ImageRegionIterator< FloatRGBAImageType >  rgbaIt( floatRgbaImage,
            floatRgbaImage->GetBufferedRegion() );

        for ( ; !priorIt.IsAtEnd(); ++priorIt, ++rgbaIt )
        {
          rgbaIt.Value()[ 0 ] += alpha * red * priorIt.Value();
          rgbaIt.Value()[ 1 ] += alpha * green * priorIt.Value();
          rgbaIt.Value()[ 2 ] += alpha * blue * priorIt.Value();
          rgbaIt.Value()[ 3 ] += alpha * priorIt.Value();
        }


      } // End loop over all labels


      // Now convert float RGBA image into uchar RGBA image, as required by the viewer
      typedef ImageViewer::RGBAImageType  RGBAImageType;
      RGBAImageType::Pointer  rgbaImage = RGBAImageType::New();
      rgbaImage->SetRegions( m_BackgroundImage->GetBufferedRegion() );
      rgbaImage->Allocate();

      itk::ImageRegionConstIterator< FloatRGBAImageType >  floatRgbaIt( floatRgbaImage,
          floatRgbaImage->GetBufferedRegion() );
      itk::ImageRegionIterator< RGBAImageType >  rgbaIt( rgbaImage,
          rgbaImage->GetBufferedRegion() );

      for ( ; !floatRgbaIt.IsAtEnd(); ++floatRgbaIt, ++rgbaIt )
      {
        for ( int i = 0; i < 4; i++ )
        {
          rgbaIt.Value()[ i ] = static_cast< unsigned char >( floatRgbaIt.Value()[ i ] * 255 + 0.5 );
        }
      }


      m_ImageViewer->SetOverlayImage( rgbaImage );
    }

  } // End test if we should show a summary image

  // Show the mesh overlaid
  if ( m_ShowMesh->value() )
  {
    m_ImageViewer->SetMesh( mesh );
    //m_ImageViewer->SetEdgeIdToHighlight( m_EdgeIdToHighlight );
  }
  else
  {
    m_ImageViewer->SetMesh( 0 );
  }

  // Redraw
  m_ImageViewer->redraw();
  Fl::check();

}



//
//
//
void
AtlasMeshViewingConsole
::SelectTriangleContainingPoint( float x, float y )
{

//   unsigned long  triangleId;
//   if ( m_ImageViewer->GetTriangleContainingPoint( x, y, triangleId ) )
//     {
//     // Set the viewer accordingly
//     m_ImageViewer->SetTriangleIdToHighlight( triangleId );
//
//     // Retrieve the point id of the nearest vertex clicked
//     unsigned long  nearestVertexPointId = m_ImageViewer->GetNearestVertexPointId( x, y );
//     std::cout << "nearestVertexPointId: " << nearestVertexPointId << std::endl;
//     m_ImageViewer->SetPointIdToShowGaussianApproximation( nearestVertexPointId );
//
//     // Redraw
//     m_ImageViewer->redraw();
//     Fl::check();
//     }

}




//
//
//
void
AtlasMeshViewingConsole
::SetSliceLocation( unsigned int  sagittalSliceNumber,
                    unsigned int  coronalSliceNumber,
                    unsigned int  axialSliceNumber )
{

  m_ImageViewer->SetSliceLocation( sagittalSliceNumber, coronalSliceNumber, axialSliceNumber );

  //this->Draw();
  // Redraw
  m_ImageViewer->redraw();
  Fl::check();


}




//
//
//
void
AtlasMeshViewingConsole
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
AtlasMeshViewingConsole
::GetScreenShot()
{

  m_ImageViewer->WriteScreenShot( "screenshot.png" );

}


//
//
//
void
AtlasMeshViewingConsole
::GetScreenShotSeries()
{

  // Screen shot at reference position
  m_MeshNumber->value( 0 );
  this->Draw();
  m_ImageViewer->WriteScreenShot( "screenshotReferenceMesh.png" );

  // Screen shot at other positions
  for ( unsigned int  meshNumber = 0; meshNumber < m_MeshCollection->GetNumberOfMeshes(); meshNumber++ )
  {
    m_MeshNumber->value( meshNumber + 1 );
    this->Draw();
    std::ostringstream  fileNameStream;
    fileNameStream << "screenshotMesh";
    if ( meshNumber < 100 )
    {
      fileNameStream << "0";
      if ( meshNumber < 10 )
      {
        fileNameStream << "0";
      }
    }
    fileNameStream << meshNumber << ".png";
    m_ImageViewer->WriteScreenShot( fileNameStream.str() );
  }

}



//
//
//
void
AtlasMeshViewingConsole
::DumpImage()
{

  typedef ImageViewer::ImageBaseType ImageBaseType;
  typedef ImageViewer::ImageType  ImageType;
  typedef ImageViewer::RGBAImageType  RGBAImageType;

  ImageBaseType::ConstPointer imageBase = m_ImageViewer->GetOverlayImage();
  if ( !imageBase )
  {
    return;
  }

  // Write out in the correct format
  const std::string  fileName = "imageDump.vtk";
  if ( dynamic_cast< const ImageType* >( imageBase.GetPointer() ) )
  {
    ImageType::ConstPointer castImage = static_cast< const ImageType* >( imageBase.GetPointer() );

    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetFileName( fileName.c_str() );
    writer->SetInput( castImage );
    writer->Update();
  }
  else if ( dynamic_cast< const RGBAImageType* >( imageBase.GetPointer() ) )
  {
    RGBAImageType::ConstPointer castImage = static_cast< const RGBAImageType* >( imageBase.GetPointer() );

    typedef itk::ImageFileWriter< RGBAImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetFileName( fileName.c_str() );
    writer->SetInput( castImage );
    writer->Update();
  }
  else
  {
    std::cerr << "Unsupported image type" << std::endl;
    return;
  }

  std::cout << "Wrote to file " << fileName << std::endl;

}



//
//
//
AtlasMeshViewingConsole::ImageType::Pointer
AtlasMeshViewingConsole
::ReadBackgroundImage( const std::string&  backgroundImageFileName )
{

  // Read the image
  typedef itk::Image< unsigned short, 3 >  InputImageType;
  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileName( backgroundImageFileName.c_str() );
  reader->Update();

  // Unset the origin and spacing because we can't deal with that right now
  const double spacing[] = { 1, 1, 1 };
  const double origin[] = { 0, 0, 0 };
  reader->GetOutput()->SetSpacing( spacing );
  reader->GetOutput()->SetOrigin( origin );


  // Convert to uchar as required by the viewer
  typedef itk::MinimumMaximumImageCalculator< InputImageType >  RangeCalculatorType;
  RangeCalculatorType::Pointer  rangeCalculator = RangeCalculatorType::New();
  rangeCalculator->SetImage( reader->GetOutput() );
  rangeCalculator->Compute();
  typedef itk::IntensityWindowingImageFilter< InputImageType, ImageType >   WindowerType;
  WindowerType::Pointer  windower = WindowerType::New();
  windower->SetInput( reader->GetOutput() );
  windower->SetWindowMinimum( rangeCalculator->GetMinimum() );
  windower->SetWindowMaximum( rangeCalculator->GetMaximum() );
  windower->SetOutputMinimum( 0 );
  windower->SetOutputMaximum( 255 );
  windower->Update();

  return windower->GetOutput();

}




//
//
//
void
AtlasMeshViewingConsole
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


  int  startSliceNumber = 0;
  int  incrementSliceNumber = 1;
  if ( m_InvertOrder->value() )
  {
    startSliceNumber = m_ImageViewer->GetMaximumImageIndex()[ directionNumber ] - 1;
    incrementSliceNumber = -1;
  }

  sliceLocation[ directionNumber ] = startSliceNumber;
  for ( int i = 0;
        i < m_ImageViewer->GetMaximumImageIndex()[ directionNumber ];
        i++, sliceLocation[ directionNumber ] += incrementSliceNumber )
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

