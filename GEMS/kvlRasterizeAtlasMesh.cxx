/**
 * @file  kvlRasterizeAtlasMesh.cxx
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
#include "kvlAtlasMeshCollection.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMGHImageIOFactory.h"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "itkIntensityWindowingImageFilter.h"




int main( int argc, char* argv[] )
{

  //
  if ( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " meshCollection templateImage [ labelNumber=0 meshNumber=0 ]" << std::endl;
    exit( -1 );
  }


  // Parse the input arguments
  const std::string  meshCollectionFileName( argv[ 1 ] );
  const std::string  templateImageFileName( argv[ 2 ] );

  int labelNumber = 0;
  if ( argc > 3 )
  {
    std::istringstream  labelNumberStream( argv[ 3 ] );
    labelNumberStream >> labelNumber;
  }

  int meshNumber = 0;
  if ( argc > 4 )
  {
    std::istringstream  meshNumberStream( argv[ 4 ] );
    meshNumberStream >> meshNumber;
  }


  // Show what we have
  std::cout << "meshCollection: " << meshCollectionFileName << std::endl;
  std::cout << "templateImage: " << templateImageFileName << std::endl;
  std::cout << "labelNumber: " << labelNumber << std::endl;
  std::cout << "meshNumber: " << meshNumber << std::endl;



  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );


  try
  {
    // Read the mesh collection and get the mesh to rasterize from it
    kvl::AtlasMeshCollection::Pointer  meshCollection = kvl::AtlasMeshCollection::New();
    if ( !meshCollection->Read( meshCollectionFileName.c_str() ) )
    {
      std::cerr << "Could not read mesh collection from " << meshCollectionFileName << std::endl;
      exit( -1 );
    }

    kvl::AtlasMesh::ConstPointer  mesh = 0;
    if ( meshNumber >= 0 )
    {
      mesh = meshCollection->GetMesh( static_cast< unsigned int >( meshNumber ) );
    }
    else
    {
      std::cout << "Using reference mesh..." << std::endl;
      mesh = meshCollection->GetReferenceMesh();
    }


    // Read the template image. We actually only care about the header information, so don't bother
    // fiddling with the correct pixel type
    typedef kvl::AtlasMeshAlphaDrawer::LabelImageType  TemplateImageType;
    typedef itk::ImageFileReader< TemplateImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( templateImageFileName.c_str() );
    reader->Update();
    TemplateImageType::ConstPointer  templateImage = reader->GetOutput();


    // Draw the alpha image
    kvl::AtlasMeshAlphaDrawer::Pointer  alphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
    alphaDrawer->SetLabelImage( templateImage );
    alphaDrawer->SetLabelNumber( static_cast< unsigned char >( labelNumber ) );
    alphaDrawer->Rasterize( mesh );


    // Convert to unsigned char pixel type
    typedef itk::Image< unsigned char, 3 >  ImageType;
    typedef itk::IntensityWindowingImageFilter< kvl::AtlasMeshAlphaDrawer::AlphaImageType,
            ImageType >   WindowerType;
    WindowerType::Pointer  windower = WindowerType::New();
    windower->SetInput( alphaDrawer->GetAlphaImage() );
    windower->SetWindowMinimum( 0 );
    windower->SetWindowMaximum( 1 );
    windower->SetOutputMinimum( 0 );
    windower->SetOutputMaximum( 255 );
    windower->Update();
    ImageType::Pointer  image = windower->GetOutput();


    // Make sure the header information is copied (should really be added to alphaDrawer's functionality...)
    image->SetOrigin( templateImage->GetOrigin() );
    image->SetSpacing( templateImage->GetSpacing() );
    image->SetDirection( templateImage->GetDirection() );


    // Write out result
    const std::string  outputFileName = "rasterized.mgz";
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetFileName( outputFileName.c_str() );
    writer->SetInput( image );
    writer->Update();

    std::cout << "Wrote output to " << outputFileName << std::endl;

  }
  catch ( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
    return -1;
  }



  return 0;
};
