/**
 * @file  kvlApplyTransform.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
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
#include "kvlRegisterer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMGHImageIOFactory.h"



typedef kvl::Registerer::ImageType  ImageType;
typedef itk::AffineTransform< double, 3 >  AffineTransformType;

// Get an image's image-to-world transfrom, mapping voxel indices into real-world coordinates
// measured in mm
AffineTransformType::Pointer  GetImageToWorldTransform( const ImageType* image )
{

  AffineTransformType::MatrixType  scale;
  AffineTransformType::OffsetType  offset;
  for ( int i = 0; i < 3; i++ )
  {
    scale[ i ][ i ] = image->GetSpacing()[ i ];
    offset[ i ] = image->GetOrigin()[ i ];
  }
  AffineTransformType::Pointer  imageToWorld = AffineTransformType::New();
  imageToWorld->SetMatrix( image->GetDirection() * scale );
  imageToWorld->SetOffset( offset );

  return imageToWorld;
};



int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 14 )
  {
    std::cerr << argv[ 0 ] << " imageFileName M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 [ invert=0 "
              "interpretAsDesiredVoxelToVoxelTransformToReference_imageFileName otherImageFileName ]" << std::endl;
    return -1;
  }

  // Parse input
  const std::string  imageFileName = argv[ 1 ];
  typedef kvl::Registerer::TransformType  TransformType;
  TransformType::OffsetType  offset;
  TransformType::MatrixType  matrix;
  for ( int i = 0; i < 12; i++ )
  {
    // Get the value
    TransformType::ScalarType  value;
    std::istringstream  valueStream( argv[ i + 2 ] );
    valueStream >> value;

    // Put the value in the correct place
    switch ( i )
    {
    case 0:
      matrix[ 0 ][ 0 ] = value;
      break;
    case 1:
      matrix[ 0 ][ 1 ] = value;
      break;
    case 2:
      matrix[ 0 ][ 2 ] = value;
      break;
    case 3:
      offset[ 0 ] = value;
      break;
    case 4:
      matrix[ 1 ][ 0 ] = value;
      break;
    case 5:
      matrix[ 1 ][ 1 ] = value;
      break;
    case 6:
      matrix[ 1 ][ 2 ] = value;
      break;
    case 7:
      offset[ 1 ] = value;
      break;
    case 8:
      matrix[ 2 ][ 0 ] = value;
      break;
    case 9:
      matrix[ 2 ][ 1 ] = value;
      break;
    case 10:
      matrix[ 2 ][ 2 ] = value;
      break;
    case 11:
      offset[ 2 ] = value;
      break;
    }

  } // End loop over all input


  bool  invert = false;
  if ( argc > 14 )
  {
    std::istringstream  invertStream( argv[ 14 ] );
    invertStream >> invert;
  }



  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  try
  {
    // Read image
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( imageFileName.c_str() );
    reader->Update();
    ImageType::Pointer  image = reader->GetOutput();


    // Construct the required transform
    typedef kvl::Registerer::TransformType  TransformType;
    TransformType::Pointer  transform = TransformType::New();

    // If this is a desired transform to be accomplished from voxel grid to voxel grid,
    // compute the netto world-to-world transform to be applied
    if ( argc > 15 )
    {
      // Read the reference image. Note: it's not actually necessary to read the pixel
      // buffer just to get access to the header inforation in ITK, but I'm too lazy
      // to figure out how to do it otherwise again (something with itkImageIOBase etc)
      const std::string  referenceImageFileName( argv[ 15 ] );
      ReaderType::Pointer  referenceReader = ReaderType::New();
      referenceReader->SetFileName( referenceImageFileName.c_str() );
      referenceReader->Update();
      ImageType::ConstPointer  referenceImage = referenceReader->GetOutput();

      // Construct the image's image-to-world transfrom
      AffineTransformType::Pointer  imageToWorld = GetImageToWorldTransform( image );
      std::cout << "imageToWorld: " << imageToWorld << std::endl;

      // Construct the reference's image-to-world transfrom
      AffineTransformType::Pointer  referenceImageToWorld = GetImageToWorldTransform( referenceImage );
      std::cout << "referenceImageToWorld: " << referenceImageToWorld << std::endl;

      // OK, so we want a transform "X" so that  desiredTransform = inv( referenceImageToWorld )  * X * imageToWorld
      // which means that X = referenceImageToWorld * desiredTransform * inv( imageToWorld )
      AffineTransformType::Pointer  desiredTransform = AffineTransformType::New();
      desiredTransform->SetMatrix( matrix );
      desiredTransform->SetOffset( offset );
      std::cout << "desiredTransform: " << desiredTransform << std::endl;

      // Invert if necessary
      if ( invert )
      {
        AffineTransformType::Pointer  inverseDesiredTransform = AffineTransformType::New();
        desiredTransform->GetInverse( inverseDesiredTransform );
        desiredTransform = inverseDesiredTransform;
      }

      AffineTransformType::Pointer  X = AffineTransformType::New();
      imageToWorld->GetInverse( X );
      X->Compose( desiredTransform );
      X->Compose( referenceImageToWorld );
      std::cout << "X: " << X << std::endl;

      // Copy into the correct format
      transform->SetMatrix( X->GetMatrix() );
      transform->SetOffset( X->GetOffset() );
    }
    else
    {
      transform->SetMatrix( matrix );
      transform->SetOffset( offset );

      // Invert if necessary
      if ( invert )
      {
        TransformType::Pointer  inverseTransform = transform->GetInverse();
        transform = inverseTransform;
      }

    }




    // Show what we have
    std::cout << "transform: " << transform << std::endl;


    // Collect all the images we have to apply the transform to
    std::vector< ImageType::Pointer >  images;
    images.push_back( image );
    std::vector< std::string >  fileNames;
    fileNames.push_back( imageFileName );
    if ( argc > 16 )
    {
      const std::string  otherImageFileName( argv[ 16 ] );
      ReaderType::Pointer  otherReader = ReaderType::New();
      otherReader->SetFileName( otherImageFileName.c_str() );
      otherReader->Update();

      images.push_back( otherReader->GetOutput() );
      fileNames.push_back( otherImageFileName );
    }

    // Apply the transform
    kvl::Registerer::ParametersType  parameters = transform->GetParameters();
    kvl::Registerer::Pointer  registerer = kvl::Registerer::New();
    registerer->SetParameters( parameters );
    for ( int imageNumber = 0; imageNumber < static_cast< int >( images.size() ); imageNumber++ )
    {
      //
      const std::string  fileName = fileNames[ imageNumber ];
      ImageType::Pointer  myImage = images[ imageNumber ];

      //
      registerer->ApplyParameters( myImage );

      // Write out
      std::ostringstream  outputFileNameStream;
      outputFileNameStream << itksys::SystemTools::GetFilenameWithoutExtension( fileName.c_str() )
                           << "_transformed.mgz";
      const std::string  outputFileName = outputFileNameStream.str();
      typedef itk::ImageFileWriter< ImageType >  WriterType;
      WriterType::Pointer  writer = WriterType::New();
      writer->SetInput( myImage );
      writer->SetFileName( outputFileName.c_str() );
      writer->Update();

      std::cout << "Wrote out " << outputFileName << std::endl;
    }

  }
  catch ( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
    return -1;
  }


  return 0;
};
