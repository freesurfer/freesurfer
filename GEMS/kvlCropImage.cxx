/**
 * @file  kvlCropImage.cxx
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
#include "kvlCroppedImageReader.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "kvlAtlasMeshCollection.h"
#include "vnl/vnl_matlab_write.h"
#include "vnl/vnl_matrix.h"
#include <vcl_fstream.h>


int main( int argc, char* argv[] )
{
  if ( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " fileName [ boundingFileName extraFraction downSamplingFactor ]" << std::endl;
    return( -1 );
  }


  kvl::CroppedImageReader::Pointer  reader = kvl::CroppedImageReader::New();
  if ( argc == 2 )
  {
    reader->Read( argv[ 1 ] );
  }
  else if ( argc == 3 )
  {
    reader->Read( argv[ 1 ], argv[ 2 ] );
  }
  else
  {
    std::istringstream  extraFractionStream( argv[ 3 ] );
    float  extraFraction;
    extraFractionStream >> extraFraction;
    reader->SetExtraFraction( extraFraction );

    if ( argc > 4 )
    {
      std::istringstream  downSamplingFactorStream( argv[ 4 ] );
      int  downSamplingFactor;
      downSamplingFactorStream >> downSamplingFactor;
      reader->SetDownSamplingFactor( downSamplingFactor );
    }

    reader->Read( argv[ 1 ], argv[ 2 ] );
  }


  // Write out image
  typedef itk::ImageFileWriter< kvl::CroppedImageReader::ImageType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( reader->GetImage() );
  writer->SetFileName( "croppedImage.img" );
  writer->Write();
  writer->SetFileName( "croppedImage.mhd" );
  writer->Write();

  // Print out transform
  //std::cout << "Transform: " << std::endl;
  //reader->GetTransform()->Print( std::cout );

  // If a boundingFileName was given, write out a transformed bounding box image
  if ( argc > 2 )
  {
    const std::string  boundingFileName( argv[ 2 ] );
    const std::string  transformedBoundingFileNameBase = "transformedBoundingBox";

    // Collect transform mapping indices in the bounding box into indices in the
    // cropped image world coordinates. Since the cropped image always has voxel
    // size 1 and offset 0 as per our design, this is just the transform mapping
    // bounding box indices into the cropped image indices
    typedef kvl::CroppedImageReader::TransformType  TransformType;
    TransformType::Pointer  transform = TransformType::New();
    transform->Compose( reader->GetTransform() );
    if ( argc > 4 )
    {
      std::istringstream  downSamplingFactorStream( argv[ 4 ] );
      int  downSamplingFactor;
      downSamplingFactorStream >> downSamplingFactor;
      transform->Scale( 1 / static_cast< float >( downSamplingFactor ) );
    }

    // OK, now we have the transform. Write the .mat file
    const std::string  matFileName = transformedBoundingFileNameBase + ".mat";
    std::cout << "Writing .mat file: " << matFileName << std::endl;

    TransformType::ParametersType  parameters = transform->GetParameters();
    vnl_matrix< double >  M( 4, 4 );
    for ( unsigned int row = 0; row < 3; row++ )
    {
      for ( unsigned int col = 0; col < 3; col++ )
      {
        M( row, col ) = parameters[ row * 3 + col ];
      }
      M( row, 3 ) = parameters[ 9 + row ];
    }
    //vnl_matlab_filewrite matlabWriter( matFileName.c_str() );
    //matlabWriter.write( matlabMatrix, "M");
    vcl_ofstream  out( matFileName.c_str() );
    vnl_matlab_write( out, (double const * const *)M.data_array(), M.rows(), M.cols(), (char const *)"M");

    // Write also the actual Analyze files
    typedef itk::Image< short, 3 >  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( boundingFileName.c_str() );
    reader->Update();
    const double spacing[] = { 1, 1, 1 };
    const double origin[] = { 0, 0, 0 };
    reader->GetOutput()->SetSpacing( spacing );
    reader->GetOutput()->SetOrigin( origin );
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( reader->GetOutput() );
    writer->SetFileName( transformedBoundingFileNameBase + ".img" );
    writer->Update();

  } // End test if boundingFileName was given

  return 0;
};

