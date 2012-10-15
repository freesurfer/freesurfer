/**
 * @file  kvlColorCodeProbabilityImages.cxx
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
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "kvlCompressionLookupTable.h"
#include "itkMGHImageIOFactory.h"


int main( int argc, char* argv[] )
{

  // Sanity check on input
  if ( argc < 3 )
  {
    std::cerr << argv[0] << " lookupTable fileName1 [ fileName2 ... ]" << std::endl;
    return -1;
  }

  // Parse input
  const std::string  lookupTableFileName( argv[ 1 ] );
  std::vector< std::string >  fileNames;
  for ( int i = 2; i < argc; i++ )
  {
    fileNames.push_back( argv[ i ] );
  }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  // Read the lookup table
  kvl::CompressionLookupTable::Pointer  compressor = kvl::CompressionLookupTable::New();
  compressor->Read( lookupTableFileName.c_str() );


  // Loop over all images
  typedef itk::RGBAPixel< unsigned char >  RGBAPixelType;
  typedef itk::Image< RGBAPixelType, 3 >  RGBAImageType;
  RGBAImageType::Pointer  rgbaImage = 0;
  for ( std::vector< std::string >::const_iterator  it = fileNames.begin();
        it != fileNames.end(); ++it )
  {
    //
    const std::string  fileName = *it;
    std::cout << "Handling " << fileName << std::endl;


    // Loop over all entries in the lookup table, and check if the freeSurfer label
    // is part of the fileName
    bool  foundLabel = false;
    unsigned char  label = 0;
    for ( kvl::CompressionLookupTable::LabelStringLookupTableType::const_iterator  labelStringIt
          = compressor->GetLabelStringLookupTable().begin();
          labelStringIt != compressor->GetLabelStringLookupTable().end(); ++labelStringIt )
    {
      if ( fileName.find( labelStringIt->second ) != std::string::npos )
      {
        label = labelStringIt->first;
        std::cout << "found that it corresponds to compressed label " << static_cast< int >( label ) << std::endl;
        foundLabel = true;
        break;
      }
    }

    //
    if ( !foundLabel )
    {
      continue;
    }


    // Determine color
    RGBAPixelType  color = compressor->GetColorLookupTable().find( label )->second;
    std::cout << "it has color [ ";
    for ( int i = 0; i < 4; i++ )
    {
      std::cout << static_cast< int >( color[ i ] ) << " ";
    }
    std::cout << "]" << std::endl;


    // Read the image
    typedef itk::Image< unsigned char, 3 >  ProbabilityImageType;
    typedef itk::ImageFileReader< ProbabilityImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( fileName.c_str() );
    reader->Update();
    ProbabilityImageType::ConstPointer  probabilityImage = reader->GetOutput();

    // Create rgba image if it doesn't exist yet
    if ( !rgbaImage )
    {
      rgbaImage = RGBAImageType::New();
      rgbaImage->SetRegions( probabilityImage->GetBufferedRegion() );
      rgbaImage->Allocate();
      RGBAImageType::PixelType  initialPixelValue( static_cast< unsigned char >( 0 ) );
      rgbaImage->FillBuffer( initialPixelValue );
    }

    // Loop over all voxels
    itk::ImageRegionConstIterator< ProbabilityImageType >   sourceIt( probabilityImage,
        probabilityImage->GetBufferedRegion() );
    itk::ImageRegionIterator< RGBAImageType >  targetIt( rgbaImage,
        rgbaImage->GetBufferedRegion() );
    const float  weights[] = { color[ 0 ] / 255.0f, color[ 1 ] / 255.0f, color[ 2 ] / 255.0f, 1.0f };
    const float  alpha = color[ 3 ] / 255.0f;
    for ( ; !sourceIt.IsAtEnd(); ++sourceIt, ++targetIt )
    {
      for ( int i = 0; i < 4; i++ )
      {
        float  newValue = targetIt.Value()[ i ] +
                          alpha * weights[ i ] * sourceIt.Value();
        if ( newValue > 255 )
        {
          newValue = 255;
        }

        targetIt.Value()[ i ] = static_cast< unsigned char >( newValue + 0.5 );
      }

    } // End loop over all voxels


  } // End loop over all images

  if ( !rgbaImage )
  {
    std::cout << "Found nothing to do" << std::endl;
    return 0;
  }


  // Write result out
  typedef itk::ImageFileWriter< RGBAImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetInput( rgbaImage );
  writer->SetFileName( "colorCodedProbabilies.vtk" );
  writer->Update();

  return 0;
};

