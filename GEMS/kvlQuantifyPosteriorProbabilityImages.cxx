/**
 * @file  kvlQuantifyPosteriorProbabilityImages.cxx
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
#include "itkImageRegionConstIterator.h"
#include "itkImageFileReader.h"
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
  std::cout << "volumeInVoxels: " << std::endl;
  for ( std::vector< std::string >::const_iterator  it = fileNames.begin();
        it != fileNames.end(); ++it )
  {
    //
    const std::string  fileName = *it;
    //std::cout << "Handling " << fileName << std::endl;


    // Loop over all entries in the lookup table, and check if the freeSurfer label
    // is part of the fileName
    bool  foundLabel = false;
    std::string  labelString;
    for ( kvl::CompressionLookupTable::LabelStringLookupTableType::const_iterator  labelStringIt
          = compressor->GetLabelStringLookupTable().begin();
          labelStringIt != compressor->GetLabelStringLookupTable().end(); ++labelStringIt )
    {
      if ( fileName.find( labelStringIt->second ) != std::string::npos )
      {
        labelString = labelStringIt->second;
        //std::cout << "found that it corresponds to label string " << labelString << std::endl;
        foundLabel = true;
        break;
      }
    }

    //
    if ( !foundLabel )
    {
      continue;
    }


    // Read the image
    typedef itk::Image< unsigned char, 3 >  ProbabilityImageType;
    typedef itk::ImageFileReader< ProbabilityImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( fileName.c_str() );
    reader->Update();
    ProbabilityImageType::ConstPointer  probabilityImage = reader->GetOutput();

    // Loop over all voxels
    itk::ImageRegionConstIterator< ProbabilityImageType >   it( probabilityImage,
        probabilityImage->GetBufferedRegion() );
    float  volumeInVoxels = 0.0f;
    for ( ; !it.IsAtEnd(); ++it )
    {
      volumeInVoxels += it.Value() / 255.0f;
    } // End loop over all voxels


    // Print out result
    std::cout << "    " << labelString << ": " << volumeInVoxels << std::endl;

  } // End loop over all images


  return 0;
};

