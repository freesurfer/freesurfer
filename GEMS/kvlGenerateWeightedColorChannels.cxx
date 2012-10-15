/**
 * @file  kvlGenerateWeightedColorChannels.cxx
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
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "kvlCompressionLookupTable.h"


int main( int argc, char* argv[] )
{

  // Sanity check on input
  if ( argc < 4 )
  {
    std::cerr << argv[0] << " templateImage meshCollection lookupTable" << std::endl;
    return -1;
  }

  // Read the input stuff from file
  typedef itk::Image< unsigned char, 3 >  TemplateImageType;
  typedef itk::ImageFileReader< TemplateImageType >  TemplateReaderType;
  TemplateReaderType::Pointer  templateReader = TemplateReaderType::New();
  templateReader->SetFileName( argv[ 1 ] );
  templateReader->Update();
  TemplateImageType::ConstPointer  templateImage = templateReader->GetOutput();

  kvl::AtlasMeshCollection::Pointer  collection = kvl::AtlasMeshCollection::New();
  collection->Read( argv[ 2 ] );

  kvl::CompressionLookupTable::Pointer  compressor = kvl::CompressionLookupTable::New();
  compressor->Read( argv[ 3 ] );


  // Allocate empty images for red, green, and blue channel
  typedef itk::Image< float, 3 >  ColorChannelImageType;
  std::vector< ColorChannelImageType::Pointer >  colorChannelImages;
  for ( int i = 0; i < 3; i++ )
  {
    ColorChannelImageType::Pointer  colorChannelImage = ColorChannelImageType::New();
    colorChannelImage->SetRegions( templateImage->GetLargestPossibleRegion() );
    colorChannelImage->Allocate();
    colorChannelImage->FillBuffer( 0.0f );

    colorChannelImages.push_back( colorChannelImage );
  }


  // Loop over all entries in the lookup table, rasterize the corresponding class in the
  // mesh, and add to each color channel. The color is weighted by the alpha channel, so
  // it's possible NOT to show some label(s) by editing the lookupTable text file
  for ( kvl::CompressionLookupTable::ColorLookupTableType::const_iterator it = compressor->GetColorLookupTable().begin();
        it != compressor->GetColorLookupTable().end(); ++it )
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
    alphaDrawer->SetLabelImage( templateImage );
    alphaDrawer->SetLabelNumber( labelNumber );
    alphaDrawer->Rasterize( collection->GetReferenceMesh() );
    kvl::AtlasMeshAlphaDrawer::AlphaImageType::ConstPointer  priorImage = alphaDrawer->GetAlphaImage();


    // Now loop over all voxels, and do the Right Thing
    itk::ImageRegionConstIterator< kvl::AtlasMeshAlphaDrawer::AlphaImageType >  priorIt( priorImage,
        priorImage->GetBufferedRegion() );
    itk::ImageRegionIterator< ColorChannelImageType >  redIt( colorChannelImages[ 0 ],
        priorImage->GetBufferedRegion() );
    itk::ImageRegionIterator< ColorChannelImageType >  greenIt( colorChannelImages[ 1 ],
        priorImage->GetBufferedRegion() );
    itk::ImageRegionIterator< ColorChannelImageType >  blueIt( colorChannelImages[ 2 ],
        priorImage->GetBufferedRegion() );

    for ( ; !priorIt.IsAtEnd(); ++priorIt, ++redIt, ++greenIt, ++blueIt )
    {
      redIt.Value() += alpha * red * priorIt.Value();
      greenIt.Value() += alpha * green * priorIt.Value();
      blueIt.Value() += alpha * blue * priorIt.Value();
    }


  } // End loop over all labels


  // Write out red, blue, and green channels
  std::vector< std::string >  fileNames;
  fileNames.push_back( "red.mhd" );
  fileNames.push_back( "green.mhd" );
  fileNames.push_back( "blue.mhd" );
  typedef itk::ImageFileWriter< ColorChannelImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  for ( int i = 0; i < 3; i++ )
  {
    writer->SetInput( colorChannelImages[ i ] );
    writer->SetFileName( fileNames[ i ].c_str() );
    writer->Write();
  }


  return 0;
};

