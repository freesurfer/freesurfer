/**
 * @file  kvlCompressionLookupTable.cxx
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
#include "kvlCompressionLookupTable.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include <fstream>


namespace kvl
{


//
//
//
CompressionLookupTable
::CompressionLookupTable()
{
}



//
//
//
CompressionLookupTable
::~CompressionLookupTable()
{
}




//
//
//
void
CompressionLookupTable
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
}



//
//
//
void
CompressionLookupTable
::Construct( const std::vector< ImageType::ConstPointer >& images )
{
  // Clear whatever we may have
  m_CompressionLookupTable.clear();
  m_LabelStringLookupTable.clear();
  m_ColorLookupTable.clear();

  // Construct a compression lookup table by looping over all images, and creating a
  // new entry for every new intensity encountered
  for ( std::vector< ImageType::ConstPointer >::const_iterator it = images.begin();
        it != images.end(); ++it )
  {

    // Loop over all pixels
    for ( itk::ImageRegionConstIterator< ImageType >  voxelIt( *it, ( *it )->GetLargestPossibleRegion() );
          !voxelIt.IsAtEnd(); ++voxelIt )
    {
      // Try to see if we've met this label before. If not, create a new lookup-table entry
      if ( m_CompressionLookupTable.find( voxelIt.Get() ) == m_CompressionLookupTable.end() )

      {
        const unsigned char  newLabel = m_CompressionLookupTable.size();
        m_CompressionLookupTable[ voxelIt.Get() ] = newLabel;
        std::cout << "    Creating a new lookup table entry: " << static_cast< int >( voxelIt.Get() )
                  << "  ->  " << static_cast< int >( newLabel ) << std::endl;
      }
    }
  }


  // Read the labels from FreeSurferColorLUT.txt if that file is found, and
  // figure out their names and RGBA colors
  const std::string  colorLUTFileName = "FreeSurferColorLUT.txt";
  std::ifstream  in( colorLUTFileName.c_str() );
  if ( in.good() )
  {
    std::cout << "Found FreeSurferColorLUT.txt; reading its relevant contents" << std::endl;

    const int size = 255;
    char buffer[ size ];
    while ( in.getline( buffer, size ) )
    {
      if ( ( buffer[0] == '#' ) || ( buffer[0] == 10 ) || ( buffer[0] == 13 ) ) // '\n' corresponds to 10 or 13, depending on platform
      {
        //std::cout << "    skipping the line: " << buffer << std::endl;
        continue;
      }

      // Parse the line
      const std::string  line = buffer;
      std::istringstream  lineStream( line );
      unsigned int  labelNumber;
      std::string  labelString;
      unsigned int  R;
      unsigned int  G;
      unsigned int  B;
      lineStream >> labelNumber >> labelString >> R >> G >> B;
      unsigned int  A = 255;
      if ( labelNumber == 0 )
      {
        A = 0;
      }


      // If this label is one that is present in our images, add it with its labelString and color to our lookup tables
      if ( m_CompressionLookupTable.find( labelNumber ) != m_CompressionLookupTable.end() )
      {
        // We can't have "/" in our strings as FLTK interprets this is a submenu in our pulldown menu...
        std::string::size_type loc = labelString.find( "/", 0 );
        if ( loc != std::string::npos )
        {
          //std::cout << "Found a / symbol. Replacing it with a - symbol" << std::endl;
          labelString.replace( loc, 1, "-" );
        }

        ColorType  color;
        color[ 0 ] = static_cast< ColorType::ComponentType >( R );
        color[ 1 ] = static_cast< ColorType::ComponentType >( G );
        color[ 2 ] = static_cast< ColorType::ComponentType >( B );
        color[ 3 ] = static_cast< ColorType::ComponentType >( A );

        std::cout << "     found that labelNumber " << labelNumber << " corresponds to labelString "
                  << labelString << " with color " << R << " " << G << " " << B << " " << A << std::endl;

        std::cout << "         Adding " << labelString << " for "
                  << static_cast< int >( ( *( m_CompressionLookupTable.find( labelNumber ) ) ).second ) << std::endl;
        m_LabelStringLookupTable[ ( *( m_CompressionLookupTable.find( labelNumber ) ) ).second ] = labelString;

        std::cout << "         Adding [";
        for ( int i = 0; i < 4; i++ )
        {
          std::cout << static_cast< int >( color[ i ] ) << " ";
        }
        std::cout << "] for "
                  << static_cast< int >( ( *( m_CompressionLookupTable.find( labelNumber ) ) ).second ) << std::endl;
        m_ColorLookupTable[ ( *( m_CompressionLookupTable.find( labelNumber ) ) ).second ] = color;
      }

    } // End loop over all lines in the colorLUT file

  }
  else
  {
    // We need to create some default lookup tables for label strings and colors
    for ( CompressedImageType::PixelType  compressedLabel = 0;
          compressedLabel < m_CompressionLookupTable.size(); compressedLabel++)
    {
      std::ostringstream  labelStringStream;
      labelStringStream << "label" << static_cast< int >( compressedLabel );
      m_LabelStringLookupTable[ compressedLabel ] = labelStringStream.str();

      ColorType  color;
      for ( int i = 0; i < 3; i++ )
      {
        color[ i ] = static_cast< ColorType::ComponentType >(
                       ( static_cast< float >( itk::NumericTraits< ColorType::ComponentType >::max() )
                         / ( m_CompressionLookupTable.size() - 1 ) )
                       * compressedLabel );
      }
      color[ 3 ] = itk::NumericTraits< ColorType::ComponentType >::max();
      m_ColorLookupTable[ compressedLabel ] = color;

    }

  }  // End test if colorLUT file exists


}



//
//
//
bool
CompressionLookupTable
::Read( const char* fileName )
{
  std::ifstream  in( fileName );
  if ( !in.is_open() )
  {
    std::cerr << "Couldn't read from file " << fileName << std::endl;
    return false;
  }


  // Clear whatever we may have
  m_CompressionLookupTable.clear();
  m_LabelStringLookupTable.clear();
  m_ColorLookupTable.clear();


  // Loop over all lines in the file
  const int size = 255;
  char buffer[ size ];
  while ( in.getline( buffer, size ) )
  {
    // Parse the line
    const std::string  line = buffer;
    std::istringstream  lineStream( line );
    unsigned int  labelNumber;
    unsigned int  compressedLabelNumber;
    std::string  labelString;
    unsigned int  R;
    unsigned int  G;
    unsigned int  B;
    unsigned int  A;
    lineStream >> labelNumber >> compressedLabelNumber >> labelString >> R >> G >> B >> A;

    ColorType  color;
    color[ 0 ] = static_cast< ColorType::ComponentType >( R );
    color[ 1 ] = static_cast< ColorType::ComponentType >( G );
    color[ 2 ] = static_cast< ColorType::ComponentType >( B );
    color[ 3 ] = static_cast< ColorType::ComponentType >( A );


    // Now fill in the entries in the lookup tables
    m_CompressionLookupTable[ labelNumber ] = compressedLabelNumber;
    m_LabelStringLookupTable[ compressedLabelNumber ] = labelString;
    m_ColorLookupTable[ compressedLabelNumber ] = color;

  }


  return true;
}



//
//
//
bool
CompressionLookupTable
::Write( const char* fileName )
{

  std::ofstream  out( fileName );
  if ( out.bad() )
  {
    std::cerr << "Can't open " << fileName << " for writing." << std::endl;
    return false;
  }

  // Loop over all entries in our lookup table and write out
  for ( CompressionLookupTableType::const_iterator  compressionIt = m_CompressionLookupTable.begin();
        compressionIt != m_CompressionLookupTable.end(); ++compressionIt )
  {
    out << ( *compressionIt ).first << "   "
        << static_cast< int >( ( *compressionIt ).second ) << "  "
        << m_LabelStringLookupTable[ ( *compressionIt ).second ]
        <<  "     ";
    for ( int i = 0; i < 4; i++ )
    {
      out << static_cast< int >( m_ColorLookupTable[ ( *compressionIt ).second ][ i ] ) << " ";
    }
    out << std::endl;
  }

  return true;

}



//
//
//
CompressionLookupTable::CompressedImageType::Pointer
CompressionLookupTable
::CompressImage( const ImageType*  image ) const
{
  // Allocate an empty image to fill in the compressed labels
  CompressedImageType::Pointer  compressedImage = CompressedImageType::New();
  compressedImage->SetRegions( image->GetLargestPossibleRegion() );
  compressedImage->SetSpacing( image->GetSpacing() );
  compressedImage->SetOrigin( image->GetOrigin() );
  compressedImage->Allocate();
  compressedImage->FillBuffer( 0 );


  // Loop over all voxels and fill in
  itk::ImageRegionConstIterator< ImageType >  imageIterator( image,
      image->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< CompressedImageType >  compressedImageIterator( compressedImage,
      image->GetLargestPossibleRegion() );
  for ( ; !imageIterator.IsAtEnd(); ++imageIterator, ++compressedImageIterator )
  {
    // Apply the lookup-table entry
    CompressionLookupTableType::const_iterator it = m_CompressionLookupTable.find( imageIterator.Get() );
    if ( it != m_CompressionLookupTable.end() )
    {
      compressedImageIterator.Set( (*it).second );
      //std::cout << "Setting pixel with intensity " << static_cast< int >( imageIterator.Get() )
      //          << " to " << static_cast< int >( (*it).second ) << std::endl;
    }
    else
    {
      // We actually didn't find the label in our lookup table!
      //std::cerr << "I couldn't find a compressed label for voxel value " << imageIterator.Get() << std::endl;
      compressedImageIterator.Set( 0 );
    }

  } // End loop over all voxels


  return compressedImage;

}



} // end namespace kvl
