#include "kvlCompressionLookupTable.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include <fstream>
#include <algorithm>   


namespace kvl
{


//
//
//
CompressionLookupTable
::CompressionLookupTable()
{
  m_NumberOfClasses = 0;
  
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
#if 0
std::unordered_map< CompressionLookupTable::ImageType::PixelType, std::vector< int > >  
#else
std::map< CompressionLookupTable::ImageType::PixelType, std::vector< int > >  
#endif
CompressionLookupTable
::ReadCollapsedLabelFile() const
{
  //
  const std::string  collapsedLabelFile( "collapsedLabels.txt" );

  // Read in list of collapsed labels
#if 0  
  std::unordered_map< ImageType::PixelType, std::vector< int > >  collapsedLabels;  
#else  
  std::map< ImageType::PixelType, std::vector< int > >  collapsedLabels;  
#endif
  std::ifstream  clfs( collapsedLabelFile.c_str() );
  if ( !( clfs.fail() ) )
    {
    std::cout << "Reading collapsedLabelFile: " << collapsedLabelFile << std::endl;

    std::string line;
    unsigned short comp, lab;
    while ( std::getline( clfs, line ) )
      {
      std::ostringstream  inputParserStream;
      inputParserStream << line;
      std::istringstream  inputStream( inputParserStream.str().c_str() );
      ImageType::PixelType  collapsedLabel;
      if ( !( inputStream >>  collapsedLabel ) )
        {
        std::cerr << "Error reading collapsed label file" << std::endl;
#if 0      
        return std::unordered_map< ImageType::PixelType, std::vector< int > >();
#else        
        return std::map< ImageType::PixelType, std::vector< int > >();
#endif        
        }
      ImageType::PixelType  realLabel;
      while ( inputStream >>  realLabel ) 
        { 
        collapsedLabels[ collapsedLabel ].push_back( realLabel );
        }
      } // End parsing line  
    } // End try to read file

  return collapsedLabels;
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
  m_NumberOfClasses = 0;
  
  // Loop over all collapsed labels, if any, adding a new class for 
  // each real label that is encountered while also pointing the
  // collapsed label to these new classes
#if 0  
  std::unordered_map< ImageType::PixelType, std::vector< int > >  
#else  
  std::map< ImageType::PixelType, std::vector< int > >  
#endif  
        collapsedLabels = this->ReadCollapsedLabelFile();
  for ( 
#if 0
        std::unordered_map< ImageType::PixelType, std::vector< int > >::const_iterator 
#else
        std::map< ImageType::PixelType, std::vector< int > >::const_iterator 
#endif
        collapsedIt = collapsedLabels.begin();
        collapsedIt != collapsedLabels.end(); ++collapsedIt )
    {
    m_CompressionLookupTable[ collapsedIt->first ] = std::vector< int >();
    for ( std::vector< int >::const_iterator  labelIt = collapsedIt->second.begin();
          labelIt != collapsedIt->second.end(); ++labelIt )
      {
      if ( m_CompressionLookupTable.find( *labelIt ) == m_CompressionLookupTable.end() )
        {
        std::cout << "Encountered new real label " << *labelIt << std::endl;
      
        const int  newClassNumber = m_NumberOfClasses;
        m_CompressionLookupTable[ *labelIt ] = std::vector< int >( 1, newClassNumber );
        m_CompressionLookupTable[ collapsedIt->first ].push_back( newClassNumber );
        m_NumberOfClasses++;
        }
        
      }  // End loop over all real labels belonging to a certain collapsed label
      
    } // End loop over all collapsed labels  
        
  
  // Also loop over all pixels of all images, and create a new entry for every new intensity encountered
  for ( std::vector< ImageType::ConstPointer >::const_iterator it = images.begin();
         it != images.end(); ++it )
    {
  
    // Loop over all pixels
    for ( itk::ImageRegionConstIterator< ImageType >  voxelIt( *it, ( *it )->GetLargestPossibleRegion() );
          !voxelIt.IsAtEnd(); ++voxelIt )
      {
      if ( m_CompressionLookupTable.find( voxelIt.Get() ) == m_CompressionLookupTable.end() )
        {
        std::cout << "Encountered new real label " << voxelIt.Get() << std::endl;
      
        const int  newClassNumber = m_NumberOfClasses;
        m_CompressionLookupTable[ voxelIt.Get() ] = std::vector< int >( 1, newClassNumber );
        m_NumberOfClasses++;
        }  
      } // End loop over all pixels
    } // End loop over all images

    
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
      if ( ( buffer[0] == 0 ) || ( buffer[0] == '#' ) || 
           ( buffer[0] == 10 ) || ( buffer[0] == 13 ) ) // '\n' corresponds to 10 or 13, depending on platform
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


      // If this label is one that is present in our images, add it with its labelString and color to our lookup tables (unless it's a virtual, collapsed label)
      if ( m_CompressionLookupTable.find( labelNumber ) != m_CompressionLookupTable.end() )
        {
        if ( collapsedLabels.find( labelNumber ) == collapsedLabels.end() )
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
          const int  classNumber =  ( m_CompressionLookupTable.find( labelNumber )->second )[0];
          std::cout << "        -> maps to classNumber: " << classNumber << std::endl;
          m_LabelStringLookupTable[ classNumber ] = labelString;
          m_ColorLookupTable[ classNumber ] = color;
          }
        }  

      } // End loop over all lines in the colorLUT file

    } // End test if we can read FreeSurferColorLUT.txt
    
  // Every class has to have a color and name -- make sure this is the case
  this->FillInMissingNamesAndColors();
  
}



//
//
//
void
CompressionLookupTable
::FillInMissingNamesAndColors()
{
  
  for ( int classNumber = 0; classNumber < m_NumberOfClasses; ++classNumber )
    {
    if ( m_LabelStringLookupTable.find( classNumber ) == m_LabelStringLookupTable.end() ) 
      {
      // Need to add name
      std::ostringstream  labelStringStream;
      labelStringStream << "class_" << static_cast< int >( classNumber );
      m_LabelStringLookupTable[ classNumber ] = labelStringStream.str();
      }
      
    if ( m_ColorLookupTable.find( classNumber ) == m_ColorLookupTable.end() )
      {
      // Need to add color
      ColorType  color;
      for ( int i = 0; i < 3; i++ )
        {
        color[ i ] = static_cast< ColorType::ComponentType >(
                        ( static_cast< float >( itk::NumericTraits< ColorType::ComponentType >::max() )
                          / ( m_NumberOfClasses - 1 ) )
                        * classNumber );
        }
      color[ 3 ] = itk::NumericTraits< ColorType::ComponentType >::max();
      m_ColorLookupTable[ classNumber ] = color;
      }
      
    } // End loop over all classes
    
}




//
//
//
void
CompressionLookupTable
::Construct( int  numberOfClasses )
{
  
  // Clear whatever we may have 
  m_CompressionLookupTable.clear();
  m_LabelStringLookupTable.clear();
  m_ColorLookupTable.clear();
  
  // 
  m_NumberOfClasses = numberOfClasses;
  for ( int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
    { 
    m_CompressionLookupTable[ classNumber ] = std::vector< int >( 1, classNumber );
    }
  
  //
  this->FillInMissingNamesAndColors();
  
  
}  



//
//
//
bool
CompressionLookupTable
::Read( const std::string& fileName )
{
  std::ifstream  in( fileName.c_str() );
  if ( !in.is_open() )
    {
    std::cerr << "Couldn't read from file " << fileName << std::endl;
    return false;
    }


  // Clear whatever we may have 
  m_CompressionLookupTable.clear();
  m_LabelStringLookupTable.clear();
  m_ColorLookupTable.clear();
  m_NumberOfClasses = 0;
    
  // Loop over all lines in the file
  const int size = 255;
  char buffer[ size ];
  while ( in.getline( buffer, size ) )
    {
    // Parse the line
    const std::string  line = buffer;
    std::istringstream  lineStream( line );
    unsigned int  labelNumber;
    unsigned int  classNumber;
    std::string  labelString;
    unsigned int  R;
    unsigned int  G;
    unsigned int  B;
    unsigned int  A;
    lineStream >> labelNumber >> classNumber >> labelString >> R >> G >> B >> A;

    ColorType  color;
    color[ 0 ] = static_cast< ColorType::ComponentType >( R );
    color[ 1 ] = static_cast< ColorType::ComponentType >( G );
    color[ 2 ] = static_cast< ColorType::ComponentType >( B );
    color[ 3 ] = static_cast< ColorType::ComponentType >( A );


    // Now fill in the entries in the lookup tables
    m_CompressionLookupTable[ labelNumber ] = std::vector< int >( 1, classNumber );
    m_LabelStringLookupTable[ classNumber ] = labelString;
    m_ColorLookupTable[ classNumber ] = color;
    m_NumberOfClasses++;
    }

  return true;
}



//
//
//
bool
CompressionLookupTable
::Write( const std::string& fileName ) const
{

  std::ofstream  out( fileName.c_str() );
  if ( out.bad() )
    {
    std::cerr << "Can't open " << fileName << " for writing." << std::endl;
    return false;
    }

  // Loop over all entries in our lookup table and write out. We only write out
  // non-collapsed labels, i.e., labels that have exactly one class number
  for ( CompressionLookupTableType::const_iterator  compressionIt = m_CompressionLookupTable.begin();
        compressionIt != m_CompressionLookupTable.end(); ++compressionIt )
    {
    if ( compressionIt->second.size() != 1 )
      {
      // Skipping
      continue;
      }  

    const int  classNumber = compressionIt->second[ 0 ];
    out << compressionIt->first << "   "
        << classNumber << "  "
        << m_LabelStringLookupTable.find( classNumber )->second
        <<  "     ";
    for ( int i = 0; i < 4; i++ )
      {
      out << static_cast< int >( m_ColorLookupTable.find( classNumber )->second[ i ] ) << " ";
      }
    out << std::endl;
    }

  return true;

}



} // end namespace kvl
