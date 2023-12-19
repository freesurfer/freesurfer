#ifndef __kvlCompressionLookupTable_h
#define __kvlCompressionLookupTable_h

#include "kvlAtlasMeshCollection.h"
#include "itkImage.h"
#include "itkRGBAPixel.h"
#if 0
  #include <unordered_map>
#endif

namespace kvl
{


class CompressionLookupTable: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef CompressionLookupTable  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CompressionLookupTable, itk::Object );

  //
  typedef itk::Image< unsigned short, 3 >  ImageType;
  
  //
  void  Construct( const std::vector< ImageType::ConstPointer >& images );
  
  //
  void  Construct( int  numberOfClasses );

  //
  bool  Read( const std::string& fileName );

  //
  bool  Write( const std::string& fileName ) const;

  //
  typedef itk::RGBAPixel< unsigned char >  ColorType;
  const std::vector< int >&  GetClassNumbers( ImageType::PixelType label ) const
    {
    return m_CompressionLookupTable.find( label )->second;
    }

  const std::vector<ImageType::PixelType> GetLabels() const
    {
      std::vector<ImageType::PixelType> labels;
      CompressionLookupTableType::const_iterator it;
      for (it = m_CompressionLookupTable.begin(); it != m_CompressionLookupTable.end(); ++it)
        labels.push_back(it->first);

      return labels;
    }
    
  const ColorType&  GetColor( int classNumber ) const
    {
    return m_ColorLookupTable.find( classNumber )->second;
    }

  const std::string&  GetLabelName( int classNumber ) const
    {
    return m_LabelStringLookupTable.find( classNumber )->second;
    }
  
  const int  GetNumberOfClasses() const
    {
    return m_NumberOfClasses;
    }  
  
  
protected :
  // Constructor
  CompressionLookupTable();
  
  // Destructor
  virtual ~CompressionLookupTable();
  
  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;

  //
  void  FillInMissingNamesAndColors();
  
private :
  CompressionLookupTable(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  //
#if 0  
  std::unordered_map< ImageType::PixelType, std::vector< int > >   ReadCollapsedLabelFile() const;
#else  
  std::map< ImageType::PixelType, std::vector< int > >   ReadCollapsedLabelFile() const;
#endif
  
  // Some typedefs
#if 0  
  typedef std::unordered_map< ImageType::PixelType, std::vector< int > >  CompressionLookupTableType;
  typedef std::unordered_map< int, std::string >  LabelStringLookupTableType;
  typedef std::unordered_map< int, ColorType >  ColorLookupTableType;
#else
  typedef std::map< ImageType::PixelType, std::vector< int > >  CompressionLookupTableType;
  typedef std::map< int, std::string >  LabelStringLookupTableType;
  typedef std::map< int, ColorType >  ColorLookupTableType;
#endif
  
  // Data members
  int  m_NumberOfClasses;
  CompressionLookupTableType  m_CompressionLookupTable;
  LabelStringLookupTableType  m_LabelStringLookupTable;
  ColorLookupTableType  m_ColorLookupTable;

};



} // end namespace kvl


#endif
