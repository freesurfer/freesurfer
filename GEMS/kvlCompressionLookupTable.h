/**
 * @file  kvlCompressionLookupTable.h
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
#ifndef __kvlCompressionLookupTable_h
#define __kvlCompressionLookupTable_h

#include "kvlAtlasMeshCollection.h"
#include "itkImage.h"
#include "itkRGBAPixel.h"


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
  typedef itk::Image< unsigned char, 3 >  CompressedImageType;
  void  Construct( const std::vector< ImageType::ConstPointer >& images );

  //
  bool  Read( const char* fileName );

  //
  bool  Write( const char* fileName );

  //
  CompressedImageType::Pointer  CompressImage( const ImageType*  image ) const;

  // Some typedefs
  typedef std::map< ImageType::PixelType, CompressedImageType::PixelType >  CompressionLookupTableType;
  typedef std::map< CompressedImageType::PixelType, std::string >  LabelStringLookupTableType;
  typedef itk::RGBAPixel< unsigned char >  ColorType;
  typedef std::map< CompressedImageType::PixelType, ColorType >  ColorLookupTableType;

  //
  const CompressionLookupTableType&  GetCompressionLookupTable() const
  {
    return m_CompressionLookupTable;
  }

  //
  const LabelStringLookupTableType&  GetLabelStringLookupTable() const
  {
    return m_LabelStringLookupTable;
  }

  //
  const ColorLookupTableType&  GetColorLookupTable() const
  {
    return m_ColorLookupTable;
  }

protected :
  // Constructor
  CompressionLookupTable();

  // Destructor
  virtual ~CompressionLookupTable();

  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;

private :
  CompressionLookupTable(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Data members
  CompressionLookupTableType  m_CompressionLookupTable;

  LabelStringLookupTableType  m_LabelStringLookupTable;

  ColorLookupTableType  m_ColorLookupTable;


};



} // end namespace kvl


#endif
