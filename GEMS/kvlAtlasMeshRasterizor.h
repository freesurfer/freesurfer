/**
 * @file  kvlAtlasMeshRasterizor.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
 *    $Revision: 1.4 $
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
#ifndef __kvlAtlasMeshRasterizor_h
#define __kvlAtlasMeshRasterizor_h

#include "kvlAtlasMesh.h"
#include "itkImage.h"



namespace kvl
{


template < class TFragmentProcessor >
class AtlasMeshRasterizor : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef AtlasMeshRasterizor  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >   Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshRasterizor, itk::Object );

  /** Some convenient typedefs. */
  typedef TFragmentProcessor   FragmentProcessorType;
  typedef itk::Image< unsigned char, 3 >  LabelImageType;

  /** Get the FragmentProcessor object. */
  FragmentProcessorType& GetFragmentProcessor()
  {
    return m_FragmentProcessor;
  }

  const FragmentProcessorType& GetFragmentProcessor() const
  {
    return m_FragmentProcessor;
  }

  /** Set the FragmentProcessor object. */
  void SetFragmentProcessor(const FragmentProcessorType& FragmentProcessor)
  {
    m_FragmentProcessor = FragmentProcessor;
    this->Modified();
  }

  /** */
  virtual void SetLabelImage( const LabelImageType*  labelImage )
  {
    m_LabelImage = labelImage;
    this->Modified();
  }

  /** */
  const LabelImageType*  GetLabelImage() const
  {
    return m_LabelImage;
  }

  /** */
  virtual void Rasterize( const AtlasMesh* mesh );


  /** */
  void RasterizeTetrahedron( const AtlasMesh::PointType& p0, const AtlasMesh::PointType& p1,  const AtlasMesh::PointType& p2, const AtlasMesh::PointType& p3 );

protected:
  AtlasMeshRasterizor();
  virtual ~AtlasMeshRasterizor() {};


  /** */
  void RasterizeTriangle( const float* vertex0, const float* vertex1, const float* vertex2,
                          const float* pisIn0, const float* pisIn1, const float* pisIn2,
                          const int zLevel );



private:
  AtlasMeshRasterizor(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  FragmentProcessorType  m_FragmentProcessor;
  LabelImageType::ConstPointer  m_LabelImage;


};


} // end namespace kvl

#include "kvlAtlasMeshRasterizor.txx"

#endif
