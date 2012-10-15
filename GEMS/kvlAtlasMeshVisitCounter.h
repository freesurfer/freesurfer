/**
 * @file  kvlAtlasMeshVisitCounter.h
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
#ifndef __kvlAtlasMeshVisitCounter_h
#define __kvlAtlasMeshVisitCounter_h

#include "kvlAtlasMeshRasterizor.h"



namespace kvl
{


namespace FragmentProcessor
{

/**
 *
 */
class CountVisit
{
public:

  typedef itk::Image< unsigned char, 3 >  ImageType;

  CountVisit()
  {
    m_Image = 0;
  }

  ~CountVisit() {};

  void AllocateImage( ImageType::SizeType  size )
  {
    m_Image = ImageType::New();
    m_Image->SetRegions( size );
    m_Image->Allocate();
    m_Image->FillBuffer( 0 );
  }

  const ImageType* GetImage() const
  {
    return m_Image;
  }

  inline void operator()( const float& pi0, const float& pi1, const float& pi2, const float& pi3 )
  {
#if 0
    std::cout << "               Visiting pixel at index " << m_Index << std::endl;
    std::cout << "                       pi0: " << pi0 << std::endl;
    std::cout << "                       pi1: " << pi1 << std::endl;
    std::cout << "                       pi2: " << pi2 << std::endl;
    std::cout << "                       pi3: " << pi3 << std::endl;
#endif


    ( m_Image->GetPixel( m_Index ) )++;

    m_Index[ 0 ]++;
  }

  inline void StartNewSpan( int x, int y, int z, const unsigned char* sourcePointer )
  {
#if 0
    std::cout << "         Starting span (" << x << ", " << y << ", " << z << ")" << std::endl;
#endif
    m_Index[ 0 ] = x;
    m_Index[ 1 ] = y;
    m_Index[ 2 ] = z;
  }

  inline bool StartNewTetrahedron( AtlasMesh::CellIdentifier cellId )
  {
    return true;
  }

  inline void SetMesh( const AtlasMesh* mesh ) {}

private:
  ImageType::Pointer  m_Image;
  ImageType::IndexType  m_Index;
};


} // End namespace FragmentProcessor



/**
 *
 */
class AtlasMeshVisitCounter: public AtlasMeshRasterizor< FragmentProcessor::CountVisit >
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshVisitCounter  Self;
  typedef AtlasMeshRasterizor< FragmentProcessor::CountVisit >  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshVisitCounter, itk::Object );

  /** Some typedefs */
  typedef Superclass::FragmentProcessorType  FragmentProcessorType;
  typedef Superclass::LabelImageType  LabelImageType;
  typedef FragmentProcessorType::ImageType  CountImageType;


  /** */
  virtual void SetLabelImage( const LabelImageType*  labelImage )
  {
    // Use the label image as a template for the count image
    this->GetFragmentProcessor().AllocateImage( labelImage->GetLargestPossibleRegion().GetSize() );

    // Invoke superclass' implementation
    Superclass::SetLabelImage( labelImage );
  }

  /** */
  const CountImageType*  GetCountImage() const
  {
    return this->GetFragmentProcessor().GetImage();
  }

protected:
  AtlasMeshVisitCounter() {};
  virtual ~AtlasMeshVisitCounter() {};

private:
  AtlasMeshVisitCounter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented


};


} // end namespace kvl

#endif
