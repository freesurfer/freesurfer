/**
 * @file  kvlAtlasMeshSummaryDrawer.h
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
#ifndef __kvlAtlasMeshSummaryDrawer_h
#define __kvlAtlasMeshSummaryDrawer_h

#include "kvlAtlasMeshRasterizor.h"


namespace kvl
{


namespace FragmentProcessor
{

/**
 *
 */
class DrawSummary
{
public:

  typedef itk::Image< float, 3 >  ImageType;

  DrawSummary()
  {
    m_Image = 0;
    m_MeanInVertex0 = 0.0f;
    m_MeanInVertex1 = 0.0f;
    m_MeanInVertex2 = 0.0f;
    m_MeanInVertex3 = 0.0f;
    m_Mesh = 0;
  }

  ~DrawSummary() {};

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
    float  mean = pi0 * m_MeanInVertex0 +
                  pi1 * m_MeanInVertex1 +
                  pi2 * m_MeanInVertex2 +
                  pi3 * m_MeanInVertex3;
    m_Image->SetPixel( m_Index, mean );

    m_Index[ 0 ]++;
  }

  inline void StartNewSpan( int x, int y, int z, const unsigned char* sourcePointer )
  {
    m_Index[ 0 ] = x;
    m_Index[ 1 ] = y;
    m_Index[ 2 ] = z;
#if 0
    if ( z == 49 )
    {
      std::cout << "         Starting span (" << x << ", " << y << ", " << z << ")" << std::endl;
    }
#endif
  }

  inline bool StartNewTetrahedron( AtlasMesh::CellIdentifier cellId )
  {
    // Cache the alpha of the specified class in each of the vertices of this tetrahedron
    AtlasMesh::CellAutoPointer  cell;
    m_Mesh->GetCell( cellId, cell );

    AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();

    const AtlasAlphasType&  alphas0 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    ++pit;
    const AtlasAlphasType&  alphas1 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    ++pit;
    const AtlasAlphasType&  alphas2 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    ++pit;
    const AtlasAlphasType&  alphas3 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;

    m_MeanInVertex0 = dot_product( alphas0, m_Means );
    m_MeanInVertex1 = dot_product( alphas1, m_Means );
    m_MeanInVertex2 = dot_product( alphas2, m_Means );
    m_MeanInVertex3 = dot_product( alphas3, m_Means );

    return true;
  }

  inline void SetMesh( const AtlasMesh* mesh )
  {
    m_Mesh = mesh;

    const int  numberOfLabels = mesh->GetPointData()->Begin().Value().m_Alphas.Size();
    AtlasAlphasType  tmp( numberOfLabels );
    for ( int labelNumber = 0; labelNumber < numberOfLabels; labelNumber++ )
    {
      tmp[ labelNumber ] = labelNumber;
    }
    m_Means = tmp;
  }

private:

  ImageType::Pointer  m_Image;
  ImageType::IndexType  m_Index;

  float  m_MeanInVertex0;
  float  m_MeanInVertex1;
  float  m_MeanInVertex2;
  float  m_MeanInVertex3;

  AtlasAlphasType  m_Means;
  AtlasMesh::ConstPointer  m_Mesh;


};


} // End namespace FragmentProcessor


/**
 *
 */
class AtlasMeshSummaryDrawer: public AtlasMeshRasterizor< FragmentProcessor::DrawSummary >
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshSummaryDrawer  Self;
  typedef AtlasMeshRasterizor< FragmentProcessor::DrawSummary >  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshSummaryDrawer, itk::Object );

  /** Some typedefs */
  typedef Superclass::FragmentProcessorType  FragmentProcessorType;
  typedef Superclass::LabelImageType  LabelImageType;
  typedef FragmentProcessorType::ImageType  SummaryImageType;


  /** */
  virtual void SetLabelImage( const LabelImageType*  labelImage )
  {
    // Use the label image as a template for the alpha image
    this->GetFragmentProcessor().AllocateImage( labelImage->GetLargestPossibleRegion().GetSize() );

    // Invoke superclass' implementation
    Superclass::SetLabelImage( labelImage );
  }

  /** */
  const SummaryImageType*  GetSummaryImage() const
  {
    return this->GetFragmentProcessor().GetImage();
  }

protected:
  AtlasMeshSummaryDrawer() {};
  virtual ~AtlasMeshSummaryDrawer() {};

private:
  AtlasMeshSummaryDrawer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented


};


} // end namespace kvl

#endif

