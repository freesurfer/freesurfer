/**
 * @file  kvlAtlasMeshLabelStatisticsCollector.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
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
#ifndef __kvlAtlasMeshLabelStatisticsCollector_h
#define __kvlAtlasMeshLabelStatisticsCollector_h

#include "kvlAtlasMeshRasterizor.h"


namespace kvl
{


namespace FragmentProcessor
{

/**
 *
 */
class CollectLabelStatistics
{
public:

  typedef itk::MapContainer< AtlasMesh::PointIdentifier , AtlasAlphasType >   StatisticsContainerType;

  CollectLabelStatistics()
  {
    m_SourcePointer = 0;
    m_Mesh = 0;
    m_LabelStatistics = 0;
    m_StatisticsInVertex0 = 0;
    m_StatisticsInVertex1 = 0;
    m_StatisticsInVertex2 = 0;
    m_StatisticsInVertex3 = 0;
    m_MinLogLikelihood = 0.0f;
  }

  ~CollectLabelStatistics() {};

  inline void operator()( const float& pi0, const float& pi1, const float& pi2,  const float& pi3 )
  {
    float  epi0 = pi0;
    float  epi1 = pi1;
    float  epi2 = pi2;
    float  epi3 = pi3;

    if ( epi0 < 0 )
    {
      epi0 = 0.0f;
    }
    if ( epi1 < 0 )
    {
      epi1 = 0.0f;
    }
    if ( epi2 < 0 )
    {
      epi2 = 0.0f;
    }
    if ( epi3 < 0 )
    {
      epi3 = 0.0f;
    }
    const float  sumOfEpis = epi0 + epi1 + epi2 + epi3;
    epi0 /= sumOfEpis;
    epi1 /= sumOfEpis;
    epi2 /= sumOfEpis;
    epi3 /= sumOfEpis;



    // Calculate the numerator of the Labelification weights W0, W1, and W2 in this point
    float  W0 = m_AlphasInVertex0[ *m_SourcePointer ] * epi0;
    float  W1 = m_AlphasInVertex1[ *m_SourcePointer ] * epi1;
    float  W2 = m_AlphasInVertex2[ *m_SourcePointer ] * epi2;
    float  W3 = m_AlphasInVertex3[ *m_SourcePointer ] * epi3;

    // Calculate the sum
#if 1
    float  denominator = W0 + W1 + W2 + W3 + 2.2204e-16 /* 0.00001 */ /* 0.001 */;
    m_MinLogLikelihood -= log( denominator );
#else
    float  denominator = W0 + W1 + W2 + W3;
    if ( denominator <= 0 )
    {
      denominator = 2.2204e-16 /* 0.00001 */ /* 0.001 */;
      // Don't contribute to the min log likelihood
    }
    else
    {
      m_MinLogLikelihood -= log( denominator );
    }
#endif

    // Normalize to obtain W0, W1, W2, and W3
    W0 /= denominator;
    W1 /= denominator;
    W2 /= denominator;
    W3 /= denominator;

    // Update the histogram entries in the three vertices accordingly
    ( *m_StatisticsInVertex0 )[ *m_SourcePointer ] += W0;
    ( *m_StatisticsInVertex1 )[ *m_SourcePointer ] += W1;
    ( *m_StatisticsInVertex2 )[ *m_SourcePointer ] += W2;
    ( *m_StatisticsInVertex3 )[ *m_SourcePointer ] += W3;

    // Move on to the next pixel
    m_SourcePointer++;
  }

  inline void StartNewSpan( int x, int y, int z, const unsigned char* sourcePointer )
  {
    m_SourcePointer = sourcePointer;
  }

  inline bool StartNewTetrahedron( AtlasMesh::CellIdentifier cellId )
  {
    // Cache relevant elements of the vertices of this triangle
    AtlasMesh::CellAutoPointer  cell;
    m_Mesh->GetCell( cellId, cell );

    AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
    m_AlphasInVertex0 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_StatisticsInVertex0 = &( m_LabelStatistics->ElementAt( *pit ) );
#if 0
    AtlasMesh::PointType  p0;
    m_Mesh->GetPoint( *pit, &p0 );
#endif
    ++pit;
    m_AlphasInVertex1 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_StatisticsInVertex1 = &( m_LabelStatistics->ElementAt( *pit ) );
#if 0
    AtlasMesh::PointType  p1;
    m_Mesh->GetPoint( *pit, &p1 );
#endif
    ++pit;
    m_AlphasInVertex2 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_StatisticsInVertex2 = &( m_LabelStatistics->ElementAt( *pit ) );
#if 0
    AtlasMesh::PointType  p2;
    m_Mesh->GetPoint( *pit, &p2 );
#endif
    ++pit;
    m_AlphasInVertex3 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_StatisticsInVertex3 = &( m_LabelStatistics->ElementAt( *pit ) );
#if 0
    AtlasMesh::PointType  p3;
    m_Mesh->GetPoint( *pit, &p3 );
#endif

#if 0
    // Check the area of this triangle
    const float twiceArea = ( p1[0] - p0[0] ) * ( p2[1] - p0[1] ) - ( p2[0] - p0[0] ) * ( p1[1] - p0[1] );

#if 1
    // Add contribution of the area prior to the gradient
    if ( twiceArea < 2*0 )
    {
      m_MinLogLikelihood = itk::NumericTraits< float >::max();
    }
#else
    {
      const float  oneOverSqrtThree = 0.577350269189626;
      const float  j11 = p1[0] - p0[0];
      const float  j21 = p1[1] - p0[1];
      const float  j12 = 2 * oneOverSqrtThree * p2[0] - oneOverSqrtThree * ( p0[0] + p1[0] );
      const float  j22 = 2 * oneOverSqrtThree * p2[1] - oneOverSqrtThree * ( p0[1] + p1[1] );

      const float  w = j11*j11 + j12*j12+ j21*j21 + j22*j22;
      const float  d = j11 * j22 - j12 * j21;
      const float  tmp = sqrt( ( w + 2*d ) * ( w - 2*d ) );
      const float  lambda1 = sqrt( ( w + tmp ) / 2 );
      const float  lambda2 = sqrt( ( w - tmp ) / 2 );

      //if ( (  lambda1 < 1 ) || ( lambda2 < 1 ) )
      if ( ( twiceArea < 0 ) || ( fabs( log( lambda1/(lambda2+0.0001) + 0.001 ) ) > 4 ) )
      {
        m_MinLogLikelihood = itk::NumericTraits< float >::max();
      }
    }
#endif
#endif

    return true;
  }

  inline void SetMesh( const AtlasMesh* mesh )
  {
    m_Mesh = mesh;

    // Create a container to hold the statistics
    unsigned int  numberOfLabeles = m_Mesh->GetPointData()->Begin().Value().m_Alphas.Size();
    AtlasAlphasType  zeroEntry( numberOfLabeles );
    zeroEntry.Fill( 0.0f );

    m_LabelStatistics = StatisticsContainerType::New();
    AtlasMesh::PointDataContainer::ConstIterator  it = m_Mesh->GetPointData()->Begin();
    while ( it != m_Mesh->GetPointData()->End() )
    {
      m_LabelStatistics->InsertElement( it.Index(), zeroEntry );

      ++it;
    }

    m_MinLogLikelihood = 0.0f;
  }

  const StatisticsContainerType* GetLabelStatistics() const
  {
    return m_LabelStatistics;
  }

  float GetMinLogLikelihood() const
  {
    return m_MinLogLikelihood;
  }

private:

  const unsigned char*  m_SourcePointer;

  AtlasAlphasType  m_AlphasInVertex0;
  AtlasAlphasType  m_AlphasInVertex1;
  AtlasAlphasType  m_AlphasInVertex2;
  AtlasAlphasType  m_AlphasInVertex3;

  AtlasAlphasType*  m_StatisticsInVertex0;
  AtlasAlphasType*  m_StatisticsInVertex1;
  AtlasAlphasType*  m_StatisticsInVertex2;
  AtlasAlphasType*  m_StatisticsInVertex3;

  AtlasMesh::ConstPointer  m_Mesh;
  StatisticsContainerType::Pointer  m_LabelStatistics;

  float  m_MinLogLikelihood;

};





} // End namespace FragmentProcessor





/**
 *
 */
class AtlasMeshLabelStatisticsCollector: public AtlasMeshRasterizor< FragmentProcessor::CollectLabelStatistics >
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshLabelStatisticsCollector  Self;
  typedef AtlasMeshRasterizor< FragmentProcessor::CollectLabelStatistics >  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshLabelStatisticsCollector, itk::Object );

  /** Some typedefs */
  typedef Superclass::FragmentProcessorType  FragmentProcessorType;
  typedef Superclass::LabelImageType  LabelImageType;
  typedef FragmentProcessorType::StatisticsContainerType  StatisticsContainerType;

  /** */
  const StatisticsContainerType*  GetLabelStatistics() const
  {
    return this->GetFragmentProcessor().GetLabelStatistics();
  }

  /** */
  float GetMinLogLikelihood() const
  {
    return this->GetFragmentProcessor().GetMinLogLikelihood();
  }

protected:
  AtlasMeshLabelStatisticsCollector() {};
  virtual ~AtlasMeshLabelStatisticsCollector() {};

private:
  AtlasMeshLabelStatisticsCollector(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented


};


} // end namespace kvl

#endif

