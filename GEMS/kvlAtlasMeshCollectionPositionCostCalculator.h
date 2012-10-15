/**
 * @file  kvlAtlasMeshCollectionPositionCostCalculator.h
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
#ifndef __kvlAtlasMeshCollectionPositionCostCalculator_h
#define __kvlAtlasMeshCollectionPositionCostCalculator_h

#include "kvlAtlasMeshRasterizor.h"
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshVertexProcessor.h"


namespace kvl
{


namespace FragmentProcessor
{

/**
 *
 */
class CalculatePosteriorPositionCurvature
{
public:

  CalculatePosteriorPositionCurvature()
  {
    m_SourcePointer = 0;

    m_X = 0;
    m_Y = 0;
    m_X0 = 0;
    m_Y0 = 0;
    m_X1 = 0;
    m_Y1 = 0;
    m_X2 = 0;
    m_Y2 = 0;
    m_TwiceArea = 0.0f;

    m_Mesh = 0;
    m_PosteriorPositionCurvature = 0;

    m_CurvatureInVertex0 = 0;
    m_CurvatureInVertex1 = 0;
    m_CurvatureInVertex2 = 0;
  }

  ~CalculatePosteriorPositionCurvature() {};

  inline void operator()( const float& pi0, const float& pi1, const float& pi2 )
  {

#if 0
    const unsigned int  centerX = static_cast< unsigned int >( ( m_X0 + m_X1 + m_X2 ) / 3 );
    const unsigned int  centerY = static_cast< unsigned int >( ( m_Y0 + m_Y1 + m_Y2 ) / 3 );
    const unsigned int  width = 0;
    if ( ( m_X < centerX-width ) || ( m_X > centerX+width ) || ( m_Y < centerY-width ) || ( m_Y > centerY+width ) )
    {
      m_X++;
      m_SourcePointer++;
      return;
    }

    std::cout << "m_X: " << m_X << ",  m_Y: " << m_Y << "  m_TwiceArea: " << m_TwiceArea << std::endl;
    std::cout << "       m_X0: " << m_X0 << ", m_Y0: " << m_Y0 << std::endl;
    std::cout << "       m_X1: " << m_X1 << ", m_Y1: " << m_Y1 << std::endl;
    std::cout << "       m_X2: " << m_X2 << ", m_Y2: " << m_Y2 << std::endl;

#endif

    float  epi0 = pi0;
    float  epi1 = pi1;
    float  epi2 = pi2;
    if ( epi0 < 0.0001 )
    {
      epi0 = 0.0001;
    }
    if ( epi1 < 0.0001 )
    {
      epi1 = 0.0001;
    }
    if ( epi2 < 0.0001 )
    {
      epi2 = 0.0001;
    }

    // Collect the data terms for each vertex
    float  dataTerm0 = m_AlphasInVertex0[ *m_SourcePointer ];
    float  dataTerm1 = m_AlphasInVertex1[ *m_SourcePointer ];
    float  dataTerm2 = m_AlphasInVertex2[ *m_SourcePointer ];

    // Calculate the likelihood
    float  likelihood = dataTerm0 * epi0 + dataTerm1 * epi1 + dataTerm2 * epi2 + 0.000001;


    // Contribution to curvature in vertex 0
    const float sqrTwiceArea = pow( m_TwiceArea, 2 ) + 0.00001;
    const float sqrLikelihood = pow( likelihood, 2 );
    Curvature  curvatureContributionToVertex0;
    curvatureContributionToVertex0.m_Curvature_dxdx =
      ( -pow( m_Y2 - m_Y1 , 2 ) +
        pow( dataTerm1 * ( m_Y - m_Y2 ) + dataTerm2 * ( m_Y1 - m_Y ), 2 ) / sqrLikelihood
      ) / sqrTwiceArea;
    curvatureContributionToVertex0.m_Curvature_dxdy =
      ( ( m_Y2 - m_Y1 ) * ( m_X2 - m_X1 ) +
        ( dataTerm1 * ( m_Y - m_Y2 ) + dataTerm2 * ( m_Y1 - m_Y ) ) *
        ( dataTerm1 * ( m_X2 - m_X ) + dataTerm2 * ( m_X - m_X1 ) ) / sqrLikelihood
      ) / sqrTwiceArea;
    curvatureContributionToVertex0.m_Curvature_dydy =
      ( -pow( m_X2 - m_X1 , 2 ) +
        pow( dataTerm1 * ( m_X2 - m_X ) + dataTerm2 * ( m_X - m_X1 ), 2 ) / sqrLikelihood
      ) / sqrTwiceArea;


    // Contribution to curvature in vertex 1
    Curvature  curvatureContributionToVertex1;
    curvatureContributionToVertex1.m_Curvature_dxdx =
      ( -pow( m_Y2 - m_Y0 , 2 ) +
        pow( dataTerm0 * ( m_Y2 - m_Y ) + dataTerm2 * ( m_Y - m_Y0 ), 2 ) / sqrLikelihood
      ) / sqrTwiceArea;
    curvatureContributionToVertex1.m_Curvature_dxdy =
      ( ( m_Y2 - m_Y0 ) * ( m_X2 - m_X0 ) +
        ( dataTerm0 * ( m_Y2 - m_Y ) + dataTerm2 * ( m_Y - m_Y0 ) ) *
        ( dataTerm0 * ( m_X - m_X2 ) + dataTerm2 * ( m_X0 - m_X ) ) / sqrLikelihood
      ) / sqrTwiceArea;
    curvatureContributionToVertex1.m_Curvature_dydy =
      ( -pow( m_X0 - m_X2 , 2 ) +
        pow( dataTerm0 * ( m_X - m_X2 ) + dataTerm2 * ( m_X0 - m_X ), 2 ) / sqrLikelihood
      ) / sqrTwiceArea;


    // Contribution to curvature in vertex 2
    Curvature  curvatureContributionToVertex2;
    curvatureContributionToVertex2.m_Curvature_dxdx =
      ( -pow( m_Y0 - m_Y1 , 2 ) +
        pow( dataTerm0 * ( m_Y - m_Y1 ) + dataTerm1 * ( m_Y0 - m_Y ), 2 ) / sqrLikelihood
      ) / sqrTwiceArea;
    curvatureContributionToVertex2.m_Curvature_dxdy =
      ( -( m_Y0 - m_Y1 ) * ( m_X1 - m_X0 )  +
        ( dataTerm0 * ( m_Y - m_Y1 ) + dataTerm1 * ( m_Y0 - m_Y ) ) *
        ( dataTerm0 * ( m_X1 - m_X ) + dataTerm1 * ( m_X - m_X0 ) ) / sqrLikelihood
      ) / sqrTwiceArea;
    curvatureContributionToVertex2.m_Curvature_dydy =
      ( -pow( m_X1 - m_X0 , 2 ) +
        pow( dataTerm0 * ( m_X1 - m_X ) + dataTerm1 * ( m_X - m_X0 ), 2 ) / sqrLikelihood
      ) / sqrTwiceArea;

#if 0
    if ( ( isnan( curvatureContributionToVertex0.m_Curvature_dxdx ) ||
           isnan( curvatureContributionToVertex0.m_Curvature_dxdy ) ||
           isnan( curvatureContributionToVertex0.m_Curvature_dydy ) ) ||
         ( isnan( curvatureContributionToVertex1.m_Curvature_dxdx ) ||
           isnan( curvatureContributionToVertex1.m_Curvature_dxdy ) ||
           isnan( curvatureContributionToVertex1.m_Curvature_dydy ) ) ||
         ( isnan( curvatureContributionToVertex2.m_Curvature_dxdx ) ||
           isnan( curvatureContributionToVertex2.m_Curvature_dxdy ) ||
           isnan( curvatureContributionToVertex2.m_Curvature_dydy ) ) )
    {
      std::cout << "m_X: " << m_X << ",  m_Y: " << m_Y << std::endl;
      std::cout << "       log( likelihood ): " << log( likelihood ) << std::endl;
      std::cout << "       *m_SourcePointer: " << static_cast< unsigned int >( *m_SourcePointer ) << std::endl;
      std::cout << "       dataTerm0: " << dataTerm0 << std::endl;
      std::cout << "       dataTerm1: " << dataTerm1 << std::endl;
      std::cout << "       dataTerm2: " << dataTerm2 << std::endl;
      std::cout << "       pi0: " << pi0 << std::endl;
      std::cout << "       pi1: " << pi1 << std::endl;
      std::cout << "       pi2: " << pi2 << std::endl;
      std::cout << "       m_X0: " << m_X0 << ", m_Y0: " << m_Y0 << std::endl;
      std::cout << "       m_X1: " << m_X1 << ", m_Y1: " << m_Y1 << std::endl;
      std::cout << "       m_X2: " << m_X2 << ", m_Y2: " << m_Y2 << std::endl;
    }
#endif

    // Add the contributions to the curvatures
    *m_CurvatureInVertex0 += curvatureContributionToVertex0;
    *m_CurvatureInVertex1 += curvatureContributionToVertex1;
    *m_CurvatureInVertex2 += curvatureContributionToVertex2;

    // Move on to the next pixel
    m_X++;
    m_SourcePointer++;
  }

  inline void StartNewSpan( int x, int y, const unsigned char* sourcePointer )
  {
    m_SourcePointer = sourcePointer;
    m_X = x;
    m_Y = y;
  }

  inline void StartNewTriangle( AtlasMesh::CellIdentifier cellId )
  {
    // Cache relevant elements of the vertices of this triangle
    AtlasMesh::CellAutoPointer  cell;
    m_Mesh->GetCell( cellId, cell );

    AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
    AtlasMesh::PointType  p;
    m_Mesh->GetPoint( *pit, &p );
    m_X0 = p[ 0 ];
    m_Y0 = p[ 1 ];
    m_AlphasInVertex0 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_CurvatureInVertex0 = &( m_PosteriorPositionCurvature->ElementAt( *pit ) );
    ++pit;
    m_Mesh->GetPoint( *pit, &p );
    m_X1 = p[ 0 ];
    m_Y1 = p[ 1 ];
    m_AlphasInVertex1 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_CurvatureInVertex1 = &( m_PosteriorPositionCurvature->ElementAt( *pit ) );
    ++pit;
    m_Mesh->GetPoint( *pit, &p );
    m_X2 = p[ 0 ];
    m_Y2 = p[ 1 ];
    m_AlphasInVertex2 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_CurvatureInVertex2 = &( m_PosteriorPositionCurvature->ElementAt( *pit ) );

    // Precalculate twice the area of this triangle
    m_TwiceArea = ( m_X1 - m_X0 ) * ( m_Y2 - m_Y0 ) - ( m_X2 - m_X0 ) * ( m_Y1 - m_Y0 );
  }

  inline void SetMesh( const AtlasMesh* mesh )
  {
    m_Mesh = mesh;

    // Create a container to hold the position curvature
    m_PosteriorPositionCurvature = AtlasPositionCurvatureContainerType::New();

    AtlasMesh::PointsContainer::ConstIterator pointIt = m_Mesh->GetPoints()->Begin();
    while ( pointIt != m_Mesh->GetPoints()->End() )
    {
      m_PosteriorPositionCurvature->InsertElement( pointIt.Index(), Curvature() );
      ++pointIt;
    }

  }

  const AtlasPositionCurvatureContainerType* GetPosteriorPositionCurvature() const
  {
    return m_PosteriorPositionCurvature;
  }

  AtlasPositionCurvatureContainerType* GetPosteriorPositionCurvature()
  {
    return m_PosteriorPositionCurvature;
  }

private:

  const unsigned char*  m_SourcePointer;
  int  m_X;
  int  m_Y;

  float  m_X0;
  float  m_Y0;
  float  m_X1;
  float  m_Y1;
  float  m_X2;
  float  m_Y2;
  float  m_TwiceArea;

  AtlasAlphasType  m_AlphasInVertex0;
  AtlasAlphasType  m_AlphasInVertex1;
  AtlasAlphasType  m_AlphasInVertex2;

  Curvature*  m_CurvatureInVertex0;
  Curvature*  m_CurvatureInVertex1;
  Curvature*  m_CurvatureInVertex2;

  AtlasMesh::ConstPointer  m_Mesh;
  AtlasPositionCurvatureContainerType::Pointer  m_PosteriorPositionCurvature;

};




} // End namespace FragmentProcessor





/**
 *
 */
class PosteriorPositionCurvatureCalculator :
  public AtlasMeshRasterizor< FragmentProcessor::CalculatePosteriorPositionCurvature >
{
public :

  /** Standard class typedefs */
  typedef PosteriorPositionCurvatureCalculator  Self;
  typedef AtlasMeshRasterizor< FragmentProcessor::CalculatePosteriorPositionCurvature >  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( PosteriorPositionCurvatureCalculator, itk::Object );

  /** Some typedefs */
  typedef Superclass::FragmentProcessorType  FragmentProcessorType;
  typedef Superclass::LabelImageType  LabelImageType;

  /** */
  const AtlasPositionCurvatureContainerType* GetPosteriorPositionCurvature() const
  {
    return this->GetFragmentProcessor().GetPosteriorPositionCurvature();
  }

  /** */
  AtlasPositionCurvatureContainerType* GetPosteriorPositionCurvature()
  {
    return this->GetFragmentProcessor().GetPosteriorPositionCurvature();
  }


protected:
  PosteriorPositionCurvatureCalculator() {};
  virtual ~PosteriorPositionCurvatureCalculator() {};

private:
  PosteriorPositionCurvatureCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};







/**
 *
 */
class AtlasMeshCollectionPositionCostCalculator: public AtlasMeshVertexProcessor
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshCollectionPositionCostCalculator  Self;
  typedef AtlasMeshVertexProcessor  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshCollectionPositionCostCalculator, AtlasMeshVertexProcessor );

  /** Some typedefs */
  typedef itk::Image< unsigned char, 2 >  LabelImageType;
  typedef Superclass::VertexNeighboringTriangleInfo  VertexNeighboringTriangleInfo;
  typedef Superclass::VertexNeighborhood VertexNeighborhood;
  typedef Superclass::VertexNeighborhoodsContainerType VertexNeighborhoodsContainerType;

  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages );

  // Get label images
  const std::vector< LabelImageType::ConstPointer >&  GetLabelImages() const
  {
    return m_LabelImages;
  }

  // Get label image
  const LabelImageType*  GetLabelImage( unsigned int labelImageNumber ) const;

  //
  float GetPositionCost()
  {
    int  numberOfProblemsDummy;
    return this->GetPositionCost( numberOfProblemsDummy );
  }

  //
  float GetPositionCost( int& numberOfProblems );


protected:
  AtlasMeshCollectionPositionCostCalculator();
  virtual ~AtlasMeshCollectionPositionCostCalculator();

private:
  AtlasMeshCollectionPositionCostCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::vector< LabelImageType::ConstPointer >  m_LabelImages;

  int  m_NumberOfLabelImages;

};






} // end namespace kvl

#endif

