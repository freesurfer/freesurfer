/**
 * @file  kvlAtlasMeshVertexProcessor.h
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
#ifndef __kvlAtlasMeshVertexProcessor_h
#define __kvlAtlasMeshVertexProcessor_h

#include "kvlAtlasMeshCollection.h"
#include "itkImage.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_matrix_fixed.h"



namespace kvl
{

namespace FragmentProcessor
{

/**
 *
 */
class CalculateDataContributions
{
public:

  CalculateDataContributions()
  {
    m_SourcePointer = 0;

    m_Cost = 0.0f;
    m_Gradient[ 0 ] = 0.0f;
    m_Gradient[ 1 ] = 0.0f;
    m_Gradient[ 2 ] = 0.0f;
  }

  ~CalculateDataContributions() {};

  inline void operator()( const float& pi0, const float& pi1, const float& pi2, const float& pi3 )
  {

    // Collect the data terms for each vertex
    float  alpha0 = m_AlphasInVertex0[ *m_SourcePointer ];
    float  alpha1 = m_AlphasInVertex1[ *m_SourcePointer ];
    float  alpha2 = m_AlphasInVertex2[ *m_SourcePointer ];
    float  alpha3 = m_AlphasInVertex3[ *m_SourcePointer ];

    // Calculate the likelihood and add it's -log to the cost
    float  likelihood = alpha0 * pi0 + alpha1 * pi1 + alpha2 * pi2 + alpha3 * pi3 + 1e-5;
    m_Cost -= log( likelihood );

#if 0
    std::cout << "\n\n=========================" << std::endl;
    std::cout << "*m_SourcePointer : " << static_cast< int >( *m_SourcePointer  ) << std::endl;
    std::cout << "alpha0: " << alpha0 << std::endl;
    std::cout << "alpha1: " << alpha1 << std::endl;
    std::cout << "alpha2: " << alpha2 << std::endl;
    std::cout << "alpha3: " << alpha3 << std::endl;
    std::cout << "pi0: " << pi0 << std::endl;
    std::cout << "pi1: " << pi1 << std::endl;
    std::cout << "pi2: " << pi2 << std::endl;
    std::cout << "pi3: " << pi3 << std::endl;
    std::cout << "likelihood: " << likelihood << std::endl;
    std::cout << "-log( likelihood ): " << -log( likelihood ) << std::endl;
    std::cout << "=========================\n\n" << std::endl;
#endif


    //
    const float  tmpX = ( m_XGradientBasis[ *m_SourcePointer ] ) / likelihood;
    const float  tmpY = ( m_YGradientBasis[ *m_SourcePointer ] ) / likelihood;
    const float  tmpZ = ( m_ZGradientBasis[ *m_SourcePointer ] ) / likelihood;

    // Add contribution to gradient in vertex 0
    AtlasPositionGradientType  gradientContributionToVertex0;
    gradientContributionToVertex0[ 0 ] = tmpX * pi0;
    gradientContributionToVertex0[ 1 ] = tmpY * pi0;
    gradientContributionToVertex0[ 2 ] = tmpZ * pi0;
    m_Gradient += gradientContributionToVertex0;

    // Move on to the next pixel
    m_SourcePointer++;
  }

  inline void StartNewSpan( int x, int y, int z, const unsigned char* sourcePointer )
  {
    m_SourcePointer = sourcePointer;
  }

  inline bool StartNewTetrahedron( AtlasMesh::CellIdentifier cellId ) {}

  bool StartNewTetrahedron( float x0, float y0, float z0,
                            float x1, float y1, float z1,
                            float x2, float y2, float z2,
                            float x3, float y3, float z3,
                            const AtlasAlphasType& alphas0,
                            const AtlasAlphasType& alphas1,
                            const AtlasAlphasType& alphas2,
                            const AtlasAlphasType& alphas3 )
  {
    //
    // Cache relevant elements of the vertices of this tetrahedron. The notation used is Y = [ p0 p1 p2 p3; 1 1 1 1 ]
    //
    const float y11 = x0;
    const float y21 = y0;
    const float y31 = z0;
    m_AlphasInVertex0 = alphas0;

    const float y12 = x1;
    const float y22 = y1;
    const float y32 = z1;
    m_AlphasInVertex1 = alphas1;

    const float y13 = x2;
    const float y23 = y2;
    const float y33 = z2;
    m_AlphasInVertex2 = alphas2;

    const float y14 = x3;
    const float y24 = y3;
    const float y34 = z3;
    m_AlphasInVertex3 = alphas3;




    //
    // Finally, precalculate some stuff that will be used over and over again, each time a voxel is visited
    //

    // Get Gamma, defined as Gamma = [ 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 1] * inv( Y )
    // where Y = [ p0 p1 p2 p3; 1 1 1 1 ]
    vnl_matrix_fixed< float, 4, 4 >  Y;
    Y.put( 0, 0, y11 );
    Y.put( 1, 0, y21 );
    Y.put( 2, 0, y31 );
    Y.put( 3, 0, 1.0f );

    Y.put( 0, 1, y12 );
    Y.put( 1, 1, y22 );
    Y.put( 2, 1, y32 );
    Y.put( 3, 1, 1.0f );

    Y.put( 0, 2, y13 );
    Y.put( 1, 2, y23 );
    Y.put( 2, 2, y33 );
    Y.put( 3, 2, 1.0f );

    Y.put( 0, 3, y14 );
    Y.put( 1, 3, y24 );
    Y.put( 2, 3, y34 );
    Y.put( 3, 3, 1.0f );

    vnl_matrix_fixed< float, 4, 4 >  invY = vnl_inverse( Y );
    const float  gamma11 = invY( 1, 0 );
    const float  gamma12 = invY( 1, 1 );
    const float  gamma13 = invY( 1, 2 );
    const float  gamma21 = invY( 2, 0 );
    const float  gamma22 = invY( 2, 1 );
    const float  gamma23 = invY( 2, 2 );
    const float  gamma31 = invY( 3, 0 );
    const float  gamma32 = invY( 3, 1 );
    const float  gamma33 = invY( 3, 2 );


    m_XGradientBasis = AtlasAlphasType( m_AlphasInVertex0.Size() );
    m_YGradientBasis = AtlasAlphasType( m_AlphasInVertex0.Size() );
    m_ZGradientBasis = AtlasAlphasType( m_AlphasInVertex0.Size() );
    for ( int labelNumber = 0; labelNumber < m_AlphasInVertex0.Size(); labelNumber++ )
    {
      //std::cout << "labelNumber: " << labelNumber << std::endl;

      const float  alpha0 = m_AlphasInVertex0[ labelNumber ];
      const float  alpha1 = m_AlphasInVertex1[ labelNumber ];
      const float  alpha2 = m_AlphasInVertex2[ labelNumber ];
      const float  alpha3 = m_AlphasInVertex3[ labelNumber ];

      m_XGradientBasis[ labelNumber ] = alpha1 * gamma11 + alpha2 * gamma21 + alpha3 * gamma31 - alpha0 * ( gamma11 + gamma21 + gamma31 );
      m_YGradientBasis[ labelNumber ] = alpha1 * gamma12 + alpha2 * gamma22 + alpha3 * gamma32 - alpha0 * ( gamma12 + gamma22 + gamma32 );
      m_ZGradientBasis[ labelNumber ] = alpha1 * gamma13 + alpha2 * gamma23 + alpha3 * gamma33 - alpha0 * ( gamma13 + gamma23 + gamma33 );
    }

    std::cout << "AtlasMeshVertexProcessor: m_XGradientBasis: " << m_XGradientBasis << std::endl;
    std::cout << "AtlasMeshVertexProcessor: m_YGradientBasis: " << m_YGradientBasis << std::endl;
    std::cout << "AtlasMeshVertexProcessor: m_ZGradientBasis: " << m_ZGradientBasis << std::endl;

    return true;
  }

  //
  float GetCost() const
  {
    return m_Cost;
  }

  //
  const AtlasPositionGradientType& GetGradient() const
  {
    return m_Gradient;
  }

  inline void SetMesh( const AtlasMesh* mesh ) {}

private:

  const unsigned char*  m_SourcePointer;

  AtlasAlphasType  m_AlphasInVertex0;
  AtlasAlphasType  m_AlphasInVertex1;
  AtlasAlphasType  m_AlphasInVertex2;
  AtlasAlphasType  m_AlphasInVertex3;

  AtlasAlphasType  m_XGradientBasis;
  AtlasAlphasType  m_YGradientBasis;
  AtlasAlphasType  m_ZGradientBasis;

  float  m_Cost;
  AtlasPositionGradientType  m_Gradient;

};


} // End namespace FragmentProcessor





class AtlasMeshVertexProcessor: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshVertexProcessor  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshVertexProcessor, itk::Object );

  // Some typedefs
  typedef itk::Image< unsigned char, 3 >  LabelImageType;

  //
  void SetMeshCollection( AtlasMeshCollection* meshCollection )
  {
    m_MeshCollection = meshCollection;
    m_VertexNeighborhoodsContainers.clear();
  }

  AtlasMeshCollection* GetMeshCollection()
  {
    return m_MeshCollection;
  }

  const AtlasMeshCollection* GetMeshCollection() const
  {
    return m_MeshCollection;
  }

  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages );

  // Get label images
  const std::vector< LabelImageType::ConstPointer >&  GetLabelImages() const
  {
    return m_LabelImages;
  }

  // Set beta
  void  SetBeta( float beta )
  {
    m_Beta = beta;
  }

  // Get beta
  float  GetBeta() const
  {
    return m_Beta;
  }

  //
  bool CalculateXstar( int meshNumber, AtlasMesh::PointIdentifier pointId,
                       float& xstar, float& ystar, float& zstar, bool verbose=false );

  //
  float  CalculateCost( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z );

  //
  AtlasPositionGradientType  CalculateGradient( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z );

  //
  Curvature  CalculateCurvature( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z, bool verbose=false );


protected :
  // Constructor
  AtlasMeshVertexProcessor();

  // Destructor
  virtual ~AtlasMeshVertexProcessor();

  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;

  //
  float  CalculatePriorCost( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z );

  //
  AtlasPositionGradientType  CalculatePriorGradient( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z );

#if 0
  //
  Curvature  CalculatePriorCurvature( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z );
#endif

  //
  struct VertexNeighboringTetrahedronInfo
  {
    AtlasMesh::CellIdentifier  m_TetrahedronId;

    const float*  m_ReferenceVolumeTimesK;

    // Elements of the inverse of the matrix [ p0 p1 p2 p3; 1 1 1 1 ]
    const float*  m_Z11;
    const float*  m_Z21;
    const float*  m_Z31;
    const float*  m_Z41;

    const float*  m_Z12;
    const float*  m_Z22;
    const float*  m_Z32;
    const float*  m_Z42;

    const float*  m_Z13;
    const float*  m_Z23;
    const float*  m_Z33;
    const float*  m_Z43;


    const AtlasMesh::PointType::ValueType*  m_X1;
    const AtlasMesh::PointType::ValueType*  m_Y1;
    const AtlasMesh::PointType::ValueType*  m_Z1;

    const AtlasMesh::PointType::ValueType*  m_X2;
    const AtlasMesh::PointType::ValueType*  m_Y2;
    const AtlasMesh::PointType::ValueType*  m_Z2;

    const AtlasMesh::PointType::ValueType*  m_X3;
    const AtlasMesh::PointType::ValueType*  m_Y3;
    const AtlasMesh::PointType::ValueType*  m_Z3;

    const AtlasAlphasType*  m_Alphas0;
    const AtlasAlphasType*  m_Alphas1;
    const AtlasAlphasType*  m_Alphas2;
    const AtlasAlphasType*  m_Alphas3;

    // Constructor
    VertexNeighboringTetrahedronInfo() : m_TetrahedronId( 0 ),
      m_ReferenceVolumeTimesK( 0 ),
      m_Z11( 0 ), m_Z21( 0 ), m_Z31( 0 ), m_Z41( 0 ),
      m_Z12( 0 ), m_Z22( 0 ), m_Z32( 0 ), m_Z42( 0 ),
      m_Z13( 0 ), m_Z23( 0 ), m_Z33( 0 ), m_Z43( 0 ),
      m_X1( 0 ), m_Y1( 0 ), m_Z1( 0 ),
      m_X2( 0 ), m_Y2( 0 ), m_Z2( 0 ),
      m_X3( 0 ), m_Y3( 0 ), m_Z3( 0 ),
      m_Alphas0( 0 ), m_Alphas1( 0 ), m_Alphas2( 0 ), m_Alphas3( 0 ) {}
  };

  typedef std::vector< VertexNeighboringTetrahedronInfo >  VertexNeighborhood;
  typedef itk::MapContainer< AtlasMesh::PointIdentifier, VertexNeighborhood >  VertexNeighborhoodsContainerType;

  //
  void CalculateVertexNeighborhoods();

  //
  const std::vector< VertexNeighborhoodsContainerType::Pointer >&  GetVertexNeighborhoodsContainers() const
  {
    return m_VertexNeighborhoodsContainers;
  }


private :
  AtlasMeshVertexProcessor(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Data members
  AtlasMeshCollection::Pointer  m_MeshCollection;
  std::vector< VertexNeighborhoodsContainerType::Pointer >  m_VertexNeighborhoodsContainers;
  float  m_Beta;
  std::vector< LabelImageType::ConstPointer >  m_LabelImages;


};



} // end namespace kvl


#endif
