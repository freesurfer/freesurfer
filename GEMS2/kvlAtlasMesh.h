#ifndef __kvlAtlasMesh_h
#define __kvlAtlasMesh_h


#include "itkMesh.h"

#ifndef USE_DYNAMIC_MESH
  #include "itkDefaultStaticMeshTraits.h"
#else
  #include "itkDefaultDynamicMeshTraits.h"
#endif  
#include "itkVector.h"


namespace kvl
{

struct ReferenceTetrahedronInfo
{
  double  m_ReferenceVolumeTimesK;

  // Elements of the inverse of the matrix [ p0 p1 p2 p3; 1 1 1 1 ]
  double  m_Z11;
  double  m_Z21;
  double  m_Z31;
  double  m_Z41;

  double  m_Z12;
  double  m_Z22;
  double  m_Z32;
  double  m_Z42;

  double  m_Z13;
  double  m_Z23;
  double  m_Z33;
  double  m_Z43;

};


typedef itk::Array< float >  AtlasAlphasType;

struct PointParameters
{
  AtlasAlphasType  m_Alphas;
  bool  m_CanChangeAlphas;
  bool  m_CanMoveX;
  bool  m_CanMoveY;
  bool  m_CanMoveZ;

};


// Some typedefs
#ifndef USE_DYNAMIC_MESH
  typedef itk::DefaultStaticMeshTraits< PointParameters, 3, 3, 
                                        double, double, ReferenceTetrahedronInfo >  AtlasMeshTraits;
#else
  typedef itk::DefaultDynamicMeshTraits< PointParameters, 3, 3, 
                                        double, double, ReferenceTetrahedronInfo >  AtlasMeshTraits;
#endif

typedef itk::Mesh< PointParameters, 3, AtlasMeshTraits >  AtlasMesh;

typedef itk::Vector< double, 3 >  AtlasPositionGradientType;
#ifndef USE_DYNAMIC_MESH
  typedef itk::VectorContainer< AtlasMesh::PointIdentifier, AtlasPositionGradientType > 
                                                                  AtlasPositionGradientContainerType;
#else
  typedef itk::MapContainer< AtlasMesh::PointIdentifier, AtlasPositionGradientType > 
                                                                  AtlasPositionGradientContainerType;
#endif
struct Curvature
  {
  float  m_Curvature_dxdx;
  float  m_Curvature_dydy;
  float  m_Curvature_dzdz;
  float  m_Curvature_dxdy;
  float  m_Curvature_dxdz;
  float  m_Curvature_dydz;

  Curvature() : m_Curvature_dxdx( 0 ), m_Curvature_dydy( 0 ), m_Curvature_dzdz( 0 ),
                m_Curvature_dxdy( 0 ), m_Curvature_dxdz( 0 ), m_Curvature_dydz( 0 ) {}
  ~Curvature() {}

  const Curvature& operator+=(const Curvature& curvature )
    {
    m_Curvature_dxdx += curvature.m_Curvature_dxdx;
    m_Curvature_dydy += curvature.m_Curvature_dydy;
    m_Curvature_dzdz += curvature.m_Curvature_dzdz;
    m_Curvature_dxdy += curvature.m_Curvature_dxdy;
    m_Curvature_dxdz += curvature.m_Curvature_dxdz;
    m_Curvature_dydz += curvature.m_Curvature_dydz;

    return *this;
    }

  const Curvature& operator-=(const Curvature& curvature )
    {
    m_Curvature_dxdx -= curvature.m_Curvature_dxdx;
    m_Curvature_dydy -= curvature.m_Curvature_dydy;
    m_Curvature_dzdz -= curvature.m_Curvature_dzdz;
    m_Curvature_dxdy -= curvature.m_Curvature_dxdy;
    m_Curvature_dxdz -= curvature.m_Curvature_dxdz;
    m_Curvature_dydz -= curvature.m_Curvature_dydz;

    return *this;
    }

  };

#ifndef USE_DYNAMIC_MESH
  typedef itk::VectorContainer< AtlasMesh::PointIdentifier, Curvature >   AtlasPositionCurvatureContainerType;
#else  
  typedef itk::MapContainer< AtlasMesh::PointIdentifier, Curvature >   AtlasPositionCurvatureContainerType;
#endif
  
} // end namespace kvl


#endif
