/**
 * @file  kvlAtlasMesh.h
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
#ifndef __kvlAtlasMesh_h
#define __kvlAtlasMesh_h

#include "itkMesh.h"
#include "itkDefaultDynamicMeshTraits.h"
#include "itkVector.h"


namespace kvl
{

struct ReferenceTetrahedronInfo
{
  float  m_ReferenceVolumeTimesK;

  // Elements of the inverse of the matrix [ p0 p1 p2 p3; 1 1 1 1 ]
  float  m_Z11;
  float  m_Z21;
  float  m_Z31;
  float  m_Z41;

  float  m_Z12;
  float  m_Z22;
  float  m_Z32;
  float  m_Z42;

  float  m_Z13;
  float  m_Z23;
  float  m_Z33;
  float  m_Z43;

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
typedef itk::DefaultDynamicMeshTraits< PointParameters, 3, 3,
        float, float, ReferenceTetrahedronInfo >  AtlasMeshTraits;


typedef itk::Mesh< PointParameters, 3, AtlasMeshTraits >  AtlasMesh;

typedef itk::Vector< float, 3 >  AtlasPositionGradientType;
typedef itk::MapContainer< AtlasMesh::PointIdentifier, AtlasPositionGradientType >
AtlasPositionGradientContainerType;
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

typedef itk::MapContainer< AtlasMesh::PointIdentifier, Curvature >   AtlasPositionCurvatureContainerType;

} // end namespace kvl


#endif
