/**
 * @file  LayerPropertiesSurface.h
 * @brief The common properties available to MRI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:39 $
 *    $Revision: 1.1 $
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

#ifndef LayerPropertiesSurface_h
#define LayerPropertiesSurface_h

#include "vtkSmartPointer.h"
#include "LayerProperties.h"

class vtkLookupTable;
class vtkRGBAColorTransferFunction;
class FSSurface;

class LayerPropertiesSurface  : public LayerProperties
{
public:
  LayerPropertiesSurface ();
  ~LayerPropertiesSurface ();

  enum CURVATURE_MAP 
    { CM_Off = 0, CM_Threshold, CM_Binary };
  
  enum SuraceRenderMode 
    { SM_Surface = 0, SM_Wireframe, SM_SurfaceAndWireframe };
  
  enum MeshColorMap
    { MC_Surface = 0, MC_Curvature, MC_Overlay, MC_Solid };

  double GetOpacity() const;
  void SetOpacity( double opacity );

  double* GetBinaryColor()
  {
    return m_dRGB;
  }
  void SetBinaryColor( double r, double g, double b );

  double* GetThresholdHighColor()
  {
    return m_dRGBThresholdHigh;
  }
  double* GetThresholdLowColor()
  {
    return m_dRGBThresholdLow;
  }
  void SetThresholdColor( double* low, double* high );

  double GetThresholdMidPoint()
  {
    return m_dThresholdMidPoint;
  }
  void SetThresholdMidPoint( double dvalue );

  double GetThresholdSlope()
  {
    return m_dThresholdSlope;
  }
  void SetThresholdSlope( double dvalue );

  double* GetEdgeColor()
  {
    return m_dRGBEdge;
  }
  void SetEdgeColor( double r, double g, double b );

  double* GetVectorColor()
  {
    return m_dRGBVector;
  }
  void SetVectorColor( double r, double g, double b );

  int GetEdgeThickness()
  {
    return m_nEdgeThickness;
  }
  void SetEdgeThickness( int nThickness );

  int GetVectorPointSize()
  {
    return m_nVectorPointSize;
  }
  void SetVectorPointSize( int nSize );

  void SetSurfaceSource( FSSurface* surf );

  vtkRGBAColorTransferFunction* GetCurvatureLUT() const;
  
  void BuildCurvatureLUT( vtkRGBAColorTransferFunction* lut, int nMap );

  int GetCurvatureMap()
  {
    return m_nCurvatureMap;
  }
  void SetCurvatureMap( int nMap );
  
  int GetSurfaceRenderMode()
  {
    return m_nSurfaceRenderMode;
  }
  void SetSurfaceRenderMode( int nMode );
  
  double* GetVertexColor()
  {
    return m_dRGBVertex;
  }
  void SetVertexColor( double r, double g, double b );
  
  bool GetShowVertices()
  {
    return m_bShowVertices;
  }
  void ShowVertices( bool bShow );
  
  int GetVertexPointSize()
  {
    return m_nVertexPointSize;
  }
  void SetVertexPointSize( int nSize );
  
  double* GetMeshColor()
  {
    return m_dRGBMesh;
  }
  void SetMeshColor( double r, double g, double b );
  
  int GetMeshColorMap()
  {
    return m_nMeshColorMap;
  }
  void SetMeshColorMap( int nMap );
  
  double* GetPosition()
  {
    return m_dPosition;
  }
  void SetPosition( double* p );

private:
  void ColorMapChanged ();

  double  m_dOpacity;
  double  m_dRGB[3];
  double  m_dRGBThresholdHigh[3];
  double  m_dRGBThresholdLow[3];
  double  m_dRGBEdge[3];
  double  m_dRGBVector[3];
  double  m_dRGBMesh[3];
  int     m_nEdgeThickness;
  int     m_nVectorPointSize;

  double  m_dThresholdMidPoint;
  double  m_dThresholdSlope;
  
  double  m_dPosition[3];
 
  int     m_nCurvatureMap;
  
  int     m_nSurfaceRenderMode;
  
  bool    m_bShowVertices;
  double  m_dRGBVertex[3];
  int     m_nVertexPointSize;

  int     m_nMeshColorMap;
  
  vtkSmartPointer<vtkRGBAColorTransferFunction> m_lutCurvature;
  
  FSSurface* m_surface;
};

#endif
