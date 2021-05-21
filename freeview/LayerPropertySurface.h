/**
 * @brief The common properties available to MRI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 *
 */

#ifndef LayerPropertySurface_h
#define LayerPropertySurface_h

#include "vtkSmartPointer.h"
#include "LayerProperty.h"
#include <QColor>
#include <QVariantMap>

class vtkLookupTable;
class vtkRGBAColorTransferFunction;
class FSSurface;

class LayerPropertySurface  : public LayerProperty
{
  Q_OBJECT
public:
  LayerPropertySurface ( QObject* parent = NULL );
  ~LayerPropertySurface ();

  enum CURVATURE_MAP
  { CM_Off = 0, CM_Threshold, CM_Binary };

  enum SuraceRenderMode
  { SM_Surface = 0, SM_Wireframe, SM_SurfaceAndWireframe };

  enum MeshColorMap
  { MC_Surface = 0, MC_Curvature, MC_Overlay, MC_Solid };

  double GetOpacity() const;

  double* GetBinaryColor()
  {
    return m_dRGB;
  }

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

  double GetThresholdSlope()
  {
    return m_dThresholdSlope;
  }

  double* GetEdgeColor()
  {
    return m_dRGBEdge;
  }

  double* GetVectorColor()
  {
    return m_dRGBVector;
  }

  int GetEdgeThickness()
  {
    return m_nEdgeThickness;
  }

  int GetVectorPointSize()
  {
    return m_nVectorPointSize;
  }

  void SetSurfaceSource( FSSurface* surf );

  vtkRGBAColorTransferFunction* GetCurvatureLUT() const;

  void BuildCurvatureLUT( vtkRGBAColorTransferFunction* lut, int nMap );

  void RebuildCurvatureLUT();

  int GetCurvatureMap()
  {
    return m_nCurvatureMap;
  }

  int GetSurfaceRenderMode()
  {
    return m_nSurfaceRenderMode;
  }

  double* GetVertexColor()
  {
    return m_dRGBVertex;
  }

  bool GetShowVertices()
  {
    return m_bShowVertices;
  }

  int GetVertexPointSize()
  {
    return m_nVertexPointSize;
  }

  double* GetMeshColor()
  {
    return m_dRGBMesh;
  }

  int GetMeshColorMap()
  {
    return m_nMeshColorMap;
  }

  double* GetPosition()
  {
    return m_dPosition;
  }
  void SetPosition( double* p );

  QVariantMap GetFullSettings();

  void RestoreFullSettings(const QVariantMap& map);

  bool GetShowOverlay()
  {
    return m_bShowOverlay;
  }

  bool GetShowAnnotation()
  {
    return m_bShowAnnotation;
  }

  bool GetUseSurfaceColorOn2D()
  {
      return m_bUseSurfaceColorOn2D;
  }

  int GetZOrderOverlay()
  {
    return m_nZOrderOverlay;
  }

  int GetZOrderLabel()
  {
    return m_nZOrderLabel;
  }

  int GetZOrderAnnotation()
  {
    return m_nZOrderAnnotation;
  }

public slots:
  void SetOpacity( double opacity );
  void SetCurvatureMap( int nMap );
  void SetSurfaceRenderMode( int nMode );
  void SetBinaryColor( double r, double g, double b );
  void SetBinaryColor( const QColor& c );
  void SetThresholdMidPoint( double dvalue );
  void SetThresholdSlope( double dvalue );
  void SetEdgeColor( double r, double g, double b );
  void SetEdgeColor( const QColor& c );
  void SetVectorColor( double r, double g, double b );
  void SetVectorColor( const QColor& c );
  void SetVertexColor( double r, double g, double b );
  void SetVertexColor( const QColor& c );
  void SetMeshColor( double r, double g, double b );
  void SetMeshColor( const QColor& c );
  void SetVertexPointSize( int nSize );
  void SetMeshColorMap( int nMap );
  void ShowVertices( bool bShow );
  void SetEdgeThickness( int nThickness );
  void SetVectorPointSize( int nSize );
  void SetShowOverlay(bool bShow);
  void SetShowAnnotation(bool bShow);
  void SetUseSurfaceColorOn2D(bool bKeep);
  void SetZOrderOverlay(int nOrder);
  void SetZOrderAnnotation(int nOrder);
  void SetZOrderLabel(int nOrder);

Q_SIGNALS:
  void OpacityChanged( double opacity );
  void EdgeThicknessChanged( int nThickness );
  void VectorPointSizeChanged( int nSize );
  void RenderModeChanged( int nMode );
  void VertexRenderChanged();
  void MeshRenderChanged();
  void ColorMapChanged();
  void PositionChanged();
  void PositionChanged(double dx, double dy, double dz);
  void OverlayChanged();
  void EdgeColorChanged();

private:
  void SetColorMapChanged();

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

  bool    m_bShowOverlay;
  bool    m_bShowAnnotation;

  bool    m_bUseSurfaceColorOn2D;

  int     m_nZOrderOverlay;
  int     m_nZOrderLabel;
  int     m_nZOrderAnnotation;

  vtkSmartPointer<vtkRGBAColorTransferFunction> m_lutCurvature;

  FSSurface* m_surface;
};

#endif
