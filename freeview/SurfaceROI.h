/**
 * @brief Surface region from a surface selection in 3D view.
 *
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

#ifndef SurfaceROI_h
#define SurfaceROI_h

#include "vtkSmartPointer.h"
#include <QObject>
#include <QColor>

class vtkRenderer;
class vtkActor;
class vtkPolyData;
class vtkPoints;
class vtkSelectPolyData;
class vtkBox;
class vtkProp;
class vtkClipPolyData;
class vtkCleanPolyData;
class RenderView3D;
class LayerSurface;

class SurfaceROI : public QObject
{
  Q_OBJECT
public:
  SurfaceROI( LayerSurface* owner );
  virtual ~SurfaceROI();

  void SetInput( vtkPolyData* polydata );

  void AddPoint( double* pt );

  void Close();

  void InitializeOutline(double * pos);

  QColor GetColor();
  void SetColor( const QColor& color );

  void Update();

  void AppendProps( vtkRenderer* renderer );

  void Show( bool bShow = true );

  int GetId()
  {
    return m_nId;
  }

  void SetId( int nId )
  {
    m_nId = nId;
  }

  vtkActor* GetActor();

  LayerSurface* GetSurface()
  {
    return m_mris;
  }

signals:
  void ColorChanged( const QColor& );

private:
  void RebuildOutline( bool bClose );

  vtkSmartPointer<vtkActor>   m_actorOutline;
  vtkSmartPointer<vtkPoints>  m_points;

  LayerSurface*   m_mris;
  QColor      m_color;
  int   m_nId;
};

#endif


