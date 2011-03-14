/**
 * @file  SurfaceRegion.h
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:59 $
 *    $Revision: 1.17 $
 *
 * Copyright (C) 2008-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef SurfaceRegion_h
#define SurfaceRegion_h

#include "RenderView.h"
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
class LayerMRI;

class SurfaceRegion : public QObject
{
    Q_OBJECT
public:
  SurfaceRegion( LayerMRI* owner );
  virtual ~SurfaceRegion();
  
  void SetInput( vtkPolyData* polydata );

  void AddPoint( double* pt );
  
  bool Close();
  
  void ResetOutline();
  
  QColor GetColor();
  void SetColor( const QColor& color );

  void Update();

  void AppendProps( vtkRenderer* renderer );

  void Show( bool bShow = true );
  
  bool HasPoint( double* pos );
  
  void Highlight( bool bHighlight = true );
  
  bool DeleteCell( RenderView3D* view, int pos_x, int pos_y );
  
  vtkActor* GetMeshActor();
  
  int GetId()
  {
    return m_nId;
  }
  
  void SetId( int nId )
  {
    m_nId = nId;
  }

  int GetGroup()
  {
    return m_nGroup;
  }

  void SetGroup( int n );

  bool Write( const QString& fn );
  
  static bool WriteHeader( FILE* fp, LayerMRI* mri_ref, int nNum = 1 );
  
  bool WriteBody( FILE* fp );
  
  bool Load( FILE* fp );

  LayerMRI* GetMRI()
  {
    return m_mri;
  }

signals:
  void ColorChanged( const QColor& );
  
private:
  void RebuildOutline( bool bClose );

  vtkSmartPointer<vtkActor>   m_actorMesh;
  vtkSmartPointer<vtkActor>   m_actorOutline;
  vtkSmartPointer<vtkPoints>  m_points;
  
  vtkSmartPointer<vtkBox>     m_clipbox;
  vtkSmartPointer<vtkClipPolyData>    m_clipperPre;
  vtkSmartPointer<vtkSelectPolyData>  m_selector;
  vtkSmartPointer<vtkCleanPolyData>   m_cleanerPost;
  
  vtkSmartPointer<vtkPolyData>        m_polydataHolder;
  
  LayerMRI*   m_mri;
  QColor      m_color;
  int   m_nId;
  int   m_nGroup;
};

#endif


