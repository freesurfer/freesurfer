/**
 * @file  SurfaceRegion.h
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/05/25 18:27:34 $
 *    $Revision: 1.6 $
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
#include "Broadcaster.h"
#include "Listener.h"
#include <wx/colour.h>

class vtkRenderer;
class vtkActor;
class vtkPolyData;
class vtkPoints;
class vtkSelectPolyData;
class vtkBox;
class vtkProp;
class vtkClipPolyData;
class RenderView3D;
class LayerMRI;

class SurfaceRegion : public Broadcaster, public Listener
{
public:
  SurfaceRegion( LayerMRI* owner );
  virtual ~SurfaceRegion();

  void SetInput( vtkPolyData* polydata );

  void AddPoint( double* pt );
  
  void Close();
  
  wxColour GetColor();
  void SetColor( const wxColour& color );

  void Update();

  void AppendActor( vtkRenderer* renderer );

  void Show( bool bShow = true );
  
  bool HasPoint( double* pos );
  
  void Highlight( bool bHighlight = true );
  
  vtkActor* GetMeshActor();
  
  int GetId()
  {
    return m_nId;
  }
  
  void SetId( int nId )
  {
    m_nId = nId;
  }

  bool Write( wxString& fn );
  
  static bool WriteHeader( FILE* fp, LayerMRI* mri_ref, int nNum = 1 );
  
  bool WriteBody( FILE* fp );
  
  bool Load( FILE* fp );
  
  
private:
  void RebuildOutline( bool bClose );

  vtkSmartPointer<vtkActor>   m_actorMesh;
  vtkSmartPointer<vtkActor>   m_actorOutline;
  vtkSmartPointer<vtkPoints>  m_points;
  
  vtkSmartPointer<vtkBox>     m_clipbox;
  vtkSmartPointer<vtkClipPolyData>    m_clipper;
  vtkSmartPointer<vtkSelectPolyData>  m_selector;
  
  LayerMRI*     m_mri;
  
  int   m_nId;
};

#endif


