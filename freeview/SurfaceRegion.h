/**
 * @file  SurfaceRegion.h
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 22:00:37 $
 *    $Revision: 1.14 $
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
class vtkCleanPolyData;
class RenderView3D;
class LayerMRI;

class SurfaceRegion : public Broadcaster, public Listener
{
public:
  SurfaceRegion( LayerMRI* owner );
  virtual ~SurfaceRegion();
  
  void SetInput( vtkPolyData* polydata );

  void AddPoint( double* pt );
  
  bool Close();
  
  void ResetOutline();
  
  wxColour GetColor();
  void SetColor( const wxColour& color );

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
  
  bool Write( wxString& fn );
  
  static bool WriteHeader( FILE* fp, LayerMRI* mri_ref, int nNum = 1 );
  
  bool WriteBody( FILE* fp );
  
  bool Load( FILE* fp );
  
  LayerMRI* GetMRI()
  {
    return m_mri;
  }
  
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
  
  LayerMRI*     m_mri;
  wxColour      m_color;
  int   m_nId;
  int   m_nGroup;
};

#endif


