/**
 * @file  SurfaceRegion.h
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/04/30 21:21:19 $
 *    $Revision: 1.1 $
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
#include <wx/colour.h>

class vtkRenderer;
class vtkActor;
class vtkPolyData;
class vtkPoints;
class vtkSelectPolyData;
class RenderView3D;

class SurfaceRegion
{
public:
  SurfaceRegion( );
  virtual ~SurfaceRegion();

  void SetInput( vtkPolyData* polydata );

  void AddPoint( double* pt );
  
  wxColour GetColor();
  void SetColor( const wxColour& color );

  void Update();

  void AppendActor( vtkRenderer* renderer );

  void Show( bool bShow = true );

private:
  void RebuildActor();

  vtkSmartPointer<vtkActor> m_actorMesh;
  vtkSmartPointer<vtkPoints>  m_points;
  
  vtkSmartPointer<vtkSelectPolyData>  m_selector;
};

#endif


