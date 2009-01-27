/**
 * @file  Cursor3D.h
 * @brief Annotation class for 2D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/01/27 18:27:24 $
 *    $Revision: 1.5 $
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

#ifndef Cursor3D_h
#define Cursor3D_h

#include "RenderView.h"
#include "vtkSmartPointer.h"
#include <wx/colour.h>

class vtkRenderer;
class vtkActor;
class RenderView3D;

class Cursor3D
{
public:
  Cursor3D( RenderView3D* view );
  virtual ~Cursor3D();

  void SetPosition( double* pos );

  double* GetPosition();
  void GetPosition( double* pos );

  void GetColor( double* rgb );
  void SetColor( double r, double g, double b );

  wxColour GetColor();
  void SetColor( const wxColour& color );

  void Update();

  void AppendActor( vtkRenderer* renderer );

  void Show( bool bShow = true );

  bool IsShown();

private:
  void RebuildActor();

  vtkSmartPointer<vtkActor> m_actorCursor;

  RenderView3D* m_view;

  double  m_dPosition[3];
};

#endif


