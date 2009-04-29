/**
 * @file  Cursor2D.h
 * @brief Annotation class for 2D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/04/29 22:53:48 $
 *    $Revision: 1.3.2.3 $
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

#ifndef Cursor2D_h
#define Cursor2D_h

#include "RenderView.h"
#include "vtkSmartPointer.h"
#include <wx/colour.h>
#include <vector>

class vtkActor2D;
class vtkRenderer;
class vtkActor;
class RenderView2D;

class Cursor2D
{
public:
  Cursor2D( RenderView2D* view );
  virtual ~Cursor2D();

  void SetPosition( double* pos, bool bConnectPrevious = false );
  void SetPosition2( double* pos);

  double* GetPosition();
  void GetPosition( double* pos );

  void SetInterpolationPoints( std::vector<double> pts );

  std::vector<double> GetInterpolationPoints()
  {
    return m_dInterpolationPoints;
  }

  void ClearInterpolationPoints()
  {
    m_dInterpolationPoints.clear();
  }

  void GetColor( double* rgb );
  void SetColor( double r, double g, double b );

  wxColour GetColor();
  void SetColor( const wxColour& color );

  int GetRadius();
  void SetRadius( int nPixels );

  void Update( bool bConnectPrevious = false );

  void AppendActor( vtkRenderer* renderer );

  void Show( bool bShow = true );

  bool IsShown();

private:
  vtkSmartPointer<vtkActor2D> m_actorCursor;

  RenderView2D* m_view;

  double  m_dPosition[3];
  double  m_dPosition2[3];
  std::vector<double> m_dInterpolationPoints;

  int   m_nRadius;
};

#endif


