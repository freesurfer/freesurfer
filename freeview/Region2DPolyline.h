/**
 * @file  Region2DPolyline.h
 * @brief Region2DPolyline data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:58 $
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

#ifndef Region2DPolyline_h
#define Region2DPolyline_h

#include "Region2D.h"
#include <vector>
#include "vtkSmartPointer.h"

class vtkActor2D;
class RenderView2D;
class vtkRenderer;
class vtkTextActor;

class Region2DPolyline : public Region2D
{
public:
  Region2DPolyline( RenderView2D* view, bool bSpline = false );
  virtual ~Region2DPolyline();

  void AddPoint( int x, int y );
  
  void Offset( int x, int y );
  
  bool Contains( int x, int y, int* nIndexOut = NULL );
  void UpdatePoint( int nIndex, int nX, int nY );
  
  void AppendProp( vtkRenderer* renderer );
  
  void Show( bool bshow = true );
  void Highlight( bool bHighlight = true );
  
  void Update();
  void UpdateStats();
  
  void UpdateSlicePosition( int nPlane, double pos );
  
  void GetWorldPoint( int nIndex, double* pt );
  
  void RemoveLastPoint();

protected:
  void UpdateWorldCoords();  
  
  vtkSmartPointer<vtkActor2D>   m_actorPolyline;
  vtkSmartPointer<vtkActor2D>   m_actorPoints;
  vtkSmartPointer<vtkTextActor> m_actorText;
  
  struct ScreenPoint {
    int pos[2];
  };
  
  struct WorldPoint {
    double pos[3];
  };
  
  std::vector<ScreenPoint>   m_screenPts;       // 2D points
  std::vector<WorldPoint>    m_worldPts;
  
  bool    m_bSpline;
};

#endif


