/**
 * @file  Region2DLine.h
 * @brief Region2DLine data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:58 $
 *    $Revision: 1.9 $
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

#ifndef Region2DLine_h
#define Region2DLine_h

#include "Region2D.h"
#include "vtkSmartPointer.h"

class vtkActor2D;
class RenderView2D;
class vtkRenderer;
class vtkTextActor;

class Region2DLine : public Region2D
{
public:
  Region2DLine( RenderView2D* view );
  virtual ~Region2DLine();

  void SetLine( int x1, int y1, int x2, int y2 );
  
  void SetPoint1( int x1, int y1 );
  void SetPoint2( int x2, int y2 );
  
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

protected:
  void UpdateWorldCoords();  
  
  vtkSmartPointer<vtkActor2D>   m_actorLine;
  vtkSmartPointer<vtkTextActor> m_actorText;
  int       m_nX1, m_nX2, m_nY1, m_nY2;  
  double    m_dPt1[3];      
  double    m_dPt2[3];

};

#endif


