/**
 * @file  Region2DRectangle.h
 * @brief Region2DRectangle data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/08/03 20:29:27 $
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

#ifndef Region2DRectangle_h
#define Region2DRectangle_h

#include "Region2D.h"
#include <wx/wx.h>
#include <string>
#include <vector>
#include "vtkSmartPointer.h"

class vtkActor2D;
class RenderView2D;
class vtkRenderer;

class Region2DRectangle : public Region2D
{
public:
  Region2DRectangle( RenderView2D* view );
  virtual ~Region2DRectangle();

  void UpdatePoint( int nIndex, int nX, int nY );
  
  void Offset( int nX, int nY );
    
  bool Contains( int nX, int nY, int* nIndexOut = NULL ); 

  void SetRect( int left, int top, int w, int h );
  
  void SetTopLeft( int left, int top );
  void SetBottomRight( int right, int bottom );
  
  void AppendProp( vtkRenderer* renderer );
  
  void Show( bool bshow );
  
  void Update();
  
  void GetWorldPoint1( double* pt );
  void GetWorldPoint2( double* pt );

  double    m_dPt0[3];       // rect in screen coordinate
  double    m_dPt1[3];
  double    m_dPt2[3];
  double    m_dPt3[3];

protected:
  void UpdateWorldCoords();  
  
  vtkSmartPointer<vtkActor2D> m_actorRect;
  int       m_nX1, m_nX2, m_nY1, m_nY2;
};

#endif


