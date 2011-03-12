/**
 * @file  Region2DRectangle.h
 * @brief Region2DRectangle data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:53 $
 *    $Revision: 1.7 $
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
#include <vector>
#include "vtkSmartPointer.h"

class vtkTextActor;
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
  
  void Highlight( bool bHighlight );
  
  void Update();
  
  void UpdateStats();
  
  void UpdateSlicePosition( int nPlane, double pos );
  
  void GetWorldPoint( int nIndex, double* pt );
  
  void SetEnableStats( bool bStats )
  {
    m_bEnableStats = bStats;
  }

protected:
  void UpdateWorldCoords();  
  int GetRange( double[3][2] );
  
  vtkSmartPointer<vtkActor2D>   m_actorRect;
  vtkSmartPointer<vtkTextActor> m_actorText;
  int       m_nX1, m_nX2, m_nY1, m_nY2;
  double    m_dPt[4][3];       // rect in world coordinate
  
  bool      m_bEnableStats;
};

#endif


