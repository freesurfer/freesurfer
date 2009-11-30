/**
 * @file  Region2D.h
 * @brief Region2D data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/11/30 21:17:20 $
 *    $Revision: 1.3 $
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

#ifndef Region2D_h
#define Region2D_h

#include <wx/wx.h>
#include "Broadcaster.h"
#include "Listener.h"

class RenderView2D;
class vtkRenderer;

class Region2D : public Broadcaster, public Listener
{
public:
  Region2D( RenderView2D* view );
  virtual ~Region2D();

  virtual void Offset( int nX, int nY ) = 0;
  
  virtual void UpdatePoint( int nIndex, int nX, int nY ) = 0;
  
  virtual bool Contains( int nX, int nY, int* nPointIndex = NULL ) = 0; 
  
  virtual void Highlight( bool bHighlight = true ) {}
  
  virtual void AppendProp( vtkRenderer* renderer ) = 0;
  
  virtual void Show( bool bshow );
  
  virtual void Update();
  
  virtual void UpdateSlicePosition( int nPlane, double pos );
  
  virtual void GetWorldPoint( int nIndex, double* pt ) = 0;

protected:
  RenderView2D* m_view;
};

#endif


