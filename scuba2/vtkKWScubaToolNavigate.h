/**
 * @file  vtkKWScubaToolNavigate.h
 * @brief Navigate tool
 *
 * Implements panning, rotating, 2DRASZ shifts, and zooming. Works in
 * 2D and 3D.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:06 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
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


#ifndef vtkKWScubaToolNavigate_h
#define vtkKWScubaToolNavigate_h

#include "vtkKWScubaTool.h"

class vtkKWScubaToolNavigate : public vtkKWScubaTool {

public:

  static vtkKWScubaToolNavigate* New ();
  vtkTypeRevisionMacro( vtkKWScubaToolNavigate, vtkKWScubaTool );

protected:

  vtkKWScubaToolNavigate ();
  ~vtkKWScubaToolNavigate ();

  // Description:
  // We don't want pick events on mouse dags so that our navigation
  // will be faster.
  virtual bool SuspendPickEvents ();

  // Description:
  // Middle mouse button scrolls. When it goes down or up, set the
  // layers' Fast Mode on or off accordingly.
  virtual void DoMouseDown ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                             vtkKWScubaLayer* iLayer, float iRAS[3] );
  virtual void DoMouseUp   ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                             vtkKWScubaLayer* iLayer, float iRAS[3] );

  // Description:
  // Tell the view to pan, scroll, or zoom.
  virtual void DoMouseDrag  ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );

};

#endif
