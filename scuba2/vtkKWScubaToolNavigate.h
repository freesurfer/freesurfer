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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.2 $
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
