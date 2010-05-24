/**
 * @file  Interactor.h
 * @brief Base Interactor class manage mouse and key input in render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/05/24 21:42:53 $
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

#ifndef Interactor_h
#define Interactor_h

#include <wx/wx.h>
#include "Broadcaster.h"

class vtkRenderer;
class RenderView;

class Interactor : public Broadcaster
{
public:
  Interactor();
  virtual ~Interactor();
  
  enum MeasureMode 
  { 
    MM_Line = 0, MM_Polyline, MM_Spline, MM_Rectangle, MM_SurfaceRegion
  };

  enum EditMode
  {
    EM_Freehand = 0, EM_Fill, EM_Polyline, EM_Livewire, EM_ColorPicker, EM_Contour
  };
  
  int GetAction();
  void SetAction( int nAction );

  // return true if to have parent interactor continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseUpEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseWheelEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseEnterEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseLeaveEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessKeyDownEvent( wxKeyEvent& event, RenderView* view );
  virtual bool ProcessKeyUpEvent( wxKeyEvent& event, RenderView* view );
  virtual void ProcessPostMouseMoveEvent( wxMouseEvent& event, RenderView* view )
  {}
  virtual void ProcessPostMouseWheelEvent( wxMouseEvent& event, RenderView* view )
  {}

protected:
  virtual void UpdateCursor( wxEvent& event, wxWindow* wnd );

  int  m_nDownPosX;
  int  m_nDownPosY;

  int  m_nAction;
};

#endif


