/**
 * @file  Interactor.h
 * @brief Base Interactor class manage mouse and key input in render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.9 $
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
    MM_Line = 0, MM_Polyline, MM_Spline, MM_Rectangle, MM_SurfaceRegion, MM_Label
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


