/**
 * @file  Interactor2DMeasure.h
 * @brief Interactor for measure tool in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.5 $
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

#ifndef Interactor2DMeasure_h
#define Interactor2DMeasure_h

#include "Interactor2D.h"
#include <wx/wx.h>
#include <vector>
#include <string>

class wxWindow;
class Region2D;

class Interactor2DMeasure : public Interactor2D
{
public:
  Interactor2DMeasure();
  virtual ~Interactor2DMeasure();
  
  // return true if to have parent Interactor2D continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseUpEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessKeyDownEvent( wxKeyEvent& event, RenderView* view );
  virtual bool ProcessKeyUpEvent( wxKeyEvent& event, RenderView* view );

protected:
  void UpdateCursor( wxEvent& event, wxWindow* wnd );
  void ProcessContextMenu( wxMouseEvent& event );

  bool        m_bEditing;
  bool        m_bDrawing;
  int         m_nPointIndex;
  Region2D*   m_region;
};

#endif


