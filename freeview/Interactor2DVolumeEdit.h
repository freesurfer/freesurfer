/**
 * @file  Interactor2DVolumeEdit.h
 * @brief Interactor for editing volume in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/07/07 22:05:04 $
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

#ifndef Interactor2DVolumeEdit_h
#define Interactor2DVolumeEdit_h

#include "Interactor2D.h"
#include <vector>
#include <string>

class wxWindow;

class Interactor2DVolumeEdit : public Interactor2D
{
public:
  Interactor2DVolumeEdit( const char* layerTypeName );
  virtual ~Interactor2DVolumeEdit();

  enum EditMode
  {
    EM_Freehand = 0, EM_Fill, EM_Polyline, EM_Livewire, EM_ColorPicker
  };

  // return true if to have parent Interactor2D continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseUpEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessKeyDownEvent( wxKeyEvent& event, RenderView* view );
  virtual bool ProcessKeyUpEvent( wxKeyEvent& event, RenderView* view );

protected:
  void UpdateCursor( wxEvent& event, wxWindow* wnd );

  bool m_bEditing;

  std::string m_strLayerTypeName;

  std::vector<double>  m_dPolylinePoints;
};

#endif


