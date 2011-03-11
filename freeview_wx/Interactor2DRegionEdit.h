/**
 * @file  Interactor2DRegionEdit.h
 * @brief Interactor2DRegionEdit to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:38 $
 *    $Revision: 1.1 $
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

#ifndef Interactor2DRegionEdit_h
#define Interactor2DRegionEdit_h

#include "Interactor2D.h"
#include <vector>

class wxWindow;

class Interactor2DRegionEdit : public Interactor2D
{
public:
  Interactor2DRegionEdit();
  virtual ~Interactor2DRegionEdit();

  enum EditMode { EM_Freehand = 0, EM_Fill, EM_Polyline };

  // return true if to have parent Interactor2DROIEdit continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseUpEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessKeyDownEvent( wxKeyEvent& event, RenderView* view );

protected:
  void UpdateCursor( wxEvent& event, wxWindow* wnd );

  bool m_bEditing;

  std::vector<double>  m_dPolylinePoints;
};

#endif


