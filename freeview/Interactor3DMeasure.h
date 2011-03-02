/**
 * @file  Interactor3DMeasure.h
 * @brief Interactor for measurement in 3D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
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

#ifndef Interactor3DMeasure_h
#define Interactor3DMeasure_h

#include "Interactor3D.h"

class Interactor3DMeasure : public Interactor3D
{
public:
  Interactor3DMeasure();
  ~Interactor3DMeasure();
  
  // return true if to have parent Interactor3D continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseUpEvent  ( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessKeyDownEvent  ( wxKeyEvent& event, RenderView* view );
  
protected:
  virtual void UpdateCursor( wxEvent& event, wxWindow* wnd );
  
  int  m_nMousePosX;
  int  m_nMousePosY;

  bool m_bSelectRegion;
};

#endif


