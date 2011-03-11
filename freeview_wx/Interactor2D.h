/**
 * @file  Interactor2D.h
 * @brief Base Interactor class to manage mouse and key input on 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:37 $
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

#ifndef Interactor2D_h
#define Interactor2D_h

#include "Interactor.h"

class Interactor2D : public Interactor
{
public:
  Interactor2D();
  virtual ~Interactor2D();

  // return true if to have parent Interactor2D continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseUpEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessKeyDownEvent( wxKeyEvent& event, RenderView* view );

  virtual void ProcessPostMouseWheelEvent( wxMouseEvent& event, RenderView* view );
  virtual void ProcessPostMouseMoveEvent( wxMouseEvent& event, RenderView* view );

protected:
  int  m_nMousePosX;
  int  m_nMousePosY;

  bool m_bWindowLevel;
  bool m_bChangeSlice;
  bool m_bMovingCursor;
  bool m_bSelecting;
};

#endif


