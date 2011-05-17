/**
 * @file  Interactor3D.h
 * @brief Base Interactor class for 3D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/05/17 14:20:14 $
 *    $Revision: 1.16 $
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
 *
 */

#ifndef Interactor3D_h
#define Interactor3D_h

#include "Interactor.h"

class SurfaceROI;

class Interactor3D : public Interactor
{
  Q_OBJECT
public:
  Interactor3D(QObject* parent = 0);
  virtual ~Interactor3D();

  // return true if to have parent Interactor3D continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent ( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseUpEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessKeyDownEvent( QKeyEvent* event, RenderView* view );

protected:
  virtual void UpdateCursor( QEvent* event, QWidget* wnd );
  bool IsInAction()
  {
    return m_bWindowLevel || m_bMoveSlice;
  }

  int  m_nMousePosX;
  int  m_nMousePosY;

  bool m_bWindowLevel;
  bool m_bMoveSlice;
  SurfaceROI*   m_surfaceROI;
};

#endif


