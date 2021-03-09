/**
 * @brief Interactor for measurement in 3D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef Interactor3DVolumeCrop_h
#define Interactor3DVolumeCrop_h

#include "Interactor3D.h"

class Interactor3DVolumeCrop : public Interactor3D
{
public:
  Interactor3DVolumeCrop(QObject* parent);
  ~Interactor3DVolumeCrop();

  // return true if to have parent Interactor3D continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseUpEvent  ( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessKeyDownEvent  ( QKeyEvent* event, RenderView* view );

protected:
  virtual void UpdateCursor( QEvent* event, QWidget* wnd );

  int  m_nMousePosX;
  int  m_nMousePosY;

  bool m_bCropping;
};

#endif


