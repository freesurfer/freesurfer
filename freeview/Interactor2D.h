/**
 * @brief Base Interactor class to manage mouse and key input on 2D render view.
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

#ifndef Interactor2D_h
#define Interactor2D_h

#include "Interactor.h"

class QString;
class Layer;
class Interactor2D : public Interactor
{
  Q_OBJECT
public:
  Interactor2D( QObject* parent = NULL );
  virtual ~Interactor2D();

  // return true if to have parent Interactor2D continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseUpEvent  ( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessKeyDownEvent  ( QKeyEvent* event, RenderView* view );

  virtual void ProcessPostMouseWheelEvent ( QWheelEvent* event, RenderView* view );
  virtual void ProcessPostMouseMoveEvent  ( QMouseEvent* event, RenderView* view );

signals:
  void Error( const QString& message, Layer* layer = NULL );
  void CursorLocationClicked();

protected:
  int  m_nMousePosX;
  int  m_nMousePosY;

  bool m_bWindowLevel;
  bool m_bChangeSlice;
  bool m_bMovingCursor;
  bool m_bSelecting;
};

#endif


