/**
 * @brief Interactor for measure tool in 2D render view.
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

#ifndef Interactor2DMeasure_h
#define Interactor2DMeasure_h

#include "Interactor2D.h"

class Region2D;

class Interactor2DMeasure : public Interactor2D
{
  Q_OBJECT
public:
  Interactor2DMeasure( QObject* parent );
  virtual ~Interactor2DMeasure();

  // return true if to have parent Interactor2D continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent  ( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseUpEvent    ( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseMoveEvent  ( QMouseEvent* event, RenderView* view );
  virtual bool ProcessKeyDownEvent    ( QKeyEvent* event, RenderView* view );
  virtual bool ProcessKeyUpEvent      ( QKeyEvent* event, RenderView* view );

signals:

protected:
  void UpdateCursor( QEvent* event, QWidget* wnd );

  bool        m_bEditing;
  bool        m_bDrawing;
  int         m_nPointIndex;
  Region2D*   m_region;
};

#endif


