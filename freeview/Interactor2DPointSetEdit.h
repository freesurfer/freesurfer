/**
 * @file  Interactor2DPointSetEdit.h
 * @brief Interactor for editing way points in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:48 $
 *    $Revision: 1.1 $
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

#ifndef Interactor2DPointSetEdit_h
#define Interactor2DPointSetEdit_h

#include "Interactor2D.h"

class Interactor2DPointSetEdit : public Interactor2D
{
public:
  Interactor2DPointSetEdit(QObject* parent);
  virtual ~Interactor2DPointSetEdit();

  // return true if to have parent Interactor2DPointSetEdit continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseUpEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessKeyDownEvent( QKeyEvent* event, RenderView* view );

protected:
  void UpdateCursor( QEvent* event, QWidget* wnd );

  bool m_bEditing;
  int  m_nCurrentIndex;
};

#endif


