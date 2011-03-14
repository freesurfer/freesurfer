/**
 * @file  Interactor2DMeasure.h
 * @brief Interactor for measure tool in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:57 $
 *    $Revision: 1.8 $
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


