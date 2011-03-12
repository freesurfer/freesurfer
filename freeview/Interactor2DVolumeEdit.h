/**
 * @file  Interactor2DVolumeEdit.h
 * @brief Interactor for editing volume in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:48 $
 *    $Revision: 1.10 $
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

class Layer;

class Interactor2DVolumeEdit : public Interactor2D
{
  Q_OBJECT
public:
  Interactor2DVolumeEdit( const QString& layerTypeName, QObject* parent );
  virtual ~Interactor2DVolumeEdit();

  // return true if to have parent Interactor2D continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseUpEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessKeyDownEvent( QKeyEvent* event, RenderView* view );
  virtual bool ProcessKeyUpEvent( QKeyEvent* event, RenderView* view );

protected:
  void UpdateCursor( QEvent* event, QWidget* wnd );
  void ProcessContextMenu( QMouseEvent* event );

  bool m_bEditing;

  QString m_strLayerTypeName;

  QList<double>  m_dPolylinePoints;
};

#endif


