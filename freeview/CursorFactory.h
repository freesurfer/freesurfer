/**
 * @file  Cursor2D.h
 * @brief Cursor creator/loader for 2D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:46 $
 *    $Revision: 1.14 $
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

#ifndef CursorFactory_h
#define CursorFactory_h

#include <QCursor>

class CursorFactory
{
public:
  CursorFactory();

  static void Initialize();

  static QCursor CursorPan;
  static QCursor CursorZoom;

  static QCursor CursorPencil;
  static QCursor CursorFill;
  static QCursor CursorPolyline;
  static QCursor CursorColorPicker;
  static QCursor CursorMeasureLine;
  static QCursor CursorMeasureRectangle;
  static QCursor CursorMeasurePolyline;
  static QCursor CursorGrab;
  static QCursor CursorContour;
  static QCursor CursorAdd;
  static QCursor CursorRemove;
};

#endif


