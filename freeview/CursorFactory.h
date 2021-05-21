/**
 * @brief Cursor creator/loader for 2D view.
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


