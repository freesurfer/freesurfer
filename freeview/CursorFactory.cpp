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

#include "CursorFactory.h"
#include <QPixmap>

QCursor CursorFactory::CursorPencil    = QCursor();
QCursor CursorFactory::CursorFill      = QCursor();
QCursor CursorFactory::CursorPolyline  = QCursor();
QCursor CursorFactory::CursorPan       = QCursor();
QCursor CursorFactory::CursorZoom      = QCursor();
QCursor CursorFactory::CursorGrab      = QCursor();
QCursor CursorFactory::CursorColorPicker  = QCursor();
QCursor CursorFactory::CursorMeasureLine  = QCursor();
QCursor CursorFactory::CursorMeasureRectangle = QCursor();
QCursor CursorFactory::CursorMeasurePolyline  = QCursor();
QCursor CursorFactory::CursorContour   = QCursor();
QCursor CursorFactory::CursorAdd       = QCursor();
QCursor CursorFactory::CursorRemove    = QCursor();

CursorFactory::CursorFactory()
{
  Initialize();
}

void CursorFactory::Initialize()
{
  CursorPencil = QCursor( QPixmap(":resource/icons/cursor_pencil.png"), 0, 0 );
  CursorFill = QCursor( QPixmap(":resource/icons/cursor_fill.png"), 2, 23 );
  CursorPolyline = QCursor( QPixmap(":resource/icons/cursor_polyline.png"), 0, 0 );
  CursorPan = QCursor( QPixmap(":resource/icons/cursor_pan.png"), 11, 11 );
  CursorZoom = QCursor( QPixmap(":resource/icons/cursor_zoom.png"), 11, 11 );
  CursorColorPicker = QCursor( QPixmap(":resource/icons/cursor_colorpicker.png"), 0, 23 );
  CursorMeasureLine = QCursor( QPixmap(":resource/icons/cursor_measure_line.png"), 0, 0 );
  CursorGrab = QCursor( QPixmap(":resource/icons/cursor_grab.png"), 11, 11 );
  CursorMeasureRectangle = QCursor( QPixmap(":resource/icons/cursor_measure_rectangle.png"), 0, 0 );
  CursorContour = QCursor( QPixmap(":resource/icons/cursor_contour.png"), 0, 0 );
  CursorMeasurePolyline = QCursor( QPixmap(":resource/icons/cursor_measure_polyline.png"), 0, 0 );
  CursorAdd = QCursor( QPixmap(":resource/icons/cursor_add.png"), 0, 0 );
  CursorRemove = QCursor( QPixmap(":resource/icons/cursor_remove.png"), 0, 0 );
}
