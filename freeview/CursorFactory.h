/**
 * @file  Cursor2D.h
 * @brief Cursor creator/loader for 2D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 22:00:36 $
 *    $Revision: 1.13 $
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
 */

#ifndef CursorFactory_h
#define CursorFactory_h

#include <wx/cursor.h>

class CursorFactory
{
public:
  CursorFactory();

  static void Initialize();

  static wxCursor CursorPan;
  static wxCursor CursorZoom;

  static wxCursor CursorPencil;
  static wxCursor CursorFill;
  static wxCursor CursorPolyline;
  static wxCursor CursorColorPicker;
  static wxCursor CursorMeasureLine;
  static wxCursor CursorMeasureRectangle;
  static wxCursor CursorMeasurePolyline;
  static wxCursor CursorGrab;
  static wxCursor CursorContour;
  static wxCursor CursorAdd;
  static wxCursor CursorRemove;
};

#endif


