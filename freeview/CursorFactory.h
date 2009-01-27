/**
 * @file  Cursor2D.h
 * @brief Annotation class for 2D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/01/27 18:27:24 $
 *    $Revision: 1.4 $
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
};

#endif


