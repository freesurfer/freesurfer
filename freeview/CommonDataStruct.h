/**
 * @file  CommonDataStruct.h
 * @brief Mix-in class for message receiving
 *
 * Simple mix-in class for use with the Broadcaster class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/04/30 21:31:05 $
 *    $Revision: 1.6 $
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


#ifndef CommonDataStruct_h
#define CommonDataStruct_h

#include <vector>
#include <wx/colour.h>

struct SettingsGeneral
{
  wxColour  BackgroundColor;
  wxColour  CursorColor;
  int       CursorStyle;  
};

struct Settings2D
{
  bool      SyncZoomFactor;
};

struct SettingsScreenshot
{
  bool HideCursor;
  bool HideCoords;
  int  Magnification;
  bool AntiAliasing;
};

struct RotationElement
{
  int     Plane;
  double  Angle;
  double  Point[3];
};

#endif
