/**
 * @file  CommonDataStruct.h
 * @brief A few common structures to hold settings
 *
 * Simple mix-in class for use with the Broadcaster class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:45 $
 *    $Revision: 1.13 $
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

#include <QString>
#include <QColor>

struct SettingsGeneral
{
  QColor    BackgroundColor;
  QColor    CursorColor;
  int       CursorStyle;  
  bool      SaveCopy;
};

struct Settings2D
{
  bool      SyncZoomFactor;
};

struct SettingsScreenshot
{
  bool HideCursor;
  bool HideCoords;
  bool AntiAliasing;
  int  Magnification;
};

struct SettingsMovieFrames
{
  QString   OutputLocation;
  QString   OutputExtension;
  double    AngleStep;
  int       StepCount;
};

struct RotationElement
{
  int     Plane;
  double  Angle;
  double  Point[3];
  int     SampleMethod;
};

#endif
