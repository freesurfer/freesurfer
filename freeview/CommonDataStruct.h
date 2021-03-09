/**
 * @brief A few common structures to hold settings
 *
 * Simple mix-in class for use with the Broadcaster class so text
 * messages with a pointer data can be sent to a list of listeners.
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
  bool HideScaleBar;
  bool AntiAliasing;
  int  Magnification;
  bool AutoTrim;
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
