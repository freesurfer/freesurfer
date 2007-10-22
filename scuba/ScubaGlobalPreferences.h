/**
 * @file  ScubaGlobalPreferences.h
 * @brief Scuba specific preferences values
 *
 * This class acts as a bridge between PreferencesManager and the rest
 * of the Scuba objects. Scuba objects can add themselves as listeners
 * to this object to get messages when preference values have been
 * changed by the user. This class also acts as the bridge between
 * prefs values and tcl.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/22 04:39:28 $
 *    $Revision: 1.21 $
 *
 * Copyright (C) 2002-2007,
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


#ifndef ScubaGlobalPreferences_h
#define ScubaGlobalPreferences_h

#include "TclCommandManager.h"
#include "Broadcaster.h"

class ScubaGlobalPreferences : public TclCommandListener, public Broadcaster {

public:

  enum PrefKey { ShowConsole,
                 AutoConfigureView, AutoReportInfoForLUT,
                 AutoDrawZeroClearForLUT,
                 ViewFlipLeftRight,
                 KeyInPlaneX, KeyInPlaneY, KeyInPlaneZ,
                 KeyMoveViewLeft, KeyMoveViewRight,
                 KeyMoveViewUp, KeyMoveViewDown,
                 KeyMoveViewIn, KeyMoveViewOut,
                 KeyZoomViewIn, KeyZoomViewOut,
                 KeyToolNavigate, KeyToolPlane, KeyToolMarker,
                 KeyToolEditVoxel,
                 KeyToolEditROI, KeyToolStraightPath, KeyToolEdgePath,
                 KeyCycleViewsInFrame, KeyShuffleLayers,
                 KeyMouseButtonOne, KeyMouseButtonTwo, KeyMouseButtonThree,
                 KeyTurnOffVisibilityInTopVisibleLayer,
                 KeyTurnOnVisibilityInTopInvisibleLayer,
                 KeyToggleVisibilityInTopmostLayer,
                 KeyToggleVisibilityInTopmostUnlockedLayer,
                 DrawCoordinateOverlay, DrawPlaneIntersections, DrawMarkers,
                 DrawPaths, SelectedTool, LockOnCursor,
                 ShowFPS, UserStructureList };

  // Gets the static reference to this class.
  static ScubaGlobalPreferences& GetPreferences();

  void SavePreferences ();

  void SetPreferencesValue ( PrefKey iKey, bool ibValue );
  void SetPreferencesValue ( PrefKey iKey, std::string isValue );

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  bool        GetPrefAsBool   ( PrefKey iKey );
  std::string GetPrefAsString ( PrefKey iKey );

protected:

  std::string GetStringForKey ( PrefKey iKey );

  ScubaGlobalPreferences ();

  void ReadPreferences ();
  void WritePreferences ();

};


#endif
