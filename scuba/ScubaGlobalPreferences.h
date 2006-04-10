#ifndef ScubaGlobalPreferences_h
#define ScubaGlobalPreferences_h

#include "TclCommandManager.h"
#include "Broadcaster.h"

// This class acts as a bridge between PreferencesManager and the rest
// of the Scuba objects. Scuba objects can add themselves as listeners
// to this object to get messages when preference values have been
// changed by the user. This class also acts as the bridge between
// prefs values and tcl.

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
		 KeyCycleViewsInFrame, KeyShuffleLayers,
		 KeyMouseButtonOne, KeyMouseButtonTwo, KeyMouseButtonThree,
		 KeyTurnOffVisibilityInTopVisibleLayer,
		 KeyTurnOnVisibilityInTopInvisibleLayer,
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
