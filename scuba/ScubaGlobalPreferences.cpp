#include "string_fixed.h"
#include "ScubaGlobalPreferences.h"
#include "PreferencesManager.h"

using namespace std;

ScubaGlobalPreferences& 
ScubaGlobalPreferences::GetPreferences() {

  static ScubaGlobalPreferences* sPreferences = NULL;
  if( NULL == sPreferences ) {
    sPreferences = new ScubaGlobalPreferences();


    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( *sPreferences, "SaveGlobalPreferences", 0, "", 
			   "Saves preferences." );
    commandMgr.AddCommand( *sPreferences, "GetPreferencesValue", 1, "key", 
			   "Return a preferences value." );
    commandMgr.AddCommand( *sPreferences, "SetPreferencesValue", 2, 
			   "key value", "Set a preferences value." );
  }

  return *sPreferences;
}

ScubaGlobalPreferences::ScubaGlobalPreferences () {

  ReadPreferences();
}

void
ScubaGlobalPreferences::SavePreferences () {

  WritePreferences();
}

TclCommandManager::TclCommandResult
ScubaGlobalPreferences::DoListenToTclCommand ( char* isCommand, int iArgc, 
					       char** iasArgv ) {

  // SaveGlobalPreferences
  if( 0 == strcmp( isCommand, "SaveGlobalPreferences" ) ) {
    SavePreferences();
  }

  // GetPreferencesValue <key>
  if( 0 == strcmp( isCommand, "GetPreferencesValue" ) ) {

    string sKey = iasArgv[1];
    if( sKey == GetStringForKey( ViewFlipLeftRight ) ||
	sKey == GetStringForKey( ShowConsole ) ||
	sKey == GetStringForKey( AutoConfigureView ) ||
	sKey == GetStringForKey( KeyInPlaneX ) ||
	sKey == GetStringForKey( KeyInPlaneY ) ||
	sKey == GetStringForKey( KeyInPlaneZ ) ||
	sKey == GetStringForKey( KeyCycleViewsInFrame ) ||
	sKey == GetStringForKey( DrawCoordinateOverlay ) ||
	sKey == GetStringForKey( DrawPlaneIntersections ) ||
	sKey == GetStringForKey( KeyMoveViewLeft ) ||
	sKey == GetStringForKey( KeyMoveViewRight ) ||
	sKey == GetStringForKey( KeyMoveViewUp ) ||
	sKey == GetStringForKey( KeyMoveViewDown ) ||
	sKey == GetStringForKey( KeyMoveViewIn ) ||
	sKey == GetStringForKey( KeyMoveViewOut ) ||
	sKey == GetStringForKey( KeyZoomViewIn ) ||
	sKey == GetStringForKey( KeyZoomViewOut ) ) {
      
      PreferencesManager& prefsMgr = PreferencesManager::GetManager();
      string sValue = prefsMgr.GetValue( sKey );
      sReturnFormat = "s";
      sReturnValues = sValue;

    } else {

      sResult = "bad value for key \"" + string(iasArgv[1]) + "\", " +
	"not recognized";
      return error;	
    }

  }

  // SetPreferencesValue <key> <value>
  if( 0 == strcmp( isCommand, "SetPreferencesValue" ) ) {

    string sKey = iasArgv[1];
    if( sKey == GetStringForKey( ViewFlipLeftRight ) ||
	sKey == GetStringForKey( ShowConsole ) ||
	sKey == GetStringForKey( AutoConfigureView ) ||
	sKey == GetStringForKey( DrawCoordinateOverlay ) ||
	sKey == GetStringForKey( DrawPlaneIntersections ) ) {

      bool bValue;

      if( 0 == strcmp( iasArgv[2], "true" ) || 
	  0 == strcmp( iasArgv[2], "1" )) {
	bValue = true;
      } else if( 0 == strcmp( iasArgv[2], "false" ) ||
		 0 == strcmp( iasArgv[2], "0" ) ) {
	bValue = false;
      } else {
	sResult = "bad value for key " + string(iasArgv[1]) + ", \"" +
	  string(iasArgv[2]) + "\", should be true, 1, false, or 0";
	return error;	
      }

      PreferencesManager& prefsMgr = PreferencesManager::GetManager();
      PreferencesManager::IntPrefValue prefValue( bValue );
      prefsMgr.SetValue( sKey, prefValue );

    } else if( sKey == GetStringForKey( KeyInPlaneX ) ||
	       sKey == GetStringForKey( KeyInPlaneY ) ||
	       sKey == GetStringForKey( KeyInPlaneZ ) ||
	       sKey == GetStringForKey( KeyCycleViewsInFrame ) ) {

      string sValue = iasArgv[2];

      PreferencesManager& prefsMgr = PreferencesManager::GetManager();
      PreferencesManager::StringPrefValue prefValue( sValue );
      prefsMgr.SetValue( sKey, prefValue );
    }

    // Send a broadcast to our listeners to notify that a prefs values
    // has changed.
    string sValue( iasArgv[2] );
    SendBroadcast( sKey, (void*)&sValue );
  }

  return ok;
}


bool 
ScubaGlobalPreferences::GetPrefAsBool ( PrefKey iKey ) {

  if( iKey == ViewFlipLeftRight  ||
      iKey == ShowConsole ||
      iKey == AutoConfigureView ||
      iKey == DrawCoordinateOverlay ||
      iKey == DrawPlaneIntersections ) {
  
    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    string sValue = prefsMgr.GetValue( GetStringForKey( iKey ) );
    PreferencesManager::IntPrefValue value( sValue );
    return value.GetValue();

  } else {
    throw runtime_error( "Not a bool key" );
  }

  return false;
}

string 
ScubaGlobalPreferences::GetPrefAsString ( PrefKey iKey ) {

  if( iKey == KeyInPlaneX ||
      iKey == KeyInPlaneY ||
      iKey == KeyInPlaneZ ||
      iKey == KeyCycleViewsInFrame ||
      iKey == KeyMoveViewLeft ||
      iKey == KeyMoveViewRight ||
      iKey == KeyMoveViewUp ||
      iKey == KeyMoveViewDown ||
      iKey == KeyMoveViewIn ||
      iKey == KeyMoveViewOut ||
      iKey == KeyZoomViewIn ||
      iKey == KeyZoomViewOut ) {
  
    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    string sValue = prefsMgr.GetValue( GetStringForKey( iKey ) );
    PreferencesManager::StringPrefValue value( sValue );
    return value.GetValue();

  } else {
    throw runtime_error( "Not a string key" );
  }

  return "";
}

string 
ScubaGlobalPreferences::GetStringForKey ( PrefKey iKey ) {

  switch( iKey ) {
  case ShowConsole:                return "ShowConsole";                 break;
  case AutoConfigureView:          return "AutoConfigureView";           break;
  case ViewFlipLeftRight:          return "ViewFlipLeftRight";           break;
  case KeyInPlaneX:                return "KeyInPlaneX";                 break;
  case KeyInPlaneY:                return "KeyInPlaneY";                 break;
  case KeyInPlaneZ:                return "KeyInPlaneZ";                 break;
  case KeyCycleViewsInFrame:       return "KeyCycleViewsInFrame";        break;
  case DrawCoordinateOverlay:      return "DrawCoordinateOverlay";       break;
  case DrawPlaneIntersections:     return "DrawPlaneIntersections";      break;
  case KeyMoveViewLeft:            return "KeyMoveViewLeft";             break;
  case KeyMoveViewRight:           return "KeyMoveViewRight";            break;
  case KeyMoveViewUp:              return "KeyMoveViewUp";               break;
  case KeyMoveViewDown:            return "KeyMoveViewDown";             break;
  case KeyMoveViewIn:              return "KeyMoveViewIn";               break;
  case KeyMoveViewOut:             return "KeyMoveViewOut";              break;
  case KeyZoomViewIn:              return "KeyZoomViewIn";               break;
  case KeyZoomViewOut:             return "KeyZoomViewOut";              break;
  default:
    throw runtime_error( "Invalid key" );
  }
  return "";
}


void
ScubaGlobalPreferences::ReadPreferences () {

  // Here we just register our prefs with default values. We won't
  // actually read them in until instructed to by the tcl code.
  PreferencesManager& prefsMgr = PreferencesManager::GetManager();

  PreferencesManager::IntPrefValue viewFlipLeftRightInYZValue( true );
  prefsMgr.RegisterValue( GetStringForKey( ViewFlipLeftRight ), 
			  "Flip the view in the right/left direction to mimic "
			  "neurological style display.", 
			  viewFlipLeftRightInYZValue );

  PreferencesManager::IntPrefValue showConsole( true );
  prefsMgr.RegisterValue( GetStringForKey( ShowConsole ), 
			  "Show the tkcon console on startup.",
			  showConsole );

  PreferencesManager::IntPrefValue autoConfigure( true );
  prefsMgr.RegisterValue( GetStringForKey( AutoConfigureView ), 
			  "Automatically set up the view when the view "
			  "configuration is changed.",
			  autoConfigure );

  PreferencesManager::StringPrefValue inPlaneX( "x" );
  prefsMgr.RegisterValue( GetStringForKey( KeyInPlaneX ), 
			  "Key to change in plane to X in the view.",
			  inPlaneX );

  PreferencesManager::StringPrefValue inPlaneY( "y" );
  prefsMgr.RegisterValue( GetStringForKey( KeyInPlaneY ), 
			  "Key to change in plane to Y in the view.",
			  inPlaneY );

  PreferencesManager::StringPrefValue inPlaneZ( "z" );
  prefsMgr.RegisterValue( GetStringForKey( KeyInPlaneZ ), 
			  "Key to change in plane to Z in the view.",
			  inPlaneZ );

  PreferencesManager::StringPrefValue cycleKey( "q" );
  prefsMgr.RegisterValue( GetStringForKey( KeyCycleViewsInFrame ), 
			  "Key to cycle view in a frame.", cycleKey );

  PreferencesManager::IntPrefValue drawCoordinateOverlay( true );
  prefsMgr.RegisterValue( GetStringForKey( DrawCoordinateOverlay ), 
			  "Draw the coordinate overlay in views.", 
			  drawCoordinateOverlay );

  PreferencesManager::IntPrefValue drawPlaneIntersections( true );
  prefsMgr.RegisterValue( GetStringForKey( DrawPlaneIntersections ), 
			  "Draw the intersections of other views with this one as lines.", 
			  drawPlaneIntersections );

  PreferencesManager::StringPrefValue moveViewLeft( "Left" );
  prefsMgr.RegisterValue( "KeyMoveViewLeft", 
			  "Key to move the view to the left.",
			  moveViewLeft );

  PreferencesManager::StringPrefValue moveViewRight( "Right" );
  prefsMgr.RegisterValue( "KeyMoveViewRight", 
			  "Key to move the view to the right.",
			  moveViewRight );

  PreferencesManager::StringPrefValue moveViewUp( "Up" );
  prefsMgr.RegisterValue( "KeyMoveViewUp", 
			  "Key to move the view up.",
			  moveViewUp );

  PreferencesManager::StringPrefValue moveViewDown( "Down" );
  prefsMgr.RegisterValue( "KeyMoveViewDown", 
			  "Key to move the view down.",
			  moveViewDown );

  PreferencesManager::StringPrefValue moveViewIn( "Prior" );
  prefsMgr.RegisterValue( "KeyMoveViewIn", 
			  "Key to move the view in in plane.",
			  moveViewIn );

  PreferencesManager::StringPrefValue moveViewOut( "Next" );
  prefsMgr.RegisterValue( "KeyMoveViewOut", 
			  "Key to move the view out in plane.",
			  moveViewOut );

  PreferencesManager::StringPrefValue zoomViewIn( "equal" );
  prefsMgr.RegisterValue( "KeyZoomViewIn", 
			  "Key to zoom the view in in plane.",
			  zoomViewIn );

  PreferencesManager::StringPrefValue zoomViewOut( "minus" );
  prefsMgr.RegisterValue( "KeyZoomViewOut", 
			  "Key to zoom the view out in plane.",
			  zoomViewOut );

}

void
ScubaGlobalPreferences::WritePreferences () {

  // Write out the file.
  PreferencesManager& prefsMgr = PreferencesManager::GetManager();
  prefsMgr.SaveFile();
}
