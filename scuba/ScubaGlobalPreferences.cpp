#include <string>
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
	sKey == GetStringForKey( KeyInPlaneX ) ||
	sKey == GetStringForKey( KeyInPlaneY ) ||
	sKey == GetStringForKey( KeyInPlaneZ ) ||
	sKey == GetStringForKey( KeyCycleViewsInFrame ) ) {
      
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
	sKey == GetStringForKey( ShowConsole ) ) {

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

    } else if( sKey == GetStringForKey( ViewFlipLeftRight ) ||
	       sKey == GetStringForKey( ShowConsole ) ) {

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
    }
  }

  return ok;
}

bool 
ScubaGlobalPreferences::GetPrefAsBool ( PrefKey iKey ) {

  if( iKey == ViewFlipLeftRight  ||
      iKey == ShowConsole ) {
  
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
ScubaGlobalPreferences::GetStringForKey ( PrefKey iKey ) {

  switch( iKey ) {
  case ShowConsole:
    return "ShowConsole";
    break;
  case ViewFlipLeftRight:
    return "ViewFlipLeftRight";
    break;
  case KeyInPlaneX:
    return "KeyInPlaneX";
    break;
  case KeyInPlaneY:
    return "KeyInPlaneY";
    break;
  case KeyInPlaneZ:
    return "KeyInPlaneZ";
    break;
  case KeyCycleViewsInFrame:
    return "KeyCycleViewsInFrame";
    break;
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
}

void
ScubaGlobalPreferences::WritePreferences () {

  // Write out the file.
  PreferencesManager& prefsMgr = PreferencesManager::GetManager();
  prefsMgr.SaveFile();
}
