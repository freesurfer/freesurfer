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

  mbViewFlipLeftRightInYZ = true;

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
    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    string sValue = prefsMgr.GetValue( sKey );
    sReturnFormat = "s";
    sReturnValues = sValue;

  }

  // SetPreferencesValue <key> <value>
  if( 0 == strcmp( isCommand, "SetPreferencesValue" ) ) {

    string sKey = iasArgv[1];
    if( sKey == "ViewFlipLeftRight" ) {

      if( 0 == strcmp( iasArgv[2], "true" ) || 
	  0 == strcmp( iasArgv[2], "1" )) {
	mbViewFlipLeftRightInYZ = true;
      } else if( 0 == strcmp( iasArgv[2], "false" ) ||
		 0 == strcmp( iasArgv[2], "0" ) ) {
	mbViewFlipLeftRightInYZ = false;
      } else {
	sResult = "bad value for key " + string(iasArgv[1]) + ", \"" +
	  string(iasArgv[2]) + "\", should be true, 1, false, or 0";
	return error;	
      }
    }

    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    string sValue = prefsMgr.GetValue( sKey );
    sReturnFormat = "s";
    sReturnValues = sValue;

  }

  return ok;
}


void
ScubaGlobalPreferences::ReadPreferences () {

  PreferencesManager& prefsMgr = PreferencesManager::GetManager();

  PreferencesManager::IntPrefValue 
    viewFlipLeftRightInYZValue( mbViewFlipLeftRightInYZ );
  prefsMgr.RegisterValue( "ViewFlipLeftRight", 
		       "Flip the view in the right/left direction to mimic "
		       "neurological style display.", 
		       viewFlipLeftRightInYZValue );
  viewFlipLeftRightInYZValue.SetFromString( prefsMgr.GetValue( "ViewFlipLeftRight" ) );
  mbViewFlipLeftRightInYZ = viewFlipLeftRightInYZValue.GetValue();

  PreferencesManager::StringPrefValue inPlaneX( "x" );
  prefsMgr.RegisterValue( "key-InPlaneX", 
			  "Key to change in plane to X in the view.",
			  inPlaneX );
  msInPlaneXKey = prefsMgr.GetValue( "key-InPlaneX" );

  PreferencesManager::StringPrefValue inPlaneY( "y" );
  prefsMgr.RegisterValue( "key-InPlaneY", 
			  "Key to change in plane to Y in the view.",
			  inPlaneY );
  msInPlaneYKey = prefsMgr.GetValue( "key-InPlaneY" );

  PreferencesManager::StringPrefValue inPlaneZ( "z" );
  prefsMgr.RegisterValue( "key-InPlaneZ", 
			  "Key to change in plane to Z in the view.",
			  inPlaneZ );
  msInPlaneZKey = prefsMgr.GetValue( "key-InPlaneZ" );

  PreferencesManager::StringPrefValue cycleKey( "q" );
  prefsMgr.RegisterValue( "key-CycleViewsInFrame", 
			  "Key to cycle view in a frame.", cycleKey );
  msCycleKey = prefsMgr.GetValue( "key-CycleViewsInFrame" );
}

void
ScubaGlobalPreferences::WritePreferences () {

  PreferencesManager& prefsMgr = PreferencesManager::GetManager();

  PreferencesManager::IntPrefValue 
    viewFlipLeftRightInYZValue( mbViewFlipLeftRightInYZ );
  prefsMgr.SetValue( "ViewFlipLeftRight", viewFlipLeftRightInYZValue );

  prefsMgr.SaveFile();
}
