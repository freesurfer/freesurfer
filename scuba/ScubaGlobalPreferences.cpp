#include "ScubaGlobalPreferences.h"
#include "PreferencesManager.h"

ScubaGlobalPreferences& 
ScubaGlobalPreferences::GetPreferences() {

  static ScubaGlobalPreferences* sPreferences = NULL;
  if( NULL == sPreferences ) {
    sPreferences = new ScubaGlobalPreferences();


    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( *sPreferences, "SaveGlobalPreferences", 0, "", 
			   "Saves preferences." );
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
}

void
ScubaGlobalPreferences::WritePreferences () {

  PreferencesManager& prefsMgr = PreferencesManager::GetManager();

  PreferencesManager::IntPrefValue 
    viewFlipLeftRightInYZValue( mbViewFlipLeftRightInYZ );
  prefsMgr.SetValue( "ViewFlipLeftRight", viewFlipLeftRightInYZValue );

  prefsMgr.SaveFile();
}
