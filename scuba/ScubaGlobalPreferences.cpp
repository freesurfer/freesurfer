/**
 * @file  ScubaGlobalPreferences.cpp
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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:37 $
 *    $Revision: 1.31 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include "string_fixed.h"
#include "ScubaGlobalPreferences.h"
#include "PreferencesManager.h"
#include "ScubaKeyCombo.h"
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h> // strcmp
#ifdef __cplusplus
}
#endif

using namespace std;

ScubaGlobalPreferences&
ScubaGlobalPreferences::GetPreferences() {

  static ScubaGlobalPreferences* sPreferences = NULL;
  if ( NULL == sPreferences ) {
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

ScubaGlobalPreferences::ScubaGlobalPreferences () :
    Broadcaster( "ScubaGlobalPreferences" ) {

  ReadPreferences();
}

void
ScubaGlobalPreferences::SavePreferences () {

  WritePreferences();
}

TclCommandManager::TclCommandResult
ScubaGlobalPreferences::DoListenToTclCommand ( char* isCommand,
    int, char** iasArgv ) {

  // SaveGlobalPreferences
  if ( 0 == strcmp( isCommand, "SaveGlobalPreferences" ) ) {
    SavePreferences();
  }

  // GetPreferencesValue <key>
  if ( 0 == strcmp( isCommand, "GetPreferencesValue" ) ) {

    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    string sKey = iasArgv[1];

    // If it's a key preference, convert it to a proper ScubaKeyCombo
    // string first. This is just to convert key values from old
    // versions of the prefs file to nice new ones.
    if ( sKey == GetStringForKey( KeyInPlaneX ) ||
         sKey == GetStringForKey( KeyInPlaneY ) ||
         sKey == GetStringForKey( KeyInPlaneZ ) ||
         sKey == GetStringForKey( KeyCycleViewsInFrame ) ||
         sKey == GetStringForKey( KeyShuffleLayers ) ||
         sKey == GetStringForKey( KeyTurnOffVisibilityInTopVisibleLayer ) ||
         sKey == GetStringForKey( KeyTurnOnVisibilityInTopInvisibleLayer ) ||
         sKey == GetStringForKey( KeyToggleVisibilityInTopmostLayer ) ||
         sKey == GetStringForKey( KeyToggleVisibilityInTopmostUnlockedLayer ) ||
         sKey == GetStringForKey( KeyMouseButtonOne ) ||
         sKey == GetStringForKey( KeyMouseButtonTwo ) ||
         sKey == GetStringForKey( KeyMouseButtonThree ) ||
         sKey == GetStringForKey( KeyMoveViewLeft ) ||
         sKey == GetStringForKey( KeyMoveViewRight ) ||
         sKey == GetStringForKey( KeyMoveViewUp ) ||
         sKey == GetStringForKey( KeyMoveViewDown ) ||
         sKey == GetStringForKey( KeyMoveViewIn ) ||
         sKey == GetStringForKey( KeyMoveViewOut ) ||
         sKey == GetStringForKey( KeyZoomViewIn ) ||
         sKey == GetStringForKey( KeyZoomViewOut ) ||
         sKey == GetStringForKey( KeyToolNavigate ) ||
         sKey == GetStringForKey( KeyToolPlane ) ||
         sKey == GetStringForKey( KeyToolMarker ) ||
         sKey == GetStringForKey( KeyToolEditVoxel  ) ||
         sKey == GetStringForKey( KeyToolEditROI ) ||
         sKey == GetStringForKey( KeyToolStraightPath ) ||
         sKey == GetStringForKey( KeyToolEdgePath )
       ) {

      string sValue = prefsMgr.GetValue( sKey );

      // Make a key combo and return the value of the ToString
      // function.
      ScubaKeyCombo* keyCombo = ScubaKeyCombo::MakeKeyCombo();
      keyCombo->SetFromString( sValue );
      sValue = keyCombo->ToString();

      // Enclose the value in quotes if it's not already in quotes.
      sReturnFormat = "s";
      stringstream ssReturnValues;
      if ( sValue[0] != '\"' )
        ssReturnValues << "\"" << sValue << "\"";
      else
        ssReturnValues << sValue;
      sReturnValues = ssReturnValues.str();

      // Another key, just return the value.
    } else if ( sKey == GetStringForKey( ViewFlipLeftRight ) ||
                sKey == GetStringForKey( ShowConsole ) ||
                sKey == GetStringForKey( AutoConfigureView ) ||
                sKey == GetStringForKey( AutoReportInfoForLUT ) ||
                sKey == GetStringForKey( AutoDrawZeroClearForLUT ) ||
                sKey == GetStringForKey( DrawCoordinateOverlay ) ||
                sKey == GetStringForKey( DrawMarkers ) ||
                sKey == GetStringForKey( DrawPaths ) ||
                sKey == GetStringForKey( DrawPlaneIntersections ) ||
                sKey == GetStringForKey( ShowFPS ) ||
                sKey == GetStringForKey( SelectedTool ) ||
                sKey == GetStringForKey( LockOnCursor ) ||
                sKey == GetStringForKey( UserStructureList )
              ) {

      // Enclose the value in quotes if it's not already in quotes.
      string sValue = prefsMgr.GetValue( sKey );
      sReturnFormat = "s";
      stringstream ssReturnValues;
      if ( sValue[0] != '\"' )
        ssReturnValues << "\"" << sValue << "\"";
      else
        ssReturnValues << sValue;
      sReturnValues = ssReturnValues.str();

    } else {

      sResult = "bad value for key \"" + string(iasArgv[1]) + "\", " +
                "not recognized";
      return error;
    }

  }

  // SetPreferencesValue <key> <value>
  if ( 0 == strcmp( isCommand, "SetPreferencesValue" ) ) {

    string sKey = iasArgv[1];
    if ( sKey == GetStringForKey( ViewFlipLeftRight ) ||
         sKey == GetStringForKey( ShowConsole ) ||
         sKey == GetStringForKey( AutoConfigureView ) ||
         sKey == GetStringForKey( AutoReportInfoForLUT ) ||
         sKey == GetStringForKey( AutoDrawZeroClearForLUT ) ||
         sKey == GetStringForKey( DrawCoordinateOverlay ) ||
         sKey == GetStringForKey( DrawMarkers ) ||
         sKey == GetStringForKey( DrawPaths ) ||
         sKey == GetStringForKey( DrawPlaneIntersections ) ||
         sKey == GetStringForKey( LockOnCursor ) ||
         sKey == GetStringForKey( ShowFPS ) ) {

      bool bValue;

      if ( 0 == strcmp( iasArgv[2], "true" ) ||
           0 == strcmp( iasArgv[2], "1" )) {
        bValue = true;
      } else if ( 0 == strcmp( iasArgv[2], "false" ) ||
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

    } else if ( sKey == GetStringForKey( KeyInPlaneX ) ||
                sKey == GetStringForKey( KeyInPlaneY ) ||
                sKey == GetStringForKey( KeyInPlaneZ ) ||
                sKey == GetStringForKey( KeyCycleViewsInFrame ) ||
                sKey == GetStringForKey( KeyShuffleLayers ) ||
                sKey ==
                GetStringForKey( KeyTurnOffVisibilityInTopVisibleLayer ) ||
                sKey ==
                GetStringForKey( KeyTurnOnVisibilityInTopInvisibleLayer ) ||
                sKey == GetStringForKey( KeyToggleVisibilityInTopmostLayer ) ||
                sKey ==
                GetStringForKey( KeyToggleVisibilityInTopmostUnlockedLayer ) ||
                sKey == GetStringForKey( KeyMouseButtonOne ) ||
                sKey == GetStringForKey( KeyMouseButtonTwo ) ||
                sKey == GetStringForKey( KeyMouseButtonThree ) ||
                sKey == GetStringForKey( KeyMoveViewIn ) ||
                sKey == GetStringForKey( KeyMoveViewOut ) ||
                sKey == GetStringForKey( KeyMoveViewLeft ) ||
                sKey == GetStringForKey( KeyMoveViewRight ) ||
                sKey == GetStringForKey( KeyMoveViewUp ) ||
                sKey == GetStringForKey( KeyMoveViewDown ) ||
                sKey == GetStringForKey( SelectedTool ) ||
                sKey == GetStringForKey( UserStructureList ) ||
                sKey == GetStringForKey( KeyToolNavigate ) ||
                sKey == GetStringForKey( KeyToolPlane ) ||
                sKey == GetStringForKey( KeyToolMarker ) ||
                sKey == GetStringForKey( KeyToolEditVoxel  ) ||
                sKey == GetStringForKey( KeyToolEditROI ) ||
                sKey == GetStringForKey( KeyToolStraightPath ) ||
                sKey == GetStringForKey( KeyToolEdgePath ) ) {

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


void
ScubaGlobalPreferences::SetPreferencesValue ( PrefKey iKey, bool ibValue ) {

  if ( iKey == ViewFlipLeftRight  ||
       iKey == ShowConsole ||
       iKey == AutoConfigureView ||
       iKey == AutoReportInfoForLUT ||
       iKey == AutoDrawZeroClearForLUT ||
       iKey == DrawCoordinateOverlay ||
       iKey == DrawMarkers ||
       iKey == DrawPaths ||
       iKey == DrawPlaneIntersections ||
       iKey == LockOnCursor ||
       iKey == ShowFPS ) {

    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    PreferencesManager::IntPrefValue value( (int) ibValue );
    prefsMgr.SetValue( GetStringForKey( iKey ), value );

    // Send a broadcast to our listeners to notify that a prefs values
    // has changed.
    stringstream ssValue;
    ssValue << (int)ibValue;
    string sValue = ssValue.str();
    SendBroadcast( GetStringForKey( iKey ), (void*)&sValue );

  } else {
    throw runtime_error( "Not a bool key" );
  }
}

void
ScubaGlobalPreferences::SetPreferencesValue ( PrefKey iKey, string isValue ) {

  if ( iKey == KeyInPlaneX ||
       iKey == KeyInPlaneY ||
       iKey == KeyInPlaneZ ||
       iKey == KeyCycleViewsInFrame ||
       iKey == KeyShuffleLayers ||
       iKey == KeyTurnOffVisibilityInTopVisibleLayer ||
       iKey == KeyTurnOnVisibilityInTopInvisibleLayer ||
       iKey == KeyToggleVisibilityInTopmostLayer ||
       iKey == KeyToggleVisibilityInTopmostUnlockedLayer ||
       iKey == KeyMouseButtonOne ||
       iKey == KeyMouseButtonTwo ||
       iKey == KeyMouseButtonThree ||
       iKey == KeyMoveViewLeft ||
       iKey == KeyMoveViewRight ||
       iKey == KeyMoveViewUp ||
       iKey == KeyMoveViewDown ||
       iKey == KeyMoveViewIn ||
       iKey == KeyMoveViewOut ||
       iKey == KeyZoomViewIn ||
       iKey == KeyZoomViewOut ||
       iKey == KeyToolNavigate ||
       iKey == KeyToolPlane ||
       iKey == KeyToolMarker ||
       iKey == KeyToolEditVoxel ||
       iKey == KeyToolEditROI ||
       iKey == KeyToolStraightPath ||
       iKey == KeyToolEdgePath ||
       iKey == SelectedTool ||
       iKey == UserStructureList ) {

    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    PreferencesManager::StringPrefValue value( isValue );
    prefsMgr.SetValue( GetStringForKey( iKey ), value );

    // Send a broadcast to our listeners to notify that a prefs values
    // has changed.
    string sValue( isValue );
    SendBroadcast( GetStringForKey( iKey ), (void*)&isValue );

  } else {
    throw runtime_error( "Not a string key" );
  }
}

bool
ScubaGlobalPreferences::GetPrefAsBool ( PrefKey iKey ) {

  if ( iKey == ViewFlipLeftRight  ||
       iKey == ShowConsole ||
       iKey == AutoConfigureView ||
       iKey == AutoReportInfoForLUT ||
       iKey == AutoDrawZeroClearForLUT ||
       iKey == DrawCoordinateOverlay ||
       iKey == DrawMarkers ||
       iKey == DrawPaths ||
       iKey == DrawPlaneIntersections ||
       iKey == LockOnCursor ||
       iKey == ShowFPS ) {

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

  // If it's a key preference, convert it to a proper ScubaKeyCombo
  // string first. This is just to convert key values from old
  // versions of the prefs file to nice new ones.
  if ( iKey == KeyInPlaneX ||
       iKey == KeyInPlaneY ||
       iKey == KeyInPlaneZ ||
       iKey == KeyCycleViewsInFrame ||
       iKey == KeyShuffleLayers ||
       iKey == KeyTurnOffVisibilityInTopVisibleLayer ||
       iKey == KeyTurnOnVisibilityInTopInvisibleLayer ||
       iKey == KeyToggleVisibilityInTopmostLayer ||
       iKey == KeyToggleVisibilityInTopmostUnlockedLayer ||
       iKey == KeyMouseButtonOne ||
       iKey == KeyMouseButtonTwo ||
       iKey == KeyMouseButtonThree ||
       iKey == KeyMoveViewLeft ||
       iKey == KeyMoveViewRight ||
       iKey == KeyMoveViewUp ||
       iKey == KeyMoveViewDown ||
       iKey == KeyMoveViewIn ||
       iKey == KeyMoveViewOut ||
       iKey == KeyZoomViewIn ||
       iKey == KeyZoomViewOut ||
       iKey == KeyToolNavigate ||
       iKey == KeyToolPlane ||
       iKey == KeyToolMarker ||
       iKey == KeyToolEditVoxel ||
       iKey == KeyToolEditROI ||
       iKey == KeyToolStraightPath ||
       iKey == KeyToolEdgePath ) {

    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    string sValue = prefsMgr.GetValue( GetStringForKey( iKey ) );

    // Make a key combo and return the value of the ToString
    // function. This is just to convert key values from old versions
    // of the prefs file to nice new ones.
    ScubaKeyCombo* keyCombo = ScubaKeyCombo::MakeKeyCombo();
    keyCombo->SetFromString( sValue );
    sValue = keyCombo->ToString();

    PreferencesManager::StringPrefValue value( sValue );
    return value.GetValue();

    // Normal string preference.
  } else if ( iKey == SelectedTool ||
              iKey == UserStructureList ) {

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

  switch ( iKey ) {
  case ShowConsole:
    return "ShowConsole";
    break;
  case AutoConfigureView:
    return "AutoConfigureView";
    break;
  case AutoReportInfoForLUT:
    return "AutoReportInfoForLUT";
    break;
  case AutoDrawZeroClearForLUT:
    return "AutoDrawZeroClearForLUT";
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
  case KeyShuffleLayers:
    return "KeyShuffleLayers";
    break;
  case KeyTurnOffVisibilityInTopVisibleLayer:
    return "KeyTurnOffVisibilityInTopVisibleLayer";
    break;
  case KeyTurnOnVisibilityInTopInvisibleLayer:
    return "KeyTurnOnVisibilityInTopInvisibleLayer";
    break;
  case KeyToggleVisibilityInTopmostLayer:
    return "KeyToggleVisibilityInTopmostLayer";
    break;
  case KeyToggleVisibilityInTopmostUnlockedLayer:
    return "KeyToggleVisibilityInTopmostUnlockedLayer";
    break;
  case KeyMouseButtonOne:
    return "KeyMouseButtonOne";
    break;
  case KeyMouseButtonTwo:
    return "KeyMouseButtonTwo";
    break;
  case KeyMouseButtonThree:
    return "KeyMouseButtonThree";
    break;
  case DrawCoordinateOverlay:
    return "DrawCoordinateOverlay";
    break;
  case DrawMarkers:
    return "DrawMarkers";
    break;
  case DrawPaths:
    return "DrawPaths";
    break;
  case DrawPlaneIntersections:
    return "DrawPlaneIntersections";
    break;
  case KeyMoveViewLeft:
    return "KeyMoveViewLeft";
    break;
  case KeyMoveViewRight:
    return "KeyMoveViewRight";
    break;
  case KeyMoveViewUp:
    return "KeyMoveViewUp";
    break;
  case KeyMoveViewDown:
    return "KeyMoveViewDown";
    break;
  case KeyMoveViewIn:
    return "KeyMoveViewIn";
    break;
  case KeyMoveViewOut:
    return "KeyMoveViewOut";
    break;
  case KeyZoomViewIn:
    return "KeyZoomViewIn";
    break;
  case KeyZoomViewOut:
    return "KeyZoomViewOut";
    break;
  case KeyToolNavigate:
    return "KeyToolNavigate";
    break;
  case KeyToolPlane:
    return "KeyToolPlane";
    break;
  case KeyToolMarker:
    return "KeyToolMarker";
    break;
  case KeyToolEditVoxel:
    return "KeyToolEditVoxel";
    break;
  case KeyToolEditROI:
    return "KeyToolEditROI";
    break;
  case KeyToolStraightPath:
    return "KeyToolStraightPath";
    break;
  case KeyToolEdgePath:
    return "KeyToolEdgePath";
    break;
  case SelectedTool:
    return "SelectedTool";
    break;
  case UserStructureList:
    return "UserStructureList";
    break;
  case LockOnCursor:
    return "LockOnCursor";
    break;
  case ShowFPS:
    return "ShowFPS";
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

  // For our default key values, create a key combo, set it, and then
  // use the ToString function to get the default string.
  ScubaKeyCombo* keyCombo = ScubaKeyCombo::MakeKeyCombo();

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

  PreferencesManager::IntPrefValue autoReportInfo( true );
  prefsMgr.RegisterValue( GetStringForKey( AutoReportInfoForLUT ),
                          "Automatically set the Report Info checkbox for "
                          "a view on if the layer in the view is an LUT "
                          "volume.",
                          autoReportInfo );

  PreferencesManager::IntPrefValue autoDrawZeroClearForLUT( true );
  prefsMgr.RegisterValue( GetStringForKey( AutoDrawZeroClearForLUT ),
                          "Automatically set the Draw Zero Clear option for "
                          "layers that have their colormaps set to LUT.",
                          autoDrawZeroClearForLUT );

  keyCombo->SetFromString( "X" );
  PreferencesManager::StringPrefValue inPlaneX( keyCombo->ToString() );
  prefsMgr.RegisterValue( GetStringForKey( KeyInPlaneX ),
                          "Key to change in plane to X in the view.",
                          inPlaneX );

  keyCombo->SetFromString( "Y" );
  PreferencesManager::StringPrefValue inPlaneY( keyCombo->ToString() );
  prefsMgr.RegisterValue( GetStringForKey( KeyInPlaneY ),
                          "Key to change in plane to Y in the view.",
                          inPlaneY );

  keyCombo->SetFromString( "Z" );
  PreferencesManager::StringPrefValue inPlaneZ( keyCombo->ToString() );
  prefsMgr.RegisterValue( GetStringForKey( KeyInPlaneZ ),
                          "Key to change in plane to Z in the view.",
                          inPlaneZ );

  keyCombo->SetFromString( "Q" );
  PreferencesManager::StringPrefValue cycleKey( keyCombo->ToString() );
  prefsMgr.RegisterValue( GetStringForKey( KeyCycleViewsInFrame ),
                          "Key to cycle view in a frame.", cycleKey );

  keyCombo->SetFromString( "C" );
  PreferencesManager::StringPrefValue shuffleKey( keyCombo->ToString() );
  prefsMgr.RegisterValue( GetStringForKey( KeyShuffleLayers ),
                          "Key to shuffle layers in a view.", shuffleKey );

  keyCombo->SetFromString( "Alt V" );
  PreferencesManager::StringPrefValue
  turnOffVisibilityKey( keyCombo->ToString() );
  prefsMgr.
  RegisterValue( GetStringForKey( KeyTurnOffVisibilityInTopVisibleLayer ),
                 "Key to turn off visibility in the topmost "
                 "visible layer in a  view.", turnOffVisibilityKey );

  keyCombo->SetFromString( "Ctrl V" );
  PreferencesManager::StringPrefValue
  turnOnVisibilityKey( keyCombo->ToString() );
  prefsMgr.
  RegisterValue( GetStringForKey( KeyTurnOnVisibilityInTopInvisibleLayer ),
                 "Key to turn on visibility in the topmost "
                 "invisible layer in a  view.", turnOnVisibilityKey );

  keyCombo->SetFromString( "Ctrl T" );
  PreferencesManager::StringPrefValue
  toggleVisibilityKey( keyCombo->ToString() );
  prefsMgr.
  RegisterValue( GetStringForKey(KeyToggleVisibilityInTopmostLayer),
                 "Key to toggle visibility in the topmost "
                 "layer in a  view.", toggleVisibilityKey );

  keyCombo->SetFromString( "Alt T" );
  PreferencesManager::StringPrefValue
  toggleUnlockedVisibilityKey( keyCombo->ToString() );
  prefsMgr.
  RegisterValue( GetStringForKey(KeyToggleVisibilityInTopmostUnlockedLayer),
                 "Key to toggle visibility in the topmost "
                 "unlocked layer in a  view.",
                 toggleUnlockedVisibilityKey );

  keyCombo->SetFromString( "W" );
  PreferencesManager::StringPrefValue mouseOneKey( keyCombo->ToString() );
  prefsMgr.RegisterValue( GetStringForKey( KeyMouseButtonOne ),
                          "Key to emulate a mouse click with button one.", mouseOneKey );

  keyCombo->SetFromString( "E" );
  PreferencesManager::StringPrefValue mouseTwoKey( keyCombo->ToString() );
  prefsMgr.RegisterValue( GetStringForKey( KeyMouseButtonTwo ),
                          "Key to emulate a mouse click with button two.", mouseTwoKey );

  keyCombo->SetFromString( "R" );
  PreferencesManager::StringPrefValue mouseThreeKey( keyCombo->ToString() );
  prefsMgr.RegisterValue( GetStringForKey( KeyMouseButtonThree ),
                          "Key to emulate a mouse click with button three.", mouseThreeKey );

  PreferencesManager::IntPrefValue drawCoordinateOverlay( true );
  prefsMgr.RegisterValue( GetStringForKey( DrawCoordinateOverlay ),
                          "Draw the coordinate overlay in views.",
                          drawCoordinateOverlay );

  PreferencesManager::IntPrefValue drawMarkers( true );
  prefsMgr.RegisterValue( GetStringForKey( DrawMarkers ),
                          "Draw markers in views.",
                          drawMarkers );

  PreferencesManager::IntPrefValue drawPaths( true );
  prefsMgr.RegisterValue( GetStringForKey( DrawPaths ),
                          "Draw paths in views.",
                          drawPaths );

  PreferencesManager::IntPrefValue drawPlaneIntersections( true );
  prefsMgr.RegisterValue( GetStringForKey( DrawPlaneIntersections ),
                          "Draw the intersections of other views with this one as lines.",
                          drawPlaneIntersections );

  keyCombo->SetFromString( "Left" );
  PreferencesManager::StringPrefValue moveViewLeft( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyMoveViewLeft",
                          "Key to move the view to the left.",
                          moveViewLeft );

  keyCombo->SetFromString( "Right" );
  PreferencesManager::StringPrefValue moveViewRight( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyMoveViewRight",
                          "Key to move the view to the right.",
                          moveViewRight );

  keyCombo->SetFromString( "Up" );
  PreferencesManager::StringPrefValue moveViewUp( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyMoveViewUp",
                          "Key to move the view up.",
                          moveViewUp );

  keyCombo->SetFromString( "Down" );
  PreferencesManager::StringPrefValue moveViewDown( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyMoveViewDown",
                          "Key to move the view down.",
                          moveViewDown );

  keyCombo->SetFromString( "PageUp" );
  PreferencesManager::StringPrefValue moveViewIn( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyMoveViewIn",
                          "Key to move the view in in plane.",
                          moveViewIn );

  keyCombo->SetFromString( "PageDown" );
  PreferencesManager::StringPrefValue moveViewOut( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyMoveViewOut",
                          "Key to move the view out in plane.",
                          moveViewOut );

  keyCombo->SetFromString( "Equal" );
  PreferencesManager::StringPrefValue zoomViewIn( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyZoomViewIn",
                          "Key to zoom the view in in plane.",
                          zoomViewIn );

  keyCombo->SetFromString( "Minus" );
  PreferencesManager::StringPrefValue zoomViewOut( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyZoomViewOut",
                          "Key to zoom the view out in plane.",
                          zoomViewOut );

  keyCombo->SetFromString( "n" );
  PreferencesManager::StringPrefValue toolNavigate( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyToolNavigate",
                          "Shortcut for the Navigate tool.",
                          toolNavigate );

  keyCombo->SetFromString( "p" );
  PreferencesManager::StringPrefValue toolPlane( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyToolPlane",
                          "Shortcut for the Plane tool.",
                          toolPlane );

  keyCombo->SetFromString( "m" );
  PreferencesManager::StringPrefValue toolMarker( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyToolMarker",
                          "Shortcut for the Marker tool.",
                          toolMarker );

  keyCombo->SetFromString( "e" );
  PreferencesManager::StringPrefValue toolEditVoxel( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyToolEditVoxel",
                          "Shortcut for the Edit Voxels tool.",
                          toolEditVoxel );

  keyCombo->SetFromString( "r" );
  PreferencesManager::StringPrefValue toolEditROI( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyToolEditROI",
                          "Shortcut for the Edit ROI tool.",
                          toolEditROI );

  keyCombo->SetFromString( "s" );
  PreferencesManager::StringPrefValue toolStraightPath( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyToolStraightPath",
                          "Shortcut for the Straigh Path tool.",
                          toolStraightPath );

  keyCombo->SetFromString( "g" );
  PreferencesManager::StringPrefValue toolEdgePath( keyCombo->ToString() );
  prefsMgr.RegisterValue( "KeyToolEdgePath",
                          "Shortcut for the Edge Path tool.",
                          toolEdgePath );

  PreferencesManager::IntPrefValue showFPS( false );
  prefsMgr.RegisterValue( "ShowFPS",
                          "Show the FPS (frames per second) for draws.",
                          showFPS );

  PreferencesManager::StringPrefValue selectedTool( "navigation" );
  prefsMgr.RegisterValue( "SelectedTool",
                          "The selected tool.",
                          selectedTool );

  PreferencesManager::StringPrefValue userStructureList( "{}" );
  prefsMgr.RegisterValue( "UserStructureList",
                          "The user's preferred structure list for editing.",
                          userStructureList );

  PreferencesManager::IntPrefValue lockOnCursor( false );
  prefsMgr.RegisterValue( "LockOnCursor",
                          "Lock On Cursor setting.",
                          lockOnCursor );

  delete keyCombo;
}

void
ScubaGlobalPreferences::WritePreferences () {

  // Write out the file.
  PreferencesManager& prefsMgr = PreferencesManager::GetManager();
  prefsMgr.SaveFile();
}
