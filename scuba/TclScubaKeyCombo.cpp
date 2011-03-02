/**
 * @file  TclScubaKeyCombo.cpp
 * @brief Tcl specific ScubaKeyCombo
 *
 * This overrides the SetFromString function so we can parse a Tk key
 * string.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.11 $
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


#include <iostream>
#include "TclScubaKeyCombo.h"
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h> // strcmp
#ifdef __cplusplus
}
#endif

using namespace std;

bool TclScubaKeyComboStaticTclListener::mbAddedTclCommands = false;

TclScubaKeyCombo::TclScubaKeyCombo () :
  ScubaKeyCombo() {}

void
TclScubaKeyCombo::SetFromString ( string const& isKey ) {

  ScubaKeyCombo::SetFromString( isKey );

  // We are initialied with a key code already; even if it didn't find
  // the right tk name, it definitely matched a last letter. So this
  // just gives the potential to override what we've already done, and
  // if not, goes with what we got.
  struct keyStringCodePair {
    string sKey;
    int code;
  };
  keyStringCodePair aKeys[] = {
                                {"bracketright", Key_BracketRight },
                                { "bracketleft", Key_BracketLeft },
                                { "asciicirum", Key_AsciiCircum },
                                { "underscore", Key_Underscore },
                                { "parenright", Key_ParenRight },
                                { "numbersign", Key_NumberSign },
                                { "apostrophe", Key_Apostrophe },
                                { "braceright", Key_BraceRight },
                                { "asciitilde", Key_AsciiTilde },
                                { "backslash", Key_Backslash },
                                { "ampersand", Key_Ampersand },
                                { "parenleft", Key_ParenLeft },
                                { "semicolon", Key_Semicolon },
                                { "braceleft", Key_BraceLeft },
                                { "backslash", Key_Backslash },
                                { "asterisk", Key_Asterisk },
                                { "quotedbl", Key_QuoteDbl },
                                { "percent", Key_Percent },
                                { "greater", Key_Greater },
                                { "exclam", Key_Exclam },
                                { "dollar", Key_Dollar },
                                { "period", Key_Period },
                                { "Prior", Key_Prior },
                                { "minus", Key_Minus },
                                { "equal", Key_Equal },
                                { "space", Key_Space },
                                { "comma", Key_Comma },
                                { "minus", Key_Minus },
                                { "slash", Key_Slash },
                                { "colon", Key_Colon },
                                { "equal", Key_Equal },
                                { "grave", Key_QuoteLeft },
                                { "less", Key_Less },
                                { "Next", Key_PageDown },
                                { "plus", Key_Plus },
                                { "bar", Key_Bar },
                                { "at", Key_At } };

  // For each key string, try a reverse find on the input string we
  // got. If we found it, check that the position we got is correct,
  // and if we got a complete word (by checking that the length of the
  // match is the same as the length of the string, or that the char
  // before the match is a space).
  for ( int nKey = 0; nKey < 37; nKey++ ) {
    size_t pos;
    pos = isKey.rfind( aKeys[nKey].sKey );
    if ( pos != string::npos &&
         pos == isKey.length () - aKeys[nKey].sKey.length() &&
         (isKey.length() == aKeys[nKey].sKey.length() ||
          isKey[pos-1] == ' ') ) {
      mKeyCode = aKeys[nKey].code;
      return;
    }
  }
}


TclScubaKeyComboStaticTclListener&
TclScubaKeyComboStaticTclListener::GetListener () {

  static TclScubaKeyComboStaticTclListener sListener;

  if ( !mbAddedTclCommands ) {

    mbAddedTclCommands = true;
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( sListener,
                           "ConvertTkInputStringToScubaKeyComboString", 5,
                           "input shift meta alt control", "Takes an input "
                           "string from Tk key text and boolean values for "
                           "modifiers and returns a ScubaKeyCombo key "
                           "combo string." );
  }

  return sListener;
}

TclCommandManager::TclCommandResult
TclScubaKeyComboStaticTclListener::DoListenToTclCommand ( char* isCommand,
							  int iArgc,
							  char** iasArgv ) {

  // ConvertTkInputStringToScubaKeyComboString <input string> <shift> <meta> <alt> <control>
  if ( 0 == strcmp( isCommand, "ConvertTkInputStringToScubaKeyComboString" ) ) {

    string sKey = iasArgv[1];
    bool bShift = TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
    bool bMeta = TclCommandManager::ConvertArgumentToBoolean( iasArgv[3] );
    bool bAlt = TclCommandManager::ConvertArgumentToBoolean( iasArgv[4] );
    bool bControl = TclCommandManager::ConvertArgumentToBoolean( iasArgv[5] );

    TclScubaKeyCombo key;
    key.SetFromString( sKey );
    key.SetShiftKeyDown( bShift );
    key.SetMetaKeyDown( bMeta );
    key.SetAltKeyDown( bAlt );
    key.SetControlKeyDown( bControl );
    stringstream ssReturnValues;
    ssReturnValues << "\"" << key.ToString() << "\"";
    sReturnValues = ssReturnValues.str();
    sReturnFormat = "s";
    return ok;
  }

  return ok;
}
