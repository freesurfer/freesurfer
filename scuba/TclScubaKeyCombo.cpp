#include "TclScubaKeyCombo.h"

using namespace std;

bool TclScubaKeyComboStaticTclListener::mbAddedTclCommands = false;

TclScubaKeyCombo::TclScubaKeyCombo () :
  ScubaKeyCombo() {
}

void
TclScubaKeyCombo::SetFromString ( string isKey ) {

  ScubaKeyCombo::SetFromString( isKey );
  
  // We are initialied with a key code already; even if it didn't find
  // the right tk name, it definitely matched a last letter. So this
  // just gives the potential to override what we've already done, and
  // if not, goes with what we got.
  if(isKey.find("bracketright")!=string::npos) mKeyCode =Key_BracketRight;
  else if(isKey.find("bracketleft")!=string::npos) mKeyCode = Key_BracketLeft;
  else if(isKey.find("asciicirum")!=string::npos) mKeyCode = Key_AsciiCircum;
  else if(isKey.find("underscore")!=string::npos) mKeyCode = Key_Underscore;
  else if(isKey.find("parenright")!=string::npos) mKeyCode = Key_ParenRight;
  else if(isKey.find("numbersign")!=string::npos) mKeyCode = Key_NumberSign;
  else if(isKey.find("apostrophe")!=string::npos) mKeyCode = Key_Apostrophe;
  else if(isKey.find("braceright")!=string::npos) mKeyCode = Key_BraceRight;
  else if(isKey.find("asciitilde")!=string::npos) mKeyCode = Key_AsciiTilde;
  else if(isKey.find("backslash")!=string::npos) mKeyCode = Key_Backslash;
  else if(isKey.find("ampersand")!=string::npos) mKeyCode = Key_Ampersand;
  else if(isKey.find("parenleft")!=string::npos) mKeyCode = Key_ParenLeft;
  else if(isKey.find("semicolon")!=string::npos) mKeyCode = Key_Semicolon;
  else if(isKey.find("braceleft")!=string::npos) mKeyCode = Key_BraceLeft;
  else if(isKey.find("backslash")!=string::npos) mKeyCode = Key_Backslash;
  else if(isKey.find("asterisk")!=string::npos) mKeyCode = Key_Asterisk;
  else if(isKey.find("quotedbl")!=string::npos) mKeyCode = Key_QuoteDbl;
  else if(isKey.find("percent")!=string::npos) mKeyCode = Key_Percent;
  else if(isKey.find("greater")!=string::npos) mKeyCode = Key_Greater;
  else if(isKey.find("exclam")!=string::npos) mKeyCode = Key_Exclam;
  else if(isKey.find("dollar")!=string::npos) mKeyCode = Key_Dollar;
  else if(isKey.find("period")!=string::npos) mKeyCode = Key_Period;
  else if(isKey.find("Prior")!=string::npos) mKeyCode = Key_PageUp;
  else if(isKey.find("minus")!=string::npos) mKeyCode = Key_Minus;
  else if(isKey.find("equal")!=string::npos) mKeyCode = Key_Equal;
  else if(isKey.find("space")!=string::npos) mKeyCode = Key_Space;
  else if(isKey.find("comma")!=string::npos) mKeyCode = Key_Comma;
  else if(isKey.find("minus")!=string::npos) mKeyCode = Key_Minus;
  else if(isKey.find("slash")!=string::npos) mKeyCode = Key_Slash;
  else if(isKey.find("colon")!=string::npos) mKeyCode = Key_Colon;
  else if(isKey.find("equal")!=string::npos) mKeyCode = Key_Equal;
  else if(isKey.find("grave")!=string::npos) mKeyCode = Key_QuoteLeft;
  else if(isKey.find("less")!=string::npos) mKeyCode = Key_Less;
  else if(isKey.find("Next")!=string::npos) mKeyCode = Key_PageDown;
  else if(isKey.find("plus")!=string::npos) mKeyCode = Key_Plus;
  else if(isKey.find("bar")!=string::npos) mKeyCode = Key_Bar;
  else if(isKey.find("at")!=string::npos) mKeyCode = Key_At;
}


TclScubaKeyComboStaticTclListener&
TclScubaKeyComboStaticTclListener::GetListener () {

  static TclScubaKeyComboStaticTclListener sListener;

  if( !mbAddedTclCommands ) {

    mbAddedTclCommands = true;
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( sListener,
			   "ConvertTkInputStringToScubaKeyComboString", 1,
			   "input", "Takes an input string from Tk key text "
			   "and returns a ScubaKeyCombo key combo string." );
  }

  return sListener;
}

TclCommandManager::TclCommandResult
TclScubaKeyComboStaticTclListener::DoListenToTclCommand ( char* isCommand, 
							  int iArgc,
							  char** iasArgv ) {
  // ConvertTkInputStringToScubaKeyComboString <input string>
  if( 0 == strcmp( isCommand, "ConvertTkInputStringToScubaKeyComboString" ) ) {

    string sKey = iasArgv[1];

    TclScubaKeyCombo key;
    key.SetFromString( sKey );
    stringstream ssReturnValues;
    ssReturnValues << "\"" << key.ToString() << "\"";
    sReturnValues = ssReturnValues.str();
    sReturnFormat = "s";
    return ok;
  }

  return ok;
}
