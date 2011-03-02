/**
 * @file  ScubaKeyCombo.cpp
 * @brief Encapsulates keypress combinations
 *
 * This is a class to encapsulate the idea of key combinations, or
 * keypresses to trigger events. These are used in a few places:
 *
 *   The main interface widget (e.g. ToglFrame) transforms keypress
 *   events into ScubaKeyCombo subclasses from their native widget
 *   type.

 *   Other widgets that check for key combos in event response
 *   functions can simply match input ScubaKeyCombos against
 *   ScubaKeyCombos stored in preferences. For example, ScubaView
 *   grabs the string pref for the MoveViewIn key, makes a
 *   ScubaKeyCombo out of it, and checks it against the current key
 *   from the InputState.
 *
 *   Preference dialog boxes read input from the user, create a
 *   ScubaKeyCombo from it, and use the ToString function to generate
 *   the string version to be stored in ScubaGlobalPreferences.
 *
 *   Each widget type should have a subclass ScubaKeyCombo
 *   (e.g. TclScubaKeyCombo) that knows how to parse input strings
 *   from that widget framework. A corresponding ScubaKeyComboFactory
 *   should be made to make those kinds of ScubaKeyCombo
 *   subclasses. The main function, knowing what kind of widget
 *   framework is being used, calls SetFactory with that subclass
 *   factory. Then, all other classes use ScubaKeyCombo::MakeKeyCombo
 *   to make their framework-specific ScubaKeyCombos.
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:37 $
 *    $Revision: 1.10 $
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
#include "ScubaKeyCombo.h"

using namespace std;

ScubaKeyComboFactory* ScubaKeyCombo::mFactory = NULL;

ScubaKeyCombo::ScubaKeyCombo () :
    mKeyCode(-1), mbShift(false), mbAlt(false), mbMeta(false), mbControl(false) {}

void
ScubaKeyCombo::SetFromString ( string const& isKey ) {

  // We're getting a string that is hopefully from a
  // ScubaKeyCombo::ToString originally, so we just look for the
  // strings we output normally in that function.

  // Look for the modifier strings first.
  if ( isKey.find ("Ctrl")!=string::npos )
    mbControl = true;
  else
    mbControl = false;
  if ( isKey.find ("Shift")!=string::npos )
    mbShift = true;
  else
    mbShift = false;
  if ( isKey.find ("Alt")!=string::npos )
    mbAlt = true;
  else
    mbAlt = false;
  if ( isKey.find ("Meta")!=string::npos )
    mbMeta = true;
  else
    mbMeta = false;

  struct keyStringCodePair {
    string sKey;
    int code;
  };
  keyStringCodePair aKeys[] = {
                                {"BracketRight", Key_BracketRight
                                },
                                {"BracketLeft", Key_BracketLeft},
                                {"AsciiCircum", Key_AsciiCircum},
                                {"Direction_L", Key_Direction_L},
                                {"Direction_R", Key_Direction_R},
                                {"NumberSign", Key_NumberSign},
                                {"Apostrophe", Key_Apostrophe},
                                {"ParenRight", Key_ParenRight},
                                {"BraceRight", Key_BraceRight},
                                {"AsciiTilde", Key_AsciiTilde},
                                {"ScrollLock", Key_ScrollLock},
                                {"Underscore", Key_Underscore},
                                {"ParenLeft", Key_ParenLeft},
                                {"Ampersand", Key_Ampersand},
                                {"Semicolon", Key_Semicolon},
                                {"Backspace", Key_Backspace},
                                {"QuoteLeft", Key_QuoteLeft},
                                {"BraceLeft", Key_BraceLeft},
                                {"Backslash", Key_Backslash},
                                {"PageDown", Key_Next},
                                {"CapsLock", Key_CapsLock},
                                {"Question", Key_Question},
                                {"Asterisk", Key_Asterisk},
                                {"QuoteDbl", Key_QuoteDbl},
                                {"Backtab", Key_Backtab},
                                {"NumLock", Key_NumLock},
                                {"Super_L", Key_Super_L},
                                {"Super_R", Key_Super_R},
                                {"Hyper_L", Key_Hyper_L},
                                {"Hyper_R", Key_Hyper_R},
                                {"Greater", Key_Greater},
                                {"Percent", Key_Percent},
                                {"Dollar", Key_Dollar},
                                {"Period", Key_Period},
                                {"Exclam", Key_Exclam},
                                {"Return", Key_Return},
                                {"Escape", Key_Escape},
                                {"Insert", Key_Insert},
                                {"Delete", Key_Delete},
                                {"SysReq", Key_SysReq},
                                {"PageUp", Key_Prior},
                                {"Comma", Key_Comma},
                                {"Minus", Key_Minus},
                                {"minus", Key_Minus},
                                {"Slash", Key_Slash},
                                {"Colon", Key_Colon},
                                {"Enter", Key_Enter},
                                {"Pause", Key_Pause},
                                {"Print", Key_Print},
                                {"Clear", Key_Clear},
                                {"Right", Key_Right},
                                {"Equal", Key_Equal},
                                {"equal", Key_Equal},
                                {"Space", Key_Space},
                                {"Plus", Key_Plus},
                                {"Less", Key_Less},
                                {"Help", Key_Help},
                                {"Home", Key_Home},
                                {"Left", Key_Left},
                                {"Down", Key_Down},
                                {"Menu", Key_Menu},
                                {"End", Key_End},
                                {"Tab", Key_Tab},
                                {"Bar", Key_Bar},
                                {"F10", Key_F10},
                                {"F11", Key_F11},
                                {"F12", Key_F12},
                                {"F13", Key_F13},
                                {"F14", Key_F14},
                                {"F15", Key_F15},
                                {"F16", Key_F16},
                                {"F17", Key_F17},
                                {"F18", Key_F18},
                                {"F19", Key_F19},
                                {"F20", Key_F20},
                                {"F21", Key_F21},
                                {"F22", Key_F22},
                                {"F23", Key_F23},
                                {"F24", Key_F24},
                                {"F25", Key_F25},
                                {"F26", Key_F26},
                                {"F27", Key_F27},
                                {"F28", Key_F28},
                                {"F29", Key_F29},
                                {"F30", Key_F30},
                                {"F31", Key_F31},
                                {"F32", Key_F32},
                                {"F33", Key_F33},
                                {"F34", Key_F34},
                                {"F35", Key_F35},
                                {"At", Key_At},
                                {"F1", Key_F1},
                                {"F2", Key_F2},
                                {"F3", Key_F3},
                                {"F4", Key_F4},
                                {"F5", Key_F5},
                                {"F6", Key_F6},
                                {"F7", Key_F7},
                                {"F8", Key_F8},
                                {"F9", Key_F9},
                                {"Up", Key_Up},
                                {"0", Key_0},
                                {"1", Key_1},
                                {"2", Key_2},
                                {"3", Key_3},
                                {"4", Key_4},
                                {"5", Key_5},
                                {"6", Key_6},
                                {"7", Key_7},
                                {"8", Key_8},
                                {"9", Key_9},
                                {"A", Key_A},
                                {"B", Key_B},
                                {"C", Key_C},
                                {"D", Key_D},
                                {"E", Key_E},
                                {"F", Key_F},
                                {"G", Key_G},
                                {"H", Key_H},
                                {"I", Key_I},
                                {"J", Key_J},
                                {"K", Key_K},
                                {"L", Key_L},
                                {"M", Key_M},
                                {"N", Key_N},
                                {"O", Key_O},
                                {"P", Key_P},
                                {"Q", Key_Q},
                                {"R", Key_R},
                                {"S", Key_S},
                                {"T", Key_T},
                                {"U", Key_U},
                                {"V", Key_V},
                                {"W", Key_W},
                                {"X", Key_X},
                                {"Y", Key_Y},
                                {"Z", Key_Z},
                                {"a", Key_A},
                                {"b", Key_B},
                                {"c", Key_C},
                                {"d", Key_D},
                                {"e", Key_E},
                                {"f", Key_F},
                                {"g", Key_G},
                                {"h", Key_H},
                                {"i", Key_I},
                                {"j", Key_J},
                                {"k", Key_K},
                                {"l", Key_L},
                                {"m", Key_M},
                                {"n", Key_N},
                                {"o", Key_O},
                                {"p", Key_P},
                                {"q", Key_Q},
                                {"r", Key_R},
                                {"s", Key_S},
                                {"t", Key_T},
                                {"u", Key_U},
                                {"v", Key_V},
                                {"w", Key_W},
                                {"x", Key_X},
                                {"y", Key_Y},
                                {"z", Key_Z}
                              };

  // For each key string, try a reverse find on the input string we
  // got. If we found it, check that the position we got is correct,
  // and if we got a complete word (by checking that the length of the
  // match is the same as the length of the string, or that the char
  // before the match is a space).
  for ( int nKey = 0; nKey < 163; nKey++ ) {
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

  //  DebugOutput( << "Couldn't find a match for key " << isKey );
  mKeyCode = Key_unknown;
}

string
ScubaKeyCombo::ToString () const {

  string sModifers = "";
  if ( mbControl ) sModifers += "Ctrl ";
  if ( mbShift ) sModifers += "Shift ";
  if ( mbAlt ) sModifers += "Alt ";
  if ( mbMeta ) sModifers += "Meta ";

  string sKey;
  switch ( mKeyCode ) {
  case Key_Escape:
    sKey = "Escape";
    break;
  case Key_Tab:
    sKey = "Tab";
    break;
  case Key_Backtab:
    sKey = "Backtab";
    break;
  case Key_Backspace:
    sKey = "Backspace";
    break;
  case Key_Return:
    sKey = "Return";
    break;
  case Key_Enter:
    sKey = "Enter";
    break;
  case Key_Insert:
    sKey = "Insert";
    break;
  case Key_Delete:
    sKey = "Delete";
    break;
  case Key_Pause:
    sKey = "Pause";
    break;
  case Key_Print:
    sKey = "Print";
    break;
  case Key_SysReq:
    sKey = "SysReq";
    break;
  case Key_Clear:
    sKey = "Clear";
    break;
  case Key_Home:
    sKey = "Home";
    break;
  case Key_End:
    sKey = "End";
    break;
  case Key_Left:
    sKey = "Left";
    break;
  case Key_Up:
    sKey = "Up";
    break;
  case Key_Right:
    sKey = "Right";
    break;
  case Key_Down:
    sKey = "Down";
    break;
  case Key_Prior:
    sKey = "PageUp";
    break;
  case Key_Next:
    sKey = "PageDown";
    break;
  case Key_Shift:
    sKey = "Shift";
    break;
  case Key_Control:
    sKey = "Control";
    break;
  case Key_Meta:
    sKey = "Meta";
    break;
  case Key_Alt:
    sKey = "Alt";
    break;
  case Key_CapsLock:
    sKey = "CapsLock";
    break;
  case Key_NumLock:
    sKey = "NumLock";
    break;
  case Key_ScrollLock:
    sKey = "ScrollLock";
    break;
  case Key_F1:
    sKey = "F1";
    break;
  case Key_F2:
    sKey = "F2";
    break;
  case Key_F3:
    sKey = "F3";
    break;
  case Key_F4:
    sKey = "F4";
    break;
  case Key_F5:
    sKey = "F5";
    break;
  case Key_F6:
    sKey = "F6";
    break;
  case Key_F7:
    sKey = "F7";
    break;
  case Key_F8:
    sKey = "F8";
    break;
  case Key_F9:
    sKey = "F9";
    break;
  case Key_F10:
    sKey = "F10";
    break;
  case Key_F11:
    sKey = "F11";
    break;
  case Key_F12:
    sKey = "F12";
    break;
  case Key_F13:
    sKey = "F13";
    break;
  case Key_F14:
    sKey = "F14";
    break;
  case Key_F15:
    sKey = "F15";
    break;
  case Key_F16:
    sKey = "F16";
    break;
  case Key_F17:
    sKey = "F17";
    break;
  case Key_F18:
    sKey = "F18";
    break;
  case Key_F19:
    sKey = "F19";
    break;
  case Key_F20:
    sKey = "F20";
    break;
  case Key_F21:
    sKey = "F21";
    break;
  case Key_F22:
    sKey = "F22";
    break;
  case Key_F23:
    sKey = "F23";
    break;
  case Key_F24:
    sKey = "F24";
    break;
  case Key_F25:
    sKey = "F25";
    break;
  case Key_F26:
    sKey = "F26";
    break;
  case Key_F27:
    sKey = "F27";
    break;
  case Key_F28:
    sKey = "F28";
    break;
  case Key_F29:
    sKey = "F29";
    break;
  case Key_F30:
    sKey = "F30";
    break;
  case Key_F31:
    sKey = "F31";
    break;
  case Key_F32:
    sKey = "F32";
    break;
  case Key_F33:
    sKey = "F33";
    break;
  case Key_F34:
    sKey = "F34";
    break;
  case Key_F35:
    sKey = "F35";
    break;
  case Key_Super_L:
    sKey = "Super_L";
    break;
  case Key_Super_R:
    sKey = "Super_R";
    break;
  case Key_Menu:
    sKey = "Menu";
    break;
  case Key_Hyper_L:
    sKey = "Hyper_L";
    break;
  case Key_Hyper_R:
    sKey = "Hyper_R";
    break;
  case Key_Help:
    sKey = "Help";
    break;
  case Key_Direction_L:
    sKey = "Direction_L";
    break;
  case Key_Direction_R:
    sKey = "Direction_R";
    break;
  case Key_Space:
    sKey = "Space";
    break;
  case Key_Exclam:
    sKey = "Exclam";
    break;
  case Key_QuoteDbl:
    sKey = "QuoteDbl";
    break;
  case Key_NumberSign:
    sKey = "NumberSign";
    break;
  case Key_Dollar:
    sKey = "Dollar";
    break;
  case Key_Percent:
    sKey = "Percent";
    break;
  case Key_Ampersand:
    sKey = "Ampersand";
    break;
  case Key_Apostrophe:
    sKey = "Apostrophe";
    break;
  case Key_ParenLeft:
    sKey = "ParenLeft";
    break;
  case Key_ParenRight:
    sKey = "ParenRight";
    break;
  case Key_Asterisk:
    sKey = "Asterisk";
    break;
  case Key_Plus:
    sKey = "Plus";
    break;
  case Key_Comma:
    sKey = "Comma";
    break;
  case Key_Minus:
    sKey = "Minus";
    break;
  case Key_Period:
    sKey = "Period";
    break;
  case Key_Slash:
    sKey = "Slash";
    break;
  case Key_0:
    sKey = "0";
    break;
  case Key_1:
    sKey = "1";
    break;
  case Key_2:
    sKey = "2";
    break;
  case Key_3:
    sKey = "3";
    break;
  case Key_4:
    sKey = "4";
    break;
  case Key_5:
    sKey = "5";
    break;
  case Key_6:
    sKey = "6";
    break;
  case Key_7:
    sKey = "7";
    break;
  case Key_8:
    sKey = "8";
    break;
  case Key_9:
    sKey = "9";
    break;
  case Key_Colon:
    sKey = "Colon";
    break;
  case Key_Semicolon:
    sKey = "Semicolon";
    break;
  case Key_Less:
    sKey = "Less";
    break;
  case Key_Equal:
    sKey = "Equal";
    break;
  case Key_Greater:
    sKey = "Greater";
    break;
  case Key_Question:
    sKey = "Question";
    break;
  case Key_At:
    sKey = "At";
    break;
  case Key_A:
    sKey = "a";
    break;
  case Key_B:
    sKey = "b";
    break;
  case Key_C:
    sKey = "c";
    break;
  case Key_D:
    sKey = "d";
    break;
  case Key_E:
    sKey = "e";
    break;
  case Key_F:
    sKey = "f";
    break;
  case Key_G:
    sKey = "g";
    break;
  case Key_H:
    sKey = "h";
    break;
  case Key_I:
    sKey = "i";
    break;
  case Key_J:
    sKey = "j";
    break;
  case Key_K:
    sKey = "k";
    break;
  case Key_L:
    sKey = "l";
    break;
  case Key_M:
    sKey = "m";
    break;
  case Key_N:
    sKey = "n";
    break;
  case Key_O:
    sKey = "o";
    break;
  case Key_P:
    sKey = "p";
    break;
  case Key_Q:
    sKey = "q";
    break;
  case Key_R:
    sKey = "r";
    break;
  case Key_S:
    sKey = "s";
    break;
  case Key_T:
    sKey = "t";
    break;
  case Key_U:
    sKey = "u";
    break;
  case Key_V:
    sKey = "v";
    break;
  case Key_W:
    sKey = "w";
    break;
  case Key_X:
    sKey = "x";
    break;
  case Key_Y:
    sKey = "y";
    break;
  case Key_Z:
    sKey = "z";
    break;
  case Key_BracketLeft:
    sKey = "BracketLeft";
    break;
  case Key_Backslash:
    sKey = "Backslash";
    break;
  case Key_BracketRight:
    sKey = "BracketRight";
    break;
  case Key_AsciiCircum:
    sKey = "AsciiCircum";
    break;
  case Key_Underscore:
    sKey = "Underscore";
    break;
  case Key_QuoteLeft:
    sKey = "QuoteLeft";
    break;
  case Key_BraceLeft:
    sKey = "BraceLeft";
    break;
  case Key_Bar:
    sKey = "Bar";
    break;
  case Key_BraceRight:
    sKey = "BraceRight";
    break;
  case Key_AsciiTilde:
    sKey = "AsciiTilde";
    break;
  case Key_unknown:
    sKey = "unknown";
    break;
  }

  return sModifers + sKey;
}

void
ScubaKeyCombo::CopyFrom ( ScubaKeyCombo const& iKey ) {

  mKeyCode = iKey.GetKeyCode();
  mbShift = iKey.IsShiftKeyDown();
  mbAlt = iKey.IsAltKeyDown();
  mbMeta = iKey.IsMetaKeyDown();
  mbControl = iKey.IsControlKeyDown();
}

void 
ScubaKeyCombo::CopyFrom ( ScubaKeyCombo const* iKey ) {

  CopyFrom( *iKey );
}

int
ScubaKeyCombo::GetKeyCode () const {
  return mKeyCode;
}

bool
ScubaKeyCombo::IsShiftKeyDown () const {
  return mbShift;
}

bool
ScubaKeyCombo::IsAltKeyDown () const {
  return mbAlt;
}

bool
ScubaKeyCombo::IsMetaKeyDown () const {
  return mbMeta;
}

bool
ScubaKeyCombo::IsControlKeyDown () const {
  return mbControl;
}

void
ScubaKeyCombo::SetKeyCode ( int iKey ) {
  mKeyCode = iKey;
}

void
ScubaKeyCombo::SetShiftKeyDown ( bool isDown ) {
  mbShift = isDown;
}

void
ScubaKeyCombo::SetAltKeyDown ( bool isDown ) {
  mbAlt = isDown;
}

void
ScubaKeyCombo::SetMetaKeyDown ( bool isDown ) {
  mbMeta = isDown;
}

void
ScubaKeyCombo::SetControlKeyDown ( bool isDown ) {
  mbControl = isDown;
}

bool
ScubaKeyCombo::IsSameAs ( ScubaKeyCombo const& iCombo ) const {

  return ( GetKeyCode() == iCombo.GetKeyCode() &&
           IsShiftKeyDown() == iCombo.IsShiftKeyDown() &&
           IsAltKeyDown() == iCombo.IsAltKeyDown() &&
           IsMetaKeyDown() == iCombo.IsMetaKeyDown() &&
           IsControlKeyDown() == iCombo.IsControlKeyDown() );
}

bool
ScubaKeyCombo::IsSameAs ( ScubaKeyCombo const* iCombo ) const {

  return IsSameAs( *iCombo );
}

void
ScubaKeyCombo::SetFactory ( ScubaKeyComboFactory* iFactory ) {

  mFactory = iFactory;
}

ScubaKeyCombo*
ScubaKeyCombo::MakeKeyCombo () {

  if ( mFactory ) {
    return mFactory->MakeKeyCombo();
  } else {
    return new ScubaKeyCombo();
  }
}

ostream&
operator <<  ( ostream& os, ScubaKeyCombo const& iKey ) {

  os.setf(ios::hex,ios::basefield);
  os << "Key " << iKey.ToString() << " (code 0x" << iKey.GetKeyCode()
  << ")";
  os.unsetf(ios::hex);
  return os;
}
