#include <iostream>
#include "ScubaKeyCombo.h"

using namespace std;

ScubaKeyCombo::ScubaKeyCombo () :
  mKeyCode(-1), mbShift(false), mbAlt(false), mbMeta(false), mbControl(false) {

}

ScubaKeyCombo::ScubaKeyCombo ( string isKey ) :
  mKeyCode(-1), mbShift(false), mbAlt(false), mbMeta(false), mbControl(false) {

  // We're getting a string that is hopefully from a
  // ScubaKeyCombo::ToString originally, so we just look for the
  // strings we output normally in that function. 

  // Look for the modifier strings first.
  if( isKey.find ("Ctrl")!=string::npos ) mbControl = true;
  if( isKey.find ("Shift")!=string::npos ) mbShift = true;
  if( isKey.find ("Alt")!=string::npos ) mbAlt = true;
  if( isKey.find ("Meta")!=string::npos ) mbMeta = true;

  // The idea here is to search for larger strings first to avoid
  // partial matches. e.g. if we searched for Up before PageUp, the
  // key could be PageUp but we would match Up first.
  if(isKey.find("BracketRight")!=string::npos) mKeyCode = Key_BracketRight;
  else if(isKey.find("BracketLeft")!=string::npos) mKeyCode = Key_BracketLeft;
  else if(isKey.find("AsciiCircum")!=string::npos) mKeyCode = Key_AsciiCircum;
  else if(isKey.find("Direction_L")!=string::npos) mKeyCode = Key_Direction_L;
  else if(isKey.find("Direction_R")!=string::npos) mKeyCode = Key_Direction_R;
  else if(isKey.find("NumberSign")!=string::npos) mKeyCode = Key_NumberSign;
  else if(isKey.find("Apostrophe")!=string::npos) mKeyCode = Key_Apostrophe;
  else if(isKey.find("ParenRight")!=string::npos) mKeyCode = Key_ParenRight;
  else if(isKey.find("BraceRight")!=string::npos) mKeyCode = Key_BraceRight;
  else if(isKey.find("AsciiTilde")!=string::npos) mKeyCode = Key_AsciiTilde;
  else if(isKey.find("ScrollLock")!=string::npos) mKeyCode = Key_ScrollLock;
  else if(isKey.find("Underscore")!=string::npos) mKeyCode = Key_Underscore;
  else if(isKey.find("ParenLeft")!=string::npos) mKeyCode = Key_ParenLeft;
  else if(isKey.find("Ampersand")!=string::npos) mKeyCode = Key_Ampersand;
  else if(isKey.find("Semicolon")!=string::npos) mKeyCode = Key_Semicolon;
  else if(isKey.find("Backspace")!=string::npos) mKeyCode = Key_Backspace;
  else if(isKey.find("QuoteLeft")!=string::npos) mKeyCode = Key_QuoteLeft;
  else if(isKey.find("BraceLeft")!=string::npos) mKeyCode = Key_BraceLeft;
  else if(isKey.find("Backslash")!=string::npos) mKeyCode = Key_Backslash;
  else if(isKey.find("PageDown")!=string::npos) mKeyCode = Key_Next;
  else if(isKey.find("CapsLock")!=string::npos) mKeyCode = Key_CapsLock;
  else if(isKey.find("Question")!=string::npos) mKeyCode = Key_Question;
  else if(isKey.find("Asterisk")!=string::npos) mKeyCode = Key_Asterisk;
  else if(isKey.find("QuoteDbl")!=string::npos) mKeyCode = Key_QuoteDbl;
  else if(isKey.find("Backtab")!=string::npos) mKeyCode = Key_Backtab;
  else if(isKey.find("NumLock")!=string::npos) mKeyCode = Key_NumLock;
  else if(isKey.find("Super_L")!=string::npos) mKeyCode = Key_Super_L;
  else if(isKey.find("Super_R")!=string::npos) mKeyCode = Key_Super_R;
  else if(isKey.find("Hyper_L")!=string::npos) mKeyCode = Key_Hyper_L;
  else if(isKey.find("Hyper_R")!=string::npos) mKeyCode = Key_Hyper_R;
  else if(isKey.find("Greater")!=string::npos) mKeyCode = Key_Greater;
  else if(isKey.find("Percent")!=string::npos) mKeyCode = Key_Percent;
  else if(isKey.find("Dollar")!=string::npos) mKeyCode = Key_Dollar;
  else if(isKey.find("Period")!=string::npos) mKeyCode = Key_Period;
  else if(isKey.find("Exclam")!=string::npos) mKeyCode = Key_Exclam;
  else if(isKey.find("Return")!=string::npos) mKeyCode = Key_Return;
  else if(isKey.find("Escape")!=string::npos) mKeyCode = Key_Escape;
  else if(isKey.find("Insert")!=string::npos) mKeyCode = Key_Insert;
  else if(isKey.find("Delete")!=string::npos) mKeyCode = Key_Delete;
  else if(isKey.find("SysReq")!=string::npos) mKeyCode = Key_SysReq;
  else if(isKey.find("PageUp")!=string::npos) mKeyCode = Key_Prior;
  else if(isKey.find("Comma")!=string::npos) mKeyCode = Key_Comma;
  else if(isKey.find("Minus")!=string::npos) mKeyCode = Key_Minus; 
  else if(isKey.find("Slash")!=string::npos) mKeyCode = Key_Slash;
  else if(isKey.find("Colon")!=string::npos) mKeyCode = Key_Colon;
  else if(isKey.find("Enter")!=string::npos) mKeyCode = Key_Enter;
  else if(isKey.find("Pause")!=string::npos) mKeyCode = Key_Pause;
  else if(isKey.find("Print")!=string::npos) mKeyCode = Key_Print;
  else if(isKey.find("Clear")!=string::npos) mKeyCode = Key_Clear;
  else if(isKey.find("Right")!=string::npos) mKeyCode = Key_Right;
  else if(isKey.find("Equal")!=string::npos) mKeyCode = Key_Equal;
  else if(isKey.find("Space")!=string::npos) mKeyCode = Key_Space;
  else if(isKey.find("Plus")!=string::npos) mKeyCode = Key_Plus;
  else if(isKey.find("Less")!=string::npos) mKeyCode = Key_Less;
  else if(isKey.find("Help")!=string::npos) mKeyCode = Key_Help;
  else if(isKey.find("Home")!=string::npos) mKeyCode = Key_Home;
  else if(isKey.find("Left")!=string::npos) mKeyCode = Key_Left;
  else if(isKey.find("Down")!=string::npos) mKeyCode = Key_Down;
  else if(isKey.find("Menu")!=string::npos) mKeyCode = Key_Menu;
  else if(isKey.find("End")!=string::npos) mKeyCode = Key_End;
  else if(isKey.find("Tab")!=string::npos) mKeyCode = Key_Tab;
  else if(isKey.find("Bar")!=string::npos) mKeyCode = Key_Bar;
  else if(isKey.find("F10")!=string::npos) mKeyCode = Key_F10;
  else if(isKey.find("F11")!=string::npos) mKeyCode = Key_F11;
  else if(isKey.find("F12")!=string::npos) mKeyCode = Key_F12;
  else if(isKey.find("F13")!=string::npos) mKeyCode = Key_F13;
  else if(isKey.find("F14")!=string::npos) mKeyCode = Key_F14;
  else if(isKey.find("F15")!=string::npos) mKeyCode = Key_F15;
  else if(isKey.find("F16")!=string::npos) mKeyCode = Key_F16;
  else if(isKey.find("F17")!=string::npos) mKeyCode = Key_F17;
  else if(isKey.find("F18")!=string::npos) mKeyCode = Key_F18;
  else if(isKey.find("F19")!=string::npos) mKeyCode = Key_F19;
  else if(isKey.find("F20")!=string::npos) mKeyCode = Key_F20;
  else if(isKey.find("F21")!=string::npos) mKeyCode = Key_F21;
  else if(isKey.find("F22")!=string::npos) mKeyCode = Key_F22;
  else if(isKey.find("F23")!=string::npos) mKeyCode = Key_F23;
  else if(isKey.find("F24")!=string::npos) mKeyCode = Key_F24;
  else if(isKey.find("F25")!=string::npos) mKeyCode = Key_F25;
  else if(isKey.find("F26")!=string::npos) mKeyCode = Key_F26;
  else if(isKey.find("F27")!=string::npos) mKeyCode = Key_F27;
  else if(isKey.find("F28")!=string::npos) mKeyCode = Key_F28;
  else if(isKey.find("F29")!=string::npos) mKeyCode = Key_F29;
  else if(isKey.find("F30")!=string::npos) mKeyCode = Key_F30;
  else if(isKey.find("F31")!=string::npos) mKeyCode = Key_F31;
  else if(isKey.find("F32")!=string::npos) mKeyCode = Key_F32;
  else if(isKey.find("F33")!=string::npos) mKeyCode = Key_F33;
  else if(isKey.find("F34")!=string::npos) mKeyCode = Key_F34;
  else if(isKey.find("F35")!=string::npos) mKeyCode = Key_F35;
  else if(isKey.find("At")!=string::npos) mKeyCode = Key_At;
  else if(isKey.find("F1")!=string::npos) mKeyCode = Key_F1;
  else if(isKey.find("F2")!=string::npos) mKeyCode = Key_F2;
  else if(isKey.find("F3")!=string::npos) mKeyCode = Key_F3;
  else if(isKey.find("F4")!=string::npos) mKeyCode = Key_F4;
  else if(isKey.find("F5")!=string::npos) mKeyCode = Key_F5;
  else if(isKey.find("F6")!=string::npos) mKeyCode = Key_F6;
  else if(isKey.find("F7")!=string::npos) mKeyCode = Key_F7;
  else if(isKey.find("F8")!=string::npos) mKeyCode = Key_F8;
  else if(isKey.find("F9")!=string::npos) mKeyCode = Key_F9;
  else if(isKey.find("Up")!=string::npos) mKeyCode = Key_Up;
  else if(isKey.find("0")!=string::npos) mKeyCode = Key_0;
  else if(isKey.find("1")!=string::npos) mKeyCode = Key_1;
  else if(isKey.find("2")!=string::npos) mKeyCode = Key_2;
  else if(isKey.find("3")!=string::npos) mKeyCode = Key_3;
  else if(isKey.find("4")!=string::npos) mKeyCode = Key_4;
  else if(isKey.find("5")!=string::npos) mKeyCode = Key_5;
  else if(isKey.find("6")!=string::npos) mKeyCode = Key_6;
  else if(isKey.find("7")!=string::npos) mKeyCode = Key_7;
  else if(isKey.find("8")!=string::npos) mKeyCode = Key_8;
  else if(isKey.find("9")!=string::npos) mKeyCode = Key_9;
  // For the single keys , make sure our match is right at the end of
  // the string by using reverse find and checking the position we
  // got.
  else if(isKey.rfind("A")==isKey.length()-1) mKeyCode = Key_A;
  else if(isKey.rfind("B")==isKey.length()-1) mKeyCode = Key_B;
  else if(isKey.rfind("C")==isKey.length()-1) mKeyCode = Key_C;
  else if(isKey.rfind("D")==isKey.length()-1) mKeyCode = Key_D;
  else if(isKey.rfind("E")==isKey.length()-1) mKeyCode = Key_E;
  else if(isKey.rfind("F")==isKey.length()-1) mKeyCode = Key_F;
  else if(isKey.rfind("G")==isKey.length()-1) mKeyCode = Key_G;
  else if(isKey.rfind("H")==isKey.length()-1) mKeyCode = Key_H;
  else if(isKey.rfind("I")==isKey.length()-1) mKeyCode = Key_I;
  else if(isKey.rfind("J")==isKey.length()-1) mKeyCode = Key_J;
  else if(isKey.rfind("K")==isKey.length()-1) mKeyCode = Key_K;
  else if(isKey.rfind("L")==isKey.length()-1) mKeyCode = Key_L;
  else if(isKey.rfind("M")==isKey.length()-1) mKeyCode = Key_M;
  else if(isKey.rfind("N")==isKey.length()-1) mKeyCode = Key_N;
  else if(isKey.rfind("O")==isKey.length()-1) mKeyCode = Key_O;
  else if(isKey.rfind("P")==isKey.length()-1) mKeyCode = Key_P;
  else if(isKey.rfind("Q")==isKey.length()-1) mKeyCode = Key_Q;
  else if(isKey.rfind("R")==isKey.length()-1) mKeyCode = Key_R;
  else if(isKey.rfind("S")==isKey.length()-1) mKeyCode = Key_S;
  else if(isKey.rfind("T")==isKey.length()-1) mKeyCode = Key_T;
  else if(isKey.rfind("U")==isKey.length()-1) mKeyCode = Key_U;
  else if(isKey.rfind("V")==isKey.length()-1) mKeyCode = Key_V;
  else if(isKey.rfind("W")==isKey.length()-1) mKeyCode = Key_W;
  else if(isKey.rfind("X")==isKey.length()-1) mKeyCode = Key_X;
  else if(isKey.rfind("Y")==isKey.length()-1) mKeyCode = Key_Y;
  else if(isKey.rfind("Z")==isKey.length()-1) mKeyCode = Key_Z;
  else if(isKey.rfind("a")==isKey.length()-1) mKeyCode = Key_A;
  else if(isKey.rfind("b")==isKey.length()-1) mKeyCode = Key_B;
  else if(isKey.rfind("c")==isKey.length()-1) mKeyCode = Key_C;
  else if(isKey.rfind("d")==isKey.length()-1) mKeyCode = Key_D;
  else if(isKey.rfind("e")==isKey.length()-1) mKeyCode = Key_E;
  else if(isKey.rfind("f")==isKey.length()-1) mKeyCode = Key_F;
  else if(isKey.rfind("g")==isKey.length()-1) mKeyCode = Key_G;
  else if(isKey.rfind("h")==isKey.length()-1) mKeyCode = Key_H;
  else if(isKey.rfind("i")==isKey.length()-1) mKeyCode = Key_I;
  else if(isKey.rfind("j")==isKey.length()-1) mKeyCode = Key_J;
  else if(isKey.rfind("k")==isKey.length()-1) mKeyCode = Key_K;
  else if(isKey.rfind("l")==isKey.length()-1) mKeyCode = Key_L;
  else if(isKey.rfind("m")==isKey.length()-1) mKeyCode = Key_M;
  else if(isKey.rfind("n")==isKey.length()-1) mKeyCode = Key_N;
  else if(isKey.rfind("o")==isKey.length()-1) mKeyCode = Key_O;
  else if(isKey.rfind("p")==isKey.length()-1) mKeyCode = Key_P;
  else if(isKey.rfind("q")==isKey.length()-1) mKeyCode = Key_Q;
  else if(isKey.rfind("r")==isKey.length()-1) mKeyCode = Key_R;
  else if(isKey.rfind("s")==isKey.length()-1) mKeyCode = Key_S;
  else if(isKey.rfind("t")==isKey.length()-1) mKeyCode = Key_T;
  else if(isKey.rfind("u")==isKey.length()-1) mKeyCode = Key_U;
  else if(isKey.rfind("v")==isKey.length()-1) mKeyCode = Key_V;
  else if(isKey.rfind("w")==isKey.length()-1) mKeyCode = Key_W;
  else if(isKey.rfind("x")==isKey.length()-1) mKeyCode = Key_X;
  else if(isKey.rfind("y")==isKey.length()-1) mKeyCode = Key_Y;
  else if(isKey.rfind("z")==isKey.length()-1) mKeyCode = Key_Z;
  else mKeyCode = Key_unknown;
}

string
ScubaKeyCombo::ToString () {
  
  string sModifers = "";
  if( mbControl ) sModifers += "Ctrl ";
  if( mbShift ) sModifers += "Shift ";
  if( mbAlt ) sModifers += "Alt ";
  if( mbMeta ) sModifers += "Meta ";

  string sKey;
  switch( mKeyCode ) {
  case Key_Escape: sKey = "Escape"; break;
  case Key_Tab: sKey = "Tab"; break;
  case Key_Backtab: sKey = "Backtab"; break;
  case Key_Backspace: sKey = "Backspace"; break;
  case Key_Return: sKey = "Return"; break;
  case Key_Enter: sKey = "Enter"; break;
  case Key_Insert: sKey = "Insert"; break;
  case Key_Delete: sKey = "Delete"; break;
  case Key_Pause: sKey = "Pause"; break;
  case Key_Print: sKey = "Print"; break;
  case Key_SysReq: sKey = "SysReq"; break;
  case Key_Clear: sKey = "Clear"; break;
  case Key_Home: sKey = "Home"; break;
  case Key_End: sKey = "End"; break;
  case Key_Left: sKey = "Left"; break;
  case Key_Up: sKey = "Up"; break;
  case Key_Right: sKey = "Right"; break;
  case Key_Down: sKey = "Down"; break;
  case Key_Prior: sKey = "PageUp"; break;
  case Key_Next: sKey = "PageDown"; break;
  case Key_Shift: sKey = "Shift"; break;
  case Key_Control: sKey = "Control"; break;
  case Key_Meta: sKey = "Meta"; break;
  case Key_Alt: sKey = "Alt"; break;
  case Key_CapsLock: sKey = "CapsLock"; break;
  case Key_NumLock: sKey = "NumLock"; break;
  case Key_ScrollLock: sKey = "ScrollLock"; break;
  case Key_F1: sKey = "F1"; break;
  case Key_F2: sKey = "F2"; break;
  case Key_F3: sKey = "F3"; break;
  case Key_F4: sKey = "F4"; break;
  case Key_F5: sKey = "F5"; break;
  case Key_F6: sKey = "F6"; break;
  case Key_F7: sKey = "F7"; break;
  case Key_F8: sKey = "F8"; break;
  case Key_F9: sKey = "F9"; break;
  case Key_F10: sKey = "F10"; break;
  case Key_F11: sKey = "F11"; break;
  case Key_F12: sKey = "F12"; break;
  case Key_F13: sKey = "F13"; break;
  case Key_F14: sKey = "F14"; break;
  case Key_F15: sKey = "F15"; break;
  case Key_F16: sKey = "F16"; break;
  case Key_F17: sKey = "F17"; break;
  case Key_F18: sKey = "F18"; break;
  case Key_F19: sKey = "F19"; break;
  case Key_F20: sKey = "F20"; break;
  case Key_F21: sKey = "F21"; break;
  case Key_F22: sKey = "F22"; break;
  case Key_F23: sKey = "F23"; break;
  case Key_F24: sKey = "F24"; break;
  case Key_F25: sKey = "F25"; break;
  case Key_F26: sKey = "F26"; break;
  case Key_F27: sKey = "F27"; break;
  case Key_F28: sKey = "F28"; break;
  case Key_F29: sKey = "F29"; break;
  case Key_F30: sKey = "F30"; break;
  case Key_F31: sKey = "F31"; break;
  case Key_F32: sKey = "F32"; break;
  case Key_F33: sKey = "F33"; break;
  case Key_F34: sKey = "F34"; break;
  case Key_F35: sKey = "F35"; break;
  case Key_Super_L: sKey = "Super_L"; break;
  case Key_Super_R: sKey = "Super_R"; break;
  case Key_Menu: sKey = "Menu"; break;
  case Key_Hyper_L: sKey = "Hyper_L"; break;
  case Key_Hyper_R: sKey = "Hyper_R"; break;
  case Key_Help: sKey = "Help"; break;
  case Key_Direction_L: sKey = "Direction_L"; break;
  case Key_Direction_R: sKey = "Direction_R"; break;
  case Key_Space: sKey = "Space"; break;
  case Key_Exclam: sKey = "Exclam"; break;
  case Key_QuoteDbl: sKey = "QuoteDbl"; break;
  case Key_NumberSign: sKey = "NumberSign"; break;
  case Key_Dollar: sKey = "Dollar"; break;
  case Key_Percent: sKey = "Percent"; break;
  case Key_Ampersand: sKey = "Ampersand"; break;
  case Key_Apostrophe: sKey = "Apostrophe"; break;
  case Key_ParenLeft: sKey = "ParenLeft"; break;
  case Key_ParenRight: sKey = "ParenRight"; break;
  case Key_Asterisk: sKey = "Asterisk"; break;
  case Key_Plus: sKey = "Plus"; break;
  case Key_Comma: sKey = "Comma"; break;
  case Key_Minus: sKey = "Minus"; break;
  case Key_Period: sKey = "Period"; break;
  case Key_Slash: sKey = "Slash"; break;
  case Key_0: sKey = "0"; break;
  case Key_1: sKey = "1"; break;
  case Key_2: sKey = "2"; break;
  case Key_3: sKey = "3"; break;
  case Key_4: sKey = "4"; break;
  case Key_5: sKey = "5"; break;
  case Key_6: sKey = "6"; break;
  case Key_7: sKey = "7"; break;
  case Key_8: sKey = "8"; break;
  case Key_9: sKey = "9"; break;
  case Key_Colon: sKey = "Colon"; break;
  case Key_Semicolon: sKey = "Semicolon"; break;
  case Key_Less: sKey = "Less"; break;
  case Key_Equal: sKey = "Equal"; break;
  case Key_Greater: sKey = "Greater"; break;
  case Key_Question: sKey = "Question"; break;
  case Key_At: sKey = "At"; break;
  case Key_A: sKey = "A"; break;
  case Key_B: sKey = "B"; break;
  case Key_C: sKey = "C"; break;
  case Key_D: sKey = "D"; break;
  case Key_E: sKey = "E"; break;
  case Key_F: sKey = "F"; break;
  case Key_G: sKey = "G"; break;
  case Key_H: sKey = "H"; break;
  case Key_I: sKey = "I"; break;
  case Key_J: sKey = "J"; break;
  case Key_K: sKey = "K"; break;
  case Key_L: sKey = "L"; break;
  case Key_M: sKey = "M"; break;
  case Key_N: sKey = "N"; break;
  case Key_O: sKey = "O"; break;
  case Key_P: sKey = "P"; break;
  case Key_Q: sKey = "Q"; break;
  case Key_R: sKey = "R"; break;
  case Key_S: sKey = "S"; break;
  case Key_T: sKey = "T"; break;
  case Key_U: sKey = "U"; break;
  case Key_V: sKey = "V"; break;
  case Key_W: sKey = "W"; break;
  case Key_X: sKey = "X"; break;
  case Key_Y: sKey = "Y"; break;
  case Key_Z: sKey = "Z"; break;
  case Key_BracketLeft: sKey = "BracketLeft"; break;
  case Key_Backslash: sKey = "Backslash"; break;
  case Key_BracketRight: sKey = "BracketRight"; break;
  case Key_AsciiCircum: sKey = "AsciiCircum"; break;
  case Key_Underscore: sKey = "Underscore"; break;
  case Key_QuoteLeft: sKey = "QuoteLeft"; break;
  case Key_BraceLeft: sKey = "BraceLeft"; break;
  case Key_Bar: sKey = "Bar"; break;
  case Key_BraceRight: sKey = "BraceRight"; break;
  case Key_AsciiTilde: sKey = "AsciiTilde"; break;
  case Key_unknown: sKey = "unknown"; break;
  }

  return sModifers + sKey;
}

bool
ScubaKeyCombo::IsSameAs ( ScubaKeyCombo& iCombo ) {

  return ( GetKeyCode() == iCombo.GetKeyCode() &&
	   IsShiftKeyDown() == iCombo.IsShiftKeyDown() &&
	   IsAltKeyDown() == iCombo.IsAltKeyDown() &&
	   IsMetaKeyDown() == iCombo.IsMetaKeyDown() &&
	   IsControlKeyDown() == iCombo.IsControlKeyDown() );
}

ostream& 
operator <<  ( ostream& os, ScubaKeyCombo iKey ) { 

  os.setf(ios::hex,ios::basefield);
  os << "Key " << iKey.ToString() << " (code 0x" << iKey.GetKeyCode()
     << ")";
  os.unsetf(ios::hex);
  return os;
}
