#include <stdlib.h>
#include "string_fixed.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "TclScubaKeyCombo.h"
extern "C" {
#include "mri.h"
}
#include "Scuba-impl.h"

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }


#define AssertTclOK(x) \
    if( TCL_OK != (x) ) { \
      stringstream ssError; \
      ssError << "Tcl_Eval returned not TCL_OK: " << endl  \
	     << "Command: " << sCommand << endl \
	     << "Result: " << iInterp->result; \
      throw runtime_error( ssError.str() ); \
    } \


using namespace std;

char* Progname = "test_TclScubaKeyCombo";



class TclScubaKeyComboTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void 
TclScubaKeyComboTester::Test ( Tcl_Interp* iInterp ) {

  try {

    struct keyTest { string sKey; int keyCode; };
    keyTest aKeyTests[] = { 
      {"Prior", ScubaKeyCombo::Key_PageUp},
      {"minus", ScubaKeyCombo::Key_Minus},
      {"equal", ScubaKeyCombo::Key_Equal},
      {"Next", ScubaKeyCombo::Key_PageDown},
      {"space", ScubaKeyCombo::Key_Space},
      {"exclam", ScubaKeyCombo::Key_Exclam},
      {"quotedbl", ScubaKeyCombo::Key_QuoteDbl},
      {"numbersign", ScubaKeyCombo::Key_NumberSign},
      {"dollar", ScubaKeyCombo::Key_Dollar},
      {"percent", ScubaKeyCombo::Key_Percent},
      {"ampersand", ScubaKeyCombo::Key_Ampersand},
      {"apostrophe", ScubaKeyCombo::Key_Apostrophe},
      {"parenleft", ScubaKeyCombo::Key_ParenLeft},
      {"parenright", ScubaKeyCombo::Key_ParenRight},
      {"asterisk", ScubaKeyCombo::Key_Asterisk},
      {"plus", ScubaKeyCombo::Key_Plus},
      {"comma", ScubaKeyCombo::Key_Comma},
      {"minus", ScubaKeyCombo::Key_Minus},
      {"period", ScubaKeyCombo::Key_Period},
      {"slash", ScubaKeyCombo::Key_Slash},
      {"colon", ScubaKeyCombo::Key_Colon},
      {"semicolon", ScubaKeyCombo::Key_Semicolon},
      {"less", ScubaKeyCombo::Key_Less},
      {"equal", ScubaKeyCombo::Key_Equal},
      {"greater", ScubaKeyCombo::Key_Greater},
      {"at", ScubaKeyCombo::Key_At},
      {"bracketleft", ScubaKeyCombo::Key_BracketLeft},
      {"backslash", ScubaKeyCombo::Key_Backslash},
      {"bracketright", ScubaKeyCombo::Key_BracketRight},
      {"asciicirum", ScubaKeyCombo::Key_AsciiCircum},
      {"underscore", ScubaKeyCombo::Key_Underscore},
      {"grave", ScubaKeyCombo::Key_QuoteLeft},
      {"braceleft", ScubaKeyCombo::Key_BraceLeft},
      {"bar", ScubaKeyCombo::Key_Bar},
      {"braceright", ScubaKeyCombo::Key_BraceRight},
      {"asciitilde", ScubaKeyCombo::Key_AsciiTilde},
      {"a", ScubaKeyCombo::Key_A},
      {"b", ScubaKeyCombo::Key_B},
      {"c", ScubaKeyCombo::Key_C},
      {"d", ScubaKeyCombo::Key_D},
      {"e", ScubaKeyCombo::Key_E},
      {"f", ScubaKeyCombo::Key_F},
      {"g", ScubaKeyCombo::Key_G},
      {"h", ScubaKeyCombo::Key_H},
      {"i", ScubaKeyCombo::Key_I},
      {"j", ScubaKeyCombo::Key_J},
      {"k", ScubaKeyCombo::Key_K},
      {"l", ScubaKeyCombo::Key_L},
      {"m", ScubaKeyCombo::Key_M},
      {"n", ScubaKeyCombo::Key_N},
      {"o", ScubaKeyCombo::Key_O},
      {"p", ScubaKeyCombo::Key_P},
      {"q", ScubaKeyCombo::Key_Q},
      {"r", ScubaKeyCombo::Key_R},
      {"s", ScubaKeyCombo::Key_S},
      {"t", ScubaKeyCombo::Key_T},
      {"u", ScubaKeyCombo::Key_U},
      {"v", ScubaKeyCombo::Key_V},
      {"w", ScubaKeyCombo::Key_W},
      {"x", ScubaKeyCombo::Key_X},
      {"y", ScubaKeyCombo::Key_Y},
      {"z", ScubaKeyCombo::Key_Z}};

    char sCommand[1024];
    int rTcl;
    for( int nKey = 0; nKey < 62; nKey++ ) {
      TclScubaKeyCombo key;
      key.SetFromString(aKeyTests[nKey].sKey);
      {
	stringstream ssError;
	ssError.setf(ios::hex,ios::basefield);
	ssError << "Failed test string " << aKeyTests[nKey].sKey 
		<< " code 0x" << aKeyTests[nKey].keyCode 
		<< " made key " << key;
	Assert( (key.GetKeyCode() == aKeyTests[nKey].keyCode), ssError.str() );
      }
      
      // Test the tcl commands.
      sprintf( sCommand, "ConvertTkInputStringToScubaKeyComboString %s",
	       aKeyTests[nKey].sKey.c_str() );
      rTcl = Tcl_Eval( iInterp, sCommand );
      AssertTclOK( rTcl );
      const char* sTclResult = Tcl_GetStringResult( iInterp );

      ScubaKeyCombo testKey;
      testKey.SetFromString( sTclResult );
      {
	stringstream ssError;
	ssError << "Failed tcl return " << aKeyTests[nKey].sKey 
		<< " returned string \"" << sTclResult
		<< "\" which made " << testKey;
	Assert( key.IsSameAs( testKey ), ssError.str() );
      }
    }

  }
  catch( exception& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed." << endl;
    exit( 1 );
  }
}



int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {

    Tcl_Interp* interp = Tcl_CreateInterp();
    Assert( interp, "Tcl_CreateInterp returned null" );
  
    int rTcl = Tcl_Init( interp );
    Assert( TCL_OK == rTcl, "Tcl_Init returned not TCL_OK" );
    
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.SetOutputStreamToCerr();
    commandMgr.Start( interp );

    TclScubaKeyComboStaticTclListener::GetListener();

    TclScubaKeyComboTester tester0;
    tester0.Test( interp );

 
  }
  catch( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}
