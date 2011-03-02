/**
 * @file  test_ScubaKeyCombo.cpp
 * @brief test ScubaKeyCombo class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.7 $
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


#include <stdlib.h>
#include "string_fixed.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "ScubaKeyCombo.h"
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

const char* Progname = "test_ScubaKeyCombo";


class ScubaKeyComboTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void
ScubaKeyComboTester::Test ( Tcl_Interp* iInterp ) {

  try {

    string asModifiers[] = {"Ctrl ", "Shift ", "Alt ", "Meta ",
                            "Ctrl Shift ", "Ctrl Alt ", "Ctrl Meta ",
                            "Shift Alt ", "Shift Meta ",
                            "Alt Meta ",
                            "Ctrl Shift Alt ", "Shift Alt Meta " };
    string asUnits[] = {"Escape", "Tab", "Backtab", "Backspace", "Return", "Enter", "Insert", "Delete", "Pause", "Print", "SysReq", "Clear", "Home", "End", "Left", "Up", "Right", "Down", "PageUp", "PageDown", "CapsLock", "NumLock", "ScrollLock", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12", "F13", "F14", "F15", "F16", "F17", "F18", "F19", "F20", "F21", "F22", "F23", "F24", "F25", "F26", "F27", "F28", "F29", "F30", "F31", "F32", "F33", "F34", "F35", "Super_L", "Super_R", "Menu", "Hyper_L", "Hyper_R", "Help", "Direction_L", "Direction_R", "Space", "Exclam", "QuoteDbl", "NumberSign", "Dollar", "Percent", "Ampersand", "Apostrophe", "ParenLeft", "ParenRight", "Asterisk", "Plus", "Comma", "Minus", "Period", "Slash", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "Colon", "Semicolon", "Less", "Equal", "Greater", "Question", "At", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "BracketLeft", "Backslash", "BracketRight", "AsciiCircum", "Underscore", "QuoteLeft", "BraceLeft", "Bar", "BraceRight", "AsciiTilde" };

    vector<string> vsKey;
    for ( int nModifier = 0; nModifier < 12; nModifier++ ) {
      for ( int nUnit = 0; nUnit < 134; nUnit++ ) {
        string sKey = asModifiers[nModifier] + asUnits[nUnit];
        vsKey.push_back( sKey );
      }
    }

    while ( vsKey.size() > 0 ) {
      string sKey = vsKey.back();
      ScubaKeyCombo* key = ScubaKeyCombo::MakeKeyCombo();
      key->SetFromString( sKey );
      {
        stringstream ssError;
        ssError << "Failed test " << *key << " == " << sKey;
        Assert( (key->ToString() == sKey), ssError.str() );
      }

      vsKey.pop_back();
      delete key;
    }

    ScubaKeyCombo* a = ScubaKeyCombo::MakeKeyCombo();
    a->SetFromString( "Ctrl a" );
    ScubaKeyCombo* a2 = ScubaKeyCombo::MakeKeyCombo();
    a2->SetFromString( "Ctrl a" );
    ScubaKeyCombo* b = ScubaKeyCombo::MakeKeyCombo();
    b->SetFromString( "Ctrl b" );
    {
      stringstream ssError;
      ssError << "Failed IsSameAs " << *a << ", " << *b;
      Assert( (!a->IsSameAs( b )), ssError.str() );
    }
    {
      stringstream ssError;
      ssError << "Failed IsSameAs " << *a << ", " << *a2;
      Assert( (a->IsSameAs( a2 )), ssError.str() );
    }

    delete a;
    delete a2;
    delete b;
  } catch ( exception& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
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

    ScubaKeyCombo::SetFactory( new ScubaKeyComboFactory() );

    ScubaKeyComboTester tester0;
    tester0.Test( interp );


  } catch ( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}
