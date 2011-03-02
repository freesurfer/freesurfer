/**
 * @file  test_TclScubaKeyCombo.cpp
 * @brief test TclScubaKeyCombo class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.8 $
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

const char* Progname = "test_TclScubaKeyCombo";



class TclScubaKeyComboTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void
TclScubaKeyComboTester::Test ( Tcl_Interp* iInterp ) {

  try {

    struct modTest {
      string sMod;
      bool bCtrl, bShift, bAlt, bMeta;
    };
    modTest aModifiers[] = {{"", false, false, false, false},
			    {"Ctrl ", true, false, false, false},
                            {"Shift ", false, true, false, false},
                            {"Alt ", false, false, true, false},
                            {"Meta ", false, false, false, true},
                            {"Ctrl Shift ", true, true, false, false},
                            {"Ctrl Alt ", true, false, true, false},
                            {"Ctrl Meta ", true, false, false, true},
                            {"Shift Alt ", false, true, true, false},
                            {"Shift Meta ", false, true, false, true},
                            {"Alt Meta ", false, false, true, true},
                            {"Ctrl Shift Alt ", true, true, true, false},
                            {"Shift Alt Meta ", false, true, true, true} };

    struct keyTest {
      string sKey;
      int keyCode;
    };
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
    for ( int nMod = 0; nMod < 12; nMod++ ) {
      for ( int nKey = 0; nKey < 62; nKey++ ) {
        string sKey = aModifiers[nMod].sMod + aKeyTests[nKey].sKey;

        ScubaKeyCombo* key = ScubaKeyCombo::MakeKeyCombo();
        key->SetFromString(sKey);
        {
          stringstream ssError;
          ssError.setf(ios::hex,ios::basefield);
          ssError << "Failed test string " << sKey
          << " code 0x" << aKeyTests[nKey].keyCode
          << " made key " << *key;
          Assert( (key->GetKeyCode() ==aKeyTests[nKey].keyCode),
                  ssError.str() );
        }
        {
          stringstream ssError;
          ssError << "Failed test shift mod " << sKey
          << " made key " << *key;
          Assert( (key->IsShiftKeyDown() == aModifiers[nMod].bShift),
                  ssError.str() );
        }
        {
          stringstream ssError;
          ssError << "Failed test ctrl mod " << sKey
          << " made key " << *key;
          Assert( (key->IsControlKeyDown() == aModifiers[nMod].bCtrl),
                  ssError.str() );
        }
        {
          stringstream ssError;
          ssError << "Failed test alt mod " << sKey
          << " made key " << *key;
          Assert( (key->IsAltKeyDown() == aModifiers[nMod].bAlt),
                  ssError.str() );
        }
        {
          stringstream ssError;
          ssError << "Failed test meta mod " << sKey
          << " made key " << *key;
          Assert( (key->IsMetaKeyDown() == aModifiers[nMod].bMeta),
                  ssError.str() );
        }

        // Test the tcl commands.
        sprintf( sCommand,
                 "ConvertTkInputStringToScubaKeyComboString %s %d %d %d %d",
                 aKeyTests[nKey].sKey.c_str(), aModifiers[nMod].bShift,
                 aModifiers[nMod].bMeta, aModifiers[nMod].bAlt,
                 aModifiers[nMod].bCtrl );
        rTcl = Tcl_Eval( iInterp, sCommand );
        AssertTclOK( rTcl );
        const char* sTclResult = Tcl_GetStringResult( iInterp );

        ScubaKeyCombo* testKey = ScubaKeyCombo::MakeKeyCombo();
        testKey->SetFromString( sTclResult );
        {
          stringstream ssError;
          ssError << "Failed tcl return " << aKeyTests[nKey].sKey
          << " shift " << aModifiers[nMod].bShift
          << " meta " << aModifiers[nMod].bMeta
          << " alt " << aModifiers[nMod].bAlt
          << " ctrl " << aModifiers[nMod].bCtrl
          << " returned string \"" << sTclResult
          << "\" which made " << *testKey
          << " which doesn't match " << *key;
          Assert( key->IsSameAs( testKey ), ssError.str() );
        }

        delete key;
        delete testKey;
      }
    }

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

    TclScubaKeyComboStaticTclListener::GetListener();
    ScubaKeyCombo::SetFactory( new TclScubaKeyComboFactory() );

    TclScubaKeyComboTester tester0;
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
