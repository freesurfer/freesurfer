#ifndef ScubaTransform_h
#define ScubaTransform_h

#include <string>
extern "C" {
#include "matrix.h"
}
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "IDTracker.h"
#include "Broadcaster.h"

class ScubaTransformStaticTclListener : public DebugReporter, public TclCommandListener {

 public:
  ScubaTransformStaticTclListener ();
  ~ScubaTransformStaticTclListener ();
  
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );
};

// Note that all 16 element matrices used in this class are in openGL
// style format:
// [ 0   4   8  12 ]
// [ 1   5   9  13 ]
// [ 2   6  10  14 ]
// [ 3   7  11  15 ]

class ScubaTransform : public DebugReporter, public IDTracker<ScubaTransform>, public TclCommandListener, public Broadcaster {

  friend class ScubaTransformTester;

 public:

  ScubaTransform();
  ~ScubaTransform();

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  void SetTransform ( float i0j0, float i1j0, float i2j0, float i3j0,
		      float i0j1, float i1j1, float i2j1, float i3j1,
		      float i0j2, float i1j2, float i2j2, float i3j2,
		      float i0j3, float i1j3, float i2j3, float i3j3 );

  void MakeIdentity ();

  void SetLabel( std::string isLabel ) { msLabel = isLabel; }
  std::string GetLabel() { return msLabel; }

  void MultiplyVector3 ( float iVector[3], float oVector[3] );
  void InvMultiplyVector3 ( float iVector[3], float oVector[3] );

 protected:

  // This is a way of changing an individual value, but note that it
  // DOES NOT call the ValuesChanged() function afterwards, which is
  // necessary to update listeners and calc in the inverse. The caller
  // should do this appropriately.
  inline void SetCR ( int iCol, int iRow, float iValue ) {
    *MATRIX_RELT(m,(iRow+1),(iCol+1)) = iValue;
  }

  inline float GetCR ( int iCol, int iRow ) {
    return *MATRIX_RELT(m,(iRow+1),(iCol+1));
  }

  void ValuesChanged ();

  void CalculateInverse ();

  static ScubaTransformStaticTclListener mStaticListener;

  std::string msLabel;

  MATRIX* m;
  MATRIX* mInv;
  MATRIX* mTmpVec4Src;
  MATRIX* mTmpVec4Dest;

};




#endif

