#ifndef ScubaTransform_h
#define ScubaTransform_h

#include "string_fixed.h"
extern "C" {
#include "matrix.h"
}
#include "Transform44.h"
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "IDTracker.h"
#include "Broadcaster.h"
#include "Point3.h"

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

class ScubaTransform : public Transform44, public IDTracker<ScubaTransform>, public TclCommandListener, public Broadcaster {

  friend class ScubaTransformTester;

 public:

  ScubaTransform();
  virtual ~ScubaTransform();

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  void SetLabel( std::string isLabel ) { msLabel = isLabel; }
  std::string GetLabel() { return msLabel; }

  Transform44& operator=(ScubaTransform& iTransform) {
    m.SetMatrix( iTransform.GetMainMatrix() );
    ValuesChanged();
    return *this;
  }

 protected:

  virtual void ValuesChanged ();

  static ScubaTransformStaticTclListener mStaticListener;

  std::string msLabel;
};

std::ostream& operator << ( std::ostream&, ScubaTransform& iTransform  );

#endif

