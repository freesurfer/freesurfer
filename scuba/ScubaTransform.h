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

class VolumeCollection;

class ScubaTransformStaticTclListener : public DebugReporter, public TclCommandListener {

 public:
  ScubaTransformStaticTclListener ();
  ~ScubaTransformStaticTclListener ();
  
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );
};

class ScubaTransform : public Transform44, public IDTracker<ScubaTransform>, public TclCommandListener, public Broadcaster {

  friend class ScubaTransformTester;

 public:

  ScubaTransform();
  virtual ~ScubaTransform();

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  void SetLabel( std::string isLabel ) { msLabel = isLabel; }
  std::string GetLabel() { return msLabel; }

  std::string GetValuesAsString();

  Transform44& operator=(ScubaTransform& iTransform) {
    m.SetMatrix( iTransform.GetMainMatrix() );
    ValuesChanged();
    return *this;
  }

  // For making this volume a registration volume.
  void TreatAsRegistration ( VolumeCollection& iSourceVolume, 
			     VolumeCollection& iDestVolume );
  void TreatAsNative ();

  bool IsRegistration () { return mIsRegistration; }
  int GetRegistrationSource ();
  int GetRegistrationDestination ();

 protected:

  virtual void ValuesChanged ();

  static ScubaTransformStaticTclListener mStaticListener;

  std::string msLabel;

  // For making this a registration transform. Save the source and
  // dest volume IDs so we can unmake it.
  bool mIsRegistration;
  VolumeCollection* mSourceVolume;
  VolumeCollection* mDestVolume;
};

std::ostream& operator << ( std::ostream&, ScubaTransform& iTransform  );

#endif

