#ifndef ScubaROI_h
#define ScubaROI_h

#include <string>
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "IDTracker.h"

class ScubaROIStaticTclListener : public DebugReporter, public TclCommandListener {

 public:
  ScubaROIStaticTclListener ();
  ~ScubaROIStaticTclListener ();

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );
};

class ScubaROI : public DebugReporter, public IDTracker<ScubaROI>, public TclCommandListener {

  friend class ScubaROITester;

 public:

  ScubaROI ();
  virtual ~ScubaROI ();

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  enum Type { Structure, Free };
  void SetType( Type iType ) { mType = iType; }
  Type GetType() { return mType; }
  
  void SetStructure( int iStructure ) { mStructure = iStructure; }
  int GetStructure() { return mStructure; }

  void SetColor( int iColor[3] );
  void GetColor( int oColor[3] );
  
  void SetLabel( std::string isLabel ) { msLabel = isLabel; }
  std::string GetLabel() { return msLabel; }
  
 protected:

  std::string msLabel;

  Type mType;

  int mStructure;
  int mColor[3];

  static ScubaROIStaticTclListener mStaticListener;
};



#endif
