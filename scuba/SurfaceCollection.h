#ifndef SurfaceCollection_H
#define SurfaceCollection_H

#include <string>
#include "DataCollection.h"
extern "C" {
#include "mrisurf.h"
}

class SurfaceCollection : public DataCollection {

 public:
  SurfaceCollection();
  virtual ~SurfaceCollection();

  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription() { return "Surface"; }

  void SetSurfaceFileName ( std::string& ifnMRIS );

  MRIS* GetMRIS();

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  virtual ScubaROI* DoNewROI ();


 protected:
  std::string mfnMRIS;
  MRIS* mMRIS;

};


#endif
