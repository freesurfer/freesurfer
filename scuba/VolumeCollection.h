#ifndef VolumeCollection_H
#define VolumeCollection_H

#include <string>
#include "DataCollection.h"
extern "C" {
#include "mri.h"
}

class VolumeCollection : public DataCollection {

  friend class VolumeCollectionTester;

 public:
  VolumeCollection ();
  virtual ~VolumeCollection ();

  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription() { return "Volume"; }

  void SetFileName ( std::string& ifnMRI );

  MRI* GetMRI ();
  
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

protected:
  std::string mfnMRI;
  MRI* mMRI;

};


#endif
