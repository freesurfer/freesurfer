#ifndef VolumeCollection_H
#define VolumeCollection_H

#include <string>
#include "DataCollection.h"
extern "C" {
#include "mri.h"
}

class VolumeCollection : public DataCollection {

 public:
  VolumeCollection( std::string& fnMRI );
  virtual ~VolumeCollection();

  MRI* GetMRI() { return mMRI; }

 protected:
  MRI* mMRI;

};


#endif
