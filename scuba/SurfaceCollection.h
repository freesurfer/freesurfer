#ifndef SurfaceCollection_H
#define SurfaceCollection_H

#include <string>
#include "DataCollection.h"
extern "C" {
#include "mrisurf.h"
}

class SurfaceCollection : public DataCollection {

 public:
  SurfaceCollection( std::string& fnMRI );
  virtual ~SurfaceCollection();

  MRIS* GetMRIS() { return mMRIS; }

 protected:
  MRIS* mMRIS;

};


#endif
