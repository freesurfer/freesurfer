#ifndef VolumeCollection_H
#define VolumeCollection_H

#include <string>
#include "DataCollection.h"
#include "ScubaROIVolume.h"
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
  void UpdateMRIValueRange ();
  float GetMRIMinValue () { return mMRIMinValue; }
  float GetMRIMaxValue () { return mMRIMaxValue; }

  void RASToMRIIndex ( float iRAS[3], int oIndex[3] );
  void RASToMRIIndex ( float iRAS[3], float oIndex[3] );

  bool IsRASInMRIBounds ( float iRAS[3] );

  float GetMRINearestValueAtRAS ( float iRAS[3] );
  float GetMRITrilinearValueAtRAS ( float iRAS[3] );
  float GetMRISincValueAtRAS ( float iRAS[3] );
  
  void SetMRIValueAtRAS ( float iRAS[3], float iValue );
  
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  virtual ScubaROI* DoNewROI ();

  void SelectRAS ( float iRAS[3] );
  void UnselectRAS ( float iRAS[3] );
  bool IsRASSelected ( float iRAS[3] );

protected:
  std::string mfnMRI;
  MRI* mMRI;

  MATRIX* mWorldToIndexMatrix;
  VECTOR* mWorldCoord;
  VECTOR* mIndexCoord;

  float mMRIMinValue, mMRIMaxValue;
};


#endif
