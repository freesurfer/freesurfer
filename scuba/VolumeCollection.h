#ifndef VolumeCollection_H
#define VolumeCollection_H

#include <string>
#include "DataCollection.h"
#include "ScubaROIVolume.h"
extern "C" {
#include "mri.h"
}
#include "Volume3.h"
#include "Point3.h"

class VolumeCollection : public DataCollection {

  friend class VolumeCollectionTester;
  friend class VolumeCollectionFlooder;

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

  float GetVoxelXSize ();
  float GetVoxelYSize ();
  float GetVoxelZSize ();

  void UpdateRASBounds ();

  void RASToMRIIndex ( float iRAS[3], int oIndex[3] );
  void RASToMRIIndex ( float iRAS[3], float oIndex[3] );
  void MRIIndexToRAS ( int iIndex[3], float oRAS[3] );
  void MRIIndexToRAS ( float iIndex[3], float oRAS[3] );

  bool IsRASInMRIBounds ( float iRAS[3] );
  bool IsMRIIndexInMRIBounds ( int iIndex[3] );

  float GetMRINearestValueAtRAS ( float iRAS[3] );
  float GetMRITrilinearValueAtRAS ( float iRAS[3] );
  float GetMRISincValueAtRAS ( float iRAS[3] );
  
  void SetMRIValueAtRAS ( float iRAS[3], float iValue );
  
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  virtual ScubaROI* DoNewROI ();

  void InitSelectionVolume ();
  void SelectRAS ( float iRAS[3] );
  void UnselectRAS ( float iRAS[3] );
  bool IsRASSelected ( float iRAS[3], int oColor[3] );
  bool IsOtherRASSelected ( float iRAS[3], int iThisROIID );

  void InitEdgeVolume ();
  void MarkRASEdge ( float iRAS[3] );
  void UnmarkRASEdge ( float iRAS[3] );
  bool IsRASEdge ( float iRAS[3] );

  void GetRASPointsInCube ( float iCenterRAS[3], int iRadius,
			    bool ibBrushX, bool ibBrushY, bool ibBrushZ,
			    std::list<Point3<float> >& oPoints );
  void GetRASPointsInSphere ( float iCenterRAS[3], int iRadius,
			      bool ibBrushX, bool ibBrushY, bool ibBrushZ,
			      std::list<Point3<float> >& oPoints );

  void WriteROIToLabel ( int iROIID, std::string ifnLabel );
  int  NewROIFromLabel ( std::string ifnLabel );
  void WriteROIsToSegmentation ( std::string ifnVolume );

protected:
  std::string mfnMRI;
  MRI* mMRI;

  MATRIX* mWorldToIndexMatrix;
  MATRIX* mIndexToWorldMatrix;
  VECTOR* mWorldCoord;
  VECTOR* mIndexCoord;

  float mMRIMinValue, mMRIMaxValue;

  Volume3<bool>* mEdgeVoxels;

  Volume3<bool>* mSelectedVoxels;

  float mMinRASBounds[3];
  float mMaxRASBounds[3];
};

class VolumeCollectionFlooder {
 public:
  VolumeCollectionFlooder ();
  virtual ~VolumeCollectionFlooder ();

  // Override these to do startup and shutdown stuff, also allow the
  // user or program to cancel the flood.
  virtual void DoBegin ();  
  virtual void DoEnd ();
  virtual bool DoStopRequested ();

  virtual bool CompareVoxel ( float iRAS[3] );
  virtual void DoVoxel ( float iRAS[3] );

  class Params {
  public:
    Params ();
    bool mbStopAtEdges;
    bool mbStopAtROIs;
    bool mb3D;
    bool mbWorkPlaneX;
    bool mbWorkPlaneY;
    bool mbWorkPlaneZ;
    int mFuzziness;
    int mMaxDistance;
    bool mbDiagonal;
  };

  void Flood ( VolumeCollection& iVolume, float iRASSeed[3], Params& iParams );
  VolumeCollection* mVolume;
  Params* mParams;
};


#endif
