#ifndef VolumeCollection_H
#define VolumeCollection_H

#include "string_fixed.h"
#include "DataCollection.h"
#include "ScubaROIVolume.h"
extern "C" {
#include "mri.h"
#ifdef X
  #undef X 
#endif
#ifdef Y
  #undef Y
#endif
#ifdef Z
  #undef Z
#endif
}
#include "Volume3.h"
#include "Point3.h"
#include "ScubaTransform.h"
#include "Broadcaster.h"

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

  float GetMRIMagnitudeMinValue ();
  float GetMRIMagnitudeMaxValue ();

  float GetVoxelXSize () { return mVoxelSize[0]; }
  float GetVoxelYSize () { return mVoxelSize[1]; }
  float GetVoxelZSize () { return mVoxelSize[2]; }

  void RASToMRIIndex ( float iRAS[3], int oIndex[3] );
  void RASToMRIIndex ( float iRAS[3], float oIndex[3] );
  void MRIIndexToRAS ( int iIndex[3], float oRAS[3] );
  void MRIIndexToRAS ( float iIndex[3], float oRAS[3] );
  void RASToDataRAS  ( float iRAS[3], float oDataRAS[3] );

  bool IsRASInMRIBounds ( float iRAS[3] );
  bool IsMRIIndexInMRIBounds ( int iIndex[3] );

  float GetMRINearestValueAtRAS ( float iRAS[3] );
  float GetMRITrilinearValueAtRAS ( float iRAS[3] );
  float GetMRISincValueAtRAS ( float iRAS[3] );
  
  void SetMRIValueAtRAS ( float iRAS[3], float iValue );
  
  void MakeMagnitudeVolume ();

  float GetMRIMagnitudeValueAtRAS ( float iRAS[3] );

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  virtual void
    DoListenToMessage ( std::string isMessage, void* iData );

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

  void GetRASPointsInCube ( float iCenterRAS[3], float iRadius,
			    bool ibBrushX, bool ibBrushY, bool ibBrushZ,
			    std::list<Point3<float> >& oPoints );
  void GetRASPointsInSphere ( float iCenterRAS[3], float iRadius,
			      bool ibBrushX, bool ibBrushY, bool ibBrushZ,
			      std::list<Point3<float> >& oPoints );

  void WriteROIToLabel ( int iROIID, std::string ifnLabel );
  int  NewROIFromLabel ( std::string ifnLabel );
  void WriteROIsToSegmentation ( std::string ifnVolume );

  virtual void SetDataToWorldTransform ( int iTransformID ) {
    DataCollection::SetDataToWorldTransform( iTransformID );
    CalcWorldToIndexTransform();
  }

protected:

  void CalcWorldToIndexTransform ();

  std::string mfnMRI;
  MRI* mMRI;
  MRI* mMagnitudeMRI;

  // We're dealing with two transforms here, data->world which is
  // defined in DataCollection, and data->index which is defined
  // here. Actually it goes world->data->index and
  // index->data->world. We composite these in the last transform.
  Transform44 mDataToIndexTransform;
  Transform44 mWorldToIndexTransform;

  Volume3<Point3<int> >* mDataToIndexCache;
  void CalcDataToIndexCache ();
  inline void DataToIndexCacheIndex ( float const iRAS[3],
				       int oCacheIndex[3] ) const;

  float mVoxelSize[3];
  float mOneOverVoxelSize[3];

  float mMRIMinValue, mMRIMaxValue;
  float mMRIMagMinValue, mMRIMagMaxValue;

  Volume3<bool>* mEdgeVoxels;

  Volume3<bool>* mSelectedVoxels;
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
