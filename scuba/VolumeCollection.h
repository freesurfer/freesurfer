/**
 * @file  VolumeCollection.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/15 20:42:21 $
 *    $Revision: 1.70 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef VolumeCollection_H
#define VolumeCollection_H

#include "string_fixed.h"
#include "DataCollection.h"
#include "ScubaROIVolume.h"
extern "C" {
#ifdef Status
#undef Status
#endif
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
#include "VectorOps.h"

class VolumeCollection;

class VolumeLocation : public DataLocation {
  friend class VolumeCollectionTester;
  friend class VolumeCollection;
public:
  VolumeLocation ( VolumeCollection& iVolume, float const iRAS[3] );
  VolumeLocation ( VolumeCollection& iVolume, int const iIndex[3] );
  VolumeLocation ( const VolumeLocation& iLoc );
  ~VolumeLocation () {}
  int*   Index      ()               {
    return mIdxi;
  }
  int    Index      ( int in ) const {
    return mIdxi[in];
  }
  float* IndexF     ()               {
    return mIdxf;
  }
  float  IndexF     ( int in ) const {
    return mIdxf[in];
  }
  void   SetFromRAS   ( float const iRAS[3] );
  void   SetFromIndex ( int const iIndex[3] );
  void   SetFrame     ( int in );
  VolumeCollection* GetVolume () const {
    return mVolume;
  }
protected:
  VolumeCollection* mVolume;
  int mFrame;
  int mIdxi[3];
  float mIdxf[3];
};

class VolumeCollection : public DataCollection {

  friend class VolumeCollectionTester;
  friend class VolumeCollectionFlooder;

public:
  VolumeCollection ();
  virtual ~VolumeCollection ();

  // Sets the frame in the location to 0.
  virtual DataLocation& MakeLocationFromRAS ( float const iRAS[3] );
  DataLocation& MakeLocationFromIndex ( int const iIndex[3] );

  // Optional way that also sets the frame.
  DataLocation& MakeLocationFromRAS ( float const iRAS[3], int iFrame );
  DataLocation& MakeLocationFromIndex ( int const iIndex[3], int iFrame );

  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription() {
    return "Volume";
  }

  // Set the file name for this ROI. This will be used to load the MRI
  // and to save it.
  void SetFileName ( std::string& ifnMRI );
  std::string GetFileName () {
    return mfnMRI;
  }

  // Creates an MRI using an existing data collection as a
  // template. If iType is -1, it will copy the type from the template
  // volume, otherwise it will use the given type. iType should be a
  // valid MRI_ type from mri.h.
  void MakeUsingTemplate ( int iCollectionID, int iType );

  // Load an MRI from the currently named volume.
  void LoadVolume ();

  // Get the MRI. If we're loading an MRI, this requires the MRI to be
  // set first.
  MRI const* GetMRI () const;

  // Saves the current volume to a file name. If none is given, uses
  // the set file name.
  void Save ();
  void Save ( std::string ifn );

  // Accessors.
  float GetMRIMinValue () {
    return mMRIMinValue;
  }
  float GetMRIMaxValue () {
    return mMRIMaxValue;
  }

  float GetMRIMagnitudeMinValue (); // (Calls MakeMagnitudeVolume if necessary)
  float GetMRIMagnitudeMaxValue ();

  float GetVoxelXSize () {
    return mVoxelSize[0];
  }
  float GetVoxelYSize () {
    return mVoxelSize[1];
  }
  float GetVoxelZSize () {
    return mVoxelSize[2];
  }

  // For multiframe volumes.
  int GetNumberOfFrames () {
    return mcFrames;
  }

  // To see whether we should use time points and conditions.
  bool InterpretFramesAsTimePoints () {
    return mbInterpretFramesAsTimePoints;
  }
  int GetNumberOfConditions () {
    return mcConditions;
  }
  int GetNumberOfTimePoints () {
    return mcTimePoints;
  }

  // Convert a time point and condition to a frame number.
  int ConvertConditionAndTimePointToFrame ( int iCondition, int iTimePoint );
  int ExtractConditionFromFrame ( int iFrame );
  int ExtractTimePointFromFrame ( int iFrame );

  virtual void GetDataRASBounds ( float oBounds[6] );
  void GetMRIIndexRange ( int oMRIIndexRange[3] );

  // Coordinate conversion.
  void RASToMRIIndex ( float const iRAS[3], int oIndex[3] );
  void RASToMRIIndex ( float const iRAS[3], float oIndex[3] );
  void MRIIndexToRAS ( int const iIndex[3], float oRAS[3] );
  void MRIIndexToRAS ( float const iIndex[3], float oRAS[3] );
  void RASToDataRAS  ( float const iRAS[3], float oDataRAS[3] );
  void DataRASToRAS  ( float const iDataRAS[3], float oRAS[3] );

  // This is based on the internal MRI's linear_transform. We check to
  // see if it's present and convert if so. IsTalTransformPresent()
  // can be used to see if there is a transform present. If there is
  // none, RASToTal will have no effect.
  bool IsTalTransformPresent();
  void RASToTal      ( float const iRAS[3], float oTal[3] );

  // Special conversion for compatibility with old edit.dat formats.
  // Calls MRIRASToTkRegRAS and MRITkRegRASToRAS
  void TkRegRASToRAS ( float const iTkRegRAS[3], float oRAS[3] );
  void RASToTkRegRAS ( float const iRAS[3], float oTkRegRAS[3] );

  // Bounds testing.
  bool IsInBounds ( VolumeLocation& iLoc ) const;

  // Calculates values.
  float GetMRINearestValue   ( VolumeLocation& iLoc ) const;
  float GetMRITrilinearValue ( VolumeLocation& iLoc ) const;
  float GetMRISincValue      ( VolumeLocation& iLoc ) const;
  float GetMRIMagnitudeValue ( VolumeLocation& iLoc );

  // Sets value.
  void SetMRIValue ( VolumeLocation& iLoc, float const iValue );

  // Returns the data type. This is a valid number from mri.h like
  // MRI_UCHAR,
  int GetDataType ();

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  virtual void
  DoListenToMessage ( std::string isMessage, void* iData );

  // Create an ROI for this data collection.
  virtual ScubaROI* DoNewROI ();


  // Select or unselect in the current ROI.
  void Select   ( VolumeLocation& iLoc );
  void Unselect ( VolumeLocation& iLoc );

  // Return whether or not an ROI is at this point, and if so, returns
  // the color. If multiple ROIs are present, blends the color.
  bool IsSelected ( VolumeLocation& iLoc, int oColor[3] );
  bool IsSelected ( VolumeLocation& iLoc );

  // Return whether or not an ROI is present other than the one passed
  // in.
  bool IsOtherRASSelected ( float iRAS[3], int iThisROIID );


  // Writes an ROI to a label file.
  void WriteROIToLabel ( int iROIID, bool ibUseRealRAS, std::string ifnLabel );

  // Creates a new ROI in this data collection from a label file.
  int  NewROIFromLabel ( bool ibUseRealRAS, std::string ifnLabel );

  // Writes all structure ROIs to a segmentation volume file.
  void WriteROIsToSegmentation ( std::string ifnVolume );


  // Sets the data to world transform, basically a user settable
  // display transform.
  virtual void SetDataToWorldTransform ( int iTransformID );

  // Returns the combined world to index transform.
  Matrix44 const& GetWorldToIndexTransform () const;

  // Enable or disable world to index transform.
  void SetUseWorldToIndexTransform ( bool ibUse );
  bool GetUseWorldToIndexTransform () {
    return mbUseDataToIndexTransform;
  }

  // Finds and returns RAS points in a square area. Each RAS will
  // translate to a unique voxel. INPUT POINTS MUST BE IN CLOCKWISE OR
  // COUNTERCLOCKWISE ORDER.
  void FindRASPointsInSquare ( float iCenter[3],
                               float iPointA[3], float iPointB[3],
                               float iPointC[3], float iPointD[3],
                               float iMaxDistance, // optional
                               std::list<Point3<float> >& oPoints );
  void FindRASPointsInCircle ( float iPointA[3], float iPointB[3],
                               float iPointC[3], float iPointD[3],
                               float iMaxDistance, // optional
                               float iCenter[3], float iRadius,
                               std::list<Point3<float> >& oPoints );

  // Finds and returns RAS points on a segment.
  void FindRASPointsOnSegment ( float iPointA[3], float iPointB[3],
                                std::list<Point3<float> >& oPoints );


  // Import and export points to a control point file.
  void ImportControlPoints ( std::string ifnControlPoints,
                             std::list<Point3<float> >& oControlPoints);
  void ExportControlPoints ( std::string ifnControlPoints,
                             std::list<Point3<float> >& iControlPoints );


  // Return a histogram from the RAS voxels passed in.
  void MakeHistogram ( std::list<Point3<float> >& iRASPoints, 
		       ScubaROIVolume* iROI,
		       int icBins,
                       float& oMinBinValue, float& oBinIncrement,
                       std::map<int,int>& oBinCounts );

  // Return a histogram for the entire volume.
  void MakeHistogram ( ScubaROIVolume* iROI, 
		       int icBins,
                       float iMinThresh, float iMaxThresh,
                       float& oMinBinValue, float& oBinIncrement,
                       std::map<int,int>& oBinCounts );

  virtual void DataChanged ();

  // For autosaving.
  void SetAutoSaveOn ( bool ibAutosave ) {
    mbAutosave = ibAutosave;
  }
  bool GetAutoSaveOn () {
    return mbAutosave;
  }

  bool IsAutosaveDirty () {
    return mbAutosaveDirty;
  }
  void AutosaveIfDirty ();

  // Given an MRI voxel index and a plane in idx coords, returns
  // whether or not the voxel intersects the plane. Tests each of the
  // voxel's edges against the plane. Returns the result of the first
  // intersection.
  VectorOps::IntersectionResult VoxelIntersectsPlane
  ( Point3<int>& iMRIIndex,
    Point3<float>& iPlaneIdx, Point3<float>& iPlaneIdxNormal,
    Point3<float>& oIntersectionIdx );

  // A longer version, this will return a list of all the
  // intersections and their results. The above function uses this
  // one. Returns true if there was an intersection.
  bool VoxelIntersectsPlane
  ( Point3<int>& iMRIIndex,
    Point3<float>& iPlaneIdx,
    Point3<float>& iPlaneIdxNormal,
    int& oNumIntersections,
    VectorOps::IntersectionResult orIntersection[12],
    Point3<float> oIntersectionIdx[12] );

  // Does something similar but tests if a segment goes through a
  // voxel. Does this by testing the segment against each face as a
  // plane.
  VectorOps::IntersectionResult VoxelIntersectsSegment
  ( Point3<int>& iMRIIndex,
    Point3<float>& iSegIdxA, Point3<float>& iSegIdxB,
    Point3<float>& oIntersectionIdx );

  // Returns a list of volume locations whose voxels have the given value.
  void GetVoxelsWithValue ( float iValue,
                            std::list<VolumeLocation>& olLocations);

  // Returns the volume of N voxels.
  float GetRASVolumeOfNVoxels ( int icVoxels );

  // Returns a best guess value increment for a GUI.
  virtual float GetPreferredValueIncrement ();

  // Returne the average value of the input locations.
  float GetAverageValue ( std::list<VolumeLocation>& ilLocations );
  float GetAverageValue ( ScubaROIVolume& iROI );
  float GetStandardDeviation ( std::list<VolumeLocation>& ilLocations, float iMean );
  float GetStandardDeviation ( ScubaROIVolume& iROI, float iMean );


  // Print out the corner RAS coordinates of a voxel.
  void PrintVoxelCornerCoords ( std::ostream& iStream,
                                Point3<int>& iMRIIdx );

  class ValueRangeFillElement {
  public:
    ValueRangeFillElement ( float iBegin, float iEnd, float iValue ) :
        mBegin(iBegin), mEnd(iEnd), mValue(iValue) {}
    float mBegin, mEnd, mValue;
  };

  // For each voxel, look at the corresponding voxel in the source
  // vol, see if it falls in one of the ranges in the list, and is in
  // an ROI if passed in, and set the value to the one specified for
  // that range.
  void DoValueRangeFill ( VolumeCollection& iSourceVol,
			  ScubaROIVolume* iROI,
			  std::vector<ValueRangeFillElement>& lElements);
protected:

  // Gets information from the MRI structure.
  void InitializeFromMRI ();

  // Look for time metadata and parse it if available.
  void TryReadingTimeMetadata ();

  // Update the value range.
  void UpdateMRIValueRange ();

  // Calculates the final world to index transfrom from the data to
  // index transform and the world to index transform.
  void CalcWorldToIndexTransform ();

  // Calculates the magnitude volume and saves it in mMagnitudeMRI.
  void MakeMagnitudeVolume ();

  // Faster way of getting values.
  float GetMRINearestValueAtIndexUnsafe ( int iIndex[3], int inFrame );

  // Faster way of checking bounds.
  inline bool IsMRIIdxInBounds ( int const iMRIIdx[3] ) const;
  inline bool IsMRIIdxInBounds ( float const iMRIIdx[3] ) const;
  inline bool IsFrameInBounds ( int const iFrame ) const;

  // Filename to use when reading or writing.
  std::string mfnMRI;

  // The actual MRI volume.
  MRI* mMRI;

  // The magnitude volume, used for the edge line tool. This is
  // generated by MakeMagnitudeVolume whenever a magnitude value is
  // set.
  MRI* mMagnitudeMRI;

  // We're dealing with two transforms here, data->world which is
  // defined in DataCollection, and data->index which is defined
  // here. Actually it goes world->data->index and
  // index->data->world. We composite these in the last transform.
  Transform44 mDataToIndexTransform;
  Transform44 mWorldToIndexTransform;

  bool mbUseDataToIndexTransform;

  // For the autosave volume.
  bool mbAutosave;
  std::string MakeAutosaveFileName ( std::string& ifn );
  void DeleteAutosave();
  bool mbAutosaveDirty;
  std::string mfnAutosave;

  // Number of frames.
  int mcFrames;

  // Time point and condition info if we're to use them. num time
  // points * num conditions = num frames.
  bool mbInterpretFramesAsTimePoints;
  int  mcConditions;
  int  mcTimePoints;

  // These are for display purposes in a time course graph. We use the
  // time resolution to draw the x axis as seconds, and use the num of
  // pre stim time points to draw a stimulus marker at the right
  // time. Only valid if mbInterpretFramesAsTimePoints is true.
  int   mcPreStimTimePoints;
  float mTimeResolution;

  // If there was a certain kind of header information available, some
  // conditions and time points should be interpreted as error
  // data. We also hold a matrix of covariances, size numTimePoints *
  // (numConditions-1) squared.
  bool    mbDataContainsErrorValues;
  typedef std::map<int,std::map<int,float> > CovarianceTable;
  CovarianceTable mCovarianceTable;

  // Voxel sizes.
  float mVoxelSize[3];

  // MRI and magnitude volume ranges.
  float mMRIMinValue, mMRIMaxValue;
  float mMRIMagMinValue, mMRIMagMaxValue;

  // The selected voxel cache.
  void InitSelectionVolume ();
  Volume3<bool>* mSelectedVoxels;

  // Bounds cache.
  bool mbBoundsCacheDirty;
  float mRASBounds[6];

  // Static transforms for doing MNI tal -> Real Tal conversion, used
  // in RASToTal.
  static Matrix44 mMatrixMNITalToTalGtz; // use when z > 0
  static Matrix44 mMatrixMNITalToTalLtz; // use when z <= 0
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

  virtual bool CompareVoxel ( float iRAS[3], int iFrame );
  virtual void DoVoxel ( float iRAS[3], int iFrame );

  class Params {
  public:
    Params ();
    int mSourceCollection;
    bool mbStopAtPaths;
    bool mbStopAtROIs;
    bool mb3D;
    bool mbWorkPlaneX;
    bool mbWorkPlaneY;
    bool mbWorkPlaneZ;
    float mViewNormal[3];
    float mFuzziness;
    float mMaxDistance;
    bool mbDiagonal;
    bool mbOnlyZero;
    enum FuzzinessType { seed, gradient };
    FuzzinessType mFuzzinessType;
  };

  class CheckPair {
  public:
    CheckPair() {}
    CheckPair( Point3<int>& iFrom, Point3<int>& iTo ) {
      mFromIndex = iFrom;
      mToIndex = iTo;
    }
    Point3<int> mFromIndex;
    Point3<int> mToIndex;
  };

  void Flood ( VolumeCollection& iVolume,
               float iRASSeed[3], int iFrame,
               Params& iParams );
  VolumeCollection* mVolume;
  Params* mParams;
};


#endif
