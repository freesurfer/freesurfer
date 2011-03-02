/**
 * @file  VolumeCollection.h
 * @brief MRI volume access and utilities
 *
 * Provides access to MRI volume data and volume ROIs. Also contains
 * functions for determining of RAS points are in a square or circular
 * area, used by tools. Can also generate histogram data given a set
 * of points. Can generate average value and standard deviation of
 * points. Can perform value-based fills that set ranges of volumes
 * to new values.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.74 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
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

  // Ctors.
  VolumeLocation ( VolumeCollection const& iVolume,
		   float const iRAS[3] );
  VolumeLocation ( VolumeCollection const& iVolume,
		   int const iIndex[3] );
  VolumeLocation ( VolumeCollection const& iVolume,
		   float const iRAS[3], int iFrame );
  VolumeLocation ( VolumeCollection const& iVolume,
		   int const iIndex[3], int iFrame );
  VolumeLocation ( VolumeLocation const& iLoc );
  ~VolumeLocation () {}

  // Accessors.
  inline int const* Index () const { return mIdxi; }
  inline int Index ( int in ) const { return mIdxi[in]; }
  inline float const* IndexF () const { return mIdxf; }
  inline float IndexF ( int in ) const { return mIdxf[in]; }
  VolumeCollection const& GetVolume () const { return *mVolume; }
  inline int GetFrame () const { return mFrame; }

  // Reinit the location from a different set of coords.
  void SetFromRAS   ( float const iRAS[3] );
  void SetFromIndex ( int const iIndex[3] );
  void SetFrame     ( int in );

  VolumeLocation& operator= ( VolumeLocation const& iLocation );

protected:
  VolumeCollection const* mVolume;
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

  // Makes specialized locations.
  VolumeLocation MakeVolumeLocationFromRAS ( float const iRAS[3] ) const;
  VolumeLocation MakeVolumeLocationFromIndex ( int const iIndex[3] ) const;
  VolumeLocation MakeVolumeLocationFromRAS ( float const iRAS[3], 
					     int iFrame ) const;
  VolumeLocation MakeVolumeLocationFromIndex ( int const iIndex[3],
					       int iFrame ) const;

  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription() const {
    return "Volume";
  }

  // Set the file name for this ROI. This will be used to load the MRI
  // and to save it, but is not actually done until LoadVolume() is
  // called or until the MRI is accessed.
  void SetFileName ( std::string const& ifnMRI );
  std::string const& GetFileName () const;

  // Creates an MRI using an existing data collection as a
  // template. If iType is -1, it will copy the type from the template
  // volume, otherwise it will use the given type. iType should be a
  // valid MRI_ type from mri.h (e.g. MRI_UCHAR).
  void MakeUsingTemplate ( int iCollectionID, int iType );

  // Load an MRI from the currently named volume.
  void LoadVolume ();

  // Get the MRI. If we're loading an MRI, this requires the MRI to be
  // set first.
  MRI const* GetMRI () const;

  // Saves the current volume to a file name. If none is given, uses
  // the set file name.
  void Save () const;
  void Save ( std::string const& ifn ) const;

  // Get the min and max values in the volume.
  float GetMRIMinValue () const;
  float GetMRIMaxValue () const;

  // Get the min and max magnitude values. If these are not computed
  // yet, this could call MakeMagnitudeVolume, so it could take a bit.
  float GetMRIMagnitudeMinValue () const; 
  float GetMRIMagnitudeMaxValue () const;

  // Get the size of the voxels.
  float GetVoxelXSize () const;
  float GetVoxelYSize () const;
  float GetVoxelZSize () const;

  // Returns the data type. This is a valid number from mri.h like
  // MRI_UCHAR,
  int GetDataType () const;

  // For multiframe volumes.
  int GetNumberOfFrames () const;

  // To see whether we should use time points and conditions.
  bool InterpretFramesAsTimePoints () const;

  // Valid if InterpretFramesAsTimePoints returns true.
  int GetNumberOfConditions () const;
  int GetNumberOfTimePoints () const;

  // Convert a time point and condition to a frame number.
  int ConvertConditionAndTimePointToFrame ( int iCondition, 
					    int iTimePoint ) const;
  int ExtractConditionFromFrame ( int iFrame ) const;
  int ExtractTimePointFromFrame ( int iFrame ) const;

  // Get bounds information.
  virtual void GetDataRASBounds ( float oBounds[6] ) const;
  void GetMRIIndexRange ( int oMRIIndexRange[3] ) const;

  // Coordinate conversion.
  void RASToMRIIndex ( float const iRAS[3], int oIndex[3] ) const;
  void RASToMRIIndex ( float const iRAS[3], float oIndex[3] ) const;
  void MRIIndexToRAS ( int const iIndex[3], float oRAS[3] ) const;
  void MRIIndexToRAS ( float const iIndex[3], float oRAS[3] ) const;
  void RASToDataRAS  ( float const iRAS[3], float oDataRAS[3] ) const;
  void DataRASToRAS  ( float const iDataRAS[3], float oRAS[3] ) const;

  // This is based on the internal MRI's linear_transform. We check to
  // see if it's present and convert if so. IsTalTransformPresent()
  // can be used to see if there is a transform present. If there is
  // none, RASToTal will have no effect.
  bool IsTalTransformPresent() const;
  void RASToTal      ( float const iRAS[3], float oTal[3] ) const;

  // Special conversion for compatibility with old edit.dat formats.
  // Calls MRIRASToTkRegRAS and MRITkRegRASToRAS
  void TkRegRASToRAS ( float const iTkRegRAS[3], float oRAS[3] ) const;
  void RASToTkRegRAS ( float const iRAS[3], float oTkRegRAS[3] ) const;

  // Bounds testing.
  bool IsInBounds ( VolumeLocation const& iLoc ) const;

  // Calculates values.
  float GetMRINearestValue   ( VolumeLocation const& iLoc ) const;
  float GetMRITrilinearValue ( VolumeLocation const& iLoc ) const;
  float GetMRISincValue      ( VolumeLocation const& iLoc ) const;
  float GetMRIMagnitudeValue ( VolumeLocation const& iLoc ) const;

  // Sets value.
  void SetMRIValue ( VolumeLocation const& iLoc, float const iValue );

  // Respond to Tcl commands.
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  // Respond to broadcasts.
  virtual void
    DoListenToMessage ( std::string isMessage, void* iData );


  // Select or unselect a voxel in the current ROI. This will work on
  // the currently selected ROI to basically add or remove a point
  // from that ROI.
  void Select   ( VolumeLocation const& iLoc );
  void Unselect ( VolumeLocation const& iLoc );

  // Return whether or not an ROI is at this point, and if so, returns
  // the color. If multiple ROIs are present, blends the color. This
  // looks in our selection volume cache.
  bool IsSelected ( VolumeLocation const& iLoc, int oColor[3] ) const;
  bool IsSelected ( VolumeLocation const& iLoc ) const;

  // Return whether or not an ROI is present other than the one passed
  // in. 
  bool IsOtherRASSelected ( float const iRAS[3], int iThisROIID ) const;


  // Writes an ROI to a label file. If ibUseRealRAS is not on, we'll
  // convert the RAS points we write to TkRegRAS points.
  void WriteROIToLabel ( int iROIID, bool ibUseRealRAS, 
			 std::string const& ifnLabel ) const;

  // Creates a new ROI in this data collection from a label file. If
  // ibUseRealRAS is not on, we'll convert the points in the label
  // from TkRegRAS to RAS.
  int NewROIFromLabel ( bool ibUseRealRAS, std::string const& ifnLabel );

  // Writes all structure ROIs to a segmentation volume file. This
  // creates a volume and whereever one of our structure ROIs has a
  // voxel selected, sets the value in the volume to the structure
  // index of the ROI. If there are voxels selected by more than one
  // ROI, it will have the value of the ROI with the higher ID.
  void WriteROIsToSegmentation ( std::string const& ifnVolume ) const;


  // Sets the data to world transform, basically a user settable
  // display transform.
  virtual void SetDataToWorldTransform ( int iTransformID );

  // Returns the combined world to index transform.
  Matrix44 const& GetWorldToIndexTransform () const;

  // Enable or disable world to index transform. If this is off, no
  // voxel->RAS transformation will be done at all and the volume will
  // be centered on 0,0,0.
  void SetUseWorldToIndexTransform ( bool ibUse );
  bool GetUseWorldToIndexTransform () const;

  // Finds and returns RAS points in a square area. Each RAS will
  // translate to a unique voxel. INPUT POINTS MUST BE IN CLOCKWISE OR
  // COUNTERCLOCKWISE ORDER.
  void FindRASPointsInSquare ( float const iCenter[3],
                               float const iPointA[3], float const iPointB[3],
                               float const iPointC[3], float const iPointD[3],
                               float iMaxDistance, // optional
                               std::list<Point3<float> >& oPoints ) const;
  void FindRASPointsInCircle ( float const iPointA[3], float const iPointB[3],
                               float const iPointC[3], float const iPointD[3],
                               float iMaxDistance, // optional
                               float const iCenter[3], float iRadius,
                               std::list<Point3<float> >& oPoints ) const;

  // Finds and returns RAS points on a segment.
  void FindRASPointsOnSegment ( float const iPointA[3], float const iPointB[3],
                                std::list<Point3<float> >& oPoints ) const;


  // Import and export points to a control point file.
  void ImportControlPoints ( std::string const& ifnControlPoints,
                             std::list<Point3<float> >& oControlPoints) const;
  void ExportControlPoints ( std::string const& ifnControlPoints,
                      std::list<Point3<float> > const& iControlPoints ) const;


  // Return a histogram from the RAS voxels passed in. If an ROI was
  // passed in as well, only count the points that are in the ROI as
  // well as the list of points we got.
  void MakeHistogram ( std::list<Point3<float> > const& iRASPoints, 
		       ScubaROIVolume const* iROI,
		       int icBins,
                       float& oMinBinValue, float& oBinIncrement,
                       std::map<int,int>& oBinCounts ) const;

  // Return a histogram for the entire volume. If an ROI was
  // passed in as well, only count the points that are in the ROI.
  void MakeHistogram ( ScubaROIVolume const* iROI, 
		       int icBins,
                       float iMinThresh, float iMaxThresh,
                       float& oMinBinValue, float& oBinIncrement,
                       std::map<int,int>& oBinCounts ) const;

  // For autosaving. This is kind of unfortunately named is the volume
  // doesn't actually decide when to autosave; it's told to by the
  // layer. But this will save the contents of the volume in the /tmp
  // directory.
  void SetAutoSaveOn ( bool ibAutosave );
  bool GetAutoSaveOn () const;
  bool IsAutosaveDirty () const;
  void AutosaveIfDirty () const;

  // Given an MRI voxel index and a plane in idx coords, returns
  // whether or not the voxel intersects the plane. Tests each of the
  // voxel's edges against the plane. Returns the result of the first
  // intersection.
  VectorOps::IntersectionResult VoxelIntersectsPlane
  ( Point3<int> const& iMRIIndex,
    Point3<float> const& iPlaneIdx, Point3<float> const& iPlaneIdxNormal,
    Point3<float>& oIntersectionIdx ) const;

  // A longer version, this will return a list of all the
  // intersections and their results. The above function uses this
  // one. Returns true if there was an intersection.
  bool VoxelIntersectsPlane
  ( Point3<int> const& iMRIIndex,
    Point3<float> const& iPlaneIdx,
    Point3<float> const& iPlaneIdxNormal,
    int& oNumIntersections,
    VectorOps::IntersectionResult orIntersection[12],
    Point3<float> oIntersectionIdx[12] ) const;

  // Does something similar but tests if a segment goes through a
  // voxel. Does this by testing the segment against each face as a
  // plane.
  VectorOps::IntersectionResult VoxelIntersectsSegment
  ( Point3<int> const& iMRIIndex,
    Point3<float> const& iSegIdxA, Point3<float> const& iSegIdxB,
    Point3<float>& oIntersectionIdx ) const;

  // Returns a list of volume locations whose voxels have the given value.
  void GetVoxelsWithValue ( float iValue,
                            std::list<VolumeLocation>& olLocations) const;

  // Returns the volume of N voxels given the voxel size.
  float GetRASVolumeOfNVoxels ( int icVoxels ) const;

  // Returns a best guess value increment for a GUI.
  virtual float GetPreferredValueIncrement () const;

  // Returne the average value of the input locations or ROI.
  float GetAverageValue ( std::list<VolumeLocation> const& ilLocations ) const;
  float GetAverageValue ( ScubaROIVolume const& iROI ) const;
  float GetStandardDeviation ( std::list<VolumeLocation> const& ilLocations, 
			       float iMean ) const;
  float GetStandardDeviation ( ScubaROIVolume const& iROI, float iMean ) const;


  // Print out the corner RAS coordinates of a voxel.
  void PrintVoxelCornerCoords ( std::ostream& iStream,
                                Point3<int> const& iMRIIdx ) const;


  // A mini class to be used with the value range fill.
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
  void DoValueRangeFill ( VolumeCollection const& iSourceVol,
			  ScubaROIVolume const* iROI,
			  std::vector<ValueRangeFillElement> const& lElements);
protected:

  // Gets information from the MRI structure and initializes all
  // internals.
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
  float GetMRINearestValueAtIndexUnsafe ( int const iIndex[3], 
					  int inFrame ) const;

  // Faster way of checking bounds.
  bool IsMRIIdxInBounds ( int const iMRIIdx[3] ) const;
  bool IsMRIIdxInBounds ( float const iMRIIdx[3] ) const;
  bool IsFrameInBounds ( int const iFrame ) const;

  // Create an ROI for this data collection.
  virtual ScubaROI* DoNewROI ();

  // Initialize our selection volume. This is volume of booleans the
  // size of our MRI that indicates whether any ROIs have that voxel
  // selected.
  void InitSelectionVolume ();

  // Called when data is changed. We mark our autosave dirty flag.
  virtual void DataChanged ();

  // Constructs an autosave file name based on a given name. Converts
  // the slashes in a file name into dots and prepends /tmp and
  // appends .mgz.
  std::string MakeAutosaveFileName ( std::string const& ifn ) const;

  // Deletes the current autosave volume.
  void DeleteAutosave();

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
  mutable bool mbAutosaveDirty;
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
  Volume3<bool>* mSelectedVoxels;

  // Bounds cache.
  mutable bool mbBoundsCacheDirty;
  mutable float mRASBounds[6];

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
