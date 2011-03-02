/**
 * @file  ScubaLayer2DMRI.h
 * @brief Layer for displaying MRI volumes
 *
 * Displays MRI volumes in grayscale, heatscale, and with an
 * LUT. Handles multiframe volumes and can interpret them as having
 * time points and conditions. Contains the logic for the Edge Path
 * and other path tools.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:38 $
 *    $Revision: 1.64 $
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


#ifndef ScubaLayer2DMRI_h
#define ScubaLayer2DMRI_h

#include "Layer.h"
#include "VolumeCollection.h"
#include "ScubaColorLUT.h"
#include "UndoManager.h"
#include "Timer.h"
#include "ShortestPathFinder.h"
#include "Path.h"
#include "Listener.h"

class EdgePathFinder;

class ScubaLayer2DMRI : public Layer {

  friend class ScubaLayer2DMRITester;

public:
  ScubaLayer2DMRI ();
  virtual ~ScubaLayer2DMRI ();

  // Associate a volume collection with this layer.
  void SetVolumeCollection ( VolumeCollection& iVolume );
  VolumeCollection* GetVolumeCollection () {
    return mVolume;
  }

  // Tell the layer to draw its contents into a GL frame buffer.
  virtual void DrawIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
                                ViewState& iViewState,
                                ScubaWindowToRASTranslator& iTranslator );

  // Tell the layer to draw its contents that need openGl commands.
  virtual void DrawIntoGL ( ViewState& iViewState,
                            ScubaWindowToRASTranslator& iTranslator );

  // Draws the MIP into the buffer.
  virtual void DrawMIPIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
                                   ViewState& iViewState,
                                   ScubaWindowToRASTranslator& iTranslator );

  // These take a volume value and a base color and outputs a final
  // color. The base is used for blending or for translucency if there
  // is no color for the value.
  void GetGrayscaleColorForValue ( float iValue,
                                   GLubyte* const iBase, int* oColor );
  void GetHeatscaleColorForValue ( float iValue,
                                   GLubyte* const iBase, int* oColor );
  void GetColorLUTColorForValue  ( float iValue,
                                   GLubyte* const iBase, int* oColor );

  // We save a cache of RAS coords for the start of each row, and
  // increments as we progress across the row. This cache is the side
  // of the buffer we get. This simply initializes the buffer
  // allocation. The values are filled in in the draw functions
  // depending on view state we get.
  void InitBufferCoordCache ( int iWidth, int iHeight );

  // We override this so we can set our cache of color values.
  virtual void SetOpacity ( float iOpacity );
  void InitColorOpacityCache ();

  // Asks the layer to describe a point of data by making InfoAtRAS
  // structs.
  virtual void GetInfoAtRAS ( float iRAS[3],
                              std::list<InfoAtRAS>& ioInfo );

  // These are the names of the reportable info values we provide.
  enum ReportableInfo { Value = 0, Index, Talairach, kcReportableInfo };
  static char* const kaReportableInfo[kcReportableInfo];

  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription () {
    return "2DMRI";
  }

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv );

  virtual void
  DoListenToMessage ( std::string isMessage, void* iData );

  // Return the volume collection.
  virtual DataCollection* GetMainDataCollection() {
    return mVolume;
  }

  // Override to change min/max range on new value.
  virtual void DataChanged();

  // Handles tool events.
  virtual void HandleTool ( float iRAS[3], ViewState& iViewState,
                            ScubaWindowToRASTranslator& iTranslator,
                            ScubaToolState& iTool, InputState& iInput );

  // Looks for command line options and sets internals from them.
  virtual void ProcessOption ( std::string isOption, std::string isValue );

  // Sets the current frame of the volume to draw.
  int GetCurrentFrame ();
  void SetCurrentFrame ( int iFrame );

  // Get and set the current time point and condition if we're using
  // those instead of straight frame numbers.
  int GetCurrentTimePoint ();
  int GetCurrentCondition ();
  void SetCurrentTimePoint ( int iTimePoint );
  void SetCurrentCondition ( int iCondition );

  // Color map determines how the volume is drawn to screen, and
  // contextually, whether it's an anatomical, functional, or
  // segmenation volume.
  enum ColorMapMethod { grayscale, heatScale, LUT };
  void SetColorMapMethod ( ColorMapMethod iMethod );
  ColorMapMethod GetColorMapMethod () {
    return mColorMapMethod;
  }
  std::string GetColorMapMethodAsString ();

  // Determines the smoothness of voxels on screen.
  enum SampleMethod { nearest, trilinear, sinc, magnitude };
  void SetSampleMethod ( SampleMethod iSampleMethod ) {
    mSampleMethod = iSampleMethod;
  }
  SampleMethod GetSampleMethod () {
    return mSampleMethod;
  }
  std::string GetSampleMethodAsString ();

  // If the color map method is LUT, this determines the LUT to use.
  void SetColorLUT ( int iLUTID );
  int GetColorLUT ();

  // If true, voxels with values of 0 are clear.
  void SetDrawZeroClear ( bool ibClearZero ) {
    mbClearZero = ibClearZero;
  }
  bool GetDrawZeroClear () {
    return mbClearZero;
  }

  // Calcs the grayscale table. cGrayscaleLUTEntries is the resolution
  // of the LUT, kMaxPixelComponentValue is a multiplier to go from a
  // 0->1 range to the e.g. uchar 0-255 range, and
  // kMaxPixelComponentValueFloat is the float version of that.
  static int const cGrayscaleLUTEntries;
  static int const kMaxPixelComponentValue;
  static float const kMaxPixelComponentValueFloat;
  void BuildGrayscaleLUT ();

  // GRAYSCALE COLOR LUT - See https://surfer.nmr.mgh.harvard.edu/fswiki/ScubaGuide_2fScubaWorkingWithData_2fScubaAnatomicalVolumes
  //
  //       |---level---|   level +- range/2
  //             |
  //              .__----               <1 black
  //             /                       1 begin visible black
  //             |                       2->3 black->white
  //             |                       3->4 white
  //             |                       4 end visible white
  //   1   2     /    3   4             >4 black
  //            .
  //   |--------          |
  //  min                max

  // This determines the levels above and below which voxels are drawn
  // as black, or no color.
  void SetMinVisibleValue ( float iValue );
  float GetMinVisibleValue () {
    return mMinVisibleValue;
  }
  void SetMaxVisibleValue ( float iValue );
  float GetMaxVisibleValue () {
    return mMaxVisibleValue;
  }
  void SetMinMaxVisibleValue ( float iMinValue, float iMaxValue );

  // These determine the visible value range. Values below the range
  // are drawn at the 'lowest' color (e.g. black) and values above are
  // drawn in the 'highest' color (e.g. white).
  void SetWindow ( float iWindow );
  float GetWindow () {
    return mWindow;
  }
  void SetLevel ( float iLevel );
  float GetLevel () {
    return mLevel;
  }

  // Brightness/contrast determines the slope of the sigmoid function
  // within the valid window/level range.
  // Brightness:  dark  1 --~~==## 0  bright
  // Contrast:    gray  0 --~~==## 30 black/white
  void SetBrightness ( float iBrightness );
  float GetBrightness () {
    return mBrightness;
  }
  void SetContrast ( float iContrast );
  float GetContrast () {
    return mContrast;
  }

  // These determine the heatscale color map. The threshold is mirrored:
  //
  // -> cyan -> blue -> trans_blue -> clear -> trans_orange -> orange -> red ->
  // -> -max -> -mid ->    -min    ->   0   ->     min      ->   mid  -> max ->
  void SetHeatScaleMinThreshold ( float iValue );
  float GetHeatScaleMinThreshold () {
    return mHeatScaleMinThreshold;
  }
  void SetHeatScaleMidThreshold ( float iValue );
  float GetHeatScaleMidThreshold () {
    return mHeatScaleMidThreshold;
  }
  void SetHeatScaleMaxThreshold ( float iValue );
  float GetHeatScaleMaxThreshold () {
    return mHeatScaleMaxThreshold;
  }
  void SetHeatScaleOffset ( float iValue );
  float GetHeatScaleOffset () {
    return mHeatScaleOffset;
  }

  void SetReverseHeatScale ( bool ib ) {
    mbReverseHeatScale = ib;
  }
  bool GetReverseHeatScale () {
    return mbReverseHeatScale;
  }
  void SetShowPositiveHeatScaleValues ( bool ib ) {
    mbShowPositiveHeatScaleValues = ib;
  }
  bool GetShowPositiveHeatScaleValues () {
    return mbShowPositiveHeatScaleValues;
  }
  void SetShowNegativeHeatScaleValues ( bool ib ) {
    mbShowNegativeHeatScaleValues = ib;
  }
  bool GetShowNegativeHeatScaleValues () {
    return mbShowNegativeHeatScaleValues;
  }

  // Sets the heatscale threshold by FDR value. If asked to use a
  // mask, this will attempt to load the given volume as a mask.
  void SetHeatScaleThresholdUsingFDR ( float       iRate,
                                       std::string ifnMask );

  // Set the volume mask of this layer.
  void SetMaskVolume ( VolumeCollection& iVolume );
  void RemoveMaskVolume ();
  int GetMaskVolume ();

  // Opacity of ROIs drawn on this layer.
  float GetROIOpacity () {
    return mROIOpacity;
  }
  void SetROIOpacity ( float iOpacity ) {
    mROIOpacity = iOpacity;
  }

  // Flag to draw the maximum intensity projection.
  void SetDrawMIP ( bool ibDrawMIP ) {
    mbDrawMIP = ibDrawMIP;
  }
  bool GetDrawMIP () {
    return mbDrawMIP;
  }

  // Stretch a path from its beginning to the end RAS point.
  void StretchPathStraight  ( Path<float>& iPath,
                              float iRASBegin[3], float iRASEnd[3] );
  void StretchPathAsEdge    ( Path<float>& iPath,
                              float iRASBegin[3], float iRASEnd[3],
                              ViewState& iViewState,
                              ScubaWindowToRASTranslator& iTranslator,
                              float iEdgeBias );


  // Select or deselect voxels along a path.
  void SelectVoxelsOnPath ( Path<float>& iPath, bool ibSelect );

  // Finds the closest path in a plane.
  Path<float>* FindClosestPathInPlane ( float iRAS[3],
                                        ViewState& iViewState );

  // Draw a path.
  void DrawRASPathIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
                               int iColor[3], ViewState& iViewState,
                               ScubaWindowToRASTranslator& iTranslator,
                               Path<float>& iRASPath );

  // The distance to travel in plane when the user navigates in and
  // out of the plane.
  virtual void GetPreferredThroughPlaneIncrements ( float oIncrements[3] );

  // A hint for an increment for a brush radius; the GUI will use this
  // granularity for drawing brushes. This is the suggested size so
  // that the user may paint one voxel, increment and paint one more
  // voxel, etc.
  virtual float GetPreferredBrushRadiusIncrement ();

  // A hint for an increment for value in this layer; the GUI will use
  // this when providing the user with a choice such as fuzziness for
  // filling. This should be a reasonable granularity based on the
  // value range of the volume, for example 1 for a volume that goes
  // from 0->255, or 0.01 for a volume that goes from -.1 to .1.
  virtual float GetPreferredValueIncrement ();

  // For filling. Sets the given flooder params based on the given
  // tool and view state.
  void SetFloodParams ( ScubaToolState& iTool, ViewState& iViewState,
                        VolumeCollectionFlooder::Params& ioParams );

protected:

  // Given an RAS point and a radius, calculates the corners of a
  // square on the view plane.
  void CalcRASSquareInViewPlane ( float iRAS[3], float iRadiusRAS,
                                  ScubaWindowToRASTranslator& iTranslator,
                                  ViewState& iViewState,
                                  float oSquareRAS[4][3] );

  // Given an RAS point and a radius, calculates the minimum update
  // rect needed to display changes done to voxels in that
  // rectangle. Adds the rect to the update list.
  void CalcAndAddUpdateSquare ( float iRAS[3], float iRadiusRAS,
                                ScubaWindowToRASTranslator& iTranslator,
                                ViewState& iViewState );

  // This is the time callback.
  virtual void DoTimer ();

  // Manages when we should do autosaves from the timer callback.
  int mTimersSinceLastAutosave;
  static int const kcTimersBetweenAutosaves;

  // The volume data set.
  VolumeCollection* mVolume;

  // General drawing settings.
  SampleMethod mSampleMethod;
  ColorMapMethod mColorMapMethod;
  bool mbClearZero;
  float mMinVisibleValue, mMaxVisibleValue;
  float mOldMinValue, mOldMaxValue;
  int mCurrentFrame;

  // For grayscale drawing.
  float mBrightness, mContrast;
  float mOriginalLevel, mOriginalWindow;
  float mWindow, mLevel;
  int mGrayscaleLUT[256]; // 0-255

  // For heatScale drawing.
  float mHeatScaleMinThreshold, mHeatScaleMidThreshold, mHeatScaleMaxThreshold;
  float mHeatScaleOffset;
  bool mbReverseHeatScale;
  bool mbShowPositiveHeatScaleValues;
  bool mbShowNegativeHeatScaleValues;

  // Our look up table.
  ScubaColorLUT* mColorLUT;

  // Optional layer mask.
  VolumeCollection* mMaskVolume;

  // ROI settings.
  float mROIOpacity;
  bool mbEditableROI;

  // Whether we are drawing MIP.
  bool mbDrawMIP;

  // For editing lines.
  Point3<float> mLastMouseUpRAS, mCurrentMouseRAS;
  bool mbDrawEditingLine;

  // For paths.
  Path<float>*   mCurrentPath;
  Point3<float>  mFirstPathRAS;
  Point3<float>  mLastPathMoveRAS;

  // Remembering the screen increments.
  int mBufferIncSize[2];
  float** mRowStartRAS;
  float** mColIncrementRAS;

  // Cache of color values for an opacity.
  GLubyte mColorTimesOpacity[256];
  GLubyte mColorTimesOneMinusOpacity[256];

  // The undo manager action ID for the current drawing operation.
  int mCurrentDrawingOperationActionID;

  // Keep a pointer to our edge finder. It's initialized when first
  // used.
  EdgePathFinder* mEdgePathFinder;
};

// Flooders ============================================================

class ScubaLayer2DMRIFloodVoxelEdit : public VolumeCollectionFlooder {
public:
  ScubaLayer2DMRIFloodVoxelEdit ( float iValue );
  ~ScubaLayer2DMRIFloodVoxelEdit () {}

  virtual void DoBegin ();
  virtual void DoEnd ();
  virtual bool DoStopRequested ();

  virtual bool CompareVoxel ( float iRAS[3], int iFrame );
  virtual void DoVoxel ( float iRAS[3], int iFrame );

  int mActionListID;
  float mValue;
};

class ScubaLayer2DMRIFloodSelect : public VolumeCollectionFlooder {
public:
  ScubaLayer2DMRIFloodSelect ( bool ibSelect );
  ~ScubaLayer2DMRIFloodSelect () {}

  virtual void DoBegin ();
  virtual void DoEnd ();
  virtual bool DoStopRequested ();

  virtual bool CompareVoxel ( float iRAS[3], int iFrame );
  virtual void DoVoxel ( float iRAS[3], int iFrame );

  int mActionListID;
  bool mbSelect;
};

// Undoers ============================================================

class UndoVoxelEditAction : public UndoAction {
public:

  UndoVoxelEditAction ( VolumeCollection* iVolume,
                        float iNewValue, float iOrigValue,
                        float iRAS[3], int iFrame );

  virtual void Undo ();
  virtual void Redo ();

protected:
  VolumeCollection* mVolume;
  float mNewValue;
  float mOrigValue;
  float mRAS[3];
  int   mFrame;
};

class UndoSelectionAction : public UndoAction {
public:

  UndoSelectionAction ( VolumeCollection* iVolume,
                        bool ibSelect, float iRAS[3] );

  virtual void Undo ();
  virtual void Redo ();

protected:
  VolumeCollection* mVolume;
  bool mbSelect;
  float mRAS[3];
};

class UndoPathAction : public UndoAction {
public:
  UndoPathAction ( Path<float>* iPath );
  virtual ~UndoPathAction ();

  virtual void Undo () {}
  virtual void Redo () {}

protected:
  Path<float>* mPath;
};

class UndoNewPathAction : public UndoPathAction {
public:
  UndoNewPathAction ( Path<float>* iPath );

  virtual void Undo ();
  virtual void Redo ();
};

class UndoDeletePathAction : public UndoPathAction {
public:
  UndoDeletePathAction ( Path<float>* iPath );
  virtual ~UndoDeletePathAction ();

  virtual void Undo ();
  virtual void Redo ();
};

// Edge Finder ============================================================

class EdgePathFinder : public ShortestPathFinder {

public:
  EdgePathFinder ( int iViewWidth, int iViewHeight,
                   ScubaWindowToRASTranslator& iTranslator,
                   VolumeCollection& iVolume );

  virtual float GetEdgeCost ( Point2<int> const& iPoint ) const;

protected:
  VolumeCollection& mVolume;
  ScubaWindowToRASTranslator& mTranslator;
};

#endif
