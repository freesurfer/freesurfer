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

class ScubaLayer2DMRI : public Layer {

  friend class ScubaLayer2DMRITester;

 public:
  ScubaLayer2DMRI ();
  virtual ~ScubaLayer2DMRI ();

  // Associate a volume collection with this layer.
  void SetVolumeCollection ( VolumeCollection& iVolume );
  VolumeCollection* GetVolumeCollection () { return mVolume; }

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
  
  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription () { return "2DMRI"; }

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv );

  virtual void
    DoListenToMessage ( std::string isMessage, void* iData );

  // Return the volume collection.
  virtual DataCollection* GetMainDataCollection() { return mVolume; }

  // Override to change min/max range on new value.
  virtual void DataChanged();

  virtual void HandleTool ( float iRAS[3], ViewState& iViewState,
			    ScubaWindowToRASTranslator& iTranslator,
			    ScubaToolState& iTool, InputState& iInput );

  virtual void ProcessOption ( std::string isOption, std::string isValue );

  enum ColorMapMethod { grayscale, heatScale, LUT };
  void SetColorMapMethod ( ColorMapMethod iMethod );
  ColorMapMethod GetColorMapMethod () { return mColorMapMethod; }
  std::string GetColorMapMethodAsString ();

  enum SampleMethod { nearest, trilinear, sinc, magnitude };
  void SetSampleMethod ( SampleMethod iSampleMethod ) {
    mSampleMethod = iSampleMethod; }
  SampleMethod GetSampleMethod () { return mSampleMethod; }
  std::string GetSampleMethodAsString ();

  void SetColorLUT ( int iLUTID );
  int GetColorLUT ();

  void SetDrawZeroClear ( bool ibClearZero ) { mbClearZero = ibClearZero; }
  bool GetDrawZeroClear () { return mbClearZero; }

  static int const cGrayscaleLUTEntries;
  static int const kMaxPixelComponentValue;  
  static float const kMaxPixelComponentValueFloat;  
  void BuildGrayscaleLUT ();

  // Brightness:  dark  1 --~~==## 0  bright
  // Contrast:    gray  0 --~~==## 30 black/white
  void SetBrightness ( float iBrightness );
  float GetBrightness () { return mBrightness; }
  void SetContrast ( float iContrast );
  float GetContrast () { return mContrast; }

  // Window / level settings.
  void SetWindow ( float iWindow );
  float GetWindow () { return mWindow; }
  void SetLevel ( float iLevel );
  float GetLevel () { return mLevel; }

  void SetMinVisibleValue ( float iValue );
  float GetMinVisibleValue () { return mMinVisibleValue; }
  void SetMaxVisibleValue ( float iValue );
  float GetMaxVisibleValue () { return mMaxVisibleValue; }
  void SetMinMaxVisibleValue ( float iMinValue, float iMaxValue );

  void SetHeatScaleMinThreshold ( float iValue );
  float GetHeatScaleMinThreshold () { return mHeatScaleMinThreshold; }
  void SetHeatScaleMidThreshold ( float iValue );
  float GetHeatScaleMidThreshold () { return mHeatScaleMidThreshold; }
  void SetHeatScaleMaxThreshold ( float iValue );
  float GetHeatScaleMaxThreshold () { return mHeatScaleMaxThreshold; }

  float GetROIOpacity () { return mROIOpacity; }
  void SetROIOpacity ( float iOpacity ) { mROIOpacity = iOpacity; }

  void SetDrawMIP ( bool ibDrawMIP ) { mbDrawMIP = ibDrawMIP; }
  bool GetDrawMIP () { return mbDrawMIP; }

  // Stretch a path from its beginning to the end RAS point.
  void StretchPathStraight  ( Path<float>& iPath,
			      float iRASBegin[3], float iRASEnd[3] );
  void StretchPathAsEdge    ( Path<float>& iPath,
			      float iRASBegin[3], float iRASEnd[3],
			      ViewState& iViewState,
			      ScubaWindowToRASTranslator& iTranslator,
			      float iStraightBias, float iEdgeBias );


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
				 
  virtual void GetPreferredThroughPlaneIncrements ( float oIncrements[3] );

  virtual float GetPreferredBrushRadiusIncrement ();

  virtual float GetPreferredValueIncrement ();

  // For filling.
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

  // For grayscale drawing.
  float mBrightness, mContrast;
  float mOriginalLevel, mOriginalWindow;
  float mWindow, mLevel;
  int mGrayscaleLUT[256]; // 0-255

  // For heatScale drawing.
  float mHeatScaleMinThreshold, mHeatScaleMidThreshold, mHeatScaleMaxThreshold;

  // Our look up table.
  ScubaColorLUT* mColorLUT;
  
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
};

// Flooders ============================================================

class ScubaLayer2DMRIFloodVoxelEdit : public VolumeCollectionFlooder {
 public:
  ScubaLayer2DMRIFloodVoxelEdit ( float iValue );
  ~ScubaLayer2DMRIFloodVoxelEdit () {}

  virtual void DoBegin ();  
  virtual void DoEnd ();
  virtual bool DoStopRequested ();

  virtual bool CompareVoxel ( float iRAS[3] );
  virtual void DoVoxel ( float iRAS[3] );

  float mValue;
};

class ScubaLayer2DMRIFloodSelect : public VolumeCollectionFlooder {
 public:
  ScubaLayer2DMRIFloodSelect ( bool ibSelect );
  ~ScubaLayer2DMRIFloodSelect () {}

  virtual void DoBegin ();  
  virtual void DoEnd ();
  virtual bool DoStopRequested ();

  virtual bool CompareVoxel ( float iRAS[3] );
  virtual void DoVoxel ( float iRAS[3] );

  bool mbSelect;
};

// Undoers ============================================================

class UndoVoxelEditAction : public UndoAction {
 public:

  UndoVoxelEditAction ( VolumeCollection* iVolume, 
			float iNewValue, float iOrigValue, float iRAS[3] );

  virtual void Undo ();
  virtual void Redo ();
  
 protected:
  VolumeCollection* mVolume;
  float mNewValue;
  float mOrigValue;
  float mRAS[3];
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
  EdgePathFinder ( int iViewWidth, int iViewHeight, int iLongestEdge,
		   ScubaWindowToRASTranslator* iTranslator,
		   VolumeCollection* iVolume );
		  
  virtual float GetEdgeCost ( Point2<int>& iPoint );

 protected:
  VolumeCollection* mVolume;
  ScubaWindowToRASTranslator* mTranslator;
};

#endif
