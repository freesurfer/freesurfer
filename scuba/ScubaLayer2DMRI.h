#ifndef ScubaLayer2DMRI_h
#define ScubaLayer2DMRI_h

#include "Layer.h"
#include "VolumeCollection.h"
#include "ScubaColorLUT.h"
#include "UndoManager.h"
#include "Timer.h"
#include "ShortestPathFinder.h"
#include "PointList3.h"

class ScubaLayer2DMRI : public Layer {

  friend class ScubaLayer2DMRITester;

 public:
  ScubaLayer2DMRI ();
  virtual ~ScubaLayer2DMRI ();

  // Associate a volume collection with this layer.
  void SetVolumeCollection ( VolumeCollection& iVolume );

  // Tell the layer to draw its contents into a GL frame buffer.
  virtual void DrawIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
				ViewState& iViewState,
				ScubaWindowToRASTranslator& iTranslator );
  
  // Asks the layer to describe a point of data by adding pairs of
  // labels and values.
  virtual void GetInfoAtRAS ( float iRAS[3],
			   std::map<std::string,std::string>& iLabelValues );
  
  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription () { return "2DMRI"; }

  virtual TclCommandResult
    DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv );

  virtual void HandleTool ( float iRAS[3], ViewState& iViewState,
			    ScubaWindowToRASTranslator& iTranslator,
			    ScubaToolState& iTool, InputState& iInput );

  enum ColorMapMethod { grayscale, heatScale, LUT };
  void SetColorMapMethod ( ColorMapMethod iMethod ) { 
    mColorMapMethod = iMethod; }
  ColorMapMethod GetColorMapMethod () { return mColorMapMethod; }

  enum SampleMethod { nearest, trilinear, sinc, magnitude };
  void SetSampleMethod ( SampleMethod iSampleMethod ) {
    mSampleMethod = iSampleMethod; }
  SampleMethod GetSampleMethod () { return mSampleMethod; }

  void SetColorLUT ( int iLUTID );

  static int const cGrayscaleLUTEntries;
  static int const kMaxPixelComponentValue;  
  void BuildGrayscaleLUT ();
  void SetBrightness ( float iBrightness ) { mBrightness = iBrightness; }
  void SetContrast ( float iContrast ) { mContrast = iContrast; }

  void SetMinVisibleValue ( float iValue ) { mMinVisibleValue = iValue; }
  float GetMinVisibleValue () { return mMinVisibleValue; }
  void SetMaxVisibleValue ( float iValue ) { mMaxVisibleValue = iValue; }
  float GetMaxVisibleValue () { return mMaxVisibleValue; }

  float GetROIOpacity () { return mROIOpacity; }
  void SetROIOpacity ( float iOpacity ) { mROIOpacity = iOpacity; }

  // Creates a new line, adds it to the list of lines, and returns it.
  PointList3<float>* NewLine ();

  // Stretch a line from its beginning to the end RAS point.
  void StretchLineStraight  ( PointList3<float>& iLine,
			      float iRASBegin[3], float iRASEnd[3] );
  void StretchLineAsEdge    ( PointList3<float>& iLine,
			      float iRASBegin[3], float iRASEnd[3],
			      ViewState& iViewState,
			      ScubaWindowToRASTranslator& iTranslator );

  // Ends a line. Adds points to volume's edge volume.
  void EndLine              ( PointList3<float>& iLine,
			      ScubaWindowToRASTranslator& iTranslator );

  // Returns the line closest to the RAS point.
  PointList3<float>*
    FindClosestLine         ( float iRAS[3],
			      ViewState& iViewState );

  void DrawRASPointListIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
				    int iColor[3], ViewState& iViewState,
				    ScubaWindowToRASTranslator& iTranslator,
				    PointList3<float>& iRASList );
				 
  virtual void GetPreferredInPlaneIncrements ( float oIncrements[3] );
 
 protected:

  VolumeCollection* mVolume;
  
  SampleMethod mSampleMethod;
  ColorMapMethod mColorMapMethod;

  float mBrightness, mContrast;
  //  std::map<int,float> mGrayscaleLUT; // 0-255
  float mGrayscaleLUT[256]; // 0-255

  ScubaColorLUT* mColorLUT;
  
  bool mbClearZero;
  float mMinVisibleValue, mMaxVisibleValue;

  float mROIOpacity;

  // For lines.
  std::list<PointList3<float>* > mLines;
  PointList3<float>*             mCurrentLine;
  Point3<float>                  mFirstLineRAS;
  Point3<float>                  mLastLineMoveRAS;

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

  Timer mFloodTimer;
  bool mbFloodDlogOpen;
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
