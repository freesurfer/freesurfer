#ifndef ScubaLayer2DMRI_h
#define ScubaLayer2DMRI_h

#include "Layer.h"
#include "VolumeCollection.h"
#include "ScubaColorLUT.h"

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

  virtual void HandleTool ( float iRAS[3],
			    ScubaToolState& iTool, InputState& iInput );

  enum ColorMapMethod { grayscale, heatScale, LUT };
  void SetColorMapMethod ( ColorMapMethod iMethod ) { 
    mColorMapMethod = iMethod; }
  ColorMapMethod GetColorMapMethod () { return mColorMapMethod; }

  enum SampleMethod { nearest, trilinear, sinc };
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

 protected:
  VolumeCollection* mVolume;
  
  SampleMethod mSampleMethod;
  ColorMapMethod mColorMapMethod;

  float mBrightness, mContrast;
  std::map<int,float> mGrayscaleLUT; // 0-255

  ScubaColorLUT* mColorLUT;
  
  bool mbClearZero;
  float mMinVisibleValue, mMaxVisibleValue;
};


#endif
