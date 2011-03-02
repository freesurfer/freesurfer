/**
 * @file  vtkKWScubaLayerCollectionMRI.h
 * @brief Implementation for MRI viewers.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.3 $
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

#ifndef vtkKWScubaLayerCollectionMRI_h
#define vtkKWScubaLayerCollectionMRI_h

#include "vtkKWScubaLayerCollection.h"
#include "ScubaCollectionPropertiesMRI.h"
#include "vtkSmartPointer.h"
extern "C" {
#include "colortab.h"
}

class vtkFSVolumeSource;
class vtkFreesurferLookupTable;
class vtkKWFrame;
class vtkKWLabel;
class vtkKWRGBATransferFunctionEditor;
class vtkKWRadioButton;
class vtkKWRange;
class vtkRGBATransferFunction;

class vtkKWScubaLayerCollectionMRI : public vtkKWScubaLayerCollection
                                     //BTX
                                     , public ScubaCollectionPropertiesMRI
				     //ETX
{

 public:

  static vtkKWScubaLayerCollectionMRI* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayerCollectionMRI, 
			vtkKWScubaLayerCollection );

  // Description:
  // Set the volume file name.
  void SetVolumeFileName ( const char* ifnVolume );

  // Description:
  // Populate a UI page with common controls for this layer type.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();

  // Settings that are shared by multiple layer types.

  // Layers will get the color table objects with these calls. These
  // impelement ScubaCollectionPropertiesMRI.
  vtkFSVolumeSource* GetSource () const;
  vtkFreesurferLookupTable* GetLUTTable () const;
  vtkRGBATransferFunction* GetGrayScaleTable () const;
  vtkRGBATransferFunction* GetHeatScaleTable () const;
  COLOR_TABLE* GetLUTCTAB () const;
  
  // Description:
  // Set the color map appropriately.
  //BTX
  void SetColorMap ( ColorMapType iType );
  ColorMapType GetColorMap () const; // Implements ScubaCollectionPropertiesMRI
  //ETX
  void SetColorMapToGrayScale ();
  void SetColorMapToHeatScale ();
  void SetColorMapToLUT ();
  void ColorMapChanged ();

  // Description.
  // Set the reslice object's interpolation.
  void SetResliceInterpolation ( int iMode );
  int GetResliceInterpolation () const; // Implements ScubaCollectionPropertiesMRI

  // Description.
  // Set the texture's interpolation (we call it smoothing to
  // differentiate from the reslice interpolation).
  void SetTextureSmoothing ( int iSmooth );
  int GetTextureSmoothing () const; // Implements ScubaCollectionPropertiesMRI

  // Description:
  // This determines the levels above and below which voxels are drawn
  // as black, or no color.
  void SetMinVisibleValue ( float iValue );
  float GetMinVisibleValue ();
  void SetMaxVisibleValue ( float iValue );
  float GetMaxVisibleValue ();
  void SetMinMaxVisibleValue ( float iMinValue, float iMaxValue );

  // Description:
  // These determine the visible value range. Values below the range
  // are drawn at the 'lowest' color (e.g. black) and values above are
  // drawn in the 'highest' color (e.g. white).
  void SetWindow ( float iWindow );
  float GetWindow ();
  void SetLevel ( float iLevel );
  float GetLevel ();
  void SetMinMaxGrayscaleWindow ( float iMin, float iMax );
  void SetMinGrayscaleWindow ( float iMin );
  void SetMaxGrayscaleWindow ( float iMax );

  // Description:
  // These determine the heatscale color map. The threshold is mirrored:
  // -> cyan -> blue -> trans_blue -> clear -> trans_orange -> orange -> red ->
  // -> -max -> -mid ->    -min    ->   0   ->     min      ->   mid  -> max ->
  void SetHeatScaleMinThreshold ( float iValue );
  float GetHeatScaleMinThreshold ();
  void SetHeatScaleMidThreshold ( float iValue );
  float GetHeatScaleMidThreshold ();
  void SetHeatScaleMaxThreshold ( float iValue );
  float GetHeatScaleMaxThreshold ();
  void SetHeatScaleOffset ( float iValue );
  float GetHeatScaleOffset ();

  // Description:
  // Various heat scale characterists.
  void SetReverseHeatScale ( bool ib );
  bool GetReverseHeatScale ();
  void SetShowPositiveHeatScaleValues ( bool ib );
  bool GetShowPositiveHeatScaleValues ();
  void SetShowNegativeHeatScaleValues ( bool ib );
  bool GetShowNegativeHeatScaleValues ();

  void EditorChangedHeatScale ();

  // Description:
  // Pop up a dialog box asking for a file name, and pass the file
  // name to LoadFreesurferLUT.
  void LoadLUTFromDlog ();

  // Description:
  // Load the LUT into the LUT color map.
  void LoadFreesurferLUT ( const char* ifn );


 protected:

  vtkKWScubaLayerCollectionMRI ();
  ~vtkKWScubaLayerCollectionMRI ();

  //BTX
  virtual vtkKWScubaLayer*
    MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode );
  //ETX

  void LoadVolumeFromFileName ();

 private:

  //BTX

  // Color tables --------------------------------------------------------
  vtkSmartPointer<vtkFreesurferLookupTable> mLUTTable;
  vtkSmartPointer<vtkRGBATransferFunction> mGrayScaleTable;
  vtkSmartPointer<vtkRGBATransferFunction> mHeatScaleTable;
  // ---------------------------------------------------------------------

  // Controls ------------------------------------------------------------
  vtkSmartPointer<vtkKWRadioButton> maRadBtnColorMap[cColorMapTypes];
  vtkSmartPointer<vtkKWRange> mRangeVisibleValue;
  vtkSmartPointer<vtkKWRange> mRangeWindowLevel;
  vtkSmartPointer<vtkKWFrame> maFrameColorMapSettings[cColorMapTypes];
  vtkSmartPointer<vtkKWRGBATransferFunctionEditor> mRGBAEditorHeatScale;
  vtkSmartPointer<vtkKWLabel> mLabelLUTFileName;
  // ---------------------------------------------------------------------

  // Color map variables -------------------------------------------------
  ColorMapType mColorMapType;

  // Applicable to all types.
  int mResliceInterpolation;
  int mTextureSmoothing;

  bool mbClearZero;
  float mMinVisibleValue, mMaxVisibleValue;

  // For grayscale drawing.
  float mMinGrayscaleWindow, mMaxGrayscaleWindow;

  // For heatScale drawing.
  float mHeatScaleMinThreshold, mHeatScaleMidThreshold, mHeatScaleMaxThreshold;
  float mHeatScaleOffset;
  bool mbReverseHeatScale;
  bool mbShowPositiveHeatScaleValues;
  bool mbShowNegativeHeatScaleValues;

  // LUT drawing.
  COLOR_TABLE* mFreeSurferCTAB;

  // ---------------------------------------------------------------------

  vtkSmartPointer<vtkFSVolumeSource> mSource;
  std::string mfnVolume;
  //ETX

};

#endif
