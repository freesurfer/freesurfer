/**
 * @file  vtkKWScubaLayerCollectionDTI.h
 * @brief Implementation for DTI viewers.
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:05 $
 *    $Revision: 1.1 $
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

#ifndef vtkKWScubaLayerCollectionDTI_h
#define vtkKWScubaLayerCollectionDTI_h

#include "vtkKWScubaLayerCollection.h"
#include "ScubaCollectionPropertiesDTI.h"

class vtkKWCheckButton;
class vtkFSVolumeSource;

class vtkKWScubaLayerCollectionDTI : public vtkKWScubaLayerCollection
                                     //BTX
                                     , public ScubaCollectionPropertiesDTI
				     //ETX
{

 public:

  static vtkKWScubaLayerCollectionDTI* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayerCollectionDTI, 
			vtkKWScubaLayerCollection );

  // Description:
  // Set the file name for the DTI DA volume. We'll use this name to
  // infer the names of the other volumes.
  void SetFAVolumeFileName ( const char* ifnVolume );

  // Description:
  // Populate a UI page with common controls for this layer type.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();

  // Settings that are shared by multiple layer types.

  vtkImageAppendComponents* GetMergedSource () const;
  vtkFSVolumeSource* GetFAVolumeSource() const;
  vtkFSVolumeSource* GetEigenValueVolumeSource() const;
  vtkFSVolumeSource* GetEigenVector1VolumeSource() const;
  vtkFSVolumeSource* GetEigenVector2VolumeSource() const;
  vtkFSVolumeSource* GetEigenVector3VolumeSource() const;

  // Description.
  void SetRenderEdges ( int ibRender );
  int GetRenderEdges () const; // Implements ScubaCollectionPropertiesDTI

  // Description.
  // Set tensor interpolation.  Could be VTK_RESLICE_LINEAR,
  // VTK_RESLICE_CUBIC, or VTK_RESLICE_NEAREST.
  void SetTensorInterpolationType ( int iMode );
  int GetTensorInterpolationType () const; // Implements ScubaCollectionPropertiesDTI

  // Description:
  // The scaling type for the tensor glyphs. Could be Uniform or FA.
  void SetTensorScalingFromInt ( int iScaling );
  //BTX
  void SetTensorScaling ( ScubaCollectionPropertiesDTI::TensorScalingType iScaling );
  ScubaCollectionPropertiesDTI::TensorScalingType GetTensorScaling () const; // Implements ScubaCollectionPropertiesDTI
  //ETX

  // Decription
  // The detail level for the tensor glyphs. Could be Least, Less, or
  // Nomral
  void SetTensorDetailFromInt ( int iDetail );
  //BTX
  void SetTensorDetail ( ScubaCollectionPropertiesDTI::TensorDetailType iDetail );
  ScubaCollectionPropertiesDTI::TensorDetailType GetTensorDetail () const; // Implements ScubaCollectionPropertiesDTI
  //ETX

 protected:

  vtkKWScubaLayerCollectionDTI ();
  ~vtkKWScubaLayerCollectionDTI ();

  //BTX
  virtual vtkKWScubaLayer*
    MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode );
  //ETX

  void LoadVolumesFromFileName ();

 private:

  //BTX

  // Common VTK pipeline elements ----------------------------------------
  vtkFSVolumeSource* mFAVolumeSource;
  vtkFSVolumeSource* mEV1VolumeSource;
  vtkFSVolumeSource* mEV2VolumeSource;
  vtkFSVolumeSource* mEV3VolumeSource;
  vtkFSVolumeSource* mEValuesVolumeSource;
  vtkImageAppendComponents* mMergedSources;
  // ---------------------------------------------------------------------

  // Controls ------------------------------------------------------------
  vtkKWCheckButton* mChkBtnRenderEdges;
  vtkKWRadioButton* mRadBtnScalingUniform;
  vtkKWRadioButton* mRadBtnScalingFA;
  vtkKWRadioButton* mRadBtnDetailLeast;
  vtkKWRadioButton* mRadBtnDetailLess;
  vtkKWRadioButton* mRadBtnDetailNormal;
  // ---------------------------------------------------------------------

  // Applicable to all types.
  int mTensorInterpolation;
  int mbRenderEdges;
  TensorScalingType mScaling;
  TensorDetailType mDetail;

  // ---------------------------------------------------------------------

  std::string mfnFAVolume;
  //ETX

};

#endif
