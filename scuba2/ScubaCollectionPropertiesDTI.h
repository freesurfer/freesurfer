/**
 * @file  ScubaCollectionPropertiesDTI.h
 * @brief The common properties available to DTI layers
 *
 * An abstract interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:04 $
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

#ifndef ScubaCollectionPropertiesDTI_h
#define ScubaCollectionPropertiesDTI_h

class vtkImageAppendComponents;
class vtkFSVolumeSource;

class ScubaCollectionPropertiesDTI {

 public:

  // Description:
  // Get a pointer to the various source volumes.
  virtual vtkImageAppendComponents* GetMergedSource () const = 0;
  virtual vtkFSVolumeSource* GetFAVolumeSource() const = 0;
  virtual vtkFSVolumeSource* GetEigenValueVolumeSource() const = 0;
  virtual vtkFSVolumeSource* GetEigenVector1VolumeSource() const = 0;
  virtual vtkFSVolumeSource* GetEigenVector2VolumeSource() const = 0;
  virtual vtkFSVolumeSource* GetEigenVector3VolumeSource() const = 0;

  // Description:
  // Get the flag for rendering edges.
  virtual int GetRenderEdges () const = 0;
  
  // Description:
  // The interpolation type for the tensor FA volume.  Returns
  // VTK_RESLICE_LINEAR, VTK_RESLICE_CUBIC, or VTK_RESLICE_NEAREST
  virtual int GetTensorInterpolationType () const = 0;

  // Description:
  // The scaling type for the tensor glyphs.
  enum TensorScalingType {
    Uniform, FA
  };
  virtual TensorScalingType GetTensorScaling () const = 0;

  // Decription
  // The detail level for the tensor glyphs.
  enum TensorDetailType {
    Least, Less, Normal
  };
  virtual TensorDetailType GetTensorDetail () const = 0;

};

#endif
