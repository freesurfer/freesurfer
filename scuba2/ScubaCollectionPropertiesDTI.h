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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.4 $
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
  virtual vtkFSVolumeSource* GetEigenVectorVolumeSource( const int iSource ) const = 0;
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

 protected:

  ScubaCollectionPropertiesDTI () {};
  virtual ~ScubaCollectionPropertiesDTI () {};

};

#endif
