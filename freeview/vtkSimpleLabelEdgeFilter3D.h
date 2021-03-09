/**
 * @brief A simple label edge filter ONLY for 2D label image (2D vtkImageData).
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 *
 */
// .NAME vtkSimpleLabelEdgeFilter3D - Simple example of an image-image filter.
// .SECTION Description
// This is an example of a simple image-image filter. It copies it's input
// to it's output (point by point). It shows how templates can be used
// to support various data types.
// .SECTION See also
// vtkSimpleImageToImageFilter

#ifndef __vtkSimpleLabelEdgeFilter3D_h
#define __vtkSimpleLabelEdgeFilter3D_h

#include "vtkSimpleImageToImageFilter.h"

#if VTK_MAJOR_VERSION > 5
#include "vtkImagingGeneralModule.h" // For export macro

class VTKIMAGINGGENERAL_EXPORT vtkSimpleLabelEdgeFilter3D : public vtkSimpleImageToImageFilter
{
public:
  static vtkSimpleLabelEdgeFilter3D *New();
  vtkTypeMacro(vtkSimpleLabelEdgeFilter3D,vtkSimpleImageToImageFilter);

protected:
  vtkSimpleLabelEdgeFilter3D() {};
  ~vtkSimpleLabelEdgeFilter3D() override {};

  void SimpleExecute(vtkImageData* input, vtkImageData* output) override;
private:
  vtkSimpleLabelEdgeFilter3D(const vtkSimpleLabelEdgeFilter3D&) = delete;  // Not implemented.
  void operator=(const vtkSimpleLabelEdgeFilter3D&) = delete;  // Not implemented.
};
#else
class vtkSimpleLabelEdgeFilter3D : public vtkSimpleImageToImageFilter
{
public:
  static vtkSimpleLabelEdgeFilter3D *New();
  vtkTypeMacro(vtkSimpleLabelEdgeFilter3D,vtkSimpleImageToImageFilter);

protected:
  vtkSimpleLabelEdgeFilter3D() {};
  ~vtkSimpleLabelEdgeFilter3D() {};

  virtual void SimpleExecute(vtkImageData* input, vtkImageData* output);
private:
  vtkSimpleLabelEdgeFilter3D(const vtkSimpleLabelEdgeFilter3D&);  // Not implemented.
  void operator=(const vtkSimpleLabelEdgeFilter3D&);  // Not implemented.
};
#endif

#endif
