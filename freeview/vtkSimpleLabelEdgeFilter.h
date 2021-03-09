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
// .NAME vtkSimpleLabelEdgeFilter - Simple example of an image-image filter.
// .SECTION Description
// This is an example of a simple image-image filter. It copies it's input
// to it's output (point by point). It shows how templates can be used
// to support various data types.
// .SECTION See also
// vtkSimpleImageToImageFilter

#ifndef __vtkSimpleLabelEdgeFilter_h
#define __vtkSimpleLabelEdgeFilter_h

#include "vtkSimpleImageToImageFilter.h"

#if VTK_MAJOR_VERSION > 5
#include "vtkImagingGeneralModule.h" // For export macro

class VTKIMAGINGGENERAL_EXPORT vtkSimpleLabelEdgeFilter : public vtkSimpleImageToImageFilter
{
public:
  static vtkSimpleLabelEdgeFilter *New();
  vtkTypeMacro(vtkSimpleLabelEdgeFilter,vtkSimpleImageToImageFilter);

protected:
  vtkSimpleLabelEdgeFilter() {};
  ~vtkSimpleLabelEdgeFilter() override {};

  void SimpleExecute(vtkImageData* input, vtkImageData* output) override;
private:
  vtkSimpleLabelEdgeFilter(const vtkSimpleLabelEdgeFilter&) = delete;  // Not implemented.
  void operator=(const vtkSimpleLabelEdgeFilter&) = delete;  // Not implemented.
};
#else
class vtkSimpleLabelEdgeFilter : public vtkSimpleImageToImageFilter
{
public:
  static vtkSimpleLabelEdgeFilter *New();
  vtkTypeMacro(vtkSimpleLabelEdgeFilter,vtkSimpleImageToImageFilter);

protected:
  vtkSimpleLabelEdgeFilter() {};
  ~vtkSimpleLabelEdgeFilter() {};

  virtual void SimpleExecute(vtkImageData* input, vtkImageData* output);
private:
  vtkSimpleLabelEdgeFilter(const vtkSimpleLabelEdgeFilter&);  // Not implemented.
  void operator=(const vtkSimpleLabelEdgeFilter&);  // Not implemented.
};
#endif

#endif
