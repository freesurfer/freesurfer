/**
 * @file  vtkSimpleLabelEdgeFilter.h
 * @brief A simple label edge filter ONLY for 2D label image (2D vtkImageData).
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:03 $
 *    $Revision: 1.2 $
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

class VTK_IMAGING_EXPORT vtkSimpleLabelEdgeFilter : public vtkSimpleImageToImageFilter
{
  public:
    static vtkSimpleLabelEdgeFilter *New();
    vtkTypeRevisionMacro(vtkSimpleLabelEdgeFilter,vtkSimpleImageToImageFilter);

  protected:

    vtkSimpleLabelEdgeFilter() {};
    ~vtkSimpleLabelEdgeFilter() {};

    virtual void SimpleExecute(vtkImageData* input, vtkImageData* output);
  private:
    vtkSimpleLabelEdgeFilter(const vtkSimpleLabelEdgeFilter&);  // Not implemented.
    void operator=(const vtkSimpleLabelEdgeFilter&);  // Not implemented.
};

#endif
