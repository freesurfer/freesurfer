/**
 * @file  vtkSimpleLabelEdgeFilter.h
 * @brief A simple label edge filter ONLY for 2D label image (2D vtkImageData).
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:55 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2008-2009,
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
