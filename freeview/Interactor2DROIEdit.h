/**
 * @file  Interactor2DROIEdit.h
 * @brief Interactor for editing ROI in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/06/17 20:41:17 $
 *    $Revision: 1.8 $
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

#ifndef Interactor2DROIEdit_h
#define Interactor2DROIEdit_h

#include "Interactor2DVolumeEdit.h"
#include <vector>

class Interactor2DROIEdit : public Interactor2DVolumeEdit
{
public:
  Interactor2DROIEdit();
  virtual ~Interactor2DROIEdit()
  {}
};

#endif


