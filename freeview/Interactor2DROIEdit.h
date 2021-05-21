/**
 * @brief Interactor for editing ROI in 2D render view.
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

#ifndef Interactor2DROIEdit_h
#define Interactor2DROIEdit_h

#include "Interactor2DVolumeEdit.h"

class Interactor2DROIEdit : public Interactor2DVolumeEdit
{
public:
  Interactor2DROIEdit( QObject* parent );
  virtual ~Interactor2DROIEdit()
  {}
};

#endif


