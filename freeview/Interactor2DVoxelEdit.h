/**
 * @file  Interactor2DVoxelEdit.h
 * @brief Interactor for editing voxel in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.11 $
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

#ifndef Interactor2DVoxelEdit_h
#define Interactor2DVoxelEdit_h

#include "Interactor2DVolumeEdit.h"

class wxWindow;

class Interactor2DVoxelEdit : public Interactor2DVolumeEdit
{
public:
  Interactor2DVoxelEdit();
  virtual ~Interactor2DVoxelEdit()
  {}
};

#endif


