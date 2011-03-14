/**
 * @file  Interactor2DVoxelEdit.h
 * @brief Interactor for editing voxel in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:57 $
 *    $Revision: 1.14 $
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

#ifndef Interactor2DVoxelEdit_h
#define Interactor2DVoxelEdit_h

#include "Interactor2DVolumeEdit.h"

class Interactor2DVoxelEdit : public Interactor2DVolumeEdit
{
public:
  Interactor2DVoxelEdit( QObject* parent );
  virtual ~Interactor2DVoxelEdit()
  {}
};

#endif


