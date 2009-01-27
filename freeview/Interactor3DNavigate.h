/**
 * @file  Interactor3DNavigate.h
 * @brief Interactor3DNavigate to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/01/27 18:43:47 $
 *    $Revision: 1.2.2.2 $
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

#ifndef Interactor3DNavigate_h
#define Interactor3DNavigate_h

#include "Interactor3D.h"

class Interactor3DNavigate : public Interactor3D
{
public:
  Interactor3DNavigate() : Interactor3D()
  {}

};

#endif


