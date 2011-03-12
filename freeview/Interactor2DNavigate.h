/**
 * @file  Interactor2DNavigate.h
 * @brief Interactor for navigating (zoom, pan, etc.) in 2D render view.
 * 
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:48 $
 *    $Revision: 1.7 $
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

#ifndef Interactor2DNavigate_h
#define Interactor2DNavigate_h

#include "Interactor2D.h"

class Interactor2DNavigate : public Interactor2D
{
public:
  Interactor2DNavigate( QObject* parent ) : Interactor2D( parent )
  {}

};

#endif


