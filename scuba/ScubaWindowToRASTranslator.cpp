/**
 * @file  ScubaWindowToRASTranslator.cpp
 * @brief Interface for a window to RAS translator
 *
 * Implemented by ScubaView, this is an interface for Layers to
 * transform between Window and RAS coords.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/22 04:39:29 $
 *    $Revision: 1.7 $
 *
 * Copyright (C) 2002-2007,
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


#include <iostream>
#include "ScubaWindowToRASTranslator.h"

using namespace std;

void
ScubaWindowToRASTranslator::TranslateWindowToRAS( int const[2], float[3] ) {
  cerr << "NO" << endl;
}

void
ScubaWindowToRASTranslator::TranslateRASToWindow( float const[3], int[2] ) {
  cerr << "NO" << endl;
}
