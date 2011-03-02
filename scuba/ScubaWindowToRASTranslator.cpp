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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.8 $
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
