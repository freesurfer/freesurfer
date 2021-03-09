/**
 * @brief Interactor for navigating in 3D render view.
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

#ifndef Interactor3DNavigate_h
#define Interactor3DNavigate_h

#include "Interactor3D.h"

class Interactor3DNavigate : public Interactor3D
{
public:
  Interactor3DNavigate( QObject* parent ) : Interactor3D( parent )
  {}

};

#endif


