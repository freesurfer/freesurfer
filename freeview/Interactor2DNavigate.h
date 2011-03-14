/**
 * @file  Interactor2DNavigate.h
 * @brief Interactor for navigating (zoom, pan, etc.) in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/14 23:44:47 $
 *    $Revision: 1.10 $
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


