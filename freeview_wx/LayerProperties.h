/**
 * @file  LayerProperties.h
 * @brief The common properties
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:39 $
 *    $Revision: 1.1 $
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

#ifndef LayerProperties_h
#define LayerProperties_h

#include "Broadcaster.h"
#include "Listener.h"

class LayerProperties : public Broadcaster, public Listener
{
public:
  LayerProperties ();
  ~LayerProperties ();
  
  void SetShowInfo( bool bShow );
  int  GetShowInfo()
  {
    return m_bShowInfo;
  }

protected:
  bool  m_bShowInfo;
};

#endif
