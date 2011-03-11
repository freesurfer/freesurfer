/**
 * @file  LayerProperties.cpp
 * @brief Implementation for generic layer properties.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
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


#include "LayerProperties.h"

LayerProperties::LayerProperties () :
    Broadcaster( "LayerProperties" ),
    Listener( "LayerProperties" )
{
  m_bShowInfo = true;
}

LayerProperties::~LayerProperties ()
{
}

void LayerProperties::SetShowInfo ( bool bShowInfo )
{
  if ( m_bShowInfo != bShowInfo )
  {
    m_bShowInfo = bShowInfo;
    
    this->SendBroadcast( "ShowInfoChanged", NULL );
  }
}
