/**
 * @brief Implementation for generic layer properties.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
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


#include "LayerProperty.h"

LayerProperty::LayerProperty ( QObject* parent ) : QObject( parent ),
  m_bShowInfo( true )
{
}

LayerProperty::~LayerProperty ()
{
}

void LayerProperty::SetShowInfo ( bool bShowInfo )
{
  if ( m_bShowInfo != bShowInfo )
  {
    m_bShowInfo = bShowInfo;

    emit ShowInfoChanged( bShowInfo );
    emit PropertyChanged();
  }
}
