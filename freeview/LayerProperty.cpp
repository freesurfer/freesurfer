/**
 * @file  LayerProperty.cpp
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
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:58 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2007 - 2009,
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
