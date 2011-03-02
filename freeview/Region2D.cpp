/**
 * @file  Region2D.cpp
 * @brief Region2D.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.5 $
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

#include "Region2D.h"
#include "RenderView2D.h"

Region2D::Region2D( RenderView2D* view ) :
    Broadcaster( "Region2D" ),
    Listener( "Region2D" ),
    m_view( view )
{
}

Region2D::~Region2D()
{}

void Region2D::UpdateStats()
{
  SendBroadcast( "RegionStatsUpdated", this );
}

void Region2D::UpdateSlicePosition( int nPlane, double pos )
{
  Update();
}

void Region2D::Show( bool bShow )
{
}

