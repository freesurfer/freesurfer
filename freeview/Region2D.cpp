/**
 * @file  Region2D.cpp
 * @brief Region2D.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:53 $
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

#include "Region2D.h"
#include "RenderView2D.h"

Region2D::Region2D( RenderView2D* view ) :
    QObject( view ),
    m_view( view )
{
}

Region2D::~Region2D()
{}

void Region2D::UpdateStats()
{
  emit StatsUpdated();
}

void Region2D::UpdateSlicePosition( int nPlane, double pos )
{
  Update();
}

void Region2D::Show( bool bShow )
{
}

