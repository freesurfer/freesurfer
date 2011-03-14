/**
 * @file  SurfaceRegionGroups.cpp
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:59 $
 *    $Revision: 1.6 $
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

#include "SurfaceRegionGroups.h"
#include "LayerMRI.h"
#include "SurfaceRegion.h"

SurfaceRegionGroups::SurfaceRegionGroups( LayerMRI* owner ) : 
    QObject( owner ),
    m_mri( owner )
{
  m_colors << Qt::blue
      << Qt::red
      << Qt::darkRed
      << Qt::green
      << Qt::darkGreen
      << Qt::blue
      << Qt::darkBlue
      << Qt::yellow
      << Qt::darkYellow
      << Qt::magenta
      << Qt::darkMagenta
      << Qt::cyan
      << Qt::darkCyan;
}

SurfaceRegionGroups::~SurfaceRegionGroups()
{}

void SurfaceRegionGroups::SetGroupColor( int n, const QColor& color )
{
  if ( m_colors.size() > n-1 )
    m_colors[n-1] = color;
  else
    m_colors.push_back( color );
  
  QList<SurfaceRegion*>& regs = m_mri->m_surfaceRegions;
  for ( int i = 0; i < regs.size(); i++ )
  {
    if ( regs[i]->GetGroup() == n )
      regs[i]->SetColor( color );
  }
}

QColor SurfaceRegionGroups::GetGroupColor( int n )
{
  if ( m_colors.size() > n-1 )
    return m_colors[n-1];
  else
    return QColor();
}

int SurfaceRegionGroups::GetGroupIdRange( SurfaceRegion* reg )
{
  QList<SurfaceRegion*>& regs = m_mri->m_surfaceRegions;
  int nMax = 0;
  for ( int i = 0; i < regs.size(); i++ )
  {
    if ( regs[i] != reg && regs[i]->GetGroup() > nMax )
      nMax = regs[i]->GetGroup();
  }
  return (nMax+1);
}
