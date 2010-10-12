/**
 * @file  SurfaceRegionGroups.cpp
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/10/12 21:22:32 $
 *    $Revision: 1.1 $
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
    Broadcaster( "SurfaceRegionGroups" ),
    Listener( "SurfaceRegionGroups" ),
    m_mri( owner )
{
  m_colors.push_back( wxColour( "blue" ) );
  m_colors.push_back( wxColour( "red" ) );
  m_colors.push_back( wxColour( "yellow" ) );
  m_colors.push_back( wxColour( "MAGENTA" ) );
  m_colors.push_back( wxColour( "CYAN" ) );
  m_colors.push_back( wxColour( "PURPLE" ) );
  m_colors.push_back( wxColour( "ORANGE" ) );
  m_colors.push_back( wxColour( "SEA GREEN" ) );
  m_colors.push_back( wxColour( "SLATE BLUE" ) );
  m_colors.push_back( wxColour( "PINK" ) );
}

SurfaceRegionGroups::~SurfaceRegionGroups()
{}

void SurfaceRegionGroups::SetGroupColor( int n, const wxColour& color )
{
  if ( (int)m_colors.size() > n-1 )
    m_colors[n-1] = color;
  else
    m_colors.push_back( color );
  
  std::vector<SurfaceRegion*>& regs = m_mri->m_surfaceRegions;
  for ( size_t i = 0; i < regs.size(); i++ )
  {
    if ( regs[i]->GetGroup() == n )
      regs[i]->SetColor( color );
  }
}

wxColour SurfaceRegionGroups::GetGroupColor( int n )
{
  if ( (int)m_colors.size() > n-1 )
    return m_colors[n-1];
  else
    return wxColour();
}

int SurfaceRegionGroups::GetGroupIdRange( SurfaceRegion* reg )
{
  std::vector<SurfaceRegion*>& regs = m_mri->m_surfaceRegions;
  int nMax = 0;
  for ( size_t i = 0; i < regs.size(); i++ )
  {
    if ( regs[i] != reg && regs[i]->GetGroup() > nMax )
      nMax = regs[i]->GetGroup();
  }
  return (nMax+1);
}
