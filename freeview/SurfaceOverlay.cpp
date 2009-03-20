/**
 * @file  SurfaceOverlay.cxx
 * @brief Implementation for surface layer properties.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/03/20 19:03:54 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2007-2009,
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


#include <assert.h>
#include "SurfaceOverlay.h"
#include "vtkLookupTable.h"
#include "vtkRGBATransferFunction.h"
#include "vtkMath.h"
#include "LayerSurface.h"
#include "SurfaceOverlayProperties.h"
#include "FSSurface.h"
#include <wx/filename.h>

using namespace std;

SurfaceOverlay::SurfaceOverlay ( LayerSurface* surf ) :
    Broadcaster( "SurfaceOverlay" ),
    Listener( "SurfaceOverlay" ),
    m_fData( NULL ),
    m_surface( surf )
{
  InitializeData();  
  
  m_properties =  new SurfaceOverlayProperties( this );
  m_properties->AddListener( this );
}

SurfaceOverlay::~SurfaceOverlay ()
{
  if ( m_fData )
    delete[] m_fData;
  
  delete m_properties;
}

void SurfaceOverlay::DoListenToMessage ( std::string const iMessage, void* iData, void* sender )
{
  if ( iMessage == "ColorMapChanged" )
  {
    this->SendBroadcast( "OverlayChanged", this );
  }
}

void SurfaceOverlay::InitializeData()
{
  if ( m_surface )
  {
    MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
    if ( m_fData )
      delete[] m_fData;
    m_nDataSize = mris->nvertices;
    m_fData = new float[ m_nDataSize ];
    if ( !m_fData )
      return;

    m_dMaxValue = m_dMinValue = mris->vertices[0].val;
    for ( int vno = 0; vno < m_nDataSize; vno++ )
    {
      m_fData[vno] = mris->vertices[vno].val;
      if ( m_dMaxValue < m_fData[vno] )
        m_dMaxValue = m_fData[vno];
      else if ( m_dMinValue > m_fData[vno] )
        m_dMinValue = m_fData[vno];
    }
  }
}

const char* SurfaceOverlay::GetName()
{
  return m_strName.c_str();
}

void SurfaceOverlay::SetName( const char* name )
{
  m_strName = name;
}

void SurfaceOverlay::MapOverlay( unsigned char* colordata )
{
  m_properties->MapOverlayColor( colordata, m_nDataSize );
}



