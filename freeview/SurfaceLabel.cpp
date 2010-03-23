/**
 * @file  SurfaceLabel.cxx
 * @brief Implementation for surface label.
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
 *    $Date: 2010/03/23 18:31:10 $
 *    $Revision: 1.1 $
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


#include <wx/wx.h>
#include <wx/ffile.h>
#include "SurfaceLabel.h"
#include "LayerSurface.h"
#include "FSSurface.h"
#include <wx/filename.h>

SurfaceLabel::SurfaceLabel ( LayerSurface* surf ) :
    Broadcaster( "SurfaceLabel" ),
    Listener( "SurfaceLabel" ),
    m_label( NULL ),
    m_surface( surf )
{
  m_rgbColor[0] = 1.0;
  m_rgbColor[1] = 1.0;
  m_rgbColor[2] = 0.0;
  
  this->AddListener( surf );
}

SurfaceLabel::~SurfaceLabel ()
{
  if ( m_label )
    ::LabelFree( &m_label );
}

void SurfaceLabel::DoListenToMessage ( std::string const iMessage, void* iData, void* sender )
{
  if ( iMessage == "ColorMapChanged" )
  {
    this->SendBroadcast( "OverlayChanged", this );
  }
}

const char* SurfaceLabel::GetName()
{
  return m_strName.c_str();
}

void SurfaceLabel::SetName( const char* name )
{
  m_strName = name;
}

bool SurfaceLabel::LoadLabel( const char* filename )
{
  if ( m_label )
    ::LabelFree( &m_label );

  char* fn = strdup( filename );
  m_label = ::LabelRead( NULL, fn );
  free( fn );

  if ( m_label == NULL )
  {
    cerr << "LabelRead failed";
    return false;
  }

  wxFFile file( wxString::FromAscii(filename) );
  wxString strg;
  if ( file.ReadAll( &strg ) )
  {
    if ( strg.Find( _("vox2ras=") ) >= 0 && 
         strg.Find( _("vox2ras=TkReg") ) < 0 )
      m_bTkReg = false;
  }

  return true;
}

void SurfaceLabel::SetColor( double r, double g, double b )
{
  m_rgbColor[0] = r;
  m_rgbColor[1] = g;
  m_rgbColor[2] = b;
  
  this->SendBroadcast( "SurfaceLabelChanged", this );
}

void SurfaceLabel::MapLabel( unsigned char* colordata, int nVertexCount )
{
  if ( !m_label )
    return;
  
  for ( int i = 0; i < m_label->n_points; i++ )
  {
    int vno = m_label->lv[i].vno;
    if ( vno < nVertexCount )
    {
      double opacity = m_label->lv[i].stat;
      if ( opacity > 1 )
        opacity = 1;
      else if ( opacity < 0 )
        opacity = 0;
      opacity = 1;    // ignore opacity for now
      colordata[vno*4]    = ( int )( colordata[vno*4]   * ( 1 - opacity ) + m_rgbColor[0] * 255 * opacity ); 
      colordata[vno*4+1]  = ( int )( colordata[vno*4+1] * ( 1 - opacity ) + m_rgbColor[1] * 255 * opacity ); 
      colordata[vno*4+2]  = ( int )( colordata[vno*4+2] * ( 1 - opacity ) + m_rgbColor[2] * 255 * opacity ); 
    }
  }
}
