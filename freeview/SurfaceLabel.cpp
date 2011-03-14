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
 *    $Author: nicks $
 *    $Date: 2011/03/14 23:44:48 $
 *    $Revision: 1.6 $
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
 *
 */


#include "SurfaceLabel.h"
#include "LayerSurface.h"
#include "FSSurface.h"
#include <QFile>

SurfaceLabel::SurfaceLabel ( LayerSurface* surf ) :
  QObject( surf ),
  m_label( NULL ),
  m_surface( surf )
{
  m_rgbColor[0] = 1.0;
  m_rgbColor[1] = 1.0;
  m_rgbColor[2] = 0.0;

  //this->AddListener( surf );
}

SurfaceLabel::~SurfaceLabel ()
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }
}

/*
void SurfaceLabel::DoListenToMessage ( std::string const iMessage, void* iData, void* sender )
{
  if ( iMessage == "ColorMapChanged" )
  {
    this->SendBroadcast( "OverlayChanged", this );
  }
}
*/

QString SurfaceLabel::GetName()
{
  return m_strName;
}

void SurfaceLabel::SetName( const QString& name )
{
  m_strName = name;
}

bool SurfaceLabel::LoadLabel( const QString& filename )
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }

  m_label = ::LabelRead( NULL, filename.toAscii().data() );

  if ( m_label == NULL )
  {
    cerr << "LabelRead failed";
    return false;
  }

  QFile file( filename );
  if ( !file.open( QIODevice::ReadOnly | QIODevice::Text ) )
  {
    return false;
  }

  QString strg;
  while (!file.atEnd())
  {
    strg += file.readLine();
  }

  if ( strg.contains( "vox2ras=", Qt::CaseInsensitive ) &&
       !strg.contains( "vox2ras=TkReg", Qt::CaseInsensitive ) )
  {
    m_bTkReg = false;
  }

  return true;
}

void SurfaceLabel::SetColor( double r, double g, double b )
{
  m_rgbColor[0] = r;
  m_rgbColor[1] = g;
  m_rgbColor[2] = b;

  emit SurfaceLabelChanged();
}

void SurfaceLabel::MapLabel( unsigned char* colordata, int nVertexCount )
{
  if ( !m_label )
  {
    return;
  }

  for ( int i = 0; i < m_label->n_points; i++ )
  {
    int vno = m_label->lv[i].vno;
    if ( vno < nVertexCount )
    {
      double opacity = m_label->lv[i].stat;
      if ( opacity > 1 )
      {
        opacity = 1;
      }
      else if ( opacity < 0 )
      {
        opacity = 0;
      }
      opacity = 1;    // ignore opacity for now
      colordata[vno*4]    = ( int )( colordata[vno*4]   * ( 1 - opacity ) + m_rgbColor[0] * 255 * opacity );
      colordata[vno*4+1]  = ( int )( colordata[vno*4+1] * ( 1 - opacity ) + m_rgbColor[1] * 255 * opacity );
      colordata[vno*4+2]  = ( int )( colordata[vno*4+2] * ( 1 - opacity ) + m_rgbColor[2] * 255 * opacity );
    }
  }
}
