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
 *    $Date: 2015/01/06 20:46:12 $
 *    $Revision: 1.14 $
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
#include "LayerPropertySurface.h"
#include "FSSurface.h"
#include "MyUtils.h"
#include <QFile>
#include <QDebug>

SurfaceLabel::SurfaceLabel ( LayerSurface* surf ) :
  QObject( surf ),
  m_label( NULL ),
  m_surface( surf ),
  m_bVisible( true ),
  m_nOutlineIndices(NULL),
  m_dThreshold(0)
{
  m_rgbColor[0] = 1.0;
  m_rgbColor[1] = 1.0;
  m_rgbColor[2] = 0.0;

  SetShowOutline(false);
  SetColor(1.0, 1.0, 0);
}

SurfaceLabel::~SurfaceLabel ()
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }
  if (m_nOutlineIndices)
    delete[] m_nOutlineIndices;
}


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

  // update vno if it is -1
  MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
  MHT* hash = MHTfillVertexTableRes(mris, NULL,
                                    CURRENT_VERTICES, 16);
  for (int i = 0; i < m_label->n_points; i++)
  {
    if (m_label->lv[i].vno < 0)
    {
      VERTEX v;
      v.x = m_label->lv[i].x;
      v.y = m_label->lv[i].y;
      v.z = m_label->lv[i].z;
      float dmin;
      int vtxno = MHTfindClosestVertexNo(hash, mris, &v, &dmin);
      if (vtxno >= 0)
        m_label->lv[i].vno = vtxno;
    }
  }

  // create outline
  m_nOutlineIndices = new int[m_label->n_points];
  UpdateOutline();

  return true;
}

void SurfaceLabel::UpdateOutline()
{
  VERTEX *v;
  MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
  MRISclearMarks(mris);
  LabelMarkSurface(m_label, mris);
  for (int n = 0 ; n < m_label->n_points ; n++)
  {
    m_nOutlineIndices[n] = 0;
    if (m_label->lv[n].vno >= 0) // && m_label->lv[n].stat > m_dThreshold)
    {
      v = &mris->vertices[m_label->lv[n].vno] ;
      if (v->ripflag)
        continue;

      for (int m = 0 ; m < v->vnum ; m++)
      {
        if (mris->vertices[v->v[m]].marked == 0)
        {
          m_nOutlineIndices[n] = 1;
          break;
        }
      }
    }
  }
}

void SurfaceLabel::SetColor( double r, double g, double b )
{
  m_rgbColor[0] = r;
  m_rgbColor[1] = g;
  m_rgbColor[2] = b;

//  emit SurfaceLabelChanged();
}

void SurfaceLabel::SetThreshold(double th)
{
  if (th != m_dThreshold)
  {
    m_dThreshold = th;
    emit SurfaceLabelChanged();
  }
}

void SurfaceLabel::MapLabel( unsigned char* colordata, int nVertexCount )
{
  if ( !m_label)
  {
    return;
  }

  for ( int i = 0; i < m_label->n_points; i++ )
  {
    int vno = m_label->lv[i].vno;
    if ( vno < nVertexCount && (!m_bShowOutline || m_nOutlineIndices[i] > 0) )
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
      if (m_label->lv[i].stat > m_dThreshold)
        opacity = 1;
      else
        opacity = 0;
      colordata[vno*4]    = ( int )( colordata[vno*4]   * ( 1 - opacity ) + m_rgbColor[0] * 255 * opacity );
      colordata[vno*4+1]  = ( int )( colordata[vno*4+1] * ( 1 - opacity ) + m_rgbColor[1] * 255 * opacity );
      colordata[vno*4+2]  = ( int )( colordata[vno*4+2] * ( 1 - opacity ) + m_rgbColor[2] * 255 * opacity );
    }
  }
}

void SurfaceLabel::SetShowOutline(bool bOutline)
{
  m_bShowOutline = bOutline;
}

void SurfaceLabel::SetVisible(bool flag)
{
  m_bVisible = flag;
  emit SurfaceLabelChanged();
}
