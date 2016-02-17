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
 *    $Author: zkaufman $
 *    $Date: 2016/02/17 20:36:46 $
 *    $Revision: 1.21 $
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
#include "vtkRGBAColorTransferFunction.h"

SurfaceLabel::SurfaceLabel ( LayerSurface* surf ) :
  QObject( surf ),
  m_label( NULL ),
  m_surface( surf ),
  m_bVisible( true ),
  m_nOutlineIndices(NULL),
  m_dThreshold(0),
  m_nColorCode(SolidColor),
  m_dHeatscaleMin(0),
  m_dHeatscaleMax(1)
{
  m_rgbColor[0] = 1.0;
  m_rgbColor[1] = 1.0;
  m_rgbColor[2] = 0.0;

  SetShowOutline(false);
  SetColor(1.0, 1.0, 0);
  m_lut = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  UpdateLut();
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

void SurfaceLabel::UpdateLut()
{
  m_lut->RemoveAllPoints();
  m_lut->AddRGBAPoint( -m_dHeatscaleMax, 0, 1, 1, 1 );
  m_lut->AddRGBAPoint( -m_dHeatscaleMin, 0, 0, 1, 1 );
  m_lut->AddRGBAPoint(  0, 0, 0, 0, 0 );
  m_lut->AddRGBAPoint(  m_dHeatscaleMin, 1, 0, 0, 1 );
  m_lut->AddRGBAPoint(  m_dHeatscaleMax, 1, 1, 0, 1 );
  m_lut->Build();
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

  if (m_label->n_points > 0)
    m_dHeatscaleMin = m_dHeatscaleMax = m_label->lv[0].stat;
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
    if (m_label->lv[i].stat < m_dHeatscaleMin)
      m_dHeatscaleMin = m_label->lv[i].stat;
    else if (m_label->lv[i].stat > m_dHeatscaleMax)
      m_dHeatscaleMax = m_label->lv[i].stat;
  }

  // create outline
  m_nOutlineIndices = new int[m_label->n_points];
  UpdateOutline();
  UpdateLut();

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

  emit SurfaceLabelChanged();
}

void SurfaceLabel::SetThreshold(double th)
{
  if (th != m_dThreshold)
  {
    m_dThreshold = th;
    emit SurfaceLabelChanged();
  }
}

void SurfaceLabel::SetColorCode(int nCode)
{
  if (m_nColorCode != nCode)
  {
    m_nColorCode = nCode;
    emit SurfaceLabelChanged();
  }
}

void SurfaceLabel::SetHeatscaleMin(double dval)
{
  if (dval != m_dHeatscaleMin)
  {
    m_dHeatscaleMin = dval;
    UpdateLut();
    emit SurfaceLabelChanged();
  }
}

void SurfaceLabel::SetHeatscaleMax(double dval)
{
  if (dval != m_dHeatscaleMax)
  {
    m_dHeatscaleMax = dval;
    UpdateLut();
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
      double opacity = 1;
      if (m_label->lv[i].stat >= m_dThreshold)
        opacity = 1;
      else
        opacity = 0;
      double rgb[4] = { m_rgbColor[0], m_rgbColor[1], m_rgbColor[2], 1 };
      if (m_nColorCode == Heatscale)
      {
        m_lut->GetColor(m_label->lv[i].stat, rgb);
      }
      colordata[vno*4]    = ( int )( colordata[vno*4]   * ( 1 - opacity ) + rgb[0] * 255 * opacity );
      colordata[vno*4+1]  = ( int )( colordata[vno*4+1] * ( 1 - opacity ) + rgb[1] * 255 * opacity );
      colordata[vno*4+2]  = ( int )( colordata[vno*4+2] * ( 1 - opacity ) + rgb[2] * 255 * opacity );
    }
  }
}

void SurfaceLabel::SetShowOutline(bool bOutline)
{
  m_bShowOutline = bOutline;
  emit SurfaceLabelChanged();
}

void SurfaceLabel::SetVisible(bool flag)
{
  m_bVisible = flag;
  emit SurfaceLabelVisibilityChanged();
}

bool SurfaceLabel::GetCentroid(double *x, double *y, double *z, int *nvo)
{
  return (LabelCentroid(m_label, m_surface->GetSourceSurface()->GetMRIS(), x, y, z, nvo) == 0);
}
