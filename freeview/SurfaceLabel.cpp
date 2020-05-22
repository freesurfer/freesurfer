/**
 * @brief Implementation for surface label.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Ruopeng Wang
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
#include "LayerMRI.h"
#include "FSVolume.h"

SurfaceLabel::SurfaceLabel ( LayerSurface* surf, bool bInitializeLabel ) :
  QObject( surf ),
  m_label( NULL ),
  m_surface( surf ),
  m_bVisible( true ),
  m_nOutlineIndices(NULL),
  m_dThreshold(0),
  m_nColorCode(SolidColor),
  m_dHeatscaleMin(0),
  m_dHeatscaleMax(1),
  m_dOpacity(1.0),
  m_bModified(false)
{
  m_rgbColor[0] = 1.0;
  m_rgbColor[1] = 1.0;
  m_rgbColor[2] = 0.0;

  SetShowOutline(false);
  SetColor(1.0, 1.0, 0);
  m_lut = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  UpdateLut();

  if (bInitializeLabel)
  {
    MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
    m_label = ::LabelAlloc( mris->nvertices, NULL, (char*)"" );
    int unassigned;
    LabelIsCompletelyUnassigned(m_label, &unassigned);
    if (unassigned)
    {
      LabelFillUnassignedVertices(mris, m_label, CURRENT_VERTICES);
    }
    m_label->vertex_label_ind = (int *)calloc(mris->nvertices, sizeof(int));
    for (int n = 0; n < mris->nvertices; n++)
      m_label->vertex_label_ind[n] = -1;

    m_nOutlineIndices = new int[mris->nvertices];
    UpdateOutline();
  }
}

SurfaceLabel::~SurfaceLabel ()
{
  if ( m_label )
    ::LabelFree( &m_label );

  foreach (LABEL* l, m_undoBuffer)
    ::LabelFree(&l);
  foreach (LABEL* l, m_redoBuffer)
    ::LabelFree(&l);

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
  int unassigned ;

  if ( m_label )
  {
    ::LabelFree( &m_label );
  }

  m_label = ::LabelRead( NULL, filename.toLatin1().data() );
  if ( m_label == NULL )
  {
    cerr << "LabelRead failed";
    return false;
  }

  LabelIsCompletelyUnassigned(m_label, &unassigned) ;
  if (unassigned)
  {
    //    LabelFree(&m_label) ;
    //    cerr << "label has not been mapped to surface";
    //    return false;
    LabelFillUnassignedVertices(m_surface->GetSourceSurface()->GetMRIS(), m_label, CURRENT_VERTICES);
    cout << "label assigned to surface";
  }

  MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
  m_label->vertex_label_ind = (int *)calloc(mris->nvertices, sizeof(int));
  for (int n = 0; n < mris->nvertices; n++)
    m_label->vertex_label_ind[n] = -1;

  for (int n = 0; n < m_label->n_points; n++)
  {
    LV* lv = &m_label->lv[n];

    if (lv->deleted) continue;
    if (lv->vno >= 0 && lv->vno <= mris->nvertices)  // vertex already assigned
    {
      m_label->vertex_label_ind[lv->vno] = n;
    }
  }

  QFile file( filename );
  if ( !file.open( QIODevice::ReadOnly | QIODevice::Text ) )
  {
    return false;
  }
  m_strFilename = filename;

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
  double max_spacing;
  int max_vno;
  MRIScomputeVertexSpacingStats(mris, NULL, NULL, &max_spacing, NULL, &max_vno, CURRENT_VERTICES);
  return true;
  MHT* hash = MHTcreateVertexTable_Resolution(mris,
                                              CURRENT_VERTICES, max_spacing/4);

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
      int vtxno = MHTfindClosestVertexNoXYZ(hash, mris, v.x,v.y,v.z, &dmin);
      if (vtxno >= 0)
        m_label->lv[i].vno = vtxno;
    }
    if (m_label->lv[i].stat < m_dHeatscaleMin)
      m_dHeatscaleMin = m_label->lv[i].stat;
    else if (m_label->lv[i].stat > m_dHeatscaleMax)
      m_dHeatscaleMax = m_label->lv[i].stat;
  }

  MHTfree(&hash);

  // create outline
  m_nOutlineIndices = new int[mris->nvertices];
  UpdateOutline();
  UpdateLut();

  return true;
}

void SurfaceLabel::UpdateOutline()
{
  VERTEX *v;
  VERTEX_TOPOLOGY* vt;
  MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
  MRISclearMarks(mris);
  LabelMarkSurface(m_label, mris);
  for (int n = 0 ; n < m_label->n_points ; n++)
  {
    m_nOutlineIndices[n] = 0;
    if (m_label->lv[n].vno >= 0 && !m_label->lv[n].deleted) // && m_label->lv[n].stat > m_dThreshold)
    {
      v = &mris->vertices[m_label->lv[n].vno];
      vt = &mris->vertices_topology[m_label->lv[n].vno];
      if (v->ripflag)
        continue;

      for (int m = 0 ; m < vt->vnum ; m++)
      {
        if (mris->vertices[vt->v[m]].marked == 0)
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

void SurfaceLabel::SetOpacity(double dval)
{
  if (dval != m_dOpacity)
  {
    m_dOpacity = dval;
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
    if (vno < nVertexCount && !m_label->lv[i].deleted && (!m_bShowOutline || m_nOutlineIndices[i] > 0) )
    {
      double opacity;
      if (m_label->lv[i].stat >= m_dThreshold)
        opacity = m_dOpacity;
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

void SurfaceLabel::Resample(LayerMRI *mri)
{
  ::LabelUnassign(m_label);
  int coord = m_surface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES;
  if (!m_label->mri_template)
    LabelInit(m_label, mri->GetSourceVolume()->GetMRI(), m_surface->GetSourceSurface()->GetMRIS(), coord);
  LABEL* label = ::LabelSampleToSurface(m_surface->GetSourceSurface()->GetMRIS(), m_label,
                                        m_label->mri_template, coord);
  if (label)
  {
    LabelCopy(label, m_label) ;
    LabelFree(&label);
    emit SurfaceLabelChanged();
  }
}

void SurfaceLabel::Dilate(int nTimes)
{
  ::LabelDilate(m_label, m_surface->GetSourceSurface()->GetMRIS(), nTimes,
                m_surface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
  emit SurfaceLabelChanged();
}

void SurfaceLabel::Erode(int nTimes)
{
  ::LabelErode(m_label, m_surface->GetSourceSurface()->GetMRIS(), nTimes);
  emit SurfaceLabelChanged();
}

void SurfaceLabel::Open(int nTimes)
{
  ::LabelErode(m_label, m_surface->GetSourceSurface()->GetMRIS(), nTimes);
  ::LabelDilate(m_label, m_surface->GetSourceSurface()->GetMRIS(), nTimes,
                m_surface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
  emit SurfaceLabelChanged();
}

void SurfaceLabel::Close(int nTimes)
{
  ::LabelDilate(m_label, m_surface->GetSourceSurface()->GetMRIS(), nTimes,
                m_surface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
  ::LabelErode(m_label, m_surface->GetSourceSurface()->GetMRIS(), nTimes);
  emit SurfaceLabelChanged();
}

bool SurfaceLabel::HasVertex(int nvo)
{
  for (int n = 0; n < m_label->n_points; n++)
  {
    if (m_label->lv[n].vno == nvo && !m_label->lv[n].deleted)
      return true;
  }
  return false;
}

void SurfaceLabel::EditVertices(const QVector<int> &verts, bool bAdd)
{
  MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
  LV *lv;
  VERTEX *v;
  double x, y, z;
  int coord = m_surface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES;
  if (bAdd)
  {
    for (int i = 0; i < verts.size(); i++)
    {
      int vno = verts[i];
      if (m_label->vertex_label_ind[vno] >= 0)
        continue;

      if (m_label->n_points >= m_label->max_points)
      {
        LabelRealloc(m_label, mris->nvertices);
//        qDebug() << "reallocated label";
      }

      int n = m_label->n_points++;  // n is the number of points before incr
      lv = &m_label->lv[n];
      v = &(mris->vertices[vno]);
      MRISgetCoords(v, coord, &x, &y, &z);

      lv->vno = vno;
      lv->x = x;
      lv->y = y;
      lv->z = z;
      lv->deleted = 0;
      m_label->vertex_label_ind[vno] = n;
    }
  }
  else
  {
    for (int i = 0; i < verts.size(); i++)
    {
      int vno = verts[i];
      int n = m_label->vertex_label_ind[vno];
      if (n < 0)
        continue;

      m_label->lv[n].deleted = 1;
      m_label->vertex_label_ind[vno] = -1;
    }

    // compress label data
    QVector<LV> list;
    for (int i = 0; i < m_label->n_points; i++)
    {
      if (!m_label->lv[i].deleted)
        list << m_label->lv[i];
    }
    for (int i = 0; i < list.size(); i++)
      m_label->lv[i] = list[i];
    m_label->n_points = list.size();
  }

  m_bModified = true;
  UpdateOutline();
  emit SurfaceLabelChanged();
}

bool SurfaceLabel::SaveToFile(const QString &filename)
{
  QString fn = filename;
  if (fn.isEmpty())
    fn = m_strFilename;

  if (fn.isEmpty())
  {
    cerr << "Can not write to empty filename\n";
    return false;
  }

  int err = ::LabelWrite( m_label, fn.toLatin1().data() );
  if ( err != 0 )
  {
    cerr << "LabelWrite failed\n";
    return false;
  }
  return true;
}

void SurfaceLabel::Undo()
{
  if (!m_undoBuffer.isEmpty())
  {
    LABEL* l = m_undoBuffer.last();
    LABEL* l2 = ::LabelCopy(m_label, NULL);
    ::LabelCopy(l, m_label);
    ::LabelFree(&l);
    m_undoBuffer.removeLast();
    m_redoBuffer << l2;
    UpdateOutline();
    emit SurfaceLabelChanged();
  }
}

void SurfaceLabel::Redo()
{
  if (!m_redoBuffer.isEmpty())
  {
    LABEL* l = m_redoBuffer.last();
    LABEL* l2 = ::LabelCopy(m_label, NULL);
    ::LabelCopy(l, m_label);
    ::LabelFree(&l);
    m_redoBuffer.removeLast();
    m_undoBuffer << l2;
    UpdateOutline();
    emit SurfaceLabelChanged();
  }
}

void SurfaceLabel::SaveForUndo()
{
  LABEL* l = ::LabelCopy(m_label, NULL);
  m_undoBuffer << l;

  // clear redo buffer
  foreach (LABEL* l, m_redoBuffer)
    ::LabelFree(&l);
  m_redoBuffer.clear();
}

void SurfaceLabel::MaskOverlay()
{
  ::LabelMaskSurface(m_label, m_surface->GetSourceSurface()->GetMRIS());
  UpdateOutline();
  emit SurfaceLabelChanged();
}
