/**
 * @brief Implementation for surface annotation.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


#include "SurfaceAnnotation.h"
#include "vtkLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkMath.h"
#include "LayerSurface.h"
#include "FSSurface.h"
#include <QFileInfo>
#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <QDebug>
#include "mri.h"

SurfaceAnnotation::SurfaceAnnotation ( LayerSurface* surf ) :
  QObject( surf ),
  m_nIndices( NULL ),
  m_nOutlineIndices(NULL),
  m_nCenterVertices( NULL ),
  m_lut( NULL ),
  m_surface( surf ),
  m_bShowOutline(false),
  m_dOpacity(1.0),
  m_nHighlightedLabel(-1),
  m_data(NULL)
{
  m_nIndexSize = surf->GetNumberOfVertices();
  m_data = new int[m_nIndexSize];
  m_nIndices = new int[m_nIndexSize];
  m_nOutlineIndices = new int[m_nIndexSize];
}

SurfaceAnnotation::~SurfaceAnnotation ()
{
  Reset();
}

void SurfaceAnnotation::Reset()
{
  if ( m_nIndices )
  {
    delete[] m_nIndices;
  }

  if ( m_nCenterVertices )
  {
    delete[] m_nCenterVertices;
  }

  if ( m_nOutlineIndices )
  {
    delete[] m_nOutlineIndices;
  }

  if (m_data)
    delete[] m_data;

  if (m_lut)
    CTABfree(&m_lut);

  m_nIndices = NULL;
  m_nOutlineIndices = NULL;
  m_lut = NULL;
  m_nCenterVertices = NULL;
  m_data = NULL;
  m_lut = NULL;

  for (int i = 0; i < m_bufferUndo.size(); i++)
    m_bufferUndo[i].Free();
  for (int i = 0; i < m_bufferRedo.size(); i++)
    m_bufferRedo[i].Free();
  m_bufferUndo.clear();
  m_bufferRedo.clear();
}

bool SurfaceAnnotation::LoadAnnotation( const QString& fn)
{
  if ( m_surface )
  {
    m_strFilename = QFileInfo(fn).canonicalFilePath();
    MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();

    int ret;
    try {
      ret = MRISreadAnnotation( mris, fn.toLatin1().data() );
    }
    catch (int return_code)
    {
      ret = return_code;
    }
    if ( ret != 0 )
    {
      cerr << "Could not load annotation from file " << qPrintable(fn) << ".\n";
      return false;
    }
    else
    {
      m_lut = CTABdeepCopy(mris->ct);
      UpdateColorList();
      for (int i = 0; i < m_nIndexSize; i++)
        m_data[i] = mris->vertices[i].annotation;
      UpdateData();
      SetSelectAllLabels();
      return true;
    }
  }
  return false;
}

bool SurfaceAnnotation::LoadFromSegmentation(const QString &fn)
{
   MRI* mri = MRIread(fn.toLatin1().data());
   if (mri->width != m_nIndexSize)
   {
       cerr << "Cannot load segmentation file. Wrong dimension size.";
       MRIfree(&mri);
       return false;
   }
   for (int i = 0; i < m_nIndexSize; i++)
   {
       int n = -1;
       switch ( mri->type )
       {
       case MRI_UCHAR:
         n = MRIseq_vox( mri, i, 0, 0, 0);
         break;
       case MRI_INT:
         n = MRIIseq_vox( mri, i, 0, 0, 0);
         break;
       case MRI_LONG:
         n = (int)MRILseq_vox( mri, i, 0, 0, 0);
         break;
       case MRI_FLOAT:
         n = (int)MRIFseq_vox( mri, i, 0, 0, 0);
         break;
       case MRI_SHORT:
         n = MRIseq_vox( mri, i, 0, 0, 0);
         break;
       default:
         break;
       }
       int annot = -1;
       if (CTABannotationAtIndex(m_lut, n, &annot) == 0)
           m_data[i] = annot;
   }
   UpdateData();
   SetSelectAllLabels();
   MRIfree(&mri);
   return true;
}

bool SurfaceAnnotation::LoadColorTable(const QString &fn)
{
    COLOR_TABLE* lut = CTABreadASCII( fn.toLatin1().data() );
    if (lut)
    {
        if (m_lut)
            CTABfree(&m_lut);
        m_lut = lut;
        UpdateColorList();
        MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
        if (mris->ct)
            CTABfree(&mris->ct);
        mris->ct = CTABdeepCopy(m_lut);
        for (int i = 0; i < m_nIndexSize; i++)
          m_data[i] = mris->vertices[i].annotation;
        UpdateData();
    }
    return (lut != 0);
}

bool SurfaceAnnotation::InitializeNewAnnotation(const QString& ct_fn)
{
  if ( m_surface )
  {
    if (ct_fn.isEmpty())
      m_lut = CTABalloc(0);
    else
      m_lut = CTABreadASCII(ct_fn.toLatin1().data());
    if (m_lut)
    {
      UpdateColorList();
      for (int i = 0; i < m_nIndexSize; i++)
        m_data[i] = -1;
      UpdateData();

      SetSelectAllLabels();
      return true;
    }
    else
    {
      cerr << "Could not load color table file" << qPrintable(ct_fn) << ".\n";
    }
  }
  return false;
}

int SurfaceAnnotation::ColorToAnnotation(const QColor &c)
{
  return CTABrgb2Annotation(c.red(), c.green(), c.blue());
}

void SurfaceAnnotation::UpdateData()
{
  MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();

  // find valid annotations
  QList<int> annotIndices;
  for ( int i = 0; i < m_nIndexSize; i++ )
  {
    int n = -1;
    if ( m_data[i] >= 0 && CTABfindAnnotation( m_lut, m_data[i], &n ) == 0 &&
         n >= 0)
    {
      m_nIndices[i] = n;
      if (!annotIndices.contains(n))
        annotIndices << n;
    }
    else
    {
      QList<int> new_ids = m_mapNewLabels.keys();
      bool bFound = false;
      foreach (int id, new_ids)
      {
        NewAnnotationLabel nl = m_mapNewLabels[id];
        if (ColorToAnnotation(nl.color) == m_data[i])
        {
          m_nIndices[i] = nl.id;
          if (!annotIndices.contains(nl.id))
            annotIndices << nl.id;
          bFound = true;
          break;
        }
      }
      if (!bFound)
        m_nIndices[i] = -1;
    }
  }
  qSort(annotIndices);
  QList<int> new_ids = m_mapNewLabels.keys();
  foreach (int id, new_ids)
  {
    if (!annotIndices.contains(id))
      m_mapNewLabels.remove(id);
  }

  int nAnnots = annotIndices.size();
  //   CTABgetNumberOfValidEntries( m_lut, &m_nAnnotations );
  if (m_nCenterVertices)
    delete[] m_nCenterVertices;
  m_nCenterVertices = new int[nAnnots];
  for ( int i = 0; i < nAnnots; i++ )
  {
    m_nCenterVertices[i] = -1;
  }

  // convert annotations to lookup table indices
  // also find the "center vertex"
  double** pts = new double*[nAnnots];
  int* vcount = new int[nAnnots];
  memset( vcount, 0, sizeof( int )*nAnnots );
  for ( int i = 0; i < nAnnots; i++ )
  {
    pts[i] = new double[3];
    memset( pts[i], 0, sizeof( double ) * 3 );
  }

  for ( int i = 0; i < m_nIndexSize; i++ )
  {
    if ( m_nIndices[i] != -1 && mris->vertices[i].ripflag == 0 )
    {
      int n = annotIndices.indexOf(m_nIndices[i]);
      vcount[n]++;
      pts[n][0] += mris->vertices[i].x;
      pts[n][1] += mris->vertices[i].y;
      pts[n][2] += mris->vertices[i].z;
      m_nCenterVertices[n] = i;
    }
  }

  for ( int i = 0; i < nAnnots; i++ )
  {
    if ( vcount[i] > 0 )
    {
      pts[i][0] /= vcount[i];
      pts[i][1] /= vcount[i];
      pts[i][2] /= vcount[i];

      int nVertex = m_surface->GetSourceSurface()->FindVertexAtSurfaceRAS( pts[i], NULL );
      if ( nVertex >=0 && m_nIndices[nVertex] == i )
      {
        m_nCenterVertices[i] = nVertex;
      }
    }
  }

  delete[] vcount;
  for ( int i = 0; i < nAnnots; i++ )
  {
    delete[] pts[i];
  }
  delete[] pts;

  // build outline indices
  memcpy(m_nOutlineIndices, m_nIndices, sizeof(int)*m_nIndexSize);
  COLOR_TABLE* old_ct = mris->ct;
  mris->ct = m_lut;
  for (int i = 0; i < nAnnots; i++)
  {
    VERTEX *v;
    VERTEX_TOPOLOGY* vt;
    MRISclearMarks(mris);
    LABEL* label = NULL;
    int nv = 0;
    for (int vno = 0; vno < m_nIndexSize; vno++)
    {
      if (m_nIndices[vno] == annotIndices[i])
        nv++;
    }
    if (nv > 0)
    {
      label = LabelAlloc(nv, NULL, "temp");
      nv = 0;
      for (int vno = 0; vno < mris->nvertices; vno++) {
        v = &mris->vertices[vno];
        if (v->ripflag)
          continue;
        if (m_nIndices[vno] == annotIndices[i])
        {
          label->lv[nv].vno = vno;
          label->n_points++;
          nv++;
        }
      }
    }

    if (label)
    {
      LabelMarkSurface(label, mris);
      for (int n = 0 ; n < label->n_points ; n++)
      {
        if (label->lv[n].vno >= 0)
        {
          m_nOutlineIndices[label->lv[n].vno] = -1;
          v = &mris->vertices[label->lv[n].vno] ;
          vt = &mris->vertices_topology[label->lv[n].vno] ;
          if (v->ripflag)
            continue;

          for (int m = 0 ; m < vt->vnum ; m++)
          {
            if (mris->vertices[vt->v[m]].marked == 0)
            {
              m_nOutlineIndices[label->lv[n].vno] = m_nIndices[label->lv[n].vno];
              break;
            }
          }
        }
      }
      LabelFree(&label);
    }
  }
  mris->ct = old_ct;

  m_listAnnotations = annotIndices;
}

void SurfaceAnnotation::GetAnnotationPoint( int nIndex, double* pt_out )
{
  m_surface->GetTargetAtVertex( m_nCenterVertices[nIndex], pt_out );
}

QString SurfaceAnnotation::GetName()
{
  return m_strName;
}

void SurfaceAnnotation::SetName( const QString& name )
{
  m_strName = name;
}

int SurfaceAnnotation::GetIndexAtVertex( int nVertex )
{
  return m_nIndices[nVertex];
}

QString SurfaceAnnotation::GetAnnotationNameAtVertex( int nVertex )
{
  return GetAnnotationNameAtIndex( m_nIndices[nVertex] );
}

QString SurfaceAnnotation::GetAnnotationNameAtIndex( int nIndex )
{
  char name[128];
  int nValid = 0;
  int nTotalCount = 0;
  if (nIndex >= UNASSIGNED_ANNOT_BASE)
    return m_mapNewLabels.value(nIndex).name;
  CTABgetNumberOfTotalEntries( m_lut, &nTotalCount );
  if ( nIndex >= 0 && nIndex < nTotalCount )
  {
    CTABisEntryValid( m_lut, nIndex, &nValid );
  }
  if ( nValid && CTABcopyName( m_lut, nIndex, name, 128 ) == 0 )
  {
    return name;
  }

  return "";
}

void SurfaceAnnotation::GetAnnotationColorAtIndex( int nIndex, int* rgb )
{
  int nValid = 0;
  int nTotalCount = 0;
  if (nIndex >= UNASSIGNED_ANNOT_BASE)
  {
    QColor c = m_mapNewLabels.value(nIndex).color;
    rgb[0] = c.red();
    rgb[1] = c.green();
    rgb[2] = c.blue();
    return;
  }
  CTABgetNumberOfTotalEntries( m_lut, &nTotalCount );
  if ( nIndex < nTotalCount )
  {
    CTABisEntryValid( m_lut, nIndex, &nValid );
  }
  if ( nValid )
  {
    CTABrgbAtIndexi( m_lut, nIndex, rgb, rgb+1, rgb+2 );
  }
}

void SurfaceAnnotation::SetShowOutline(bool bOutline)
{
  this->m_bShowOutline = bOutline;
}

void SurfaceAnnotation::MapAnnotationColor( unsigned char* colordata )
{
  int c[4];
  int* indices = (m_bShowOutline ? m_nOutlineIndices : m_nIndices);
  for ( int i = 0; i < m_nIndexSize; i++ )
  {
    if (indices[i] >= 0 && m_listVisibleLabels.contains(indices[i]))
    {
      if (indices[i] >= UNASSIGNED_ANNOT_BASE)
      {
        NewAnnotationLabel nl = m_mapNewLabels[indices[i]];
        colordata[i*4] = ( int )( colordata[i*4] * ( 1 - m_dOpacity ) + nl.color.red() * m_dOpacity );
        colordata[i*4+1] = ( int )( colordata[i*4+1] * ( 1 - m_dOpacity ) + nl.color.green() * m_dOpacity );
        colordata[i*4+2] = ( int )( colordata[i*4+2] * ( 1 - m_dOpacity ) + nl.color.blue() * m_dOpacity );
      }
      else
      {
        char name[128] = {0};
        if ( CTABcopyName( m_lut, indices[i], name, 128 ) == 0 && QString(name).toLower() != "unknown" )
        {
          if (CTABrgbAtIndexi( m_lut, indices[i], c, c+1, c+2 ) == 0) // no error & no black color
          {
            colordata[i*4] = ( int )( colordata[i*4] * ( 1 - m_dOpacity ) + c[0] * m_dOpacity );
            colordata[i*4+1] = ( int )( colordata[i*4+1] * ( 1 - m_dOpacity ) + c[1] * m_dOpacity );
            colordata[i*4+2] = ( int )( colordata[i*4+2] * ( 1 - m_dOpacity ) + c[2] * m_dOpacity );
          }
        }
      }
    }
  }

  if (m_nHighlightedLabel >= 0)
  {
    int* indices = m_nOutlineIndices;
    c[0] = 255; c[1] = 255; c[2] = 255;
    for (int i = 0; i < m_nIndexSize; i++)
    {
      if (indices[i] == m_nHighlightedLabel) // no error & no black color
      {
        colordata[i*4] = ( int )( colordata[i*4] * ( 1 - m_dOpacity ) + c[0] * m_dOpacity );
        colordata[i*4+1] = ( int )( colordata[i*4+1] * ( 1 - m_dOpacity ) + c[1] * m_dOpacity );
        colordata[i*4+2] = ( int )( colordata[i*4+2] * ( 1 - m_dOpacity ) + c[2] * m_dOpacity );
      }
    }
  }
}

void SurfaceAnnotation::SetSelectAllLabels()
{
  m_listVisibleLabels.clear();
  COLOR_TABLE* ct = GetColorTable();
  if (ct)
  {
    int cEntries;
    CTABgetNumberOfTotalEntries( ct, &cEntries );
    for ( int nEntry = 0; nEntry < cEntries; nEntry++ )
    {
      int bValid;
      CTABisEntryValid( ct, nEntry, &bValid );
      if ( bValid )
      {
        m_listVisibleLabels << nEntry;
      }
    }
    m_listVisibleLabels << m_mapNewLabels.keys();
  }
  SetModified();
}

void SurfaceAnnotation::SetUnselectAllLabels()
{
  m_listVisibleLabels.clear();
  SetModified();
}

void SurfaceAnnotation::SetSelectLabel(int nVal, bool bSelected)
{
  if (bSelected && !m_listVisibleLabels.contains(nVal))
    m_listVisibleLabels << nVal;
  else if (!bSelected)
    m_listVisibleLabels.removeOne(nVal);

  SetModified();
}

void SurfaceAnnotation::SetHighlightedLabel(int n)
{
  m_nHighlightedLabel = n;
  SetModified();
}

void SurfaceAnnotation::EditLabel(const QVector<int> &verts, int fill_index, const QVariantMap &options)
{
  SaveForUndo();
  int fill_annot = -1;
  if (options["CreateLabel"].toBool())
  {
    NewAnnotationLabel nl;
    if (m_mapNewLabels.isEmpty())
      nl.id = UNASSIGNED_ANNOT_BASE;
    else
    {
      QList<int> ids = m_mapNewLabels.keys();
      nl.id = m_mapNewLabels[ids.last()].id+1;
    }
    nl.color = GenerateNewColor();
    nl.name = QString("Unnamed %1").arg(nl.id-UNASSIGNED_ANNOT_BASE);
    fill_index = nl.id;
    if (!m_listVisibleLabels.contains(fill_index))
      m_listVisibleLabels << fill_index;
    m_mapNewLabels[nl.id] = nl;
    fill_annot = ColorToAnnotation(nl.color);
    m_listColors << nl.color;
  }
  setProperty("current_fill_index", fill_index);

  if (options.value("RemoveFromLabel").toBool())
  {
    foreach (int i, verts)
    {
      m_data[i] = -1;
    }
  }
  else
  {
    if (fill_index >= 0 && fill_index < UNASSIGNED_ANNOT_BASE)
      CTABannotationAtIndex(m_lut, fill_index, &fill_annot);
    foreach (int i, verts)
    {
      m_data[i] = fill_annot;
    }
  }

  UpdateData();
  SetModified();
}

QColor SurfaceAnnotation::GenerateNewColor()
{
  QColor c;
  while (!c.isValid() || (c.red()+c.green()+c.blue())/3 < 50 ||
         m_listColors.contains(c))
  {
    c = QColor(qrand()%256, qrand()%256, qrand()%256);
  }
  return c;
}

void SurfaceAnnotation::UpdateColorTable(int nIndex, const QString &name, const QColor &color)
{
  CTE* cte;
  if (nIndex < m_lut->nentries)
  {
    cte = m_lut->entries[nIndex];
    if (!cte)
    {
      cte = (CTE *)malloc(sizeof(CTE));
      m_lut->entries[nIndex] = cte;
    }
  }
  else
  {
    COLOR_TABLE* ct = CTABalloc(nIndex+1);
    for (int i = 0; i < m_lut->nentries; i++)
    {
      if (m_lut->entries[i])
        memmove(ct->entries[i], m_lut->entries[i], sizeof(COLOR_TABLE_ENTRY));
      else
      {
        free(ct->entries[i]);
        ct->entries[i] = NULL;
      }
    }
    for (int i = m_lut->nentries; i < nIndex; i++)
    {
      free(ct->entries[i]);
      ct->entries[i] = NULL;
    }
    CTABfree(&m_lut);
    m_lut = ct;
    cte = m_lut->entries[nIndex];
  }
  memset(cte->name, 0, STRLEN);
  sprintf(cte->name, "%s", qPrintable(name));
  cte->ri = color.red();
  cte->gi = color.green();
  cte->bi = color.blue();
  cte->rf = color.redF();
  cte->gf = color.greenF();
  cte->bf = color.blueF();
}

void SurfaceAnnotation::ReassignNewLabel(int nId, int ctab_id, const QString& name, const QColor& color)
{
  SaveForUndo();

  int nValid = 0;
  if (ctab_id >= 0 && ctab_id < m_lut->nentries)
    CTABisEntryValid( m_lut, ctab_id, &nValid );
  if (!nValid)
  {
    if (name.isEmpty() || !color.isValid())
    {
      cerr << "Invalid index in color table: " << ctab_id << endl;
      return;
    }
    UpdateColorTable(ctab_id, name, color);
  }

  int src_annot;  
  if (nId >= UNASSIGNED_ANNOT_BASE)
  {
    QColor color = m_mapNewLabels[nId].color;
    m_mapNewLabels.remove(nId);
    src_annot = ColorToAnnotation(color);
  }
  else
  {
     CTABannotationAtIndex(m_lut, nId, &src_annot);
  }
  m_listVisibleLabels.removeOne(nId);
  if (!m_listVisibleLabels.contains(ctab_id))
    m_listVisibleLabels << ctab_id;
  int annot = -1;
  CTABannotationAtIndex(m_lut, ctab_id, &annot);
  for (int i = 0; i < m_nIndexSize; i++)
  {
    if (m_data[i] == src_annot)
      m_data[i] = annot;
  }
  if (m_nHighlightedLabel == nId)
    m_nHighlightedLabel = ctab_id;
  UpdateData();
  SetModified();
}

void SurfaceAnnotation::DeleteLabel(int nId)
{
  SaveForUndo();
  int src_annot;
  if (nId >= UNASSIGNED_ANNOT_BASE)
  {
    QColor color = m_mapNewLabels[nId].color;
    m_mapNewLabels.remove(nId);
    src_annot = ColorToAnnotation(color);
  }
  else
  {
     CTABannotationAtIndex(m_lut, nId, &src_annot);
  }
  m_listVisibleLabels.removeOne(nId);
  for (int i = 0; i < m_nIndexSize; i++)
  {
    if (m_data[i] == src_annot)
      m_data[i] = -1;
  }
  if (m_nHighlightedLabel == nId)
    m_nHighlightedLabel = -1;
  UpdateData();
  SetModified();
}

void SurfaceAnnotation::UpdateColorList()
{
  m_listColors.clear();
  int nTotalCount = 0;
  CTABgetNumberOfTotalEntries( m_lut, &nTotalCount );
  int nValid = 0;
  for ( int i = 0; i < nTotalCount; i++ )
  {
    CTABisEntryValid( m_lut, i, &nValid );
    if ( nValid )
    {
      int nr, ng, nb;
      CTABrgbAtIndexi( m_lut, i, &nr, &ng, &nb );
      m_listColors << QColor( nr, ng, nb );
    }
  }
}

void SurfaceAnnotation::UpdateLabelInfo(int i, const QString &name, const QColor &color)
{
  SaveForUndo();

  QColor old_color;
  if (i < UNASSIGNED_ANNOT_BASE)
  {
    if (!name.isEmpty())
    {
      memset(m_lut->entries[i]->name, 0, STRLEN);
      strncpy(m_lut->entries[i]->name, name.toLatin1().data(), name.size());
    }
    if (color.isValid())
    {
      COLOR_TABLE_ENTRY* cte = m_lut->entries[i];
      old_color = QColor(cte->ri, cte->gi, cte->bi);
      m_listColors.removeOne(old_color);
      cte->ri = color.red();
      cte->gi = color.green();
      cte->bi = color.blue();
      cte->rf = color.redF();
      cte->gf = color.greenF();
      cte->bf = color.blueF();
      m_listColors << color;
    }
  }
  else if (m_mapNewLabels.contains(i))
  {
    NewAnnotationLabel nl = m_mapNewLabels[i];
    if (!name.isEmpty())
      nl.name = name;
    if (color.isValid())
    {
      old_color = nl.color;
      m_listColors.removeOne(nl.color);
      nl.color = color;
      m_listColors << color;
    }
    m_mapNewLabels[i] = nl;
  }

  if (color.isValid())
  {
    int old_annot = ColorToAnnotation(old_color);
    int new_annot = ColorToAnnotation(color);
    for (int i = 0; i < m_nIndexSize; i++)
    {
        if (m_data[i] == old_annot)
          m_data[i] = new_annot;
    }
  }

  if (color.isValid())
    SetModified();
}

bool SurfaceAnnotation::SaveToFile(const QString &fn)
{
  MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
  if (mris->ct)
    CTABfree(&mris->ct);

  mris->ct = CTABdeepCopy(m_lut);
  for (int i = 0; i < m_nIndexSize; i++)
    mris->vertices[i].annotation = m_data[i];

  bool ret = (::MRISwriteAnnotation(mris, qPrintable(fn)) == 0);
  if (ret)
    m_strFilename = fn;
  return ret;
}

bool SurfaceAnnotation::HasUndo()
{
  return m_bufferUndo.size() > 0;
}

bool SurfaceAnnotation::HasRedo()
{
  return m_bufferRedo.size() > 0;
}

AnnotUndoRedoBufferItem SurfaceAnnotation::SaveCurrentUndoRedoBuffer()
{
  AnnotUndoRedoBufferItem item;
  item.m_data = new int[m_nIndexSize];
  memcpy(item.m_data, m_data, sizeof(int)*m_nIndexSize);
  item.m_listVisibleLabels = m_listVisibleLabels;
  item.m_mapNewLabels = m_mapNewLabels;
  item.m_ctab = CTABdeepCopy(m_lut);
  return item;
}

void SurfaceAnnotation::RestoreFromUndoRedoBuffer(const AnnotUndoRedoBufferItem &item)
{
  memcpy(m_data, item.m_data, sizeof(int)*m_nIndexSize);
  m_listVisibleLabels = item.m_listVisibleLabels;
  m_mapNewLabels = item.m_mapNewLabels;
  CTABfree(&m_lut);
  m_lut = CTABdeepCopy(item.m_ctab);
  UpdateData();
  SetModified();
}

void SurfaceAnnotation::Undo()
{
  if ( m_bufferUndo.size() > 0 )
  {
    AnnotUndoRedoBufferItem item = m_bufferUndo[m_bufferUndo.size()-1];
    m_bufferUndo.pop_back();

    m_bufferRedo.push_back( SaveCurrentUndoRedoBuffer() );

    RestoreFromUndoRedoBuffer( item );
    item.Free();

    SetModified();
//    emit ActorUpdated();
    emit Modified();
  }
}

void SurfaceAnnotation::Redo()
{
  if ( m_bufferRedo.size() > 0 )
  {
    AnnotUndoRedoBufferItem item = m_bufferRedo[m_bufferRedo.size()-1];
    m_bufferRedo.pop_back();

    m_bufferUndo.push_back( SaveCurrentUndoRedoBuffer() );

    RestoreFromUndoRedoBuffer( item );
    item.Free();

    SetModified();
//    emit ActorUpdated();
    emit Modified();
  }
}

void SurfaceAnnotation::SaveForUndo()
{
  m_bufferUndo.push_back( SaveCurrentUndoRedoBuffer() );

  // clear redo buffer
  for ( size_t i = 0; i < m_bufferRedo.size(); i++ )
  {
    m_bufferRedo[i].Free();
  }
  m_bufferRedo.clear();
}

void SurfaceAnnotation::CleanUpColorTable()
{
  SaveForUndo();

  for (int i = 0; i < m_lut->nentries; i++)
  {
    CTE* cte = m_lut->entries[i];
    if (cte)
    {
      if (!m_listAnnotations.contains(i))
      {
        free(cte);
        m_lut->entries[i] = NULL;
      }
    }
  }
}
